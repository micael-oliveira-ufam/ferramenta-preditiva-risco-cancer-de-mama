# O modelo pode ser usado para inferência preditiva?

## Análise de aptidão, limites de confiabilidade e documentação da plataforma

---

## 1. Resposta curta

**Sim para ordenar risco; não para estimar risco individual; nunca para decisão clínica.**

Três usos precisam ser separados, porque a resposta é diferente para cada um:

| Uso pretendido | Apto? | Base da conclusão |
|---|---|---|
| **Ordenar pacientes por risco** (triagem para estudo, priorização de seguimento em pesquisa) | **Sim** | Log-rank entre tercis p = 1,8 × 10⁻³⁶; sobrevida mediana de 217,6 contra 85,7 meses entre os extremos |
| **Estimar a probabilidade individual de óbito em 10 anos** | **Com ressalvas** | Calibração adequada após correção isotônica (Hosmer-Lemeshow p = 0,66 no modelo clínico), mas sem validação externa |
| **Apoiar decisão clínica sobre uma paciente** | **Não** | Sem validação externa, sem randomização, coorte histórica pré-trastuzumabe, AUC de 0,72 |
| **Substituir plataformas validadas** (Oncotype, MammaPrint) | **Não** | Nenhuma comparação foi feita; os genes sequer superaram variáveis clínicas de rotina |

---

## 2. Uma correção importante em relação ao que reportei antes

Ao preparar a exportação dos modelos, detectei **vazamento de seleção de variáveis** na primeira
versão: a assinatura de 106 genes era escolhida usando toda a base e só depois submetida à
validação cruzada. Isso inflou artificialmente os resultados:

| Métrica | Com vazamento (incorreto) | Sem vazamento (correto) |
|---|---|---|
| AUC do modelo de genes | 0,727 | **0,639** |
| AUC do modelo combinado | 0,768 | **0,692** |
| Basal — AUC | 0,710 (IC 0,632–0,783) → "confiável" | **0,577 (IC 0,491–0,659) → não confiável** |

A versão corrigida faz a seleção de genes **e** a calibração dentro de cada partição de treino.
Os números da plataforma são os da coluna da direita, e voltam a concordar com o benchmark
anterior de onze algoritmos (modelo de genes ≈ 0,65). Registro isso explicitamente porque o
erro, se tivesse passado, apareceria justamente onde é mais perigoso: declarando confiável a
predição em tumores basais.

---

## 3. Em que aspectos a predição é confiável

### 3.1 Confiável — a ordenação de risco

Os tercis de risco previsto separam a sobrevida real de forma inequívoca:

| Grupo | n | Óbitos em 10 anos observados | Sobrevida global mediana |
|---|---|---|---|
| Baixo | 535 | 26,2% | 217,6 meses |
| Médio | 505 | 45,7% | 142,4 meses |
| Alto | 520 | 63,6% | 85,7 meses |

**Log-rank p = 1,8 × 10⁻³⁶.** Uma diferença de 132 meses entre os extremos. Para qualquer
tarefa que dependa de *ordem* — selecionar pacientes de maior risco para um protocolo, definir
estratos de um estudo, priorizar amostras para sequenciamento — o modelo cumpre a função.

### 3.2 Confiável, após correção — a calibração

Antes da correção isotônica, o modelo atribuía 6,8% de risco a um grupo onde morriam 25,1%
(Hosmer-Lemeshow p = 5,7 × 10⁻¹⁶³). Após a correção, ajustada dentro da validação cruzada:

| Modelo | Brier | Hosmer-Lemeshow |
|---|---|---|
| Clínico | 0,211 | **p = 0,66** |
| Genes | 0,233 | p = 0,25 |
| Combinado | 0,223 | p = 0,20 |

Os três passam no teste de calibração. Tabela observada contra prevista, no modelo em produção:

| Risco previsto | n | Previsto | Observado |
|---|---|---|---|
| ≤ 20% | 161 | 14,8% | 19,3% |
| 20–30% | 139 | 23,1% | 26,6% |
| 30–40% | 273 | 36,2% | 31,9% |
| 40–50% | 418 | 45,4% | 46,2% |
| 50–60% | 148 | 54,2% | 54,1% |
| 60–70% | 280 | 62,9% | 62,9% |
| > 70% | 141 | 76,9% | 69,5% |

A concordância é boa na faixa intermediária e piora nos extremos — o modelo ainda subestima
risco nos casos que julga muito favoráveis e superestima nos que julga muito graves.

### 3.3 Não confiável — dois estratos específicos

Testando cada estrato com IC 95% por bootstrap (1.500 reamostragens):

| Estrato | n | AUC | IC 95% | Veredito |
|---|---|---|---|---|
| Normal-like | 111 | 0,719 | 0,613 – 0,813 | confiável |
| claudin-low | 150 | 0,714 | 0,621 – 0,793 | confiável |
| LumA | 544 | 0,678 | 0,630 – 0,727 | confiável |
| LumB | 394 | 0,661 | 0,608 – 0,717 | confiável |
| Her2 | 184 | 0,608 | 0,522 – 0,688 | limítrofe |
| **Basal** | 177 | **0,577** | **0,491 – 0,659** | **não confiável** |
| Estágio I | 390 | 0,651 | 0,591 – 0,710 | confiável |
| Estágio II | 668 | 0,643 | 0,599 – 0,685 | confiável |
| **Estágio III–IV** | 109 | **0,506** | **0,387 – 0,626** | **não confiável** |

No subtipo Basal e no estágio III–IV o intervalo cruza 0,5: a predição não se distingue
estatisticamente de um sorteio. **A plataforma bloqueia a leitura do número nesses casos**, em
vez de exibi-lo com uma ressalva discreta.

### 3.4 Não confiável — o caso que mais importaria

Os óbitos que o modelo não antecipa são de pacientes **mais jovens (61,3 contra 69,1 anos),
com tumores menores (20 contra 30 mm), sem linfonodos comprometidos e com índice de Nottingham
favorável** (todos FDR < 0,001). É exatamente o nicho onde um marcador molecular teria valor
clínico — e é onde o modelo falha. Isso não é um defeito de implementação: é a consequência
direta de os genes não acrescentarem informação sobre as variáveis clínicas.

---

## 4. Barreiras estruturais ao uso preditivo real

1. **Sem validação externa.** Toda a estimativa é interna. Um modelo pode manter desempenho no
   próprio conjunto e desabar em outra população — é a regra, não a exceção.
2. **Coorte histórica.** Pacientes tratadas entre os anos 1980 e 2000, majoritariamente antes do
   trastuzumabe. O prognóstico do subtipo Her2 na base **não corresponde** ao da prática atual.
3. **Desenho observacional.** Sem randomização de tratamento; o modelo aprende também os padrões
   de indicação terapêutica da época.
4. **Dependência de escala.** Os z-scores são calculados sobre a própria coorte. Aplicar a uma
   paciente nova exige recalibrar contra uma referência — não basta padronizar isoladamente.
5. **Plataforma de medição.** Microarranjo Illumina HT-12; dados de RNA-seq não são
   intercambiáveis sem harmonização.
6. **AUC de 0,72 é modesto.** Comparável ao que o estadiamento clínico isolado já entrega.

---

## 5. A plataforma construída

Aplicação estática de página única, com todo o cálculo no navegador — nenhum dado do usuário
trafega pela rede. Os modelos foram exportados como coeficientes (61 KB de JSON), e a predição
é reimplementada em JavaScript: padronização, produto escalar, sigmoide e interpolação da curva
isotônica. Verifiquei que a implementação em JS reproduz os valores do scikit-learn.

**Decisões de projeto que decorrem da análise acima:**

- O **modelo clínico é o padrão**, não o genômico — seria incoerente destacar o modelo que a
  própria análise mostrou inferior. Ao selecionar o modelo de genes ou o combinado, a interface
  informa que ele discrimina pior e por quê.
- O **elemento central da tela é o painel de confiabilidade**, não o número do risco. Ele mostra
  a barra de AUC com IC 95% para cada subtipo, com a linha do acaso marcada, e emite um veredito
  explícito. Em Basal, a mensagem é literal: ignore o número acima.
- O **classificador de subtipo informa a própria confiança** e avisa quando ela está abaixo de
  0,6 — faixa em que acerta pouco mais que a metade das vezes.
- **Faixa de aviso permanente** no topo, com o alcance e os limites do uso.
- O **desfecho real dos casos-exemplo fica oculto** até o usuário pedir, para que a predição seja
  feita sem ancoragem.

### Dados de entrada de exemplo

Seis pacientes reais da coorte, com os 489 valores de expressão preenchidos, escolhidos para
cobrir condições clínicas distintas (`exemplos_entrada.csv` traz os mesmos casos no formato da
planilha de treino):

| Caso | Perfil | Desfecho real |
|---|---|---|
| Luminal A, estágio I, axila negativa | 56,8 anos, grau 1, 17 mm, NPI 2,03 | viva, 153,6 meses |
| Luminal B, estágio II, axila positiva | 77,0 anos, grau 3, 40 mm, 8 linfonodos, NPI 6,08 | óbito, 41,4 meses |
| Basal, grau 3 | subtipo em que o modelo não é confiável | óbito, 8,1 meses |
| HER2-enriquecido, tumor grande | coorte pré-trastuzumabe | óbito, 36,3 meses |
| Claudin-low, evolução favorável | assinatura imune/estromal | viva, 140,5 meses |
| **Luminal A jovem, aparência favorável** | **< 55 anos, 22 mm, axila negativa — e óbito** | **óbito, 32,6 meses** |

O último caso está incluído de propósito: é o modo de falha da seção 3.4. O modelo clínico
atribui a ele **11,0%** de risco, e a paciente morreu em 32,6 meses. O caso demonstra, na própria
interface, o limite que nenhuma métrica agregada comunica tão bem.

### Publicação no Cloudflare

Deixei prontos `wrangler.toml` (Worker com assets), `_headers` (CSP restritiva, `connect-src 'none'`,
sem enquadramento em iframe) e o passo a passo no `README.md`. O caminho mais curto:

```bash
npm install -g wrangler
wrangler login
cd plataforma_metabric
wrangler pages project create painel-metabric --production-branch main
wrangler pages deploy . --project-name painel-metabric
```

O `wrangler` devolve o link externo `https://painel-metabric.pages.dev` ao final.

**Não consegui executar a publicação por você:** este ambiente só alcança um conjunto restrito
de domínios (PyPI, npm, GitHub), e a API do Cloudflare não está entre eles. A autenticação
também exigiria sua conta. Os comandos acima rodam em qualquer máquina com Node instalado.

Dado que o material é de pesquisa e trata de prognóstico oncológico, sugiro publicar **atrás do
Cloudflare Access** (Zero Trust → Access → Applications), restringindo ao domínio institucional
da Fiocruz ou da UFAM. Uma URL `.pages.dev` aberta é indexável, e um painel de risco de câncer
sem controle de acesso tende a ser lido fora do contexto de pesquisa para o qual foi feito.

---

## 6. O que faria o modelo apto a mais do que ordenar risco

1. **Validação externa em TCGA-BRCA**, com atenção à diferença de plataforma. É o passo que
   separa um exercício metodológico de um modelo utilizável.
2. **Restringir ao contexto onde há sinal** — luminais em estágio II — e rederivar ali.
3. **Modelos de risco competitivo** (Fine-Gray), dada a alta proporção de óbitos por outras causas.
4. **Random Survival Forest**, que usa o tempo até o evento em vez de dicotomizar em 10 anos.
5. **Coorte contemporânea**, com terapias atuais, para que o prognóstico do Her2 faça sentido.
