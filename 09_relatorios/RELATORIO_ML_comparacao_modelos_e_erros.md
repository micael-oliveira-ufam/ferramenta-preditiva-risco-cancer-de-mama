# Comparação de algoritmos de aprendizado de máquina e anatomia dos erros

## Extensão do relatório METABRIC — 11 algoritmos, 36 combinações modelo × conjunto de variáveis, com validação estatística formal

**Coorte:** 1.897 pacientes (classificação de subtipo) e 1.560 pacientes com desfecho de 10 anos determinado (predição de sobrevida)
**Protocolo:** validação cruzada estratificada 5-fold repetida 3 vezes (15 partições por modelo), sementes fixas
**Testes empregados:** t pareado corrigido de Nadeau-Bengio, McNemar, DeLong, qui-quadrado, Mann-Whitney, binomial, Hosmer-Lemeshow, log-rank

---

## 1. Por que comparar vários algoritmos, e por que testar a comparação

Um único modelo com bom desempenho não diz se o resultado vem **do algoritmo** ou **dos dados**. Se onze algoritmos com pressupostos radicalmente diferentes — lineares, baseados em distância, em árvores, em margens, em redes neurais — chegam a desempenhos semelhantes, o teto observado é uma propriedade da informação disponível, não uma limitação de escolha técnica. Essa é a pergunta que esta seção responde.

Igualmente importante: **diferenças de desempenho entre modelos precisam de teste estatístico**. Comparar duas acurácias médias e declarar vencedor o maior número é um erro comum e sério. As estimativas vêm das mesmas partições de dados, são correlacionadas entre si, e um teste t comum aplicado a folds de validação cruzada infla drasticamente a taxa de falsos positivos. Por isso foi usado o **t corrigido de Nadeau-Bengio**, que ajusta a variância pela sobreposição entre conjuntos de treino, e o **McNemar**, que compara os modelos paciente a paciente em vez de média a média.

**Glossário rápido dos testes usados:**

| Teste | O que responde | Como ler |
|---|---|---|
| **t de Nadeau-Bengio** | Dois modelos diferem em desempenho médio? | Corrige a dependência entre folds; p < 0,05 = diferença real |
| **McNemar** | Dois modelos erram nos *mesmos* pacientes? | Compara acertos exclusivos de cada um; p alto = empate técnico |
| **DeLong** | Duas curvas ROC diferem? | Teste específico para AUC pareado no mesmo conjunto |
| **Qui-quadrado** | O erro se concentra em algum grupo? | p baixo = o erro não é aleatório entre grupos |
| **Mann-Whitney** | Casos errados diferem dos acertados em alguma variável? | Compara distribuições sem assumir normalidade |
| **Hosmer-Lemeshow** | O risco previsto corresponde ao observado? | p baixo = modelo mal calibrado |
| **Log-rank** | Grupos de risco previstos têm sobrevida real distinta? | p baixo = a estratificação funciona na prática |

---

## 2. Tarefa 1 — Classificação dos seis subtipos moleculares

### 2.1 Desempenho dos 11 algoritmos

Média de 15 partições (5-fold × 3 repetições), a partir dos 489 genes:

| Modelo | Acurácia | DP | Acurácia balanceada | F1-macro |
|---|---|---|---|---|
| **Gradient Boosting** | **0,7872** | 0,019 | 0,7256 | 0,7430 |
| SVM RBF | 0,7760 | 0,016 | 0,7241 | 0,7385 |
| Random Forest | 0,7663 | 0,021 | 0,6780 | 0,7021 |
| Rede neural (MLP) | 0,7591 | 0,017 | 0,7166 | 0,7266 |
| Regressão logística (L2) | 0,7582 | 0,021 | 0,7183 | 0,7265 |
| Extra Trees | 0,7563 | 0,023 | 0,6595 | 0,6818 |
| LDA | 0,7412 | 0,017 | 0,6895 | 0,7044 |
| SVM linear | 0,7401 | 0,020 | 0,7107 | 0,7127 |
| Naive Bayes | 0,7062 | 0,021 | 0,6968 | 0,6808 |
| k-NN (k = 15) | 0,6375 | 0,017 | 0,5584 | 0,5923 |
| Baseline (classe majoritária) | 0,3579 | 0,001 | 0,1667 | 0,0879 |

**Leitura da acurácia balanceada.** A acurácia simples é otimista quando as classes são desbalanceadas: o LumA sozinho representa 35,8% da coorte, e um modelo que sempre respondesse "LumA" já acertaria 35,8%. A acurácia balanceada é a média do acerto **dentro de cada classe**, e revela algo que a acurácia simples esconde: o Random Forest e o Extra Trees caem muito mais (0,678 e 0,660) do que a regressão logística (0,718), apesar de terem acurácia simples maior. Métodos baseados em árvores, sem reponderação de classes, favorecem as classes grandes.

### 2.2 O ranking é real? Os testes dizem que só em parte

**t pareado corrigido de Nadeau-Bengio**, tomando o Gradient Boosting como referência:

| Comparação | Diferença de acurácia | t | p corrigido | Conclusão |
|---|---|---|---|---|
| vs. SVM RBF | +0,0112 | 1,19 | **0,254** | **empate estatístico** |
| vs. Random Forest | +0,0209 | 1,79 | **0,095** | **empate estatístico** |
| vs. Rede neural (MLP) | +0,0281 | 3,09 | 0,008 | GB superior |
| vs. Regressão logística | +0,0290 | 2,69 | 0,018 | GB superior |
| vs. Extra Trees | +0,0309 | 2,69 | 0,018 | GB superior |
| vs. LDA | +0,0460 | 3,80 | 0,002 | GB superior |
| vs. SVM linear | +0,0471 | 3,25 | 0,006 | GB superior |
| vs. Naive Bayes | +0,0810 | 5,84 | 4 × 10⁻⁵ | GB superior |
| vs. k-NN | +0,1497 | 13,18 | < 10⁻⁸ | GB superior |
| vs. Baseline | +0,4293 | 39,84 | < 10⁻¹⁵ | GB superior |

**McNemar entre os dois primeiros colocados** (Gradient Boosting vs. SVM RBF, sobre as 1.897 predições fora da amostra): 130 pacientes que só o Gradient Boosting acerta contra 114 que só o SVM RBF acerta; χ² = 0,922; **p = 0,337**.

A conclusão é clara e vale ser dita sem rodeios: **o "vencedor" do ranking não é estatisticamente melhor que o segundo nem que o terceiro colocado**. A vantagem de 1,1 ponto percentual sobre o SVM RBF está dentro da flutuação amostral. Escolher o Gradient Boosting como modelo final é uma decisão legítima por outros critérios (velocidade, robustez a hiperparâmetros), mas não por superioridade demonstrada.

Por outro lado, três conclusões **são** estatisticamente sólidas: (i) todos os métodos superam massivamente o baseline; (ii) o k-NN é claramente inadequado para este problema — em 489 dimensões, a noção de "vizinho próximo" perde sentido, um efeito conhecido como maldição da dimensionalidade; (iii) modelos lineares (logística, LDA, SVM linear) ficam 3–5 pontos atrás dos não-lineares, indicando que existe estrutura de interação entre genes que os modelos lineares não capturam.

### 2.3 Onde exatamente o modelo erra

Matriz de confusão do Gradient Boosting (linhas = subtipo real, colunas = predito):

| real \ predito | LumA | LumB | Her2 | Basal | claudin-low | Normal |
|---|---|---|---|---|---|---|
| **LumA** | **592** | 57 | 17 | 1 | 1 | 11 |
| **LumB** | 67 | **377** | 12 | 1 | 4 | 0 |
| **Her2** | 18 | 32 | **159** | 4 | 4 | 3 |
| **Basal** | 5 | 3 | 16 | **156** | 12 | 7 |
| **claudin-low** | 13 | 3 | 3 | 20 | **156** | 3 |
| **Normal** | **67** | 1 | 7 | 5 | 11 | **49** |

Métricas por classe:

| Subtipo | Precisão | Recall | F1 | n |
|---|---|---|---|---|
| LumA | 0,777 | **0,872** | 0,822 | 679 |
| LumB | 0,797 | 0,818 | 0,807 | 461 |
| Basal | 0,834 | 0,784 | 0,808 | 199 |
| claudin-low | 0,830 | 0,788 | 0,808 | 198 |
| Her2 | 0,743 | 0,723 | 0,733 | 220 |
| **Normal** | 0,671 | **0,350** | **0,460** | 140 |

**Prova de que o erro não é aleatório:** qui-quadrado entre subtipo real e acerto/erro: **χ² = 195,3; p = 2,9 × 10⁻⁴⁰**. A taxa de acerto varia de 87,2% (LumA) a 35,0% (Normal).

**O erro tem uma geometria.** Todos os erros relevantes ocorrem entre pares **biologicamente adjacentes**: LumA ↔ LumB (124 erros nas duas direções), Her2 → LumB/LumA (50), Basal → Her2 (16), claudin-low → Basal (20), Normal → LumA (67). O que praticamente **não** acontece é igualmente informativo: apenas 5 casos Basal foram chamados de LumA e apenas 1 caso LumA foi chamado de Basal. **O modelo nunca confunde os extremos do espectro** — confunde vizinhos, exatamente onde a própria fronteira biológica é difusa.

### 2.4 Prova de que a maior falha é do rótulo, não do algoritmo

O subtipo "Normal-like" concentra o fracasso (recall 0,350; 67 dos 140 casos vão para LumA). Duas evidências independentes indicam que o problema está na natureza desse rótulo:

**Evidência 1 — a celularidade tumoral prediz o erro.** Testando se características da amostra se associam ao acerto (com correção FDR):

| Variável | Teste | p | FDR | Associado ao erro? |
|---|---|---|---|---|
| **Celularidade** | qui-quadrado | 0,0020 | **0,0111** | **sim** |
| **Integrative cluster** | qui-quadrado | 3 × 10⁻⁵ | **0,0003** | **sim** |
| Status de ER | qui-quadrado | 0,0176 | 0,064 | limítrofe |
| Tamanho tumoral | Mann-Whitney | 0,050 | 0,137 | não |
| Grau histológico | qui-quadrado | 0,394 | 0,620 | não |
| Status HER2 | qui-quadrado | 1,000 | 1,000 | não |
| Idade ao diagnóstico | Mann-Whitney | 0,512 | 0,701 | não |
| Contagem de mutações | Mann-Whitney | 0,573 | 0,701 | não |
| Índice de Nottingham | Mann-Whitney | 0,751 | 0,826 | não |

Taxa de acerto por celularidade: **alta 82,1%** (n = 936) → **moderada 76,6%** (n = 709) → **baixa 72,7%** (n = 198).

O gradiente é monótono e a interpretação é direta: quanto menor a fração de células tumorais na amostra, maior a contaminação por tecido mamário normal, e mais o perfil de expressão se afasta do tumor que se pretende classificar. **O modelo não está errando o cálculo — está recebendo, em parte, o transcriptoma do tecido errado.**

**Evidência 2 — nada do que é clinicamente relevante prediz o erro.** Grau, HER2, idade, NPI e carga mutacional têm FDR > 0,13. Se o modelo estivesse falhando por incapacidade de aprender, esperaríamos erro concentrado em tumores clinicamente difíceis. Não é o caso: **o erro segue a qualidade e a biologia da amostra, não a complexidade clínica do caso.**

### 2.5 O modelo sabe quando não sabe

Uma propriedade decisiva para uso prático: a confiança do modelo é informativa?

Confiança (probabilidade máxima) mediana nos acertos: **0,581**; nos erros: **0,442** — Mann-Whitney **p = 4,2 × 10⁻⁶⁰**.

| Faixa de confiança | n | Acerto observado | Confiança média |
|---|---|---|---|
| ≤ 0,4 | 324 | **47,5%** | 0,341 |
| 0,4 – 0,5 | 449 | 68,6% | 0,452 |
| 0,5 – 0,6 | 414 | 79,0% | 0,548 |
| 0,6 – 0,7 | 350 | 89,7% | 0,648 |
| 0,7 – 0,8 | 239 | 94,6% | 0,747 |
| 0,8 – 0,9 | 106 | 98,1% | 0,846 |
| > 0,9 | 15 | **100,0%** | 0,925 |

A relação é monótona e a calibração é conservadora (o acerto real supera a confiança declarada em todas as faixas). Isso habilita uma estratégia operacional concreta: **classificar automaticamente os casos com confiança acima de 0,7 (acurácia ≥ 94,6%, 360 pacientes) e encaminhar à revisão os 324 casos abaixo de 0,4**, onde o modelo acerta pouco mais que uma moeda.

---

## 3. Tarefa 2 — Predição de óbito em até 10 anos

Coorte de 1.560 pacientes com desfecho de 10 anos determinado (702 óbitos, 45,0%), evitando censura informativa. Doze algoritmos × três conjuntos de variáveis = 36 combinações.

### 3.1 Ranking (AUC média em 15 partições)

**Melhores de cada conjunto:**

| Conjunto | Melhor modelo | AUC | Brier |
|---|---|---|---|
| **Clínico (6 variáveis)** | **Rede neural (MLP)** | **0,7356** | 0,207 |
| **Combinado** | Gradient Boosting | 0,7329 | 0,234 |
| **Genes (489)** | Random Forest | **0,6561** | 0,232 |

Ranking completo dos dez primeiros:

| Posição | Conjunto \| Modelo | AUC |
|---|---|---|
| 1 | clínico \| Rede neural (MLP) | 0,7356 |
| 2 | combinado \| Gradient Boosting | 0,7329 |
| 3 | clínico \| SVM RBF | 0,7271 |
| 4 | clínico \| Random Forest | 0,7258 |
| 5 | combinado \| Logística LASSO | 0,7226 |
| 6 | clínico \| Extra Trees | 0,7210 |
| 7 | clínico \| Logística L2 | 0,7165 |
| 8 | combinado \| Random Forest | 0,7155 |
| 9 | clínico \| LDA | 0,7154 |
| 10 | clínico \| SVM linear | 0,7146 |

**Oito das dez primeiras posições usam apenas as 6 variáveis clínicas.** O melhor modelo genômico (0,656) fica na 21ª posição.

### 3.2 A prova estatística: os genes perdem em todos os algoritmos

Comparação pareada genes × clínico **dentro de cada algoritmo** (t de Nadeau-Bengio, correção FDR) — isso isola o efeito do conjunto de variáveis, eliminando o algoritmo como explicação:

| Algoritmo | AUC genes | AUC clínico | Diferença | FDR |
|---|---|---|---|---|
| Rede neural (MLP) | 0,614 | 0,736 | −0,122 | 0,0003 |
| k-NN | 0,602 | 0,714 | −0,112 | 0,0003 |
| Logística L2 | 0,604 | 0,717 | −0,112 | 0,0003 |
| SVM linear | 0,606 | 0,715 | −0,109 | 0,0005 |
| LDA | 0,607 | 0,715 | −0,109 | 0,0005 |
| SVM RBF | 0,655 | 0,727 | −0,072 | 0,0021 |
| Extra Trees | 0,651 | 0,721 | −0,070 | 0,0005 |
| Random Forest | 0,656 | 0,726 | −0,070 | 0,0016 |
| Naive Bayes | 0,630 | 0,693 | −0,064 | 0,0011 |
| Gradient Boosting | 0,638 | 0,700 | −0,061 | 0,0025 |
| Logística LASSO | 0,654 | 0,714 | −0,060 | 0,0006 |

**Onze de onze algoritmos, todos com FDR < 0,003, na mesma direção.** Não existe interpretação alternativa: a desvantagem dos 489 genes frente a 6 variáveis clínicas não é artefato de escolha de modelo, de linearidade, de regularização ou de hiperparâmetro. É uma propriedade dos dados.

**Teste de DeLong** (comparação formal de curvas ROC nas predições fora da amostra), tendo como referência o Gradient Boosting clínico (AUC 0,702):

| Modelo | AUC | Diferença | z | p | FDR |
|---|---|---|---|---|---|
| clínico \| SVM RBF | 0,731 | +0,029 | 2,72 | 0,007 | 0,013 |
| clínico \| Random Forest | 0,730 | +0,028 | 3,81 | 0,0001 | 0,0008 |
| combinado \| Gradient Boosting | 0,734 | +0,032 | 2,59 | 0,010 | 0,014 |
| genes \| Random Forest | 0,656 | −0,046 | −2,70 | 0,007 | 0,013 |
| genes \| SVM RBF | 0,653 | −0,049 | −2,90 | 0,004 | 0,010 |
| genes \| Gradient Boosting | 0,643 | −0,059 | −3,43 | 0,0006 | 0,002 |
| genes \| Logística L2 | 0,604 | −0,099 | −5,45 | < 10⁻⁷ | < 10⁻⁶ |

**Observação central:** acrescentar 489 genes ao modelo clínico **não melhora significativamente** o desempenho (combinado vs. clínico: +0,004, p = 0,82 no teste de Nadeau-Bengio) e, em vários algoritmos, **piora** — a logística combinada cai de 0,717 para 0,658. Isso é dilução de sinal: 489 preditores adicionais, quase todos sem informação prognóstica independente, aumentam a variância do modelo sem acrescentar sinal.

### 3.3 Onde o modelo acerta e onde falha, com intervalos de confiança

AUC por estrato do melhor modelo combinado, com IC 95% por bootstrap (1.500 reamostragens):

| Estrato | n | Eventos | AUC | IC 95% | Melhor que acaso? |
|---|---|---|---|---|---|
| **Coorte inteira** | 1.560 | 702 | 0,734 | 0,710 – 0,759 | sim |
| **LumA** | 544 | 183 | **0,772** | 0,728 – 0,815 | sim |
| claudin-low | 150 | 59 | 0,762 | 0,675 – 0,838 | sim |
| Normal | 111 | 52 | 0,759 | 0,656 – 0,849 | sim |
| Her2 | 184 | 109 | 0,712 | 0,636 – 0,785 | sim |
| LumB | 394 | 204 | 0,671 | 0,620 – 0,722 | sim |
| **Basal** | 177 | 95 | **0,572** | **0,488 – 0,655** | **NÃO** |
| Estágio I | 390 | 113 | 0,679 | 0,622 – 0,737 | sim |
| Estágio II | 668 | 324 | 0,711 | 0,669 – 0,748 | sim |
| Estágio III–IV | 109 | 78 | 0,650 | 0,533 – 0,759 | sim |

**O subtipo Basal é o único estrato onde o intervalo de confiança inclui 0,5.** Essa é a demonstração formal — e não apenas a observação de um número baixo — de que o modelo **não funciona** em tumores basais. O IC 95% de 0,488 a 0,655 significa que, nesse grupo, não se pode rejeitar a hipótese de que a predição equivale ao acaso.

Taxa de acerto por subtipo (χ² = 27,48; **p = 4,6 × 10⁻⁵**): LumA 73,3%, claudin-low 68,7%, Her2 67,4%, LumB 61,4%, Normal 61,3%, **Basal 55,4%**. Já por estágio, a diferença **não** é significativa (χ² = 4,55; p = 0,103) — o modelo é razoavelmente estável entre estágios, e instável entre subtipos.

### 3.4 Que tipo de erro o modelo comete

Usando a mediana do risco previsto como ponto de corte (acurácia global 66,3%):

| Resultado | n |
|---|---|
| Acerto: baixo risco / paciente viva | 556 |
| Acerto: alto risco / óbito | 478 |
| **Falso positivo** (alarme sem óbito) | **302** |
| **Falso negativo** (óbito não previsto) | **224** |

**Teste binomial:** o modelo erra de forma assimétrica, com mais falsos positivos que falsos negativos (302 vs. 224; **p = 0,0008**). Em contexto de estratificação de risco, esse é o viés menos danoso — o modelo peca por cautela.

**O perfil dos erros, com teste formal.** Comparando os óbitos **não previstos** com os óbitos **previstos** (Mann-Whitney, todos com FDR < 0,001):

| Variável | Óbito não previsto | Óbito previsto | FDR |
|---|---|---|---|
| Idade ao diagnóstico | 61,3 anos | 69,1 anos | < 0,001 |
| Tamanho tumoral | 20 mm | 30 mm | < 0,001 |
| Índice de Nottingham | 4,04 | 5,04 | < 0,001 |
| Linfonodos positivos | 0 | 2 | < 0,001 |
| Tempo até o óbito | 70,1 meses | 46,0 meses | < 0,001 |

Aqui está o retrato mais nítido da falha: **os óbitos que o modelo não antecipa são de pacientes mais jovens, com tumores menores, sem linfonodos comprometidos e índice prognóstico favorável — que morrem mais tarde dentro da janela de 10 anos.** Em outras palavras, o modelo captura bem a morte que "parece" provável pelos critérios clínicos clássicos e falha justamente nos casos que mais precisariam de um marcador molecular: pacientes com aparência clínica de baixo risco e evolução desfavorável. Esse é o nicho não atendido, e ele coincide exatamente com o resultado de que os genes não acrescentam sinal — se acrescentassem, seria aqui.

A composição do erro também depende do subtipo (χ² = 188,6; **p = 4,3 × 10⁻³²**). O grupo Normal-like tem a maior proporção de falsos negativos (25,2%), e Basal e LumB as maiores proporções de falsos positivos (28,8% e 24,9%).

### 3.5 Calibração: o risco previsto é confiável como número?

| Faixa de risco previsto | n | Risco médio previsto | Óbito observado | Diferença |
|---|---|---|---|---|
| ≤ 0,2 | 609 | 0,068 | **0,251** | −0,183 |
| 0,2 – 0,3 | 120 | 0,253 | 0,367 | −0,114 |
| 0,3 – 0,4 | 87 | 0,344 | 0,506 | −0,162 |
| 0,4 – 0,5 | 110 | 0,452 | 0,436 | +0,016 |
| 0,5 – 0,6 | 105 | 0,545 | 0,486 | +0,060 |
| 0,6 – 0,7 | 91 | 0,649 | 0,571 | +0,077 |
| > 0,7 | 438 | 0,882 | **0,708** | +0,174 |

**Hosmer-Lemeshow: χ² = 779,4; gl = 8; p = 5,7 × 10⁻¹⁶³.** A calibração é **claramente ruim**, e de forma sistemática: o modelo é excessivamente confiante nos dois extremos. Entre as pacientes a quem atribui risco de 6,8%, morrem 25,1%; entre aquelas a quem atribui 88,2%, morrem 70,8%.

Isso tem consequência prática direta e merece destaque: **a ordenação de risco funciona, mas o número absoluto não pode ser comunicado como probabilidade.** Dizer a uma paciente "seu risco é de 7%" quando a taxa real no grupo é de 25% seria um erro clínico grave. Antes de qualquer uso, o modelo exigiria recalibração (Platt scaling ou regressão isotônica) — e é justamente por isso que uma métrica de discriminação como o AUC, isoladamente, é insuficiente para avaliar um modelo de risco.

### 3.6 A ordenação, essa sim, funciona

Dividindo a coorte em tercis pelo risco previsto e observando a sobrevida **real** (não a dicotomizada em 10 anos):

| Tercil | n | Óbitos | Sobrevida global mediana |
|---|---|---|---|
| Baixo | 520 | 246 | **227,9 meses** |
| Médio | 520 | 358 | 144,7 meses |
| Alto | 520 | 433 | **74,5 meses** |

**Log-rank: p = 5,2 × 10⁻⁶¹.** Uma diferença de **153 meses** — mais de doze anos — na sobrevida mediana entre os tercis extremos. Ou seja: apesar da calibração ruim e do AUC modesto, o modelo **ordena** as pacientes de forma extremamente significativa. Discriminação e calibração são propriedades distintas, e este modelo tem a primeira sem a segunda.

---

## 4. Síntese — o que ficou demonstrado

**Sobre a comparação de algoritmos.** Onze algoritmos com pressupostos distintos convergem para o mesmo teto: ~0,79 de acurácia na classificação de subtipos e ~0,73 de AUC na predição de sobrevida. As diferenças entre os três primeiros colocados **não são estatisticamente significativas** (Nadeau-Bengio p ≥ 0,095; McNemar p = 0,337). A escolha do algoritmo não é o gargalo — a informação disponível é.

**Sobre onde o modelo acerta.** Classifica com alta confiabilidade os subtipos com identidade molecular definida (LumA recall 0,872; Basal e claudin-low F1 = 0,808) e ordena o risco de sobrevida de forma extremamente robusta (log-rank p = 5 × 10⁻⁶¹; 153 meses de diferença entre tercis extremos). Sabe reconhecer sua própria incerteza (confiança correlacionada com acerto, p = 4 × 10⁻⁶⁰; 100% de acerto acima de 0,9 de confiança).

**Sobre onde falha, e por quê.** (i) O grupo Normal-like, com recall de 0,350, cuja falha se explica pela celularidade tumoral (FDR = 0,011) e não por características clínicas (todas com FDR > 0,13) — falha do rótulo, não do algoritmo. (ii) O subtipo Basal na predição de sobrevida, único estrato cujo IC 95% inclui o acaso (0,488–0,655). (iii) A calibração absoluta (Hosmer-Lemeshow p = 6 × 10⁻¹⁶³), que impede comunicar o risco como probabilidade. (iv) Os óbitos de pacientes clinicamente favoráveis — mais jovens, tumores menores, sem linfonodos comprometidos (todos FDR < 0,001) —, que são exatamente os casos onde um marcador molecular seria mais valioso.

**Sobre a expressão gênica como preditor de desfecho.** Onze de onze algoritmos mostram os 489 genes inferiores às 6 variáveis clínicas, todos com FDR < 0,003 e na mesma direção; adicionar os genes ao modelo clínico não traz ganho (p = 0,82) e frequentemente prejudica. Este é o achado mais robusto de toda a análise, precisamente por ser independente de qualquer escolha metodológica.

---

## 5. Arquivos gerados

| Arquivo | Conteúdo |
|---|---|
| M01, M02 | Desempenho por partição e resumo dos 11 modelos (subtipo) |
| M03 | Teste t corrigido de Nadeau-Bengio (subtipo) |
| M04, M05 | Matriz de confusão e métricas por classe |
| M06 | Predição individual por paciente (subtipo) |
| M07 | Testes de associação entre características da amostra e erro |
| M08 | Tabela de calibração da confiança |
| M10, M11 | Desempenho por partição e resumo das 36 combinações (sobrevida) |
| M12 | Teste t corrigido pareado (sobrevida) |
| M13 | Comparação genes × clínico dentro de cada algoritmo |
| M14 | Predições fora da amostra dos 12 melhores modelos |
| M15 | Teste de DeLong |
| M16 | Predição individual por paciente (sobrevida), com tipo de erro e tercil |
| M17 | AUC por subtipo e estágio com IC 95% bootstrap |
| M18, M21 | Perfil clínico dos erros e comparação óbito previsto × não previsto |
| M19 | Medianas por tipo de resultado |
| M20 | Tabela de calibração e Hosmer-Lemeshow |

Scripts: `ml_tarefa1_subtipo.py`, `ml_tarefa2_sobrevida.py`, `ml_tarefa2_final.py`. Semente 42 em todas as etapas; resultados integralmente reproduzíveis.
