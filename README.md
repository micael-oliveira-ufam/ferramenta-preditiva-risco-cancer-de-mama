# Pipeline reprodutível — METABRIC: expressão gênica, sobrevida e modelos preditivos

Pacote completo: todo o código de treinamento, as análises de expressão gênica, o cálculo do AUC
de cada método, as figuras de validação, as tabelas de genes associados à sobrevida e a plataforma
de inferência que roda no navegador.

**Coorte:** METABRIC (Curtis et al., *Nature* 2012; Pereira et al., *Nat Commun* 2016)
**Semente fixa 42 em todas as etapas aleatórias** — os números deste pacote são reproduzidos integralmente.

---

## 1. Estrutura do pacote

| Pasta | Conteúdo |
|---|---|
| `01_pipeline_R/` | `metabric_pipeline.R` — pipeline principal em R, 12 módulos, autossuficiente |
| `02_analises_complementares/` | Estágio tumoral, interação gene × tratamento, Random Forest complementar |
| `03_benchmark_ml/` | Benchmark de 11 algoritmos com testes estatísticos formais |
| `04_treino_e_exportacao/` | Treino final dos 16 modelos e exportação para o navegador |
| `05_figuras/` | 8 figuras de validação + o script que as gera |
| `06_tabelas/` | Genes de sobrevida e métricas dos 16 modelos |
| `07_plataforma/` | Aplicação web de inferência (roda offline, abrindo `index.html`) |
| `08_resultados_brutos/` | Todas as tabelas intermediárias de cada etapa |

---

## 2. Como reproduzir do zero

### Requisitos

```bash
# R (pipeline principal) — o script instala os pacotes que faltarem
R --version                       # ≥ 4.3

# Python (análises complementares, ML, exportação)
pip install pandas numpy scipy scikit-learn statsmodels lifelines xgboost matplotlib
```

### Ordem de execução

```bash
# Etapa 1 — pipeline principal (~7 min). Baixa os dados e confere o MD5.
cd 01_pipeline_R && Rscript metabric_pipeline.R

# Etapa 2 — estágio, farmacogenômica e RF complementar (~40 min)
cd ../02_analises_complementares
python3 analise_complementar.py     # Cox por estágio nos 489 genes
python3 parte_b_farmaco.py          # interação gene × tratamento (1.467 modelos)
python3 parte_c_rf.py               # RF de subtipo e de mortalidade em 10 anos

# Etapa 3 — benchmark de algoritmos com testes estatísticos (~30 min)
cd ../03_benchmark_ml
python3 ml_tarefa1_subtipo.py       # 11 algoritmos, classificação de subtipo
python3 ml_tarefa2_sobrevida.py     # 12 algoritmos × 3 conjuntos, sobrevida
python3 ml_tarefa2_final.py         # análise de erro, calibração, log-rank

# Etapa 4 — treino final e exportação para o navegador (~5 min)
cd ../04_treino_e_exportacao
python3 treina_exporta.py           # 16 modelos → JSON executável em JS

# Etapa 5 — figuras e tabelas finais (~1 min)
cd ../05_figuras && python3 gera_figuras.py
```

O arquivo de dados (`METABRIC_RNA_Mutation.csv`, MD5 `c619471beb87af0f9fd4e4a40058654d`) é baixado
automaticamente pela etapa 1 e reaproveitado pelas demais; espera-se encontrá-lo em `dados/`.

**Cache entre execuções.** A etapa 4 grava um arquivo por modelo em `saidas_multi/parcial/`. Se a
execução for interrompida, basta rodar de novo: os modelos já concluídos são pulados.

**Atenção ao paralelismo.** Em máquinas com poucos núcleos, `n_jobs=-1` causa contenção de threads
e chega a multiplicar o tempo por 40. Os scripts já vêm com `n_jobs=1`; rode com
`OMP_NUM_THREADS=1` se notar lentidão anômala.

---

## 3. Como o AUC de cada método é obtido

Este é o ponto onde a maioria dos trabalhos superestima o desempenho. O procedimento usado aqui,
implementado em `04_treino_e_exportacao/treina_exporta.py`:

1. **Definição da coorte.** Só entram pacientes cujo desfecho de 10 anos é determinado —
   quem morreu antes de 120 meses ou quem foi seguido por pelo menos 120 meses (n = 1.560,
   702 óbitos). Isso evita censura informativa.

2. **Partição externa.** `StratifiedKFold(5, shuffle=True, random_state=42)`.

3. **Busca de hiperparâmetros DENTRO de cada dobra.** `RandomizedSearchCV` com validação
   interna de 3 dobras (`random_state=43`). Escolher hiperparâmetros com os dados de teste à
   vista infla o AUC; aqui o teste externo nunca participa da escolha.

4. **Seleção de variáveis dentro da dobra.** Quando há seleção L1, ela é um passo do `Pipeline`,
   não uma etapa prévia. *Este ponto custou uma correção no meio do projeto:* a primeira versão
   selecionava os genes com toda a base antes de validar, o que elevou artificialmente o AUC do
   modelo genômico de **0,639 para 0,727** — e, pior, fazia o subtipo Basal parecer confiável.

5. **Predições fora da amostra.** `cross_val_predict` devolve, para cada paciente, a probabilidade
   estimada por um modelo que nunca o viu no treino.

6. **AUC e IC 95%.** Calculados sobre essas predições com o método de **DeLong**, que estima a
   variância da AUC a partir dos midranks — e não por bootstrap ingênuo.

7. **Calibração.** Regressão isotônica ajustada dentro da validação, avaliada por escore de
   Brier e teste de Hosmer-Lemeshow.

8. **Verificação da reimplementação.** Cada modelo exportado é reexecutado em Python puro
   (percorrendo as árvores nó a nó, como o JavaScript faz) e comparado ao scikit-learn.
   **Erro máximo: 3×10⁻⁶.**

---

## 4. Resultados: desempenho dos 16 modelos

Tabela completa em `06_tabelas/G04_metricas_16_modelos.csv`; figuras `F01` e `F02`.

| Conjunto | Algoritmo | Var. | AUC | IC 95% | Brier cal. | H-L p |
|---|---|---|---|---|---|---|
| **Combinado** | **XGBoost** | 500 | **0,7488** | 0,725–0,773 | 0,198 | 0,999 |
| Clínico | Rede neural (MLP) | 11 | 0,7379 | 0,713–0,762 | 0,202 | 0,95 |
| Clínico | XGBoost | 11 | 0,7352 | 0,711–0,760 | 0,203 | 0,43 |
| Clínico | Random Forest | 11 | 0,7334 | 0,709–0,758 | 0,203 | 0,99 |
| Clínico | Gradient Boosting | 11 | 0,7328 | 0,708–0,758 | 0,204 | 0,52 |
| Combinado | Extra Trees | 500 | 0,7253 | 0,701–0,750 | 0,206 | 0,998 |
| Clínico | Regressão logística | 11 | 0,7190 | 0,694–0,744 | 0,208 | 0,65 |
| Clínico | Extra Trees | 11 | 0,7182 | 0,693–0,743 | 0,208 | 0,94 |
| Combinado | Random Forest | 500 | 0,7047 | 0,679–0,730 | 0,213 | 0,84 |
| Combinado | Regressão logística | 500 | 0,6914 | 0,665–0,717 | 0,217 | 0,26 |
| Combinado | Rede neural (MLP) | 500 | 0,6702 | 0,644–0,697 | 0,224 | 0,97 |
| Genes | XGBoost | 489 | 0,6632 | 0,636–0,690 | 0,223 | 0,99 |
| Genes | Random Forest | 489 | 0,6511 | 0,624–0,678 | 0,226 | 0,88 |
| Genes | Extra Trees | 489 | 0,6488 | 0,622–0,676 | 0,227 | 0,996 |
| Genes | Regressão logística | 489 | 0,6314 | 0,604–0,659 | 0,231 | 0,83 |
| Genes | Rede neural (MLP) | 489 | 0,6194 | 0,592–0,647 | 0,234 | 0,94 |

**Leitura honesta do ranking.** Os intervalos dos treze primeiros se sobrepõem quase inteiramente.
O XGBoost combinado é o único modelo que supera o melhor puramente clínico, mas seu IC
(0,725–0,773) cobre o do MLP clínico (0,713–0,762): a diferença de 1,1 ponto **não é conclusiva**.
O teto de ~0,75 não cede a troca de algoritmo — é limite da informação nos dados.

### Acurácia do classificador de subtipo (figura `F08`)

Gradient Boosting, 6 classes, validação cruzada: **acurácia 78,7%**, F1-macro 0,743. Recall por
classe: LumA 0,872 · LumB 0,818 · Basal 0,784 · claudin-low 0,788 · Her2 0,723 · **Normal 0,350**.
O erro concentra-se entre pares biologicamente adjacentes e associa-se à **celularidade tumoral**
(FDR = 0,011), não a características clínicas (todas com FDR > 0,13) — é falha do rótulo, não do
algoritmo. Testes em `08_resultados_brutos/benchmark_ml/M03`, `M05`, `M07`.

---

## 5. Genes associados à sobrevida

Cox univariado nos 489 genes, HR por **+1 desvio-padrão de expressão** (z-score), correção
Benjamini-Hochberg. **125 genes significativos (FDR < 0,05): 61 protetores e 64 de risco.**
Tabelas em `06_tabelas/G01`–`G03`; figura `F06`.

### Genes associados a MAIOR sobrevida (HR < 1)

| Gene | HR (OS) | IC 95% | HR (DSS) |
|---|---|---|---|
| **IGF1** | 0,804 | 0,750–0,862 | 0,820 |
| **STAT5A** | 0,807 | 0,759–0,858 | 0,753 |
| **SPRY2** | 0,810 | 0,759–0,864 | 0,864 |
| **CDKN2C** | 0,815 | 0,760–0,874 | 0,874 |
| **FLT3** | 0,820 | 0,767–0,876 | 0,792 |
| **ABCB1** | 0,822 | 0,770–0,877 | 0,875 |
| LAMA2 | 0,836 | 0,788–0,887 | 0,855 |
| STAT5B | 0,838 | 0,790–0,890 | 0,744 |
| PDGFRA | 0,844 | 0,796–0,895 | 0,945 |
| CCND2 | 0,846 | 0,796–0,899 | 0,898 |
| BCL2 | 0,849 | 0,801–0,900 | 0,723 |
| KIT | 0,849 | 0,796–0,905 | 0,868 |

### Genes associados a MENOR sobrevida (HR > 1)

| Gene | HR (OS) | IC 95% | HR (DSS) |
|---|---|---|---|
| **GSK3B** | 1,225 | 1,159–1,295 | 1,326 |
| **AURKA** | 1,191 | 1,123–1,262 | **1,458** |
| **VEGFA** | 1,165 | 1,100–1,233 | 1,209 |
| **FANCD2** | 1,163 | 1,096–1,233 | 1,323 |
| KRAS | 1,153 | 1,086–1,224 | 1,199 |
| E2F7 | 1,148 | 1,084–1,215 | 1,276 |
| TUBB4B | 1,144 | 1,078–1,214 | 1,227 |
| SLC19A1 | 1,142 | 1,078–1,209 | 1,239 |
| MMP11 | 1,134 | 1,067–1,206 | 1,162 |
| CTCF | 1,132 | 1,068–1,200 | 1,166 |
| RPS6KB2 | 1,131 | 1,066–1,200 | — |
| MAML1 | 1,128 | 1,062–1,197 | — |

Os HR são **maiores na sobrevida específica do câncer (DSS)** do que na global, como esperado ao
remover o ruído da mortalidade competitiva — *AURKA* passa de 1,191 para 1,458.

### Convergência entre três métodos independentes

O núcleo **GSK3B, AURKA, STAT5A, FLT3, BCL2** aparece simultaneamente no Cox univariado, na
seleção estável do LASSO (GSK3B em 72% das 25 repetições, STAT5A em 60%) e no topo da importância
do Random Forest para o desfecho (figura `F07`). Três famílias metodológicas com pressupostos
distintos chegam ao mesmo conjunto.

### O sinal é dependente do contexto

Testando dentro de cada estrato (figura `F07`, tabela `C01`):

- **Subtipo LumA:** 12 de 12 genes significativos (*PDGFRA* HR 0,695, p = 5,4×10⁻¹¹)
- **Subtipos Basal e Her2:** **nenhum** gene sobrevive à correção
- **Estágio II:** 71 genes significativos — *ACVR1B* 1,279, *STAT5A* 0,787, *IGF1* 0,762
- **Estágio I:** apenas *NRIP1* (0,729) · **Estágio III–IV:** apenas *DNAH11* (1,523)

Um painel prognóstico derivado da coorte inteira é, na prática, um painel para doença luminal.

---

## 6. Figuras de validação

| Figura | O que mostra |
|---|---|
| `F01_curvas_roc.png` | Curvas ROC dos 6 algoritmos, em cada um dos 3 conjuntos |
| `F02_auc_ic95.png` | AUC com IC 95% de DeLong para os 16 modelos, ordenados |
| `F03_calibracao.png` | Risco previsto × observado após correção isotônica (H-L p = 0,66) |
| `F04_confiabilidade_estratos.png` | AUC por subtipo e por estágio, com a linha do acaso marcada |
| `F05_tercis_risco.png` | Desfecho e sobrevida real por tercil (log-rank p = 1,8×10⁻³⁶) |
| `F06_genes_forest.png` | Forest plot dos 24 genes mais fortes, protetores e de risco |
| `F07_genes_ml_e_estagio.png` | Importância no Random Forest e genes do estágio II |
| `F08_matriz_confusao_subtipo.png` | Matriz de confusão da classificação dos 6 subtipos |

---

## 7. Pipeline de inferência

O modelo treinado é exportado como coeficientes e arrays de nós de árvore em
`07_plataforma/modelos.json`, e a predição é reimplementada em JavaScript puro. Para usar em
outra linguagem, o algoritmo é o seguinte:

```
1. Monte o vetor x na ordem exata de modelo["variaveis"]
   - nomes iniciados por "sub_"        → 1 se for o subtipo do paciente, 0 caso contrário
   - nomes de variáveis clínicas       → o valor informado (NUNCA zero quando ausente)
   - demais nomes                      → o z-score do gene (0 = média da coorte)

2. Conforme modelo["tipo"]:
   linear        z = intercepto + Σ ((x[i] − media[i]) / escala[i]) × coef[i]
                 p = 1 / (1 + e^(−z))
   floresta      p = média, sobre as árvores, do valor da folha alcançada
   boosting_sk   z = base + lr × Σ (valor da folha de cada árvore);  p = sigmoide(z)
   boosting_xgb  z = logit(base) + Σ (valor da folha de cada árvore); p = sigmoide(z)
   mlp           padronize, aplique W/b com ReLU camada a camada, depois a saída logística

   Percurso da árvore: comece no nó 0; enquanto l[i] ≠ −1,
   vá para l[i] se x[f[i]] < t[i], senão para r[i]. Devolva v[i].
   ATENÇÃO: a folha é identificada por l[i] = −1. No sklearn, f[i] vale −2 nas folhas,
   e usar esse campo como teste de parada gera laço infinito.

3. Calibre: interpole linearmente o p bruto na curva modelo["calibracao"] (x → y).
   O valor calibrado é a probabilidade final.
```

Para o **XGBoost**, `base` é a prevalência real (≈ 0,45), não 0,5 — ela é lida de
`booster.save_config()`. Usar 0,5 produz erro de 4,8×10⁻² na probabilidade.

**Abrir a plataforma:** basta abrir `07_plataforma/index.html` no navegador; funciona offline,
sem servidor. Para publicar: `wrangler pages deploy 07_plataforma`.

---

## 8. Limites que o pipeline não resolve

1. **Sem validação externa.** Toda estimativa é interna. Validação em TCGA-BRCA é o passo seguinte
   obrigatório antes de qualquer afirmação de utilidade.
2. **Coorte histórica**, majoritariamente anterior ao trastuzumabe — o prognóstico do subtipo Her2
   não corresponde à prática atual.
3. **Basal e estágio III–IV:** o IC 95% da AUC inclui 0,5. Nesses estratos a predição não se
   distingue do acaso, e a plataforma bloqueia a leitura do número.
4. **Discriminação moderada (≈ 0,75)** e ausência de randomização de tratamento.
5. **Os óbitos não previstos** são de pacientes mais jovens, com tumores menores e sem linfonodos
   comprometidos (todos FDR < 0,001) — exatamente o nicho onde um marcador molecular seria mais
   valioso.

Uso restrito a pesquisa e ensino. Não é dispositivo médico.
