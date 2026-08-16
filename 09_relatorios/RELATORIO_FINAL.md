# Expressão gênica, prognóstico e aprendizado de máquina no câncer de mama

## Relatório final consolidado — coorte METABRIC (n = 1.897)

**Instituição:** Faculdade de Ciências Farmacêuticas, Universidade Federal do Amazonas
**Coorte:** METABRIC — Curtis et al., *Nature* 2012; Pereira et al., *Nat Commun* 2016
**Dados:** `METABRIC_RNA_Mutation.csv`, MD5 `c619471beb87af0f9fd4e4a40058654d`
**Ambiente:** R 4.3.3 (22 pacotes) e Python 3.12 (scikit-learn 1.8.0, lifelines 0.30.3, XGBoost 3.4.1)
**Semente fixa 42** em todas as etapas aleatórias — todos os números deste relatório são reproduzíveis.

---

## Sumário

1. [Sumário executivo](#1-sumário-executivo)
2. [Como ler os números](#2-como-ler-os-números)
3. [Dados, coorte e integridade](#3-dados-coorte-e-integridade)
4. [Métodos](#4-métodos)
5. [Resultados — caracterização clínica e molecular](#5-resultados--caracterização-clínica-e-molecular)
6. [Resultados — genes associados à sobrevida](#6-resultados--genes-associados-à-sobrevida)
7. [Resultados — aprendizado de máquina](#7-resultados--aprendizado-de-máquina)
8. [Resultados — onde o modelo falha](#8-resultados--onde-o-modelo-falha)
9. [Camada farmacogenômica](#9-camada-farmacogenômica)
10. [Plataforma de inferência](#10-plataforma-de-inferência)
11. [Erros cometidos e corrigidos](#11-erros-cometidos-e-corrigidos)
12. [Limitações](#12-limitações)
13. [Próximos passos](#13-próximos-passos)
14. [Inventário do pacote](#14-inventário-do-pacote)
15. [Referências](#15-referências)

---

## 1. Sumário executivo

Este trabalho analisou 1.897 pacientes da coorte METABRIC combinando expressão de 489 genes,
status mutacional de 173 genes e 31 variáveis clínicas, com três objetivos: identificar genes
associados à sobrevida, comparar algoritmos de aprendizado de máquina para prever desfecho, e
construir uma plataforma de inferência com limites de confiabilidade explícitos.

**Os dez achados centrais:**

1. **125 genes associam-se à sobrevida global** (FDR < 0,05): 61 protetores e 64 de risco.
   Núcleo de pior prognóstico: *GSK3B*, *AURKA*, *VEGFA*, *FANCD2*. Núcleo de melhor prognóstico:
   *IGF1*, *STAT5A*, *SPRY2*, *CDKN2C*, *FLT3*, *ABCB1*.

2. **Três famílias metodológicas convergem** para o mesmo núcleo (*GSK3B*, *AURKA*, *STAT5A*,
   *FLT3*, *BCL2*): Cox univariado, seleção estável do LASSO em 25 repetições e importância do
   Random Forest — apesar de pressupostos radicalmente distintos.

3. **O sinal gênico é dependente do contexto.** No subtipo LumA, 12 de 12 genes testados são
   significativos; nos subtipos Basal e Her2, **nenhum**. No estágio II, 71 genes; nos estágios I
   e III–IV, apenas um em cada. Um painel derivado da coorte inteira é, na prática, um painel
   para doença luminal.

4. **568 associações gene–subtipo** com assinaturas biologicamente coerentes, incluindo um
   controle positivo interno perfeito (*ERBB2* superexpresso exatamente no subtipo Her2).

5. **A classificação dos subtipos por expressão atinge 78,7% de acurácia** (Gradient Boosting),
   com erro concentrado entre pares biologicamente adjacentes.

6. **A maior falha do classificador é do rótulo, não do algoritmo.** O grupo "Normal-like" tem
   recall de 0,350, e o erro se associa à celularidade tumoral (FDR = 0,011) — não a grau, HER2,
   idade ou carga mutacional (todos FDR > 0,13).

7. **Dezesseis modelos preditivos foram treinados** (6 algoritmos × 3 conjuntos de variáveis) com
   busca de hiperparâmetros aninhada. O melhor é XGBoost sobre variáveis clínicas + expressão
   (AUC 0,749), mas seu intervalo de confiança se sobrepõe ao do segundo colocado.

8. **A expressão gênica isolada perde para 6 variáveis clínicas** em onze de onze algoritmos
   testados, todos com FDR < 0,003 e na mesma direção — o achado mais robusto do trabalho,
   precisamente por ser independente de escolha metodológica.

9. **A ordenação de risco funciona; a probabilidade absoluta exigia correção.** Os tercis de
   risco separam sobrevida mediana de 217,6 a 85,7 meses (log-rank p = 1,8 × 10⁻³⁶), mas o
   modelo era grosseiramente descalibrado (Hosmer-Lemeshow p = 5,7 × 10⁻¹⁶³) até a aplicação de
   regressão isotônica ajustada dentro da validação (p = 0,66 depois).

10. **Existe sinal farmacogenômico para hormonioterapia:** 14 genes modificam significativamente
    o efeito prognóstico conforme a paciente recebeu ou não bloqueio endócrino, e todos resistem
    ao ajuste por status de ER e subtipo molecular.

---

## 2. Como ler os números

| Termo | Significado neste relatório |
|---|---|
| **z-score de expressão** | Expressão padronizada: 0 é a média da coorte, +1 é um desvio-padrão acima |
| **HR (hazard ratio)** | Risco relativo de óbito por **+1 desvio-padrão** de expressão. HR 1,20 = 20% mais risco; 0,80 = 20% menos |
| **IC 95%** | Faixa compatível com os dados. Se não cruza 1 (HR) ou 0,5 (AUC), o efeito é detectável |
| **FDR** | Valor-p corrigido para múltiplos testes (Benjamini-Hochberg) |
| **Δz** | Diferença de expressão média entre um subtipo e os demais, em desvios-padrão |
| **AUC** | Probabilidade de ordenar corretamente dois pacientes quanto ao risco. 0,5 = acaso |
| **C-index** | Análogo da AUC para dados de sobrevida censurados |
| **Escore de Brier** | Erro quadrático médio da probabilidade prevista. Menor é melhor |
| **Hosmer-Lemeshow** | Teste de calibração. **p alto = boa calibração** (ao contrário do usual) |
| **OS / DSS** | Sobrevida global (qualquer causa) / específica do câncer |
| **Erro out-of-bag** | Estimativa honesta de erro do Random Forest, com os casos que cada árvore não viu |

---

## 3. Dados, coorte e integridade

| Item | Valor |
|---|---|
| Pacientes no arquivo bruto | 1.904 |
| **Coorte analítica** | **1.897** |
| Colunas | 693 (31 clínicas + 489 expressão + 173 mutação) |
| Plataforma | Illumina HT-12 v3 (microarranjo); sequenciamento alvo para mutações |
| Seguimento mediano | 115,6 meses |
| Óbitos por qualquer causa | 1.098 |
| Óbitos atribuídos ao câncer | 619 |
| IDs duplicados / ausentes | 0 / 0 |

**Exclusões:** 1 sarcoma mamário, 6 casos sem classificação de subtipo (`NC`).

**Auditoria da codificação do desfecho — etapa crítica.** Este conjunto de dados é frequentemente
analisado com o evento invertido, o que produziria conclusões exatamente opostas. A tabulação
cruzada (`T02`) confirma que `overall_survival = 1` significa **paciente viva**; o evento de óbito
foi definido como `1 − overall_survival`.

**Coorte de modelagem preditiva.** Para os modelos de óbito em 10 anos, só entram pacientes cujo
desfecho é determinado — quem morreu antes de 120 meses ou foi seguido por pelo menos 120 meses:
**n = 1.560, 702 óbitos (45,0%)**. Isso evita censura informativa.

**Estágio tumoral:** disponível em 1.400 pacientes (73,8%). Agrupamento: estágios 0–1 como **I**
(n = 479), estágio 2 como **II** (n = 797), estágios 3–4 como **III–IV** (n = 124).

---

## 4. Métodos

### 4.1 Pipeline principal em R — 12 módulos

Parâmetros declarados *a priori* em um único objeto de configuração, evitando escolhas post-hoc:
FDR < 0,01 e |Δz| ≥ 0,50 para expressão diferencial; FDR < 0,05 para modelos de Cox; partição
70/30; 500 árvores; 10 dobras na validação cruzada; frequência mínima de 2% para mutações.

| Módulo | Função |
|---|---|
| 0–2 | Ambiente, ingestão, integridade, auditoria do desfecho, matriz 489 × 1.897 |
| 3 | Caracterização clínica por subtipo (qui-quadrado, Kruskal-Wallis, correção BH) |
| 4 | Expressão diferencial com **limma** — modelos lineares por gene, moderação bayesiana empírica |
| 5 | **Random Forest** de classificação dos 6 subtipos, erro OOB, importância |
| 6 | **Cox univariado** dos 489 genes em OS e DSS; teste de riscos proporcionais; Cox ajustado |
| 7 | Kaplan-Meier e **Cox estratificado dentro de cada subtipo** |
| 8 | **Assinatura LASSO-Cox** com validação em partição independente |
| 8b | **Validação repetida em 25 partições** e estabilidade de seleção |
| 9 | Mutações somáticas: heterogeneidade e associação com sobrevida |
| 10 | Integração das camadas e escore de prioridade por subtipo |

### 4.2 Análises complementares em Python

| Análise | Desenho |
|---|---|
| Cox por estágio | 489 modelos univariados dentro de cada estrato, correção BH separada |
| Interação gene × tratamento | 1.467 modelos com termo de interação, ajustados por idade, grau e tamanho |
| Random Forest de desfecho | Óbito em 10 anos, validação cruzada 5-fold estratificada |
| Benchmark de 11 algoritmos | 5-fold × 3 repetições, testes de Nadeau-Bengio, McNemar e DeLong |
| Treino final multi-algoritmo | 16 modelos com busca de hiperparâmetros aninhada |

### 4.3 Como o AUC de cada método é obtido

Este é o ponto onde a maioria dos trabalhos superestima o desempenho. O procedimento:

1. **Partição externa:** `StratifiedKFold(5, shuffle=True, random_state=42)`.
2. **Busca de hiperparâmetros DENTRO de cada dobra:** `RandomizedSearchCV` com 3 dobras internas.
   Escolher hiperparâmetros com o teste à vista infla a AUC.
3. **Seleção de variáveis dentro da dobra:** quando há seleção L1, ela é um passo do `Pipeline`,
   não uma etapa prévia.
4. **Predições fora da amostra:** `cross_val_predict` devolve, para cada paciente, a probabilidade
   estimada por um modelo que nunca o viu no treino.
5. **AUC e IC 95%:** método de **DeLong**, que estima a variância a partir dos midranks.
6. **Calibração:** regressão isotônica ajustada dentro da validação; avaliação por Brier e
   Hosmer-Lemeshow.
7. **Verificação da reimplementação:** cada modelo exportado é reexecutado em Python puro,
   percorrendo as árvores nó a nó como o JavaScript faz, e comparado ao scikit-learn.
   **Erro máximo: 3 × 10⁻⁶.**

---

## 5. Resultados — caracterização clínica e molecular

### 5.1 Perfil clínico por subtipo

| Subtipo | n (%) | Idade | Tamanho (mm) | G3 (%) | N+ (%) | ER+ (%) | HER2+ (%) | NPI | Óbitos (%) |
|---|---|---|---|---|---|---|---|---|---|
| LumA | 679 (35,8) | 63,2 | 20,5 | 25,9 | 42,3 | 99,6 | 3,1 | 3,08 | 53,6 |
| LumB | 461 (24,3) | 66,2 | 25,0 | 55,6 | 51,4 | 100,0 | 9,1 | 4,05 | 65,7 |
| Her2 | 220 (11,6) | 58,8 | 25,0 | 73,7 | 57,3 | 42,3 | 56,8 | 4,08 | 70,5 |
| Basal | 199 (10,5) | 54,5 | 25,0 | 90,5 | 53,8 | 13,6 | 10,1 | 5,02 | 55,8 |
| claudin-low | 198 (10,4) | 58,5 | 20,0 | 67,9 | 47,5 | 38,4 | 7,6 | 4,05 | 44,9 |
| Normal | 140 (7,4) | 57,8 | 23,0 | 34,1 | 40,0 | 85,7 | 9,3 | 4,03 | 54,3 |

As 14 variáveis clínicas testadas diferiram entre subtipos após correção BH. As mais
discriminantes: status de ER (FDR = 4,0 × 10⁻²³³), PR (3,9 × 10⁻¹⁰⁶), HER2 (2,8 × 10⁻⁹⁷) e grau
histológico (5,1 × 10⁻⁸⁰).

### 5.2 Sobrevida por subtipo e por estágio

| Subtipo | n | Óbitos | Mediana de OS |   | Estágio | n | Óbitos | Mediana de OS |
|---|---|---|---|---|---|---|---|---|
| claudin-low | 198 | 89 | 219,2 meses |  | I | 479 | 216 | 227,9 meses |
| LumA | 679 | 364 | 186,6 meses |  | II | 797 | 480 | 140,8 meses |
| Normal | 140 | 76 | 158,5 meses |  | III–IV | 124 | 94 | 55,5 meses |
| Basal | 199 | 111 | 130,9 meses |  | | | | |
| LumB | 461 | 303 | 123,0 meses |  | | | | |
| Her2 | 220 | 155 | 104,0 meses |  | | | | |

Log-rank por subtipo: **p = 3,74 × 10⁻¹⁰**. Por estágio: **p = 1,21 × 10⁻²⁵** — o estágio clínico
separa a coorte com força muito superior à do subtipo molecular.

O METABRIC é uma coorte majoritariamente **pré-trastuzumabe**, o que explica o Her2 aparecer com a
pior sobrevida mediana. O Basal, apesar do perfil clínico mais agressivo, não ocupa a última
posição — reflexo do padrão conhecido de risco precoce alto seguido de platô, enquanto os
luminais mantêm risco tardio persistente. Esse cruzamento de curvas é também a razão pela qual o
pressuposto de riscos proporcionais é problemático nesta coorte.

### 5.3 Genes que definem cada subtipo

**568 associações gene–subtipo** com FDR < 0,01 **e** |Δz| ≥ 0,50 (figuras `F02`, `F03`).

| Subtipo | Marcadores principais (Δz) |
|---|---|
| **Basal** | ↑ *CCNE1* (+1,67), *CDKN2A* (+1,66), *CHEK1* (+1,57), *E2F3* (+1,51) · ↓ *TGFB3* (−1,41) |
| **Her2** | ↑ ***ERBB2*** (+1,57), *ARRDC1* (+1,12), *GSK3B* (+1,04) · ↓ *SMAD4*, *MYC*, *BCL2* |
| **LumA** | ↑ *GATA3* (+1,13), *MAPT* (+1,12), *BCL2* (+0,91) · ↓ *AURKA*, *CCNE1*, *E2F2* |
| **LumB** | ↑ *GATA3* (+1,00) · ↓ *EGFR* (−1,06), *NOTCH1*, *FOXO1*, *SPRY2* |
| **claudin-low** | ↑ *FOLR2* (+1,57), *KLRG1* (+1,56), *CSF1R* (+1,37) · ↓ *ERBB3*, *RAB25* |
| **Normal** | ↑ *NR2F1*, *ABCB1*, *SPRY2* · ↓ *CDK1*, *CHEK1*, *AURKA* |

A superexpressão de *ERBB2* exatamente no subtipo Her2 funciona como **controle positivo interno**:
se a pipeline não recuperasse esse achado, haveria erro de processamento em algum ponto.

### 5.4 Mutações somáticas

| Gene | Global | LumA | LumB | Her2 | Basal | claudin-low |
|---|---|---|---|---|---|---|
| PIK3CA | 41,8% | **57,4%** | 34,9% | 41,4% | 16,1% | 24,2% |
| TP53 | 34,7% | 11,9% | 24,3% | 70,0% | **88,4%** | 52,0% |
| MUC16 | 17,2% | 16,6% | 16,1% | 27,3% | 23,1% | 9,6% |
| GATA3 | 12,1% | 19,6% | 13,9% | 7,7% | 0,0% | 2,0% |
| MAP3K1 | 10,3% | 16,2% | 8,7% | 6,8% | 4,5% | 4,0% |

O eixo *TP53*-mutado / *PIK3CA*-selvagem (Basal, Her2) contra *PIK3CA*-mutado / *TP53*-selvagem
(LumA) reproduz a dicotomia estabelecida na literatura.

Associação com sobrevida: *GATA3* mutado HR 0,590 (FDR 3,3 × 10⁻⁵), *CBFB* 0,501, *TP53* 1,289.
**Ressalva:** essas mutações concentram-se em subtipos específicos, de modo que os HR são
fortemente confundidos pelo subtipo.

---

## 6. Resultados — genes associados à sobrevida

Cox univariado nos 489 genes, HR por +1 desvio-padrão de expressão, correção Benjamini-Hochberg.
**125 genes significativos (FDR < 0,05): 61 protetores e 64 de risco.** Figura `F06`.

### 6.1 Genes associados a MAIOR sobrevida (HR < 1)

| Gene | HR (OS) | IC 95% | HR (DSS) | Comentário |
|---|---|---|---|---|
| **IGF1** | 0,804 | 0,750–0,862 | 0,820 | fator de crescimento; efeito protetor consistente |
| **STAT5A** | 0,807 | 0,759–0,858 | 0,753 | selecionado em 60% das repetições do LASSO |
| **SPRY2** | 0,810 | 0,759–0,864 | 0,864 | regulador negativo de RTK |
| **CDKN2C** | 0,815 | 0,760–0,874 | 0,874 | inibidor de CDK |
| **FLT3** | 0,820 | 0,767–0,876 | 0,792 | 3º no ranking do Random Forest |
| **ABCB1** | 0,822 | 0,770–0,877 | 0,875 | transportador de efluxo |
| LAMA2 | 0,836 | 0,788–0,887 | 0,855 | matriz extracelular |
| STAT5B | 0,838 | 0,790–0,890 | 0,744 | par de *STAT5A* |
| PDGFRA | 0,844 | 0,796–0,895 | 0,945 | HR 0,695 dentro do LumA (p = 5,4 × 10⁻¹¹) |
| CCND2 | 0,846 | 0,796–0,899 | 0,898 | ciclina D2 |
| BCL2 | 0,849 | 0,801–0,900 | 0,723 | marcador de diferenciação luminal |
| KIT | 0,849 | 0,796–0,905 | 0,868 | receptor tirosina-quinase |

### 6.2 Genes associados a MENOR sobrevida (HR > 1)

| Gene | HR (OS) | IC 95% | HR (DSS) | Comentário |
|---|---|---|---|---|
| **GSK3B** | 1,225 | 1,159–1,295 | 1,326 | selecionado em **72%** das repetições do LASSO |
| **AURKA** | 1,191 | 1,123–1,262 | **1,458** | 1º no Random Forest; proliferação |
| **VEGFA** | 1,165 | 1,100–1,233 | 1,209 | angiogênese |
| **FANCD2** | 1,163 | 1,096–1,233 | 1,323 | reparo de DNA |
| KRAS | 1,153 | 1,086–1,224 | 1,199 | oncogene clássico |
| E2F7 | 1,148 | 1,084–1,215 | 1,276 | ciclo celular |
| TUBB4B | 1,144 | 1,078–1,214 | 1,227 | citoesqueleto |
| SLC19A1 | 1,142 | 1,078–1,209 | 1,239 | transportador de folato reduzido |
| MMP11 | 1,134 | 1,067–1,206 | 1,162 | remodelamento de matriz |
| CTCF | 1,132 | 1,068–1,200 | 1,166 | organização da cromatina |
| RPS6KB2 | 1,131 | 1,066–1,200 | — | via mTOR |
| MAML1 | 1,128 | 1,062–1,197 | — | via Notch |

**Os HR são maiores na sobrevida específica do câncer**, como esperado ao remover o ruído da
mortalidade competitiva: *AURKA* passa de 1,191 (OS) para 1,458 (DSS).

### 6.3 Após ajuste clínico

Ajustando por idade, grau, tamanho, linfonodos, subtipo e os três tratamentos, **10 dos 20 genes
mais fortes mantêm significância**: *STAT5A* (0,867), *GSK3B* (1,133), *VEGFA* (1,125), *BCL2*
(0,885), *CDKN2C* (0,896), *CIR1* (0,911), *STAT5B* (0,912), *IGF1* (0,910), *AURKA* (1,107),
*FLT3* (0,923). A atenuação (*GSK3B* de 1,225 para 1,133) mostra que parte do efeito univariado é
mediada pelo subtipo e pelo grau — mas não toda.

**Alerta metodológico:** entre os 20 genes mais fortes, **12 violam o pressuposto de riscos
proporcionais**. Com seguimento mediano de quase 10 anos, os HR devem ser lidos como efeitos
médios ao longo do seguimento, não como razões constantes.

### 6.4 Onde o sinal está — por subtipo e por estágio

| Subtipo | n / eventos | Genes significativos |   | Estágio | n / eventos | Genes significativos |
|---|---|---|---|---|---|---|
| **LumA** | 679 / 364 | **12 de 12** |  | I | 479 / 216 | **1** (*NRIP1* 0,729) |
| LumB | 461 / 303 | 4 |  | **II** | 797 / 480 | **71** |
| Normal | 140 / 76 | 4 |  | III–IV | 124 / 94 | **1** (*DNAH11* 1,523) |
| claudin-low | 198 / 89 | 0 |  | | | |
| **Her2** | 220 / 155 | **0** |  | | | |
| **Basal** | 199 / 111 | **0** |  | | | |

No estágio II destacam-se *ACVR1B* (1,279), *MEN1* (1,242), *MAML1* (1,236) e *AR* (1,217) como
risco; *PDGFRA* (0,791), *STAT5A* (0,787), *IGF1* (0,762) e *CDKN2C* (0,777) como proteção.

**A informação transcricional é útil onde existe heterogeneidade de risco a resolver.** No estágio
I o prognóstico já é bom para quase todas; no III–IV já é ruim para quase todas. É no estágio II —
o maior e mais heterogêneo — que a expressão gênica efetivamente estratifica.

### 6.5 Integração: o eixo dominante é a proliferação

Cruzando marcadores de subtipo com genes prognósticos, 223 pares são simultaneamente as duas
coisas. A simetria entre Basal e LumA — os mesmos genes (*AURKA*, *FANCD2*, *E2F2*, *BCL2*) com
direções invertidas — indica que o eixo prognóstico dominante não é um conjunto de genes
independentes, mas um **gradiente único de proliferação**. Isso explica por que a assinatura
multigênica acrescenta tão pouco a um modelo que já contém subtipo e grau: ela mede, por um
caminho mais longo e instável, algo que o grau histológico já mede.

---

## 7. Resultados — aprendizado de máquina

### 7.1 Classificação dos subtipos moleculares

Onze algoritmos, 5-fold × 3 repetições, a partir dos 489 genes:

| Modelo | Acurácia | Acurácia balanceada | F1-macro |
|---|---|---|---|
| **Gradient Boosting** | **0,7872** | 0,7256 | 0,7430 |
| SVM RBF | 0,7760 | 0,7241 | 0,7385 |
| Random Forest | 0,7663 | 0,6780 | 0,7021 |
| Rede neural (MLP) | 0,7591 | 0,7166 | 0,7266 |
| Regressão logística | 0,7582 | 0,7183 | 0,7265 |
| Extra Trees | 0,7563 | 0,6595 | 0,6818 |
| LDA | 0,7412 | 0,6895 | 0,7044 |
| SVM linear | 0,7401 | 0,7107 | 0,7127 |
| Naive Bayes | 0,7062 | 0,6968 | 0,6808 |
| k-NN (k = 15) | 0,6375 | 0,5584 | 0,5923 |
| Baseline | 0,3579 | 0,1667 | 0,0879 |

**O ranking é real? Só em parte.** Pelo t pareado corrigido de Nadeau-Bengio, o Gradient Boosting
**não é estatisticamente superior** ao SVM RBF (p = 0,254) nem ao Random Forest (p = 0,095).
McNemar sobre as 1.897 predições: 130 acertos exclusivos do primeiro contra 114 do segundo,
χ² = 0,922, **p = 0,337**. A vantagem de 1,1 ponto está dentro da flutuação amostral.

São conclusões sólidas: todos superam massivamente o baseline; o k-NN é inadequado (maldição da
dimensionalidade em 489 dimensões); e os modelos lineares ficam 3–5 pontos atrás dos não-lineares,
indicando estrutura de interação entre genes.

**Genes mais informativos:** *GATA3*, *AURKA*, *EGFR*, *CDK1*, *CDC25A*, *MAPT*, *CHEK1*, *IGF1R*,
*BCL2*, *E2F2*. O Random Forest, método multivariado e não-linear, chega aos mesmos genes que o
limma identificou gene a gene.

### 7.2 Predição de óbito em 10 anos — os 16 modelos

Busca de hiperparâmetros aninhada, IC 95% de DeLong (figuras `F01`, `F02`):

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

### 7.3 Como funciona cada algoritmo

**Regressão logística.** Fronteira linear: cada variável recebe um peso, os pesos são somados e a
soma vira probabilidade pela função logística. É o modelo mais transparente e serve de referência:
se um método complexo não a supera, a estrutura dos dados é essencialmente linear.

**Random Forest.** Centenas de árvores de decisão, cada uma treinada num reamostreio dos pacientes
e enxergando apenas um sorteio das variáveis a cada divisão. A predição é a média. Essa dupla
aleatoriedade descorrelaciona os erros: cada árvore erra de um jeito diferente e a média cancela
boa parte do ruído.

**Extra Trees.** Variante em que os pontos de corte são sorteados em vez de otimizados. Mais
aleatoriedade aumenta o viés de cada árvore mas reduz ainda mais a variância do conjunto — por
isso é o segundo melhor no conjunto combinado, com 489 variáveis de sinal fraco.

**Gradient Boosting.** As árvores são construídas em sequência: cada nova árvore corrige o erro
residual das anteriores e entra com peso pequeno (taxa de aprendizado). Ganha precisão, mas é mais
sensível a sobreajuste que a floresta.

**XGBoost.** Gradient boosting com freios: regularização explícita dos pesos das folhas, poda por
ganho mínimo, amostragem de linhas e colunas a cada rodada e busca de cortes por histograma. É o
que permite rodar boosting em 500 colunas sem decorar o treino — e por isso lidera.

**Rede neural (MLP).** Camada oculta com ativação não-linear e saída logística, aprendendo
combinações que nenhuma variável expressa isoladamente. Lidera no conjunto clínico e é a **pior**
nos genes: com 489 entradas e 1.560 casos, decora em vez de generalizar.

### 7.4 A expressão gênica perde para 6 variáveis clínicas

Comparação pareada genes × clínico **dentro de cada algoritmo** (t de Nadeau-Bengio, correção FDR):

| Algoritmo | AUC genes | AUC clínico | Diferença | FDR |
|---|---|---|---|---|
| Rede neural (MLP) | 0,614 | 0,736 | −0,122 | 0,0003 |
| k-NN | 0,602 | 0,714 | −0,112 | 0,0003 |
| Regressão logística | 0,604 | 0,717 | −0,112 | 0,0003 |
| SVM linear | 0,606 | 0,715 | −0,109 | 0,0005 |
| LDA | 0,607 | 0,715 | −0,109 | 0,0005 |
| SVM RBF | 0,655 | 0,727 | −0,072 | 0,0021 |
| Extra Trees | 0,651 | 0,721 | −0,070 | 0,0005 |
| Random Forest | 0,656 | 0,726 | −0,070 | 0,0016 |
| Naive Bayes | 0,630 | 0,693 | −0,064 | 0,0011 |
| Gradient Boosting | 0,638 | 0,700 | −0,061 | 0,0025 |
| Regressão logística LASSO | 0,654 | 0,714 | −0,060 | 0,0006 |

**Onze de onze algoritmos, todos com FDR < 0,003, na mesma direção.** Não existe interpretação
alternativa: a desvantagem dos 489 genes frente a 6 variáveis clínicas não é artefato de escolha
de modelo, de linearidade, de regularização ou de hiperparâmetro. É propriedade dos dados.

**Qualificação importante.** O XGBoost sobre o conjunto combinado (0,7488) é o **único** modelo a
superar o melhor puramente clínico. Isso qualifica — sem anular — a conclusão acima: a regressão
logística não extraía nada dos genes além do clínico, mas um método capaz de capturar interações
extrai **algum** sinal. Duas ressalvas: o IC se sobrepõe amplamente ao do MLP clínico, e o ganho
custa 489 variáveis adicionais.

### 7.5 A assinatura multigênica que não funcionou

O LASSO-Cox selecionou 65 genes em `lambda.min` e 12 em `lambda.1se`. Em 25 partições independentes:

| Modelo | C-index médio | DP | Mín | Máx |
|---|---|---|---|---|
| Assinatura gênica | 0,5549 | 0,0354 | 0,500 | 0,611 |
| Clínico | 0,6553 | 0,0139 | 0,627 | 0,689 |
| Combinado | 0,6556 | 0,0142 | 0,634 | 0,691 |

Comparações pareadas (Wilcoxon): combinado *vs.* clínico → +0,0003, **p = 0,70**; assinatura
*vs.* clínico → −0,1004, **p = 6 × 10⁻⁸**.

Na amostra de derivação, o teste de razão de verossimilhança sugeria ganho espetacular
(χ² = 161,5; p = 5,3 × 10⁻³⁷). Fora da amostra, esse ganho **desaparece por completo** — a
assinatura digital do sobreajuste. Relatar apenas o valor-p intra-amostral produziria a conclusão
oposta e falsa.

**Estabilidade de seleção:** apenas quatro genes foram escolhidos com consistência — *GSK3B* (72%),
*STAT5A* (60%), *FLT3* (32%), *SPRY2* (32%); todos os demais em ≤ 12%.

### 7.6 Calibração — discriminação e calibração são coisas distintas

Antes da correção isotônica, o modelo atribuía 6,8% de risco a um grupo onde morriam 25,1%:
**Hosmer-Lemeshow χ² = 779,4; p = 5,7 × 10⁻¹⁶³**. Depois da correção ajustada dentro da validação
(figura `F03`):

| Modelo | Brier | Hosmer-Lemeshow |
|---|---|---|
| Clínico | 0,211 | **p = 0,66** |
| Genes | 0,233 | p = 0,25 |
| Combinado | 0,223 | p = 0,20 |

**A ordenação de risco, essa sim, sempre funcionou** (figura `F05`):

| Tercil | n | Óbitos em 10 anos | Sobrevida global mediana |
|---|---|---|---|
| Baixo | 535 | 26,2% | 217,6 meses |
| Médio | 505 | 45,7% | 142,4 meses |
| Alto | 520 | 63,6% | 85,7 meses |

**Log-rank p = 1,8 × 10⁻³⁶** — 132 meses de diferença entre os extremos.

---

## 8. Resultados — onde o modelo falha

### 8.1 A falha do classificador é do rótulo

Matriz de confusão (figura `F08`), recall por classe: LumA 0,872 · LumB 0,818 · Basal 0,784 ·
claudin-low 0,788 · Her2 0,723 · **Normal 0,350**.

Qui-quadrado entre subtipo real e acerto: **χ² = 195,3; p = 2,9 × 10⁻⁴⁰** — o erro não é aleatório.

**O erro tem uma geometria.** Todos os erros relevantes ocorrem entre pares **biologicamente
adjacentes**: LumA ↔ LumB (124 erros), Normal → LumA (67), Her2 → LumB/LumA (50), claudin-low →
Basal (20). O que praticamente **não** acontece é igualmente informativo: apenas 5 casos Basal
foram chamados de LumA. **O modelo nunca confunde os extremos do espectro.**

**Duas evidências de que a falha é do rótulo:**

| Variável | Teste | FDR | Associada ao erro? |
|---|---|---|---|
| **Celularidade** | qui-quadrado | **0,0111** | **sim** |
| **Integrative cluster** | qui-quadrado | **0,0003** | **sim** |
| Status de ER | qui-quadrado | 0,064 | limítrofe |
| Tamanho tumoral | Mann-Whitney | 0,137 | não |
| Grau histológico | qui-quadrado | 0,620 | não |
| Status HER2 | qui-quadrado | 1,000 | não |
| Idade | Mann-Whitney | 0,701 | não |
| Carga mutacional | Mann-Whitney | 0,701 | não |

Taxa de acerto por celularidade: **alta 82,1% → moderada 76,6% → baixa 72,7%**. Quanto menor a
fração de células tumorais, maior a contaminação por tecido normal. **O modelo não está errando o
cálculo — está recebendo, em parte, o transcriptoma do tecido errado.**

### 8.2 O modelo sabe quando não sabe

Confiança mediana nos acertos: **0,581**; nos erros: **0,442** — Mann-Whitney **p = 4,2 × 10⁻⁶⁰**.

| Faixa de confiança | n | Acerto observado |
|---|---|---|
| ≤ 0,4 | 324 | **47,5%** |
| 0,4 – 0,5 | 449 | 68,6% |
| 0,5 – 0,6 | 414 | 79,0% |
| 0,6 – 0,7 | 350 | 89,7% |
| 0,7 – 0,8 | 239 | 94,6% |
| 0,8 – 0,9 | 106 | 98,1% |
| > 0,9 | 15 | **100,0%** |

Isso habilita uma estratégia operacional: classificar automaticamente os casos com confiança acima
de 0,7 e encaminhar à revisão os 324 casos abaixo de 0,4, onde o modelo acerta pouco mais que uma
moeda.

### 8.3 Onde a predição de sobrevida não é confiável

AUC por estrato com IC 95% bootstrap (figura `F04`):

| Estrato | n | Eventos | AUC | IC 95% | Veredito |
|---|---|---|---|---|---|
| Normal-like | 111 | 52 | 0,719 | 0,613–0,813 | confiável |
| claudin-low | 150 | 59 | 0,714 | 0,621–0,793 | confiável |
| LumA | 544 | 183 | 0,678 | 0,630–0,727 | confiável |
| LumB | 394 | 204 | 0,661 | 0,608–0,717 | confiável |
| Her2 | 184 | 109 | 0,608 | 0,522–0,688 | limítrofe |
| **Basal** | 177 | 95 | **0,577** | **0,491–0,659** | **não confiável** |
| Estágio I | 390 | 113 | 0,651 | 0,591–0,710 | confiável |
| Estágio II | 668 | 324 | 0,643 | 0,599–0,685 | confiável |
| **Estágio III–IV** | 109 | 78 | **0,506** | **0,387–0,626** | **não confiável** |

No subtipo Basal e no estágio III–IV o intervalo cruza 0,5: a predição não se distingue de um
sorteio. A plataforma **bloqueia a leitura do número** nesses casos.

### 8.4 O erro que mais importa

O modelo erra assimetricamente: 302 falsos positivos contra 224 falsos negativos (teste binomial
**p = 0,0008**). Em estratificação de risco, esse é o viés menos danoso.

Comparando os óbitos **não previstos** com os **previstos** (Mann-Whitney, todos FDR < 0,001):

| Variável | Óbito não previsto | Óbito previsto |
|---|---|---|
| Idade ao diagnóstico | 61,3 anos | 69,1 anos |
| Tamanho tumoral | 20 mm | 30 mm |
| Índice de Nottingham | 4,04 | 5,04 |
| Linfonodos positivos | 0 | 2 |
| Tempo até o óbito | 70,1 meses | 46,0 meses |

**Os óbitos que o modelo não antecipa são de pacientes mais jovens, com tumores menores, sem
linfonodos comprometidos e índice prognóstico favorável.** O modelo captura bem a morte que
"parece" provável pelos critérios clínicos clássicos e falha justamente nos casos que mais
precisariam de um marcador molecular. Esse é o nicho não atendido, e ele coincide exatamente com
o resultado de que os genes não acrescentam sinal — se acrescentassem, seria aqui.

---

## 9. Camada farmacogenômica

Para cada um dos 489 genes e cada um dos três tratamentos registrados, foi ajustado um modelo de
Cox contendo o gene, o tratamento e o **termo de interação gene × tratamento** (1.467 modelos),
com ajuste por idade, grau e tamanho tumoral.

| Tratamento | n (tratadas / não tratadas) | Interações com FDR < 0,05 |
|---|---|---|
| **Hormonioterapia** | 1.170 / 727 | **14** |
| Quimioterapia | 396 / 1.501 | 0 (menor FDR: *AKT1S1*, 0,071) |
| Radioterapia | 1.137 / 760 | 0 |

### Os 14 genes com interação significativa

O HR de referência é o efeito do gene em quem **não** recebeu hormonioterapia; o HR de interação é
o fator multiplicativo aplicado a quem recebeu.

| Gene | HR sem terapia | HR de interação | FDR | Efeito resultante | Leitura |
|---|---|---|---|---|---|
| **EIF4E** | 1,184 | 0,762 | 0,012 | ≈ 0,90 | risco anulado sob terapia |
| **TWIST1** | 1,205 | 0,779 | 0,013 | ≈ 0,94 | risco anulado |
| **DNAH2** | 0,802 | 1,294 | 0,013 | ≈ 1,04 | proteção anulada |
| **SETDB1** | 0,881 | 1,275 | 0,017 | ≈ 1,12 | proteção invertida |
| **PDGFRA** | 1,148 | 0,784 | 0,022 | ≈ 0,90 | risco anulado |
| **EPCAM** | 1,139 | 0,789 | 0,023 | ≈ 0,90 | risco anulado |
| **MYC** | 1,139 | 0,804 | 0,030 | ≈ 0,92 | risco anulado |
| **MTOR** | 0,891 | 1,243 | 0,030 | ≈ 1,11 | proteção invertida |
| **MAP3K13** | 0,931 | 1,238 | 0,030 | ≈ 1,15 | proteção invertida |
| **ACKR3** | 1,209 | 0,801 | 0,031 | ≈ 0,97 | risco anulado |
| **BAP1** | 0,870 | 1,230 | 0,032 | ≈ 1,07 | proteção invertida |
| **DTX3** | 0,869 | 1,213 | 0,034 | ≈ 1,05 | proteção invertida |
| **CASP6** | 1,108 | 0,809 | 0,034 | ≈ 0,90 | risco anulado |
| **MMP16** | 1,176 | 0,811 | 0,042 | ≈ 0,95 | risco anulado |

**Teste de confundimento — a verificação decisiva.** A hormonioterapia é prescrita conforme o
status de receptor de estrogênio, então uma "interação" poderia refletir apenas a diferença entre
biologias ER+ e ER−. Reanalisando com ajuste adicional por **status de ER, subtipo molecular e
quimioterapia**, os 14 permaneceram significativos (FDR entre 0,0009 e 0,013), com HR
praticamente inalterados (*EIF4E* 0,762 → 0,790; *MYC* 0,804 → 0,794; *MTOR* 1,243 → 1,223).

**Interpretação.** Emergem dois grupos coerentes com a biologia da resistência endócrina:

- Genes de **tradução, proliferação e transição epitélio-mesenquimal** (*EIF4E*, *MYC*, *TWIST1*,
  *EPCAM*, *PDGFRA*, *MMP16*, *ACKR3*): alta expressão marca pior prognóstico em pacientes não
  tratadas, e a desvantagem desaparece sob terapia — compatível com esses tumores serem os que
  mais se beneficiam do bloqueio endócrino.
- Genes de **remodelamento de cromatina e sinalização mTOR** (*MTOR*, *SETDB1*, *MAP3K13*, *BAP1*,
  *DTX3*): alta expressão é neutra ou favorável sem terapia e passa a marcar **pior** desfecho sob
  hormonioterapia — perfil candidato a **resistência endócrina**, biologicamente plausível dado
  que a via PI3K/AKT/mTOR é o mecanismo de escape endócrino mais estabelecido no luminal.

**Nomenclatura.** Trata-se de farmacogenômica **somática** — expressão tumoral como modificadora
de efeito do tratamento — e não da farmacogenômica clássica de polimorfismos germinativos
(*CYP2D6*, *DPYD*, *UGT1A1*), para a qual o METABRIC não tem dados.

---

## 10. Plataforma de inferência

Aplicação estática de página única, com todo o cálculo no navegador — nenhum dado do usuário
trafega pela rede. Os 16 modelos foram exportados como coeficientes e arrays compactos de nós de
árvore, e a predição reimplementada em JavaScript puro.

### Algoritmo de inferência (para reimplementar em qualquer linguagem)

```
1. Monte o vetor x na ordem exata de modelo["variaveis"]:
   - nomes iniciados por "sub_"  → 1 se for o subtipo do paciente, 0 caso contrário
   - variáveis clínicas          → o valor informado (NUNCA zero quando ausente)
   - demais nomes                → o z-score do gene (0 = média da coorte)

2. Conforme modelo["tipo"]:
   linear        z = intercepto + Σ ((x[i] − media[i]) / escala[i]) × coef[i];  p = sigmoide(z)
   floresta      p = média, sobre as árvores, do valor da folha alcançada
   boosting_sk   z = base + lr × Σ (folha de cada árvore);      p = sigmoide(z)
   boosting_xgb  z = logit(base) + Σ (folha de cada árvore);    p = sigmoide(z)
   mlp           padronize; aplique W/b com ReLU camada a camada; saída logística

   Percurso da árvore: comece no nó 0; enquanto l[i] ≠ −1, vá para l[i] se x[f[i]] < t[i],
   senão para r[i]. Devolva v[i].
   ATENÇÃO: a folha é identificada por l[i] = −1. No sklearn, f[i] vale −2 nas folhas.

3. Calibre: interpole linearmente o p bruto na curva modelo["calibracao"].
```

Para o **XGBoost**, `base` é a prevalência real (≈ 0,45), lida de `booster.save_config()`.

### Decisões de projeto derivadas da análise

- **Estados de indisponibilidade explícitos.** Sem expressão carregada ou com campo clínico em
  branco, a tela mostra "Aguardando dados" e explica o que falta — em vez de calcular sobre um
  paciente médio inexistente. Campo em branco nunca vira zero.
- **O painel de confiabilidade é o elemento central**, não o número do risco. Exibe a barra de AUC
  com IC 95% por subtipo, com a linha do acaso marcada, e emite veredito explícito. Em Basal, a
  mensagem é literal: ignore o número acima.
- **Seletor de conjunto e de algoritmo**, com AUC e IC visíveis no momento da escolha, e curvas
  ROC sobrepostas.
- **O classificador de subtipo informa a própria confiança** e avisa abaixo de 0,6.

---

## 11. Erros cometidos e corrigidos

Registro explícito, porque cada um afetaria as conclusões:

**1. Vazamento de seleção de variáveis.** A primeira versão da exportação selecionava a assinatura
de 106 genes usando toda a base e só depois validava. Efeito:

| Métrica | Com vazamento | Corrigido |
|---|---|---|
| AUC do modelo de genes | 0,727 | **0,639** |
| AUC do modelo combinado | 0,768 | **0,692** |
| Basal | 0,710 → "confiável" | **0,577 (IC 0,491–0,659) → não confiável** |

Se tivesse passado, a plataforma declararia confiável justamente a predição em tumores basais.

**2. Laço infinito na travessia de árvores.** As folhas do sklearn usam `feature = −2`, não `−1`;
o teste de parada nunca disparava. Diagnosticado com `py-spy` após 13 minutos de processo travado.
A condição correta é `children_left == −1`.

**3. XGBoost com erro de 4,8 × 10⁻² na reimplementação.** O `base_score` real é a prevalência
(≈ 0,45), não 0,5, e os limiares precisam da precisão de float32 do booster. Após corrigir ambos,
o erro caiu para 8,4 × 10⁻⁸.

**4. Sete incoerências na interface auditada.** Modelo padrão contradizendo a própria nota da tela;
veredito de confiabilidade de um modelo exibido para outro; predição de 40,4% exibida sem dado
algum carregado; campo clínico vazio virando zero (risco caía de 40,3% para 15,7%); subtipo
indefinido caindo silenciosamente em LumA; afirmação falsa sobre "106 genes + PAM50"; e perda da
CSP restritiva na migração. Todas corrigidas na versão final.

---

## 12. Limitações

1. **Sem validação externa.** Toda a validação foi interna. Confirmação em TCGA-BRCA ou coorte
   própria é o passo obrigatório seguinte.
2. **Cobertura gênica parcial:** 489 genes de painel curado, não o transcriptoma completo.
3. **Riscos proporcionais violados** em 12 dos 20 genes mais fortes.
4. **Coorte histórica**, majoritariamente pré-trastuzumabe, com esquemas adjuvantes heterogêneos e
   sem registro de duração ou adesão à hormonioterapia.
5. **Tratamento não randomizado:** interação gene × tratamento em dados observacionais é vulnerável
   a confundimento por indicação.
6. **Estágio ausente em 26,2%** e apenas 124 casos em III–IV.
7. **Mortalidade competitiva:** a OS inclui 480 óbitos por outras causas; o tratamento formal
   exigiria modelos de risco competitivo (Fine-Gray).
8. **Grupo "Normal-like"** provavelmente reflete contaminação por tecido normal.
9. **z-scores calculados sobre a coorte inteira:** aplicar a uma paciente nova exige recalibração.
10. **Sem dados de genótipo germinativo:** a camada farmacogenômica é somática, não clássica.
11. **Discriminação moderada (≈ 0,75)**, comparável ao que o estadiamento clínico já entrega.

---

## 13. Próximos passos

Em ordem de retorno esperado:

1. **Restringir aos tumores luminais em estágio II** e rederivar a assinatura nesse contexto — é
   onde todo o sinal se concentra.
2. **Validar externamente** em TCGA-BRCA, com atenção à diferença de plataforma (microarranjo
   vs. RNA-seq).
3. **Aprofundar o eixo mTOR/*SETDB1*** como marcador de resistência endócrina, cruzando com as
   mutações de *PIK3CA* já disponíveis nesta coorte.
4. **Obter a matriz de expressão completa** (~24.000 genes) e repetir os módulos de expressão
   diferencial, Cox e assinatura.
5. **Modelos de risco competitivo** (Fine-Gray) e **Random Survival Forest**, que usa o tempo até
   o evento em vez de dicotomizar em 10 anos.
6. **Modelagem tempo-dependente** para os genes que violam riscos proporcionais.

---

## 14. Inventário do pacote

| Pasta | Conteúdo |
|---|---|
| `00_dados_entrada/` | `METABRIC_RNA_Mutation.csv` (8,4 MB) e `MD5SUM.txt` |
| `01_pipeline_R/` | `metabric_pipeline.R` — 12 módulos, autossuficiente |
| `02_analises_complementares/` | Estágio, farmacogenômica, Random Forest complementar (Python) |
| `03_benchmark_ml/` | Benchmark de 11 algoritmos com testes estatísticos (Python) |
| `04_treino_e_exportacao/` | Treino dos 16 modelos e exportação para o navegador (Python) |
| `05_figuras/` | 8 figuras de validação + script gerador |
| `06_tabelas/` | Genes de sobrevida (G01–G03) e métricas dos 16 modelos (G04) |
| `07_plataforma/` | Aplicação web de inferência — abra `index.html` |
| `08_resultados_brutos/` | 33 tabelas do R, 13 figuras do R, todas as saídas intermediárias |
| `09_relatorios/` | Este relatório e os relatórios temáticos anteriores |
| `EXECUTAR_TUDO.sh` | Executa o pipeline completo na ordem correta |
| `requisitos.txt` | Dependências de Python e R |

### Figuras de validação

| Figura | Conteúdo |
|---|---|
| `F01_curvas_roc.png` | Curvas ROC dos 6 algoritmos, nos 3 conjuntos |
| `F02_auc_ic95.png` | AUC com IC 95% de DeLong para os 16 modelos |
| `F03_calibracao.png` | Risco previsto × observado após correção isotônica |
| `F04_confiabilidade_estratos.png` | AUC por subtipo e estágio, com linha do acaso |
| `F05_tercis_risco.png` | Desfecho e sobrevida real por tercil |
| `F06_genes_forest.png` | Forest plot dos 24 genes mais fortes |
| `F07_genes_ml_e_estagio.png` | Importância no Random Forest e genes do estágio II |
| `F08_matriz_confusao_subtipo.png` | Matriz de confusão dos 6 subtipos |

### Bibliotecas e finalidade

**R:** data.table (leitura de alto desempenho), dplyr/tidyr/purrr/stringr/forcats (manipulação),
**limma** (expressão diferencial com moderação bayesiana empírica), **survival** (Kaplan-Meier,
log-rank, Cox, `cox.zph`), survminer (curvas com tabela de risco), **glmnet** (Cox penalizado por
LASSO), **randomForest** (classificação multiclasse), matrixStats, broom, ggplot2, ggrepel,
scales, RColorBrewer, viridis, pheatmap, patchwork, cowplot, gridExtra.

**Python:** **scikit-learn** (Random Forest, Extra Trees, Gradient Boosting, MLP, SVM, validação
cruzada, calibração isotônica, métricas), **XGBoost** (boosting regularizado), **lifelines**
(modelos de Cox, Kaplan-Meier, log-rank), statsmodels (Benjamini-Hochberg, McNemar), scipy
(testes estatísticos), pandas/numpy, matplotlib (figuras).

---

## 15. Referências

- Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*. 2012;486(7403):346-352.
- Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. *Nat Commun*. 2016;7:11479.
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47.
- Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. *J Stat Softw*. 2011;39(5):1-13.
- Breiman L. Random Forests. *Mach Learn*. 2001;45(1):5-32.
- Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *KDD*. 2016:785-794.
- DeLong ER, DeLong DM, Clarke-Pearson DL. Comparing the areas under two or more correlated receiver operating characteristic curves. *Biometrics*. 1988;44(3):837-845.
- Nadeau C, Bengio Y. Inference for the generalization error. *Mach Learn*. 2003;52(3):239-281.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*. 1995;57(1):289-300.
- Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. Springer; 2000.
- Meinshausen N, Bühlmann P. Stability selection. *J R Stat Soc B*. 2010;72(4):417-473.
- Niculescu-Mizil A, Caruana R. Predicting good probabilities with supervised learning. *ICML*. 2005:625-632.

---

*Uso restrito a pesquisa e ensino. Não é dispositivo médico e não deve orientar decisão clínica individual.*
