# Farmacogenômica e prognóstico molecular no câncer de mama

## Relatório integrado da coorte METABRIC (n = 1.897): subtipos moleculares, estágio clínico, aprendizado de máquina e modificação de efeito por tratamento

**Coorte:** METABRIC / cBioPortal `brca_metabric` · **Pacientes analisadas:** 1.897 · **Seguimento mediano:** 115,6 meses
**Camadas de dados:** 489 genes de expressão (z-score), 173 genes com status mutacional, 31 variáveis clínicas, desfechos de sobrevida
**Ambiente:** R 4.3.3 (pipeline principal, 12 módulos) + Python 3 / lifelines 0.30.3 / scikit-learn 1.8.0 (análises complementares) · Semente fixa 42 em todas as etapas aleatórias

---

## 1. Sumário executivo

Este relatório reúne, em um único documento, tudo o que a análise produziu sobre a coorte METABRIC — dados **reais**, sem qualquer simulação ou imputação por modelo. A investigação foi conduzida em duas etapas: um pipeline reprodutível em R (33 tabelas, 13 figuras) e um conjunto de análises complementares em Python que acrescentou três camadas ausentes na primeira rodada — **estratificação por estágio tumoral**, **modelos de aprendizado de máquina para desfecho de sobrevida** e **teste formal de interação gene × tratamento**, que é o que autoriza a leitura propriamente farmacogenômica.

Os oito achados centrais:

1. **Os subtipos moleculares diferem fortemente em prognóstico** (log-rank p = 3,7 × 10⁻¹⁰): sobrevida global mediana de 219,2 meses no claudin-low a 104,0 meses no Her2.
2. **O estágio clínico separa a coorte com muito mais força que qualquer gene isolado** (log-rank p = 1,2 × 10⁻²⁵): mediana de 227,9 meses no estágio I, 140,8 no II e 55,5 no III–IV.
3. **568 associações gene–subtipo** significativas, com assinaturas biologicamente coerentes e um controle positivo interno perfeito (*ERBB2* superexpresso no subtipo Her2).
4. **Random Forest classifica os 6 subtipos com 76,1% de acurácia** (erro out-of-bag 23,9%) a partir apenas da expressão gênica, com *GATA3*, *EGFR*, *CDK1*, *AURKA* e *MAPT* como genes mais informativos. O erro concentra-se no grupo "Normal-like" (80,7%), o que é informativo sobre a natureza desse grupo, não um defeito do modelo.
5. **125 genes associam-se à sobrevida global e 137 à sobrevida específica por câncer.** Núcleo robusto de **pior** sobrevida: *GSK3B*, *AURKA*, *VEGFA*, *FANCD2*, *E2F2/E2F7*, *CCNE1*. Núcleo de **melhor** sobrevida: *STAT5A*, *SPRY2*, *IGF1*, *PDGFRA*, *ABCB1*, *BCL2*, *CDKN2C*, *STAT5B*, *FLT3*.
6. **O sinal prognóstico transcricional é dependente do contexto**: concentra-se nos tumores **luminais** e no **estágio II**, e praticamente desaparece em Basal, Her2, no estágio I e no estágio III–IV. Isso vale tanto para o Cox estratificado quanto para o Random Forest (AUC 0,666 nos luminais contra 0,420 no Basal).
7. **Nenhuma combinação de genes superou um modelo clínico simples fora da amostra.** A assinatura LASSO-Cox obteve C-index médio de 0,555 contra 0,655 do modelo clínico em 25 partições independentes; no Random Forest de mortalidade em 10 anos, AUC de 0,653 (genes) contra 0,730 (clínico). Este resultado negativo é relatado com destaque porque contradiz a leitura ingênua dos testes feitos apenas na amostra de derivação.
8. **Camada farmacogenômica:** 14 genes apresentam **interação significativa com hormonioterapia** (FDR < 0,05), e todos os 14 permanecem significativos após ajuste por status de ER e subtipo molecular — ou seja, não são simples reflexo de quem recebe a terapia. Para *EIF4E*, *TWIST1*, *MYC*, *PDGFRA*, *EPCAM*, *CASP6*, *ACKR3* e *MMP16*, a alta expressão é fator de risco em quem **não** recebeu hormonioterapia e deixa de sê-lo em quem recebeu; para *MTOR*, *SETDB1*, *MAP3K13*, *BAP1*, *DNAH2* e *DTX3*, o padrão se inverte. Nenhuma interação sobreviveu à correção para quimioterapia ou radioterapia.

---

## 2. Como ler os números deste relatório

Um glossário curto, para que qualquer seção possa ser lida sem consultar material externo.

| Termo | O que significa aqui |
|---|---|
| **z-score de expressão** | Expressão do gene padronizada: 0 é a média da coorte, +1 é um desvio-padrão acima. Permite comparar genes com escalas diferentes. |
| **HR (hazard ratio)** | Risco relativo de óbito por **aumento de 1 desvio-padrão** na expressão do gene. HR = 1,20 → 20% mais risco; HR = 0,80 → 20% menos risco. HR = 1 significa ausência de efeito. |
| **IC 95%** | Faixa de valores compatíveis com os dados. Se o intervalo não cruza 1, o efeito é estatisticamente detectável. |
| **FDR** | Valor-p corrigido para múltiplos testes (Benjamini-Hochberg). Testar 489 genes gera falsos positivos por acaso; o FDR controla essa proporção. |
| **Δz** | Diferença de expressão média entre um subtipo e todos os demais, em desvios-padrão. É um tamanho de efeito, não apenas um p-valor. |
| **C-index** | Probabilidade de o modelo ordenar corretamente duas pacientes quanto ao risco. 0,5 = acaso; 1,0 = perfeito. Na prática oncológica, 0,65–0,70 é um modelo útil. |
| **AUC** | Análogo do C-index para classificação binária (aqui: óbito em até 10 anos). Mesma escala de leitura. |
| **Erro out-of-bag (OOB)** | Estimativa honesta de erro do Random Forest, calculada em cada árvore com os pacientes que ela não viu no treino. |
| **OS / DSS** | Sobrevida global (óbito por qualquer causa) / sobrevida específica do câncer (óbitos por outras causas são censurados, não contados como evento). |
| **Interação gene × tratamento** | Teste de se o efeito prognóstico de um gene **muda** conforme a paciente recebeu ou não determinado tratamento. É o desenho que sustenta uma afirmação farmacogenômica. |

---

## 3. Dados, coorte e integridade

| Item | Valor |
|---|---|
| Estudo | METABRIC (Curtis et al., *Nature* 2012; Pereira et al., *Nat Commun* 2016) |
| Arquivo | `METABRIC_RNA_Mutation.csv`, MD5 `c619471beb87af0f9fd4e4a40058654d` (verificado) |
| Pacientes no arquivo bruto | 1.904 · **Coorte analítica: 1.897** |
| Colunas | 693 (31 clínicas + 489 expressão + 173 mutação) |
| IDs duplicados / ausentes na matriz de expressão | 0 / 0 |
| Plataforma | Illumina HT-12 v3 (microarranjo); painel de sequenciamento alvo para mutações |
| Óbitos por qualquer causa | 1.098 |
| Óbitos atribuídos ao câncer de mama | 619 |

**Auditoria da codificação do desfecho.** Este conjunto de dados é frequentemente analisado com o evento invertido, o que produziria conclusões exatamente opostas. A tabulação cruzada confirma que `overall_survival = 1` significa **paciente viva**; o evento de óbito foi definido como `1 − overall_survival`.

**Exclusões:** 1 caso de sarcoma mamário e 6 casos sem classificação de subtipo (`NC`). Nenhum caso perdido por tempo de seguimento inválido.

**Estágio tumoral:** disponível em 1.400 pacientes (73,8%); 497 sem registro. Estágios 0 e 1 foram agrupados como **I** (n = 479), estágio 2 como **II** (n = 797), estágios 3 e 4 como **III–IV** (n = 124), pela raridade dos extremos.

---

## 4. Métodos

### 4.1 Pipeline principal em R (12 módulos)

Parâmetros declarados *a priori* em um único objeto de configuração, evitando escolhas post-hoc: FDR < 0,01 e |Δz| ≥ 0,50 para expressão diferencial; FDR < 0,05 para modelos de Cox; partição 70/30; 500 árvores no Random Forest; 10 dobras na validação cruzada; frequência mínima de 2% para testar mutações.

| Módulo | O que faz |
|---|---|
| 0–2 | Ambiente, ingestão, verificação de integridade, auditoria do desfecho, montagem da matriz 489 × 1.897 |
| 3 | Caracterização clínica por subtipo (qui-quadrado e Kruskal-Wallis, correção BH) |
| 4 | Expressão diferencial com **limma** — modelo linear por gene, contrastes um-vs-demais, moderação bayesiana empírica |
| 5 | **Random Forest** de classificação dos 6 subtipos, erro OOB e ranking de importância |
| 6 | **Cox univariado** para os 489 genes, em OS e DSS; teste de riscos proporcionais e Cox ajustado por 8 covariáveis clínicas nos 20 genes mais fortes |
| 7 | Kaplan-Meier e **Cox estratificado dentro de cada subtipo** |
| 8 | **Assinatura multigênica LASSO-Cox** com validação em partição independente |
| 8b | **Validação repetida em 25 partições** e estabilidade de seleção dos genes |
| 9 | Mutações somáticas: heterogeneidade entre subtipos (Fisher com Monte Carlo) e associação com sobrevida |
| 10 | Integração das camadas e escore de prioridade por subtipo |

### 4.2 Análises complementares em Python

| Análise | Desenho |
|---|---|
| **Cox por estágio** | 489 modelos de Cox univariados dentro de cada grupo de estágio (I, II, III–IV), correção BH separada por estágio → tabela `C01` |
| **Interação gene × tratamento** | 1.467 modelos de Cox (489 genes × 3 tratamentos) com os termos gene, tratamento e **gene × tratamento**, ajustados por idade, grau e tamanho tumoral → `C02`; HR por estrato de tratamento → `C03`; reanálise dos genes significativos com ajuste adicional por status de ER, subtipo e quimioterapia → `C08` |
| **Random Forest de subtipo** | Replicação independente do módulo 5 em outra linguagem e outra implementação → `C04` |
| **Random Forest de mortalidade em 10 anos** | Classificação de óbito em até 120 meses (n = 1.642, 735 eventos), validação cruzada 5-fold estratificada, comparação entre modelo de genes, clínico e combinado → `C05`, `C06`; desempenho dentro de cada subtipo e estágio → `C07` |

Todas as bibliotecas utilizadas, com versão e finalidade específica, estão listadas na seção 9.

---

## 5. Resultados

### 5.1 Perfil clínico por subtipo

| Subtipo | n (%) | Idade mediana | Tamanho (mm) | G3 (%) | N+ (%) | ER+ (%) | HER2+ (%) | NPI | Óbitos (%) |
|---|---|---|---|---|---|---|---|---|---|
| LumA | 679 (35,8) | 63,2 | 20,5 | 25,9 | 42,3 | 99,6 | 3,1 | 3,08 | 53,6 |
| LumB | 461 (24,3) | 66,2 | 25,0 | 55,6 | 51,4 | 100,0 | 9,1 | 4,05 | 65,7 |
| Her2 | 220 (11,6) | 58,8 | 25,0 | 73,7 | 57,3 | 42,3 | 56,8 | 4,08 | 70,5 |
| Basal | 199 (10,5) | 54,5 | 25,0 | 90,5 | 53,8 | 13,6 | 10,1 | 5,02 | 55,8 |
| claudin-low | 198 (10,4) | 58,5 | 20,0 | 67,9 | 47,5 | 38,4 | 7,6 | 4,05 | 44,9 |
| Normal | 140 (7,4) | 57,8 | 23,0 | 34,1 | 40,0 | 85,7 | 9,3 | 4,03 | 54,3 |

As 14 variáveis clínicas testadas diferiram entre subtipos após correção BH. As mais discriminantes: status de ER (FDR = 4,0 × 10⁻²³³), status de PR (3,9 × 10⁻¹⁰⁶), HER2 (2,8 × 10⁻⁹⁷) e grau histológico (5,1 × 10⁻⁸⁰).

### 5.2 Sobrevida por subtipo

| Subtipo | n | Óbitos | Mediana de OS (meses) | IC 95% |
|---|---|---|---|---|
| claudin-low | 198 | 89 | 219,2 | 194,1 – 238,1 |
| LumA | 679 | 364 | 186,6 | 169,0 – 198,1 |
| Normal | 140 | 76 | 158,5 | 125,8 – 203,5 |
| Basal | 199 | 111 | 130,9 | 83,4 – 206,6 |
| LumB | 461 | 303 | 123,0 | 114,9 – 143,1 |
| Her2 | 220 | 155 | 104,0 | 88,9 – 142,4 |

Log-rank global: **p = 3,74 × 10⁻¹⁰**.

O METABRIC é uma coorte majoritariamente **pré-trastuzumabe**, o que explica o Her2 aparecer com a pior sobrevida mediana. O Basal, apesar do perfil clínico mais agressivo (90,5% G3, maior NPI), não ocupa a última posição em sobrevida global de longo prazo — reflexo do padrão conhecido de risco precoce alto seguido de platô, enquanto os luminais mantêm risco tardio persistente. Esse cruzamento de curvas é também a razão pela qual o pressuposto de riscos proporcionais é problemático nesta coorte (seção 5.6).

### 5.3 Sobrevida por estágio tumoral

| Estágio | n | Óbitos | Mediana de OS (meses) |
|---|---|---|---|
| I | 479 | 216 | 227,9 |
| II | 797 | 480 | 140,8 |
| III–IV | 124 | 94 | 55,5 |

Log-rank global: **p = 1,21 × 10⁻²⁵** — uma separação com ordem de grandeza muito superior à obtida por qualquer gene individual, e superior à do próprio subtipo molecular. Distribuição conjunta estágio × subtipo (tabela `C01` e seção 5.9): o estágio I concentra LumA (211/479) e o estágio III–IV é pequeno e heterogêneo, o que limita o poder estatístico nesse estrato.

### 5.4 Genes que definem cada subtipo

**568 associações gene–subtipo** com FDR < 0,01 **e** |Δz| ≥ 0,50.

| Subtipo | Marcadores principais (Δz) |
|---|---|
| **Basal** | ↑ *CCNE1* (+1,67), *CDKN2A* (+1,66), *CHEK1* (+1,57), *MAP2* (+1,53), *E2F3* (+1,51), *CDC25A* (+1,45) · ↓ *TGFB3* (−1,41) — desregulação do ciclo G1/S e resposta a dano de DNA |
| **Her2** | ↑ *ERBB2* (+1,57), *ARRDC1* (+1,12), *GSK3B* (+1,04), *AKT1* (+0,99) · ↓ *SMAD4* (−1,11), *MYC* (−0,98), *BCL2* (−0,96) |
| **LumA** | ↑ *GATA3* (+1,13), *MAPT* (+1,12), *BCL2* (+0,91) · ↓ *AURKA* (−0,96), *CCNE1* (−0,94), *E2F2* (−0,93) — baixa proliferação, diferenciação luminal preservada |
| **LumB** | ↑ *GATA3* (+1,00) · ↓ *EGFR* (−1,06), *NOTCH1* (−0,94), *FOXO1* (−0,90), *SPRY2* (−0,90) — identidade luminal com perda de reguladores negativos de crescimento |
| **claudin-low** | ↑ *FOLR2* (+1,57), *KLRG1* (+1,56), *CSF1R* (+1,37), *TGFBR2* (+1,27) · ↓ *ERBB3* (−1,51), *RAB25* (−1,42) — assinatura imune/estromal com perda de identidade epitelial |
| **Normal** | ↑ *NR2F1* (+1,05), *ABCB1* (+1,02), *SPRY2* (+0,98) · ↓ *CDK1* (−1,04), *CHEK1* (−0,97), *AURKA* (−0,94) |

A superexpressão de *ERBB2* exatamente no subtipo Her2 funciona como **controle positivo interno**: se a pipeline não recuperasse esse achado, haveria erro de processamento em algum ponto.

### 5.5 Aprendizado de máquina I — Random Forest de classificação dos subtipos

Modelo de 500 árvores treinado com os 489 genes para prever o subtipo molecular.

| Implementação | Erro OOB | Acurácia |
|---|---|---|
| R (`randomForest` 4.7-1.1) | 23,77% | 76,2% |
| Python (`scikit-learn` 1.8.0) | 23,88% | 76,1% |

A concordância entre duas linguagens e duas implementações independentes é uma verificação de robustez do resultado, não uma repetição redundante.

**Erro por classe (matriz de confusão OOB):**

| Subtipo | Erro | Para onde vão os erros |
|---|---|---|
| LumA | 11,3% | 66 → LumB |
| LumB | 18,0% | 74 → LumA |
| Basal | 20,6% | 22 → Her2 |
| claudin-low | 29,8% | 33 → Basal, 18 → LumA |
| Her2 | 35,5% | 47 → LumB, 29 → LumA |
| **Normal** | **80,7%** | **93 dos 140 → LumA** |

**Genes mais importantes** (concordantes entre as duas implementações): *GATA3*, *AURKA*, *EGFR*, *CDK1*, *CDC25A*, *MAPT*, *CHEK1*, *IGF1R*, *BCL2*, *E2F2*, *CCNE1*, *TGFBR2*, *ERBB3*, *CCNB1*, *E2F7*, *ERBB2*.

Dois pontos interpretativos importam. Primeiro, o erro de 80,7% no grupo "Normal-like" é **informação, não falha**: quase todos os casos migram para LumA, o que é consistente com a literatura que questiona se esse grupo é uma entidade biológica ou um artefato de baixa celularidade tumoral com contaminação por tecido mamário normal. Segundo, o Random Forest é um método multivariado e não-linear que **chega aos mesmos genes** que o limma identificou gene a gene — convergência entre métodos com pressupostos diferentes.

### 5.6 Genes associados à sobrevida na coorte completa

**125 genes** com FDR < 0,05 para sobrevida global e **137** para sobrevida específica do câncer. Os 15 mais fortes em OS (HR por +1 desvio-padrão de expressão):

| Gene | HR | IC 95% | FDR | Direção |
|---|---|---|---|---|
| GSK3B | 1,225 | 1,16–1,29 | 3,8 × 10⁻¹⁰ | **pior** sobrevida |
| STAT5A | 0,807 | 0,76–0,86 | 1,7 × 10⁻⁹ | melhor |
| SPRY2 | 0,810 | 0,76–0,86 | 2,6 × 10⁻⁸ | melhor |
| IGF1 | 0,804 | 0,75–0,86 | 1,0 × 10⁻⁷ | melhor |
| ABCB1 | 0,822 | 0,77–0,88 | 2,6 × 10⁻⁷ | melhor |
| AURKA | 1,191 | 1,12–1,26 | 2,6 × 10⁻⁷ | **pior** |
| LAMA2 | 0,836 | 0,79–0,89 | 2,6 × 10⁻⁷ | melhor |
| FLT3 | 0,820 | 0,77–0,88 | 2,6 × 10⁻⁷ | melhor |
| STAT5B | 0,838 | 0,79–0,89 | 4,9 × 10⁻⁷ | melhor |
| CDKN2C | 0,815 | 0,76–0,87 | 5,0 × 10⁻⁷ | melhor |
| PDGFRA | 0,844 | 0,80–0,89 | 6,2 × 10⁻⁷ | melhor |
| RPS6 | 0,855 | 0,81–0,90 | 1,6 × 10⁻⁶ | melhor |
| CCND2 | 0,846 | 0,80–0,90 | 2,3 × 10⁻⁶ | melhor |
| BCL2 | 0,849 | 0,80–0,90 | 5,0 × 10⁻⁶ | melhor |
| VEGFA | 1,165 | 1,10–1,23 | 5,8 × 10⁻⁶ | **pior** |

Na sobrevida específica do câncer os efeitos são **maiores em magnitude**, como esperado ao remover o ruído da mortalidade competitiva: *AURKA* HR 1,458 (FDR 1,5 × 10⁻²⁰), *BCL2* 0,723, *MAPT* 0,735, *GSK3B* 1,326, *FANCD2* 1,323, *E2F2* 1,310, *CCNE1* 1,271.

**Após ajuste por idade, grau, tamanho, linfonodos, subtipo e os três tratamentos**, 10 dos 20 genes testados mantiveram significância: *STAT5A* (0,867), *GSK3B* (1,133), *VEGFA* (1,125), *BCL2* (0,885), *CDKN2C* (0,896), *CIR1* (0,911), *STAT5B* (0,912), *IGF1* (0,910), *AURKA* (1,107), *FLT3* (0,923). A atenuação dos HR (por exemplo, *GSK3B* de 1,225 para 1,133) mostra que parte do efeito univariado é mediada pelo subtipo e pelo grau — mas não toda, o que caracteriza informação prognóstica parcialmente independente.

**Alerta metodológico.** Entre os 20 genes mais fortes, **12 violam o pressuposto de riscos proporcionais**. Com seguimento mediano de quase 10 anos, isso é esperado: o efeito de vários genes é forte nos primeiros anos e se dilui depois. Os HR devem ser lidos como **efeitos médios ao longo do seguimento**, não como razões de risco constantes.

### 5.7 Onde o sinal prognóstico está — por subtipo

Testando os 12 genes mais prognósticos **dentro** de cada subtipo:

| Subtipo | n / eventos | Genes significativos | Destaques |
|---|---|---|---|
| **LumA** | 679 / 364 | **12 de 12** | *PDGFRA* 0,695 (p = 5,4 × 10⁻¹¹), *SPRY2* 0,692, *IGF1* 0,727, *ABCB1* 0,741, *CDKN2C* 0,746, *STAT5A* 0,757, *GSK3B* 1,273 |
| **LumB** | 461 / 303 | 4 | *GSK3B* 1,304, *IGF1* 0,807, *STAT5A* 0,829, *FLT3* 0,798 |
| **Normal** | 140 / 76 | 4 | *AURKA* 1,835 (p = 5,8 × 10⁻⁴), *GSK3B* 1,388, *STAT5B* 0,674, *FLT3* 0,648 |
| **claudin-low** | 198 / 89 | 0 | efeitos fracos, nenhum sobrevive ao FDR |
| **Her2** | 220 / 155 | **0** | — |
| **Basal** | 199 / 111 | **0** | — |

Nos subtipos Her2 e Basal, o próprio subtipo já determina a maior parte do risco e a expressão desses genes acrescenta pouco. A consequência prática é direta: **um painel prognóstico derivado da coorte inteira será, na prática, um painel para doença luminal.**

### 5.8 Onde o sinal prognóstico está — por estágio

Cox univariado dos 489 genes dentro de cada estrato de estágio (análise complementar `C01`):

| Estágio | n / eventos | Genes com FDR < 0,05 | Principais |
|---|---|---|---|
| **I** | 479 / 216 | **1** | *NRIP1* 0,729 (FDR 0,0097). Abaixo do limiar, mas coerentes: *STAT5A* 0,786, *CDKN2C* 0,775, *MAPT* 0,808 (proteção); *FN1* 1,260, *MMP11* 1,276, *E2F7* 1,269 (risco) |
| **II** | 797 / 480 | **71** | Risco: *ACVR1B* 1,279, *MEN1* 1,242, *MAML1* 1,236, *AR* 1,217, *VEGFA* 1,181, *SMARCC2* 1,212. Proteção: *PDGFRA* 0,791, *STAT5A* 0,787, *IGF1* 0,762, *CCND2* 0,790, *CDKN2C* 0,777, *PTPN22* 0,780, *JAK2* 0,800, *LAMA2* 0,804 |
| **III–IV** | 124 / 94 | **1** | *DNAH11* 1,523 (FDR 0,039); nenhum outro sobrevive à correção |

O padrão espelha o encontrado por subtipo: **a informação transcricional é útil onde existe heterogeneidade de risco a resolver**. No estágio I o prognóstico já é bom para quase todas as pacientes; no III–IV já é ruim para quase todas (e n = 124 limita o poder); é no **estágio II** — o maior grupo e o mais heterogêneo — que a expressão gênica efetivamente estratifica.

Note também a convergência entre os dois eixos: *STAT5A*, *PDGFRA*, *IGF1* e *CDKN2C* aparecem como protetores tanto no LumA quanto no estágio II, e *VEGFA* como fator de risco em ambos.

### 5.9 Aprendizado de máquina II — Random Forest para o desfecho de sobrevida

Classificação de **óbito em até 10 anos** (n = 1.642, 735 eventos), validação cruzada 5-fold estratificada.

| Modelo | Variáveis | AUC |
|---|---|---|
| Expressão gênica | 489 genes | **0,653** |
| Clínico | idade, grau, tamanho, linfonodos, NPI, subtipo | **0,730** |
| Combinado | genes + clínico | 0,715 |

**Genes mais importantes no modelo de desfecho** (índice de Gini): *AURKA*, *STAT5A*, *FLT3*, *DIRAS3*, *STAT5B*, *GSK3B*, *MAPT*, *NCOA3*, *BCL2*, *MAP2K4*, *RUNX1*, *IGF1R*, *E2F2*, *ABCB1*, *VEGFA*, *IGF1*.

Esse ranking é notável por dois motivos. Primeiro, ele é **quase disjunto** do ranking de importância para classificar subtipos (*GATA3*, *EGFR*, *CDK1*, *CDC25A*): os genes que melhor definem *o que* o tumor é não são os mesmos que melhor predizem *o que vai acontecer* com a paciente — com a exceção instrutiva de *AURKA*, que lidera as duas listas. Segundo, ele **reproduz de forma independente** o núcleo identificado pelo Cox e pela seleção estável do LASSO (*AURKA*, *STAT5A*, *GSK3B*, *FLT3*, *BCL2*), apesar de o Random Forest não assumir linearidade nem proporcionalidade de riscos. Três famílias metodológicas distintas convergem para o mesmo conjunto de genes.

**Desempenho por estrato — a mesma dependência de contexto:**

| Estrato | n / eventos | AUC (só genes) |
|---|---|---|
| Luminais (LumA + LumB) | 986 / 410 | **0,666** |
| LumA | 574 / 196 | 0,648 |
| claudin-low | 158 / 61 | 0,596 |
| LumB | 412 / 214 | 0,580 |
| Her2 | 196 / 113 | 0,558 |
| Basal | 180 / 97 | **0,420** (abaixo do acaso) |
| Estágio I | 408 / 116 | 0,612 |
| Estágio II | 691 / 338 | 0,601 |
| Estágio III–IV | 113 / 79 | 0,549 |

O AUC de 0,420 no Basal não indica um sinal invertido explorável — indica **ausência de sinal** com flutuação amostral em um estrato pequeno. A leitura correta é que, dentro do subtipo Basal, a expressão desses 489 genes não carrega informação prognóstica utilizável.

### 5.10 Assinatura multigênica: o que ela ensina ao não funcionar

O LASSO-Cox selecionou 65 genes em `lambda.min` e 12 em `lambda.1se` (*FLT3*, *CDKN2C*, *SMAD6*, *STAT5A*, *GSK3B*, *MMP25*, *VEGFA*, *ABCB1*, *SPRY2*, *IGF1*, *CTCF*, *PDGFRA*).

**Desempenho em 25 partições independentes:**

| Modelo | C-index médio | DP | Mín | Máx |
|---|---|---|---|---|
| Assinatura gênica | 0,5549 | 0,0354 | 0,500 | 0,611 |
| Clínico | 0,6553 | 0,0139 | 0,627 | 0,689 |
| Combinado | 0,6556 | 0,0142 | 0,634 | 0,691 |

Comparações pareadas (Wilcoxon): combinado *vs.* clínico → **+0,0003, p = 0,70** (nenhum ganho); assinatura *vs.* clínico → **−0,1004, p = 6 × 10⁻⁸** (a assinatura é claramente inferior).

Na amostra de derivação, o teste de razão de verossimilhança sugeria que a assinatura acrescentava informação de forma espetacular (χ² = 161,5; p = 5,3 × 10⁻³⁷). Fora da amostra, esse ganho **desaparece por completo**. Essa discrepância é a assinatura digital do sobreajuste: com 489 candidatos e cerca de 1.100 eventos, o LASSO captura estrutura específica do conjunto de treino. Some-se a isso que o subtipo molecular já está no modelo clínico e ele próprio é um resumo da expressão gênica — a assinatura, em boa parte, apenas reexpressa informação já contida no subtipo e no grau. Relatar apenas o valor-p intra-amostral produziria a conclusão oposta e falsa.

**Estabilidade de seleção nas 25 repetições:** apenas quatro genes foram escolhidos com alguma consistência — **GSK3B (72%)**, **STAT5A (60%)**, *FLT3* (32%) e *SPRY2* (32%); todos os demais ficaram em ≤ 12%. A lista de 65 genes é em grande parte instável, mas o núcleo *GSK3B* / *STAT5A* é reprodutível e coincide com o Cox univariado, o Cox ajustado e o Random Forest de desfecho.

### 5.11 Mutações somáticas

Dos 173 genes com dado mutacional, 71 alcançaram frequência ≥ 2%.

| Gene | Global | LumA | LumB | Her2 | Basal | claudin-low | Normal |
|---|---|---|---|---|---|---|---|
| PIK3CA | 41,8% | **57,4%** | 34,9% | 41,4% | 16,1% | 24,2% | 50,0% |
| TP53 | 34,7% | 11,9% | 24,3% | 70,0% | **88,4%** | 52,0% | 22,9% |
| MUC16 | 17,2% | 16,6% | 16,1% | 27,3% | 23,1% | 9,6% | 10,0% |
| GATA3 | 12,1% | 19,6% | 13,9% | 7,7% | 0,0% | 2,0% | 8,6% |
| MAP3K1 | 10,3% | 16,2% | 8,7% | 6,8% | 4,5% | 4,0% | 10,0% |
| CDH1 | 9,1% | 12,7% | 9,8% | 5,0% | 2,0% | 4,0% | 12,9% |

O eixo *TP53*-mutado / *PIK3CA*-selvagem (Basal, Her2) contra *PIK3CA*-mutado / *TP53*-selvagem (LumA) reproduz com precisão a dicotomia estabelecida na literatura — outra validação implícita da pipeline.

Associação com sobrevida global: *GATA3* mutado HR 0,590 (FDR 3,3 × 10⁻⁵), *CBFB* mutado 0,501 (1,3 × 10⁻³), *TP53* mutado 1,289 (1,3 × 10⁻³). **Ressalva:** essas mutações concentram-se em subtipos específicos, de modo que os HR são fortemente confundidos pelo subtipo e não devem ser lidos como efeito causal independente.

### 5.12 Camada farmacogenômica — modificação do efeito prognóstico pelo tratamento

Esta é a análise que sustenta a leitura farmacogenômica do trabalho. Para cada um dos 489 genes e cada um dos três tratamentos registrados, foi ajustado um modelo de Cox contendo o gene, o tratamento e o **termo de interação gene × tratamento**, com ajuste por idade, grau e tamanho tumoral (1.467 modelos, correção BH separada por tratamento).

| Tratamento | n (tratadas / não tratadas) | Interações com FDR < 0,05 |
|---|---|---|
| **Hormonioterapia** | 1.170 / 727 | **14** |
| Quimioterapia | 396 / 1.501 | 0 (menor FDR: *AKT1S1*, 0,071) |
| Radioterapia | 1.137 / 760 | 0 (menor FDR: *PRPS2*, 0,272) |

**Os 14 genes com interação significativa com hormonioterapia.** O HR de referência é o efeito do gene em quem **não** recebeu hormonioterapia; o HR de interação é o fator multiplicativo aplicado a quem recebeu; a última coluna é o efeito resultante nas tratadas.

| Gene | HR sem hormonioterapia | HR de interação | FDR | HR resultante nas tratadas | Leitura |
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

**Teste de confundimento — a verificação decisiva.** A hormonioterapia não é atribuída ao acaso: ela é prescrita conforme o status de receptor de estrogênio. Uma "interação" poderia, portanto, ser apenas o reflexo de que tumores ER+ e ER− têm biologias diferentes. Para separar as duas explicações, os 14 genes foram reanalisados com ajuste adicional por **status de ER, subtipo molecular e quimioterapia** (tabela `C08`):

**Os 14 permaneceram significativos** (FDR entre 0,0009 e 0,013), com HR de interação praticamente inalterados — por exemplo, *EIF4E* de 0,762 para 0,790; *MYC* de 0,804 para 0,794; *MTOR* de 1,243 para 1,223. Ou seja, a modificação de efeito **não é explicada** por quem recebe a terapia nem pelo subtipo do tumor.

**Interpretação.** Emergem dois grupos coerentes com a biologia da resistência endócrina:

- Genes ligados à **tradução, proliferação e transição epitélio-mesenquimal** (*EIF4E*, *MYC*, *TWIST1*, *EPCAM*, *PDGFRA*, *MMP16*, *ACKR3*): sua alta expressão marca pior prognóstico em pacientes não tratadas com hormonioterapia, e essa desvantagem deixa de existir sob a terapia — padrão compatível com esses tumores serem justamente os que mais se beneficiam do bloqueio endócrino.
- Genes de **remodelamento de cromatina e sinalização mTOR** (*MTOR*, *SETDB1*, *MAP3K13*, *BAP1*, *DTX3*): sua alta expressão é neutra ou favorável sem terapia e passa a marcar **pior** desfecho sob hormonioterapia — um perfil candidato a **resistência endócrina**, e biologicamente plausível dado que a via PI3K/AKT/mTOR é o mecanismo de escape endócrino mais bem estabelecido no câncer de mama luminal (e alvo clínico de inibidores de mTOR).

**Ausência de sinal para quimioterapia e radioterapia.** Nenhuma interação sobreviveu à correção. Para quimioterapia isso é em boa parte uma questão de poder: apenas 396 pacientes tratadas, com 212 eventos. Vale registrar como sinal exploratório, sem valor confirmatório, que os menores FDR apontaram *AKT1S1*, *SLC19A1*, *SF3B1* e *ERBB2* — nomes plausíveis (*SLC19A1* é o transportador de folato reduzido, via de entrada do metotrexato, componente do esquema CMF usado nessa coorte histórica).

### 5.13 Integração das camadas — genes prioritários por subtipo

Cruzando marcadores de subtipo com genes prognósticos, 223 pares gene–subtipo são simultaneamente as duas coisas. Escore de prioridade = |Δz| × |log HR| × −log₁₀(FDR).

| Subtipo | Genes prioritários | Leitura |
|---|---|---|
| **Basal** | *AURKA* (↑, HR 1,19), *FANCD2* (↑, 1,16), *VEGFA* (↑, 1,16), *IGF1* (↓, 0,80), *LAMA2* (↓), *BCL2* (↓) | Coerentes com **pior** prognóstico: superexpressa genes de risco e perde genes protetores |
| **Her2** | *GSK3B* (↑, 1,22), *AURKA* (↑), *SPRY2* (↓), *BCL2* (↓), *STAT5B* (↓), *KIT* (↓) | Pior prognóstico, dominado pela superexpressão de *GSK3B* |
| **LumB** | *SPRY2* (↓), *ABCB1* (↓), *STAT5A* (↓), *PDGFRA* (↓), *CCND2* (↓) | Pior prognóstico por **perda de genes protetores**, não por ganho de genes de risco |
| **LumA** | *AURKA* (↓), *BCL2* (↑), *FANCD2* (↓), *E2F7* (↓), *E2F2* (↓) | Perfil de **melhor** prognóstico — o espelho exato do Basal |
| **claudin-low** | *STAT5A* (↑), *ABCB1* (↑), *IGF1* (↑), *CDKN2C* (↑), *CCND2* (↑) | Melhor prognóstico, consistente com a maior sobrevida mediana |
| **Normal** | *SPRY2* (↑), *ABCB1* (↑), *LAMA2* (↑), *AURKA* (↓), *KIT* (↑) | Melhor prognóstico |

A simetria entre Basal e LumA — os mesmos genes (*AURKA*, *FANCD2*, *E2F2*, *BCL2*) com direções invertidas — indica que o eixo prognóstico dominante na coorte não é um conjunto de genes independentes, mas um **gradiente único de proliferação**. Isso também explica por que a assinatura LASSO acrescenta tão pouco a um modelo que já contém subtipo e grau: ela está medindo, por um caminho mais longo e mais instável, algo que o grau histológico já mede.

---

## 6. Respostas diretas às perguntas da pesquisa

**Quais genes definem cada subtipo molecular?**
*ERBB2* (Her2); *CCNE1*, *CDKN2A*, *CHEK1*, *E2F3*, *CDC25A* (Basal); *GATA3*, *MAPT*, *BCL2* (LumA); perda de *EGFR*, *NOTCH1*, *FOXO1*, *SPRY2* (LumB); *FOLR2*, *KLRG1*, *CSF1R*, *TGFBR2* com perda de *ERBB3* e *RAB25* (claudin-low). No plano mutacional: *TP53* (Basal, Her2) contra *PIK3CA*, *GATA3*, *MAP3K1*, *CBFB* (LumA).

**Quais genes se associam a MENOR sobrevida?**
*GSK3B* (HR 1,23 bruto; 1,13 ajustado; selecionado em 72% das repetições do LASSO), *AURKA* (1,19 em OS; 1,46 em mortalidade específica; primeiro no ranking do Random Forest de desfecho), *VEGFA* (1,17; 1,13 ajustado), *FANCD2*, *E2F2*, *E2F7*, *CCNE1*. Por estágio, no estágio II acrescentam-se *ACVR1B* (1,28), *MEN1* (1,24), *MAML1* (1,24) e *AR* (1,22).

**Quais genes se associam a MAIOR sobrevida?**
*STAT5A* (HR 0,81; segundo no ranking do Random Forest; selecionado em 60% das repetições), *SPRY2* (0,81), *IGF1* (0,80), *ABCB1* (0,82), *FLT3* (0,82), *CDKN2C* (0,82), *STAT5B* (0,84), *PDGFRA* (0,84), *BCL2* (0,85), *CCND2* (0,85), *LAMA2* (0,84). No estágio I, *NRIP1* (0,73) é o único a atingir significância.

**Em que contexto esses genes funcionam como marcadores?**
Nos tumores **luminais** (12 de 12 genes significativos no LumA; AUC 0,666 no Random Forest) e no **estágio II** (71 genes significativos). Nos subtipos Her2 e Basal, e nos estágios I e III–IV, praticamente nenhum gene sobrevive à correção para múltiplos testes. O contexto molecular e clínico não é um detalhe do resultado — é parte do resultado.

**Existe sinal farmacogenômico?**
Sim, e apenas para hormonioterapia: 14 genes modificam significativamente o efeito prognóstico conforme a paciente tenha recebido ou não o bloqueio endócrino, e todos resistem ao ajuste por status de ER e subtipo. O grupo *MTOR* / *SETDB1* / *MAP3K13* / *BAP1* é o candidato mais interessante a marcador de resistência endócrina. Para quimioterapia e radioterapia, nenhuma interação sobreviveu à correção — para quimioterapia, sobretudo por limitação de poder (396 tratadas).

**Com que grau de confiança?**
**Alto** para a existência das associações: n = 1.897, FDR rigoroso, efeitos replicados entre OS e DSS, entre R e Python, e entre três famílias metodológicas (modelos lineares, Cox, Random Forest). **Baixo** para utilidade preditiva incremental: nenhuma combinação de genes superou, fora da amostra, um modelo clínico simples — nem por LASSO-Cox (0,555 contra 0,655) nem por Random Forest (0,653 contra 0,730). **Intermediário e exploratório** para a camada farmacogenômica: os achados são estatisticamente sólidos e resistentes a confundimento medido, mas derivam de uma coorte observacional histórica, sem randomização de tratamento e sem validação externa.

---

## 7. Limitações

1. **Cobertura gênica parcial.** 489 genes de um painel curado, não o transcriptoma completo (~24.000 genes). Genes prognósticos fora do painel não podem ser detectados, e análises de enriquecimento de vias ficariam enviesadas pela composição do painel.
2. **Riscos proporcionais violados** em 12 dos 20 genes mais fortes; os HR são efeitos médios ao longo de quase 10 anos de seguimento.
3. **Ausência de validação externa.** Toda a validação foi interna (25 partições, validação cruzada 5-fold). Confirmação em TCGA-BRCA ou coorte própria é o passo obrigatório seguinte.
4. **Coorte histórica.** METABRIC é majoritariamente pré-trastuzumabe, com esquemas adjuvantes heterogêneos e sem registro de duração ou adesão à hormonioterapia (apenas sim/não). Os achados farmacogenômicos refletem a prática de uma era.
5. **Tratamento não randomizado.** Interação gene × tratamento em dados observacionais é vulnerável a confundimento por indicação. O ajuste por ER e subtipo mitiga a explicação mais óbvia, mas não substitui um ensaio.
6. **Estágio ausente em 26,2% das pacientes** e apenas 124 casos em III–IV, o que limita fortemente o poder nesse estrato.
7. **Mortalidade competitiva.** A OS inclui 480 óbitos por outras causas em coorte de mediana etária elevada; a análise de DSS mitiga, mas o tratamento formal exigiria modelos de risco competitivo (Fine-Gray).
8. **Grupo "Normal-like"** provavelmente reflete contaminação por tecido normal (80,7% de erro de classificação); resultados nesse estrato pedem cautela.
9. **z-scores calculados sobre a coorte inteira** implicam que o valor de cada amostra depende da composição do estudo, limitando a aplicação direta do escore a uma paciente individual sem recalibração.
10. **Ausência de dados de genótipo germinativo.** Farmacogenômica clássica estuda polimorfismos em genes de metabolização (*CYP2D6* para tamoxifeno, *DPYD*, *UGT1A1*). O que este trabalho analisa é **expressão tumoral como modificadora de efeito do tratamento** — uma farmacogenômica somática, que deve ser descrita como tal.

---

## 8. Próximos passos sugeridos, em ordem de retorno esperado

1. **Restringir o escopo aos tumores luminais em estágio II** e rederivar a assinatura nesse contexto. É onde todo o sinal se concentra e a hipótese mais promissora para superar o modelo clínico.
2. **Validar externamente** os 14 genes de interação com hormonioterapia em TCGA-BRCA, com atenção à diferença de plataforma (microarranjo vs. RNA-seq).
3. **Aprofundar o eixo mTOR/*SETDB1*** como marcador de resistência endócrina, incluindo o cruzamento com mutações de *PIK3CA* já disponíveis nesta mesma coorte.
4. **Obter a matriz de expressão completa** (~24.000 genes) do cBioPortal e repetir os módulos de expressão diferencial, Cox e assinatura.
5. **Modelagem tempo-dependente** para os 12 genes que violam riscos proporcionais, ou restrição da janela de análise (OS em 5 anos).
6. **Modelos de risco competitivo** (Fine-Gray) dada a alta proporção de óbitos por outras causas.
7. **Random Survival Forest e gradient boosting de sobrevida**, que usam o tempo até o evento diretamente em vez de dicotomizar em 10 anos.

---

## 9. Bibliotecas utilizadas e finalidade de cada uma

### R (pipeline principal)

| Pacote | Versão | Finalidade nesta análise |
|---|---|---|
| data.table | 1.14.10 | Leitura de alto desempenho do CSV de 8,1 MB e escrita das 33 tabelas |
| dplyr | 1.1.4 | Filtros, agrupamentos e sumarizações |
| tidyr | 1.3.1 | Conversões wide↔long para heatmaps e boxplots |
| purrr | 1.0.2 | Iteração funcional sobre os 489 genes nos modelos de Cox |
| stringr | 1.5.1 | Normalização de nomes de genes e rótulos |
| forcats | 1.0.0 | Ordenação de fatores em gráficos |
| **limma** | 3.58.1 | Núcleo da expressão diferencial: modelos lineares por gene com moderação bayesiana empírica da variância |
| **survival** | 3.5-8 | Objetos `Surv`, Kaplan-Meier, log-rank, modelos de Cox, teste `cox.zph`, concordância |
| survminer | 0.4.9 | Curvas de sobrevida com tabela de números sob risco |
| **glmnet** | 4.1-8 | Cox penalizado por LASSO com validação cruzada |
| **randomForest** | 4.7-1.1 | Classificação multiclasse dos subtipos e importância dos genes |
| matrixStats | 1.2.0 | Estatísticas por linha da matriz de expressão no controle de qualidade |
| broom | 1.0.5 | Padronização de saídas de modelos |
| ggplot2 | 3.4.4 | Motor gráfico das 13 figuras |
| ggrepel | 0.9.5 | Rótulos sem sobreposição nos volcano plots |
| scales | 1.3.0 | Formatação de eixos |
| RColorBrewer | 1.1-3 | Paleta divergente dos heatmaps |
| viridis | 0.6.5 | Escalas contínuas perceptualmente uniformes |
| pheatmap | 1.0.12 | Heatmap com clusterização hierárquica |
| patchwork | 1.2.0 | Composição de painéis |
| cowplot | 1.1.3 | Tema gráfico consistente |
| gridExtra | 2.3 | Manipulação de objetos grid na exportação |

### Python (análises complementares)

| Pacote | Versão | Finalidade nesta análise |
|---|---|---|
| **lifelines** | 0.30.3 | Modelos de Cox por estágio, termos de interação gene × tratamento, Kaplan-Meier e log-rank por estágio |
| **scikit-learn** | 1.8.0 | Random Forest (subtipo e mortalidade em 10 anos), validação cruzada estratificada, AUC |
| statsmodels | — | Correção de Benjamini-Hochberg para múltiplos testes |
| pandas / numpy | 3.0.2 / — | Manipulação da matriz de dados e operações numéricas |

---

## 10. Inventário completo dos arquivos

### 10.1 Pipeline em R — 33 tabelas (`saidas/tabelas/`)

| Arquivo | Conteúdo |
|---|---|
| T01–T04 | Integridade dos dados, auditoria da codificação do evento, fluxograma de exclusões, QC dos 489 genes |
| T05–T06 | Caracterização clínica por subtipo e testes de associação |
| T07 | **2.934 linhas** — expressão diferencial completa (489 genes × 6 subtipos) |
| T08–T09 | Resumo da expressão diferencial e top 15 genes por subtipo |
| T10–T11 | Matriz de confusão e ranking de importância do Random Forest |
| T12–T13 | Cox univariado dos 489 genes para OS e para DSS |
| T14–T15 | Teste de riscos proporcionais e Cox ajustado por 8 covariáveis |
| T16–T17 | Medianas de sobrevida por subtipo e log-rank dos 4 genes principais |
| T18 | Cox estratificado: 12 genes × 6 subtipos |
| T19, T19b, T20–T22 | Assinatura LASSO (65 e 12 genes), C-index dos 5 modelos, razão de verossimilhança, tercis de risco |
| T23–T24 | Frequência mutacional por subtipo e Cox das mutações |
| T25–T26 | **223 pares** gene–subtipo integrados e genes prioritários por subtipo |
| T27–T28 | Resumo da execução e versões das 22 bibliotecas |
| T29–T32 | Validação repetida em 25 partições, resumo, comparações pareadas, estabilidade de seleção |

### 10.2 Pipeline em R — 13 figuras (`saidas/figuras/`)

F01 distribuição dos subtipos · F02 volcano por subtipo · F03 heatmap dos genes discriminantes · F04 importância no Random Forest · F05 forest plot do Cox univariado · F06 Kaplan-Meier por subtipo · F07 Kaplan-Meier dos 4 genes principais · F08 heatmap do Cox estratificado · F09 grupos de risco na validação · F10 coeficientes do LASSO · F11 mutações por subtipo · F12 integração subtipo × prognóstico · F13 estabilidade e validação repetida

### 10.3 Análises complementares (`saidas_complementares/`)

| Arquivo | Conteúdo |
|---|---|
| C01_cox_por_estagio_todos_genes.csv | 489 genes × 3 estágios: HR, IC 95%, p, FDR |
| C02_interacao_gene_tratamento.csv | 1.467 modelos de interação gene × tratamento |
| C03_cox_por_estrato_tratamento.csv | HR dos genes-núcleo dentro de cada estrato de tratamento |
| C04_rf_subtipo_importancia.csv | Importância dos 489 genes na classificação dos subtipos (Python) |
| C05_rf_mortalidade10a_importancia.csv | Importância dos 489 genes na predição de óbito em 10 anos |
| C06_rf_auc_mortalidade10a.csv | AUC dos modelos de genes, clínico e combinado |
| C07_rf_auc_por_estrato.csv | AUC por subtipo e por estágio |
| C08_interacao_hormonio_ajustada_ER_subtipo.csv | Reanálise dos 14 genes com ajuste por ER, subtipo e quimioterapia |

### 10.4 Código e reprodutibilidade

- `metabric_pipeline.R` — script único, autossuficiente: detecta o próprio diretório, instala dependências ausentes, baixa os dados, verifica o MD5 e executa os 12 módulos. Execução: `Rscript metabric_pipeline.R` (~6,6 min).
- `analise_complementar.py`, `parte_b_farmaco.py`, `parte_c_rf.py` — análises de estágio, farmacogenômica e aprendizado de máquina.
- `saidas/objetos_analise.rds` — todos os objetos da análise, para reanálise sem reexecutar.
- `saidas/logs/execucao.log` e `sessionInfo.txt` — ambiente exato para reprodução.

Todas as etapas aleatórias usam semente fixa (42), de modo que os números deste relatório são reproduzidos integralmente.

---

## 11. Referências

- Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*. 2012;486(7403):346-352.
- Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. *Nat Commun*. 2016;7:11479.
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47.
- Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. *J Stat Softw*. 2011;39(5):1-13.
- Breiman L. Random Forests. *Mach Learn*. 2001;45(1):5-32.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*. 1995;57(1):289-300.
- Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. Springer; 2000.
- Meinshausen N, Bühlmann P. Stability selection. *J R Stat Soc B*. 2010;72(4):417-473.
