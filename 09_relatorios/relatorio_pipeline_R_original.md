# Genes associados a subtipos moleculares e a pior prognóstico no câncer de mama

## Relatório final da pipeline analítica em R — coorte METABRIC (n = 1.897)

**Data da execução:** 06/08/2026 · **Ambiente:** R 4.3.3 (2024-02-29), Ubuntu 24.04 · **Tempo total de execução:** 6,57 min · **Semente global:** 42

---

## 1. Sumário executivo

Foi construída e executada uma pipeline reprodutível em R, de ponta a ponta, sobre **dados reais e não simulados** da coorte METABRIC (Molecular Taxonomy of Breast Cancer International Consortium). A análise integrou expressão gênica (489 genes, z-scores de microarranjo), status mutacional (173 genes) e desfechos de sobrevida de 1.897 pacientes com seguimento mediano de 115,6 meses (1.098 óbitos por qualquer causa; 619 óbitos atribuídos ao câncer de mama).

Achados centrais:

1. **Os subtipos moleculares diferem fortemente em prognóstico** (log-rank global p = 3,7 × 10⁻¹⁰). A sobrevida global mediana variou de 219,2 meses (claudin-low) a 104,0 meses (Her2).
2. **568 associações gene–subtipo** atingiram significância (FDR < 0,01 e |Δz| ≥ 0,50), com assinaturas biologicamente coerentes: *ERBB2* no Her2; *CCNE1*, *CDKN2A*, *CHEK1*, *E2F3* no Basal; *GATA3*, *MAPT*, *BCL2* no LumA; *FOLR2*, *KLRG1*, *CSF1R*, *CSF1R*/*TGFBR2* (assinatura imune-estromal) com perda de *ERBB3*/*RAB25* no claudin-low.
3. **125 genes** associaram-se independentemente à sobrevida global e **137** à sobrevida específica por câncer (FDR < 0,05). Os efeitos mais robustos, mantidos após ajuste por idade, grau, tamanho, linfonodos, subtipo e tratamentos, foram *STAT5A* (HR ajustado 0,867), *GSK3B* (1,133), *VEGFA* (1,125), *BCL2* (0,885) e *AURKA* (1,107).
4. **A assinatura multigênica derivada por LASSO-Cox não superou o modelo clínico**: em 25 partições independentes, C-index médio de 0,555 (assinatura) versus 0,655 (clínico), e o modelo combinado não acrescentou discriminação sobre o clínico (diferença média +0,0003; Wilcoxon pareado p = 0,70). Esse resultado negativo é relatado explicitamente porque contradiz a leitura ingênua do teste de razão de verossimilhança feito na amostra de derivação (χ² = 161,5; p = 5,3 × 10⁻³⁷) — a discrepância é a assinatura digital do sobreajuste.
5. **Mutações somáticas segregam nitidamente por subtipo**: *TP53* em 88,4% dos Basal e 70,0% dos Her2 contra 11,9% dos LumA; *PIK3CA* em 57,4% dos LumA contra 16,1% dos Basal (ambos FDR = 5,1 × 10⁻⁴). Mutação em *GATA3* associou-se a melhor sobrevida (HR 0,590; FDR 3,3 × 10⁻⁵) e em *TP53* a pior (HR 1,289; FDR 1,3 × 10⁻³).

---

## 2. Dados: origem, proveniência e integridade

| Item | Valor |
|---|---|
| Estudo | METABRIC / cBioPortal `brca_metabric` |
| Publicações originais | Curtis et al., *Nature* 486:346-352 (2012); Pereira et al., *Nat Commun* 7:11479 (2016) |
| Arquivo processado | `METABRIC_RNA_Mutation.csv` (recorte curado do cBioPortal) |
| Checksum MD5 verificado na pipeline | `c619471beb87af0f9fd4e4a40058654d` |
| Pacientes no arquivo bruto | 1.904 |
| Colunas | 693 (31 clínicas + 489 expressão + 173 mutação) |
| IDs duplicados | 0 |
| Valores ausentes na matriz de expressão | 0 |
| Faixa dos z-scores | −7,2432 a +20,3950 |
| Plataforma de expressão | Illumina HT-12 v3 (microarranjo), z-scores calculados pelo cBioPortal |
| Plataforma mutacional | sequenciamento alvo de painel de genes drivers |

**Natureza dos dados.** Nenhum valor foi simulado, imputado por modelo ou gerado sinteticamente. O único tratamento aplicado a valores ausentes foi a exclusão por caso (*complete-case*) nos modelos que exigiam covariáveis clínicas.

### 2.1 Auditoria da codificação do desfecho (etapa crítica)

Esse conjunto de dados é frequentemente analisado com o evento invertido. A pipeline audita a codificação antes de qualquer modelagem (tabela `T02`):

| `overall_survival` | Died of Disease | Died of Other Causes | Living |
|---|---|---|---|
| **0** | 622 | 480 | 0 |
| **1** | 0 | 0 | 801 |

Ou seja, `overall_survival = 1` significa **paciente vivo**. Definições adotadas:

- **OS (sobrevida global):** evento = `1 − overall_survival` → 1.098 óbitos.
- **DSS (sobrevida específica do câncer):** evento = 1 se `death_from_cancer == "Died of Disease"`; óbitos por outras causas são **censurados** (não tratados como evento) → 619 eventos. Um paciente sem informação de causa permaneceu como `NA` e foi excluído apenas das análises de DSS.

### 2.2 Fluxograma de exclusões

| Etapa | Excluídos | Restantes |
|---|---|---|
| Coorte inicial | — | 1.904 |
| Histologia não mamária (Breast Sarcoma) | 1 | 1.903 |
| Subtipo PAM50/claudin-low = "NC" (não classificado) | 6 | 1.897 |
| Tempo de seguimento ausente ou ≤ 0 | 0 | 1.897 |
| **Coorte analítica final** | — | **1.897** |

---

## 3. Ambiente computacional e bibliotecas

O script instala automaticamente qualquer biblioteca ausente (CRAN + Bioconductor) antes de iniciar a análise — ver seção 9. As versões efetivamente utilizadas nesta execução ficam registradas em `T28_bibliotecas_versoes.csv` e em `saidas/logs/sessionInfo.txt`.

| Pacote | Versão | Finalidade específica nesta pipeline |
|---|---|---|
| **data.table** | 1.14.10 | Leitura de alto desempenho (`fread`) do CSV de 8,1 MB e escrita de todas as 33 tabelas de saída (`fwrite`) |
| **dplyr** | 1.1.4 | Manipulação tabular: filtros, agrupamentos e sumarizações da caracterização clínica e das tabelas de integração |
| **tidyr** | 1.3.1 | Conversão wide↔long (`pivot_longer`) para os heatmaps de mutação e para o boxplot de C-index |
| **purrr** | 1.0.2 | Iteração funcional (`map_dfr`) sobre os 489 genes nos modelos de Cox e sobre os contrastes de subtipo |
| **stringr** | 1.5.1 | Normalização de nomes de genes e rótulos de figuras |
| **forcats** | 1.0.0 | Ordenação de fatores em gráficos (`fct_reorder`) e definição da ordem dos subtipos |
| **limma** | 3.58.1 | Núcleo da expressão diferencial: modelos lineares por gene com moderação bayesiana empírica da variância (`lmFit`, `contrasts.fit`, `eBayes`) |
| **survival** | 3.5-8 | Objetos `Surv`, estimador de Kaplan-Meier, teste log-rank (`survdiff`), modelos de Cox (`coxph`), teste de riscos proporcionais (`cox.zph`) e cálculo de concordância (`concordance`) |
| **survminer** | 0.4.9 | Curvas de sobrevida com tabela de números sob risco e valor-p do log-rank (`ggsurvplot`) |
| **glmnet** | 4.1-8 | Regressão de Cox penalizada por LASSO com validação cruzada (`cv.glmnet`, `family = "cox"`) para derivação da assinatura multigênica |
| **randomForest** | 4.7-1.1 | Classificação multiclasse dos subtipos e ranqueamento multivariado de importância dos genes (Gini e queda de acurácia) |
| **matrixStats** | 1.2.0 | Estatísticas por linha da matriz de expressão (desvio-padrão, MAD, mínimo, máximo) no controle de qualidade |
| **broom** | 1.0.5 | Padronização de saídas de modelos em `data.frame` |
| **ggplot2** | 3.4.4 | Motor gráfico de todas as 13 figuras |
| **ggrepel** | 0.9.5 | Rótulos de genes sem sobreposição nos volcano plots e no gráfico de integração |
| **scales** | 1.3.0 | Formatação de eixos (escala log₂ do forest plot) |
| **RColorBrewer** | 1.1-3 | Paleta divergente RdBu do heatmap e paletas categóricas dos subtipos |
| **viridis** | 0.6.5 | Escalas contínuas perceptualmente uniformes (importância do RF, frequências de mutação) |
| **pheatmap** | 1.0.12 | Heatmap com clusterização hierárquica dos genes discriminantes |
| **patchwork** | 1.2.0 | Composição de painéis múltiplos em uma única figura |
| **cowplot** | 1.1.3 | Tema gráfico consistente e alinhamento de painéis |
| **gridExtra** | 2.3 | Manipulação de objetos grid na exportação das figuras |

---

## 4. Metodologia

A pipeline está organizada em 12 módulos sequenciais no arquivo `metabric_pipeline.R` (764 linhas, integralmente comentado em português). Parâmetros analíticos foram declarados *a priori* em um único objeto `PAR`, evitando escolhas post-hoc: FDR < 0,01 e |Δz| ≥ 0,50 para expressão diferencial; FDR < 0,05 para Cox; 70/30 para partição; 500 árvores no Random Forest; 10 dobras na validação cruzada; frequência mínima de 2% para testar mutações.

**Módulo 0 — Ambiente.** Carregamento de bibliotecas, semente fixa (42), definição de diretórios e funções auxiliares de exportação.

**Módulo 1 — Ingestão e integridade.** Leitura do CSV, identificação automática dos blocos de colunas (clínicas, expressão, `*_mut`), verificação de duplicatas, ausentes, faixa de valores e checksum MD5.

**Módulo 2 — QC e pré-processamento.** Auditoria da codificação do desfecho (seção 2.1), aplicação dos critérios de exclusão, padronização de covariáveis clínicas, montagem da matriz genes × amostras (489 × 1.897) e binarização das colunas mutacionais (`0` = selvagem; qualquer variante = mutado).

**Módulo 3 — Caracterização clínica.** Estatísticas descritivas por subtipo; associação com variáveis categóricas por qui-quadrado de Pearson e com variáveis contínuas por Kruskal-Wallis (não paramétrico, sem assumir normalidade), com correção de Benjamini-Hochberg (BH).

**Módulo 4 — Expressão diferencial (limma).** Modelo linear sem intercepto (`~ 0 + subtipo`) com contrastes **um-vs-demais** para cada um dos 6 subtipos, isto é, cada subtipo contra a média dos outros cinco. Moderação bayesiana empírica com `trend = TRUE` e `robust = TRUE`; ajuste BH. Como a expressão já está em z-score, o coeficiente do contraste (Δz) é diretamente interpretável como diferença em desvios-padrão — um tamanho de efeito, não apenas um p-valor.

**Módulo 5 — Random Forest.** Classificação dos 6 subtipos a partir dos 489 genes (500 árvores), com erro out-of-bag como estimativa honesta de generalização e ranqueamento de importância (queda média de Gini e de acurácia). Serve como verificação multivariada e não-linear do que o limma detecta gene a gene.

**Módulo 6 — Cox univariado gene a gene.** 489 modelos de Cox para OS e 489 para DSS, com a expressão em z-score, de modo que o HR expressa o risco por **aumento de 1 desvio-padrão** de expressão. Correção BH. Para os 20 genes mais fortes: teste formal do pressuposto de riscos proporcionais (`cox.zph`) e reanálise em modelo multivariável ajustado por idade, grau, tamanho tumoral, status linfonodal, subtipo, quimioterapia, hormonioterapia e radioterapia.

**Módulo 7 — Kaplan-Meier e análise estratificada.** Curvas por subtipo e por expressão dicotomizada na mediana (com a ressalva metodológica de que dicotomizar variáveis contínuas reduz poder — usada aqui apenas para visualização, enquanto a inferência permanece nos modelos contínuos). Cox dentro de cada subtipo para os 12 genes mais prognósticos, testando se o efeito é global ou dependente do contexto molecular.

**Módulo 8 — Assinatura multigênica (LASSO-Cox).** Partição estratificada 70/30 preservando a distribuição conjunta de subtipo e evento; `cv.glmnet` com 10 dobras; extração dos coeficientes em `lambda.min` e `lambda.1se`; escore de risco linear; comparação de C-index entre modelo gênico, clínico e combinado; teste de razão de verossimilhança; estratificação em tercis de risco definidos **no conjunto de derivação** e aplicados ao de validação.

**Módulo 8b — Estabilidade e validação repetida.** Como uma única partição 70/30 produz estimativas instáveis, todo o procedimento (partição → `cv.glmnet` → avaliação) foi repetido **25 vezes** com sementes distintas, registrando a distribuição do C-index dos três modelos, comparações pareadas por Wilcoxon e a frequência de seleção de cada gene (*stability selection*).

**Módulo 9 — Mutações somáticas.** Dos 173 genes, 71 atingiram frequência ≥ 2% e foram testados: teste exato de Fisher (com p simulado por Monte Carlo, B = 20.000, para tabelas 2×6) para heterogeneidade entre subtipos, e Cox univariado para associação com OS; ambos com correção BH.

**Módulo 10 — Integração.** Cruzamento das camadas: genes simultaneamente marcadores de subtipo (FDR < 0,01, |Δz| ≥ 0,50) e prognósticos (FDR < 0,05), classificados quanto à **coerência direcional** (se o gene está superexpresso no subtipo e sua alta expressão prediz maior risco, é marcador de pior prognóstico naquele subtipo). Escore de prioridade = |Δz| × |log HR| × −log₁₀(FDR_OS), com anotação da posição no ranking do Random Forest e do coeficiente LASSO.

---

## 5. Resultados

### 5.1 Caracterização clínica por subtipo (T05, T06, F01)

| Subtipo | n (%) | Idade mediana | Tamanho (mm) | G3 (%) | N+ (%) | ER+ (%) | HER2+ (%) | NPI | Óbitos (%) | Óbitos por câncer |
|---|---|---|---|---|---|---|---|---|---|---|
| LumA | 679 (35,8) | 63,2 | 20,5 | 25,9 | 42,3 | 99,6 | 3,1 | 3,08 | 53,6 | 146 |
| LumB | 461 (24,3) | 66,2 | 25,0 | 55,6 | 51,4 | 100,0 | 9,1 | 4,05 | 65,7 | 181 |
| Her2 | 220 (11,6) | 58,8 | 25,0 | 73,7 | 57,3 | 42,3 | 56,8 | 4,08 | 70,5 | 107 |
| Basal | 199 (10,5) | 54,5 | 25,0 | 90,5 | 53,8 | 13,6 | 10,1 | 5,02 | 55,8 | 81 |
| claudin-low | 198 (10,4) | 58,5 | 20,0 | 67,9 | 47,5 | 38,4 | 7,6 | 4,05 | 44,9 | 56 |
| Normal | 140 (7,4) | 57,8 | 23,0 | 34,1 | 40,0 | 85,7 | 9,3 | 4,03 | 54,3 | 48 |

Todas as 14 variáveis clínicas testadas diferiram significativamente entre subtipos após correção BH (T06). As associações mais fortes: status de ER (χ² = 1.093,9; FDR = 4,0 × 10⁻²³³), grau histológico (χ² = 403,8; FDR = 5,1 × 10⁻⁸⁰), PR (FDR = 3,9 × 10⁻¹⁰⁶), HER2 (FDR = 2,8 × 10⁻⁹⁷) e NPI (Kruskal-Wallis, FDR = 8,2 × 10⁻⁴⁶). O subtipo Basal concentra pacientes mais jovens (mediana 54,5 anos), tumores G3 (90,5%) e o maior NPI (5,02).

### 5.2 Sobrevida global por subtipo (T16, F06)

| Subtipo | n | Óbitos | Mediana de OS (meses) | IC 95% |
|---|---|---|---|---|
| claudin-low | 198 | 89 | 219,2 | 194,1 – 238,1 |
| LumA | 679 | 364 | 186,6 | 169,0 – 198,1 |
| Normal | 140 | 76 | 158,5 | 125,8 – 203,5 |
| Basal | 199 | 111 | 130,9 | 83,4 – 206,6 |
| LumB | 461 | 303 | 123,0 | 114,9 – 143,1 |
| Her2 | 220 | 155 | 104,0 | 88,9 – 142,4 |

Teste log-rank global: χ² correspondente a **p = 3,74 × 10⁻¹⁰**.

Observação relevante: no METABRIC, coorte majoritariamente pré-trastuzumabe, o subtipo **Her2 apresenta a pior sobrevida mediana**, e o Basal — apesar do perfil clínico mais agressivo — não se destaca como o pior em sobrevida global de longo prazo, refletindo o padrão bem descrito de risco precoce elevado seguido de platô, enquanto os luminais mantêm risco tardio persistente. Isso também explica por que o pressuposto de riscos proporcionais é problemático nesta coorte (seção 5.5).

### 5.3 Expressão diferencial por subtipo (T07–T09, F02, F03)

**568 associações gene–subtipo** significativas (FDR < 0,01 **e** |Δz| ≥ 0,50):

| Subtipo | Significativos por FDR | Com efeito ≥ 0,50 | Superexpressos | Subexpressos |
|---|---|---|---|---|
| Basal | 300 | 151 | 89 | 62 |
| claudin-low | 269 | 107 | 62 | 45 |
| LumB | 308 | 92 | 33 | 59 |
| Normal | 201 | 78 | 38 | 40 |
| LumA | 287 | 77 | 30 | 47 |
| Her2 | 248 | 63 | 36 | 27 |

Principais marcadores por subtipo (Δz; todos com FDR ≪ 0,01):

- **Basal** — *CCNE1* (+1,67), *CDKN2A* (+1,66), *CHEK1* (+1,57), *MAP2* (+1,53), *E2F3* (+1,51), *TTYH1* (+1,49), *CDC25A* (+1,45); *TGFB3* (−1,41). Assinatura clássica de desregulação do ciclo celular G1/S e resposta a dano de DNA.
- **Her2** — *ERBB2* (+1,57), *ARRDC1* (+1,12), *GSK3B* (+1,04), *MMP15* (+1,02), *AKT1* (+0,99); *SMAD4* (−1,11), *MYC* (−0,98), *BCL2* (−0,96). A superexpressão de *ERBB2* funciona como controle positivo interno da pipeline.
- **LumA** — *GATA3* (+1,13), *MAPT* (+1,12), *APH1B* (+0,94), *BCL2* (+0,91); *AURKA* (−0,96), *CCNE1* (−0,94), *E2F2* (−0,93), *CHEK1* (−0,92). Perfil de baixa proliferação com diferenciação luminal preservada.
- **LumB** — *GATA3* (+1,00); *EGFR* (−1,06), *LAMB3* (−0,94), *NOTCH1* (−0,94), *PLAGL1* (−0,92), *FOXO1* (−0,90), *SPRY2* (−0,90), *MMP7* (−0,88). Identidade luminal com perda de reguladores negativos de crescimento.
- **claudin-low** — *FOLR2* (+1,57), *KLRG1* (+1,56), *CSF1R* (+1,37), *TGFBR2* (+1,27), *CCND2* (+1,26), *ACVRL1* (+1,19); *ERBB3* (−1,51), *RAB25* (−1,42). Assinatura imune/estromal (macrófagos, linfócitos) com perda de identidade epitelial — exatamente o fenótipo mesenquimal que define esse grupo.
- **Normal** — *NR2F1* (+1,05), *ABCB1* (+1,02), *SPRY2* (+0,98), *LAMA2* (+0,98); *CDK1* (−1,04), *CHEK1* (−0,97), *AURKA* (−0,94).

### 5.4 Classificação por Random Forest (T10, T11, F04)

Erro out-of-bag global: **23,77%** (acurácia ≈ 76%) na discriminação de 6 classes a partir de 489 genes. Erro por classe:

| Subtipo | Erro de classificação |
|---|---|
| LumA | 11,3% |
| LumB | 18,0% |
| Basal | 20,6% |
| claudin-low | 29,8% |
| Her2 | 35,5% |
| **Normal** | **80,7%** |

O erro de 80,7% no grupo "Normal-like" é um resultado informativo, não um defeito do modelo: 93 dos 140 casos são classificados como LumA. Isso é consistente com a literatura que questiona se o Normal-like é uma entidade biológica ou um artefato de baixa celularidade tumoral com contaminação por tecido mamário normal. A confusão residual Her2 → LumA/LumB (76 casos) reflete a sobreposição entre HER2 clínico e subtipo intrínseco.

### 5.5 Genes associados à sobrevida (T12–T15, F05)

**125 genes** com FDR < 0,05 para OS e **137** para DSS. Top 15 para sobrevida global (HR por +1 DP de expressão):

| Gene | HR | IC 95% | FDR |
|---|---|---|---|
| GSK3B | 1,225 | 1,16–1,29 | 3,8 × 10⁻¹⁰ |
| STAT5A | 0,807 | 0,76–0,86 | 1,7 × 10⁻⁹ |
| SPRY2 | 0,810 | 0,76–0,86 | 2,6 × 10⁻⁸ |
| IGF1 | 0,804 | 0,75–0,86 | 1,0 × 10⁻⁷ |
| ABCB1 | 0,822 | 0,77–0,88 | 2,6 × 10⁻⁷ |
| AURKA | 1,191 | 1,12–1,26 | 2,6 × 10⁻⁷ |
| LAMA2 | 0,836 | 0,79–0,89 | 2,6 × 10⁻⁷ |
| FLT3 | 0,820 | 0,77–0,88 | 2,6 × 10⁻⁷ |
| STAT5B | 0,838 | 0,79–0,89 | 4,9 × 10⁻⁷ |
| CDKN2C | 0,815 | 0,76–0,87 | 5,0 × 10⁻⁷ |
| PDGFRA | 0,844 | 0,80–0,89 | 6,2 × 10⁻⁷ |
| RPS6 | 0,855 | 0,81–0,90 | 1,6 × 10⁻⁶ |
| CCND2 | 0,846 | 0,80–0,90 | 2,3 × 10⁻⁶ |
| BCL2 | 0,849 | 0,80–0,90 | 5,0 × 10⁻⁶ |
| VEGFA | 1,165 | 1,10–1,23 | 5,8 × 10⁻⁶ |

Para o desfecho específico do câncer (DSS), os efeitos são **maiores em magnitude**, como esperado ao remover o ruído da mortalidade competitiva: *AURKA* HR 1,458 (FDR 1,5 × 10⁻²⁰), *BCL2* 0,723, *MAPT* 0,735, *GSK3B* 1,326, *FANCD2* 1,323, *E2F2* 1,310, *CCNE1* 1,271.

**Ajuste por covariáveis clínicas** (idade, grau, tamanho, linfonodos, subtipo, quimio, hormônio, radioterapia). Dos 20 genes testados, 10 mantiveram significância após correção BH:

| Gene | HR ajustado | IC 95% | FDR |
|---|---|---|---|
| STAT5A | 0,867 | 0,81–0,93 | 5,8 × 10⁻⁴ |
| GSK3B | 1,133 | 1,06–1,21 | 2,5 × 10⁻³ |
| VEGFA | 1,125 | 1,05–1,20 | 2,5 × 10⁻³ |
| BCL2 | 0,885 | 0,82–0,96 | 0,015 |
| CDKN2C | 0,896 | 0,83–0,97 | 0,018 |
| CIR1 | 0,911 | 0,85–0,97 | 0,018 |
| STAT5B | 0,912 | 0,85–0,98 | 0,021 |
| IGF1 | 0,910 | 0,84–0,98 | 0,038 |
| AURKA | 1,107 | 1,02–1,20 | 0,038 |
| FLT3 | 0,923 | 0,86–0,99 | 0,038 |

A atenuação dos HR (por exemplo, *GSK3B* de 1,225 para 1,133) mostra que parte do efeito univariado é mediada pelo subtipo e pelo grau — mas não toda, o que caracteriza informação prognóstica parcialmente independente.

**Alerta metodológico (T14).** Entre os 20 genes mais prognósticos, **12 violaram o pressuposto de riscos proporcionais** (`cox.zph`, p < 0,05) e apenas 8 o atenderam. Com seguimento mediano de quase 10 anos, isso é esperado: o efeito de vários genes é forte nos primeiros anos e se dilui depois. Consequência prática: os HR reportados devem ser lidos como **efeitos médios ao longo do seguimento**, não como razões de risco constantes. Modelagem com efeitos dependentes do tempo ou restrição da janela de análise (por exemplo, OS em 5 anos) é a extensão natural desta pipeline.

### 5.6 Kaplan-Meier por expressão gênica (T17, F07)

| Gene | Mediana OS — expressão baixa | Mediana OS — expressão alta | Log-rank p |
|---|---|---|---|
| GSK3B | 186,4 meses | 124,8 meses | 2,8 × 10⁻¹⁰ |
| STAT5A | 124,8 meses | 187,9 meses | 2,5 × 10⁻¹⁰ |
| SPRY2 | 139,3 meses | 174,6 meses | 7,6 × 10⁻⁶ |
| IGF1 | 137,9 meses | 176,3 meses | 1,4 × 10⁻⁵ |

Alta expressão de *GSK3B* está associada a uma perda de aproximadamente **61 meses** de sobrevida mediana; alta expressão de *STAT5A*, a um ganho de aproximadamente **63 meses**.

### 5.7 Efeito prognóstico dentro de cada subtipo (T18, F08)

Este é um dos resultados mais informativos da análise. Testando os 12 genes mais prognósticos dentro de cada subtipo:

- **LumA** (679 pacientes, 364 eventos) — **todos os 12 genes** foram significativos, e os efeitos são maiores do que na coorte inteira: *PDGFRA* HR 0,695 (p = 5,4 × 10⁻¹¹), *SPRY2* 0,692, *CDKN2C* 0,746, *STAT5A* 0,757, *ABCB1* 0,741, *GSK3B* 1,273.
- **LumB** (461, 303 eventos) — 4 genes significativos: *GSK3B* 1,304, *FLT3* 0,798, *STAT5A* 0,829, *IGF1* 0,807.
- **Normal** (140, 76 eventos) — *AURKA* HR 1,835 (p = 5,8 × 10⁻⁴), *GSK3B* 1,388, *STAT5B* 0,674, *FLT3* 0,648.
- **claudin-low** (198, 89 eventos) — efeitos fracos, nenhum sobrevivendo à correção FDR.
- **Her2 e Basal** — **nenhum gene significativo**.

Interpretação: a informação prognóstica transcricional se concentra nos tumores **luminais**, onde existe heterogeneidade de risco a ser resolvida. Nos subtipos Her2 e Basal, o próprio subtipo já determina a maior parte do risco e a expressão desses genes acrescenta pouco — o que tem consequência direta para o desenho de painéis prognósticos: um painel derivado da coorte inteira será, na prática, um painel para doença luminal.

### 5.8 Assinatura multigênica e sua validação (T19–T22, T29–T32, F09, F10, F13)

O LASSO-Cox selecionou **65 genes** em `lambda.min` e **12 genes** em `lambda.1se`. Os maiores coeficientes em `lambda.min`: *HSD3B7* (−0,115), *CDKN2A* (−0,103), *CDKN2C* (−0,093), *FLT3* (−0,084), *MMP25* (−0,083), *SMAD6* (+0,080), *GSK3B* (+0,075). A assinatura parcimoniosa (`lambda.1se`) é: *FLT3*, *CDKN2C*, *SMAD6*, *STAT5A*, *GSK3B*, *MMP25*, *VEGFA*, *ABCB1*, *SPRY2*, *IGF1*, *CTCF*, *PDGFRA*.

**Desempenho em partição única (T20):**

| Modelo | Genes | C-index derivação | C-index validação |
|---|---|---|---|
| Assinatura LASSO (`lambda.min`) | 65 | 0,683 | 0,588 |
| Assinatura LASSO (`lambda.1se`) | 12 | 0,623 | 0,585 |
| Clínico | — | 0,657 | **0,661** |
| Combinado (`lambda.min`) | 65 | 0,711 | 0,651 |
| Combinado (`lambda.1se`) | 12 | 0,675 | 0,667 |

**Desempenho em 25 partições independentes (T30, T31):**

| Modelo | C-index médio | DP | Mín | Máx |
|---|---|---|---|---|
| Assinatura gênica | 0,5549 | 0,0354 | 0,500 | 0,611 |
| Clínico | 0,6553 | 0,0139 | 0,627 | 0,689 |
| Combinado | 0,6556 | 0,0142 | 0,634 | 0,691 |

Comparações pareadas (Wilcoxon): combinado *vs.* clínico → diferença média **+0,0003, p = 0,70** (sem ganho); assinatura *vs.* clínico → **−0,1004, p = 6 × 10⁻⁸** (a assinatura é inferior).

**Este é o resultado mais importante do relatório em termos metodológicos.** Na amostra de derivação, o teste de razão de verossimilhança indica que a assinatura acrescenta informação de forma espetacular (χ² = 161,5; 1 gl; p = 5,3 × 10⁻³⁷). Fora da amostra, esse ganho **desaparece por completo**. A explicação: com 489 candidatos e ~1.100 eventos, o LASSO captura estrutura específica da amostra de treino; além disso, o subtipo molecular já está no modelo clínico e ele próprio é um resumo da expressão gênica — a assinatura, em grande parte, apenas reexpressa informação já contida no subtipo e no grau. Relatar apenas o valor-p intra-amostral produziria a conclusão oposta e falsa.

A estratificação em tercis reproduz o mesmo padrão: na derivação a separação é ampla (33,8% / 60,5% / 79,6% de óbitos), mas na validação ela se comprime (48,5% / 61,1% / 64,7%).

**Estabilidade de seleção (T32, F13B).** Nas 25 repetições, apenas quatro genes foram selecionados com alguma consistência: **GSK3B (72%)**, **STAT5A (60%)**, *FLT3* (32%) e *SPRY2* (32%). Todos os demais ficaram em ≤ 12%. Ou seja, a lista de 65 genes é em grande parte instável — mas o núcleo *GSK3B* / *STAT5A* é reprodutível e coincide com os achados do Cox univariado e ajustado. Esse é o subconjunto que merece investigação subsequente.

### 5.9 Mutações somáticas (T23, T24, F11)

Dos 173 genes com dado mutacional, 71 alcançaram frequência ≥ 2%. Distribuição por subtipo dos principais (todos com FDR ≤ 5,1 × 10⁻⁴):

| Gene | Global | LumA | LumB | Her2 | Basal | claudin-low | Normal |
|---|---|---|---|---|---|---|---|
| PIK3CA | 41,8% | 57,4% | 34,9% | 41,4% | 16,1% | 24,2% | 50,0% |
| TP53 | 34,7% | 11,9% | 24,3% | 70,0% | **88,4%** | 52,0% | 22,9% |
| MUC16 | 17,2% | 16,6% | 16,1% | 27,3% | 23,1% | 9,6% | 10,0% |
| GATA3 | 12,1% | 19,6% | 13,9% | 7,7% | 0,0% | 2,0% | 8,6% |
| MAP3K1 | 10,3% | 16,2% | 8,7% | 6,8% | 4,5% | 4,0% | 10,0% |
| CDH1 | 9,1% | 12,7% | 9,8% | 5,0% | 2,0% | 4,0% | 12,9% |
| CBFB | 4,8% | 8,4% | 3,7% | 0,9% | 0,5% | 2,5% | 6,4% |
| RB1 | 2,6% | 1,5% | 2,4% | 3,2% | 7,0% | 3,5% | 0,0% |

O eixo *TP53*-mutado/*PIK3CA*-selvagem (Basal, Her2) versus *PIK3CA*-mutado/*TP53*-selvagem (LumA) reproduz com precisão a dicotomia estabelecida na literatura de câncer de mama — uma validação externa implícita da pipeline.

Associação com sobrevida global (Cox univariado):

| Gene mutado | Frequência | HR | IC 95% | FDR |
|---|---|---|---|---|
| GATA3 | 12,1% | 0,590 | 0,48–0,72 | 3,3 × 10⁻⁵ |
| CBFB | 4,8% | 0,501 | 0,36–0,70 | 1,3 × 10⁻³ |
| TP53 | 34,7% | 1,289 | 1,14–1,46 | 1,3 × 10⁻³ |

Nenhum outro gene mutado sobreviveu à correção BH. Ressalva importante: mutação em *GATA3* e *CBFB* concentra-se em tumores LumA, e mutação em *TP53* em Basal/Her2 — de modo que esses HR são **fortemente confundidos pelo subtipo** e não devem ser lidos como efeito causal independente sem modelagem estratificada.

### 5.10 Integração: subtipo × prognóstico (T25, T26, F12)

O cruzamento produziu **223 pares gene–subtipo** que são simultaneamente marcadores de subtipo e prognósticos. Os prioritários por subtipo (escore = |Δz| × |log HR| × −log₁₀ FDR):

| Subtipo | Genes prioritários | Leitura |
|---|---|---|
| **Basal** | *AURKA* (↑, HR 1,19), *FANCD2* (↑, HR 1,16), *IGF1* (↓, HR 0,80), *LAMA2* (↓, HR 0,84), *BCL2* (↓, HR 0,85) | Todos coerentes com **pior prognóstico**: o Basal superexpressa genes de risco e perde genes protetores |
| **Her2** | *GSK3B* (↑, HR 1,23), *SPRY2* (↓), *LAMA2* (↓), *BCL2* (↓), *STAT5B* (↓) | Padrão de pior prognóstico dominado pela superexpressão de *GSK3B* |
| **LumB** | *SPRY2* (↓), *ABCB1* (↓), *STAT5A* (↓), *PDGFRA* (↓), *CCND2* (↓) | Pior prognóstico por **perda de genes protetores**, não por ganho de genes de risco |
| **LumA** | *AURKA* (↓), *BCL2* (↑), *FANCD2* (↓), *E2F7* (↓), *E2F2* (↓) | Perfil de **melhor prognóstico** — o espelho exato do Basal |
| **claudin-low** | *STAT5A* (↑), *ABCB1* (↑), *IGF1* (↑), *CDKN2C* (↑), *CCND2* (↑) | Melhor prognóstico, consistente com a maior sobrevida mediana observada |
| **Normal** | *SPRY2* (↑), *ABCB1* (↑), *LAMA2* (↑), *AURKA* (↓), *KIT* (↑) | Melhor prognóstico |

Note a simetria entre Basal e LumA: os mesmos genes (*AURKA*, *FANCD2*, *E2F2*, *BCL2*) aparecem nos dois, com direções invertidas. Isso indica que o eixo prognóstico dominante na coorte não é um conjunto de genes independentes, mas um **gradiente único de proliferação** — o que também explica por que a assinatura LASSO acrescenta tão pouco ao modelo que já contém subtipo e grau.

---

## 6. Resposta direta à pergunta da análise

**Quais genes estão mais associados a determinados subtipos?**
*ERBB2* (Her2); *CCNE1*, *CDKN2A*, *CHEK1*, *E2F3*, *CDC25A* (Basal); *GATA3*, *MAPT*, *BCL2* (LumA); *FOLR2*, *KLRG1*, *CSF1R*, *TGFBR2* com perda de *ERBB3*/*RAB25* (claudin-low); perda de *EGFR*, *NOTCH1*, *FOXO1*, *SPRY2* (LumB). No plano mutacional: *TP53* (Basal, Her2), *PIK3CA*/*GATA3*/*MAP3K1*/*CBFB* (LumA).

**Quais genes estão mais associados a pior prognóstico?**
Consistentes em Cox univariado, Cox ajustado e seleção estável por LASSO: **GSK3B** (HR 1,23 bruto; 1,13 ajustado; selecionado em 72% das repetições), **AURKA** (1,19 bruto; 1,46 para mortalidade específica), **VEGFA** (1,17; 1,13 ajustado), *FANCD2*, *E2F7*, *E2F2*, *CCNE1*. Como marcadores de **melhor** prognóstico: **STAT5A** (HR 0,81; selecionado em 60% das repetições), *SPRY2*, *IGF1*, *ABCB1*, *BCL2*, *STAT5B*, *PDGFRA*, *CDKN2C*.

**Com que grau de confiança?**
Alta para a existência das associações (FDR rigoroso, n = 1.897, efeitos replicados entre OS e DSS e entre métodos independentes). **Baixa** para a utilidade preditiva incremental: nenhuma combinação desses genes superou, fora da amostra, um modelo clínico simples baseado em idade, grau, tamanho, linfonodos e subtipo.

---

## 7. Limitações

1. **Cobertura gênica parcial.** Foram analisados 489 genes de um painel curado, não o transcriptoma completo (~24.000 genes) do METABRIC — este está disponível no cBioPortal apenas via Git LFS/S3, fora do alcance da rede deste ambiente. Genes prognósticos fora do painel não podem ser detectados, e as análises de enriquecimento de vias ficam enviesadas pela composição do painel.
2. **Pressuposto de riscos proporcionais violado** em 12 dos 20 genes mais fortes; os HR representam efeitos médios no seguimento.
3. **Ausência de validação externa independente.** Toda a validação foi interna (partições repetidas). Confirmação em TCGA-BRCA, GSE ou coorte própria é o passo seguinte obrigatório.
4. **Confundimento por tratamento.** A coorte é histórica (pré-trastuzumabe, esquemas adjuvantes heterogêneos); prognóstico e resposta terapêutica não são separáveis aqui.
5. **Dicotomização pela mediana** nas curvas KM tem finalidade descritiva; a inferência formal usa a expressão contínua.
6. **Mortalidade competitiva.** A OS inclui 480 óbitos por outras causas em uma coorte de mediana etária elevada; a análise de DSS mitiga, mas o tratamento formal exigiria modelos de risco competitivo (Fine-Gray).
7. **Grupo "Normal-like"** provavelmente reflete contaminação por tecido normal (80,7% de erro de classificação), e resultados nesse estrato devem ser interpretados com cautela.
8. **z-scores calculados sobre a coorte inteira** implicam que os valores de cada amostra dependem da composição do estudo, limitando a aplicação direta do escore a um paciente individual sem recalibração.

---

## 8. Dicionário completo das saídas geradas

### 8.1 Tabelas (`saidas/tabelas/`, 33 arquivos CSV)

| Arquivo | Conteúdo |
|---|---|
| T01_integridade_dados | Verificações de integridade e checksum MD5 do arquivo de origem |
| T02_auditoria_codificacao_evento | Tabulação cruzada que define a codificação do evento de óbito |
| T03_fluxograma_exclusoes | Fluxograma de exclusões com n em cada etapa |
| T04_qc_expressao_por_gene | Média, DP, MAD, mínimo e máximo dos 489 genes |
| T05_caracteristicas_clinicas_por_subtipo | 17 variáveis clínicas resumidas por subtipo |
| T06_testes_associacao_clinica_subtipo | Qui-quadrado e Kruskal-Wallis com FDR |
| T07_expressao_diferencial_completa_todos_subtipos | **2.934 linhas** — resultado completo do limma (489 genes × 6 subtipos) |
| T08_resumo_expressao_diferencial | Contagem de genes significativos por subtipo |
| T09_top15_genes_por_subtipo | 15 genes de maior efeito em cada subtipo |
| T10_randomforest_matriz_confusao | Matriz de confusão out-of-bag e erro por classe |
| T11_randomforest_importancia_genes | Ranking completo dos 489 genes por importância |
| T12_cox_univariado_OS_todos_genes | Cox de sobrevida global para os 489 genes (HR, IC 95%, p, FDR, concordância) |
| T13_cox_univariado_DSS_todos_genes | Idem para sobrevida específica do câncer |
| T14_teste_riscos_proporcionais_top20 | Teste `cox.zph` dos 20 genes mais prognósticos |
| T15_cox_multivariado_ajustado_top20 | Cox ajustado por 8 covariáveis clínicas |
| T16_mediana_sobrevida_por_subtipo | Medianas de OS com IC 95% |
| T17_logrank_top4_genes | Log-rank dos 4 genes principais dicotomizados |
| T18_cox_estratificado_por_subtipo | 12 genes × 6 subtipos: HR e p dentro de cada estrato |
| T19_assinatura_lasso_cox_coeficientes | 65 genes e coeficientes em `lambda.min` |
| T19b_assinatura_lasso_lambda1se | 12 genes da assinatura parcimoniosa |
| T20_desempenho_modelos_C_index | C-index dos 5 modelos em partição única |
| T21_teste_razao_verossimilhanca | LRT clínico vs. clínico+assinatura |
| T22_grupos_de_risco_por_conjunto | Distribuição de óbitos por tercil de risco |
| T23_mutacoes_por_subtipo | Frequência mutacional de 71 genes por subtipo, Fisher + FDR |
| T24_cox_mutacoes_OS | Cox de sobrevida para cada gene mutado |
| T25_integracao_subtipo_prognostico | **223 pares** gene–subtipo com HR, coerência direcional e anotações |
| T26_genes_prioritarios_por_subtipo | Top 10 por subtipo pelo escore de prioridade |
| T27_resumo_execucao | Metadados e métricas-chave da execução |
| T28_bibliotecas_versoes | 22 pacotes com versão exata |
| T29_validacao_repetida_25_particoes | C-index das 3 abordagens em cada uma das 25 repetições |
| T30_resumo_validacao_repetida | Média, DP, mínimo e máximo do C-index |
| T31_comparacao_pareada_modelos | Testes de Wilcoxon pareados entre modelos |
| T32_estabilidade_selecao_lasso | Frequência de seleção de cada gene nas 25 repetições |

### 8.2 Figuras (`saidas/figuras/`, 13 arquivos PNG a 150–170 dpi)

| Arquivo | Conteúdo |
|---|---|
| F01_distribuicao_subtipos | Distribuição dos subtipos e mortalidade bruta |
| F02_volcano_por_subtipo | Seis volcano plots com genes rotulados |
| F03_heatmap_genes_discriminantes | Heatmap com clusterização dos genes mais discriminantes |
| F04_importancia_randomforest | Top 30 genes por importância no Random Forest |
| F05_forest_cox_univariado_OS | Forest plot dos 25 genes mais prognósticos |
| F06_km_subtipos | Kaplan-Meier por subtipo com tabela de risco |
| F07_km_top4_genes | Kaplan-Meier dos 4 genes principais |
| F08_heatmap_cox_estratificado | HR por gene dentro de cada subtipo |
| F09_km_grupos_risco_validacao | Validação independente dos tercis de risco |
| F10_coeficientes_assinatura_lasso | Coeficientes dos 65 genes da assinatura |
| F11_mutacoes_por_subtipo | Frequência mutacional por subtipo |
| F12_integracao_subtipo_prognostico | Δz vs. log₂(HR) por subtipo |
| F13_estabilidade_validacao_repetida | Distribuição do C-index e estabilidade de seleção |

### 8.3 Demais artefatos

- `metabric_pipeline.R` — script único, 764 linhas, comentado, executável de ponta a ponta
- `saidas/objetos_analise.rds` — todos os objetos da análise para reanálise sem reexecutar
- `saidas/logs/execucao.log` — log completo com marcação de tempo por módulo
- `saidas/logs/sessionInfo.txt` — ambiente exato para reprodução

---

## 9. Como reproduzir

O script é **autossuficiente**: não exige edição de caminhos, instalação prévia de pacotes nem download manual de dados. Em qualquer computador com R ≥ 4.0 (Linux, macOS ou Windows):

```bash
Rscript metabric_pipeline.R
```

O que ele faz sozinho, nesta ordem:

1. **Detecta o próprio diretório** (via `--file=` do Rscript, via `ofile` quando usado com `source()`, ou o diretório atual em sessão interativa) e cria a estrutura `dados/` e `saidas/{tabelas,figuras,logs}`.
2. **Configura o locale UTF-8** testando os nomes usados em Linux, macOS e Windows — sem isso os acentos quebram a renderização das figuras do `survminer`.
3. **Instala as dependências ausentes**: 21 pacotes do CRAN (`https://cloud.r-project.org`) e o `limma` via `BiocManager`. Pacotes já presentes são preservados. Se a instalação falhar, a mensagem de erro sugere o comando `apt-get` equivalente para Debian/Ubuntu, que evita compilação.
4. **Baixa o arquivo de entrada** (`METABRIC_RNA_Mutation.csv`, 8,1 MB) de `raw.githubusercontent.com` e **verifica o MD5** contra a referência `c619471beb87af0f9fd4e4a40058654d`. O resultado da verificação e a URL efetivamente usada ficam registrados em `T01_integridade_dados.csv`.
5. **Executa toda a análise** e grava 33 tabelas, 13 figuras, os logs e o `.rds` com os objetos.

Comportamento em execuções repetidas e em falhas — testado explicitamente:

| Situação | Comportamento verificado |
|---|---|
| Diretório limpo, só com o script | Baixa os dados e executa do início ao fim |
| Arquivo já presente e íntegro | Confirma o MD5 e pula o download |
| Arquivo local corrompido ou truncado | Detecta pela validação estrutural, descarta e rebaixa |
| MD5 divergente mas estrutura válida | Emite aviso e prossegue (tolera atualização do espelho) |
| Todas as fontes indisponíveis | Falha com instrução explícita de download manual e o caminho exato de destino |
| Sem rede, arquivo posto manualmente em `dados/` | Valida e executa normalmente, sem tentar baixar |

Variáveis de ambiente opcionais:

```bash
METABRIC_DIR=/caminho/alternativo Rscript metabric_pipeline.R   # outro diretório de trabalho
METABRIC_NREP=3 Rscript metabric_pipeline.R                     # teste rápido (~2,6 min)
```

**Tempo de execução:** ~6,6 min com o padrão de 25 repetições, em 1 vCPU e 3 GB de RAM. O Módulo 8b responde por cerca de 4,7 min desse total; com `METABRIC_NREP=3` a execução completa cai para ~2,6 min (útil para validar o ambiente antes da rodada definitiva). Todas as etapas com aleatoriedade — Random Forest, validação cruzada do LASSO e partições — usam sementes fixas, de modo que os números deste relatório são reproduzidos bit a bit.

**Validação estrutural do arquivo baixado.** Antes de qualquer análise, o script confere a presença das colunas críticas (`patient_id`, `pam50_+_claudin-low_subtype`, `overall_survival`, `overall_survival_months`, `death_from_cancer`, `tp53`, `erbb2`), o número mínimo de colunas e o tamanho do arquivo — o que evita que uma página de erro HTML salva como CSV seja processada silenciosamente como se fossem dados.

## 10. Próximos passos sugeridos

1. **Obter a matriz de expressão completa** (~24.000 genes) do cBioPortal em rede sem restrição, repetindo os Módulos 4, 6 e 8 — a cobertura atual de 489 genes é a limitação mais séria.
2. **Modelagem tempo-dependente** (`tt()` no `coxph` ou modelos de tempo de falha acelerado) para os 12 genes que violam riscos proporcionais.
3. **Restringir o escopo aos tumores luminais**, onde toda a informação prognóstica transcricional se concentra (seção 5.7), e rederivar a assinatura nesse contexto — é a hipótese mais promissora para superar o modelo clínico.
4. **Validação externa** em TCGA-BRCA, com atenção à diferença de plataforma (microarranjo vs. RNA-seq).
5. **Análise de enriquecimento funcional** dos 125 genes prognósticos, com correção pelo viés de composição do painel.
6. **Modelos de risco competitivo** (Fine-Gray) dada a alta proporção de óbitos por outras causas.

---

## 11. Referências

- Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*. 2012;486(7403):346-352.
- Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. *Nat Commun*. 2016;7:11479.
- Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47.
- Simon N, Friedman J, Hastie T, Tibshirani R. Regularization paths for Cox's proportional hazards model via coordinate descent. *J Stat Softw*. 2011;39(5):1-13.
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. *J R Stat Soc B*. 1995;57(1):289-300.
- Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model*. Springer; 2000.
- Meinshausen N, Bühlmann P. Stability selection. *J R Stat Soc B*. 2010;72(4):417-473.
