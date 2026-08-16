# Resumo para submissão — Congresso AMAZOOmicas

**Eixo temático sugerido:** Bioinformática / Genômica aplicada à saúde · **Formato:** pôster ou comunicação oral

---

## APRENDIZADO DE MÁQUINA APLICADO A DADOS TRANSCRIPTÔMICOS IDENTIFICA GENES ASSOCIADOS À SOBREVIDA EM SUBTIPOS E ESTÁGIOS DO CÂNCER DE MAMA: ANÁLISE DE 1.897 PACIENTES DA COORTE METABRIC

**Autores:** Micael Davi Lima de Oliveira¹; [coautores]; Pritesh J. Lalwani¹

¹ Laboratório de Imunologia e Doenças Infecciosas (LABIDI), Instituto Leônidas e Maria Deane, Fiocruz Amazônia, Manaus, AM, Brasil

---

### Resumo

**Introdução.** O câncer de mama é uma doença molecularmente heterogênea, e a decisão terapêutica depende de estratificar corretamente o risco de cada paciente. Métodos de aprendizado de máquina permitem interrogar centenas de genes simultaneamente, sem assumir linearidade ou efeitos constantes no tempo, mas seu ganho real sobre variáveis clínicas de rotina raramente é testado fora da amostra de derivação.

**Objetivo.** Identificar, por aprendizado de máquina e modelos de sobrevida, os genes associados a maior e menor sobrevida global em diferentes subtipos moleculares e estágios clínicos do câncer de mama, e quantificar o ganho preditivo desses genes sobre um modelo clínico.

**Métodos.** Análise reprodutível de dados reais da coorte METABRIC (n = 1.897; 489 genes de expressão em z-score; 1.098 óbitos; seguimento mediano de 115,6 meses), em pipeline automatizado em R e Python. Foram aplicados: **Random Forest** (500 árvores) para classificação dos seis subtipos moleculares, com erro *out-of-bag* e ranking de importância por índice de Gini; **Random Forest** para predição de óbito em até 10 anos, com validação cruzada estratificada 5-fold e comparação entre modelo transcriptômico, clínico e combinado; **modelos de Cox** univariados para os 489 genes, estratificados por subtipo e por estágio tumoral; e **LASSO-Cox** com validação repetida em 25 partições independentes e análise de estabilidade de seleção. Correção de Benjamini-Hochberg em todas as etapas; sementes fixas garantindo reprodutibilidade integral.

**Resultados.** O Random Forest classificou os seis subtipos com **76,1% de acurácia** (erro OOB 23,9%), destacando *GATA3*, *EGFR*, *CDK1*, *AURKA* e *MAPT* como genes mais informativos; o erro concentrou-se no grupo Normal-like (80,7%, com 93/140 casos reclassificados como LumA), reforçando a hipótese de contaminação por tecido normal. Na predição de óbito em 10 anos, os genes mais importantes foram **AURKA, STAT5A, FLT3, DIRAS3, STAT5B e GSK3B** — conjunto quase disjunto do que define os subtipos, indicando que os genes que caracterizam o tumor não são os que predizem o desfecho. Cento e vinte e cinco genes associaram-se à sobrevida global (FDR < 0,05): **menor sobrevida** para *GSK3B* (HR 1,23), *AURKA* (1,19), *VEGFA* (1,17), *FANCD2* e *E2F2*; **maior sobrevida** para *STAT5A* (HR 0,81), *SPRY2* (0,81), *IGF1* (0,80), *ABCB1* (0,82), *FLT3* (0,82), *CDKN2C* (0,82) e *PDGFRA* (0,84). O sinal mostrou-se **fortemente dependente do contexto**: no LumA, 12 de 12 genes testados foram significativos (*PDGFRA* HR 0,70; p = 5,4 × 10⁻¹¹), contra nenhum em Basal e Her2; por estágio, 71 genes atingiram significância no estágio II (*ACVR1B* 1,28; *STAT5A* 0,79; *IGF1* 0,76) contra apenas um no estágio I (*NRIP1* 0,73) e um no III–IV (*DNAH11* 1,52). O desempenho do Random Forest reproduziu esse gradiente (AUC 0,666 nos luminais; 0,558 no Her2; 0,420 no Basal). Contudo, o modelo transcriptômico **não superou o clínico**: AUC 0,653 contra 0,730 (combinado 0,715), e C-index de 0,555 contra 0,655 em 25 partições (p = 6 × 10⁻⁸). Apenas *GSK3B* (72%) e *STAT5A* (60%) mostraram seleção estável pelo LASSO.

**Conclusão.** O aprendizado de máquina recuperou de forma robusta e convergente entre três famílias metodológicas os genes associados a maior e menor sobrevida, mas mostrou que essa informação **se concentra nos tumores luminais e no estágio II** e não acrescenta capacidade preditiva sobre variáveis clínicas de rotina. O achado orienta o desenho de painéis prognósticos: em vez de assinaturas derivadas de coortes inteiras, o esforço deve concentrar-se nos estratos onde ainda há heterogeneidade de risco a resolver.

**Palavras-chave:** câncer de mama; aprendizado de máquina; Random Forest; análise de sobrevida; METABRIC.

---

### Notas para ajuste antes da submissão

- **Extensão:** o corpo do resumo (Introdução → Conclusão) tem cerca de **480 palavras**. Congressos costumam limitar a 250, 300 ou 500 — se o limite for 250–300, os cortes naturais, nesta ordem, são: (1) a frase sobre o grupo Normal-like, (2) a lista completa de genes protetores, reduzida aos quatro primeiros, (3) os detalhes da estabilidade de seleção pelo LASSO.
- **Título alternativo mais curto:** "Random Forest e modelos de sobrevida identificam genes prognósticos contexto-dependentes no câncer de mama (METABRIC, n = 1.897)".
- **Autoria:** conferir a ordem e completar os coautores conforme a contribuição de cada um; verificar o limite de coautores do evento.
- **Referências:** a maioria dos congressos não aceita referências no corpo do resumo; nenhuma foi incluída.
- **Camada farmacogenômica:** deliberadamente fora deste resumo, que foi centrado em aprendizado de máquina conforme solicitado. Se o eixo temático do AMAZOOmicas favorecer farmacogenômica, há material para um **segundo resumo** com os 14 genes que modificam o efeito da hormonioterapia (*EIF4E*, *TWIST1*, *MYC*, *MTOR*, *SETDB1* e outros), todos resistentes ao ajuste por status de ER e subtipo — o que daria dois trabalhos independentes a partir da mesma análise.
