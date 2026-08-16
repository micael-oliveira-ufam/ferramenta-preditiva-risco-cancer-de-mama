###############################################################################
# PIPELINE METABRIC — Genes associados a subtipos moleculares e prognóstico
# em câncer de mama
#
# Coorte: METABRIC (Molecular Taxonomy of Breast Cancer International
#         Consortium; Curtis et al., Nature 2012; Pereira et al., Nat Commun 2016)
# Dados:  arquivo curado a partir do cBioPortal (estudo brca_metabric),
#         1.904 pacientes, 31 variáveis clínicas, 489 genes com expressão
#         (z-scores de microarranjo Illumina HT-12 v3) e 173 genes com
#         status mutacional (painel de sequenciamento alvo).
#
# ---------------------------------------------------------------------------
# EXECUÇÃO — o script é autossuficiente em qualquer computador com R >= 4.0:
#
#   Rscript metabric_pipeline.R
#
# Ele (1) detecta sozinho o próprio diretório, (2) instala as dependências
# ausentes a partir do CRAN e do Bioconductor, (3) BAIXA o arquivo de dados
# de entrada e verifica o checksum MD5, e (4) executa toda a análise,
# gravando tabelas, figuras e logs em ./saidas/.
#
# Nada precisa ser editado. Variáveis de ambiente opcionais:
#   METABRIC_DIR   — diretório de trabalho alternativo
#   METABRIC_NREP  — número de repetições da validação (padrão 25; use 3 para
#                    um teste rápido de ponta a ponta)
#
# Requisito de rede: acesso HTTPS a raw.githubusercontent.com (dados) e aos
# repositórios CRAN/Bioconductor (apenas se faltarem pacotes). Em ambientes
# sem rede, coloque manualmente o arquivo METABRIC_RNA_Mutation.csv em
# ./dados/ — o script detecta, valida e prossegue sem baixar nada.
#
# ---------------------------------------------------------------------------
# Estrutura:
#   MÓDULO 0 — Ambiente, dependências e reprodutibilidade
#   MÓDULO 1 — Aquisição automática dos dados e verificação de integridade
#   MÓDULO 2 — Controle de qualidade e pré-processamento
#   MÓDULO 3 — Caracterização clínico-descritiva por subtipo
#   MÓDULO 4 — Expressão diferencial por subtipo (limma)
#   MÓDULO 5 — Classificação supervisionada e importância de genes (Random Forest)
#   MÓDULO 6 — Sobrevida univariada gene a gene (Cox proporcional de riscos)
#   MÓDULO 7 — Kaplan-Meier (subtipos e genes de maior efeito)
#   MÓDULO 8 — Assinatura prognóstica multigênica (LASSO-Cox + validação)
#   MÓDULO 8b — Estabilidade e validação repetida
#   MÓDULO 9 — Mutações somáticas por subtipo e prognóstico
#   MÓDULO 10 — Integração final e exportação
#
# Autor: pipeline gerado para análise reprodutível
###############################################################################

## =========================================================================
## MÓDULO 0 — AMBIENTE, DEPENDÊNCIAS E REPRODUTIBILIDADE
## =========================================================================
t_inicio <- Sys.time()

options(stringsAsFactors = FALSE, scipen = 6,
        timeout = max(900, getOption("timeout")))   # downloads grandes em redes lentas

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# --- 0.1 Locale UTF-8 (necessário para acentos nas figuras; multiplataforma) --
# Tenta os nomes usados em Linux, macOS e Windows, nesta ordem, sem falhar.
invisible(suppressWarnings({
  for (loc in c("C.UTF-8", "en_US.UTF-8", "pt_BR.UTF-8", "C.utf8",
                "Portuguese_Brazil.utf8", "English_United States.utf8")) {
    if (!identical(try(Sys.setlocale("LC_ALL", loc), silent = TRUE), "")) break
  }
}))

# --- 0.2 Detecção automática do diretório do script --------------------------
# Funciona via Rscript (--file=), via source() e em sessão interativa.
detectar_diretorio <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  alvo <- grep("^--file=", args, value = TRUE)
  if (length(alvo) > 0) return(dirname(normalizePath(sub("^--file=", "", alvo[1]))))
  caminho <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(caminho)) return(dirname(normalizePath(caminho)))
  getwd()
}
DIR <- Sys.getenv("METABRIC_DIR", unset = detectar_diretorio())

DIR_DADOS <- file.path(DIR, "dados")
DIR_TAB   <- file.path(DIR, "saidas", "tabelas")
DIR_FIG   <- file.path(DIR, "saidas", "figuras")
DIR_LOG   <- file.path(DIR, "saidas", "logs")
for (d in c(DIR_DADOS, DIR_TAB, DIR_FIG, DIR_LOG))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
msg("Diretório de trabalho: ", DIR)

# --- 0.3 Instalação automática das dependências ------------------------------
# Bibliotecas e finalidade de cada uma:
#   data.table   — leitura rápida (fread) e escrita das tabelas de saída
#   dplyr        — manipulação tabular (filter/mutate/summarise)
#   tidyr        — reformatação wide<->long (pivot_longer/pivot_wider)
#   purrr        — iteração funcional (map_dfr) nos loops de modelagem
#   stringr      — manipulação de strings (nomes de genes, rótulos)
#   forcats      — reordenação de fatores em gráficos
#   survival     — Surv, Kaplan-Meier, log-rank, Cox, cox.zph, concordance
#   survminer    — curvas de sobrevida publicáveis com tabela de risco
#   glmnet       — Cox penalizado por LASSO (assinatura multigênica)
#   randomForest — classificação multiclasse e importância de variáveis
#   matrixStats  — estatísticas por linha da matriz de expressão (QC)
#   broom        — padronização da saída de modelos em data.frame
#   ggplot2      — motor gráfico de todas as figuras
#   ggrepel      — rótulos não sobrepostos nos volcano plots
#   scales       — formatação de eixos e escalas
#   RColorBrewer — paletas categóricas e divergentes
#   viridis      — paletas contínuas perceptualmente uniformes
#   pheatmap     — heatmaps com clusterização hierárquica
#   patchwork    — composição de múltiplos painéis
#   cowplot      — tema e alinhamento de painéis
#   gridExtra    — manipulação de objetos grid na exportação
#   limma        — (Bioconductor) modelos lineares + Bayes empírico para
#                  expressão diferencial
PKGS_CRAN <- c("data.table", "dplyr", "tidyr", "purrr", "stringr", "forcats",
               "survival", "survminer", "glmnet", "randomForest", "matrixStats",
               "broom", "ggplot2", "ggrepel", "scales", "RColorBrewer",
               "viridis", "pheatmap", "patchwork", "cowplot", "gridExtra")
PKGS_BIOC <- c("limma")

instalar_dependencias <- function(cran, bioc) {
  presente <- function(p) requireNamespace(p, quietly = TRUE)
  repos <- getOption("repos")
  if (is.null(repos[["CRAN"]]) || repos[["CRAN"]] %in% c("@CRAN@", ""))
    repos <- c(CRAN = "https://cloud.r-project.org")

  falta_cran <- cran[!vapply(cran, presente, logical(1))]
  if (length(falta_cran) > 0) {
    msg("Instalando do CRAN: ", paste(falta_cran, collapse = ", "))
    install.packages(falta_cran, repos = repos, dependencies = TRUE)
  }

  falta_bioc <- bioc[!vapply(bioc, presente, logical(1))]
  if (length(falta_bioc) > 0) {
    if (!presente("BiocManager")) install.packages("BiocManager", repos = repos)
    msg("Instalando do Bioconductor: ", paste(falta_bioc, collapse = ", "))
    BiocManager::install(falta_bioc, ask = FALSE, update = FALSE)
  }

  ausentes <- c(cran, bioc)[!vapply(c(cran, bioc), presente, logical(1))]
  if (length(ausentes) > 0) {
    stop("Não foi possível instalar: ", paste(ausentes, collapse = ", "),
         "\nEm Debian/Ubuntu, uma alternativa sem compilação é:\n  sudo apt-get install -y ",
         paste0("r-cran-", tolower(setdiff(ausentes, "limma"))," ", collapse = ""),
         if ("limma" %in% ausentes) "r-bioc-limma" else "", call. = FALSE)
  }
  invisible(TRUE)
}
instalar_dependencias(PKGS_CRAN, PKGS_BIOC)

suppressPackageStartupMessages({
  library(data.table); library(dplyr);     library(tidyr);        library(purrr)
  library(stringr);    library(forcats);   library(limma);        library(survival)
  library(survminer);  library(glmnet);    library(randomForest); library(matrixStats)
  library(broom);      library(ggplot2);   library(ggrepel);      library(scales)
  library(RColorBrewer); library(viridis); library(pheatmap);     library(patchwork)
  library(cowplot);    library(gridExtra)
})

set.seed(42)            # semente global: reprodutibilidade de RF, LASSO e partições

# --- 0.4 Parâmetros analíticos globais (declarados a priori) -----------------
PAR <- list(
  fdr_de        = 0.01,   # limiar de FDR para expressão diferencial
  lfc_de        = 0.50,   # diferença mínima em unidades de z-score (efeito)
  fdr_cox       = 0.05,   # limiar de FDR para Cox univariado
  prop_treino   = 0.70,   # proporção da coorte no conjunto de derivação
  n_arvores     = 500,    # árvores no Random Forest
  n_folds_cv    = 10,     # dobras da validação cruzada do LASSO-Cox
  freq_min_mut  = 0.02    # frequência mínima de mutação para teste (2%)
)

salvar_tab <- function(x, nome) {
  fwrite(as.data.frame(x), file.path(DIR_TAB, paste0(nome, ".csv")))
  invisible(x)
}
salvar_fig <- function(plot, nome, w = 9, h = 6, dpi = 150) {
  ggsave(file.path(DIR_FIG, paste0(nome, ".png")), plot, width = w, height = h,
         dpi = dpi, bg = "white")
  invisible(NULL)
}

msg("MÓDULO 0 concluído — ", R.version.string, " | ", length(c(PKGS_CRAN, PKGS_BIOC)),
    " pacotes disponíveis")

## =========================================================================
## MÓDULO 1 — AQUISIÇÃO AUTOMÁTICA DOS DADOS E VERIFICAÇÃO DE INTEGRIDADE
## =========================================================================
# O arquivo é um recorte curado do estudo brca_metabric do cBioPortal.
# A matriz de expressão completa do METABRIC (~24.000 genes) está no
# cBioPortal apenas via Git LFS/S3; este recorte de 489 genes é a versão
# publicamente redistribuída e com download direto.
ARQUIVO <- file.path(DIR_DADOS, "METABRIC_RNA_Mutation.csv")
MD5_REF <- "c619471beb87af0f9fd4e4a40058654d"
FONTES  <- c(
  "https://raw.githubusercontent.com/thomas-smithh/metabric-RNA-mutation/master/Data/METABRIC_RNA_Mutation.csv"
)

validar_estrutura <- function(caminho) {
  if (!file.exists(caminho) || file.size(caminho) < 5e6) return(FALSE)
  cab <- tryCatch(names(fread(caminho, nrows = 0)), error = function(e) character(0))
  obrigatorias <- c("patient_id", "pam50_+_claudin-low_subtype", "overall_survival",
                    "overall_survival_months", "death_from_cancer", "tp53", "erbb2")
  all(obrigatorias %in% cab) && length(cab) > 600
}

baixar_dados <- function(destino, urls, md5_ref) {
  if (file.exists(destino)) {
    md5 <- unname(tools::md5sum(destino))
    if (identical(md5, md5_ref)) {
      msg("Arquivo já presente e verificado por MD5 — download dispensado")
      return(invisible("cache"))
    }
    if (validar_estrutura(destino)) {
      warning("MD5 diferente do de referência (", md5, "), mas a estrutura do ",
              "arquivo é válida. Prosseguindo com o arquivo local.", call. = FALSE)
      return(invisible("cache_md5_divergente"))
    }
    msg("Arquivo local inválido — será rebaixado")
    file.remove(destino)
  }
  for (u in urls) {
    msg("Baixando dados de: ", u)
    ok <- tryCatch({
      utils::download.file(u, destfile = destino, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) { msg("  falhou: ", conditionMessage(e)); FALSE },
       warning = function(w) { msg("  aviso: ", conditionMessage(w)); FALSE })
    if (isTRUE(ok) && validar_estrutura(destino)) {
      msg("Download concluído (", round(file.size(destino) / 1024^2, 1), " MB)")
      return(invisible(u))
    }
    if (file.exists(destino)) file.remove(destino)
  }
  stop("Não foi possível obter o arquivo de dados automaticamente.\n",
       "Verifique o acesso HTTPS a raw.githubusercontent.com ou baixe o arquivo\n",
       "METABRIC_RNA_Mutation.csv manualmente (Kaggle: 'Breast Cancer Gene\n",
       "Expression Profiles (METABRIC)') e coloque-o em:\n  ", destino, call. = FALSE)
}
origem <- baixar_dados(ARQUIVO, FONTES, MD5_REF)

bruto <- fread(ARQUIVO, na.strings = c("", "NA"), showProgress = FALSE)
msg("MÓDULO 1 — arquivo lido: ", nrow(bruto), " pacientes x ", ncol(bruto), " colunas")

# Definição dos blocos de colunas (layout do arquivo cBioPortal curado)
cols <- names(bruto)
col_clin <- cols[1:31]                              # variáveis clínico-patológicas
col_mut  <- grep("_mut$", cols, value = TRUE)       # status mutacional (173 genes)
col_expr <- setdiff(cols[32:length(cols)], col_mut) # expressão (z-scores)
stopifnot(length(col_expr) > 400, length(col_mut) > 100)

integridade <- data.frame(
  item = c("n_pacientes", "n_colunas_total", "n_variaveis_clinicas",
           "n_genes_expressao", "n_genes_mutacao", "ids_duplicados",
           "NAs_matriz_expressao", "min_zscore", "max_zscore",
           "md5_arquivo_origem", "md5_confere_com_referencia", "origem_do_arquivo"),
  valor = c(nrow(bruto), ncol(bruto), length(col_clin), length(col_expr),
            length(col_mut), sum(duplicated(bruto$patient_id)),
            sum(is.na(bruto[, ..col_expr])),
            round(min(as.matrix(bruto[, ..col_expr]), na.rm = TRUE), 4),
            round(max(as.matrix(bruto[, ..col_expr]), na.rm = TRUE), 4),
            unname(tools::md5sum(ARQUIVO)),
            identical(unname(tools::md5sum(ARQUIVO)), MD5_REF),
            as.character(origem))
)
salvar_tab(integridade, "T01_integridade_dados")
print(integridade)

## =========================================================================
## MÓDULO 2 — CONTROLE DE QUALIDADE E PRÉ-PROCESSAMENTO
## =========================================================================
dados <- as.data.frame(bruto)

# --- 2.1 Desfechos de sobrevida ---------------------------------------------
# Auditoria da codificação: no arquivo, overall_survival = 1 corresponde a
# "Living"; portanto o EVENTO (óbito) é 1 - overall_survival.
audit_evento <- table(dados$overall_survival, dados$death_from_cancer, useNA = "ifany")
salvar_tab(as.data.frame.matrix(audit_evento), "T02_auditoria_codificacao_evento")
print(audit_evento)

dados$os_meses   <- as.numeric(dados$overall_survival_months)
dados$os_evento  <- ifelse(dados$overall_survival == 1, 0, 1)   # 1 = óbito por qualquer causa
dados$dss_evento <- ifelse(is.na(dados$death_from_cancer), NA_real_,
                    ifelse(dados$death_from_cancer == "Died of Disease", 1, 0))  # específico da doença
# (1 paciente sem informação de causa do óbito permanece NA e é excluído apenas das análises de DSS)

# --- 2.2 Subtipo molecular ---------------------------------------------------
dados$subtipo <- dados$`pam50_+_claudin-low_subtype`
dados$subtipo[dados$subtipo == "NC"] <- NA   # "NC" = não classificado

# --- 2.3 Critérios de exclusão (documentados) --------------------------------
n0 <- nrow(dados)
exc <- data.frame(etapa = character(), n_excluidos = integer(), n_restante = integer())
reg <- function(etapa, antes, depois) {
  exc <<- rbind(exc, data.frame(etapa = etapa, n_excluidos = antes - depois, n_restante = depois))
}

d1 <- dados[dados$cancer_type == "Breast Cancer", ]; reg("Histologia não mamária (Breast Sarcoma)", n0, nrow(d1))
d2 <- d1[!is.na(d1$subtipo), ];                      reg("Subtipo PAM50/claudin-low = NC", nrow(d1), nrow(d2))
d3 <- d2[!is.na(d2$os_meses) & d2$os_meses > 0, ];   reg("Tempo de seguimento ausente ou <= 0", nrow(d2), nrow(d3))
dt <- d3
salvar_tab(exc, "T03_fluxograma_exclusoes")
print(exc)

dt$subtipo <- factor(dt$subtipo,
                     levels = c("LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal"))

# --- 2.4 Covariáveis clínicas padronizadas -----------------------------------
dt$idade      <- as.numeric(dt$age_at_diagnosis)
dt$grau       <- factor(dt$neoplasm_histologic_grade, levels = c(1, 2, 3),
                        labels = c("G1", "G2", "G3"))
dt$tamanho    <- as.numeric(dt$tumor_size)
dt$linfonodos <- as.numeric(dt$lymph_nodes_examined_positive)
dt$ln_pos     <- factor(ifelse(dt$linfonodos > 0, "N+", "N0"), levels = c("N0", "N+"))
dt$npi        <- as.numeric(dt$nottingham_prognostic_index)
dt$er         <- factor(dt$er_status,   levels = c("Negative", "Positive"))
dt$pr         <- factor(dt$pr_status,   levels = c("Negative", "Positive"))
dt$her2       <- factor(dt$her2_status, levels = c("Negative", "Positive"))
dt$quimio     <- factor(ifelse(dt$chemotherapy   == 1, "Sim", "Não"), levels = c("Não", "Sim"))
dt$hormonio   <- factor(ifelse(dt$hormone_therapy == 1, "Sim", "Não"), levels = c("Não", "Sim"))
dt$radio      <- factor(ifelse(dt$radio_therapy   == 1, "Sim", "Não"), levels = c("Não", "Sim"))
dt$menopausa  <- factor(dt$inferred_menopausal_state)
dt$cluster_int<- factor(dt$integrative_cluster)

# --- 2.5 Matriz de expressão -------------------------------------------------
# Genes x amostras (formato exigido pelo limma). Valores já são z-scores
# calculados pelo cBioPortal em relação a todas as amostras do estudo.
EXPR <- t(as.matrix(dt[, col_expr]))
colnames(EXPR) <- dt$patient_id
mode(EXPR) <- "numeric"

qc_genes <- data.frame(
  gene    = rownames(EXPR),
  media   = rowMeans(EXPR),
  dp      = rowSds(EXPR),
  mad     = rowMads(EXPR),
  minimo  = rowMins(EXPR),
  maximo  = rowMaxs(EXPR)
)
salvar_tab(qc_genes[order(-qc_genes$dp), ], "T04_qc_expressao_por_gene")

# --- 2.6 Matriz de mutações (binarização) ------------------------------------
# Colunas *_mut trazem a alteração observada; "0" indica ausência de mutação.
MUT <- as.data.frame(lapply(dt[, col_mut], function(x) as.integer(!(x == "0" | is.na(x)))))
rownames(MUT) <- dt$patient_id
colnames(MUT) <- sub("_mut$", "", col_mut)

msg("MÓDULO 2 concluído — coorte analítica: ", nrow(dt), " pacientes | ",
    nrow(EXPR), " genes | ", ncol(MUT), " genes com status mutacional")

## =========================================================================
## MÓDULO 3 — CARACTERIZAÇÃO CLÍNICO-DESCRITIVA POR SUBTIPO
## =========================================================================
tab_subtipo <- dt %>%
  group_by(subtipo) %>%
  summarise(
    n                    = n(),
    perc                 = round(100 * n() / nrow(dt), 1),
    idade_mediana        = round(median(idade, na.rm = TRUE), 1),
    idade_iqr            = paste0(round(quantile(idade, .25, na.rm = TRUE), 1), "-",
                                  round(quantile(idade, .75, na.rm = TRUE), 1)),
    tamanho_mediano_mm   = round(median(tamanho, na.rm = TRUE), 1),
    perc_G3              = round(100 * mean(grau == "G3", na.rm = TRUE), 1),
    perc_N_positivo      = round(100 * mean(ln_pos == "N+", na.rm = TRUE), 1),
    perc_ER_positivo     = round(100 * mean(er == "Positive", na.rm = TRUE), 1),
    perc_HER2_positivo   = round(100 * mean(her2 == "Positive", na.rm = TRUE), 1),
    perc_quimioterapia   = round(100 * mean(quimio == "Sim"), 1),
    perc_hormonioterapia = round(100 * mean(hormonio == "Sim"), 1),
    npi_mediano          = round(median(npi, na.rm = TRUE), 2),
    seguimento_mediano   = round(median(os_meses), 1),
    obitos_totais        = sum(os_evento),
    perc_obitos          = round(100 * mean(os_evento), 1),
    obitos_por_cancer    = sum(dss_evento, na.rm = TRUE),
    .groups = "drop"
  )
salvar_tab(tab_subtipo, "T05_caracteristicas_clinicas_por_subtipo")
print(as.data.frame(tab_subtipo))

# Testes de associação subtipo x variáveis categóricas
vars_cat <- c("grau", "ln_pos", "er", "pr", "her2", "quimio", "hormonio", "radio", "menopausa")
testes_cat <- map_dfr(vars_cat, function(v) {
  tb <- table(dt[[v]], dt$subtipo)
  ct <- suppressWarnings(chisq.test(tb))
  data.frame(variavel = v, teste = "Qui-quadrado de Pearson",
             estatistica = round(unname(ct$statistic), 2), gl = unname(ct$parameter),
             valor_p = ct$p.value)
})
# Testes para variáveis contínuas (Kruskal-Wallis: não assume normalidade)
vars_num <- c("idade", "tamanho", "linfonodos", "npi", "os_meses")
testes_num <- map_dfr(vars_num, function(v) {
  kt <- kruskal.test(dt[[v]] ~ dt$subtipo)
  data.frame(variavel = v, teste = "Kruskal-Wallis",
             estatistica = round(unname(kt$statistic), 2), gl = unname(kt$parameter),
             valor_p = kt$p.value)
})
testes_clin <- rbind(testes_cat, testes_num)
testes_clin$valor_p_fdr <- p.adjust(testes_clin$valor_p, method = "BH")
salvar_tab(testes_clin, "T06_testes_associacao_clinica_subtipo")
print(testes_clin)

# Figura 1 — distribuição dos subtipos e mortalidade bruta
f1a <- ggplot(tab_subtipo, aes(x = fct_reorder(subtipo, -n), y = n, fill = subtipo)) +
  geom_col() + geom_text(aes(label = paste0(n, "\n(", perc, "%)")), vjust = -0.15, size = 3) +
  scale_fill_brewer(palette = "Set2") + expand_limits(y = max(tab_subtipo$n) * 1.18) +
  labs(title = "A. Distribuição dos subtipos moleculares", x = NULL, y = "Pacientes (n)") +
  theme_cowplot(11) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

f1b <- ggplot(tab_subtipo, aes(x = fct_reorder(subtipo, -perc_obitos), y = perc_obitos, fill = subtipo)) +
  geom_col() + geom_text(aes(label = paste0(perc_obitos, "%")), vjust = -0.3, size = 3) +
  scale_fill_brewer(palette = "Set2") + expand_limits(y = max(tab_subtipo$perc_obitos) * 1.15) +
  labs(title = "B. Mortalidade global bruta por subtipo", x = NULL, y = "Óbitos (%)") +
  theme_cowplot(11) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

salvar_fig(f1a + f1b, "F01_distribuicao_subtipos", w = 11, h = 5)
msg("MÓDULO 3 concluído — caracterização clínica exportada")

## =========================================================================
## MÓDULO 4 — EXPRESSÃO DIFERENCIAL POR SUBTIPO (limma)
## =========================================================================
# Estratégia: modelo linear sem intercepto (~0 + subtipo) e contrastes
# "um subtipo versus a média dos demais" (one-vs-rest), com moderação
# bayesiana empírica da variância (eBayes) e correção de Benjamini-Hochberg.
# Como a expressão já está em z-score, o coeficiente do contraste é a
# diferença média em desvios-padrão (interpretação direta do tamanho de efeito).

grupo  <- factor(make.names(as.character(dt$subtipo)))
design <- model.matrix(~ 0 + grupo)
colnames(design) <- levels(grupo)

fit <- lmFit(EXPR, design)

lvs <- levels(grupo)
contrastes <- sapply(lvs, function(g) {
  outros <- setdiff(lvs, g)
  paste0(g, " - (", paste(outros, collapse = " + "), ")/", length(outros))
})
cm  <- makeContrasts(contrasts = contrastes, levels = design)
colnames(cm) <- lvs
fit2 <- eBayes(contrasts.fit(fit, cm), trend = TRUE, robust = TRUE)

de_todos <- map_dfr(lvs, function(g) {
  tt <- topTable(fit2, coef = g, number = Inf, adjust.method = "BH", sort.by = "P")
  data.frame(subtipo = g, gene = rownames(tt), diff_z = tt$logFC, media_expr = tt$AveExpr,
             estatistica_t = tt$t, valor_p = tt$P.Value, fdr = tt$adj.P.Val, B = tt$B)
})
de_todos$direcao <- ifelse(de_todos$diff_z > 0, "Superexpresso", "Subexpresso")
de_todos$significativo <- de_todos$fdr < PAR$fdr_de & abs(de_todos$diff_z) >= PAR$lfc_de
salvar_tab(de_todos, "T07_expressao_diferencial_completa_todos_subtipos")

resumo_de <- de_todos %>% group_by(subtipo) %>%
  summarise(genes_testados = n(),
            signif_FDR = sum(fdr < PAR$fdr_de),
            signif_FDR_e_efeito = sum(significativo),
            super = sum(significativo & diff_z > 0),
            sub   = sum(significativo & diff_z < 0), .groups = "drop")
salvar_tab(resumo_de, "T08_resumo_expressao_diferencial")
print(as.data.frame(resumo_de))

top_de <- de_todos %>% filter(significativo) %>% group_by(subtipo) %>%
  slice_max(order_by = abs(diff_z), n = 15) %>% arrange(subtipo, -abs(diff_z)) %>% ungroup()
salvar_tab(top_de, "T09_top15_genes_por_subtipo")

# Figura 2 — volcano plots por subtipo
volc <- de_todos %>% mutate(neglog10p = -log10(valor_p))
rot <- volc %>% group_by(subtipo) %>% slice_max(order_by = abs(diff_z) * neglog10p, n = 8) %>% ungroup()
f2 <- ggplot(volc, aes(x = diff_z, y = neglog10p)) +
  geom_point(aes(color = ifelse(significativo, direcao, "Não significativo")), alpha = .6, size = 1) +
  geom_vline(xintercept = c(-PAR$lfc_de, PAR$lfc_de), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(max(volc$valor_p[volc$fdr < PAR$fdr_de])), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = rot, aes(label = toupper(gene)), size = 2.6, max.overlaps = 20, segment.alpha = .4) +
  scale_color_manual(values = c("Superexpresso" = "#D7263D", "Subexpresso" = "#1B98E0",
                                "Não significativo" = "grey75"), name = NULL) +
  facet_wrap(~ subtipo, scales = "free_y", ncol = 3) +
  labs(title = "Expressão diferencial por subtipo (contraste um-vs-demais, limma)",
       subtitle = paste0("Linhas tracejadas: |Δz| ≥ ", PAR$lfc_de, " e FDR < ", PAR$fdr_de),
       x = "Diferença média de expressão (Δ z-score)", y = expression(-log[10](p))) +
  theme_cowplot(10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey92"))
salvar_fig(f2, "F02_volcano_por_subtipo", w = 12, h = 8)

# Figura 3 — heatmap dos genes mais discriminantes (médias por subtipo)
genes_hm <- unique(de_todos %>% filter(significativo) %>% group_by(subtipo) %>%
                     slice_max(order_by = abs(diff_z), n = 10) %>% pull(gene))
mat_hm <- sapply(levels(dt$subtipo), function(s) rowMeans(EXPR[genes_hm, dt$subtipo == s, drop = FALSE]))
rownames(mat_hm) <- toupper(genes_hm)
png(file.path(DIR_FIG, "F03_heatmap_genes_discriminantes.png"), width = 1500, height = 2100, res = 170)
pheatmap(mat_hm, cluster_cols = TRUE, cluster_rows = TRUE, fontsize_row = 7, fontsize_col = 10,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), border_color = NA,
         main = "Expressão média (z-score) por subtipo — top 10 genes/subtipo")
dev.off()
msg("MÓDULO 4 concluído — ", sum(de_todos$significativo), " associações gene-subtipo significativas")

## =========================================================================
## MÓDULO 5 — CLASSIFICAÇÃO SUPERVISIONADA E IMPORTÂNCIA (RANDOM FOREST)
## =========================================================================
# Objetivo: ranquear genes por poder discriminante multivariado e não-linear,
# de forma complementar ao limma (que é univariado por gene).
X_rf <- as.data.frame(t(EXPR))
y_rf <- droplevels(dt$subtipo)
set.seed(42)
rf <- randomForest(x = X_rf, y = y_rf, ntree = PAR$n_arvores, importance = TRUE)

rf_conf <- as.data.frame.matrix(rf$confusion)
salvar_tab(cbind(subtipo_real = rownames(rf_conf), rf_conf), "T10_randomforest_matriz_confusao")

imp <- as.data.frame(importance(rf))
imp$gene <- rownames(imp)
imp_rank <- imp[order(-imp$MeanDecreaseGini), c("gene", "MeanDecreaseAccuracy", "MeanDecreaseGini")]
imp_rank$posicao <- seq_len(nrow(imp_rank))
salvar_tab(imp_rank, "T11_randomforest_importancia_genes")

erro_oob <- round(100 * rf$err.rate[PAR$n_arvores, "OOB"], 2)
msg("MÓDULO 5 — erro OOB do Random Forest: ", erro_oob, "%")

f4 <- ggplot(head(imp_rank, 30), aes(x = MeanDecreaseGini, y = fct_reorder(toupper(gene), MeanDecreaseGini))) +
  geom_segment(aes(xend = 0, yend = fct_reorder(toupper(gene), MeanDecreaseGini)), color = "grey70") +
  geom_point(aes(color = MeanDecreaseAccuracy), size = 3) +
  scale_color_viridis_c(option = "C", name = "Queda média\nde acurácia") +
  labs(title = "Top 30 genes discriminantes de subtipo (Random Forest)",
       subtitle = paste0("Erro out-of-bag = ", erro_oob, "% | ", PAR$n_arvores, " árvores"),
       x = "Queda média do índice de Gini", y = NULL) +
  theme_cowplot(10)
salvar_fig(f4, "F04_importancia_randomforest", w = 8, h = 8)

## =========================================================================
## MÓDULO 6 — SOBREVIDA UNIVARIADA GENE A GENE (COX)
## =========================================================================
# Modelo de riscos proporcionais de Cox para cada gene, com a expressão em
# z-score: HR representa o risco por aumento de 1 desvio-padrão na expressão.
# Dois desfechos: sobrevida global (OS) e sobrevida específica por câncer (DSS).

cox_univar <- function(tempo, evento, rotulo) {
  surv_obj <- Surv(tempo, evento)
  map_dfr(rownames(EXPR), function(g) {
    m <- coxph(surv_obj ~ EXPR[g, ])
    s <- summary(m)
    data.frame(desfecho = rotulo, gene = g,
               HR = unname(s$coefficients[1, "exp(coef)"]),
               IC95_inf = unname(s$conf.int[1, "lower .95"]),
               IC95_sup = unname(s$conf.int[1, "upper .95"]),
               z = unname(s$coefficients[1, "z"]),
               valor_p = unname(s$coefficients[1, "Pr(>|z|)"]),
               concordancia = unname(s$concordance[1]))
  })
}

cox_os  <- cox_univar(dt$os_meses, dt$os_evento,  "Sobrevida global (OS)")
cox_dss <- cox_univar(dt$os_meses, dt$dss_evento, "Sobrevida específica (DSS)")
cox_os$fdr  <- p.adjust(cox_os$valor_p,  method = "BH")
cox_dss$fdr <- p.adjust(cox_dss$valor_p, method = "BH")
cox_os  <- cox_os[order(cox_os$valor_p), ]
cox_dss <- cox_dss[order(cox_dss$valor_p), ]
salvar_tab(cox_os,  "T12_cox_univariado_OS_todos_genes")
salvar_tab(cox_dss, "T13_cox_univariado_DSS_todos_genes")

msg("MÓDULO 6 — genes com FDR<", PAR$fdr_cox, ": OS = ", sum(cox_os$fdr < PAR$fdr_cox),
    " | DSS = ", sum(cox_dss$fdr < PAR$fdr_cox))

# Verificação do pressuposto de riscos proporcionais para os 20 principais genes
top20_os <- head(cox_os$gene, 20)
ph_test <- map_dfr(top20_os, function(g) {
  m <- coxph(Surv(dt$os_meses, dt$os_evento) ~ EXPR[g, ])
  z <- cox.zph(m)
  data.frame(gene = g, chisq_zph = unname(z$table[1, "chisq"]), p_zph = unname(z$table[1, "p"]),
             pressuposto = ifelse(z$table[1, "p"] > 0.05, "Atendido", "Violado"))
})
salvar_tab(ph_test, "T14_teste_riscos_proporcionais_top20")

# Modelo de Cox ajustado por covariáveis clínicas para os 20 principais genes
clin_df <- data.frame(tempo = dt$os_meses, evento = dt$os_evento, idade = dt$idade,
                      grau = dt$grau, tamanho = dt$tamanho, ln_pos = dt$ln_pos,
                      subtipo = dt$subtipo, quimio = dt$quimio, hormonio = dt$hormonio,
                      radio = dt$radio)
cox_ajustado <- map_dfr(top20_os, function(g) {
  df <- cbind(clin_df, gene_z = EXPR[g, ])
  m <- coxph(Surv(tempo, evento) ~ gene_z + idade + grau + tamanho + ln_pos + subtipo +
               quimio + hormonio + radio, data = df)
  s <- summary(m)
  data.frame(gene = g, n_modelo = s$n, eventos = s$nevent,
             HR_ajustado = unname(s$coefficients["gene_z", "exp(coef)"]),
             IC95_inf = unname(s$conf.int["gene_z", "lower .95"]),
             IC95_sup = unname(s$conf.int["gene_z", "upper .95"]),
             valor_p = unname(s$coefficients["gene_z", "Pr(>|z|)"]),
             C_index_modelo = unname(s$concordance[1]))
})
cox_ajustado$fdr <- p.adjust(cox_ajustado$valor_p, method = "BH")
cox_ajustado <- cox_ajustado[order(cox_ajustado$valor_p), ]
salvar_tab(cox_ajustado, "T15_cox_multivariado_ajustado_top20")

# Figura 5 — forest plot dos 25 genes prognósticos mais fortes (OS)
fp <- head(cox_os, 25) %>% mutate(gene = toupper(gene),
                                  efeito = ifelse(HR > 1, "Risco aumentado", "Risco reduzido"))
f5 <- ggplot(fp, aes(x = HR, y = fct_reorder(gene, HR), color = efeito)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = IC95_inf, xmax = IC95_sup), height = .25) +
  geom_point(size = 2.4) +
  scale_x_continuous(trans = "log2", breaks = c(0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5)) +
  scale_color_manual(values = c("Risco aumentado" = "#D7263D", "Risco reduzido" = "#1B98E0"), name = NULL) +
  labs(title = "Genes com maior associação à sobrevida global (Cox univariado)",
       subtitle = "HR por aumento de 1 desvio-padrão (z-score) na expressão; barras = IC 95%",
       x = "Hazard Ratio (escala log2)", y = NULL) +
  theme_cowplot(10) + theme(legend.position = "bottom")
salvar_fig(f5, "F05_forest_cox_univariado_OS", w = 8.5, h = 8)

## =========================================================================
## MÓDULO 7 — KAPLAN-MEIER (SUBTIPOS E GENES DE MAIOR EFEITO)
## =========================================================================
km_sub <- survfit(Surv(os_meses, os_evento) ~ subtipo, data = dt)
lr_sub <- survdiff(Surv(os_meses, os_evento) ~ subtipo, data = dt)
p_lr_sub <- pchisq(lr_sub$chisq, df = length(lr_sub$n) - 1, lower.tail = FALSE)

mediana_sub <- as.data.frame(summary(km_sub)$table)
mediana_sub$subtipo <- sub("subtipo=", "", rownames(mediana_sub))
salvar_tab(mediana_sub[, c("subtipo", "records", "events", "median", "0.95LCL", "0.95UCL")],
           "T16_mediana_sobrevida_por_subtipo")

g_km_sub <- ggsurvplot(km_sub, data = dt, pval = TRUE, risk.table = TRUE, conf.int = FALSE,
                       xlim = c(0, 300), break.time.by = 60, palette = "Set2",
                       legend.title = "Subtipo", legend.labs = levels(dt$subtipo),
                       xlab = "Tempo (meses)", ylab = "Probabilidade de sobrevida global",
                       title = "Sobrevida global por subtipo molecular (METABRIC)",
                       risk.table.height = .28, ggtheme = theme_cowplot(10))
png(file.path(DIR_FIG, "F06_km_subtipos.png"), width = 1700, height = 1500, res = 170)
print(g_km_sub); dev.off()

# KM para os 4 genes prognósticos mais fortes (dicotomizados pela mediana)
genes_km <- head(cox_os$gene, 4)
dt_km <- dt
for (g in genes_km) {
  dt_km[[paste0("grp_", g)]] <- factor(ifelse(EXPR[g, ] > median(EXPR[g, ]), "Alta", "Baixa"),
                                       levels = c("Baixa", "Alta"))
}
km_genes <- map_dfr(genes_km, function(g) {
  f <- as.formula(paste0("Surv(os_meses, os_evento) ~ grp_", g))
  sd <- survdiff(f, data = dt_km)
  fit <- survfit(f, data = dt_km)
  med <- summary(fit)$table[, "median"]
  data.frame(gene = g, mediana_baixa = med[1], mediana_alta = med[2],
             qui_quadrado_logrank = unname(sd$chisq),
             valor_p_logrank = pchisq(sd$chisq, df = 1, lower.tail = FALSE))
})
salvar_tab(km_genes, "T17_logrank_top4_genes")

plots_km <- lapply(genes_km, function(g) {
  f <- as.formula(paste0("Surv(os_meses, os_evento) ~ grp_", g))
  ggsurvplot(surv_fit(f, data = dt_km), data = dt_km, pval = TRUE, conf.int = TRUE,
             xlim = c(0, 300), break.time.by = 100, palette = c("#1B98E0", "#D7263D"),
             legend.title = toupper(g), legend.labs = c("Expressão baixa", "Expressão alta"),
             xlab = "Tempo (meses)", ylab = "Sobrevida global", ggtheme = theme_cowplot(9))$plot
})
salvar_fig(wrap_plots(plots_km, ncol = 2) +
             plot_annotation(title = "Sobrevida global segundo expressão gênica (corte pela mediana)"),
           "F07_km_top4_genes", w = 11, h = 8)

# Análise estratificada: efeito prognóstico dos top 12 genes DENTRO de cada subtipo
top12 <- head(cox_os$gene, 12)
cox_por_subtipo <- map_dfr(levels(dt$subtipo), function(s) {
  idx <- which(dt$subtipo == s)
  if (sum(dt$os_evento[idx]) < 15) return(NULL)
  map_dfr(top12, function(g) {
    m <- try(coxph(Surv(dt$os_meses[idx], dt$os_evento[idx]) ~ EXPR[g, idx]), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    s2 <- summary(m)
    data.frame(subtipo = s, gene = g, n = length(idx), eventos = sum(dt$os_evento[idx]),
               HR = unname(s2$coefficients[1, "exp(coef)"]),
               valor_p = unname(s2$coefficients[1, "Pr(>|z|)"]))
  })
})
cox_por_subtipo$fdr <- p.adjust(cox_por_subtipo$valor_p, method = "BH")
salvar_tab(cox_por_subtipo, "T18_cox_estratificado_por_subtipo")

f8 <- ggplot(cox_por_subtipo, aes(x = subtipo, y = toupper(gene), fill = log2(HR))) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(valor_p < 0.001, "***", ifelse(valor_p < 0.01, "**",
                        ifelse(valor_p < 0.05, "*", "")))), size = 4, vjust = .75) +
  scale_fill_gradient2(low = "#1B98E0", mid = "white", high = "#D7263D", midpoint = 0,
                       name = "log2(HR)") +
  labs(title = "Efeito prognóstico dos principais genes dentro de cada subtipo",
       subtitle = "Cox univariado por estrato; * p<0,05  ** p<0,01  *** p<0,001",
       x = NULL, y = NULL) +
  theme_cowplot(10) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
salvar_fig(f8, "F08_heatmap_cox_estratificado", w = 8, h = 7)
msg("MÓDULO 7 concluído — log-rank global entre subtipos: p = ", format.pval(p_lr_sub, digits = 3))

## =========================================================================
## MÓDULO 8 — ASSINATURA PROGNÓSTICA MULTIGÊNICA (LASSO-COX)
## =========================================================================
# Partição estratificada (70% derivação / 30% validação) preservando a
# distribuição conjunta de subtipo e evento.
set.seed(42)
estrato <- interaction(dt$subtipo, dt$os_evento, drop = TRUE)
idx_treino <- unlist(lapply(split(seq_len(nrow(dt)), estrato), function(i)
  sample(i, size = max(1, floor(PAR$prop_treino * length(i))))))
idx_treino <- sort(idx_treino)
idx_teste  <- setdiff(seq_len(nrow(dt)), idx_treino)

X <- t(EXPR)
y_tr <- Surv(dt$os_meses[idx_treino], dt$os_evento[idx_treino])
y_te <- Surv(dt$os_meses[idx_teste],  dt$os_evento[idx_teste])

set.seed(42)
cvfit <- cv.glmnet(X[idx_treino, ], y_tr, family = "cox", alpha = 1,
                   nfolds = PAR$n_folds_cv, standardize = TRUE)
coefs <- coef(cvfit, s = "lambda.min")
assinatura <- data.frame(gene = rownames(coefs)[which(coefs != 0)],
                         coeficiente = as.numeric(coefs[which(coefs != 0)]))
assinatura$efeito <- ifelse(assinatura$coeficiente > 0, "Risco (↑ expressão = pior)",
                            "Proteção (↑ expressão = melhor)")
assinatura$HR_por_DP <- exp(assinatura$coeficiente)
assinatura <- assinatura[order(-abs(assinatura$coeficiente)), ]
salvar_tab(assinatura, "T19_assinatura_lasso_cox_coeficientes")

# Sensibilidade: modelo mais parcimonioso (lambda.1se)
coefs_1se <- coef(cvfit, s = "lambda.1se")
assinatura_1se <- data.frame(gene = rownames(coefs_1se)[which(coefs_1se != 0)],
                             coeficiente = as.numeric(coefs_1se[which(coefs_1se != 0)]))
assinatura_1se$HR_por_DP <- exp(assinatura_1se$coeficiente)
assinatura_1se <- assinatura_1se[order(-abs(assinatura_1se$coeficiente)), ]
salvar_tab(assinatura_1se, "T19b_assinatura_lasso_lambda1se")
escore_1se <- as.numeric(predict(cvfit, newx = X, s = "lambda.1se", type = "link"))

escore <- as.numeric(predict(cvfit, newx = X, s = "lambda.min", type = "link"))
dt$escore_risco <- escore

# Desempenho discriminativo (C-index) — gene, clínico e combinado
d_tr <- data.frame(tempo = dt$os_meses[idx_treino], evento = dt$os_evento[idx_treino],
                   escore = escore[idx_treino], idade = dt$idade[idx_treino],
                   grau = dt$grau[idx_treino], tamanho = dt$tamanho[idx_treino],
                   ln_pos = dt$ln_pos[idx_treino], subtipo = dt$subtipo[idx_treino])
d_te <- data.frame(tempo = dt$os_meses[idx_teste], evento = dt$os_evento[idx_teste],
                   escore = escore[idx_teste], idade = dt$idade[idx_teste],
                   grau = dt$grau[idx_teste], tamanho = dt$tamanho[idx_teste],
                   ln_pos = dt$ln_pos[idx_teste], subtipo = dt$subtipo[idx_teste])

m_gene  <- coxph(Surv(tempo, evento) ~ escore, data = d_tr)
m_clin  <- coxph(Surv(tempo, evento) ~ idade + grau + tamanho + ln_pos + subtipo, data = d_tr)
m_comb  <- coxph(Surv(tempo, evento) ~ escore + idade + grau + tamanho + ln_pos + subtipo, data = d_tr)

c_index <- function(modelo, dados) {
  lp <- predict(modelo, newdata = dados, type = "lp")
  ok <- complete.cases(lp, dados$tempo, dados$evento)
  unname(concordance(Surv(dados$tempo[ok], dados$evento[ok]) ~ lp[ok], reverse = TRUE)$concordance)
}
d_tr$escore_1se <- escore_1se[idx_treino]; d_te$escore_1se <- escore_1se[idx_teste]
m_gene1se <- coxph(Surv(tempo, evento) ~ escore_1se, data = d_tr)
m_comb1se <- coxph(Surv(tempo, evento) ~ escore_1se + idade + grau + tamanho + ln_pos + subtipo, data = d_tr)
desempenho <- data.frame(
  modelo = c("Assinatura gênica LASSO (lambda.min)", "Assinatura gênica LASSO (lambda.1se)",
             "Clínico (idade, grau, tamanho, linfonodo, subtipo)",
             "Combinado (clínico + assinatura lambda.min)",
             "Combinado (clínico + assinatura lambda.1se)"),
  n_genes = c(nrow(assinatura), nrow(assinatura_1se), 0, nrow(assinatura), nrow(assinatura_1se)),
  C_index_treino = c(c_index(m_gene, d_tr), c_index(m_gene1se, d_tr), c_index(m_clin, d_tr),
                     c_index(m_comb, d_tr), c_index(m_comb1se, d_tr)),
  C_index_validacao = c(c_index(m_gene, d_te), c_index(m_gene1se, d_te), c_index(m_clin, d_te),
                        c_index(m_comb, d_te), c_index(m_comb1se, d_te))
)
desempenho[, 3:4] <- round(desempenho[, 3:4], 4)
salvar_tab(desempenho, "T20_desempenho_modelos_C_index")
print(desempenho)

# Teste de verossimilhança: a assinatura acrescenta informação ao modelo clínico?
lrt <- anova(m_clin, m_comb, test = "LRT")
salvar_tab(data.frame(comparacao = "Clínico vs Clínico+Assinatura",
                      qui_quadrado = round(lrt$Chisq[2], 3), gl = lrt$Df[2],
                      valor_p = lrt$`Pr(>|Chi|)`[2]), "T21_teste_razao_verossimilhanca")

# Estratificação de risco (tercis definidos no conjunto de derivação)
cortes <- quantile(escore[idx_treino], probs = c(1/3, 2/3))
dt$grupo_risco <- cut(escore, breaks = c(-Inf, cortes, Inf),
                      labels = c("Baixo risco", "Risco intermediário", "Alto risco"))
dt$conjunto <- ifelse(seq_len(nrow(dt)) %in% idx_treino, "Derivação", "Validação")

km_risco_val <- survfit(Surv(os_meses, os_evento) ~ grupo_risco, data = dt[dt$conjunto == "Validação", ])
sd_val <- survdiff(Surv(os_meses, os_evento) ~ grupo_risco, data = dt[dt$conjunto == "Validação", ])
p_val <- pchisq(sd_val$chisq, df = 2, lower.tail = FALSE)

tab_risco <- dt %>% group_by(conjunto, grupo_risco) %>%
  summarise(n = n(), obitos = sum(os_evento), perc_obitos = round(100 * mean(os_evento), 1),
            .groups = "drop")
salvar_tab(tab_risco, "T22_grupos_de_risco_por_conjunto")

g_km_risco <- ggsurvplot(km_risco_val, data = dt[dt$conjunto == "Validação", ], pval = TRUE,
                         risk.table = TRUE, conf.int = TRUE, xlim = c(0, 300), break.time.by = 60,
                         palette = c("#2E933C", "#F2A65A", "#D7263D"),
                         legend.title = "Grupo", legend.labs = levels(dt$grupo_risco),
                         xlab = "Tempo (meses)", ylab = "Sobrevida global",
                         title = "Validação independente da assinatura LASSO-Cox (30% da coorte)",
                         risk.table.height = .28, ggtheme = theme_cowplot(10))
png(file.path(DIR_FIG, "F09_km_grupos_risco_validacao.png"), width = 1700, height = 1500, res = 170)
print(g_km_risco); dev.off()

f10 <- ggplot(assinatura, aes(x = coeficiente, y = fct_reorder(toupper(gene), coeficiente),
                              fill = coeficiente > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "#D7263D", "FALSE" = "#1B98E0"),
                    labels = c("TRUE" = "Risco", "FALSE" = "Proteção"), name = NULL) +
  labs(title = "Assinatura prognóstica selecionada por LASSO-Cox",
       subtitle = paste0(nrow(assinatura), " genes com coeficiente não nulo em lambda.min"),
       x = "Coeficiente de Cox penalizado", y = NULL) +
  theme_cowplot(9)
salvar_fig(f10, "F10_coeficientes_assinatura_lasso", w = 8, h = max(6, nrow(assinatura) * 0.16))
msg("MÓDULO 8 concluído — ", nrow(assinatura), " genes na assinatura; C-index validação = ",
    desempenho$C_index_validacao[1])

## =========================================================================
## MÓDULO 8b — ESTABILIDADE E VALIDAÇÃO REPETIDA (25 PARTIÇÕES ALEATÓRIAS)
## =========================================================================
# Uma única partição 70/30 produz estimativas instáveis de C-index. Repetimos
# todo o procedimento (partição estratificada -> cv.glmnet -> avaliação) 25
# vezes, registrando (i) a distribuição do C-index dos três modelos e (ii) a
# frequência com que cada gene é selecionado pelo LASSO (stability selection).
N_REP <- as.integer(Sys.getenv("METABRIC_NREP", unset = "25"))
sel_freq <- setNames(integer(nrow(EXPR)), rownames(EXPR))
rep_perf <- vector("list", N_REP)

for (r in seq_len(N_REP)) {
  set.seed(1000 + r)
  itr <- sort(unlist(lapply(split(seq_len(nrow(dt)), estrato), function(i)
    sample(i, size = max(1, floor(PAR$prop_treino * length(i)))))))
  ite <- setdiff(seq_len(nrow(dt)), itr)
  cvr <- cv.glmnet(X[itr, ], Surv(dt$os_meses[itr], dt$os_evento[itr]),
                   family = "cox", alpha = 1, nfolds = PAR$n_folds_cv, standardize = TRUE)
  cf  <- coef(cvr, s = "lambda.1se")
  gsel <- rownames(cf)[which(cf != 0)]
  sel_freq[gsel] <- sel_freq[gsel] + 1L
  esc <- as.numeric(predict(cvr, newx = X, s = "lambda.1se", type = "link"))
  dtr <- data.frame(tempo = dt$os_meses[itr], evento = dt$os_evento[itr], escore = esc[itr],
                    idade = dt$idade[itr], grau = dt$grau[itr], tamanho = dt$tamanho[itr],
                    ln_pos = dt$ln_pos[itr], subtipo = dt$subtipo[itr])
  dte <- data.frame(tempo = dt$os_meses[ite], evento = dt$os_evento[ite], escore = esc[ite],
                    idade = dt$idade[ite], grau = dt$grau[ite], tamanho = dt$tamanho[ite],
                    ln_pos = dt$ln_pos[ite], subtipo = dt$subtipo[ite])
  mg <- coxph(Surv(tempo, evento) ~ escore, data = dtr)
  mc <- coxph(Surv(tempo, evento) ~ idade + grau + tamanho + ln_pos + subtipo, data = dtr)
  mk <- coxph(Surv(tempo, evento) ~ escore + idade + grau + tamanho + ln_pos + subtipo, data = dtr)
  rep_perf[[r]] <- data.frame(repeticao = r, n_genes = length(gsel),
                              C_gene = c_index(mg, dte), C_clinico = c_index(mc, dte),
                              C_combinado = c_index(mk, dte))
}
rep_perf <- bind_rows(rep_perf)
salvar_tab(rep_perf, "T29_validacao_repetida_25_particoes")

resumo_rep <- data.frame(
  modelo = c("Assinatura gênica (LASSO lambda.1se)", "Clínico", "Combinado"),
  C_index_medio = round(c(mean(rep_perf$C_gene), mean(rep_perf$C_clinico), mean(rep_perf$C_combinado)), 4),
  desvio_padrao = round(c(sd(rep_perf$C_gene), sd(rep_perf$C_clinico), sd(rep_perf$C_combinado)), 4),
  minimo = round(c(min(rep_perf$C_gene), min(rep_perf$C_clinico), min(rep_perf$C_combinado)), 4),
  maximo = round(c(max(rep_perf$C_gene), max(rep_perf$C_clinico), max(rep_perf$C_combinado)), 4)
)
salvar_tab(resumo_rep, "T30_resumo_validacao_repetida")
print(resumo_rep)

teste_pareado <- data.frame(
  comparacao = c("Combinado vs Clínico", "Combinado vs Assinatura", "Assinatura vs Clínico"),
  diferenca_media = round(c(mean(rep_perf$C_combinado - rep_perf$C_clinico),
                            mean(rep_perf$C_combinado - rep_perf$C_gene),
                            mean(rep_perf$C_gene - rep_perf$C_clinico)), 4),
  valor_p_wilcoxon = c(wilcox.test(rep_perf$C_combinado, rep_perf$C_clinico, paired = TRUE)$p.value,
                       wilcox.test(rep_perf$C_combinado, rep_perf$C_gene, paired = TRUE)$p.value,
                       wilcox.test(rep_perf$C_gene, rep_perf$C_clinico, paired = TRUE)$p.value)
)
salvar_tab(teste_pareado, "T31_comparacao_pareada_modelos")

estab <- data.frame(gene = names(sel_freq), vezes_selecionado = as.integer(sel_freq),
                    freq_selecao_perc = round(100 * as.numeric(sel_freq) / N_REP, 1))
estab <- estab[order(-estab$vezes_selecionado), ]
estab <- left_join(estab, cox_os[, c("gene", "HR", "fdr")], by = "gene")
salvar_tab(estab, "T32_estabilidade_selecao_lasso")

f13a <- rep_perf %>% select(-repeticao, -n_genes) %>%
  pivot_longer(everything(), names_to = "modelo", values_to = "C") %>%
  mutate(modelo = recode(modelo, C_gene = "Assinatura", C_clinico = "Clínico", C_combinado = "Combinado")) %>%
  ggplot(aes(x = modelo, y = C, fill = modelo)) +
  geom_boxplot(alpha = .8, outlier.size = .8) + geom_jitter(width = .12, size = .8, alpha = .5) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = paste0("A. C-index em ", N_REP, " partições independentes"), x = NULL, y = "C-index (conjunto de validação)") +
  theme_cowplot(10) + theme(legend.position = "none")

f13b <- estab %>% filter(vezes_selecionado > 0) %>% head(25) %>%
  ggplot(aes(x = freq_selecao_perc, y = fct_reorder(toupper(gene), freq_selecao_perc),
             fill = ifelse(HR > 1, "HR > 1", "HR < 1"))) +
  geom_col() + scale_fill_manual(values = c("HR > 1" = "#D7263D", "HR < 1" = "#1B98E0"), name = NULL) +
  labs(title = paste0("B. Estabilidade de seleção (LASSO, ", N_REP, " repetições)"),
       x = "Frequência de seleção (%)", y = NULL) +
  theme_cowplot(9)
salvar_fig(f13a + f13b + plot_layout(widths = c(1, 1.4)), "F13_estabilidade_validacao_repetida", w = 12, h = 7)
msg("MÓDULO 8b concluído — C-index médio (", N_REP, " partições): assinatura = ",
    resumo_rep$C_index_medio[1], " | clínico = ", resumo_rep$C_index_medio[2],
    " | combinado = ", resumo_rep$C_index_medio[3])

## =========================================================================
## MÓDULO 9 — MUTAÇÕES SOMÁTICAS POR SUBTIPO E PROGNÓSTICO
## =========================================================================
freq_mut <- colMeans(MUT)
genes_mut <- names(freq_mut)[freq_mut >= PAR$freq_min_mut]

mut_subtipo <- map_dfr(genes_mut, function(g) {
  tb <- table(MUT[[g]], dt$subtipo)
  ft <- fisher.test(tb, simulate.p.value = TRUE, B = 20000)
  prop <- round(100 * tapply(MUT[[g]], dt$subtipo, mean), 1)
  data.frame(gene = g, freq_global_perc = round(100 * freq_mut[g], 1),
              t(prop), valor_p = ft$p.value)
})
mut_subtipo$fdr <- p.adjust(mut_subtipo$valor_p, method = "BH")
mut_subtipo <- mut_subtipo[order(mut_subtipo$valor_p), ]
salvar_tab(mut_subtipo, "T23_mutacoes_por_subtipo")

cox_mut <- map_dfr(genes_mut, function(g) {
  m <- coxph(Surv(dt$os_meses, dt$os_evento) ~ MUT[[g]])
  s <- summary(m)
  data.frame(gene = g, freq_perc = round(100 * freq_mut[g], 1),
             HR = unname(s$coefficients[1, "exp(coef)"]),
             IC95_inf = unname(s$conf.int[1, "lower .95"]),
             IC95_sup = unname(s$conf.int[1, "upper .95"]),
             valor_p = unname(s$coefficients[1, "Pr(>|z|)"]))
})
cox_mut$fdr <- p.adjust(cox_mut$valor_p, method = "BH")
cox_mut <- cox_mut[order(cox_mut$valor_p), ]
salvar_tab(cox_mut, "T24_cox_mutacoes_OS")

mut_long <- mut_subtipo %>% head(20) %>%
  select(gene, all_of(make.names(levels(dt$subtipo)))) %>%
  pivot_longer(-gene, names_to = "subtipo", values_to = "freq")
f11 <- ggplot(mut_long, aes(x = factor(subtipo, levels = make.names(levels(dt$subtipo))),
                            y = fct_reorder(toupper(gene), freq), fill = freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0(freq, "%")), size = 2.6) +
  scale_fill_viridis_c(option = "magma", direction = -1, name = "Frequência (%)") +
  labs(title = "Frequência de mutação somática por subtipo molecular",
       subtitle = "20 genes com maior heterogeneidade entre subtipos (teste exato de Fisher)",
       x = NULL, y = NULL) +
  theme_cowplot(10) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
salvar_fig(f11, "F11_mutacoes_por_subtipo", w = 8.5, h = 7)
msg("MÓDULO 9 concluído — ", length(genes_mut), " genes com frequência >= ",
    100 * PAR$freq_min_mut, "% testados")

## =========================================================================
## MÓDULO 10 — INTEGRAÇÃO FINAL E EXPORTAÇÃO
## =========================================================================
# Cruzamento: genes que (i) marcam um subtipo específico e (ii) carregam
# informação prognóstica independente na coorte completa.
integra <- de_todos %>%
  filter(significativo) %>%
  select(subtipo, gene, diff_z, fdr_subtipo = fdr) %>%
  inner_join(cox_os %>% select(gene, HR_OS = HR, IC95_inf, IC95_sup, fdr_OS = fdr), by = "gene") %>%
  filter(fdr_OS < PAR$fdr_cox) %>%
  mutate(
    direcao_subtipo = ifelse(diff_z > 0, "Superexpresso", "Subexpresso"),
    direcao_risco   = ifelse(HR_OS > 1, "Maior risco de óbito", "Menor risco de óbito"),
    coerencia = ifelse((diff_z > 0 & HR_OS > 1) | (diff_z < 0 & HR_OS < 1),
                       "Marcador de pior prognóstico no subtipo",
                       "Marcador de melhor prognóstico no subtipo"),
    escore_prioridade = abs(diff_z) * abs(log(HR_OS)) * -log10(fdr_OS)
  ) %>%
  left_join(imp_rank %>% select(gene, posicao_RF = posicao), by = "gene") %>%
  left_join(assinatura %>% select(gene, coef_LASSO = coeficiente), by = "gene") %>%
  arrange(subtipo, desc(escore_prioridade))
salvar_tab(integra, "T25_integracao_subtipo_prognostico")

prioridade <- integra %>% group_by(subtipo) %>% slice_max(escore_prioridade, n = 10) %>% ungroup()
salvar_tab(prioridade, "T26_genes_prioritarios_por_subtipo")

f12 <- ggplot(integra, aes(x = diff_z, y = log2(HR_OS))) +
  geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
  geom_point(aes(color = subtipo, size = -log10(fdr_OS)), alpha = .7) +
  geom_text_repel(data = prioridade %>% group_by(subtipo) %>% slice_max(escore_prioridade, n = 4),
                  aes(label = toupper(gene)), size = 2.5, max.overlaps = 25) +
  scale_color_brewer(palette = "Set2", name = "Subtipo") +
  scale_size_continuous(name = expression(-log[10](FDR[OS]))) +
  facet_wrap(~ subtipo, ncol = 3) +
  labs(title = "Integração: especificidade de subtipo × efeito prognóstico",
       subtitle = "Quadrantes superiores-direito e inferiores-esquerdo indicam coerência entre superexpressão e risco",
       x = "Δ z-score (subtipo vs demais)", y = "log2(HR) para sobrevida global") +
  theme_cowplot(10) + theme(legend.position = "bottom", strip.background = element_rect(fill = "grey92"))
salvar_fig(f12, "F12_integracao_subtipo_prognostico", w = 12, h = 8)

# --- Metadados da execução ---------------------------------------------------
resumo_exec <- data.frame(
  parametro = c("Data/hora da execução", "Versão do R", "Pacientes na coorte analítica",
                "Genes de expressão analisados", "Genes com dado mutacional testados",
                "Óbitos por qualquer causa", "Óbitos por câncer de mama",
                "Seguimento mediano (meses)", "Associações gene-subtipo significativas",
                "Genes prognósticos (FDR<0,05, OS)", "Genes na assinatura LASSO-Cox",
                "C-index validação (assinatura)", "C-index validação (clínico)",
                "C-index validação (combinado)", paste0("C-index médio ", N_REP, " partições (assinatura)"),
                paste0("C-index médio ", N_REP, " partições (clínico)"), paste0("C-index médio ", N_REP, " partições (combinado)"),
                "Erro OOB Random Forest (%)",
                "Tempo total de execução (min)"),
  valor = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R.version.string, nrow(dt),
            nrow(EXPR), length(genes_mut), sum(dt$os_evento), sum(dt$dss_evento, na.rm = TRUE),
            round(median(dt$os_meses), 1), sum(de_todos$significativo),
            sum(cox_os$fdr < PAR$fdr_cox), nrow(assinatura),
            desempenho$C_index_validacao[1], desempenho$C_index_validacao[3],
            desempenho$C_index_validacao[4], resumo_rep$C_index_medio[1],
            resumo_rep$C_index_medio[2], resumo_rep$C_index_medio[3], erro_oob,
            round(as.numeric(difftime(Sys.time(), t_inicio, units = "mins")), 2))
)
salvar_tab(resumo_exec, "T27_resumo_execucao")
print(resumo_exec)

bibs <- data.frame(
  pacote = c("data.table", "dplyr", "tidyr", "purrr", "stringr", "forcats", "limma",
             "survival", "survminer", "glmnet", "randomForest", "matrixStats", "broom",
             "ggplot2", "ggrepel", "scales", "RColorBrewer", "viridis", "pheatmap",
             "patchwork", "cowplot", "gridExtra"),
  versao = sapply(c("data.table", "dplyr", "tidyr", "purrr", "stringr", "forcats", "limma",
                    "survival", "survminer", "glmnet", "randomForest", "matrixStats", "broom",
                    "ggplot2", "ggrepel", "scales", "RColorBrewer", "viridis", "pheatmap",
                    "patchwork", "cowplot", "gridExtra"),
                  function(p) as.character(packageVersion(p)))
)
salvar_tab(bibs, "T28_bibliotecas_versoes")

writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "sessionInfo.txt"))
saveRDS(list(dt = dt, EXPR = EXPR, de = de_todos, cox_os = cox_os, cox_dss = cox_dss,
             assinatura = assinatura, integra = integra, estabilidade = estab,
             desempenho = desempenho, rep_perf = rep_perf),
        file.path(DIR, "saidas", "objetos_analise.rds"))

msg("PIPELINE CONCLUÍDA em ", round(as.numeric(difftime(Sys.time(), t_inicio, units = "mins")), 2), " min")
