#!/usr/bin/env python3
"""Retomada do trecho final da Tarefa 2: tipologia de erro, calibracao, log-rank."""
import numpy as np, pandas as pd, warnings
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED = 42
OUT = "saidas_ml"

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})
d = d[(d["cancer_type"] != "Breast Sarcoma") & (d["subtipo"] != "NC") &
      (d["overall_survival_months"] > 0)].reset_index(drop=True)
d["evento"] = 1 - d["overall_survival"]; d["tempo"] = d["overall_survival_months"]
clin_vars = ["age_at_diagnosis", "neoplasm_histologic_grade", "tumor_size",
             "lymph_nodes_examined_positive", "nottingham_prognostic_index"]
d10 = d[(d["tempo"] >= 120) | (d["evento"] == 1)].dropna(subset=clin_vars).reset_index(drop=True)
d10["obito10"] = ((d10["evento"] == 1) & (d10["tempo"] < 120)).astype(int)

oof = pd.read_csv(f"{OUT}/M14_predicoes_oof_sobrevida.csv")
assert len(oof) == len(d10)
d10["risco"] = oof["combinado | Gradient Boosting"].values
d10["risco_genes"] = oof["genes (489) | Random Forest"].values
d10["risco_clin"] = oof["clinico (6) | Random Forest"].values
d10["pred"] = (d10["risco"] > d10["risco"].median()).astype(int)
d10["acertou"] = (d10["pred"] == d10["obito10"]).astype(int)
d10["estagio"] = d10["tumor_stage"].map({0:"I",1:"I",2:"II",3:"III-IV",4:"III-IV"})
print("n =", len(d10), "| obitos em 10a =", int(d10.obito10.sum()),
      "| acuracia (corte mediana) =", round(d10.acertou.mean(), 3))

# ---- AUC por estrato com IC bootstrap, para os tres modelos ----
def auc_ic(yy, pp, B=2000, seed=SEED):
    rng = np.random.default_rng(seed); idx = np.arange(len(yy)); vals = []
    for _ in range(B):
        s = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(yy[s])) < 2: continue
        vals.append(roc_auc_score(yy[s], pp[s]))
    return roc_auc_score(yy, pp), np.percentile(vals, 2.5), np.percentile(vals, 97.5)

linhas = []
grupos = [("coorte inteira", d10)]
for v in ["LumA","LumB","Her2","Basal","claudin-low","Normal"]:
    grupos.append((f"subtipo {v}", d10[d10.subtipo == v]))
for v in ["I","II","III-IV"]:
    grupos.append((f"estagio {v}", d10[d10.estagio == v]))
for nome, s in grupos:
    if len(s) < 50 or s["obito10"].nunique() < 2:
        linhas.append(dict(estrato=nome, n=len(s), eventos=int(s.obito10.sum()))); continue
    a, lo, hi = auc_ic(s.obito10.values, s.risco.values, B=1500)
    ag = roc_auc_score(s.obito10, s.risco_genes)
    ac = roc_auc_score(s.obito10, s.risco_clin)
    linhas.append(dict(estrato=nome, n=len(s), eventos=int(s.obito10.sum()),
                       AUC_combinado=a, IC_inf=lo, IC_sup=hi,
                       AUC_genes=ag, AUC_clinico=ac, acuracia=s.acertou.mean(),
                       significativo=("sim" if lo > 0.5 else "NAO")))
E = pd.DataFrame(linhas)
E.to_csv(f"{OUT}/M17_auc_por_estrato.csv", index=False)
print("\n=== Onde acerta e onde falha (AUC com IC 95% bootstrap) ===")
print(E.round(3).to_string(index=False))

for col in ["subtipo", "estagio"]:
    s = d10.dropna(subset=[col])
    c, pv, _, _ = stats.chi2_contingency(pd.crosstab(s[col], s["acertou"]))
    print(f"\nA taxa de acerto depende de {col}? chi2={c:.2f} p={pv:.3g}")
    print(s.groupby(col)["acertou"].agg(["mean","size"]).round(3).to_string())

# ---- tipologia de erro ----
d10["tipo"] = np.select(
    [(d10.pred==1)&(d10.obito10==1), (d10.pred==0)&(d10.obito10==0),
     (d10.pred==1)&(d10.obito10==0), (d10.pred==0)&(d10.obito10==1)],
    ["acerto: alto risco / obito", "acerto: baixo risco / viva",
     "FALSO POSITIVO (alarme sem obito)", "FALSO NEGATIVO (obito nao previsto)"],
    default="")
print("\n=== Tipos de resultado ===")
print(d10["tipo"].value_counts().to_string())
erros = d10[d10.acertou == 0]
print(f"\nEntre os {len(erros)} erros: {(erros.tipo.str.startswith('FALSO POSITIVO')).sum()} falsos positivos, "
      f"{(erros.tipo.str.startswith('FALSO NEGATIVO')).sum()} falsos negativos")
fp = int((erros.tipo.str.startswith('FALSO POSITIVO')).sum()); fn = len(erros) - fp
b = stats.binomtest(fp, len(erros), 0.5)
print(f"Teste binomial (o modelo erra mais para um lado?): p = {b.pvalue:.3g}")

print("\nComposicao dos erros por subtipo (proporcao dentro do subtipo):")
print(pd.crosstab(d10.subtipo, d10.tipo, normalize="index").round(3).to_string())
d10.groupby("tipo")[["age_at_diagnosis","tumor_size","nottingham_prognostic_index",
                     "lymph_nodes_examined_positive","tempo","risco"]].median().round(2) \
   .to_csv(f"{OUT}/M19_tipos_de_erro.csv")
print("\nMedianas por tipo de resultado:")
print(d10.groupby("tipo")[["age_at_diagnosis","tumor_size","nottingham_prognostic_index",
                           "lymph_nodes_examined_positive","tempo","risco"]].median().round(2).to_string())

# falso negativo vs acerto de obito: o que diferencia?
fn_df = d10[d10.tipo.str.startswith("FALSO NEGATIVO")]
tp_df = d10[d10.tipo.str.startswith("acerto: alto risco")]
print("\n=== O que distingue um obito NAO previsto de um obito previsto? ===")
tt = []
for v in ["age_at_diagnosis","tumor_size","nottingham_prognostic_index",
          "lymph_nodes_examined_positive","tempo"]:
    u, pv = stats.mannwhitneyu(fn_df[v].dropna(), tp_df[v].dropna())
    tt.append(dict(variavel=v, mediana_nao_previsto=fn_df[v].median(),
                   mediana_previsto=tp_df[v].median(), p=pv))
from statsmodels.stats.multitest import multipletests
TT = pd.DataFrame(tt); TT["fdr"] = multipletests(TT.p, method="fdr_bh")[1]
print(TT.round(4).to_string(index=False))
TT.to_csv(f"{OUT}/M21_falsos_negativos_vs_acertos.csv", index=False)
c, pv, _, _ = stats.chi2_contingency(pd.crosstab(d10.subtipo, d10.tipo))
print(f"\nTipo de erro depende do subtipo? chi2={c:.2f} p={pv:.3g}")

# ---- calibracao ----
faixas = pd.cut(d10["risco"], [0,.2,.3,.4,.5,.6,.7,1.0])
cal = d10.groupby(faixas, observed=True).agg(n=("obito10","size"), risco_previsto=("risco","mean"),
                                             obito_observado=("obito10","mean"))
cal["diferenca"] = cal["risco_previsto"] - cal["obito_observado"]
cal.to_csv(f"{OUT}/M20_calibracao_sobrevida.csv")
print("\n=== Calibracao: o risco previsto corresponde ao observado? ===")
print(cal.round(3).to_string())
# Hosmer-Lemeshow
g = pd.qcut(d10["risco"], 10, duplicates="drop")
hl_tab = d10.groupby(g, observed=True).agg(n=("obito10","size"), obs=("obito10","sum"), esp=("risco","sum"))
hl = (((hl_tab.obs - hl_tab.esp)**2) / (hl_tab.esp * (1 - hl_tab.esp/hl_tab.n))).sum()
print(f"Hosmer-Lemeshow: chi2 = {hl:.2f}, gl = {len(hl_tab)-2}, p = {stats.chi2.sf(hl, len(hl_tab)-2):.3g}")

# ---- tercis de risco e sobrevida real ----
d10["tercil"] = pd.qcut(d10["risco"], 3, labels=["baixo","medio","alto"])
lr = multivariate_logrank_test(d10["tempo"], d10["tercil"], d10["evento"])
print(f"\n=== Tercis de risco previsto vs sobrevida real — log-rank p = {lr.p_value:.3g} ===")
for t in ["baixo","medio","alto"]:
    s = d10[d10.tercil == t]; k = KaplanMeierFitter().fit(s["tempo"], s["evento"])
    print(f"{t}: n={len(s)} obitos={int(s.evento.sum())} mediana OS={k.median_survival_time_:.1f} meses")
d10[["patient_id","subtipo","estagio","obito10","risco","pred","acertou","tipo","tercil"]] \
  .to_csv(f"{OUT}/M16_predicoes_paciente_sobrevida.csv", index=False)
print("\nOK-FINAL")
