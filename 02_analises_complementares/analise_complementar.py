#!/usr/bin/env python3
"""Analises complementares a pipeline METABRIC em R:
   (A) Cox por estagio tumoral (489 genes)
   (B) Interacao gene x tratamento (leitura farmacogenomica)
   (C) Random Forest: subtipo (replicacao) e mortalidade em 10 anos
"""
import numpy as np, pandas as pd, warnings
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)
OUT = "saidas_complementares"
import os; os.makedirs(OUT, exist_ok=True)

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})

# ---- coorte analitica (mesmos criterios do pipeline em R) ----
d = d[d["cancer_type"] != "Breast Sarcoma"]
d = d[d["subtipo"] != "NC"]
d = d[d["overall_survival_months"] > 0]
print("coorte:", d.shape)

# blocos de colunas
clin_end = list(d.columns).index("death_from_cancer")
expr_cols = [c for c in d.columns[clin_end+1:] if not c.endswith("_mut")]
print("genes expressao:", len(expr_cols))

d["evento_os"] = 1 - d["overall_survival"]
d["tempo"] = d["overall_survival_months"]
d["evento_dss"] = np.where(d["death_from_cancer"] == "Died of Disease", 1, 0)

# ---- grupos de estagio ----
def grupo_estagio(s):
    if pd.isna(s): return np.nan
    if s in (0, 1): return "I"
    if s == 2: return "II"
    if s in (3, 4): return "III-IV"
    return np.nan
d["estagio_grp"] = d["tumor_stage"].apply(grupo_estagio)
print(d["estagio_grp"].value_counts(dropna=False))

def cox_uni(df, gene, extra=None):
    cols = ["tempo", "evento_os", gene] + (extra or [])
    sub = df[cols].dropna()
    if sub[gene].std() == 0 or sub["evento_os"].sum() < 10: return None
    cph = CoxPHFitter()
    try:
        cph.fit(sub, "tempo", "evento_os")
    except Exception:
        return None
    s = cph.summary
    return dict(HR=float(np.exp(s.loc[gene, "coef"])),
                ic_inf=float(np.exp(s.loc[gene, "coef lower 95%"])),
                ic_sup=float(np.exp(s.loc[gene, "coef upper 95%"])),
                p=float(s.loc[gene, "p"]), n=len(sub), eventos=int(sub["evento_os"].sum()))

# ================= (A) Cox por estagio =================
linhas = []
for est in ["I", "II", "III-IV"]:
    sub = d[d["estagio_grp"] == est]
    for g in expr_cols:
        r = cox_uni(sub, g)
        if r: linhas.append(dict(estagio=est, gene=g, **r))
A = pd.DataFrame(linhas)
A["fdr"] = np.nan
for est in A["estagio"].unique():
    m = A["estagio"] == est
    A.loc[m, "fdr"] = multipletests(A.loc[m, "p"], method="fdr_bh")[1]
A = A.sort_values(["estagio", "fdr"])
A.to_csv(f"{OUT}/C01_cox_por_estagio_todos_genes.csv", index=False)
print("\n=== A: significativos (FDR<0.05) por estagio ===")
print(A[A.fdr < 0.05].groupby("estagio").size())
for est in ["I", "II", "III-IV"]:
    s = A[(A.estagio == est) & (A.fdr < 0.05)].head(12)
    print(f"\n--- {est} (n={s.n.max() if len(s) else 'NA'}) top ---")
    print(s[["gene", "HR", "ic_inf", "ic_sup", "fdr", "eventos"]].to_string(index=False))

# ================= (B) Interacao gene x tratamento =================
inter = []
for trat in ["hormone_therapy", "chemotherapy", "radio_therapy"]:
    for g in expr_cols:
        sub = d[["tempo", "evento_os", g, trat, "age_at_diagnosis",
                 "neoplasm_histologic_grade", "tumor_size"]].dropna()
        sub = sub.rename(columns={g: "gene", trat: "trat"})
        sub["gene_x_trat"] = sub["gene"] * sub["trat"]
        cph = CoxPHFitter()
        try:
            cph.fit(sub, "tempo", "evento_os")
            s = cph.summary
            inter.append(dict(tratamento=trat, gene=g,
                              HR_interacao=float(np.exp(s.loc["gene_x_trat", "coef"])),
                              p_interacao=float(s.loc["gene_x_trat", "p"]),
                              HR_gene_sem_trat=float(np.exp(s.loc["gene", "coef"])),
                              n=len(sub)))
        except Exception:
            continue
B = pd.DataFrame(inter)
B["fdr"] = np.nan
for t in B["tratamento"].unique():
    m = B["tratamento"] == t
    B.loc[m, "fdr"] = multipletests(B.loc[m, "p_interacao"], method="fdr_bh")[1]
B = B.sort_values(["tratamento", "fdr"])
B.to_csv(f"{OUT}/C02_interacao_gene_tratamento.csv", index=False)
print("\n=== B: interacoes gene x tratamento (FDR<0.05) ===")
print(B[B.fdr < 0.05].groupby("tratamento").size())
for t in B["tratamento"].unique():
    print(f"\n--- {t} top 10 ---")
    print(B[B.tratamento == t].head(10)[["gene", "HR_interacao", "p_interacao", "fdr"]].to_string(index=False))

# HR do gene dentro de cada estrato de tratamento, para os genes-nucleo
nucleo = ["gsk3b", "stat5a", "aurka", "vegfa", "bcl2", "spry2", "igf1", "abcb1",
          "flt3", "cdkn2c", "pdgfra", "erbb2", "esr1", "mapt", "ar"]
nucleo = [g for g in nucleo if g in expr_cols]
estr = []
for trat in ["hormone_therapy", "chemotherapy"]:
    for val, rot in [(1, "tratado"), (0, "nao tratado")]:
        sub = d[d[trat] == val]
        for g in nucleo:
            r = cox_uni(sub, g)
            if r: estr.append(dict(tratamento=trat, estrato=rot, gene=g, **r))
E = pd.DataFrame(estr)
E["fdr"] = multipletests(E["p"], method="fdr_bh")[1]
E.to_csv(f"{OUT}/C03_cox_por_estrato_tratamento.csv", index=False)
print("\n=== B2: HR por estrato de tratamento (genes-nucleo) ===")
print(E.pivot_table(index="gene", columns=["tratamento", "estrato"], values="HR").round(3).to_string())

# ================= (C) Random Forest =================
X = d[expr_cols].values
y = d["subtipo"].values
rf = RandomForestClassifier(n_estimators=500, oob_score=True, random_state=SEED, n_jobs=-1)
rf.fit(X, y)
print("\n=== C: RF subtipo — acuracia OOB:", round(rf.oob_score_, 4), " erro OOB:", round(1-rf.oob_score_, 4))
imp = pd.DataFrame({"gene": expr_cols, "importancia_gini": rf.feature_importances_}) \
        .sort_values("importancia_gini", ascending=False).reset_index(drop=True)
imp["posicao"] = imp.index + 1
imp.to_csv(f"{OUT}/C04_rf_subtipo_importancia.csv", index=False)
print(imp.head(20).to_string(index=False))

# RF de mortalidade em 10 anos (120 meses)
d10 = d[(d["tempo"] >= 120) | (d["evento_os"] == 1)].copy()
d10["obito_10a"] = ((d10["evento_os"] == 1) & (d10["tempo"] < 120)).astype(int)
Xg = d10[expr_cols].values
yg = d10["obito_10a"].values
print("\nRF mortalidade 10a: n =", len(yg), " eventos =", yg.sum())
cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
rf2 = RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1,
                             class_weight="balanced", min_samples_leaf=5)
p_gen = cross_val_predict(rf2, Xg, yg, cv=cv, method="predict_proba")[:, 1]
auc_gen = roc_auc_score(yg, p_gen)
clin_vars = ["age_at_diagnosis", "neoplasm_histologic_grade", "tumor_size",
             "lymph_nodes_examined_positive", "nottingham_prognostic_index"]
dc = d10.dropna(subset=clin_vars)
Xc = pd.get_dummies(dc[clin_vars + ["subtipo"]], columns=["subtipo"]).values.astype(float)
yc = dc["obito_10a"].values
p_clin = cross_val_predict(rf2, Xc, yc, cv=cv, method="predict_proba")[:, 1]
auc_clin = roc_auc_score(yc, p_clin)
Xcomb = np.hstack([dc[expr_cols].values, Xc])
p_comb = cross_val_predict(rf2, Xcomb, yc, cv=cv, method="predict_proba")[:, 1]
auc_comb = roc_auc_score(yc, p_comb)
print(f"AUC (5-fold CV) — genes: {auc_gen:.3f} | clinico: {auc_clin:.3f} | combinado: {auc_comb:.3f}")
rf2.fit(Xg, yg)
imp2 = pd.DataFrame({"gene": expr_cols, "importancia_gini": rf2.feature_importances_}) \
         .sort_values("importancia_gini", ascending=False).reset_index(drop=True)
imp2["posicao"] = imp2.index + 1
imp2.to_csv(f"{OUT}/C05_rf_mortalidade10a_importancia.csv", index=False)
print(imp2.head(20).to_string(index=False))
pd.DataFrame([{"modelo": "genes (489)", "AUC": auc_gen, "n": len(yg)},
              {"modelo": "clinico", "AUC": auc_clin, "n": len(yc)},
              {"modelo": "combinado", "AUC": auc_comb, "n": len(yc)}]) \
  .to_csv(f"{OUT}/C06_rf_auc_mortalidade10a.csv", index=False)

# RF por subtipo (luminais vs nao luminais) para mortalidade
for grp, nome in [(["LumA", "LumB"], "luminais"), (["Basal", "Her2"], "basal_her2")]:
    sg = d10[d10["subtipo"].isin(grp)]
    if sg["obito_10a"].sum() < 30: continue
    p = cross_val_predict(rf2, sg[expr_cols].values, sg["obito_10a"].values,
                          cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                          method="predict_proba")[:, 1]
    print(f"AUC RF (genes) em {nome}: {roc_auc_score(sg['obito_10a'], p):.3f} "
          f"(n={len(sg)}, eventos={int(sg['obito_10a'].sum())})")
print("\nOK")
