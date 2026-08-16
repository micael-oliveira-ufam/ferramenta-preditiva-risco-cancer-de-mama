#!/usr/bin/env python3
"""(C) Random Forest: subtipo (replicacao) e mortalidade em 10 anos."""
import numpy as np, pandas as pd, warnings, os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
SEED = 42; np.random.seed(SEED)
OUT = "saidas_complementares"; os.makedirs(OUT, exist_ok=True)

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})
d = d[(d["cancer_type"] != "Breast Sarcoma") & (d["subtipo"] != "NC") &
      (d["overall_survival_months"] > 0)]
clin_end = list(d.columns).index("death_from_cancer")
expr_cols = [c for c in d.columns[clin_end+1:] if not c.endswith("_mut")]
d["evento_os"] = 1 - d["overall_survival"]
d["tempo"] = d["overall_survival_months"]

X = d[expr_cols].values; y = d["subtipo"].values
rf = RandomForestClassifier(n_estimators=500, oob_score=True, random_state=SEED, n_jobs=-1)
rf.fit(X, y)
print("RF subtipo — acuracia OOB:", round(rf.oob_score_, 4), "| erro OOB:", round(1-rf.oob_score_, 4))
imp = pd.DataFrame({"gene": expr_cols, "importancia_gini": rf.feature_importances_}) \
        .sort_values("importancia_gini", ascending=False).reset_index(drop=True)
imp["posicao"] = imp.index + 1
imp.to_csv(f"{OUT}/C04_rf_subtipo_importancia.csv", index=False)
print(imp.head(25).to_string(index=False))

d10 = d[(d["tempo"] >= 120) | (d["evento_os"] == 1)].copy()
d10["obito_10a"] = ((d10["evento_os"] == 1) & (d10["tempo"] < 120)).astype(int)
Xg = d10[expr_cols].values; yg = d10["obito_10a"].values
print("\nRF mortalidade 10a: n =", len(yg), "| eventos =", int(yg.sum()))
cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
rf2 = RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1,
                             class_weight="balanced", min_samples_leaf=5)
auc_gen = roc_auc_score(yg, cross_val_predict(rf2, Xg, yg, cv=cv, method="predict_proba")[:, 1])
clin_vars = ["age_at_diagnosis", "neoplasm_histologic_grade", "tumor_size",
             "lymph_nodes_examined_positive", "nottingham_prognostic_index"]
dc = d10.dropna(subset=clin_vars)
Xc = pd.get_dummies(dc[clin_vars + ["subtipo"]], columns=["subtipo"]).values.astype(float)
yc = dc["obito_10a"].values
auc_clin = roc_auc_score(yc, cross_val_predict(rf2, Xc, yc, cv=cv, method="predict_proba")[:, 1])
Xcomb = np.hstack([dc[expr_cols].values, Xc])
auc_comb = roc_auc_score(yc, cross_val_predict(rf2, Xcomb, yc, cv=cv, method="predict_proba")[:, 1])
print(f"AUC (5-fold CV) — genes: {auc_gen:.3f} | clinico: {auc_clin:.3f} | combinado: {auc_comb:.3f} (n clin={len(yc)})")
rf2.fit(Xg, yg)
imp2 = pd.DataFrame({"gene": expr_cols, "importancia_gini": rf2.feature_importances_}) \
         .sort_values("importancia_gini", ascending=False).reset_index(drop=True)
imp2["posicao"] = imp2.index + 1
imp2.to_csv(f"{OUT}/C05_rf_mortalidade10a_importancia.csv", index=False)
print(imp2.head(25).to_string(index=False))
pd.DataFrame([{"modelo": "genes (489)", "AUC": auc_gen, "n": len(yg)},
              {"modelo": "clinico", "AUC": auc_clin, "n": len(yc)},
              {"modelo": "combinado (genes+clinico)", "AUC": auc_comb, "n": len(yc)}]) \
  .to_csv(f"{OUT}/C06_rf_auc_mortalidade10a.csv", index=False)

linhas = []
for grp, nome in [(["LumA", "LumB"], "luminais (LumA+LumB)"), (["Basal", "Her2"], "Basal+Her2"),
                  (["LumA"], "LumA"), (["LumB"], "LumB"), (["Basal"], "Basal"), (["Her2"], "Her2"),
                  (["claudin-low"], "claudin-low")]:
    sg = d10[d10["subtipo"].isin(grp)]
    if sg["obito_10a"].sum() < 25 or (len(sg) - sg["obito_10a"].sum()) < 25: 
        print("pulado (poucos eventos):", nome); continue
    p = cross_val_predict(rf2, sg[expr_cols].values, sg["obito_10a"].values,
                          cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                          method="predict_proba")[:, 1]
    a = roc_auc_score(sg["obito_10a"], p)
    print(f"AUC RF (genes) em {nome}: {a:.3f} (n={len(sg)}, eventos={int(sg['obito_10a'].sum())})")
    linhas.append(dict(estrato=nome, AUC=a, n=len(sg), eventos=int(sg["obito_10a"].sum())))

# por estagio
def gest(s):
    if pd.isna(s): return np.nan
    if s in (0, 1): return "I"
    if s == 2: return "II"
    if s in (3, 4): return "III-IV"
    return np.nan
d10["estagio_grp"] = d10["tumor_stage"].apply(gest)
for est in ["I", "II", "III-IV"]:
    sg = d10[d10["estagio_grp"] == est]
    if sg["obito_10a"].sum() < 25 or (len(sg) - sg["obito_10a"].sum()) < 25:
        print("pulado (poucos eventos): estagio", est); continue
    p = cross_val_predict(rf2, sg[expr_cols].values, sg["obito_10a"].values,
                          cv=StratifiedKFold(5, shuffle=True, random_state=SEED),
                          method="predict_proba")[:, 1]
    a = roc_auc_score(sg["obito_10a"], p)
    print(f"AUC RF (genes) no estagio {est}: {a:.3f} (n={len(sg)}, eventos={int(sg['obito_10a'].sum())})")
    linhas.append(dict(estrato=f"estagio {est}", AUC=a, n=len(sg), eventos=int(sg["obito_10a"].sum())))
pd.DataFrame(linhas).to_csv(f"{OUT}/C07_rf_auc_por_estrato.csv", index=False)
print("OK-C")
