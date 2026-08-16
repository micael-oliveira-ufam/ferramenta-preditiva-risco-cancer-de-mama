#!/usr/bin/env python3
"""Treino multi-algoritmo com ajuste aninhado de hiperparametros, curvas ROC
e exportacao de modelos executaveis no navegador."""
import numpy as np, pandas as pd, json, warnings, os, time
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV, cross_val_predict
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss, average_precision_score
from scipy import stats
import xgboost as xgb
warnings.filterwarnings("ignore")
SEED = 42; np.random.seed(SEED)
OUT = "saidas_multi"; os.makedirs(OUT, exist_ok=True)

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})
d = d[(d["cancer_type"] != "Breast Sarcoma") & (d["subtipo"] != "NC") &
      (d["overall_survival_months"] > 0)].reset_index(drop=True)
clin_end = list(d.columns).index("death_from_cancer")
expr = [c for c in d.columns[clin_end+1:] if not c.endswith("_mut")]
d["evento"] = 1 - d["overall_survival"]; d["tempo"] = d["overall_survival_months"]
clin_vars = ["age_at_diagnosis", "neoplasm_histologic_grade", "tumor_size",
             "lymph_nodes_examined_positive", "nottingham_prognostic_index"]
d10 = d[(d["tempo"] >= 120) | (d["evento"] == 1)].dropna(subset=clin_vars).reset_index(drop=True)
d10["obito10"] = ((d10["evento"] == 1) & (d10["tempo"] < 120)).astype(int)
y = d10["obito10"].values
classes = sorted(d10["subtipo"].unique())
Cl = pd.concat([d10[clin_vars],
                pd.get_dummies(d10["subtipo"]).reindex(columns=classes, fill_value=0)
                  .add_prefix("sub_").astype(float)], axis=1)
clin_feats = list(Cl.columns)
CONJ = {"clinico": (Cl.values.astype(float), clin_feats),
        "genes": (d10[expr].values.astype(float), expr),
        "combinado": (np.hstack([d10[expr].values, Cl.values]).astype(float), expr + clin_feats)}
print(f"n={len(y)} eventos={y.sum()} ({y.mean():.1%})", flush=True)

cv_out = StratifiedKFold(5, shuffle=True, random_state=SEED)
cv_in = StratifiedKFold(3, shuffle=True, random_state=SEED + 1)

def alg(nome, n_feat):
    """Retorna (estimador, grade). Grades pequenas: busca aninhada e cara."""
    if nome == "logistica":
        return (Pipeline([("sc", StandardScaler()),
                          ("m", LogisticRegression(max_iter=4000, random_state=SEED))]),
                {"m__C": [0.01, 0.05, 0.1, 0.5, 1.0]})
    if nome == "random_forest":
        return (RandomForestClassifier(random_state=SEED, n_jobs=-1, class_weight="balanced"),
                {"n_estimators": [150], "max_leaf_nodes": [16, 24, 32],
                 "min_samples_leaf": [3, 8, 15],
                 "max_features": ["sqrt", 0.3]})
    if nome == "extra_trees":
        return (ExtraTreesClassifier(random_state=SEED, n_jobs=-1, class_weight="balanced"),
                {"n_estimators": [150], "max_leaf_nodes": [16, 32],
                 "min_samples_leaf": [3, 8], "max_features": ["sqrt", 0.3]})
    if nome == "gradient_boosting":
        return (GradientBoostingClassifier(random_state=SEED),
                {"n_estimators": [120, 200], "learning_rate": [0.03, 0.08],
                 "max_depth": [2, 3], "subsample": [0.8, 1.0]})
    if nome == "xgboost":
        return (xgb.XGBClassifier(random_state=SEED, n_jobs=-1, eval_metric="logloss",
                                  tree_method="hist"),
                {"n_estimators": [200, 350], "learning_rate": [0.03, 0.06],
                 "max_depth": [2, 3], "subsample": [0.8], "colsample_bytree": [0.6, 0.9],
                 "reg_lambda": [1.0, 5.0], "min_child_weight": [3, 8]})
    if nome == "svm_rbf":
        return (Pipeline([("sc", StandardScaler()),
                          ("m", SVC(kernel="rbf", probability=True, random_state=SEED))]),
                {"m__C": [0.5, 1, 5], "m__gamma": ["scale", 0.001]})
    if nome == "mlp":
        return (Pipeline([("sc", StandardScaler()),
                          ("m", MLPClassifier(max_iter=900, random_state=SEED))]),
                {"m__hidden_layer_sizes": [(32,), (64,), (64, 32)],
                 "m__alpha": [0.3, 1.0, 3.0]})
    raise ValueError(nome)

ALGS = ["logistica", "random_forest", "extra_trees", "gradient_boosting", "xgboost", "svm_rbf", "mlp"]

def delong_var(y, p):
    """Variancia de DeLong para uma AUC (para o IC 95%)."""
    y = np.asarray(y); pos = p[y == 1]; neg = p[y == 0]
    m, n = len(pos), len(neg)
    def midrank(x):
        s = np.argsort(x); xs = x[s]; N = len(x); T = np.zeros(N); i = 0
        while i < N:
            j = i
            while j < N - 1 and xs[j+1] == xs[i]: j += 1
            T[i:j+1] = 0.5*(i+j) + 1; i = j+1
        o = np.empty(N); o[s] = T; return o
    tx, ty, tz = midrank(pos), midrank(neg), midrank(np.r_[pos, neg])
    auc = (tz[:m].sum() - m*(m+1)/2) / (m*n)
    v01 = (tz[:m] - tx) / n
    v10 = 1 - (tz[m:] - ty) / m
    var = np.var(v01, ddof=1)/m + np.var(v10, ddof=1)/n
    return float(auc), float(var)

resultados = []
export = {}
t0 = time.time()
for cj, (Xd, feats) in CONJ.items():
    for nome in ALGS:
        est, grade = alg(nome, Xd.shape[1])
        busca = RandomizedSearchCV(est, grade, n_iter=6, cv=cv_in, scoring="roc_auc",
                                   random_state=SEED, n_jobs=-1, refit=True)
        # AUC honesto: busca de hiperparametros dentro de cada dobra externa
        oof = cross_val_predict(clone(busca), Xd, y, cv=cv_out, method="predict_proba", n_jobs=1)[:, 1]
        auc, var = delong_var(y, oof)
        lo, hi = auc - 1.96*np.sqrt(var), auc + 1.96*np.sqrt(var)
        fpr, tpr, _ = roc_curve(y, oof)
        idx = np.unique(np.linspace(0, len(fpr)-1, 120).astype(int))
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.01, y_max=0.99).fit(oof, y)
        oof_cal = iso.predict(oof)
        g = pd.qcut(pd.Series(oof_cal).rank(method="first"), 10)
        t = pd.DataFrame({"p": oof_cal, "y": y}).groupby(g, observed=True).agg(
            n=("y", "size"), obs=("y", "sum"), esp=("p", "sum"))
        hlst = float((((t.obs - t.esp)**2) / (t.esp*(1 - t.esp/t.n))).sum())
        hlp = float(stats.chi2.sf(hlst, len(t) - 2))
        resultados.append(dict(conjunto=cj, algoritmo=nome, auc=auc, ic_inf=lo, ic_sup=hi,
                               ap=float(average_precision_score(y, oof)),
                               brier=float(brier_score_loss(y, oof)),
                               brier_cal=float(brier_score_loss(y, oof_cal)),
                               hosmer_p=hlp))
        export[f"{cj}|{nome}"] = dict(
            roc=dict(fpr=[round(float(v), 4) for v in fpr[idx]],
                     tpr=[round(float(v), 4) for v in tpr[idx]]),
            oof=oof)
        print(f"{cj:10s} {nome:18s} AUC={auc:.4f} [{lo:.3f}-{hi:.3f}] HL p={hlp:.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)

R = pd.DataFrame(resultados).sort_values("auc", ascending=False)
R.to_csv(f"{OUT}/X01_desempenho_algoritmos.csv", index=False)
print("\n=== RANKING ===")
print(R.round(4).to_string(index=False), flush=True)
np.save(f"{OUT}/oof_todos.npy", np.array([export[k]["oof"] for k in export]))
with open(f"{OUT}/chaves.json", "w") as f:
    json.dump(list(export), f)
with open(f"{OUT}/roc.json", "w") as f:
    json.dump({k: v["roc"] for k, v in export.items()}, f, separators=(",", ":"))
print("\nOK-BENCH")
