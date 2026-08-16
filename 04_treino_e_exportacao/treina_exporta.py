#!/usr/bin/env python3
"""Passagem unica: busca de hiperparametros, AUC honesto, curva ROC e exportacao
para execucao no navegador. Salva incrementalmente para sobreviver a interrupcoes."""
import numpy as np, pandas as pd, json, warnings, os, time, sys
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, GradientBoostingClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import StratifiedKFold, RandomizedSearchCV, cross_val_predict
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, roc_curve, brier_score_loss
from scipy import stats
import xgboost as xgb
warnings.filterwarnings("ignore")
SEED = 42; np.random.seed(SEED)
OUT = "saidas_multi"; os.makedirs(OUT, exist_ok=True)
PARCIAL = f"{OUT}/parcial"; os.makedirs(PARCIAL, exist_ok=True)

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
cv_out = StratifiedKFold(5, shuffle=True, random_state=SEED)
cv_in = StratifiedKFold(3, shuffle=True, random_state=SEED + 1)

def grade(nome):
    if nome == "logistica":
        return (Pipeline([("sc", StandardScaler()),
                          ("m", LogisticRegression(max_iter=4000, random_state=SEED))]),
                {"m__C": [0.01, 0.05, 0.1, 0.5, 1.0]})
    if nome == "random_forest":
        return (RandomForestClassifier(random_state=SEED, n_jobs=1, class_weight="balanced"),
                {"n_estimators": [150], "max_leaf_nodes": [16, 24, 32],
                 "min_samples_leaf": [3, 8, 15], "max_features": ["sqrt", 0.3]})
    if nome == "extra_trees":
        return (ExtraTreesClassifier(random_state=SEED, n_jobs=1, class_weight="balanced"),
                {"n_estimators": [150], "max_leaf_nodes": [16, 32],
                 "min_samples_leaf": [3, 8], "max_features": ["sqrt", 0.3]})
    if nome == "gradient_boosting":
        return (GradientBoostingClassifier(random_state=SEED),
                {"n_estimators": [120, 200], "learning_rate": [0.03, 0.08],
                 "max_depth": [2, 3], "subsample": [0.8, 1.0]})
    if nome == "xgboost":
        return (xgb.XGBClassifier(random_state=SEED, n_jobs=1, eval_metric="logloss",
                                  tree_method="hist"),
                {"n_estimators": [200, 350], "learning_rate": [0.03, 0.06],
                 "max_depth": [2, 3], "subsample": [0.8], "colsample_bytree": [0.6, 0.9],
                 "reg_lambda": [1.0, 5.0], "min_child_weight": [3, 8]})
    if nome == "mlp":
        return (Pipeline([("sc", StandardScaler()),
                          ("m", MLPClassifier(max_iter=900, random_state=SEED))]),
                {"m__hidden_layer_sizes": [(32,), (64,)], "m__alpha": [0.3, 1.0, 3.0]})
    raise ValueError(nome)

# gradient_boosting classico e proibitivo em 489 colunas; xgboost cobre a familia de boosting
TAREFAS = []
for cj in ["clinico", "genes", "combinado"]:
    for a in ["logistica", "random_forest", "extra_trees", "gradient_boosting", "xgboost", "mlp"]:
        if a == "gradient_boosting" and cj != "clinico":
            continue
        TAREFAS.append((cj, a))

def arvore_sk(tree, prob=True):
    t = tree.tree_
    p = {"f": [(-1 if v < 0 else int(v)) for v in t.feature],
         "t": [round(float(v), 6) for v in t.threshold],
         "l": t.children_left.tolist(), "r": t.children_right.tolist()}
    if prob:
        v = t.value[:, 0, :]
        p["v"] = [round(float(r[1]/r.sum()), 6) if r.sum() > 0 else 0.0 for r in v]
    else:
        p["v"] = [round(float(t.value[i][0][0]), 6) for i in range(t.node_count)]
    return p

def arvores_xgb(modelo):
    saida = []
    for dj in modelo.get_booster().get_dump(dump_format="json"):
        raiz = json.loads(dj)
        F, T, L, R, V = [], [], [], [], []
        def visita(n):
            i = len(F)
            F.append(-1); T.append(0.0); L.append(-1); R.append(-1); V.append(0.0)
            if "leaf" in n:
                V[i] = float(np.float32(n["leaf"])); return i
            s = str(n["split"])
            F[i] = int(s[1:]) if s.startswith("f") else int(s)
            T[i] = float(np.float32(n["split_condition"]))
            filhos = {c["nodeid"]: c for c in n["children"]}
            L[i] = visita(filhos[n["yes"]])
            R[i] = visita(filhos[n["no"]])
            return i
        visita(raiz)
        saida.append({"f": F, "t": T, "l": L, "r": R, "v": V})
    return saida

def delong(y, p):
    pos = p[y == 1]; neg = p[y == 0]; m, n = len(pos), len(neg)
    def midrank(x):
        s = np.argsort(x); xs = x[s]; N = len(x); T = np.zeros(N); i = 0
        while i < N:
            j = i
            while j < N-1 and xs[j+1] == xs[i]: j += 1
            T[i:j+1] = 0.5*(i+j)+1; i = j+1
        o = np.empty(N); o[s] = T; return o
    tx, ty, tz = midrank(pos), midrank(neg), midrank(np.r_[pos, neg])
    auc = (tz[:m].sum() - m*(m+1)/2)/(m*n)
    v01 = (tz[:m]-tx)/n; v10 = 1-(tz[m:]-ty)/m
    return float(auc), float(np.var(v01, ddof=1)/m + np.var(v10, ddof=1)/n)

def confere(p, X, ref):
    def desce(a, x):
        i = 0
        while a["l"][i] != -1:
            i = a["l"][i] if x[a["f"][i]] < a["t"][i] else a["r"][i]
        return a["v"][i]
    out = []
    for x in X:
        if p["tipo"] == "linear":
            z = p["intercepto"] + sum(((x[i]-p["media"][i])/p["escala"][i])*p["coef"][i]
                                      for i in range(len(x)))
            out.append(1/(1+np.exp(-z)))
        elif p["tipo"] == "floresta":
            out.append(float(np.mean([desce(a, x) for a in p["arvores"]])))
        elif p["tipo"] == "boosting_sk":
            z = p["base"] + p["lr"]*sum(desce(a, x) for a in p["arvores"])
            out.append(1/(1+np.exp(-z)))
        elif p["tipo"] == "boosting_xgb":
            z = np.log(p["base"]/(1-p["base"])) + sum(desce(a, x) for a in p["arvores"])
            out.append(1/(1+np.exp(-z)))
        elif p["tipo"] == "mlp":
            h = (np.array(x)-np.array(p["media"]))/np.array(p["escala"])
            for W, b in zip(p["W"], p["b"]):
                h = np.maximum(np.array(W).T @ h + np.array(b), 0)
            z = float(np.array(p["Wout"]) @ h + p["bout"])
            out.append(1/(1+np.exp(-z)))
    return float(np.max(np.abs(np.array(out)-ref)))

t0 = time.time()
for cj, nome in TAREFAS:
    chave = f"{cj}|{nome}"
    alvo = f"{PARCIAL}/{cj}__{nome}.json"
    if os.path.exists(alvo):
        print(f"[pular] {chave}", flush=True); continue
    est, g = grade(nome)
    busca = RandomizedSearchCV(est, g, n_iter=4, cv=cv_in, scoring="roc_auc",
                               random_state=SEED, n_jobs=1, refit=True)
    busca.fit(CONJ[cj][0], y)
    melhor = busca.best_estimator_
    Xd, feats = CONJ[cj]
    oof = cross_val_predict(melhor, Xd, y, cv=cv_out, method="predict_proba")[:, 1]
    auc, var = delong(y, oof)
    lo, hi = auc-1.96*np.sqrt(var), auc+1.96*np.sqrt(var)
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.01, y_max=0.99).fit(oof, y)
    cal = iso.predict(oof)
    q = pd.qcut(pd.Series(cal).rank(method="first"), 10)
    tt = pd.DataFrame({"p": cal, "y": y}).groupby(q, observed=True).agg(
        n=("y", "size"), obs=("y", "sum"), esp=("p", "sum"))
    hlst = float((((tt.obs-tt.esp)**2)/(tt.esp*(1-tt.esp/tt.n))).sum())
    hlp = float(stats.chi2.sf(hlst, len(tt)-2))
    fpr, tpr, _ = roc_curve(y, oof)
    idx = np.unique(np.linspace(0, len(fpr)-1, 110).astype(int))
    xs = np.linspace(0.001, 0.999, 160)
    p = {"conjunto": cj, "algoritmo": nome, "variaveis": feats,
         "hiper": {k: str(v) for k, v in busca.best_params_.items()},
         "auc": round(auc, 4), "ic_inf": round(lo, 4), "ic_sup": round(hi, 4),
         "brier": round(float(brier_score_loss(y, oof)), 4),
         "brier_calibrado": round(float(brier_score_loss(y, cal)), 4),
         "hosmer_p": round(hlp, 4),
         "roc": {"fpr": [round(float(v), 4) for v in fpr[idx]],
                 "tpr": [round(float(v), 4) for v in tpr[idx]]},
         "calibracao": {"x": [round(float(v), 5) for v in xs],
                        "y": [round(float(v), 5) for v in iso.predict(xs)]}}
    if nome == "logistica":
        sc = melhor.named_steps["sc"]; lr = melhor.named_steps["m"]
        p.update(tipo="linear", media=[round(float(v), 6) for v in sc.mean_],
                 escala=[round(float(v), 6) for v in sc.scale_],
                 coef=[round(float(v), 6) for v in lr.coef_[0]],
                 intercepto=round(float(lr.intercept_[0]), 6))
    elif nome in ("random_forest", "extra_trees"):
        p.update(tipo="floresta", arvores=[arvore_sk(t) for t in melhor.estimators_])
    elif nome == "gradient_boosting":
        pri = float(melhor.init_.class_prior_[1])
        p.update(tipo="boosting_sk", base=round(float(np.log(pri/(1-pri))), 6),
                 lr=float(melhor.learning_rate),
                 arvores=[arvore_sk(e[0], prob=False) for e in melhor.estimators_])
    elif nome == "xgboost":
        cfg = json.loads(melhor.get_booster().save_config())
        bs = float(str(cfg["learner"]["learner_model_param"]["base_score"]).strip("[]"))
        p.update(tipo="boosting_xgb", base=round(bs, 8),
                 arvores=arvores_xgb(melhor))
    elif nome == "mlp":
        sc = melhor.named_steps["sc"]; m = melhor.named_steps["m"]
        p.update(tipo="mlp", media=[round(float(v), 6) for v in sc.mean_],
                 escala=[round(float(v), 6) for v in sc.scale_],
                 W=[[[round(float(v), 6) for v in linha] for linha in W]
                    for W in m.coefs_[:-1]],
                 b=[[round(float(v), 6) for v in b] for b in m.intercepts_[:-1]],
                 Wout=[round(float(v), 6) for v in m.coefs_[-1][:, 0]],
                 bout=round(float(m.intercepts_[-1][0]), 6))
    Xv = Xd[:12].astype(np.float32).astype(np.float64)
    erro = confere(p, Xv, melhor.predict_proba(Xd[:12])[:, 1])
    p["erro_verificacao"] = float(erro)
    with open(alvo, "w") as f:
        json.dump(p, f, separators=(",", ":"))
    print(f"{chave:28s} AUC={auc:.4f} [{lo:.3f}-{hi:.3f}] HLp={hlp:.3f} "
          f"| {os.path.getsize(alvo)/1024:6.0f} KB | erro JS={erro:.1e} "
          f"| {time.time()-t0:.0f}s", flush=True)
print("\nOK-TODOS")
