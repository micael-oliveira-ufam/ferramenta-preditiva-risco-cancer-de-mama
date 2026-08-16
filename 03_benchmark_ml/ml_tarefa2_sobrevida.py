#!/usr/bin/env python3
"""Benchmark multi-modelo — TAREFA 2: predicao de obito em ate 10 anos.
Compara algoritmos e conjuntos de variaveis, com testes estatisticos e analise de erro."""
import numpy as np, pandas as pd, warnings, os, time
from sklearn.base import clone
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss, average_precision_score
from scipy import stats
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings("ignore")
SEED = 42; np.random.seed(SEED)
OUT = "saidas_ml"; os.makedirs(OUT, exist_ok=True)

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})
d = d[(d["cancer_type"] != "Breast Sarcoma") & (d["subtipo"] != "NC") &
      (d["overall_survival_months"] > 0)].reset_index(drop=True)
clin_end = list(d.columns).index("death_from_cancer")
expr = [c for c in d.columns[clin_end+1:] if not c.endswith("_mut")]
d["evento"] = 1 - d["overall_survival"]; d["tempo"] = d["overall_survival_months"]

# coorte com desfecho de 10 anos determinado (evita censura informativa)
clin_vars = ["age_at_diagnosis", "neoplasm_histologic_grade", "tumor_size",
             "lymph_nodes_examined_positive", "nottingham_prognostic_index"]
d10 = d[(d["tempo"] >= 120) | (d["evento"] == 1)].dropna(subset=clin_vars).reset_index(drop=True)
d10["obito10"] = ((d10["evento"] == 1) & (d10["tempo"] < 120)).astype(int)
y = d10["obito10"].values
Xg = d10[expr].values
Xc = pd.get_dummies(d10[clin_vars + ["subtipo"]], columns=["subtipo"]).astype(float).values
Xk = np.hstack([Xg, Xc])
print(f"coorte: n={len(y)} obitos10a={y.sum()} ({y.mean():.1%})", flush=True)

MODELOS = {
    "Baseline (prevalencia)": DummyClassifier(strategy="prior"),
    "Regressao logistica (L2)": make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, C=0.05, random_state=SEED)),
    "Regressao logistica (LASSO)": make_pipeline(StandardScaler(), LogisticRegression(max_iter=5000, C=0.05, penalty="l1", solver="liblinear", random_state=SEED)),
    "SVM linear": make_pipeline(StandardScaler(), SVC(kernel="linear", C=0.01, probability=True, random_state=SEED)),
    "SVM RBF": make_pipeline(StandardScaler(), SVC(kernel="rbf", C=1, probability=True, random_state=SEED)),
    "Random Forest": RandomForestClassifier(n_estimators=500, min_samples_leaf=5, class_weight="balanced", random_state=SEED, n_jobs=-1),
    "Extra Trees": ExtraTreesClassifier(n_estimators=500, min_samples_leaf=5, class_weight="balanced", random_state=SEED, n_jobs=-1),
    "Gradient Boosting": HistGradientBoostingClassifier(max_iter=250, learning_rate=0.05, random_state=SEED),
    "k-NN (k=25)": make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=25)),
    "LDA": make_pipeline(StandardScaler(), LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")),
    "Naive Bayes": make_pipeline(StandardScaler(), GaussianNB()),
    "Rede neural (MLP)": make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(64,), alpha=1.0, max_iter=800, random_state=SEED)),
}
CONJUNTOS = {"genes (489)": Xg, "clinico (6)": Xc, "combinado": Xk}

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=SEED)
folds = list(rskf.split(Xg, y))
linhas = []; t0 = time.time()
for cj, Xd in CONJUNTOS.items():
    for nome, mod in MODELOS.items():
        for i, (tr, te) in enumerate(folds):
            m = clone(mod); m.fit(Xd[tr], y[tr])
            p = m.predict_proba(Xd[te])[:, 1]
            linhas.append(dict(conjunto=cj, modelo=nome, fold=i,
                               auc=roc_auc_score(y[te], p),
                               ap=average_precision_score(y[te], p),
                               brier=brier_score_loss(y[te], p)))
        print(f"{cj} | {nome}: ok ({time.time()-t0:.0f}s)", flush=True)
F = pd.DataFrame(linhas)
F.to_csv(f"{OUT}/M10_desempenho_por_fold_sobrevida.csv", index=False)
resumo = F.groupby(["conjunto", "modelo"]).agg(
    AUC_media=("auc", "mean"), AUC_dp=("auc", "std"),
    AP_media=("ap", "mean"), Brier_medio=("brier", "mean")).sort_values("AUC_media", ascending=False)
resumo.to_csv(f"{OUT}/M11_resumo_modelos_sobrevida.csv")
print("\n=== RESUMO (15 folds) ===")
print(resumo.round(4).to_string(), flush=True)

# teste t corrigido de Nadeau-Bengio
def nb(a, b, n_te, n_tr):
    dif = a - b; n = len(dif)
    if dif.std(ddof=1) == 0: return 0.0, 1.0
    t = dif.mean() / np.sqrt(dif.var(ddof=1) * (1/n + n_te/n_tr))
    return t, 2*stats.t.sf(abs(t), n-1)
n_te = len(y)/5; n_tr = len(y)*4/5
F["chave"] = F["conjunto"] + " | " + F["modelo"]
piv = F.pivot(index="fold", columns="chave", values="auc")
ref = resumo.index[0][0] + " | " + resumo.index[0][1]
comp = []
for k in piv.columns:
    if k == ref: continue
    t, p = nb(piv[ref].values, piv[k].values, n_te, n_tr)
    comp.append(dict(referencia=ref, comparado=k, dif_AUC=piv[ref].mean()-piv[k].mean(), t=t, p=p))
C = pd.DataFrame(comp)
C["fdr"] = multipletests(C["p"], method="fdr_bh")[1]
C = C.sort_values("dif_AUC")
C.to_csv(f"{OUT}/M12_teste_pareado_sobrevida.csv", index=False)
print(f"\n=== Teste t corrigido — referencia: {ref} ===")
print(C.round(5).to_string(index=False), flush=True)

# comparacao direta genes x clinico dentro do MESMO algoritmo
print("\n=== Genes vs clinico, dentro de cada algoritmo (teste pareado) ===")
gc = []
for nome in MODELOS:
    a = F[(F.conjunto == "genes (489)") & (F.modelo == nome)].sort_values("fold")["auc"].values
    b = F[(F.conjunto == "clinico (6)") & (F.modelo == nome)].sort_values("fold")["auc"].values
    t, p = nb(a, b, n_te, n_tr)
    gc.append(dict(modelo=nome, AUC_genes=a.mean(), AUC_clinico=b.mean(), dif=a.mean()-b.mean(), p=p))
GC = pd.DataFrame(gc); GC["fdr"] = multipletests(GC["p"], method="fdr_bh")[1]
GC.to_csv(f"{OUT}/M13_genes_vs_clinico.csv", index=False)
print(GC.round(4).to_string(index=False), flush=True)

# ---- teste de DeLong (AUC pareado nas predicoes fora da amostra) ----
def delong(y, p1, p2):
    y = np.asarray(y); pos = y == 1; neg = y == 0
    m, n = pos.sum(), neg.sum()
    def midrank(x):
        s = np.argsort(x); xs = x[s]; N = len(x); T = np.zeros(N); i = 0
        while i < N:
            j = i
            while j < N-1 and xs[j+1] == xs[i]: j += 1
            T[i:j+1] = 0.5*(i+j) + 1; i = j+1
        out = np.empty(N); out[s] = T; return out
    V01 = np.empty((2, m)); V10 = np.empty((2, n)); aucs = np.empty(2)
    for k, p in enumerate([p1, p2]):
        px, nx = p[pos], p[neg]
        tx, ty, tz = midrank(px), midrank(nx), midrank(np.r_[px, nx])
        aucs[k] = (tz[:m].sum() - m*(m+1)/2) / (m*n)
        V01[k] = (tz[:m] - tx) / n
        V10[k] = 1 - (tz[m:] - ty) / m
    S = np.cov(V01)/m + np.cov(V10)/n
    var = S[0,0] + S[1,1] - 2*S[0,1]
    z = (aucs[0]-aucs[1]) / np.sqrt(var) if var > 0 else 0.0
    return aucs[0], aucs[1], z, 2*stats.norm.sf(abs(z))

cv1 = StratifiedKFold(5, shuffle=True, random_state=SEED)
melhores = ["Gradient Boosting", "Random Forest", "Regressao logistica (L2)", "SVM RBF"]
oof = {}
for cj, Xd in CONJUNTOS.items():
    for nome in melhores:
        oof[f"{cj} | {nome}"] = cross_val_predict(clone(MODELOS[nome]), Xd, y, cv=cv1,
                                                  method="predict_proba")[:, 1]
np.save(f"{OUT}/oof_predicoes.npy", np.array([oof[k] for k in oof]))
pd.DataFrame(oof).assign(y=y).to_csv(f"{OUT}/M14_predicoes_oof_sobrevida.csv", index=False)
dl = []
chaves = list(oof)
base = "clinico (6) | Gradient Boosting"
for k in chaves:
    if k == base: continue
    a1, a2, z, p = delong(y, oof[k], oof[base])
    dl.append(dict(modelo=k, AUC=a1, AUC_referencia=a2, dif=a1-a2, z=z, p_DeLong=p))
D = pd.DataFrame(dl); D["fdr"] = multipletests(D["p_DeLong"], method="fdr_bh")[1]
D.to_csv(f"{OUT}/M15_teste_delong.csv", index=False)
print(f"\n=== Teste de DeLong (AUC pareado) — referencia: {base} ===")
print(D.round(4).to_string(index=False), flush=True)

# ---- analise de erro do melhor modelo combinado ----
mkey = resumo.index[0][0] + " | " + resumo.index[0][1]
mkey = mkey if mkey in oof else "combinado | Gradient Boosting"
p = oof[mkey]
d10["risco"] = p
lim = np.median(p)
d10["pred"] = (p > lim).astype(int)
d10["acertou"] = (d10["pred"] == d10["obito10"]).astype(int)
d10[["patient_id", "subtipo", "tumor_stage", "obito10", "risco", "pred", "acertou"]] \
    .to_csv(f"{OUT}/M16_predicoes_paciente_sobrevida.csv", index=False)
print(f"\n=== Analise de erro — {mkey} (corte na mediana) ===")
print("acuracia:", round(d10["acertou"].mean(), 3))

# AUC por estrato + IC bootstrap
def auc_ic(yy, pp, B=2000, seed=SEED):
    rng = np.random.default_rng(seed); a = roc_auc_score(yy, pp); vals = []
    idx = np.arange(len(yy))
    for _ in range(B):
        s = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(yy[s])) < 2: continue
        vals.append(roc_auc_score(yy[s], pp[s]))
    return a, np.percentile(vals, 2.5), np.percentile(vals, 97.5)

estratos = []
for col, vals in [("subtipo", ["LumA","LumB","Her2","Basal","claudin-low","Normal"]),
                  ("estagio", ["I","II","III-IV"])]:
    if col == "estagio":
        d10["estagio"] = d10["tumor_stage"].map({0:"I",1:"I",2:"II",3:"III-IV",4:"III-IV"})
    for v in vals:
        s = d10[d10[col] == v]
        if s["obito10"].nunique() < 2 or len(s) < 60: 
            estratos.append(dict(estrato=f"{col}={v}", n=len(s), eventos=int(s["obito10"].sum()), AUC=np.nan)); continue
        a, lo, hi = auc_ic(s["obito10"].values, s["risco"].values, B=1000)
        estratos.append(dict(estrato=f"{col}={v}", n=len(s), eventos=int(s["obito10"].sum()),
                             AUC=a, IC_inf=lo, IC_sup=hi, acuracia=s["acertou"].mean()))
E = pd.DataFrame(estratos)
E.to_csv(f"{OUT}/M17_auc_por_estrato.csv", index=False)
print("\n=== Onde o modelo acerta e onde falha (AUC com IC 95% bootstrap) ===")
print(E.round(3).to_string(index=False))

# o erro depende do subtipo / estagio?
for col in ["subtipo", "estagio"]:
    s = d10.dropna(subset=[col])
    tab = pd.crosstab(s[col], s["acertou"])
    c, pv, _, _ = stats.chi2_contingency(tab)
    print(f"\nErro depende de {col}? chi2={c:.1f} p={pv:.3g}")
    print(s.groupby(col)["acertou"].agg(["mean","size"]).round(3).to_string())

# quem sao os erros: perfil clinico
perf = []
for var in ["age_at_diagnosis", "tumor_size", "nottingham_prognostic_index",
            "lymph_nodes_examined_positive", "mutation_count", "tempo"]:
    s = d10.dropna(subset=[var])
    u, pv = stats.mannwhitneyu(s.loc[s.acertou==1, var], s.loc[s.acertou==0, var])
    perf.append(dict(variavel=var, mediana_acerto=s.loc[s.acertou==1, var].median(),
                     mediana_erro=s.loc[s.acertou==0, var].median(), U=u, p=pv))
P = pd.DataFrame(perf); P["fdr"] = multipletests(P["p"], method="fdr_bh")[1]
P.to_csv(f"{OUT}/M18_perfil_erros_sobrevida.csv", index=False)
print("\n=== Perfil clinico dos casos errados ===")
print(P.round(4).to_string(index=False))

# falsos positivos vs falsos negativos
d10["tipo_erro"] = np.select(
    [(d10.pred==1)&(d10.obito10==1), (d10.pred==0)&(d10.obito10==0),
     (d10.pred==1)&(d10.obito10==0), (d10.pred==0)&(d10.obito10==1)],
    ["acerto (alto risco, obito)", "acerto (baixo risco, viva)",
     "falso positivo", "falso negativo"], default="")
print("\n=== Tipos de erro ===")
print(d10["tipo_erro"].value_counts().to_string())
print("\nDistribuicao por subtipo:")
print(pd.crosstab(d10["subtipo"], d10["tipo_erro"], normalize="index").round(3).to_string())
d10.groupby("tipo_erro")[["age_at_diagnosis","tumor_size","nottingham_prognostic_index","tempo"]] \
   .median().round(2).to_csv(f"{OUT}/M19_tipos_de_erro.csv")
print("\nMedianas por tipo de erro:")
print(d10.groupby("tipo_erro")[["age_at_diagnosis","tumor_size","nottingham_prognostic_index","tempo"]].median().round(2).to_string())

# calibracao
faixas = pd.cut(d10["risco"], [0,.2,.3,.4,.5,.6,.7,1.0])
cal = d10.groupby(faixas, observed=True).agg(n=("obito10","size"), risco_medio=("risco","mean"),
                                             obito_observado=("obito10","mean"))
cal.to_csv(f"{OUT}/M20_calibracao_sobrevida.csv")
print("\n=== Calibracao (risco previsto vs observado) ===")
print(cal.round(3).to_string())

# separacao de sobrevida por tercil de risco (log-rank)
from lifelines.statistics import multivariate_logrank_test
from lifelines import KaplanMeierFitter
d10["tercil"] = pd.qcut(d10["risco"], 3, labels=["baixo","medio","alto"])
lr = multivariate_logrank_test(d10["tempo"], d10["tercil"], d10["evento"])
print(f"\n=== Tercis de risco previsto — log-rank p = {lr.p_value:.3g} ===")
for t in ["baixo","medio","alto"]:
    s = d10[d10.tercil==t]; k = KaplanMeierFitter().fit(s["tempo"], s["evento"])
    print(f"{t}: n={len(s)} obitos={int(s['evento'].sum())} mediana OS={k.median_survival_time_:.1f} meses")
print("\nOK-TAREFA2", flush=True)
