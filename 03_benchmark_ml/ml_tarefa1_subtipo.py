#!/usr/bin/env python3
"""Benchmark multi-modelo — TAREFA 1: classificacao dos 6 subtipos moleculares.
Compara 10 algoritmos, com testes estatisticos formais e analise de erro."""
import numpy as np, pandas as pd, warnings, os, json, time
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
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, f1_score,
                             confusion_matrix, classification_report, brier_score_loss)
from scipy import stats
warnings.filterwarnings("ignore")
SEED = 42; np.random.seed(SEED)
OUT = "saidas_ml"; os.makedirs(OUT, exist_ok=True)

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})
d = d[(d["cancer_type"] != "Breast Sarcoma") & (d["subtipo"] != "NC") &
      (d["overall_survival_months"] > 0)].reset_index(drop=True)
clin_end = list(d.columns).index("death_from_cancer")
expr = [c for c in d.columns[clin_end+1:] if not c.endswith("_mut")]
X = d[expr].values; y = d["subtipo"].values
print("X:", X.shape, "| classes:", pd.Series(y).value_counts().to_dict(), flush=True)

MODELOS = {
    "Baseline (classe majoritaria)": DummyClassifier(strategy="most_frequent"),
    "Regressao logistica (L2)": make_pipeline(StandardScaler(), LogisticRegression(max_iter=3000, C=0.1, random_state=SEED)),
    "SVM linear": make_pipeline(StandardScaler(), SVC(kernel="linear", C=0.1, probability=False, random_state=SEED)),
    "SVM RBF": make_pipeline(StandardScaler(), SVC(kernel="rbf", C=10, gamma="scale", random_state=SEED)),
    "Random Forest": RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1),
    "Extra Trees": ExtraTreesClassifier(n_estimators=500, random_state=SEED, n_jobs=-1),
    "Gradient Boosting": HistGradientBoostingClassifier(max_iter=300, random_state=SEED),
    "k-NN (k=15)": make_pipeline(StandardScaler(), KNeighborsClassifier(n_neighbors=15)),
    "LDA": make_pipeline(StandardScaler(), LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")),
    "Naive Bayes": make_pipeline(StandardScaler(), GaussianNB()),
    "Rede neural (MLP)": make_pipeline(StandardScaler(), MLPClassifier(hidden_layer_sizes=(128,), max_iter=800, random_state=SEED)),
}

rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=3, random_state=SEED)
folds = list(rskf.split(X, y))
res_fold = []
t0 = time.time()
for nome, mod in MODELOS.items():
    for i, (tr, te) in enumerate(folds):
        from sklearn.base import clone
        m = clone(mod); m.fit(X[tr], y[tr]); p = m.predict(X[te])
        res_fold.append(dict(modelo=nome, fold=i,
                             acuracia=accuracy_score(y[te], p),
                             acuracia_balanceada=balanced_accuracy_score(y[te], p),
                             f1_macro=f1_score(y[te], p, average="macro")))
    print(f"{nome}: ok ({time.time()-t0:.0f}s)", flush=True)
F = pd.DataFrame(res_fold)
F.to_csv(f"{OUT}/M01_desempenho_por_fold_subtipo.csv", index=False)

resumo = F.groupby("modelo").agg(
    acuracia_media=("acuracia", "mean"), acuracia_dp=("acuracia", "std"),
    balanceada_media=("acuracia_balanceada", "mean"), balanceada_dp=("acuracia_balanceada", "std"),
    f1_macro_media=("f1_macro", "mean"), f1_macro_dp=("f1_macro", "std")
).sort_values("acuracia_media", ascending=False)
resumo.to_csv(f"{OUT}/M02_resumo_modelos_subtipo.csv")
print("\n=== RESUMO (15 folds: 5-fold x 3 repeticoes) ===")
print(resumo.round(4).to_string(), flush=True)

# ---- teste t corrigido de Nadeau-Bengio (comparacao pareada entre modelos) ----
def nb_ttest(a, b, n_test, n_train):
    dif = a - b; n = len(dif)
    if dif.std(ddof=1) == 0: return 0.0, 1.0
    var = dif.var(ddof=1)
    t = dif.mean() / np.sqrt(var * (1/n + n_test/n_train))
    p = 2 * stats.t.sf(abs(t), n - 1)
    return t, p

n_te = len(y)/5; n_tr = len(y)*4/5
piv = F.pivot(index="fold", columns="modelo", values="acuracia")
melhor = resumo.index[0]
comp = []
for m in piv.columns:
    if m == melhor: continue
    t, p = nb_ttest(piv[melhor].values, piv[m].values, n_te, n_tr)
    comp.append(dict(referencia=melhor, comparado=m,
                     dif_media=piv[melhor].mean()-piv[m].mean(), t=t, p_corrigido=p))
C = pd.DataFrame(comp).sort_values("dif_media")
C.to_csv(f"{OUT}/M03_teste_pareado_subtipo.csv", index=False)
print(f"\n=== Teste t corrigido (Nadeau-Bengio) — referencia: {melhor} ===")
print(C.round(5).to_string(index=False), flush=True)

# ---- predicoes fora da amostra do melhor modelo + McNemar contra o 2o ----
from sklearn.base import clone
cv1 = StratifiedKFold(5, shuffle=True, random_state=SEED)
pred = {}
for nome in list(resumo.index[:3]):
    pred[nome] = cross_val_predict(clone(MODELOS[nome]), X, y, cv=cv1, n_jobs=1)
segundo = resumo.index[1]
a = pred[melhor] == y; b = pred[segundo] == y
n01 = int(((~a) & b).sum()); n10 = int((a & (~b)).sum())
from statsmodels.stats.contingency_tables import mcnemar
mc = mcnemar([[int((a & b).sum()), n10], [n01, int(((~a) & (~b)).sum())]], exact=False, correction=True)
print(f"\n=== McNemar: {melhor} vs {segundo} ===")
print(f"so {melhor} acerta: {n10} | so {segundo} acerta: {n01} | chi2={mc.statistic:.3f} p={mc.pvalue:.3g}", flush=True)

# ---- analise de erro do melhor modelo ----
yp = pred[melhor]
classes = ["LumA", "LumB", "Her2", "Basal", "claudin-low", "Normal"]
cm = pd.DataFrame(confusion_matrix(y, yp, labels=classes), index=classes, columns=classes)
cm.to_csv(f"{OUT}/M04_matriz_confusao_melhor.csv")
print(f"\n=== Matriz de confusao — {melhor} (linhas = real) ===")
print(cm.to_string())
rep = pd.DataFrame(classification_report(y, yp, labels=classes, output_dict=True)).T
rep.to_csv(f"{OUT}/M05_metricas_por_classe.csv")
print("\n=== Metricas por classe ===")
print(rep.round(3).to_string())

d["acertou"] = (yp == y).astype(int); d["pred"] = yp
d[["patient_id", "subtipo", "pred", "acertou"]].to_csv(f"{OUT}/M06_predicoes_paciente.csv", index=False)

# teste: erro depende do subtipo real?
tab = pd.crosstab(d["subtipo"], d["acertou"])
chi2, p, _, _ = stats.chi2_contingency(tab)
print(f"\nErro depende do subtipo real? chi2={chi2:.1f} p={p:.3g}")

# teste: erro associado a celularidade / grau / pureza?
testes = []
for var in ["cellularity", "neoplasm_histologic_grade", "er_status", "her2_status",
            "tumor_other_histologic_subtype", "integrative_cluster"]:
    sub = d.dropna(subset=[var])
    t = pd.crosstab(sub[var], sub["acertou"])
    if t.shape[0] < 2: continue
    c, pv, _, _ = stats.chi2_contingency(t)
    testes.append(dict(variavel=var, teste="qui-quadrado", estatistica=c, p=pv, n=len(sub)))
for var in ["age_at_diagnosis", "tumor_size", "mutation_count", "nottingham_prognostic_index",
            "overall_survival_months"]:
    sub = d.dropna(subset=[var])
    u, pv = stats.mannwhitneyu(sub.loc[sub.acertou == 1, var], sub.loc[sub.acertou == 0, var])
    testes.append(dict(variavel=var, teste="Mann-Whitney", estatistica=u, p=pv, n=len(sub),
                       mediana_acerto=sub.loc[sub.acertou == 1, var].median(),
                       mediana_erro=sub.loc[sub.acertou == 0, var].median()))
T = pd.DataFrame(testes)
from statsmodels.stats.multitest import multipletests
T["fdr"] = multipletests(T["p"], method="fdr_bh")[1]
T.to_csv(f"{OUT}/M07_perfil_dos_erros.csv", index=False)
print("\n=== O que caracteriza os casos que o modelo erra? ===")
print(T.round(5).to_string(index=False))

# celularidade detalhada
if "cellularity" in d.columns:
    print("\nTaxa de acerto por celularidade:")
    print(d.groupby("cellularity")["acertou"].agg(["mean", "size"]).round(3).to_string())
print("\nTaxa de acerto por subtipo:")
print(d.groupby("subtipo")["acertou"].agg(["mean", "size"]).round(3).to_string())

# confianca da predicao (modelo probabilistico) — usa RF sempre
rfp = cross_val_predict(RandomForestClassifier(n_estimators=500, random_state=SEED, n_jobs=-1),
                        X, y, cv=cv1, method="predict_proba")
conf = rfp.max(axis=1)
rf_pred = np.array(sorted(np.unique(y)))[rfp.argmax(axis=1)]
ok = rf_pred == y
u, pv = stats.mannwhitneyu(conf[ok], conf[~ok])
print(f"\nConfianca do RF — acertos mediana {np.median(conf[ok]):.3f} vs erros {np.median(conf[~ok]):.3f} "
      f"(Mann-Whitney p={pv:.3g})")
faixas = pd.cut(conf, [0, .4, .5, .6, .7, .8, .9, 1.0])
cal = pd.DataFrame({"conf": conf, "ok": ok, "faixa": faixas}).groupby("faixa", observed=True).agg(
    n=("ok", "size"), acerto_observado=("ok", "mean"), confianca_media=("conf", "mean"))
cal.to_csv(f"{OUT}/M08_calibracao_confianca.csv")
print("\n=== Calibracao: o modelo sabe quando nao sabe? ===")
print(cal.round(3).to_string())
print("\nOK-TAREFA1", flush=True)
