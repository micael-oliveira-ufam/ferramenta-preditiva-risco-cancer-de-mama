#!/usr/bin/env python3
"""Treino final SEM vazamento: selecao de genes e calibracao dentro da validacao cruzada."""
import numpy as np, pandas as pd, json, warnings, os
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import SelectFromModel
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy import stats
warnings.filterwarnings("ignore")
SEED = 42; np.random.seed(SEED)
OUT = "plataforma"; os.makedirs(OUT, exist_ok=True)

d = pd.read_csv("dados/METABRIC_RNA_Mutation.csv", low_memory=False)
d = d.rename(columns={"pam50_+_claudin-low_subtype": "subtipo"})
d = d[(d["cancer_type"] != "Breast Sarcoma") & (d["subtipo"] != "NC") &
      (d["overall_survival_months"] > 0)].reset_index(drop=True)
clin_end = list(d.columns).index("death_from_cancer")
expr = [c for c in d.columns[clin_end+1:] if not c.endswith("_mut")]
d["evento"] = 1 - d["overall_survival"]; d["tempo"] = d["overall_survival_months"]
cv_ext = StratifiedKFold(5, shuffle=True, random_state=SEED)
cv_int = StratifiedKFold(5, shuffle=True, random_state=SEED + 1)

# ---------- MODELO 1: subtipo ----------
X = d[expr].values; y = d["subtipo"].values
pipe1 = Pipeline([("sc", StandardScaler()),
                  ("lr", LogisticRegression(max_iter=5000, C=0.1, random_state=SEED))])
oof1 = cross_val_predict(clone(pipe1), X, y, cv=cv_ext, method="predict_proba")
pipe1.fit(X, y); classes = list(pipe1.named_steps["lr"].classes_)
pred1 = np.array(classes)[oof1.argmax(1)]
acc1 = float((pred1 == y).mean())
print(f"[subtipo] acuracia fora da amostra: {acc1:.4f}")
conf = oof1.max(1)
tab_conf = pd.DataFrame({"c": conf, "ok": pred1 == y,
                         "f": pd.cut(conf, [0, .4, .5, .6, .7, .8, 1.01])}) \
    .groupby("f", observed=True).agg(n=("ok", "size"), acerto=("ok", "mean"))
print(tab_conf.round(3).to_string())
acc_classe = {c: float((pred1[y == c] == c).mean()) for c in classes}
print("acerto por classe:", {k: round(v, 3) for k, v in acc_classe.items()})

# ---------- MODELO 2: obito em 10 anos ----------
clin_vars = ["age_at_diagnosis", "neoplasm_histologic_grade", "tumor_size",
             "lymph_nodes_examined_positive", "nottingham_prognostic_index"]
d10 = d[(d["tempo"] >= 120) | (d["evento"] == 1)].dropna(subset=clin_vars).reset_index(drop=True)
d10["obito10"] = ((d10["evento"] == 1) & (d10["tempo"] < 120)).astype(int)
yv = d10["obito10"].values
Cl = pd.concat([d10[clin_vars],
                pd.get_dummies(d10["subtipo"]).reindex(columns=classes, fill_value=0)
                  .add_prefix("sub_").astype(float)], axis=1)
clin_feats = list(Cl.columns); Xc = Cl.values.astype(float)
Xg = d10[expr].values
Xk = np.hstack([Xg, Xc])
n_clin = Xc.shape[1]

def faz_pipe(seleciona):
    passos = [("sc", StandardScaler())]
    if seleciona:
        passos.append(("sel", SelectFromModel(
            LogisticRegression(penalty="l1", solver="liblinear", C=0.05,
                               max_iter=5000, random_state=SEED))))
    passos.append(("lr", LogisticRegression(max_iter=5000, C=0.05, random_state=SEED)))
    return Pipeline(passos)

def hl(p, yy):
    g = pd.qcut(pd.Series(p).rank(method="first"), 10)
    t = pd.DataFrame({"p": p, "y": yy}).groupby(g, observed=True).agg(
        n=("y", "size"), obs=("y", "sum"), esp=("p", "sum"))
    s = (((t.obs - t.esp) ** 2) / (t.esp * (1 - t.esp / t.n))).sum()
    return float(s), float(stats.chi2.sf(s, len(t) - 2))

def avalia(Xd, nome, seleciona):
    """AUC/Brier honestos: selecao E calibracao dentro de cada fold externo."""
    aucs, briers, cal_oof = [], [], np.zeros(len(yv))
    for tr, te in cv_ext.split(Xd, yv):
        p = clone(faz_pipe(seleciona))
        p.fit(Xd[tr], yv[tr])
        # calibracao ajustada SO no treino, via validacao interna
        pin = cross_val_predict(clone(p), Xd[tr], yv[tr], cv=cv_int, method="predict_proba")[:, 1]
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.01, y_max=0.99).fit(pin, yv[tr])
        praw = p.predict_proba(Xd[te])[:, 1]
        cal_oof[te] = iso.predict(praw)
        aucs.append(roc_auc_score(yv[te], praw))
        briers.append(brier_score_loss(yv[te], praw))
    auc = float(np.mean(aucs)); br = float(np.mean(briers))
    br_cal = brier_score_loss(yv, cal_oof)
    _, p_hl = hl(cal_oof, yv)
    print(f"{nome}: AUC={auc:.4f} (dp {np.std(aucs):.4f}) | Brier {br:.4f} -> {br_cal:.4f} "
          f"| Hosmer-Lemeshow calibrado p={p_hl:.3g}")
    return dict(auc=auc, auc_dp=float(np.std(aucs)), brier=br, brier_cal=float(br_cal),
                hl_p_cal=p_hl, cal_oof=cal_oof)

print("\n[sobrevida] desempenho honesto (selecao dentro da CV):")
AV = {"clinico": avalia(Xc, "clinico", False),
      "genes": avalia(Xg, "genes", True),
      "combinado": avalia(Xk, "combinado", True)}

# modelo final (treinado em tudo) + curva de calibracao a partir de predicoes fora da amostra
FINAL = {}
for nome, Xd, sel in [("clinico", Xc, False), ("genes", Xg, True), ("combinado", Xk, True)]:
    p = clone(faz_pipe(sel)); p.fit(Xd, yv)
    oof = cross_val_predict(clone(p), Xd, yv, cv=cv_ext, method="predict_proba")[:, 1]
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.01, y_max=0.99).fit(oof, yv)
    if sel:
        mask = p.named_steps["sel"].get_support()
        nomes = [(expr + clin_feats)[i] if nome == "combinado" else expr[i]
                 for i in np.where(mask)[0]]
        sc = p.named_steps["sc"]
        media = sc.mean_[mask]; escala = sc.scale_[mask]
    else:
        nomes = clin_feats; sc = p.named_steps["sc"]
        media = sc.mean_; escala = sc.scale_
    FINAL[nome] = dict(nomes=nomes, media=media, escala=escala,
                       coef=p.named_steps["lr"].coef_[0],
                       inter=float(p.named_steps["lr"].intercept_[0]), iso=iso)
    print(f"  {nome}: {len(nomes)} variaveis no modelo final")

# confiabilidade por subtipo, usando as predicoes calibradas honestas
d10["risco"] = AV["combinado"]["cal_oof"]
conf_estrato = {}
for s in classes:
    sub = d10[d10.subtipo == s]
    if sub.obito10.nunique() < 2: continue
    rng = np.random.default_rng(SEED); vals = []
    yy = sub.obito10.values; pp = sub.risco.values; idx = np.arange(len(yy))
    for _ in range(1500):
        b = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(yy[b])) < 2: continue
        vals.append(roc_auc_score(yy[b], pp[b]))
    a = float(roc_auc_score(yy, pp)); lo = float(np.percentile(vals, 2.5))
    conf_estrato[s] = dict(auc=a, ic_inf=lo, ic_sup=float(np.percentile(vals, 97.5)),
                           n=int(len(sub)), eventos=int(sub.obito10.sum()),
                           confiavel=bool(lo > 0.5))
print("\nConfiabilidade por subtipo (combinado, calibrado, sem vazamento):")
for k, v in conf_estrato.items():
    print(f"  {k}: AUC {v['auc']:.3f} [{v['ic_inf']:.3f}-{v['ic_sup']:.3f}] n={v['n']} "
          f"-> {'confiavel' if v['confiavel'] else 'NAO CONFIAVEL'}")

# estagio
d10["estagio"] = d10["tumor_stage"].map({0: "I", 1: "I", 2: "II", 3: "III-IV", 4: "III-IV"})
conf_estagio = {}
for s in ["I", "II", "III-IV"]:
    sub = d10[d10.estagio == s]
    if len(sub) < 50 or sub.obito10.nunique() < 2: continue
    rng = np.random.default_rng(SEED); vals = []
    yy = sub.obito10.values; pp = sub.risco.values; idx = np.arange(len(yy))
    for _ in range(1500):
        b = rng.choice(idx, len(idx), replace=True)
        if len(np.unique(yy[b])) < 2: continue
        vals.append(roc_auc_score(yy[b], pp[b]))
    lo = float(np.percentile(vals, 2.5))
    conf_estagio[s] = dict(auc=float(roc_auc_score(yy, pp)), ic_inf=lo,
                           ic_sup=float(np.percentile(vals, 97.5)), n=int(len(sub)),
                           eventos=int(sub.obito10.sum()), confiavel=bool(lo > 0.5))
print("\nPor estagio:")
for k, v in conf_estagio.items():
    print(f"  {k}: AUC {v['auc']:.3f} [{v['ic_inf']:.3f}-{v['ic_sup']:.3f}] n={v['n']} "
          f"-> {'confiavel' if v['confiavel'] else 'NAO CONFIAVEL'}")

# tabela de calibracao final
faixas = pd.cut(d10["risco"], [0, .2, .3, .4, .5, .6, .7, 1.01])
cal_tab = d10.groupby(faixas, observed=True).agg(
    n=("obito10", "size"), previsto=("risco", "mean"), observado=("obito10", "mean"))
print("\nCalibracao apos correcao isotonica:")
print(cal_tab.round(3).to_string())

from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
q1, q2 = np.percentile(d10["risco"], [33.333, 66.667])
d10["tercil"] = pd.cut(d10["risco"], [-1, q1, q2, 2], labels=["baixo", "medio", "alto"])
lr = multivariate_logrank_test(d10["tempo"], d10["tercil"], d10["evento"])
tercis = {}
for t in ["baixo", "medio", "alto"]:
    s = d10[d10.tercil == t]; k = KaplanMeierFitter().fit(s["tempo"], s["evento"])
    tercis[t] = dict(n=int(len(s)), obitos=int(s.evento.sum()),
                     mediana_os=float(k.median_survival_time_),
                     obito10_obs=float(s.obito10.mean()))
print(f"\nTercis (log-rank p={lr.p_value:.2e}):",
      {k: round(v['mediana_os'], 1) for k, v in tercis.items()})

# ---------- exportacao ----------
def iso_pts(iso):
    xs = np.linspace(0.001, 0.999, 200)
    return dict(x=[round(float(v), 5) for v in xs],
                y=[round(float(v), 5) for v in iso.predict(xs)])

sc1 = pipe1.named_steps["sc"]; lr1 = pipe1.named_steps["lr"]
modelo = {
    "meta": {"coorte": "METABRIC", "n_subtipo": int(len(y)), "n_sobrevida": int(len(yv)),
             "eventos_10a": int(yv.sum()), "semente": SEED,
             "aviso": "Uso exclusivo em pesquisa e ensino. Nao e dispositivo medico."},
    "subtipo": {"genes": expr,
                "media": [round(float(v), 6) for v in sc1.mean_],
                "escala": [round(float(v), 6) for v in sc1.scale_],
                "classes": classes,
                "coef": [[round(float(v), 6) for v in r] for r in lr1.coef_],
                "intercepto": [round(float(v), 6) for v in lr1.intercept_],
                "acuracia_oof": round(acc1, 4),
                "acerto_por_classe": {k: round(v, 4) for k, v in acc_classe.items()},
                "calibracao_confianca": [{"faixa": str(i), "n": int(r.n),
                                          "acerto": round(float(r.acerto), 3)}
                                         for i, r in tab_conf.iterrows()]},
    "sobrevida": {"variaveis_clinicas": clin_vars, "classes_subtipo": classes,
                  "modelos": {}, "confiabilidade_por_subtipo": conf_estrato,
                  "confiabilidade_por_estagio": conf_estagio,
                  "calibracao": [{"faixa": str(i), "n": int(r.n),
                                  "previsto": round(float(r.previsto), 3),
                                  "observado": round(float(r.observado), 3)}
                                 for i, r in cal_tab.iterrows()],
                  "tercis": {"corte_baixo": round(float(q1), 4),
                             "corte_alto": round(float(q2), 4),
                             "logrank_p": float(lr.p_value), "grupos": tercis}}
}
for nome, f in FINAL.items():
    modelo["sobrevida"]["modelos"][nome] = {
        "variaveis": f["nomes"],
        "media": [round(float(v), 6) for v in f["media"]],
        "escala": [round(float(v), 6) for v in f["escala"]],
        "coef": [round(float(v), 6) for v in f["coef"]],
        "intercepto": round(f["inter"], 6),
        "auc": round(AV[nome]["auc"], 4), "auc_dp": round(AV[nome]["auc_dp"], 4),
        "brier": round(AV[nome]["brier"], 4),
        "brier_calibrado": round(AV[nome]["brier_cal"], 4),
        "hosmer_p_calibrado": AV[nome]["hl_p_cal"],
        "calibracao": iso_pts(f["iso"])}
with open(f"{OUT}/modelos.json", "w") as fh:
    json.dump(modelo, fh, separators=(",", ":"))
print(f"\nmodelos.json: {os.path.getsize(f'{OUT}/modelos.json')/1024:.0f} KB")

# ---------- pacientes-exemplo ----------
def monta(r, rot, desc):
    return dict(rotulo=rot, descricao=desc,
                clinico={"age_at_diagnosis": float(r.age_at_diagnosis),
                         "neoplasm_histologic_grade": float(r.neoplasm_histologic_grade),
                         "tumor_size": float(r.tumor_size),
                         "lymph_nodes_examined_positive": float(r.lymph_nodes_examined_positive),
                         "nottingham_prognostic_index": float(r.nottingham_prognostic_index),
                         "subtipo": str(r.subtipo)},
                estagio=(None if pd.isna(r.tumor_stage) else int(r.tumor_stage)),
                er=str(r.er_status), her2=str(r.her2_status),
                desfecho={"obito_10a": int(r.obito10), "tempo_meses": round(float(r.tempo), 1),
                          "vital": ("obito" if r.evento == 1 else "viva")},
                expressao={g: round(float(r[g]), 4) for g in expr})

base = d10.dropna(subset=["neoplasm_histologic_grade", "tumor_size"])
cands = [
    ((base.subtipo == "LumA") & (base.tumor_stage == 1) & (base.neoplasm_histologic_grade == 1)
     & (base.lymph_nodes_examined_positive == 0) & (base.obito10 == 0),
     "Luminal A, estágio I, axila negativa",
     "Baixo risco clássico: tumor bem diferenciado, sem linfonodos, sobrevida longa."),
    ((base.subtipo == "LumB") & (base.tumor_stage == 2) & (base.lymph_nodes_examined_positive >= 2)
     & (base.obito10 == 1), "Luminal B, estágio II, axila positiva",
     "Risco intermediário-alto: proliferação elevada com comprometimento nodal."),
    ((base.subtipo == "Basal") & (base.neoplasm_histologic_grade == 3) & (base.obito10 == 1),
     "Basal, grau 3", "Subtipo de menor confiabilidade do modelo — veja o aviso na tela."),
    ((base.subtipo == "Her2") & (base.obito10 == 1) & (base.tumor_size >= 25),
     "HER2-enriquecido, tumor grande",
     "Coorte anterior ao trastuzumabe: o desfecho não reflete a terapia atual."),
    ((base.subtipo == "claudin-low") & (base.obito10 == 0), "Claudin-low, evolução favorável",
     "Assinatura imune e estromal, com perda de identidade epitelial."),
    ((base.subtipo == "LumA") & (base.obito10 == 1) & (base.age_at_diagnosis < 55)
     & (base.lymph_nodes_examined_positive == 0) & (base.tumor_size <= 22),
     "Luminal A jovem, aparência favorável, óbito em 10 anos",
     "O caso difícil: perfil clínico tranquilizador e desfecho ruim. É onde o modelo mais falha."),
]
exemplos = []
for cond, rot, desc in cands:
    s = base[cond]
    if len(s) == 0: print("sem candidato:", rot); continue
    e = monta(s.iloc[0], rot, desc); exemplos.append(e)
    print(f"exemplo: {rot} -> {e['desfecho']}")
with open(f"{OUT}/exemplos.json", "w") as fh:
    json.dump(exemplos, fh, separators=(",", ":"))
linhas = [{"rotulo": e["rotulo"], **e["clinico"], "estagio": e["estagio"],
           "er_status": e["er"], "her2_status": e["her2"],
           "desfecho_obito_10a": e["desfecho"]["obito_10a"],
           "seguimento_meses": e["desfecho"]["tempo_meses"], **e["expressao"]} for e in exemplos]
pd.DataFrame(linhas).to_csv(f"{OUT}/exemplos_entrada.csv", index=False)
print(f"exemplos.json e exemplos_entrada.csv gravados ({len(exemplos)} pacientes)")
print("\nOK-EXPORT")
