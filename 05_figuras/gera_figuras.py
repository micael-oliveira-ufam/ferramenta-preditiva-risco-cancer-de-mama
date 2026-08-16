#!/usr/bin/env python3
"""Figuras de validacao e tabela consolidada de genes de sobrevida.
Consome os resultados ja gravados pelas etapas anteriores do pipeline."""
import json, os, glob
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "figuras"; os.makedirs(OUT, exist_ok=True)
TAB = "tabelas_finais"; os.makedirs(TAB, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
AZUL, TEAL, AMBAR, VERM, ROXO, CINZA = "#0284C7", "#0D9488", "#D97706", "#DC2626", "#7C3AED", "#64748B"
CORES = {"logistica": CINZA, "random_forest": TEAL, "extra_trees": AZUL,
         "gradient_boosting": AMBAR, "xgboost": VERM, "mlp": ROXO}
NOMES = {"logistica": "Regressão logística", "random_forest": "Random Forest",
         "extra_trees": "Extra Trees", "gradient_boosting": "Gradient Boosting",
         "xgboost": "XGBoost", "mlp": "Rede neural (MLP)"}
CONJ = {"clinico": "Clínico (11 var.)", "genes": "Expressão gênica (489 var.)",
        "combinado": "Combinado (500 var.)"}

M = json.load(open("plataforma2/modelos.json"))
ALG, ROC, SB = M["algoritmos"], M["roc"], M["sobrevida_base"]

# ---------- F1: curvas ROC por conjunto ----------
fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.9))
for ax, cj in zip(axes, ["clinico", "genes", "combinado"]):
    ax.plot([0, 1], [0, 1], ls="--", lw=.9, color="#94A3B8")
    chaves = sorted([c for c in ROC if c.startswith(cj + "|")],
                    key=lambda c: -ALG[c]["auc"])
    for ch in chaves:
        k = ch.split("|")[1]; r = ROC[ch]
        ax.plot(r["fpr"], r["tpr"], color=CORES[k], lw=1.6,
                label=f"{NOMES[k]} — {ALG[ch]['auc']:.3f}")
    ax.set_title(CONJ[cj], fontsize=10, fontweight="bold")
    ax.set_xlabel("1 − especificidade"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=6.6, loc="lower right", frameon=False)
axes[0].set_ylabel("sensibilidade")
fig.suptitle("Curvas ROC — validação cruzada 5-fold, óbito em 10 anos (n = 1.560)",
             fontsize=11, fontweight="bold", y=1.03)
fig.savefig(f"{OUT}/F01_curvas_roc.png"); plt.close(fig)

# ---------- F2: AUC com IC 95% (todos os 16 modelos) ----------
linhas = sorted(ALG.items(), key=lambda kv: kv[1]["auc"])
fig, ax = plt.subplots(figsize=(7.2, 5.6))
for i, (ch, m) in enumerate(linhas):
    k = ch.split("|")[1]
    ax.plot([m["ic_inf"], m["ic_sup"]], [i, i], color=CORES[k], lw=2.4, alpha=.45,
            solid_capstyle="round")
    ax.plot(m["auc"], i, "o", color=CORES[k], ms=5.5)
ax.axvline(.5, color="#94A3B8", ls="--", lw=1)
ax.set_yticks(range(len(linhas)))
ax.set_yticklabels([f"{CONJ[m['conjunto']].split(' (')[0]} · {NOMES[ch.split('|')[1]]}"
                    for ch, m in linhas], fontsize=8)
ax.set_xlabel("AUC (IC 95%, DeLong)"); ax.set_xlim(.45, .82)
ax.set_title("Desempenho dos 16 modelos — predição de óbito em 10 anos",
             fontsize=10.5, fontweight="bold")
fig.savefig(f"{OUT}/F02_auc_ic95.png"); plt.close(fig)

# ---------- F3: calibração ----------
cal = pd.DataFrame(SB["calibracao"])
fig, ax = plt.subplots(figsize=(4.6, 4.3))
ax.plot([0, 1], [0, 1], ls="--", lw=.9, color="#94A3B8", label="calibração perfeita")
ax.plot(cal["previsto"], cal["observado"], "o-", color=AZUL, lw=1.6, ms=6,
        label="após correção isotônica")
for _, r in cal.iterrows():
    ax.annotate(f"n={r['n']}", (r["previsto"], r["observado"]), fontsize=6,
                xytext=(4, -8), textcoords="offset points", color=CINZA)
ax.set_xlabel("risco previsto"); ax.set_ylabel("óbitos observados")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title("Calibração do modelo em produção\nHosmer-Lemeshow p = 0,66",
             fontsize=10, fontweight="bold")
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
fig.savefig(f"{OUT}/F03_calibracao.png"); plt.close(fig)

# ---------- F4: confiabilidade por subtipo e estagio ----------
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6),
                         gridspec_kw={"width_ratios": [1.35, 1]})
for ax, chave, titulo in [(axes[0], "confiabilidade_por_subtipo", "Por subtipo molecular"),
                          (axes[1], "confiabilidade_por_estagio", "Por estágio tumoral")]:
    d = SB[chave]; ks = list(d)
    for i, k in enumerate(ks):
        o = d[k]; cor = TEAL if o["confiavel"] else VERM
        ax.plot([o["ic_inf"], o["ic_sup"]], [i, i], color=cor, lw=6, alpha=.32,
                solid_capstyle="round")
        ax.plot(o["auc"], i, "o", color=cor, ms=6)
        ax.text(.83, i, f"n={o['n']}", va="center", fontsize=7, color=CINZA)
    ax.axvline(.5, color="#475569", ls="--", lw=1.1)
    ax.set_yticks(range(len(ks))); ax.set_yticklabels(ks, fontsize=8.5)
    ax.set_xlim(.42, .9); ax.set_xlabel("AUC (IC 95% bootstrap)")
    ax.set_title(titulo, fontsize=10, fontweight="bold")
fig.legend(handles=[Patch(color=TEAL, alpha=.5, label="IC não cruza o acaso"),
                    Patch(color=VERM, alpha=.5, label="IC inclui 0,5 — não confiável")],
           fontsize=7.5, frameon=False, loc="lower center", ncol=2, bbox_to_anchor=(.5, -.11))
fig.suptitle("Onde a predição é confiável — modelo combinado calibrado",
             fontsize=10.5, fontweight="bold", y=1.02)
fig.savefig(f"{OUT}/F04_confiabilidade_estratos.png"); plt.close(fig)

# ---------- F5: tercis de risco ----------
t = SB["tercis"]["grupos"]
fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
gs = ["baixo", "medio", "alto"]; cs = [TEAL, AMBAR, VERM]
axes[0].bar(range(3), [t[g]["obito10_obs"]*100 for g in gs], color=cs, width=.62)
for i, g in enumerate(gs):
    axes[0].text(i, t[g]["obito10_obs"]*100+1.5, f"{t[g]['obito10_obs']*100:.1f}%",
                 ha="center", fontsize=8.5, fontweight="bold")
axes[0].set_xticks(range(3)); axes[0].set_xticklabels(["Baixo", "Médio", "Alto"])
axes[0].set_ylabel("óbitos em 10 anos (%)"); axes[0].set_ylim(0, 75)
axes[0].set_title("Desfecho observado por tercil", fontsize=9.5, fontweight="bold")
axes[1].bar(range(3), [t[g]["mediana_os"] for g in gs], color=cs, width=.62)
for i, g in enumerate(gs):
    axes[1].text(i, t[g]["mediana_os"]+5, f"{t[g]['mediana_os']:.0f}", ha="center",
                 fontsize=8.5, fontweight="bold")
axes[1].set_xticks(range(3)); axes[1].set_xticklabels(["Baixo", "Médio", "Alto"])
axes[1].set_ylabel("sobrevida global mediana (meses)"); axes[1].set_ylim(0, 250)
axes[1].set_title("Sobrevida real por tercil de risco previsto", fontsize=9.5, fontweight="bold")
fig.suptitle(f"Estratificação em tercis — log-rank p = {SB['tercis']['logrank_p']:.1e}",
             fontsize=10.5, fontweight="bold", y=1.04)
fig.savefig(f"{OUT}/F05_tercis_risco.png"); plt.close(fig)

# ---------- genes de sobrevida ----------
os_ = pd.read_csv("saidas/tabelas/T12_cox_univariado_OS_todos_genes.csv")
dss = pd.read_csv("saidas/tabelas/T13_cox_univariado_DSS_todos_genes.csv")
cols = {c.lower(): c for c in os_.columns}
def col(df, *cands):
    for c in cands:
        for orig in df.columns:
            if orig.lower() == c: return orig
    raise KeyError(cands)
g_os = col(os_, "gene"); hr_os = col(os_, "hr"); fdr_os = col(os_, "fdr", "fdr_bh")
li = col(os_, "ic95_inf"); ls = col(os_, "ic95_sup")
g_d = col(dss, "gene"); hr_d = col(dss, "hr"); fdr_d = col(dss, "fdr", "fdr_bh")

base = os_[[g_os, hr_os, li, ls, fdr_os]].copy()
base.columns = ["gene", "HR_os", "ic_inf", "ic_sup", "fdr_os"]
base = base.merge(dss[[g_d, hr_d, fdr_d]].rename(
    columns={g_d: "gene", hr_d: "HR_dss", fdr_d: "fdr_dss"}), on="gene", how="left")
base["direcao"] = np.where(base.HR_os < 1, "maior sobrevida (protetor)", "menor sobrevida (risco)")
base["significativo_os"] = base.fdr_os < 0.05
sig = base[base.significativo_os].copy()
sig["forca"] = np.abs(np.log(sig.HR_os))
sig = sig.sort_values("forca", ascending=False)
sig.to_csv(f"{TAB}/G01_genes_sobrevida_significativos.csv", index=False)
prot = sig[sig.HR_os < 1].head(20); risco = sig[sig.HR_os > 1].head(20)
prot.to_csv(f"{TAB}/G02_genes_maior_sobrevida.csv", index=False)
risco.to_csv(f"{TAB}/G03_genes_menor_sobrevida.csv", index=False)
print(f"genes significativos (FDR<0,05): {len(sig)} "
      f"({(sig.HR_os<1).sum()} protetores, {(sig.HR_os>1).sum()} de risco)")

# ---------- F6: forest plot dos genes ----------
sel = pd.concat([risco.head(12).sort_values("HR_os"),
                 prot.head(12).sort_values("HR_os")])
fig, ax = plt.subplots(figsize=(6.6, 7.4))
for i, (_, r) in enumerate(sel.iterrows()):
    cor = VERM if r.HR_os > 1 else TEAL
    ax.plot([r.ic_inf, r.ic_sup], [i, i], color=cor, lw=2, alpha=.5, solid_capstyle="round")
    ax.plot(r.HR_os, i, "o", color=cor, ms=5)
ax.axvline(1, color="#475569", ls="--", lw=1.1)
ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel.gene.str.upper(), fontsize=7.6)
ax.set_xlabel("Hazard ratio por +1 desvio-padrão de expressão (sobrevida global)")
ax.set_xlim(0.68, 1.42)
ax.set_title("Genes associados à sobrevida\n(Cox univariado, 489 genes, FDR < 0,05)",
             fontsize=10.5, fontweight="bold")
ax.legend(handles=[Patch(color=VERM, alpha=.6, label="HR > 1 — menor sobrevida"),
                   Patch(color=TEAL, alpha=.6, label="HR < 1 — maior sobrevida")],
          fontsize=7.5, frameon=False, loc="lower right")
fig.savefig(f"{OUT}/F06_genes_forest.png"); plt.close(fig)

# ---------- F7: importancia no RF de desfecho + genes por estagio ----------
imp = pd.read_csv("saidas_complementares/C05_rf_mortalidade10a_importancia.csv").head(15)
est = pd.read_csv("saidas_complementares/C01_cox_por_estagio_todos_genes.csv")
est2 = est[(est.estagio == "II") & (est.fdr < 0.05)].sort_values("p").head(12)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
axes[0].barh(range(len(imp))[::-1], imp.importancia_gini, color=AZUL, alpha=.85)
axes[0].set_yticks(range(len(imp))[::-1]); axes[0].set_yticklabels(imp.gene.str.upper(), fontsize=7.5)
axes[0].set_xlabel("importância (índice de Gini)")
axes[0].set_title("Random Forest — predição de óbito em 10 anos", fontsize=9.5, fontweight="bold")
c2 = [VERM if h > 1 else TEAL for h in est2.HR]
axes[1].barh(range(len(est2))[::-1], est2.HR - 1, left=1, color=c2, alpha=.85)
axes[1].axvline(1, color="#475569", lw=1)
axes[1].set_yticks(range(len(est2))[::-1]); axes[1].set_yticklabels(est2.gene.str.upper(), fontsize=7.5)
axes[1].set_xlabel("hazard ratio"); axes[1].set_xlim(0.7, 1.4)
axes[1].set_title("Estágio II — onde o sinal gênico se concentra", fontsize=9.5, fontweight="bold")
fig.savefig(f"{OUT}/F07_genes_ml_e_estagio.png"); plt.close(fig)

# ---------- F8: matriz de confusao do classificador de subtipo ----------
cm = pd.read_csv("saidas_ml/M04_matriz_confusao_melhor.csv", index_col=0)
prop = cm.div(cm.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(5.4, 4.6))
im = ax.imshow(prop.values, cmap="Blues", vmin=0, vmax=1)
ax.set_xticks(range(len(cm.columns))); ax.set_xticklabels(cm.columns, rotation=45, ha="right", fontsize=8)
ax.set_yticks(range(len(cm.index))); ax.set_yticklabels(cm.index, fontsize=8)
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        v = cm.values[i, j]
        if v:
            ax.text(j, i, str(v), ha="center", va="center", fontsize=7.5,
                    color="white" if prop.values[i, j] > .55 else "#0F172A")
ax.set_xlabel("subtipo predito"); ax.set_ylabel("subtipo real")
ax.set_title("Classificação dos subtipos — Gradient Boosting\nacurácia 78,7% (erro OOB 23,9% no RF)",
             fontsize=9.5, fontweight="bold")
fig.colorbar(im, ax=ax, shrink=.8, label="proporção da linha")
fig.savefig(f"{OUT}/F08_matriz_confusao_subtipo.png"); plt.close(fig)

# ---------- tabela final de metricas ----------
met = pd.DataFrame([{"conjunto": m["conjunto"], "algoritmo": m["algoritmo"],
                     "n_variaveis": len(m["variaveis"]), "AUC": m["auc"],
                     "IC95_inf": m["ic_inf"], "IC95_sup": m["ic_sup"],
                     "Brier": m["brier"], "Brier_calibrado": m["brier_calibrado"],
                     "Hosmer_Lemeshow_p": m["hosmer_p"],
                     "erro_reimplementacao_js": m["erro_verificacao"],
                     "hiperparametros": json.dumps(m["hiper"], ensure_ascii=False)}
                    for m in ALG.values()]).sort_values("AUC", ascending=False)
met.to_csv(f"{TAB}/G04_metricas_16_modelos.csv", index=False)
print(met[["conjunto", "algoritmo", "AUC", "IC95_inf", "IC95_sup"]].head(5).to_string(index=False))
print("\nfiguras:", sorted(os.listdir(OUT)))
