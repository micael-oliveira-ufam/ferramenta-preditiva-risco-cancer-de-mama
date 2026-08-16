#!/usr/bin/env python3
"""(B) Interacao gene x tratamento — leitura farmacogenomica."""
import numpy as np, pandas as pd, warnings, os, time
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
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

t0 = time.time()
inter = []
for trat in ["hormone_therapy", "chemotherapy", "radio_therapy"]:
    base = d[["tempo", "evento_os", trat, "age_at_diagnosis",
              "neoplasm_histologic_grade", "tumor_size"] + expr_cols].dropna(
              subset=["tempo", "evento_os", trat, "age_at_diagnosis",
                      "neoplasm_histologic_grade", "tumor_size"])
    for i, g in enumerate(expr_cols):
        sub = base[["tempo", "evento_os", g, trat, "age_at_diagnosis",
                    "neoplasm_histologic_grade", "tumor_size"]].copy()
        sub = sub.rename(columns={g: "gene", trat: "trat"})
        sub["gene_x_trat"] = sub["gene"] * sub["trat"]
        try:
            cph = CoxPHFitter().fit(sub, "tempo", "evento_os")
            s = cph.summary
            inter.append(dict(tratamento=trat, gene=g,
                              HR_interacao=float(np.exp(s.loc["gene_x_trat", "coef"])),
                              p_interacao=float(s.loc["gene_x_trat", "p"]),
                              HR_gene_referencia=float(np.exp(s.loc["gene", "coef"])),
                              n=len(sub), eventos=int(sub["evento_os"].sum())))
        except Exception:
            continue
    print(f"{trat}: {len(inter)} ajustes acumulados, {time.time()-t0:.0f}s", flush=True)

B = pd.DataFrame(inter)
B["fdr"] = np.nan
for t in B["tratamento"].unique():
    m = B["tratamento"] == t
    B.loc[m, "fdr"] = multipletests(B.loc[m, "p_interacao"], method="fdr_bh")[1]
B = B.sort_values(["tratamento", "fdr"])
B.to_csv(f"{OUT}/C02_interacao_gene_tratamento.csv", index=False)
print("\n=== interacoes com FDR<0.05 ===", flush=True)
print(B[B.fdr < 0.05].groupby("tratamento").size())
for t in B["tratamento"].unique():
    print(f"\n--- {t}: top 12 ---")
    print(B[B.tratamento == t].head(12)[["gene", "HR_interacao", "HR_gene_referencia",
                                         "p_interacao", "fdr"]].to_string(index=False))

def cox_uni(df, gene):
    sub = df[["tempo", "evento_os", gene]].dropna()
    if sub["evento_os"].sum() < 10: return None
    try:
        cph = CoxPHFitter().fit(sub, "tempo", "evento_os")
    except Exception:
        return None
    s = cph.summary
    return dict(HR=float(np.exp(s.loc[gene, "coef"])),
                ic_inf=float(np.exp(s.loc[gene, "coef lower 95%"])),
                ic_sup=float(np.exp(s.loc[gene, "coef upper 95%"])),
                p=float(s.loc[gene, "p"]), n=len(sub), eventos=int(sub["evento_os"].sum()))

nucleo = [g for g in ["gsk3b","stat5a","aurka","vegfa","bcl2","spry2","igf1","abcb1",
                      "flt3","cdkn2c","pdgfra","erbb2","esr1","mapt","ar","top2a","tubb4a",
                      "nrip1","ncoa3","jak2"] if g in expr_cols]
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
print("\n=== HR por estrato de tratamento (genes-nucleo) ===")
print(E.pivot_table(index="gene", columns=["tratamento","estrato"], values="HR").round(3).to_string())
print("\nn por estrato:")
print(E.groupby(["tratamento","estrato"])[["n","eventos"]].max().to_string())
print("OK-B", flush=True)
