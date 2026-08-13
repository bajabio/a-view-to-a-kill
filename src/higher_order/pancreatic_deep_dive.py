import numpy as np, pandas as pd
from scipy import stats
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
lines=M.index; panc=(lin.reindex(lines)=="Pancreas").values
G=M.values.astype(np.float32); cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1])
genes=np.array(M.columns)
# 1) genome-wide pancreatic-selective dependencies (t-test panc vs rest)
t=np.array([stats.ttest_ind(G[panc,j],G[~panc,j],equal_var=False).statistic for j in range(G.shape[1])])
diff=G[panc].mean(0)-G[~panc].mean(0); lm=G[panc].mean(0)
d=pd.DataFrame({"gene":genes,"t":t,"diff":diff,"panc_mean":lm})
sel=d[(d["diff"]<-0.1)&(d.panc_mean<-0.4)].sort_values("t").head(15)
print("[A] Genome-wide PANCREATIC-SELECTIVE dependencies (more essential in pancreatic vs other cancers):")
print(sel[["gene","panc_mean","diff","t"]].to_string(index=False,float_format=lambda v:f"{v:.2f}"))

# 2) passenger-paralog scan: strong paralog SL where the PARTNER is frequently lost in pancreatic
ml=pd.read_csv("paralog_ml_v2.csv")   # A=partner lost, B=target, diff = B more essential when A lost
E=pd.read_csv("uni_expr.csv",index_col=0); El=E.reindex(lines)
ecol=set(El.columns); glo={g:E[g].quantile(0.25) for g in El.columns}   # per-gene global low threshold
def panc_low_frac(A):
    if A not in ecol: return np.nan
    v=El[A].values[panc]; return np.nanmean(v<glo[A])
gidx={g:i for i,g in enumerate(genes)}
ml=ml[ml.A.isin(ecol)&ml.B.isin(gidx)]
ml["panc_partner_loss"]=ml.A.map(panc_low_frac)
ml["B_panc_dep"]=ml.B.map(lambda b: G[panc,gidx[b]].mean())
cand=ml[(ml["diff"]<-0.1)&(ml.panc_partner_loss>=0.25)&(ml.B_panc_dep<-0.4)].copy()
cand=cand.sort_values("diff").drop_duplicates("B")
print("\n[B] Genome-wide PANCREATIC passenger-paralog targets (partner often lost in pancreatic, target essential):")
print(f"{'target(B)':10}{'partner(A)':11}{'A lost in panc':>14}{'B dep(panc)':>12}{'SL diff':>9}")
for _,r in cand.head(14).iterrows():
    print(f"{r.B:10}{r.A:11}{r.panc_partner_loss*100:12.0f}% {r.B_panc_dep:11.2f}{r['diff']:9.2f}")
cand.to_csv("/home/jmilton/ml/ppset/pancreatic_paralog_targets.csv",index=False)
