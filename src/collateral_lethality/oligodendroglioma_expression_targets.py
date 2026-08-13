# Expression-based (dosage-sensitive) collateral targets for heterozygous 1p/19q loss.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
from scipy.stats import pearsonr
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
E=pd.read_csv("uni_expr.csv",index_col=0)
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
p=pd.read_csv("paralog_1p19q.csv")
lines=M.index.intersection(E.index); Gv=M.loc[lines]; Ev=E.loc[lines]; cns=(lin.reindex(lines)=="CNS/Brain").values
gi={g:i for i,g in enumerate(Gv.columns)}; Gm=Gv.values.astype(float); cm=np.nanmean(Gm,0); ii=np.where(np.isnan(Gm)); Gm[ii]=np.take(cm,ii[1])
rows=[]
for A,sub in p.groupby("A"):
    if A not in Ev.columns: continue
    e=Ev[A].values; thr=np.nanpercentile(e,25); low=e<=thr
    if np.nansum(low)<15: continue
    for _,r in sub.iterrows():
        B=r.B
        if B not in gi: continue
        bi=gi[B]; bl=Gm[low,bi].mean(); bh=Gm[~low,bi].mean(); diff=bl-bh
        rr,_=pearsonr(np.nan_to_num(e),Gm[:,bi])
        cns_low=cns&low; b_cns=Gm[cns_low,bi].mean() if cns_low.sum()>=3 else np.nan
        rows.append((A,str(r.locA),B,round(bl,2),round(diff,3),round(rr,2),int(low.sum()),round(b_cns,2) if not np.isnan(b_cns) else np.nan))
d=pd.DataFrame(rows,columns=["partner","locus","target","B_dep_lowA","diff","corr","n_lowA","B_dep_CNSlowA"])
d=d[(d["diff"]<-0.15)&(d.B_dep_lowA<-0.5)].sort_values("diff").drop_duplicates("target")
d.to_csv("/home/jmilton/ml/ppset/oligo_1p19q_expression_targets.csv",index=False)
print("EXPRESSION-actionable 1p/19q collateral targets (paralog essential when partner expression is low):")
print(d.head(16).to_string(index=False))
