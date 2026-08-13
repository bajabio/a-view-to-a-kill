# Per-cancer genetic profiles and passenger-deletion target opportunities (e.g. VPS4A, PRMT5).
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
L=pd.read_csv("lof_matrix.csv").set_index("ModelID").astype(bool)
vps=pd.read_csv("vps_expr.csv").set_index("ModelID")
lines=M.index; ln=lin.reindex(lines)
vps4b=vps["VPS4B"].reindex(lines); vps4b_low=(vps4b<vps4b.quantile(0.25))
def eff(g): return M[g] if g in M.columns else pd.Series(np.nan,index=lines)
for L_name in ["Breast","Prostate","Pancreas"]:
    m=(ln==L_name).values; n=m.sum()
    print(f"\n================  {L_name.upper()}  (n={n} cell lines)  ================")
    print("Genetic profile (fraction with loss-of-function):")
    for g in ["TP53","PTEN","CDKN2A","MTAP","SMAD4","RB1","ARID1A","NF1"]:
        if g in L.columns:
            fr=L[g].reindex(lines).fillna(False).values[m].mean(); print(f"   {g:7} {fr*100:4.0f}%")
    vlf=vps4b_low.reindex(lines).fillna(False).astype(bool).values[m].mean(); print(f"   VPS4B-low (18q) {vlf*100:4.0f}%")
    print("Passenger-deletion target opportunities:")
    mtap=L["MTAP"].reindex(lines).fillna(False).astype(bool).values
    for tgt in ["PRMT5","WDR77"]:
        e=eff(tgt).values; sub=m&mtap
        if sub.sum()>=3: print(f"   {tgt:6} in MTAP-lost {L_name} lines: {np.nanmean(e[sub]):+.2f} (n={sub.sum()}) vs MTAP-intact {np.nanmean(e[m&~mtap]):+.2f}")
    e=eff("VPS4A").values; sub=m&vps4b_low.reindex(lines).fillna(False).astype(bool).values
    if sub.sum()>=3: print(f"   VPS4A in VPS4B-low {L_name} lines: {np.nanmean(e[sub]):+.2f} (n={sub.sum()}) vs VPS4B-intact {np.nanmean(e[m&~vps4b_low.values]):+.2f}")
    # lineage-selective dependencies
    G=M.values.astype(np.float32); cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1])
    diff=G[m].mean(0)-G[~m].mean(0); linmean=G[m].mean(0)
    d=pd.DataFrame({"gene":M.columns,"diff":diff,"linmean":linmean})
    sel=d[(d["diff"]<-0.15)&(d.linmean<-0.5)].sort_values("diff").head(8)
    print(f"   Top {L_name}-selective dependencies:", ", ".join(f"{r.gene}({r['diff']:+.2f})" for _,r in sel.iterrows()))
