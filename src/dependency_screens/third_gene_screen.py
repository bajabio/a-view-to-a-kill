# Screens a two-gene-loss background genome-wide for a third gene that becomes selectively essential.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
print("loading CRISPR..."); M=pd.read_csv("gene_effect_full.csv",index_col=0)
M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
L=pd.read_csv("lof_matrix.csv").set_index("ModelID")
genes=list(M.columns); G=M.values.astype(np.float32)
cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1]); lines=M.index

def design(v):
    d=pd.get_dummies(v,dummy_na=True).astype(float); d.insert(0,"_int",1.0); return d.values
def resid(Y,Ld): b,_,_,_=np.linalg.lstsq(Ld,Y,rcond=None); return Y-Ld@b

def screen(a,b,label):
    flag=(L[a]&L[b]).reindex(lines).fillna(False).astype(float).values
    Ld=design(lin.reindex(lines).values)
    Xr=resid(G,Ld); yr=resid(flag.reshape(-1,1),Ld).ravel()
    r=(Xr.T@yr)/np.sqrt((Xr**2).sum(0)*(yr**2).sum()+1e-12)
    df=len(lines)-np.linalg.matrix_rank(Ld)-1
    t=r*np.sqrt(df/np.clip(1-r**2,1e-9,None)); p=2*stats.t.sf(np.abs(t),df)
    res=pd.DataFrame({"gene":genes,"t":t,"r":r,"p":p})
    res["fdr"]=multipletests(res.p,method="fdr_bh")[1]
    res["eff_double"]=G[flag==1].mean(0); res["eff_rest"]=G[flag==0].mean(0)
    n=int(flag.sum())
    sl=res[(res.r<0)&(res.eff_double<-0.3)].sort_values("t")
    print(f"\n===== {label}  (n double-loss = {n}) =====")
    print("top genes MORE essential in the double-loss background:")
    print(sl.head(10)[["gene","t","fdr","eff_double","eff_rest"]].to_string(index=False,float_format=lambda v:f"{v:.3f}"))
    return res,flag

# validation
res,flag=screen("MTAP","CDKN2A","MTAP + CDKN2A (9p21 co-deletion)")
for g in ["PRMT5","MAT2A","RIOK1","MTAP"]:
    if g in set(res.gene):
        row=res[res.gene==g].iloc[0]; rank=int((res.t<row.t).sum())+1
        print(f"   validation {g}: t={row.t:+.2f} fdr={row.fdr:.1e} rank #{rank}")
# specificity for PRMT5: needs both, or driven by MTAP?
L2=L.reindex(M.index).fillna(False)
prmt5=M["PRMT5"].values
for lbl,mask in [("both MTAP&CDKN2A",(L2.MTAP&L2.CDKN2A)),("MTAP only",(L2.MTAP&~L2.CDKN2A)),
                 ("CDKN2A only",(~L2.MTAP&L2.CDKN2A)),("neither",(~L2.MTAP&~L2.CDKN2A))]:
    m=mask.values; print(f"   PRMT5 effect | {lbl:18}: {prmt5[m].mean():+.3f}  (n={m.sum()})")

screen("STK11","KEAP1","STK11 + KEAP1 (co-mutated lung)")
screen("TP53","RB1","TP53 + RB1")
