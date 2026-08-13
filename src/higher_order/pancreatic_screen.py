import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
print("loading..."); M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
L=pd.read_csv("lof_matrix.csv").set_index("ModelID").reindex(M.index).fillna(False).astype(bool)
genes=np.array(M.columns); G=M.values.astype(np.float32); cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1]); lines=M.index
panc=(lin.reindex(lines)=="Pancreas").values

print(f"\npancreatic cell lines in panel: {panc.sum()}")
for g in ["SMAD4","CDKN2A","TP53","KRAS"]:
    print(f"  {g} LoF: pancreatic {int(L[g].values[panc].sum())}/{panc.sum()}  | pan-cancer {int(L[g].sum())}/{len(L)}")

def design(v): d=pd.get_dummies(v,dummy_na=True).astype(float); d.insert(0,"_int",1.0); return d.values
def resid(Y,Ld): b,_,_,_=np.linalg.lstsq(Ld,Y,rcond=None); return Y-Ld@b
Ld=design(lin.reindex(lines).values)
def screen(flag):
    Xr=resid(G,Ld); yr=resid(flag.reshape(-1,1).astype(float),Ld).ravel()
    r=(Xr.T@yr)/np.sqrt((Xr**2).sum(0)*(yr**2).sum()+1e-12)
    df=len(lines)-np.linalg.matrix_rank(Ld)-1; t=r*np.sqrt(df/np.clip(1-r**2,1e-9,None)); p=2*stats.t.sf(np.abs(t),df)
    o=pd.DataFrame({"gene":genes,"t":t,"p":p}); o["fdr"]=multipletests(o.p,method="fdr_bh")[1]
    o["eff_flag"]=G[flag].mean(0); o["eff_rest"]=G[~flag].mean(0); return o

# ---- validation: SMAD4 single-loss recovers known partners? ----
s=screen(L["SMAD4"].values)
print("\n[VALIDATION] SMAD4-loss (pan-cancer): ranks of literature partners")
for g in ["RAB10","KLF5","AURKA","BRD2","BRD4","GPX4"]:
    if g in set(s.gene):
        r=s[s.gene==g].iloc[0]; rank=int((s.t<r.t).sum())+1
        print(f"   {g:6} t={r.t:+.2f} fdr={r.fdr:.2g} rank#{rank}  (eff loss={r.eff_flag:+.2f} vs {r.eff_rest:+.2f})")

# ---- PDAC two-gene LoF backgrounds ----
for a,b in [("SMAD4","CDKN2A"),("SMAD4","TP53"),("CDKN2A","TP53")]:
    flag=(L[a]&L[b]).values
    o=screen(flag); sl=o[(o.t<0)&(o.eff_flag<-0.3)].sort_values("t")
    npd=int((flag&panc).sum())
    print(f"\n===== {a}+{b} lost  (n={int(flag.sum())} pan-cancer, {npd} pancreatic) =====")
    print(sl.head(10)[["gene","t","fdr","eff_flag","eff_rest"]].to_string(index=False,float_format=lambda v:f"{v:.3f}"))
    o.to_csv(f"pdac_{a}_{b}.csv",index=False)
