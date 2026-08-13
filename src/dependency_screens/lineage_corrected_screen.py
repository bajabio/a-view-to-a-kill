# Lineage-corrected differential-dependency screen (Frisch-Waugh-Lovell) that recovers known synthetic-lethal partners.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats
from statsmodels.stats.multitest import multipletests
print("loading CRISPR..."); M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
mut=pd.read_csv("bg_mut.csv").set_index("ModelID"); expr=pd.read_csv("bg_expr.csv").set_index("ModelID")
genes=np.array(M.columns); G=M.values.astype(np.float32); cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1]); lines=M.index
gi={g:k for k,g in enumerate(genes)}

# ---------- Fig 2: MMP14 knockout effect vs PTEN expression ----------
d=pd.DataFrame({"MMP14":M["MMP14"]}).join(expr[["PTEN"]].rename(columns={"PTEN":"PTEN_expr"}),how="inner").dropna()
q20=d.PTEN_expr.quantile(.20)
fig,ax=plt.subplots(figsize=(9,5.4))
ax.axhspan(-0.5,d.MMP14.min()-0.3,color="#c0392b",alpha=.06); ax.axhline(-0.5,color="#c0392b",lw=1.3,ls="--")
ax.text(d.PTEN_expr.max(),-0.55,"below here = MMP14 is a dependency",color="#c0392b",fontsize=9,ha="right",va="top")
ax.axvspan(d.PTEN_expr.min()-0.3,q20,color="#0f8f74",alpha=.08)
ax.scatter(d.PTEN_expr,d.MMP14,s=16,c="#1f4e79",alpha=.55,edgecolor="none")
ax.text(q20,d.MMP14.max(),"PTEN-low\n(bottom 20%)",color="#0a5f4d",fontsize=9,ha="right",va="top")
ax.set_xlabel("PTEN expression  (log TPM)  ->  PTEN-deficient at left",fontsize=10)
ax.set_ylabel("MMP14 CRISPR gene effect\n(more negative = MMP14 knockout hurts)",fontsize=10)
ax.set_title("If MMP14-PTEN were synthetic-lethal, PTEN-low cells (left)\nshould drop below the dashed line. None do.",fontsize=12.5,fontweight="bold")
ax.set_ylim(d.MMP14.min()-0.25,d.MMP14.max()+0.15); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig("sl_stratified.png",dpi=150,bbox_inches="tight"); plt.close()

# ---------- Fig 3: lineage-corrected screens (mutation-defined loss) ----------
def design(v): dd=pd.get_dummies(v,dummy_na=True).astype(float); dd.insert(0,"_int",1.0); return dd.values
def resid(Y,Ld): b,_,_,_=np.linalg.lstsq(Ld,Y,rcond=None); return Y-Ld@b
Ld=design(lin.reindex(lines).values)
def screen(tgt):
    s=(mut[tgt].reindex(lines).fillna(0).values>0).astype(float)
    Xr=resid(G,Ld); yr=resid(s.reshape(-1,1),Ld).ravel()
    r=(Xr.T@yr)/np.sqrt((Xr**2).sum(0)*(yr**2).sum()+1e-12)
    df=len(lines)-np.linalg.matrix_rank(Ld)-1; t=r*np.sqrt(df/np.clip(1-r**2,1e-9,None)); p=2*stats.t.sf(np.abs(t),df)
    out=pd.DataFrame({"gene":genes,"r":r,"t":t,"p":p,"nlp":-np.log10(np.clip(p,1e-300,None))}); out["fdr"]=multipletests(out.p,method="fdr_bh")[1]
    out["n_def"]=int(s.sum()); return out
panels=[("PTEN","PIK3CB"),("SMARCA4","SMARCA2"),("ARID1A","ARID1B")]; extra={"ARID1A":["WRN","ARID2"]}
fig,axes=plt.subplots(1,3,figsize=(15,5.2))
for ax,(tgt,partner) in zip(axes,panels):
    dd=screen(tgt); sig=dd[(dd.r<0)&(dd.fdr<0.10)]
    ax.scatter(dd.r,dd.nlp,s=7,c="#c9d2cf",edgecolor="none"); ax.scatter(sig.r,sig.nlp,s=14,c="#0f8f74",edgecolor="none")
    pr=dd[dd.gene==partner].iloc[0]; ax.scatter([pr.r],[pr.nlp],s=70,c="#c0392b",zorder=6,edgecolor="white",linewidth=1)
    rank=int((dd.t<pr.t).sum())+1
    ax.annotate(f"{partner}\n(known SL, rank #{rank})",(pr.r,pr.nlp),color="#c0392b",fontsize=9,fontweight="bold",xytext=(8,-6),textcoords="offset points")
    for gx in extra.get(tgt,[]):
        rr=dd[dd.gene==gx].iloc[0]; ax.scatter([rr.r],[rr.nlp],s=34,c="#0a5f4d",zorder=5); ax.annotate(gx,(rr.r,rr.nlp),fontsize=8,color="#0a5f4d",xytext=(6,2),textcoords="offset points")
    ax.axvline(0,color="#999",lw=.9,ls="--"); ax.set_title(f"{tgt}-deficient cells\n({int(dd.n_def.iloc[0])} mutant lines)",fontsize=12,fontweight="bold")
    ax.set_xlabel("partial correlation\n<- more required when "+tgt+" lost",fontsize=9); ax.spines[["top","right"]].set_visible(False)
axes[0].set_ylabel("-log10 p",fontsize=10)
fig.suptitle("Lineage-corrected synthetic-lethal screens recover known partners at the top",fontsize=13.5,fontweight="bold",y=1.02)
fig.tight_layout(); fig.savefig("sl_v2_panels.png",dpi=150,bbox_inches="tight"); plt.close()
print("regen_b done: sl_stratified.png, sl_v2_panels.png")
