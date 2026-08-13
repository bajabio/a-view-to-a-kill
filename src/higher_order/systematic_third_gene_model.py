import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
L=pd.read_csv("lof_matrix.csv").set_index("ModelID").reindex(M.index).fillna(False).astype(bool)
genes=np.array(M.columns); G=M.values.astype(np.float32); cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1]); lines=M.index
D=pd.get_dummies(lin.reindex(lines),dummy_na=True).astype(float); D.insert(0,"_int",1.0); Dv=D.values
b,_,_,_=np.linalg.lstsq(Dv,G,rcond=None); Xr=G-Dv@b            # residualize dependency on lineage ONCE
dfree=len(lines)-np.linalg.matrix_rank(Dv)-1
ts=[g for g in ["ARID1A","BAP1","CDKN2A","KEAP1","MTAP","NF1","PBRM1","PTEN","RB1","SMAD4","SMARCA4","STK11","TP53","VHL"] if g in L.columns]
gi={g:i for i,g in enumerate(genes)}
rows=[]
for a_i in range(len(ts)):
    for b_i in range(a_i+1,len(ts)):
        A,B=ts[a_i],ts[b_i]; both=(L[A]&L[B]).values; n=int(both.sum())
        if n<20: continue
        yr=both.astype(float)-Dv@np.linalg.lstsq(Dv,both.astype(float),rcond=None)[0]
        r=(Xr.T@yr)/np.sqrt((Xr**2).sum(0)*(yr**2).sum()+1e-12); t=r*np.sqrt(dfree/np.clip(1-r**2,1e-9,None))
        p=2*stats.t.sf(np.abs(t),dfree); fdr=multipletests(p,method="fdr_bh")[1]
        Aonly=(L[A]&~L[B]).values; Bonly=(~L[A]&L[B]).values
        for j in np.argsort(t)[:6]:
            if t[j]>-3 or G[both,j].mean()>-0.4: continue
            ed=G[both,j].mean(); ea=G[Aonly,j].mean() if Aonly.sum()>=5 else np.nan; eb=G[Bonly,j].mean() if Bonly.sum()>=5 else np.nan
            worse=np.nanmin([ea,eb]); syn=ed-worse
            drv=A if (not np.isnan(ea) and ea<=np.nanmin([eb,0])) else B
            genuine = (syn<-0.1) and (t[j]<-3.5)
            rows.append((f"{A}+{B}",genes[j],round(float(t[j]),2),round(float(fdr[j]),4),round(ed,2),round(syn,3),
                         "genuine 3-way" if genuine else f"driven by {drv}",n))
d=pd.DataFrame(rows,columns=["background","target","t","fdr","eff_double","synergy","interpretation","n_double"])
d=d[d.fdr<0.1].sort_values("t").drop_duplicates(["target"]).reset_index(drop=True)
d.to_csv("/home/jmilton/ml/ppset/thirdgene_model_targets.csv",index=False)
print(f"THIRD-GENE MODEL: {len(d)} distinct targets (FDR<0.1) across {d.background.nunique()} two-loss backgrounds\n")
print(d.head(22).to_string(index=False))
print("\ngenuine 3-way count:", (d.interpretation=="genuine 3-way").sum(), "| passenger/single-driven:", (d.interpretation!="genuine 3-way").sum())
