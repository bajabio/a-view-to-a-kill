import numpy as np, pandas as pd
print("loading...")
G=pd.read_csv("gene_effect_full.csv",index_col=0); G.columns=[c.split(" (")[0] for c in G.columns]
E=pd.read_csv("uni_expr.csv",index_col=0)
lines=G.index.intersection(E.index); G=G.loc[lines]; E=E.loc[lines]
pairs=pd.read_csv("paralog_pairs.csv")
gcol={g:i for i,g in enumerate(G.columns)}; ecol={g:i for i,g in enumerate(E.columns)}
Gv=G.values.astype(np.float32); cm=np.nanmean(Gv,0); ii=np.where(np.isnan(Gv)); Gv[ii]=np.take(cm,ii[1])
Ev=E.values.astype(np.float32)
n=len(lines)
# per-gene essentiality features + expression-low mask
ess_mean=np.nanmean(Gv,0); ess_min=np.nanmin(Gv,0); ess_var=np.nanvar(Gv,0)
# standardized for correlations
Gz=(Gv-Gv.mean(0))/(Gv.std(0)+1e-6)
Elow_thr=np.nanpercentile(np.where(np.isnan(Ev),np.inf,Ev),15,axis=0)
Ez=(np.nan_to_num(Ev-np.nanmean(Ev,0)))/(np.nanstd(Ev,0)+1e-6)
rows=[]
pairs=pairs[pairs.A.isin(gcol)&pairs.B.isin(gcol)&pairs.A.isin(ecol)&pairs.B.isin(ecol)]
for A,sub in pairs.groupby("A"):
    ai=gcol[A]; aei=ecol[A]
    low=Ev[:,aei] < Elow_thr[aei]
    nlow=int(np.nansum(low))
    if nlow<10 or nlow>n-10: continue
    Bs=sub.B.values; Bi=[gcol[b] for b in Bs]; Bei=[ecol[b] for b in Bs]
    Gc=Gv[:,Bi]; bl=Gc[low].mean(0); bh=Gc[~low].mean(0); diff=bl-bh
    codep=(Gz[:,ai:ai+1]*Gz[:,Bi]).mean(0)
    coexpr=(Ez[:,aei:aei+1]*Ez[:,Bei]).mean(0)
    lab=((diff<-0.15)&(bl<-0.5)).astype(int)
    for k,b in enumerate(Bs):
        rows.append((A,b,sub.maxid.values[k],sub.id_BA.values[k],sub.famA.values[k],
                     ess_mean[Bi[k]],ess_min[Bi[k]],ess_var[Bi[k]],ess_mean[ai],nlow/n,
                     float(codep[k]),float(coexpr[k]),float(diff[k]),lab[k]))
d=pd.DataFrame(rows,columns=["A","B","maxid","minid","famA","B_ess_mean","B_ess_min","B_ess_var",
      "A_ess_mean","A_loss_freq","codep","coexpr","diff","SL"])
d.to_csv("paralog_ml.csv",index=False)
print("labeled ordered pairs:",len(d),"| SL positives:",int(d.SL.sum()),f"({d.SL.mean()*100:.1f}%)")
for a,b in [("VPS4A","VPS4B"),("SMARCA4","SMARCA2"),("MAGOH","MAGOHB")]:
    r=d[(d.A==a)&(d.B==b)]
    if len(r): print(f"  {a}->{b}: SL_label={int(r.SL.iloc[0])} diff={r['diff'].iloc[0]:+.2f} maxid={r.maxid.iloc[0]:.0f}")
