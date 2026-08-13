import numpy as np, pandas as pd, pickle
print("loading...")
G=pd.read_csv("gene_effect_full.csv",index_col=0); G.columns=[c.split(" (")[0] for c in G.columns]
CN=pd.read_csv("uni_cn.csv",index_col=0); E=pd.read_csv("uni_expr.csv",index_col=0)
lines=G.index.intersection(CN.index).intersection(E.index)
G=G.loc[lines]; CN=CN.loc[lines]; E=E.loc[lines]; n=len(lines); print("common lines:",n)
gtex=pd.read_csv("gtex_feat.csv").set_index("gene")
bpp=pickle.load(open("bioplex_pairs.pkl","rb"))
pairs=pd.read_csv("paralog_pairs.csv")
# ---- loss = ploidy-normalized CN loss OR very low expression ----
cn=CN.values.astype(np.float32); ploidy=np.nanmedian(cn,axis=1,keepdims=True); cn_norm=cn/ (ploidy+1e-6)
ev=E.values.astype(np.float32); e_thr=np.nanpercentile(np.where(np.isnan(ev),np.inf,ev),5,axis=0)
CNcol={g:i for i,g in enumerate(CN.columns)}; Ecol={g:i for i,g in enumerate(E.columns)}
gcol={g:i for i,g in enumerate(G.columns)}
Gv=G.values.astype(np.float32); cm=np.nanmean(Gv,0); ii=np.where(np.isnan(Gv)); Gv[ii]=np.take(cm,ii[1])
ess_mean=Gv.mean(0); ess_min=Gv.min(0); ess_var=Gv.var(0)
Gz=(Gv-Gv.mean(0))/(Gv.std(0)+1e-6); Ez=np.nan_to_num((ev-np.nanmean(ev,0))/(np.nanstd(ev,0)+1e-6))
def lossmask(gene):
    m=np.zeros(n,bool)
    if gene in CNcol: m|= (cn_norm[:,CNcol[gene]]<0.75)
    if gene in Ecol:  m|= (ev[:,Ecol[gene]] < e_thr[Ecol[gene]])
    return m
pairs=pairs[pairs.A.isin(gcol)&pairs.B.isin(gcol)]
rows=[]
for A,sub in pairs.groupby("A"):
    if A not in gcol: continue
    low=lossmask(A); nlow=int(low.sum())
    if nlow<10 or nlow>n-10: continue
    ai=gcol[A]; Bs=sub.B.values; Bi=[gcol[b] for b in Bs]
    Gc=Gv[:,Bi]; bl=Gc[low].mean(0); bh=Gc[~low].mean(0); diff=bl-bh
    codep=(Gz[:,ai:ai+1]*Gz[:,Bi]).mean(0)
    aei=Ecol.get(A); coexpr=(Ez[:,aei:aei+1]*Ez[:,[Ecol.get(b,0) for b in Bs]]).mean(0) if aei is not None else np.zeros(len(Bs))
    for k,b in enumerate(Bs):
        gt=gtex.loc[b] if b in gtex.index else None
        rows.append((A,b,sub.maxid.values[k],sub.id_BA.values[k],sub.famA.values[k],
            ess_mean[Bi[k]],ess_min[Bi[k]],ess_var[Bi[k]],ess_mean[ai],nlow/n,
            float(codep[k]),float(coexpr[k]),
            float(gt.gtex_median) if gt is not None else np.nan,
            float(gt.gtex_breadth) if gt is not None else np.nan,
            int(frozenset((A,b)) in bpp),
            float(diff[k])))
d=pd.DataFrame(rows,columns=["A","B","maxid","minid","famA","B_ess_mean","B_ess_min","B_ess_var",
   "A_ess_mean","A_loss_freq","codep","coexpr","B_gtex_med","B_gtex_breadth","bioplex","diff"])
d.to_csv("paralog_ml_v2.csv",index=False)
print("pairs:",len(d))
for a,b in [("VPS4B","VPS4A"),("SMARCA2","SMARCA4"),("MAGOHB","MAGOH")]:
    r=d[(d.A==a)&(d.B==b)]
    if len(r): print(f"  {a}->{b}: diff={r['diff'].iloc[0]:+.3f} (top {(d['diff']<r['diff'].iloc[0]).mean()*100:.2f}%) bioplex={int(r.bioplex.iloc[0])} gtex_breadth={r.B_gtex_breadth.iloc[0]:.0f}")
