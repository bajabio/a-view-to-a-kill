# Paralog model, step 5: test generalizable features vs essentiality on 13 independently validated synthetic lethals.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
d=pd.read_csv("paralog_ml_v2.csv"); d["SL"]=(d["diff"]<-0.15).astype(int)
val=pd.read_csv("validated_SLs.txt")
vset=set(frozenset((a,b)) for a,b in zip(val.A1,val.A2)); vgenes=set(val.A1)|set(val.A2)
ess=["B_ess_mean","B_ess_min","B_ess_var","A_ess_mean"]
gen=["maxid","minid","famA","A_loss_freq","codep","coexpr","B_gtex_med","B_gtex_breadth","bioplex"]  # NO essentiality
hold=d.A.isin(vgenes)|d.B.isin(vgenes)     # remove all validated-gene pairs from training
def fit_pred(cols):
    c=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.05,max_depth=4,l2_regularization=1.0,class_weight="balanced",random_state=0).fit(d[cols][~hold],d.SL[~hold])
    return c.predict_proba(d[cols])[:,1]
d["pe"]=fit_pred(ess); d["pg"]=fit_pred(gen)
print("Recovery of 13 INDEPENDENT validated SL pairs (held out of training):")
print(f"{'pair':20} {'essentiality-only':>18} {'generalizable-only':>20}")
rowspe=[]; rowspg=[]
for a,b in [tuple(x) for x in val[['A1','A2']].values]:
    both=d[((d.A==a)&(d.B==b))|((d.A==b)&(d.B==a))]
    if not len(both): print(f"{a}/{b:12} not in paralog set"); continue
    pe=both.pe.max(); pg=both.pg.max()
    pce=(d.pe<pe).mean()*100; pcg=(d.pg<pg).mean()*100
    rowspe.append(pce); rowspg.append(pcg)
    print(f"{a+'/'+b:20} top {100-pce:5.2f}%          top {100-pcg:5.2f}%")
print(f"\nMEDIAN percentile of validated SLs:")
print(f"  essentiality-only : top {100-np.median(rowspe):.1f}%")
print(f"  generalizable-only: top {100-np.median(rowspg):.1f}%   (identity/complex/coexpr/GTEx, NO essentiality)")
