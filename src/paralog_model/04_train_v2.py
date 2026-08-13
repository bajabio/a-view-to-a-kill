# Paralog model, step 4: train v2, compare against v1, and blind-test the validated pairs.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
rng=np.random.default_rng(0)
d=pd.read_csv("paralog_ml_v2.csv"); d["SL"]=(d["diff"]<-0.15).astype(int)
feats=["maxid","minid","famA","B_ess_mean","B_ess_min","B_ess_var","A_ess_mean","A_loss_freq",
       "codep","coexpr","B_gtex_med","B_gtex_breadth","bioplex"]
print(f"pairs={len(d)}  SL+={d.SL.sum()} ({d.SL.mean()*100:.2f}%)")
X=d[feats].values; y=d.SL.values
genes=sorted(set(d.A)|set(d.B)); gf={g:i%5 for i,g in enumerate(rng.permutation(genes))}
fa=d.A.map(gf).values; fb=d.B.map(gf).values
def cv(cols):
    Xf=d[cols].values; pr=np.full(len(y),np.nan)
    for k in range(5):
        te=(fa==k)&(fb==k); tr=(fa!=k)&(fb!=k)
        if te.sum()<5 or y[tr].sum()<5: continue
        c=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.06,max_depth=4,l2_regularization=1.0,class_weight="balanced",random_state=0).fit(Xf[tr],y[tr])
        pr[te]=c.predict_proba(Xf[te])[:,1]
    q=~np.isnan(pr); return roc_auc_score(y[q],pr[q]),average_precision_score(y[q],pr[q])
a,ap=cv(feats); print(f"\nV2 GENE-HELD-OUT: AUC={a:.3f}  AP={ap:.3f}   (V1 was AUC 0.984 / AP 0.163)")
for nm,cs in [("sequence identity only",["maxid","minid"]),
              ("essentiality only",["B_ess_mean","B_ess_min","B_ess_var","A_ess_mean"]),
              ("+GTEx+BioPlex (no ess)",["maxid","minid","famA","A_loss_freq","codep","coexpr","B_gtex_med","B_gtex_breadth","bioplex"]),
              ("full v2",feats)]:
    aa,_=cv(cs); print(f"  ablation [{nm:26}] AUC={aa:.3f}")
# blind test both directions
tg={"VPS4A","VPS4B","SMARCA2","SMARCA4","MAGOH","MAGOHB"}; hold=d.A.isin(tg)|d.B.isin(tg)
clf=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.05,max_depth=4,l2_regularization=1.0,class_weight="balanced",random_state=0).fit(X[~hold],y[~hold])
d["pred"]=clf.predict_proba(X)[:,1]
print("\nBLIND TEST (test genes fully held out) — best direction of each pair:")
for g1,g2 in [("VPS4A","VPS4B"),("SMARCA4","SMARCA2"),("MAGOH","MAGOHB")]:
    both=d[((d.A==g1)&(d.B==g2))|((d.A==g2)&(d.B==g1))]
    if len(both):
        r=both.loc[both.pred.idxmax()]; pct=(d.pred<r.pred).mean()*100
        print(f"  {r.A}->{r.B}: p(SL)={r.pred:.3f}  top {100-pct:.2f}%  (diff={r['diff']:+.2f})")
from sklearn.inspection import permutation_importance
pi=permutation_importance(clf,X[~hold],y[~hold],scoring="roc_auc",n_repeats=5,random_state=0)
print("\ntop features (AUC drop):")
for nm,v in sorted(zip(feats,pi.importances_mean),key=lambda z:-z[1])[:7]: print(f"  {nm:16} {v:.3f}")
d.sort_values("pred",ascending=False).to_csv("/home/jmilton/ml/ppset/paralog_SL_predictions_v2.csv",index=False)
