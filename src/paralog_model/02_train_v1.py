import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
rng=np.random.default_rng(0)
d=pd.read_csv("paralog_ml.csv")
d["SL"]=(d["diff"]<-0.15).astype(int)
feats=["maxid","minid","famA","B_ess_mean","B_ess_min","B_ess_var","A_ess_mean","A_loss_freq","codep","coexpr"]
print(f"pairs={len(d)}  SL positives={d.SL.sum()} ({d.SL.mean()*100:.2f}%)")
X=d[feats].values; y=d.SL.values

# gene-held-out CV (both genes of a test pair unseen)
genes=sorted(set(d.A)|set(d.B)); gf={g:i%5 for i,g in enumerate(rng.permutation(genes))}
fa=d.A.map(gf).values; fb=d.B.map(gf).values
proba=np.full(len(y),np.nan)
for k in range(5):
    te=(fa==k)&(fb==k); tr=(fa!=k)&(fb!=k)
    if te.sum()<5 or y[tr].sum()<5: continue
    clf=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.06,max_depth=4,l2_regularization=1.0,
        class_weight="balanced",random_state=0).fit(X[tr],y[tr])
    proba[te]=clf.predict_proba(X[te])[:,1]
m=~np.isnan(proba)&np.isin(True,True)
mm=~np.isnan(proba)
print(f"\nGENE-HELD-OUT (both genes unseen): AUC={roc_auc_score(y[mm],proba[mm]):.3f}  AP={average_precision_score(y[mm],proba[mm]):.3f}  (base {y[mm].mean()*100:.2f}%)")

# ablation
for name,fs in [("identity only",["maxid","minid"]),("essentiality only",["B_ess_mean","B_ess_min","B_ess_var","A_ess_mean"]),("full",feats)]:
    Xf=d[fs].values; pr=np.full(len(y),np.nan)
    for k in range(5):
        te=(fa==k)&(fb==k); tr=(fa!=k)&(fb!=k)
        if te.sum()<5 or y[tr].sum()<5: continue
        c=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.06,max_depth=4,l2_regularization=1.0,class_weight="balanced",random_state=0).fit(Xf[tr],y[tr])
        pr[te]=c.predict_proba(Xf[te])[:,1]
    q=~np.isnan(pr); print(f"  ablation [{name:18}] AUC={roc_auc_score(y[q],pr[q]):.3f}")

# blind test: hold out all pairs touching the test genes, predict them
testgenes={"VPS4A","VPS4B","SMARCA2","SMARCA4","MAGOH","MAGOHB"}
hold=d.A.isin(testgenes)|d.B.isin(testgenes)
clf=HistGradientBoostingClassifier(max_iter=400,learning_rate=0.05,max_depth=4,l2_regularization=1.0,class_weight="balanced",random_state=0).fit(X[~hold],y[~hold])
d["pred"]=clf.predict_proba(X)[:,1]
print("\nBLIND TEST (all test-gene pairs excluded from training), predicted SL probability + genome-wide percentile:")
for a,b in [("VPS4B","VPS4A"),("SMARCA2","SMARCA4"),("MAGOHB","MAGOH")]:
    r=d[(d.A==a)&(d.B==b)]
    if len(r):
        p=r.pred.iloc[0]; pct=(d.pred<p).mean()*100
        print(f"  {a:8}-> {b:8}  p(SL)={p:.3f}  (higher than {pct:.1f}% of all paralog pairs)  true diff={r['diff'].iloc[0]:+.2f}")
# feature importance
from sklearn.inspection import permutation_importance
pi=permutation_importance(clf,X[~hold],y[~hold],scoring="roc_auc",n_repeats=5,random_state=0)
print("\nfeature importance (AUC drop):")
for nm,v in sorted(zip(feats,pi.importances_mean),key=lambda z:-z[1])[:6]: print(f"  {nm:14} {v:.3f}")
