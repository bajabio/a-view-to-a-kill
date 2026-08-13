import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import pearsonr
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
E=pd.read_csv("uni_expr.csv",index_col=0); lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
lines=M.index.intersection(E.index)
d=pd.DataFrame({"dep":M["ATP1A1"].reindex(lines),"part":E["ATP1A3"].reindex(lines),
                "cns":(lin.reindex(lines)=="CNS/Brain").values}).dropna()
r,p=pearsonr(d.part,d.dep)
fig,ax=plt.subplots(figsize=(8.6,5.3)); thr=d.part.quantile(0.25)
ax.axvspan(d.part.min()-0.3,thr,color="#0f8f74",alpha=.07); ax.axhline(-0.5,color="#c0392b",lw=1.2,ls="--")
ax.text(d.part.max(),-0.55,"ATP1A1 is a dependency",color="#c0392b",fontsize=9,ha="right",va="top")
oth=d[~d.cns]; cns=d[d.cns]
ax.scatter(oth.part,oth.dep,s=13,c="#c3ccc9",edgecolor="none",label="other cancers")
ax.scatter(cns.part,cns.dep,s=34,c="#0a5f4d",edgecolor="white",linewidth=.4,label="CNS / brain")
ax.text(thr,d.dep.max(),"ATP1A3 low\n(19q loss reduces dosage)",color="#0a5f4d",fontsize=9,ha="right",va="top")
ax.set_xlabel("ATP1A3 expression (log TPM)  ->  reduced by 19q loss at left",fontsize=10)
ax.set_ylabel("ATP1A1 knockout effect\n(more negative = essential)",fontsize=10)
ax.set_title(f"Dosage-sensitive collateral target: ATP1A1 becomes essential as its\n19q paralog ATP1A3 is reduced (r={r:+.2f}, p={p:.0e}) — no homozygous deletion needed",fontsize=11.5,fontweight="bold")
ax.legend(fontsize=9,frameon=False,loc="lower right"); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig("atp1a1_target.png",dpi=150,bbox_inches="tight"); print("saved; r=",round(r,2))
