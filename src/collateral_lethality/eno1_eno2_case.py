# The ENO1/ENO2 collateral-lethality case: ENO2 is essential only where ENO1 is homozygously deleted.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=pd.read_csv("eno_cn.csv")
fig,ax=plt.subplots(figsize=(8.6,5.2))
ax.axvspan(0,0.3,color="#c0392b",alpha=.08); ax.axhline(-0.5,color="#c0392b",lw=1.2,ls="--")
ax.text(d.ENO1_norm.max(),-0.55,"ENO2 is a dependency",color="#c0392b",fontsize=9,ha="right",va="top")
deep=d[d.ENO1_norm<0.3]; rest=d[d.ENO1_norm>=0.3]
ax.scatter(rest.ENO1_norm,rest.ENO2_dep,s=14,c="#c3ccc9",edgecolor="none")
ax.scatter(deep.ENO1_norm,deep.ENO2_dep,s=90,c="#c0392b",edgecolor="white",linewidth=1,zorder=5,label=f"ENO1 homozygously deleted (n={len(deep)})")
ax.text(0.15,-0.9,"ENO1-deleted\n(1p36 loss)",color="#8f2318",fontsize=9,ha="center",fontweight="bold")
ax.set_xlabel("ENO1 copy number (relative to ploidy)  ->  1p36 deleted at left",fontsize=10)
ax.set_ylabel("ENO2 knockout effect\n(more negative = essential)",fontsize=10)
ax.set_title("Collateral lethality: ENO2 becomes essential only when\nits paralog ENO1 is homozygously deleted",fontsize=12.5,fontweight="bold")
ax.set_xlim(-0.05,d.ENO1_norm.quantile(0.99)); ax.legend(fontsize=9,frameon=False,loc="lower right"); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig("eno_case.png",dpi=150,bbox_inches="tight"); print("saved; deep-del ENO2:",round(deep.ENO2_dep.mean(),2))
