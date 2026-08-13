import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
E=pd.read_csv("uni_expr.csv",index_col=0); lines=M.index; ln=lin.reindex(lines)
El=E.reindex(lines)
# paralog targets: (target, partner-whose-loss-creates-vulnerability)
targets=[("VPS4A","VPS4B"),("PRMT5","MTAP"),("COPG1","COPG2"),("MAGOH","MAGOHB"),
         ("NXT1","NXT2"),("INTS6","INTS6L"),("PPP2CA","PPP2CB")]
partners=[p for _,p in targets]
thr={p:El[p].quantile(0.25) for p in partners if p in El.columns}
counts=ln.value_counts(); keep=[c for c in counts.index if counts[c]>=12 and isinstance(c,str)]
order=["Pancreas","Breast","Lung","Bowel","Ovary/Fallopian Tube","Esophagus/Stomach","Bladder/Urinary Tract",
       "Kidney","CNS/Brain","Skin","Head and Neck","Uterus","Prostate","Lymphoid","Myeloid","Bone","Soft Tissue","Biliary Tract"]
cols=[c for c in order if c in keep]
mat=np.full((len(targets),len(cols)),np.nan)
for i,(t,p) in enumerate(targets):
    if p not in thr: continue
    for j,c in enumerate(cols):
        m=(ln==c).values; mat[i,j]=np.nanmean(El[p].values[m]<thr[p])*100
labels=[f"{t}  (via {p} loss)" for t,p in targets]
fig,ax=plt.subplots(figsize=(14,5.6))
im=ax.imshow(mat,cmap="YlGnBu",aspect="auto",vmin=0,vmax=55)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([f"{c}\n(n={int((ln==c).sum())})" for c in cols],fontsize=8.5,rotation=40,ha="right")
ax.set_yticks(range(len(targets))); ax.set_yticklabels(labels,fontsize=10,fontfamily="monospace")
for i in range(len(targets)):
    for j in range(len(cols)):
        if not np.isnan(mat[i,j]): ax.text(j,i,f"{mat[i,j]:.0f}",ha="center",va="center",fontsize=8,color="white" if mat[i,j]>32 else "#333")
cb=fig.colorbar(im,ax=ax,fraction=0.02,pad=0.01); cb.set_label("% of cancer addressable\n(partner lost)",fontsize=9)
ax.set_title("Passenger-deletion paralog targets across cancer types: percent of each cancer addressable",fontsize=13,fontweight="bold")
fig.tight_layout(); fig.savefig("cross_cancer_targets.png",dpi=150,bbox_inches="tight")
print("saved. Pancreas column:")
for (t,p),v in zip(targets,mat[:,cols.index("Pancreas")]): print(f"  {t:8} {v:.0f}%")
