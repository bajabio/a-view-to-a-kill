import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
M=pd.read_csv("gene_effect_full.csv",index_col=0)
lin=pd.read_csv("Model.csv").set_index("ModelID")["OncotreeLineage"]
L=pd.read_csv("lof_matrix.csv").set_index("ModelID").reindex(M.index).fillna(False).astype(bool)
ln=lin.reindex(M.index)
# background -> third-gene target (from the model)
apps=[("CDKN2A","MTAP","PRMT5 / WDR77"),("RB1","TP53","E2F3 / CDK2"),("SMAD4","TP53","VPS4A"),
      ("BAP1","PBRM1","PSMD6"),("ARID1A","SMAD4","KLF5"),("BAP1","KEAP1","ERBB2"),
      ("PTEN","RB1","YEATS2"),("SMARCA4","TP53","HAUS8")]
counts=ln.value_counts()
order=["Pancreas","Lung","Breast","Bowel","Ovary/Fallopian Tube","Esophagus/Stomach","Kidney",
       "Bladder/Urinary Tract","CNS/Brain","Skin","Head and Neck","Uterus","Prostate","Biliary Tract","Myeloid","Lymphoid"]
cols=[c for c in order if counts.get(c,0)>=12]
mat=np.zeros((len(apps),len(cols)))
for i,(a,b,tg) in enumerate(apps):
    dbl=(L[a]&L[b]).values
    for j,c in enumerate(cols): m=(ln==c).values; mat[i,j]=dbl[m].mean()*100
rows=[f"{a}+{b}  →  {tg}" for a,b,tg in apps]
fig,ax=plt.subplots(figsize=(13.5,5))
im=ax.imshow(mat,cmap="PuBuGn",aspect="auto",vmin=0,vmax=40)
ax.set_xticks(range(len(cols))); ax.set_xticklabels([f"{c}\n(n={int((ln==c).sum())})" for c in cols],fontsize=8,rotation=40,ha="right")
ax.set_yticks(range(len(apps))); ax.set_yticklabels(rows,fontsize=9.5,fontfamily="monospace")
for i in range(len(apps)):
    for j in range(len(cols)):
        ax.text(j,i,f"{mat[i,j]:.0f}",ha="center",va="center",fontsize=8,color="white" if mat[i,j]>24 else "#333")
cb=fig.colorbar(im,ax=ax,fraction=0.02,pad=0.01); cb.set_label("% of cancer with the\ntwo-loss background",fontsize=9)
ax.set_title("Third-gene targets across cancer types: percent carrying each two-loss background",fontsize=13,fontweight="bold")
fig.tight_layout(); fig.savefig("thirdgene_apps.png",dpi=150,bbox_inches="tight"); cp=1
print("saved. top per cancer:")
for j,c in enumerate(cols):
    i=mat[:,j].argmax(); print(f"  {c:22} {rows[i].split('→')[1].strip():14} ({mat[i,j]:.0f}%)")
