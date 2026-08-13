# Decomposes a two-gene-loss dependency into its single-loss parts to tell a real interaction from a co-deleted passenger.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
L=pd.read_csv("lof_matrix.csv").set_index("ModelID").reindex(M.index).fillna(False).astype(bool)
G=M.values.astype(np.float32); cm=np.nanmean(G,0); ii=np.where(np.isnan(G)); G[ii]=np.take(cm,ii[1])
genes=np.array(M.columns)
def groups(a,b):
    A,B=L[a].values.astype(bool),L[b].values.astype(bool)
    return {"both":A&B,"Aonly":A&~B,"Bonly":~A&B,"neither":~A&~B}
def synergy(a,b,label,min_single=10):
    g=groups(a,b)
    if g["both"].sum()<20 or g["Aonly"].sum()<min_single or g["Bonly"].sum()<min_single:
        print(f"{label}: insufficient single-loss lines"); return
    ed=G[g["both"]].mean(0); ea=G[g["Aonly"]].mean(0); eb=G[g["Bonly"]].mean(0); en=G[g["neither"]].mean(0)
    worse_single=np.minimum(ea,eb)                       # the more-essential single background
    syn=ed-worse_single                                  # <0 => double worse than either single (true 3-way)
    real=(ed<-0.4)                                       # must be a genuine dependency in double
    d=pd.DataFrame({"gene":genes,"eff_both":ed,"eff_Aonly":ea,"eff_Bonly":eb,"eff_neither":en,"synergy":syn})
    d=d[real].sort_values("synergy")
    print(f"\n===== {label}: genes requiring loss of BOTH (synergy = both - worse single) =====")
    print(f"  n: both={g['both'].sum()} {a}only={g['Aonly'].sum()} {b}only={g['Bonly'].sum()}")
    print(d.head(8).to_string(index=False,float_format=lambda v:f"{v:.3f}"))
synergy("TP53","RB1","TP53 + RB1")
synergy("MTAP","CDKN2A","MTAP + CDKN2A")
synergy("STK11","KEAP1","STK11 + KEAP1")
synergy("TP53","PTEN","TP53 + PTEN")
synergy("ARID1A","TP53","ARID1A + TP53")
