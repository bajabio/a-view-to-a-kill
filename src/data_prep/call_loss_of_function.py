# Calls loss of function per gene per cell line (mutation or deletion/low expression) and counts two-gene-loss backgrounds.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import numpy as np, pandas as pd
expr=pd.read_csv("bg_expr.csv").set_index("ModelID")
mut =pd.read_csv("bg_mut.csv").set_index("ModelID")
# cell lines present in CRISPR gene effect
ge_ids=pd.read_csv("gene_effect_full.csv",usecols=[0]).iloc[:,0].tolist()
genes=[c for c in expr.columns if c in mut.columns]
# LoF = damaging mutation OR expression below gene's 15th percentile
lof={}
for g in genes:
    e=expr[g]; low=e < e.quantile(0.15)
    m=(mut[g].reindex(expr.index).fillna(0)>0)
    lof[g]=(low | m)
L=pd.DataFrame(lof)
L=L.reindex(ge_ids).dropna(how="all")           # restrict to CRISPR lines
L=L.fillna(False)
print(f"cell lines (with CRISPR + omics): {len(L)}\n")
print("single-gene LoF frequency:")
print((L.sum().sort_values(ascending=False)).to_string())
print("\ndouble-LoF backgrounds (n cell lines with BOTH lost):")
pairs=[("MTAP","CDKN2A"),("STK11","KEAP1"),("TP53","CDKN2A"),("TP53","RB1"),
       ("CDKN2A","PTEN"),("SMAD4","TP53"),("ARID1A","TP53"),("NF1","TP53"),
       ("CDKN2A","CDKN2B"),("TP53","PTEN"),("VHL","PBRM1")]
for a,b in pairs:
    n=int((L[a]&L[b]).sum())
    print(f"  {a:7} + {b:7}: {n:4d}")
L.to_csv("lof_matrix.csv")
