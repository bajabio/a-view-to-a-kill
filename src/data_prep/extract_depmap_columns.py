# Streams a large DepMap matrix and extracts the requested gene columns.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import sys,csv
csv.field_size_limit(10**7)
want=set(open("bg_genes.txt").read().split()); out=sys.argv[1]
r=csv.reader(sys.stdin); h=next(r); sym=[c.split(" (")[0] for c in h]
idx={s:i for i,s in enumerate(sym) if s in want}; keep=[0]+[idx[s] for s in sorted(idx)]
cols=["ModelID"]+[sym[i] for i in keep[1:]]
with open(out,"w",newline="") as o:
    w=csv.writer(o); w.writerow(cols); n=0
    for row in r:
        if row: w.writerow([row[i] for i in keep]); n+=1
print(out,"->",cols[1:],"| rows",n)
