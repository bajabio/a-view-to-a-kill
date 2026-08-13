# Streams the DepMap copy-number matrix and keeps only the genes of interest.
#
# Baja Bio · Bioinformatics — synthetic-lethality target-discovery toolkit
# Primary author: Jeff Milton
# Public dependency/genomic data per docs/DATA_SOURCES.md; set input paths at the top for your setup.

import sys,csv
csv.field_size_limit(10**7); want=set(open('universe.txt').read().split())
r=csv.reader(sys.stdin); h=next(r); sym=[c.split(' (')[0] for c in h]
idx={s:i for i,s in enumerate(sym) if s in want}; keep=[0]+[idx[s] for s in sorted(idx)]
cols=['ModelID']+[sym[i] for i in keep[1:]]
w=csv.writer(open('uni_cn.csv','w',newline='')); w.writerow(cols)
for row in r:
    if row: w.writerow([row[i] for i in keep])
print('CN genes:',len(cols)-1)
