import json,urllib.request,urllib.parse,time
genes=open("qpcr_genes.txt").read().split()
# mygene querymany POST
def post(url,data):
    req=urllib.request.Request(url,data=urllib.parse.urlencode(data).encode(),
        headers={"Content-Type":"application/x-www-form-urlencoded"})
    return json.load(urllib.request.urlopen(req,timeout=60))
DEEP={559292,4932,4896,3702,44689,6239,7227,7955,6669}  # yeast,fungi,plant,dicty,worm,fly,fish,daphnia
res=post("https://mygene.info/v3/query",{"q":",".join(genes),"scopes":"symbol",
         "fields":"pantherdb.ortholog.taxid,symbol","species":"human"})
rows=[]
for r in res:
    sym=r.get("symbol") or r.get("query")
    pj=r.get("pantherdb") or {}
    orth=pj.get("ortholog") or []
    if isinstance(orth,dict): orth=[orth]
    taxa={o.get("taxid") for o in orth if o.get("taxid")}
    rows.append((r.get("query"),sym,len(taxa),int(bool(taxa & {559292,4932,4896})),int(bool(taxa & DEEP))))
import pandas as pd
df=pd.DataFrame(rows,columns=["query","symbol","n_ortho_species","yeast_ortholog","deep_ortholog"]).drop_duplicates("query")
df.to_csv("conservation.csv",index=False)
print("genes with conservation:",len(df),"| have any ortholog:",(df.n_ortho_species>0).sum())
print("breadth quartiles:",df.n_ortho_species.quantile([0,.25,.5,.75,1]).tolist())
print("yeast-conserved:",int(df.yeast_ortholog.sum()),"| deep-conserved:",int(df.deep_ortholog.sum()))
print(df.sort_values("n_ortho_species",ascending=False).head(6).to_string(index=False))
print(df.sort_values("n_ortho_species").head(6).to_string(index=False))
