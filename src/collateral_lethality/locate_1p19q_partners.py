import numpy as np, pandas as pd, json, urllib.request, urllib.parse
# candidate paralog pairs where target B is a selective dependency
M=pd.read_csv("gene_effect_full.csv",index_col=0); M.columns=[c.split(" (")[0] for c in M.columns]
E=pd.read_csv("uni_expr.csv",index_col=0)
p=pd.read_csv("paralogs.tsv",sep="\t",header=None,names=["A","B","idBA","idAB"]).dropna()
p["maxid"]=p[["idBA","idAB"]].max(axis=1); p=p[p.maxid>=20]
dep=set(M.columns); ex=set(E.columns)
essmin=M.min(0)
p=p[p.A.isin(ex)&p.B.isin(dep)]
p["Bmin"]=p.B.map(lambda b: essmin[b])
p=p[p.Bmin<-0.5]                         # target is a real dependency somewhere
partners=sorted(set(p.A))
print("candidate partner genes to locate:",len(partners))
# mygene map_location (cytoband)
def post(q):
    req=urllib.request.Request("https://mygene.info/v3/query",
        data=urllib.parse.urlencode({"q":",".join(q),"scopes":"symbol","fields":"map_location","species":"human"}).encode(),
        headers={"Content-Type":"application/x-www-form-urlencoded"})
    return json.load(urllib.request.urlopen(req,timeout=90))
loc={}
for i in range(0,len(partners),700):
    for r in post(partners[i:i+700]):
        ml=r.get("map_location"); 
        if ml: loc[r.get("query")]=ml
pd.Series(loc).to_csv("partner_loc.csv")
on1p=[g for g,m in loc.items() if str(m).startswith("1p")]
on19q=[g for g,m in loc.items() if str(m).startswith("19q")]
print(f"partners on 1p: {len(on1p)} | on 19q: {len(on19q)}")
p["locA"]=p.A.map(loc)
p1919=p[p.locA.astype(str).str.match(r'1p|19q')].copy()
p1919.to_csv("paralog_1p19q.csv",index=False)
print("paralog pairs with partner on 1p/19q (target is a dependency):",len(p1919))
print("example partners:",on1p[:8],"|",on19q[:8])
