"""Run from gesui chat workspace with PR95 and the archived inputs (see gesui_paired.py).
FR archive from gesui Actions 43373639391. Fixed original hydrological folds.
"""
from pathlib import Path
import hashlib,sys
import pandas as pd,numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
out=Path(__file__).resolve().parent
repo=Path("work/gesui")
sys.path.insert(0,str((repo/"analysis/20261115_hns_family_taxon").resolve()))
from analyse import load_reconstruction,FAMS,COMMON
orig=Path("work/original")
d,t,B,s,known,fam,raw,_,_=load_reconstruction(orig/"H7_454_HNS_panel.csv",orig/"GESUI_303_HOST_NUTRIENT_SCORES.csv",Path("work/facets"))
mask=known&(fam=="Lycaenidae");n=B[:,mask];yes=n>0
ok=(n.sum(1)>=10)&(yes.sum(1)>=2)&(n.sum(1)>=.7*raw[:,FAMS.index("Lycaenidae")])
d["host_eq"]=np.where(ok,yes@s[mask]/np.maximum(yes.sum(1),1),np.nan)
d["host_rec"]=np.where(ok,n@s[mask]/np.maximum(n.sum(1),1),np.nan)
d["family_share"]=d.lycaenidae_fivefamily_index_share_logodds
path=next(Path("work/fr73").rglob("FR73*.csv"))
assert hashlib.sha256(path.read_bytes()).hexdigest()=="f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b"
fr=pd.read_csv(path,dtype={"HYBAS_ID_L07":str,"HYBAS_ID_L06":str})
d=d.merge(fr,left_on=["h","h6"],right_on=["HYBAS_ID_L07","HYBAS_ID_L06"],validate="many_to_one")
for pat in ["*ELIGIBILITY*.csv","*MEANS*.csv"]:
    q=pd.read_csv(next(Path("work/plant").glob(pat)),dtype={"h":str})
    d=d.merge(q,on=["h","period"],suffixes=("","_source"),validate="one_to_one")
assert len(d)==146 and d.full_source_trait_eligible.all()
late=d.period=="2020-2022";d["late"]=late.astype(float)
for label,field in [("TP","logTP"),("EQR","EQR")]:
    d[label]=np.where(late,d["mean_original_2022_"+field],d["mean_original_2018_"+field])
d["plant"]=d.Plant_N_incidence
d["scored"]=np.log1p(d.n_scored_species)
base=COMMON+["scored","NUTRI_coverage_all_Plantae_speciesKeys","log1p_original_same_site_2018_22_pair_n","late","W"]
# Recursive observational regressions. Common standardized scales fixed for block refits.
z=d.dropna(subset=["host_eq"]).reset_index(drop=True)
assert (len(z), z.h.nunique(), z.h6.nunique())==(139,72,24)
controls=[c for c in base if c!="W" and z[c].std()>1e-10]
equations={"TP":["W"],"EQR":["W","TP"],"plant":["W","TP"],"host_eq":["W","TP","EQR","plant"]}
columns=list(dict.fromkeys(controls+["W","TP","EQR","plant","host_eq"]))
x=(z[columns]-z[columns].mean())/z[columns].std(ddof=0)
groups=pd.factorize(z.h6,sort=True)[0]
G=len(set(groups)); size=np.bincount(groups)
models=[]
for y,terms in equations.items():
    X=np.column_stack([np.ones(len(z)),x[controls+terms]])
    models.append((y,terms,X,x[y].to_numpy()))
def fit(multiplicity):
    w=np.sqrt(multiplicity[groups]/size[groups])
    coefficients={}
    for y,terms,X,Y in models:
        b,_,rank,_=np.linalg.lstsq(X*w[:,None],Y*w,rcond=None)
        if rank!=X.shape[1]:return None
        for k,v in zip(terms,b[-len(terms):]):coefficients[k+">"+y]=v
    for path in ["W>TP>plant>host_eq","W>TP>EQR>host_eq","TP>plant>host_eq","W>plant>host_eq"]:
        nodes=path.split(">")
        coefficients[path]=np.prod([coefficients[a+">"+b] for a,b in zip(nodes,nodes[1:])])
    return coefficients
point=fit(np.ones(G));assert point is not None
rng=np.random.default_rng(1791);boot=[]
for i in range(3000):
    b=fit(np.bincount(rng.integers(G,size=G),minlength=G))
    if b is not None:boot.append(b)
assert len(boot)>2850
ci=pd.DataFrame(boot).quantile([.025,.975])
res=pd.DataFrame({"standardized":point,"low":ci.iloc[0],"high":ci.iloc[1],"valid_bootstraps":len(boot)})
res.to_csv(out/"fr_interface_paths.csv",index_label="path")
print(res.round(3))
# Test whether adding the two water indicators predicts an unseen catchment.
for response,background,extra in [("plant",base,["TP","EQR"]),("host_eq",base+["plant"],["TP","EQR"])]:
    err=[]
    for cols in [background,background+extra]:
        y=z[response].to_numpy();pred=np.full(len(z),np.nan)
        for f in sorted(z.original_H6_source_fold.unique()):
            test=z.original_H6_source_fold==f;train=~test
            scaler=StandardScaler().fit(z.loc[train,cols]);g=z.loc[train,"h6"]
            wt=1/g.map(g.value_counts()).to_numpy(float)
            model=Ridge(alpha=10).fit(scaler.transform(z.loc[train,cols]),y[train],sample_weight=wt/wt.mean())
            pred[test]=model.predict(scaler.transform(z.loc[test,cols]))
        err.append(pd.Series((pred-y)**2).groupby(z.h6).mean().to_numpy())
    a,b=err;draw=rng.integers(G,size=(4000,G));ci=np.quantile(100*(1-b[draw].mean(1)/a[draw].mean(1)),[.025,.975])
    print(response,"water joint heldout gain",100*(1-b.mean()/a.mean()),ci)
