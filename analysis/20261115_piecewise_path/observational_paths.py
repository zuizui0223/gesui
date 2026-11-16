#!/usr/bin/env python3
"""EPA NRSA 2013–14 same-visit observational piecewise path decomposition.
This is statistical projection, NOT a causal SEM or WWTP-effluent effect.
Reuses original 1094 native visit UIDs reconstructed by existing source script.
"""
import argparse,collections,hashlib,json,os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

NUT=["logTN","logTP"]
GEO=["lat","lon","year","season_sin","season_cos","stream_order","logcond"]
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","logwidth"]
MED=["logchla","logafdm"]
CATS=["state","eco9"]
Y="ept_fraction"

def load(path):
    df=pd.read_csv(path,dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
    req=["UID","huc8","state","eco9",Y]+NUT+GEO+PHAB+MED
    if set(req)-set(df.columns):raise ValueError("Missing original source fields")
    if (len(df),df.huc8.nunique(),df.state.nunique())!=(1094,693,47):
        raise ValueError("Source native UID/HUC8/STATE cohort drift")
    if df.UID.duplicated().any() or not df[Y].between(0,1).all():raise ValueError("Invalid original EPT/UID")
    if not np.isfinite(df[NUT+GEO+PHAB+MED+[Y]].to_numpy(float)).all():
        raise ValueError("Non-finite source measurements")
    return df.sort_values("UID").reset_index(drop=True)

def matrix(df):
    z={}
    for name in NUT+GEO+PHAB+MED:
        a=df[name].to_numpy(float)
        if a.std(ddof=0)<1e-12:raise ValueError("Constant numeric "+name)
        z[name]=(a-a.mean())/a.std(ddof=0)
    names=["intercept"]+NUT+GEO+PHAB
    col=[np.ones(len(df))]+[z[x] for x in names[1:]]
    for cat in CATS:
        t=df[cat].astype(str).to_numpy()
        for level in sorted(set(t))[1:]:
            names.append(cat+"_"+level);col.append((t==level).astype(float))
    X=np.column_stack(col);M=np.column_stack([z[x] for x in MED])
    Z=np.column_stack([X,M]);y=df[Y].to_numpy(float)
    if np.linalg.matrix_rank(Z)!=Z.shape[1]:
        raise ValueError("Rank-deficient exact original model")
    h=df.huc8.astype(str).to_numpy()
    c=collections.Counter(h);w=np.array([1./c[k] for k in h])
    return X,Z,M,y,names,h,w

def solve(X,Z,M,y,names,w):
    sw=np.sqrt(w)
    b0=np.linalg.lstsq(sw[:,None]*X,sw[:,None]*np.column_stack([M,y]),rcond=None)[0]
    b1=np.linalg.lstsq(sw[:,None]*Z,sw*y,rcond=None)[0]
    out={}
    for nutrient in NUT:
        i=names.index(nutrient)
        a=b0[i,:2];b=b1[-2:]
        product=a*b
        total=float(b0[i,2]);direct=float(b1[i])
        if abs(total-direct-product.sum())>1e-7:
            raise ValueError("Statistical coefficient identity failed")
        out[nutrient]=dict(total=total,direct=direct,a_Chla=float(a[0]),a_AFDM=float(a[1]),
                           b_Chla=float(b[0]),b_AFDM=float(b[1]),product_Chla=float(product[0]),
                           product_AFDM=float(product[1]),product_sum=float(product.sum()))
    return out

def bootstrap(X,Z,M,y,names,h,B=500):
    unique,inv=np.unique(h,return_inverse=True)
    n=len(unique);cnt=np.bincount(inv)
    rng=np.random.default_rng(20261115)
    distributions={v:collections.defaultdict(list) for v in NUT}
    for iteration in range(B):
        sel=rng.integers(0,n,size=n)
        weights=np.bincount(sel,minlength=n)[inv]/cnt[inv]
        result=solve(X,Z,M,y,names,weights)
        for k,fields in result.items():
            for name,value in fields.items():distributions[k][name].append(value)
    return {nutrient:{name:list(map(float,np.quantile(values,[.025,.975])))
                      for name,values in fields.items()}
            for nutrient,fields in distributions.items()}

def cluster_intervals(X,Z,M,y,names,h,w):
    import statsmodels.api as sm
    fits={
      "Chla":sm.WLS(M[:,0],X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":h},use_t=False),
      "AFDM":sm.WLS(M[:,1],X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":h},use_t=False),
      "EPT_total":sm.WLS(y,X,weights=w).fit(cov_type="cluster",cov_kwds={"groups":h},use_t=False),
      "EPT_direct":sm.WLS(y,Z,weights=w).fit(cov_type="cluster",cov_kwds={"groups":h},use_t=False)}
    result={}
    for k in NUT:
        j=names.index(k);result[k]={}
        for stage,fit in fits.items():
            result[k][stage]={"beta":float(fit.params[j]),
                               "HUC8_cluster_normal95":list(map(float,fit.conf_int()[j]))}
    for stage,j in [("Chla",-2),("AFDM",-1)]:
        fit=fits["EPT_direct"]
        result["EPT_via_"+stage]={"beta":float(fit.params[j]),
                                 "HUC8_cluster_normal95":list(map(float,fit.conf_int()[j]))}
    return result

def features(tr,te,cols):
    x=[tr[cols].to_numpy(float)];v=[te[cols].to_numpy(float)]
    for cat in CATS:
        a=tr[cat].astype(str).to_numpy();b=te[cat].astype(str).to_numpy()
        for level in sorted(set(a))[1:]:
            x.append((a==level).astype(float)[:,None])
            v.append((b==level).astype(float)[:,None])
    return np.column_stack(x),np.column_stack(v)

def heldout(df):
    h=df.huc8.astype(str).to_numpy()
    fold=np.array([int(hashlib.sha256(k.encode()).hexdigest()[:8],16)%5 for k in h])
    stages={"geo_nutrients":GEO+NUT,
            "plus_habitat":GEO+NUT+PHAB,
            "plus_habitat_resources":GEO+NUT+PHAB+MED}
    yy=df[Y].to_numpy(float)
    pred={k:np.full(len(df),np.nan) for k in stages}
    for f in range(5):
        i=fold!=f;j=~i
        if i.sum()<500 or j.sum()<80:raise ValueError("HUC8 fold source support changed")
        for name,cols in stages.items():
            a,b=features(df.loc[i],df.loc[j],cols)
            ss=StandardScaler().fit(a)
            fit=Ridge(alpha=10,solver="svd").fit(ss.transform(a),yy[i])
            pred[name][j]=fit.predict(ss.transform(b))
    if any(not np.isfinite(p).all() for p in pred.values()):
        raise ValueError("No source-heldout prediction")
    groups,inv=np.unique(h,return_inverse=True)
    err={k:np.bincount(inv,weights=(yy-p)**2)/np.bincount(inv) for k,p in pred.items()}
    sel=np.random.default_rng(20261115).integers(0,len(groups),size=(2000,len(groups)))
    def gain(base,added):
        x=err[base];y=err[added]
        bootstrap=100*(x[sel].mean(axis=1)-y[sel].mean(axis=1))/x[sel].mean(axis=1)
        return {"HUC8_equal_MSE_gain_percent":float(100*(x.mean()-y.mean())/x.mean()),
                "fixed_error_HUC8_95":list(map(float,np.quantile(bootstrap,[.025,.975])))}
    return {"MSE_HUC8_equal":{k:float(x.mean()) for k,x in err.items()},
            "habitat_after_nutrient":gain("geo_nutrients","plus_habitat"),
            "resources_after_nutrient_and_habitat":gain("plus_habitat","plus_habitat_resources"),
            "all_folds_HUC8_disjoint":True}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path)
    ap.add_argument("--output",type=Path)
    ap.add_argument("--bootstrap",type=int,default=500)
    a=ap.parse_args()
    root=Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1314_independent_numeric"
    panel=a.panel or root/"AUDIT_ONLY_NRSA1314_AREAL_AFDM_SOURCE_PANEL.csv"
    out=a.output or root/"NRSA1314_OBSERVATIONAL_PIECEWISE_PATH.json"
    df=load(panel);X,Z,M,y,names,h,w=matrix(df)
    b=solve(X,Z,M,y,names,w)
    r={"schema":"GESUI_NRSA1314_OBSERVATIONAL_PATH_V1",
       "observed_UIDs":len(df),"HUC8":df.huc8.nunique(),"states":df.state.nunique(),
       "original_same_visit_nonrandomized":True,"ambient_TN_TP_not_sewage_attributed":True,
       "biofilm_standing_stock_not_primary_production":True,
       "EPT_fraction_not_emergence_or_absolute_biomass":True,
       "statistical_products_not_causal_mediation":True,
       "posthoc_exploratory":True,"multiple_testing_adjusted":False,
       "path_model":"TN/TP+geo+PHab -> Chla and AFDM; TN/TP+geo+PHab+Chla+AFDM -> EPT fraction",
       "same_original_data_exact_association_decomposition":b,
       "HUC8_cluster_bootstrap_product_95":bootstrap(X,Z,M,y,names,h,a.bootstrap),
       "HUC8_cluster_normal_intervals":cluster_intervals(X,Z,M,y,names,h,w),
       "HUC8_disjoint_heldout":heldout(df)}
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False))
    print(json.dumps({"paths":b,"out_of_HUC8_predictive":r["HUC8_disjoint_heldout"]},
                     ensure_ascii=False,indent=2,allow_nan=False))

if __name__=="__main__": main()
