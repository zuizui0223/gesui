#!/usr/bin/env python3
"""Independent numeric FWL+augmented Ridge QA for original FR73 aquatic/land.

Original analysis: statsmodels cluster-WLS + normal-equation Ridge10
true-H6 fivefold. Independent audit: weighted Frisch–Waugh–Lovell residual
slope + augmented penalized least-squares (np.linalg.lstsq), checked against
4 frozen source-derived reported outputs. Additional 24-H6 leave-one-parent-
out sensitivity, all four tests no sign filtering. No independent field data,
no chemical sewage origin or ecological causal path.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_fr73_direct_TP_EIV_HNS"))
from same73_H7_water_land import source_join,BASE_EQR,BASE_PLANT,BASE_BUTTERFLY

SCHEMA="GESUI_20261116_FR73_FR_TP_NUTRI_HNS_INDEPENDENT_NUMERIC_QA_V1"
REFERENCE={
 "aquatic_EQR_on_measured_TP":(-.08455287908609084,24.68300304921941),
 "plant_NUTRI_on_measured_TP":(.020289920222520463,-2.978770289691385),
 "butterfly_HNS_on_measured_TP_after_plantN":(-.02305627778135466,-3.8987276606566685),
 "plant_NUTRI_on_butterfly_HNS_after_TP":(.14996914406426515,30.20856833781452)
}
CASES=[
 ("aquatic_EQR_on_measured_TP","actual_French_I2M2_EQR",BASE_EQR,"source_log_measured_TP"),
 ("plant_NUTRI_on_measured_TP","Plant_N_incidence",BASE_PLANT,"source_log_measured_TP"),
 ("butterfly_HNS_on_measured_TP_after_plantN","primary_incidence_mean",BASE_BUTTERFLY+["Plant_N_incidence"],"source_log_measured_TP"),
 ("plant_NUTRI_on_butterfly_HNS_after_TP","primary_incidence_mean",BASE_BUTTERFLY+["source_log_measured_TP"],"Plant_N_incidence")
]

def groupweights(groups):
    counts=collections.Counter(groups)
    return np.array([len(groups)/(len(counts)*counts[x]) for x in groups],float)

def mat(train,test,numeric):
    X=[np.ones(len(train),float)];Y=[np.ones(len(test),float)]
    for col in numeric:
        a=train[col].to_numpy(float);b=test[col].to_numpy(float)
        mean=a.mean();sd=a.std(ddof=0)
        if sd<1e-12:
            if len(train)==len(test):raise RuntimeError("Original WLS full source design has constant covariate")
            sd=1.
        X.append((a-mean)/sd);Y.append((b-mean)/sd)
    for col in ("country_method","period"):
        a=train[col].astype(str);b=test[col].astype(str)
        for k in sorted(a.unique().tolist())[1:]:
            X.append((a==k).to_numpy(float))
            Y.append((b==k).to_numpy(float))
    return np.column_stack(X),np.column_stack(Y)

def fwl(df,y,base,x):
    allfeatures=base+[x]
    design,_=mat(df,df,allfeatures)
    j=1+len(base)
    XX=np.delete(design,j,axis=1)
    xsrc=design[:,j]
    resp=df[y].to_numpy(float)
    w=np.sqrt(groupweights(df.h6.astype(str).to_numpy()))
    xr=xsrc-XX@np.linalg.lstsq(XX*w[:,None],xsrc*w,rcond=None)[0]
    yr=resp-XX@np.linalg.lstsq(XX*w[:,None],resp*w,rcond=None)[0]
    return float((w*w*xr*yr).sum()/(w*w*xr*xr).sum())

def hashfold(k):
    return int(hashlib.sha256(str(k).encode()).hexdigest()[:8],16)%5

def original_ridge_cv(df,y,features,mode):
    groups=df.h6.astype(str).to_numpy()
    folds=np.array([hashfold(g) for g in groups]) if mode=="5fold" else groups
    result=np.full(len(df),np.nan,float)
    for fold in sorted(set(folds)):
        tr=np.flatnonzero(folds!=fold);te=np.flatnonzero(folds==fold)
        if not len(tr) or not len(te) or set(groups[tr])&set(groups[te]):
            raise RuntimeError("Original source true H6 CV fold invalid")
        train=df.iloc[tr].reset_index(drop=True)
        test=df.iloc[te].reset_index(drop=True)
        X,T=mat(train,test,features)
        w=groupweights(train.h6.astype(str).to_numpy())
        penalty=np.diag([0.]+[np.sqrt(10.)]*(X.shape[1]-1))
        aa=np.vstack([np.sqrt(w)[:,None]*X,penalty])
        bb=np.concatenate([np.sqrt(w)*train[y].to_numpy(float),np.zeros(X.shape[1])])
        coef=np.linalg.lstsq(aa,bb,rcond=None)[0]
        result[te]=T@coef
    if not np.isfinite(result).all():raise RuntimeError("Missing original H6 out of fold prediction")
    return result

def gain(df,y,base,x,mode):
    ysrc=df[y].to_numpy(float)
    g=df.h6.astype(str).to_numpy()
    a=(ysrc-original_ridge_cv(df,y,base,mode))**2
    b=(ysrc-original_ridge_cv(df,y,base+[x],mode))**2
    unique=sorted(set(g))
    m0=np.array([a[g==key].mean() for key in unique])
    m1=np.array([b[g==key].mean() for key in unique])
    delta=float(100*(m0.mean()-m1.mean())/m0.mean())
    rng=np.random.default_rng(20261116)
    ix=rng.integers(0,len(unique),size=(2000,len(unique)))
    boot=100*(m0[ix].mean(axis=1)-m1[ix].mean(axis=1))/m0[ix].mean(axis=1)
    return {"gain_percent":delta,
        "H6_fixed_outofold_squared_error_resampling95":np.quantile(boot,[.025,.975]).tolist(),
        "n_true_H6":len(unique),"folds":len(set(np.array([hashfold(z) for z in g]) if mode=="5fold" else g))}

def main(a):
    state,df=source_join(a.fr,a.eiv)
    if df is None or len(df)!=146 or df.h.nunique()!=73 or df.h6.nunique()!=24:
        raise RuntimeError("Source-authentic 73 French H7 2018+22 EIV cohort failed")
    r=[]
    for name,y,base,x in CASES:
        beta=fwl(df,y,base,x)
        five=gain(df,y,base,x,"5fold")
        loho=gain(df,y,base,x,"leave_one_true_H6_out")
        b0,g0=REFERENCE[name]
        epsb=abs(beta-b0);epsg=abs(five["gain_percent"]-g0)
        if epsb>1e-9 or epsg>1e-7:
            raise RuntimeError("Original main source arithmetic independent verification FAIL "+name+
                  ": "+str((epsb,epsg,beta,five["gain_percent"])))
        r.append({"source_case":name,"original_beta":b0,"independent_FWL_beta":beta,
            "beta_abs_error":epsb,"original_H6_fivefold_gain":g0,
            "independent_augmented_lstsq_H6_fivefold_gain":five,
            "fivefold_gain_error_pp":epsg,
            "strict_24_true_H6_leave_one_parent_out":loho,
            "source_matched_French_H7periods":len(df),"true_H6":int(df.h6.nunique()),
            "numerical_identity_pass":True})
    a.output.mkdir(parents=True,exist_ok=True)
    out={"schema":SCHEMA,"status":"4_4_INDEPENDENT_NUMERIC_QA_PASS",
       "QA_method":"weighted Frisch-Waugh-Lovell + augmented least squares penalty",
       "original_method":"statsmodels cluster WLS + manual normal-equation Ridge10",
       "original_263_physical_field_stations_not_independent_GBIF_observations":True,
       "not_independent_field_replication":True,
       "not_source_specific_wwtp_or_causal":True,
       "cases":r}
    (a.output/"FR73_TP_EQR_NUTRI_HNS_INDEPENDENT_FWL_RIDGE_NUMERIC_QA.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_FR73_INDEPENDENT_SOURCE_NUMERIC_QA="+json.dumps({
        "status":out["status"],"cases":[{"name":q["source_case"],
        "original_beta":q["original_beta"],"QA_beta":q["independent_FWL_beta"],
        "original_5fold_gain":q["original_H6_fivefold_gain"],
        "QA_5fold_gain":q["independent_augmented_lstsq_H6_fivefold_gain"]["gain_percent"],
        "24H6_LOHO_gain":q["strict_24_true_H6_leave_one_parent_out"]["gain_percent"],
        "24H6_LOHO_interval":q["strict_24_true_H6_leave_one_parent_out"]["H6_fixed_outofold_squared_error_resampling95"]}
        for q in r]},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--fr",required=True,type=Path)
    p.add_argument("--eiv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
