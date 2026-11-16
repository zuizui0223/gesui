#!/usr/bin/env python3
"""Independent matrix-level numeric QA of same-H7 water–land bridge.

No reuse of original statsmodels, H6_model_if_gate.design_full, H6
Ridge normal-equation predictor or cv_h6 in this independent verifier.
Verify primary conditional β via weighted FWL residualization; verify H6
group-disjoint Ridge10 5-fold predictions via augmented least squares.
Compare frozen authentic original EIV, GBIF, EQR, WWTP-source panels.
Not a new ecological replication or a causal test.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_aquatic_plant_n_bridge"))
from same_H7_EQR_plantN import source_cohort,EQR,P_N,EQR_BG
from source_report_direction_controls import load_original,UP,DOWN
from paired_h7_EQR_plantN_temporal import pair_source,CORE

ALPHA=10.
CATS=("country_method","period")

def matrix(train,test,numeric):
    X=[np.ones(len(train),dtype=float)]
    T=[np.ones(len(test),dtype=float)]
    for var in numeric:
        x=train[var].to_numpy(dtype=float)
        t=test[var].to_numpy(dtype=float)
        mu=x.mean();sd=x.std(ddof=0)
        if sd<1e-12:sd=1.
        X.append((x-mu)/sd)
        T.append((t-mu)/sd)
    for col in CATS:
        g=train[col].astype(str)
        v=test[col].astype(str)
        for level in sorted(g.unique().tolist())[1:]:
            X.append((g==level).to_numpy(float))
            T.append((v==level).to_numpy(float))
    return np.column_stack(X),np.column_stack(T)

def groupweights(d):
    h=d.h6.astype(str)
    counts=h.value_counts().to_dict()
    return np.asarray([len(d)/(h.nunique()*counts[z]) for z in h],float)

def fwl_beta(d,response,base,x):
    mat,_=matrix(d,d,base+[x])
    w=groupweights(d)
    xj=mat[:,-1] if len(CATS)==0 else mat[:,len(base)+1]
    X0=np.delete(mat,len(base)+1,axis=1)
    y=d[response].to_numpy(float)
    r=np.sqrt(w)
    xc=xj-X0@np.linalg.lstsq(X0*r[:,None],xj*r,rcond=None)[0]
    yc=y-X0@np.linalg.lstsq(X0*r[:,None],y*r,rcond=None)[0]
    if (xc*w*xc).sum()<1e-10:raise RuntimeError("No residual focal variation")
    return float((xc*w*yc).sum()/(xc*w*xc).sum())

def fold(group):
    return int(hashlib.sha256(str(group).encode()).hexdigest()[:8],16)%5

def pred_aug(train,test,response,numeric):
    X,T=matrix(train,test,numeric)
    w=groupweights(train)
    pen=np.diag([0.]+[np.sqrt(ALPHA)]*(X.shape[1]-1))
    A=np.vstack([X*np.sqrt(w)[:,None],pen])
    b=np.concatenate([train[response].to_numpy(float)*np.sqrt(w),np.zeros(X.shape[1])])
    coef=np.linalg.lstsq(A,b,rcond=None)[0]
    return T@coef

def gain(d,response,base,x):
    h=d.h6.astype(str).to_numpy()
    y=d[response].to_numpy(float)
    folds=np.asarray([fold(z) for z in h],int)
    p0=np.full(len(d),np.nan);p1=np.full(len(d),np.nan)
    for f in range(5):
        i=np.flatnonzero(folds!=f);j=np.flatnonzero(folds==f)
        if not len(i) or not len(j) or len(i)<30 or set(h[i])&set(h[j]):
            raise RuntimeError("Actual H6 holdout source invalid")
        a=d.iloc[i].reset_index(drop=True);b=d.iloc[j].reset_index(drop=True)
        p0[j]=pred_aug(a,b,response,base)
        p1[j]=pred_aug(a,b,response,base+[x])
    if not np.isfinite(p0).all() or not np.isfinite(p1).all():
        raise RuntimeError("Original heldout prediction null")
    aa=[];bb=[]
    for group in sorted(set(h)):
        sel=h==group
        aa.append(((y[sel]-p0[sel])**2).mean())
        bb.append(((y[sel]-p1[sel])**2).mean())
    return float(100*(np.mean(aa)-np.mean(bb))/np.mean(aa))

def main(a):
    d=load_original(source_cohort(a.source))
    original,pairs,ok=pair_source(a.source)
    if not ok or len(pairs)!=170 or pairs.h6.nunique()!=53:
        raise RuntimeError("Original same-method paired 170 H7 / 53 H6 changed")
    source_cases=[
        ("EQR_plantN_no_report",d,EQR,EQR_BG,P_N,-.03445700952886715,2.837935847960906),
        ("EQR_plantN_with_UPDOWN_flags",d,EQR,EQR_BG+[UP,DOWN],P_N,-.033912458074206134,1.3790440807410311),
        ("EQR_plantN_with_UPDOWN_reported_N",d,EQR,EQR_BG+[UP,DOWN,"log1p_upstream_N_reported_equiv","log1p_downstream_N_reported_equiv"],P_N,-.03334852849171598,.9158337505521202),
        ("paired_EQR_delta_plantN",pairs,"d_EQR",CORE,"d_plant_N",.012090635823906692,-2.029868787231651)]
    rows=[]
    for name,dset,y,base,x,ref_b,ref_g in source_cases:
        b=fwl_beta(dset,y,base,x)
        g=gain(dset,y,base,x)
        berr=abs(b-ref_b);gerr=abs(g-ref_g)
        rows.append({"name":name,"rows":len(dset),"H6":int(dset.h6.nunique()),
                     "independent_FWL_beta":b,"published_original_beta":ref_b,
                     "absolute_beta_error":berr,
                     "independent_augmented_lstsq_H6_out_gain":g,
                     "published_original_H6_out_gain":ref_g,
                     "absolute_gain_error_pp":gerr,
                     "numerical_QA_pass":bool(berr<1e-9 and gerr<1e-7)})
        if berr>1e-9 or gerr>1e-7:
            raise RuntimeError("Independent original aquatic–terrestrial numeric replay mismatch "+name)
    report={"schema":"GESUI_20261116_INDEPENDENT_EQR_PLANTN_WATER_LAND_NUMERIC_QA_V1",
      "status":"PASS_ALL_4" if len(rows)==4 else "STOP",
      "all_original_source_sha_EIV_and_direction_verified_in_input":True,
      "not_independent_ecological_replication":True,
      "not_wastewater_causation":True,
      "main_original_solver":"statsmodels WLS and weighted normal-equation Ridge10",
      "QA_solver":"weighted Frisch-Waugh-Lovell and augmented least squares Ridge10",
      "results":rows}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"GESUI_WATER_LAND_NUMERIC_QA.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_WATER_LAND_NUMERIC_QA="+json.dumps(report,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
