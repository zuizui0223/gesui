#!/usr/bin/env python3
"""Independent source-matched weighted FWL and augmented-lstsq H6 Ridge QA.

Rechecks 2 original-source primary responses x 3 report-source exposures x
two UP/DOWN directions. External data join uses immutable original source,
but coefficient and holdout numerical algorithms are independent from the
original statsmodels/normal-equation routines. Observational, not causal.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from conditional_models import load,OUTCOMES,TYPES

RIDGE_ALPHA=10.0

def weights(d):
    h=d.h6.astype(str)
    return len(d)/(h.nunique()*h.map(h.value_counts()).to_numpy(float))

def design(tr,te,columns):
    X=[np.ones(len(tr))];T=[np.ones(len(te))]
    for k in columns:
        x=tr[k].to_numpy(float);u=te[k].to_numpy(float)
        if not np.isfinite(x).all() or not np.isfinite(u).all():
            raise RuntimeError("Nonfinite original source field "+k)
        mu=float(x.mean());sd=float(x.std(ddof=0))
        if sd<1e-12:sd=1.0
        X.append((x-mu)/sd);T.append((u-mu)/sd)
    for cat in ("country_method","period"):
        train=tr[cat].astype(str)
        test=te[cat].astype(str)
        for name in sorted(train.unique())[1:]:
            X.append((train==name).to_numpy(float))
            T.append((test==name).to_numpy(float))
    return np.column_stack(X),np.column_stack(T)

def full_coeff(d,target,columns):
    X,_=design(d,d,columns)
    y=d[target].to_numpy(float)
    w=weights(d);root=np.sqrt(w)
    z=root[:,None]*X
    if np.linalg.matrix_rank(z)!=z.shape[1]:
        raise RuntimeError("Original complete H6 source design is not full rank")
    j=len(columns) # intercept, numeric columns, then categorical indicators
    Z0=np.delete(z,j,axis=1);z1=z[:,j]
    yw=root*y
    coeffx=np.linalg.lstsq(Z0,z1,rcond=None)[0]
    rx=z1-Z0@coeffx
    coeffy=np.linalg.lstsq(Z0,yw,rcond=None)[0]
    ry=yw-Z0@coeffy
    return float((rx@ry)/(rx@rx))

def ridge_predict(tr,te,target,columns):
    A,B=design(tr,te,columns)
    y=tr[target].to_numpy(float)
    root=np.sqrt(weights(tr))
    penalty=np.eye(A.shape[1])*np.sqrt(RIDGE_ALPHA)
    penalty[0,0]=0
    augX=np.vstack([A*root[:,None],penalty])
    augy=np.concatenate([y*root,np.zeros(A.shape[1])])
    coeff=np.linalg.lstsq(augX,augy,rcond=None)[0]
    return B@coeff

def cv_gain(d,target,base,added):
    hs=d.h6.astype(str).to_numpy()
    folds=np.array([int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5
                    for h in hs])
    pred0=np.full(len(d),np.nan);pred1=np.full(len(d),np.nan)
    for f in range(5):
        tr=d.iloc[np.flatnonzero(folds!=f)].reset_index(drop=True)
        te=d.iloc[np.flatnonzero(folds==f)].reset_index(drop=True)
        idx=np.flatnonzero(folds==f)
        if not len(idx) or set(tr.h6)&set(te.h6):
            raise RuntimeError("Geographic H6 blocked fold leakage")
        pred0[idx]=ridge_predict(tr,te,target,base)
        pred1[idx]=ridge_predict(tr,te,target,base+[added])
    if not np.isfinite(pred0).all() or not np.isfinite(pred1).all():
        raise RuntimeError("Unfilled H6 out-of-fold predictions")
    y=d[target].to_numpy(float)
    keys=sorted(set(hs))
    r0=np.array([np.mean((y[hs==h]-pred0[hs==h])**2) for h in keys])
    r1=np.array([np.mean((y[hs==h]-pred1[hs==h])**2) for h in keys])
    return float(100*(r0.mean()-r1.mean())/r0.mean())

def main(a):
    saved=json.loads(Path("results/20261116_directional_reported_np/UPDOWN_NP_PLANT_BUTTERFLY_VERIFIED_RESULTS.json").read_text())
    if saved.get("source_job_conclusion")!="success" or saved.get("source_run_id")!=43790119834:
        raise RuntimeError("Original source saved report was not successful exact published run")
    d=load(a)
    result=saved["result"]["outcomes"]
    primary=["plant_species_equal_NUTRI","butterfly_species_equal_HNS_adjusted_for_plantN"]
    n=0;max_b=0.;max_g=0.;checks=[]
    for outcome in primary:
        target,bg=OUTCOMES[outcome]
        for kind,(up,down) in TYPES.items():
            for side,other in (("up",down),("down",up)):
                base=bg if kind=="report_indicator" else bg+[
                   "upstream_report_positive","downstream_report_positive"]
                base=base+[other]
                added=up if side=="up" else down
                expected=result[outcome][kind]
                b=full_coeff(d,target,base+[added])
                gain=cv_gain(d,target,base,added)
                original_b=float(expected[side+"_beta"])
                original_gain=float(expected[side+"_gain"])
                delta_b=abs(b-original_b);delta_g=abs(gain-original_gain)
                max_b=max(max_b,delta_b);max_g=max(max_g,delta_g)
                if delta_b>1e-9 or delta_g>1e-7:
                    raise RuntimeError(f"Independent source reproduction mismatch: {outcome}/{kind}/{side}: beta {b} versus {original_b}; gain {gain} versus {original_gain}")
                checks.append({"target":outcome,"report_source":kind,"source_direction":side,
                               "verified_H6_weighted_beta":b,"verified_H6_oof_added_gain":gain,
                               "beta_abs_difference":delta_b,"gain_abs_difference":delta_g})
                n+=1
    if n!=12:raise RuntimeError("Did not verify all preselected 12 source contrasts")
    qa={"schema":"GESUI_20261116_ORIGINAL_REPORTED_NP_UPDOWN_INDEPENDENT_QA_V1",
        "n_original_H7period":len(d),"n_true_H6":int(d.h6.nunique()),
        "n_verified_comparisons":n,
        "max_absolute_beta_reproduction_difference":max_b,
        "max_absolute_H6_gain_percent_reproduction_difference":max_g,
        "source_numeric_algorithm":"independent weighted FWL plus augmented-lstsq Ridge10",
        "original_H6_unseen_parent_5fold_hash_preserved":True,
        "train_only_normalization":True,
        "not_2018_2022_receiver_water_TNTP":True,
        "not_causal_wwtp_effect":True,
        "all_tests_pass":True,
        "all_comparisons":checks}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"INDEPENDENT_SOURCE_12_CONTRASTS_QA.json").write_text(
        json.dumps(qa,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_INDEPENDENT_REPORTED_UPDOWN_NP_QA="+json.dumps(
        {k:qa[k] for k in ["n_original_H7period","n_true_H6",
                         "n_verified_comparisons",
                         "max_absolute_beta_reproduction_difference",
                         "max_absolute_H6_gain_percent_reproduction_difference",
                         "all_tests_pass"]},sort_keys=True),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--eiv",required=True,type=Path)
    ap.add_argument("--gate",required=True,type=Path)
    ap.add_argument("--output",required=True,type=Path)
    main(ap.parse_args())
