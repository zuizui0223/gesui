#!/usr/bin/env python3
"""Does original upstream-vs-downstream reported WWTP N/P explain the EQR–plantN spatial relation?

Source-authentic same H7×period 451 rows, SAME native H6 folds and existing
observational W/geography/effort/country-method controls. Report status of 0
means *no matched positive reported outlet*, NOT effluent-free. Reported load
divided by model flow is potential concentration equivalent, NOT measured
receiving water chemistry, root-zone nutrient or causal treatment assignment.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_aquatic_plant_n_bridge"))
from same_H7_EQR_plantN import source_cohort,EQR,P_N,EQR_BG
sys.path.insert(0,str(ROOT/"analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import source_cluster_regression,cv_h6

DIR_PATH=ROOT/"results/20261115_np_direction/H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv"
BLOB_SHA="702131afc7d1a7bfb15e80357a506754cae7a70a"
UP="upstream_report_positive"
DOWN="downstream_report_positive"

def gitblob(path):
    d=path.read_bytes()
    return hashlib.sha1(f"blob {len(d)}\0".encode()+d).hexdigest()

def load_original(d):
    if gitblob(DIR_PATH)!=BLOB_SHA:
        raise RuntimeError("Original PR92/108 byte-identical derived source report index missing")
    rep=pd.read_csv(DIR_PATH,dtype={"h":str,"h6":str})
    if len(rep)!=454 or rep.duplicated(["h","period"]).any():
        raise RuntimeError("Original 454 annual source index grain changed")
    fields=["h","period","h6","W",EQR,UP,DOWN,
        "upstream_N_equiv_mgL","downstream_N_equiv_mgL",
        "upstream_P_equiv_mgL","downstream_P_equiv_mgL"]
    if not set(fields).issubset(rep):
        raise RuntimeError("Missing original directional WWTP source columns")
    check=d.merge(rep[fields],on=["h","period"],how="left",validate="one_to_one",
                  suffixes=("","_r"),indicator=True)
    if len(check)!=451 or (check._merge!="both").any():
        raise RuntimeError("Different true EIV trait and original WWTP reported source grain")
    if not (check.h6==check.h6_r).all():
        raise RuntimeError("Official parent H6 mismatch")
    if np.max(np.abs(check.W-check.W_r))>1e-10 or np.max(np.abs(check[EQR]-check[EQR+"_r"]))>1e-10:
        raise RuntimeError("Mismatch in original aquatic EQR or modeled W shared response source")
    check=check.drop(columns=["_merge","h6_r","W_r",EQR+"_r"])
    if check[UP].sum()!=56 or check[DOWN].sum()!=61 or check.h6.nunique()!=90:
        raise RuntimeError("Known source report opportunities unexpected")
    if set(check[UP].unique())!={0,1} or set(check[DOWN].unique())!={0,1}:
        raise RuntimeError("Report flags malformed")
    for n in ["N","P"]:
        for side in ["upstream","downstream"]:
            col=f"{side}_{n}_equiv_mgL"
            if not np.isfinite(check[col]).all() or (check[col]<0).any():
                raise RuntimeError("Original reported-source index invalid")
            check[f"log1p_{side}_{n}_reported_equiv"]=np.log1p(check[col].astype(float))
    return check

def fit(d,target,extra,base):
    val=source_cluster_regression(d,target,base+[extra],extra)
    cv=cv_h6(d,target,base,extra)
    return {"beta_per_1SD":val["adjusted_beta_per_1SD_focal"],
            "beta_H6_descriptive95":val["H6_cluster_descriptive95"],
            "H6_out_added_MSE_gain_percent":cv["H6_equal_added_MSE_gain_pct"],
            "H6_out_gain95_fixed_error":cv["fixed_out_of_fold_error_H6_resampling95_pct"],
            "source_predictor":extra,
            "n_H7period":len(d),"n_H6":d.h6.nunique(),
            "model_noncausal":True}

def main(a):
    d=load_original(source_cohort(a.source))
    a.output.mkdir(parents=True,exist_ok=True)
    cells={}
    for u in (0,1):
        for v in (0,1):
            cells[f"{u}_{v}"]=int(((d[UP]==u)&(d[DOWN]==v)).sum())
    if cells!={"0_0":369,"1_0":21,"0_1":26,"1_1":35}:
        raise RuntimeError("Original frozen source 2x2 report status drift")
    if d[UP].sum()!=56 or d.loc[d[UP]==1,"h6"].nunique()!=31:
        raise RuntimeError("Original EEA source report-positive true H6 identity drift")
    # Existing comparable tested original slope from PR110 must be reproduced;
    # no optimistic inference if source clone changed.
    orig=fit(d,EQR,P_N,EQR_BG)
    if abs(orig["beta_per_1SD"]+0.03445700952886715)>1e-9:
        raise RuntimeError("Original same-H7 EQR–plantN beta could not reproduce")
    source_presence=EQR_BG+[UP,DOWN]
    adjusted_presence=fit(d,EQR,P_N,source_presence)
    result={"schema":"GESUI_EU451_EQR_PLANT_NUTRI_UPDOWN_WWTP_REPORTED_SOURCE_ADVERSARIAL_V1",
        "source_original_454_report_index_gitblob":BLOB_SHA,
        "true_source_rows_eligible":len(d),
        "original_true_H6":int(d.h6.nunique()),
        "positive_report_status_2x2":cells,
        "upstream_report_positive_H7period":int(d[UP].sum()),
        "upstream_report_positive_true_H6":int(d.loc[d[UP]==1,"h6"].nunique()),
        "source_status_00_not_verified_zero_sewage":True,
        "full_451_original_EQR_plantN_no_directional_report_controls":orig,
        "full_451_original_EQR_plantN_after_UP_DOWN_status":adjusted_presence,
        "reported_status_increment_EQR_after_plantN":{
            "UP_after_DOWN":fit(d,EQR,UP,EQR_BG+[P_N,DOWN]),
            "DOWN_after_UP":fit(d,EQR,DOWN,EQR_BG+[P_N,UP])},
        "source_intensity_controls":{},
        "all_original_geographic_source_and_effort_controls_retained":True,
        "modeled_W_is_static":True,
        "annual_reported_NP_div_by_model_flow_not_measured_receiving_water":True,
        "not_ecological_causal_mediation_or_predator_energy_flux":True,
        "exploratory_posthoc_falsifier_only":True}
    for nutrient in ("N","P"):
        up=f"log1p_upstream_{nutrient}_reported_equiv"
        down=f"log1p_downstream_{nutrient}_reported_equiv"
        controls=EQR_BG+[UP,DOWN,up,down]
        result["source_intensity_controls"][nutrient]={
            "EQR_on_plantN_after_direction_and_reported_"+nutrient:
                fit(d,EQR,P_N,controls),
            "reported_UP_"+nutrient+"_after_report_status_DOWN_intensity_and_plantN":
                fit(d,EQR,up,EQR_BG+[P_N,UP,DOWN,down])}
    file=a.output/"GESUI_EQR_PLANTN_DIRECTIONAL_NP_SOURCE_CONTROL.json"
    file.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    s={"n":len(d),"H6":d.h6.nunique(),"cells":cells,
       "original_no_report":orig,
       "conditioned_both_report_flags":adjusted_presence,
       "report_status_EQR_after_plantN":result["reported_status_increment_EQR_after_plantN"],
       "reported_N_EQR_on_plantN":result["source_intensity_controls"]["N"]["EQR_on_plantN_after_direction_and_reported_N"],
       "reported_P_EQR_on_plantN":result["source_intensity_controls"]["P"]["EQR_on_plantN_after_direction_and_reported_P"],
       "source_UP_N_increment_EQR_after_plantN":result["source_intensity_controls"]["N"]["reported_UP_N_after_report_status_DOWN_intensity_and_plantN"]}
    print("GESUI_EQR_PLANTN_DIRECTIONAL_NP_CONTROL="+json.dumps(s,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
