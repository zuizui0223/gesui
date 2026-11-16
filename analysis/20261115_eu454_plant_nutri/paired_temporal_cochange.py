#!/usr/bin/env python3
"""Original within-H7 two-period DELTA plant EIV NUTRI vs butterfly HNS.

Source-corrected 451 H7-period plant scores (2000s source tags unchanged);
paired H7 only. W is STATIC, hence cannot explain temporal change. This
is a posthoc noncausal cochange stress test, not treatment-response or
environmentally mediated food-web evidence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
from h6_model_if_gate import (ORIGINAL_PANEL_SHA,FACETS_SHA,EIV_SHA,
                              source_cluster_regression,cv_h6)

ROOT=Path(__file__).resolve().parents[2]
DIFF_CONTROLS=["d_nPlant","d_nButterfly","d_plant_coverage",
               "d_primary_record_coverage_raw","d_primary_host_completeness"]
START="2017-2019"
END="2020-2022"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ci_h6_mean(r,c,seed=20261115):
    # One equally weighted value per genuine H6 source, not a pseudoindependent
    # second-period H7. H6 resampling fixed values is descriptive only.
    dd=r.groupby("h6")[[c[0],c[1]]].mean().reset_index()
    a=dd[c[0]].to_numpy(float);b=dd[c[1]].to_numpy(float)
    if a.std()<1e-12 or b.std()<1e-12:
        return {"status":"UNINFORMATIVE_CONSTANT_H6"}
    rho=float(np.corrcoef(a,b)[0,1])
    rng=np.random.default_rng(seed)
    arr=[]
    for _ in range(2000):
        take=rng.integers(0,len(a),len(a))
        if np.std(a[take])<1e-10 or np.std(b[take])<1e-10:continue
        arr.append(float(np.corrcoef(a[take],b[take])[0,1]))
    return {"H6_groupmean_Pearson_r":rho,
            "H6_fixed_difference_resampling95":np.quantile(arr,[.025,.975]).tolist() if len(arr)>100 else None,
            "original_true_H6":int(dd.h6.nunique()),
            "not_a_treatment_effect_or_independent_replication":True}

def main(args):
    gate=json.loads((args.source/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not gate.get("hard_source_support_pass")
        or gate.get("source_facets_SHA256")!=FACETS_SHA
        or gate.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA):
        raise RuntimeError("True original v2 hybrid source gate required")
    panel=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if sha(panel)!=ORIGINAL_PANEL_SHA:raise RuntimeError("True original 454-panel SHA mismatch")
    full=pd.read_csv(panel,dtype={"h":str,"h6":str})
    source=pd.read_csv(args.source/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",dtype={"h":str})
    elig=pd.read_csv(args.source/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",dtype={"h":str})
    if len(full)!=454 or len(source)!=451 or len(elig)!=454:
        raise RuntimeError("Source panel/source trait quality cohort mismatch")
    sub=elig.loc[elig.full_source_trait_eligible.astype(str)=="True",
                 ["h","period","NUTRI_coverage_all_Plantae_speciesKeys"]]
    panel=full.merge(source,on=["h","period"],validate="one_to_one",how="inner")
    panel=panel.merge(sub,on=["h","period"],validate="one_to_one",how="inner")
    if len(panel)!=451:raise RuntimeError("Source 451 H7period scores not exact")
    if set(panel.period.unique())!={START,END}:raise RuntimeError("Source survey periods unexpectedly changed")
    rows=[]
    for h,g in panel.groupby("h",sort=True):
        if len(g)!=2:continue
        if set(g.period)!={START,END}:raise RuntimeError("Duplicate period within original H7")
        x=g.loc[g.period==START].iloc[0];y=g.loc[g.period==END].iloc[0]
        if x.h6!=y.h6 or x.country_holdout!=y.country_holdout:
            raise RuntimeError("H7 changes H6 parent or original holdout-country")
        row={"h":str(h),"h6":str(x.h6),"country_method":str(x.country_holdout),
             "period":"delta_2020_2022_minus_2017_2019"}
        keys={
          "d_plant_N_incidence":"Plant_N_incidence",
          "d_plant_N_record":"Plant_N_record",
          "d_butterfly_HNS_incidence":"primary_incidence_mean",
          "d_butterfly_HNS_record":"primary_record_mean",
          "d_nPlant":"nPlant",
          "d_nButterfly":"nButterfly",
          "d_plant_coverage":"NUTRI_coverage_all_Plantae_speciesKeys",
          "d_primary_record_coverage_raw":"primary_record_coverage_raw",
          "d_primary_host_completeness":"primary_host_completeness"}
        for name,orig in keys.items():row[name]=float(y[orig]-x[orig])
        rows.append(row)
    d=pd.DataFrame(rows)
    if len(d)<100 or d.h6.nunique()<40:
        raise RuntimeError("Insufficient authentic original paired H7/H6 cochange source support")
    # Main WLS regression covariates require otherwise unrelated W and
    # Plant_N_incidence fields only for generic diagnostic metadata.
    # These have no role in the temporal fit; set them to genuine
    # time-difference equivalents and avoid model-level interpretation.
    d["W"]=0.
    d["Plant_N_incidence"]=d["d_plant_N_incidence"]
    targetpairs=[
        ("d_plant_N_incidence","d_butterfly_HNS_incidence"),
        ("d_plant_N_incidence","d_butterfly_HNS_record"),
        ("d_plant_N_record","d_butterfly_HNS_incidence")]
    out={"schema":"GESUI_EU454_TWO_PERIOD_H7_PLANT_BUTTERFLY_NUTRI_COCHANGE_V1",
         "original_454_panel_SHA256":ORIGINAL_PANEL_SHA,
         "official_NUTRI_SHA256":EIV_SHA,
         "original_original_H7period_source_eligible":451,
         "n_exact_twoperiod_H7":len(d),
         "n_true_H6":int(d.h6.nunique()),
         "H7_original_periods":[START,END],
         "posthoc_observational_not_causal":True,
         "W_exposure_static_by_design_does_not_explain_Delta":True,
         "no_actual_sewage_or_local_host_consumption_measured":True,
         "source_effort_difference_adjusted_covariates":DIFF_CONTROLS,
         "cochanges":{}}
    for independent,dependent in targetpairs:
        label=independent+"__"+dependent
        # Generic clustered fit metadata uses the same blocked source regressions,
        # but W is static zero; temporarily omit W/Plant metadata using dummy
        # only to retain interface, never as regression covariates.
        d["Plant_N_incidence"]=d[independent]
        # source_cluster_regression requires W.std only for diagnostic
        # metadata, not W as a regressor; it's intentionally 0.
        fit=source_cluster_regression(d,dependent,DIFF_CONTROLS+[independent],independent)
        cv=cv_h6(d,dependent,DIFF_CONTROLS,independent)
        out["cochanges"][label]={"adjusted_source_effort_cluster_WLS":fit,
                                 "H6_out_DeltaPlantN_added_value":cv,
                                 "H6_groupmean_raw_pearson":ci_h6_mean(d,[independent,dependent])}
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"EU454_TWO_PERIOD_H7_PLANT_HNS_COVARIATION.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    summary={"n_exact_twoperiod_H7":len(d),"n_true_H6":int(d.h6.nunique()),
             "cochanges":{}}
    for k,v in out["cochanges"].items():
        f=v["adjusted_source_effort_cluster_WLS"];cv=v["H6_out_DeltaPlantN_added_value"];b=v["H6_groupmean_raw_pearson"]
        summary["cochanges"][k]={"beta":f["adjusted_beta_per_1SD_focal"],
                                 "95":f["H6_cluster_descriptive95"],
                                 "h6_added_gain":cv["H6_equal_added_MSE_gain_pct"],
                                 "gain95":cv["fixed_out_of_fold_error_H6_resampling95_pct"],
                                 "H6_mean_r":b.get("H6_groupmean_Pearson_r"),
                                 "H6_mean_r95":b.get("H6_fixed_difference_resampling95")}
    print("GESUI_EU454_TEMPORAL_PLANT_N_BUTTERFLY_COCCHANGE="+json.dumps(
        summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())

