#!/usr/bin/env python3
"""2021 EPA NWCA standardized plant Shannon vs exact-feature DMR MONITORING nearness.

Original wetland 2021 H_ALL/FQAI/TNTP records opened ONLY AFTER the historical
2021 calendar DMR source/permit/feature gate passed, source 21N/54F in 12
state-ecoregions; observed N/P in an EXO 10km away is NOT a receptor dose.
Both all-source 2026 ECHO registered geo and historical 2021 DMR sampling are
explicitly posthoc observational, not causal wastewater-biodiversity effects.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

from gesui_us_NWCA2021_actual_site_plant_Shannon_FQAI_TP_TN_after_source_gate import (
    FILES,select_and_join,valid_rows)

DIR=Path("results/us_NWCA2021_plant_TNTP_20261110")
SAM=DIR/"NWCA2021_CAL2021_ACTUAL_DMR_TNTP_NEARFAR_STATEECO_SOURCE_ONLY_GATE.json"
OLD=DIR/"EPA_NWCA2021_NPDES_WETLAND_HALL_FQAI_AND_ACTUAL_NP_STRICT_NEAR10_FAR20_RESULTS.json"
RESPONSES=("HALL","HNAT","HALIEN","FQAIALL","FQAINAT")

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def source_id(s,visit):
    return hashlib.sha256((str(s)+"|"+str(int(visit))).encode()).hexdigest()

def model(d,outcome,weighted=False):
    q=d.loc[pd.to_numeric(d[outcome],errors="coerce").notna()].copy()
    if len(q)<50 or q.NEAR.sum()<10 or (q.NEAR==0).sum()<30 or q.PSTLCODE.nunique()<8:
        return {"STUDY_STOP":"Sparse measured plant community support in original DMR source quality cohort",
                "n":len(q),"near":int(q.NEAR.sum()),"far":int((q.NEAR==0).sum()),
                "states":int(q.PSTLCODE.nunique())}
    formula=outcome+" ~ NEAR + C(STATE_ECO) + C(WETCLSGRP) + C(HYDCLS) + LATANALYS"
    if weighted:
        w=np.asarray(q.WGTTPCORE,dtype=float)
        if np.any(w<=0) or not np.isfinite(w).all():raise ValueError("Invalid EPA NWCA survey design weight")
        w=w/np.median(w)
        fit=smf.wls(formula,data=q,weights=w).fit(
            cov_type="cluster",cov_kwds={"groups":q.PSTLCODE})
    else:
        fit=smf.ols(formula,data=q).fit(
            cov_type="cluster",cov_kwds={"groups":q.PSTLCODE})
    b=float(fit.params["NEAR"])
    ci=list(map(float,fit.conf_int().loc["NEAR"]))
    return {
       "n_firstvisit_study_wetlands":len(q),
       "n_near_2021_DMR_monitored_same_external_feature":int(q.NEAR.sum()),
       "n_control_far20km":int((q.NEAR==0).sum()),
       "states_with_both_original_source_groups":int(q.PSTLCODE.nunique()),
       "state_ecoregion_original_source_comparison_cells":int(q.STATE_ECO.nunique()),
       "original_wetland_response":outcome,
       "formula":formula,
       "fit_kind":"original_NWCA_2021_survey_area_weighted_WLS" if weighted else "site_equal_OLS",
       "beta_near_DMR_reported_POTW_minus_far_in_original_response_units":b,
       "state_cluster95":ci,
       "nominal_p_exploratory_not_preregistered":float(fit.pvalues["NEAR"]),
       "adjusted_R2":float(fit.rsquared_adj),
       "no_causal_NP_or_real_wetland_discharge_interpretation":True
    }

def measure(d):
    d=d.copy()
    d["STATE_ECO"]=d.PSTLCODE.astype("string").str.strip()+"|"+d.AGECO9.astype("string").str.strip()
    for c in ("LATANALYS","WGTTPCORE",*RESPONSES):
        d[c]=pd.to_numeric(d[c],errors="coerce")
    def val(y):
        a=d.loc[d.NEAR==1,y].dropna().to_numpy(dtype=float)
        b=d.loc[d.NEAR==0,y].dropna().to_numpy(dtype=float)
        return {
            "original_near_records":len(a),
            "original_far_records":len(b),
            "original_near_raw_mean":float(np.mean(a)) if len(a) else None,
            "original_far_raw_mean":float(np.mean(b)) if len(b) else None,
            "near_minus_far_raw_unadjusted":float(np.mean(a)-np.mean(b)) if len(a) and len(b) else None
        }
    body={
      "source_original_plots":len(d),
      "near":int((d.NEAR==1).sum()),"far":int((d.NEAR==0).sum()),
      "states":int(d.PSTLCODE.nunique()),"state_ecoregion":int(d.STATE_ECO.nunique()),
      "original_actual_2021_standardized_vegetation_raw_group_summary":{
         y:val(y) for y in RESPONSES},
      "original_2021_standardized_plants_vs_2021_reported_near_POTW_feature_models":{
         y:model(d,y,weighted=False) for y in RESPONSES},
      "original_2021_all_plant_Shannon_NWCA_area_weighted_sensitivity":
         model(d,"HALL",weighted=True),
      "area_weight_is_not_full_EPA_national_complex_stratified_survey_effect_estimator":True,
      "posthoc_DMR_presence_does_not_verify_receiving_wetland_or_sewage_nutrient_dose":True
    }
    weights=d.WGTTPCORE.to_numpy(dtype=float)
    body["source_NWCA_2021_survey_weights_dispersion"]={
        "median_area_weight":float(np.median(weights)),
        "Kish_effective_n":float(np.sum(weights)**2/np.sum(weights**2)),
        "max_weight_fraction_of_total":float(weights.max()/weights.sum()),
        "five_largest_fraction_of_total":float(np.sum(np.sort(weights)[-5:])/np.sum(weights))
    }
    return body

def run(root,directory):
    source=json.loads((root/SAM).read_text())
    old=json.loads((root/OLD).read_text())
    if old["original_STATE_ECO_source_only_comparison_population"]["source_matched_state_ecoregion_near"]!=103:
        raise ValueError("Original EPA 27-state plant plot subset changed")
    if not source["SOURCE_ONLY_no_2021_NWCA_plant_H_ALL_H_NAT_FQAI_or_measured_wetland_water_values_reopened"]:
        raise ValueError("Actual DMR source sites were selected using plant outcome")
    paired,_,audit=select_and_join(root,{k:directory/v for k,v in FILES.items()})
    if len(paired)!=568:raise ValueError("Original 568 official NWCA plant+water source join did not survive")
    paired["site_SHA"]=[source_id(a,b) for a,b in zip(paired.SITEID,paired.VISITNO)]
    if paired.site_SHA.duplicated().any():raise ValueError("Original NWCA site/visit duplicated")
    original_sites=set(paired.site_SHA)
    output={
        "schema":"gesui_US_NWCA2021_2021calendar_ECHO_NP_EXO_DMR_reported_veg_Shannon_12state_21near54far_posthoc_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_2021_DMR_exact_feature_NP_monitoring_sites_pre_response_SHA256":sha(root/SAM),
        "original_27state_240plant_NPDES2026_broad_registered_outfall_source_result_SHA256":sha(root/OLD),
        "source_568_wetland_probability_2021_firstvisits_verifed":True,
        "all_source_2021_DMR_reports_not_direct_wetland_effluent_receiving_site":True,
        "all_2026_NPDES_registered_outfalls_historical_2021_date_plausible_only":True,
        "all_new_HALL_vs_DMR_hypothesis_explicitly_POSTHOC_AFTER_old_weighted_positive_and_equal_negative":True,
        "measured_2021_wetland_TNTP_not_used_to_select_sources_or_favorable_outcomes":True
    }
    for k,v in source["original_source_DMR_2021_ANY_NP_and_BOTH_NP_strict_same_state_ecoregion_controls"].items():
        if not v["may_descriptively_fit_2021_HALL_with_monitored_POTW_proximity_ONLY_NONCAUSAL"]:
            output[k]={"STUDY_STOP":"Previously frozen identical original state/ecoregion-source support failed"}
            continue
        near=set(v["original_NWCA2021_NEAR_source_site_SHA256s"])
        far=set(v["original_NWCA2021_FAR_source_site_SHA256s"])
        if not near.isdisjoint(far) or not near.issubset(original_sites) or not far.issubset(original_sites):
            raise ValueError("Frozen DMR 2021 site sets not compatible with original 568 wetland biological survey")
        t=paired.loc[paired.site_SHA.isin(near|far)].copy()
        if len(t)!=len(near)+len(far):raise ValueError("Original source-defined DMR near/far cohort not all available")
        t["NEAR"]=t.site_SHA.isin(near).astype(int)
        t["STATE_ECO"]=t.PSTLCODE.astype("string").str.strip()+"|"+t.AGECO9.astype("string").str.strip()
        if len(t.STATE_ECO.unique())!=v["state_ecoregion_cells_with_NEAR_AND_FAR"]:
            raise ValueError("Old source-frozen state ecology cells changed before new plant outcome")
        # Fit original all 568 vegetation first-visit response on fixed source
        # locations, independent of TN/TP water sample quality.
        output[k]=measure(t)
        # Actual NWCA water TN/TP are SECONDARY, source measured at wetland,
        # not source-derived WWTP outfall concentrations; do not mediate.
        actual=valid_rows(t)
        if len(actual)>0:
            near_actual=actual.loc[actual.NEAR==1]
            far_actual=actual.loc[actual.NEAR==0]
            output[k]["original_2021_actual_wetland_water_TNTP_quality_supported_sites"]={
                "n_source_only_wetlands_with_measured_TNTP_and_HALL":len(actual),
                "near":len(near_actual),"far":len(far_actual),
                "median_wetland_TN_mg_N_L_NEAR":float(near_actual.NTLRESULT.median()) if len(near_actual) else None,
                "median_wetland_TN_mg_N_L_FAR":float(far_actual.NTLRESULT.median()) if len(far_actual) else None,
                "median_wetland_TP_mg_P_L_NEAR":float(near_actual.PTLRESULT.median()/1000) if len(near_actual) else None,
                "median_wetland_TP_mg_P_L_FAR":float(far_actual.PTLRESULT.median()/1000) if len(far_actual) else None,
                "2021_wetland_nutrient_concentration_NOT_established_as_2021_source_outfall_discharge":True
            }
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nwca-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root,a.nwca_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
       k:({"source":{z:v[z] for z in ("source_original_plots","near","far","states","state_ecoregion")},
           "HALL":v["original_2021_standardized_plants_vs_2021_reported_near_POTW_feature_models"]["HALL"],
           "HALL_design_weighted":v["original_2021_all_plant_Shannon_NWCA_area_weighted_sensitivity"],
           "wetland_TNTP":v.get("original_2021_actual_wetland_water_TNTP_quality_supported_sites")}
          if "source_original_plots" in v else v)
       for k,v in out.items() if k.startswith("DMR_actual2021_")
    },indent=2))

if __name__=="__main__":main()
