#!/usr/bin/env python3
"""Independent SOURCE/SAMPLE + NumPy OLS/WLS QA: TRUE External Outfall only.

Same original 568 EPA NWCA2021 first-visit probable wetlands, 2026 EPA ECHO
EXO-only source group frozen before corrected vegetation H_ALL response.
This secondary implementation does not call original statsmodels fit code.
External outfall geographic proximity does NOT establish 2021 N/P exposure.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

from gesui_us_NWCA2021_actual_site_plant_Shannon_FQAI_TP_TN_after_source_gate import (
    FILES,select_and_join,valid_rows)
from audit_gesui_NWCA2021_POTW_near10_HALL_designweights_independent import (
    independent_fit,qd)

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SRC=D/"EPA_NWCA2021_NPDES_EXTERNAL_OUTFALL_ONLY_568SITE_SOURCE_GATE.json"
BIA=D/"EPA_NWCA2021_EXOONLY_239SITE_PLANTS_TNTP_EXPOSURE_SENSITIVITY.json"
PREV=D/"EPA_NWCA2021_NPDES_WETLAND_HALL_FQAI_AND_ACTUAL_NP_STRICT_NEAR10_FAR20_RESULTS.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sid(a,b):return hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()
def numeric_eq(a,b,where,tol=1e-6):
    if a is None or b is None or not np.isfinite(a) or not np.isfinite(b) or abs(a-b)>tol:
        raise ValueError(where+" expected "+str(b)+" got "+str(a))

def run(root,directory):
    geo=json.loads((root/SRC).read_text())
    out=json.loads((root/BIA).read_text())
    original=json.loads((root/PREV).read_text())
    if out["external_source_only_SHA256"]!=sha(root/SRC):
        raise ValueError("Biological source model did not use frozen 568 External Outfall source")
    if out["old_original_EPA_source_type_unrestricted_GEO_SHA256"]!=geo["original_EPA_NWCA568_SOURCE_only_SHA256"]:
        raise ValueError("Old NPDES geometry not identical in new EXO study")
    if not geo["frozen_before_corrected_source_type_biological_comparison_GATES"]["may_open_new_externalOutfall_source_sensitive_2021_biodiversity_effect"]:
        raise ValueError("The original source-only EXO comparison support gate did not pass")
    files={k:directory/p for k,p in FILES.items()}
    observations,_,base=select_and_join(root,files)
    if len(observations)!=568:raise ValueError("Original US EPA NWCA source no longer 568 sites")
    source={x["NWCA2021_source_site_visit_SHA256"]:x
        for x in geo["all_original_568_site_hashes_and_external_OUTFALL_source_receipts_no_plant_or_nutrient_outcomes"]}
    if len(source)!=568:raise ValueError("Original EXO source site IDs duplicated")
    obs=observations.copy()
    obs["SITE_SHA"]=[sid(a,b) for a,b in zip(obs.SITEID,obs.VISITNO)]
    if set(obs.SITE_SHA)!=set(source):raise ValueError("New source group HALL changed original wetland site population")
    obs["SOURCE_NEAR"]=obs.SITE_SHA.map(lambda k:int(source[k]["outfall_rows_within_10km"]>0))
    obs["SOURCE_FAR"]=obs.SITE_SHA.map(lambda k:int(
        source[k]["nearest_outfall_geodesic_km_within25km_else_null"] is None or
        source[k]["nearest_outfall_geodesic_km_within25km_else_null"]>=20))
    if obs.SOURCE_NEAR.sum()!=124 or obs.SOURCE_FAR.sum()!=314:raise ValueError("EXO 124 source near and 314 far not reproduced")
    obs["STATE_ECO"]=obs.PSTLCODE.astype("string").str.strip()+"|"+obs.AGECO9.astype("string").str.strip()
    strata=set(obs.loc[obs.SOURCE_NEAR==1,"STATE_ECO"])&set(obs.loc[obs.SOURCE_FAR==1,"STATE_ECO"])
    d=obs.loc[obs.STATE_ECO.isin(strata)&
          (obs.SOURCE_NEAR.eq(1)|obs.SOURCE_FAR.eq(1))].copy()
    d["NEAR"]=d.SOURCE_NEAR.astype(int)
    if Counter(d.NEAR)!=Counter({1:99,0:140}) or len(d)!=239 or len(strata)!=32 or d.PSTLCODE.nunique()!=27:
        raise ValueError("Frozen 239 original EXO near/far 32-source strata not reproduced independently")
    if d.SITE_SHA.duplicated().any():raise ValueError("Duplicate original NWCA source plot")
    for col in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        d[col]=pd.to_numeric(d[col],errors="coerce")
    o=out["same_original_standardized_wetlands_2021_plant_models"]
    tests={}
    for resp in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        numeric=independent_fit(d,resp,None)
        initial=o[resp+"_site_equal"]["source_outfall_near10_beta_and_state_cluster95"]
        numeric_eq(numeric["independent_numpy_near10_minus_far20_coefficient"],
                   initial["beta_near10_minus_far20_in_reported_response_units"],resp+" site β")
        for q,r in zip(numeric["independent_state_cluster_sandwich95"],initial["original_US_state_cluster95"]):
            numeric_eq(q,r,resp+" site state CI")
        tests[resp+"_site_equal"]=numeric
    weighted=independent_fit(d,"HALL",d.WGTTPCORE.astype(float))
    i=o["HALL_2021_original_design_weighted"]["source_outfall_near10_beta_and_state_cluster95"]
    numeric_eq(weighted["independent_numpy_near10_minus_far20_coefficient"],
               i["beta_near10_minus_far20_in_reported_response_units"],"H_ALL acreage weighted β")
    for q,r in zip(weighted["independent_state_cluster_sandwich95"],i["original_US_state_cluster95"]):
        numeric_eq(q,r,"H_ALL acreage weighted state CI")
    tests["HALL_original_design_weighted"]=weighted
    weight=d.WGTTPCORE.astype(float)
    clipped={}
    for pct in (.95,.99):
        limit=float(np.quantile(weight,pct))
        tmp=independent_fit(d,"HALL",weight.clip(upper=limit))
        clipped["weight_capped_at_"+str(int(pct*100))+"pct"]={
            "cap_WGT_TP_CORE":limit,
            "site_count_above_cap":int((weight>limit).sum()),
            "original_239_site_cluster_coefficient":tmp,
            "not_official_NWCA_design_estimator":True
        }
    loostate={}
    for state in sorted(d.PSTLCODE.astype(str).unique()):
        sub=d.loc[d.PSTLCODE.astype(str)!=state]
        w=independent_fit(sub,"HALL",sub.WGTTPCORE.astype(float))
        v=independent_fit(sub,"HALL",None)
        loostate[str(state)]={"weighted_beta":w["independent_numpy_near10_minus_far20_coefficient"],
                              "site_equal_beta":v["independent_numpy_near10_minus_far20_coefficient"]}
    # Prior true H_ALL with same-source measured TN/TP; audit positive TP
    # source association separately from actual WWTP discharge causal claim.
    real=valid_rows(obs)
    if len(real)!=499:raise ValueError("Original water quality 499 source site population changed")
    act=real.loc[real.STATE_ECO.isin(strata)&
          (real.SOURCE_NEAR.eq(1)|real.SOURCE_FAR.eq(1))].copy()
    act["NEAR"]=act.SOURCE_NEAR.astype(int)
    if len(act)!=204 or Counter(act.NEAR)!=Counter({1:81,0:123}):
        raise ValueError("Actual original water quality matched sample changed after EXO source")
    act["log10_TP"]=np.log10(act.PTLRESULT/1000.)
    act["log10_TN"]=np.log10(act.NTLRESULT)
    chemical={}
    refsource=out["original_measured_2021_wetland_N_and_P_sample_by_externalOutfall_proximity"]
    for y,k in (("log10_TN","TN_geodesic_near10_vs_far20"),("log10_TP","TP_geodesic_near10_vs_far20")):
        v=independent_fit(act,y,None)
        r=refsource[k]["source_outfall_near10_beta_and_state_cluster95"]
        numeric_eq(v["independent_numpy_near10_minus_far20_coefficient"],
                   r["beta_near10_minus_far20_in_reported_response_units"],y+" real nutrient β")
        for a,b in zip(v["independent_state_cluster_sandwich95"],r["original_US_state_cluster95"]):
            numeric_eq(a,b,y+" real nutrient CI")
        chemical[y]=v
    return {
        "schema":"gesui_NWCA2021_239plots_exact_EXTERNAL_OUTFALL_ONLY_NP_Shannon_independent_numpy_USstate_QA_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_EXO_source_only_SHA256":sha(root/SRC),
        "source_original_2021_239plant_sensitivity_model_SHA256":sha(root/BIA),
        "original_previous_2021_240plant_source_model_SHA256":sha(root/PREV),
        "original_568_site_biodiversity_source_unchanged":True,
        "US_state_and_ecoregion_source_overlap_239_exact":True,
        "matched_US_states":27,"matched_STATE_ECO_strata":32,
        "source_near10":99,"source_far20":140,
        "independent_2021_all_five_NWCA_plant_endpoints_and_HALL_original_design_weight_beta_CI":tests,
        "independent_original_2021_TN_TP_exact_204site_observational_concentrations":chemical,
        "original_area_weight_WGT_TP_CORE_distribution":{
            "all239":qd(weight.to_numpy()),
            "near99":qd(d.loc[d.NEAR.eq(1),"WGTTPCORE"].to_numpy()),
            "far140":qd(d.loc[d.NEAR.eq(0),"WGTTPCORE"].to_numpy())
        },
        "descriptive_truncated_survey_area_weight_not_original_national_area_estimator":clipped,
        "leave_original_one_US_state_out_2021_HALL":{
            "n_states_dropped":len(loostate),
            "site_equal_negative_count":sum(v["site_equal_beta"]<0 for v in loostate.values()),
            "area_weighted_positive_count":sum(v["weighted_beta"]>0 for v in loostate.values()),
            "site_beta_range":[min(v["site_equal_beta"] for v in loostate.values()),
                               max(v["site_equal_beta"] for v in loostate.values())],
            "weighted_beta_range":[min(v["weighted_beta"] for v in loostate.values()),
                                   max(v["weighted_beta"] for v in loostate.values())]
        },
        "actual2021_nutrient_not_confirmed_wastewater_source":True,
        "observational_source_outfall_2026_registered_not_historical_wetland_dose":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    res=run(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(res,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "main":{k:{"beta":v["independent_numpy_near10_minus_far20_coefficient"],
                    "ci":v["independent_state_cluster_sandwich95"]} for k,v in res["independent_2021_all_five_NWCA_plant_endpoints_and_HALL_original_design_weight_beta_CI"].items()},
        "TNTP":{k:{"beta":v["independent_numpy_near10_minus_far20_coefficient"],
                      "ci":v["independent_state_cluster_sandwich95"]} for k,v in res["independent_original_2021_TN_TP_exact_204site_observational_concentrations"].items()},
        "weight_source":res["original_area_weight_WGT_TP_CORE_distribution"],
        "weighted_trim":res["descriptive_truncated_survey_area_weight_not_original_national_area_estimator"],
        "state":res["leave_original_one_US_state_out_2021_HALL"]
    },indent=2))
if __name__=="__main__":
    main()
