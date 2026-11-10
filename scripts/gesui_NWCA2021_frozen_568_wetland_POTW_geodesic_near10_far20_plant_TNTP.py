#!/usr/bin/env python3
"""EPA NWCA2021 standardized plants × strict NPDES POTW feature proximity.

STRICT post-source 2026 EPA registered outfall geography already frozen on
original 568 probability wetlands. 2021 H_ALL/H_NAT/H_ALIEN/FQAI and actual
TN/TP only opened AFTER source n128 near10/137 comparable far matched gate.
This is source-distance association, NOT hydrologically confirmed WWTP dose.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import spearmanr

from gesui_us_NWCA2021_actual_site_plant_Shannon_FQAI_TP_TN_after_source_gate import (
    FILES,select_and_join,valid_rows)

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
EXPOSURE=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
PREV=D/"EPA_NWCA2021_REAL_PLANT_SHANNON_TNTP_WETLANDS_STATE_ECOREGION_QA.json"
PRIMARY_Y="HALL"
SECONDARY_Y=("HNAT","HALIEN","FQAIALL","FQAINAT")
MIN_NEAR_SOURCE=60
MIN_FAR_SOURCE=60
MIN_STATES_SOURCE=10
MIN_CELLS_SOURCE=15
MIN_WATER_NEAR=50
MIN_WATER_FAR=50
MIN_WATER_STATES=10

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def site_hash(a,b):
    return hashlib.sha256((str(a)+"|"+str(int(b))).encode()).hexdigest()

def fit(d,y,regressors,design_weight=False):
    usable=d.loc[pd.to_numeric(d[y],errors="coerce").notna()].copy()
    if len(usable)<80:return None
    # SOURCE STRATA fixed before opening any new biodiversity endpoint.
    formula=y+"~"+regressors+"+C(STATE_ECO)+C(WETCLSGRP)+C(HYDCLS)+LATANALYS"
    if design_weight:
        w=usable.WGTTPCORE.to_numpy(dtype=float)
        if np.any(w<=0) or not np.isfinite(w).all():
            raise ValueError("Probability design weight invalid on source-selected sites")
        w=w/np.median(w)
        m=smf.wls(formula,data=usable,weights=w).fit(cov_type="cluster",cov_kwds={"groups":usable.PSTLCODE})
    else:
        m=smf.ols(formula,data=usable).fit(cov_type="cluster",cov_kwds={"groups":usable.PSTLCODE})
    results={
        "model_formula":formula,
        "original_same_state_ecoregion_matched_sites_n":len(usable),
        "independent_original_US_state_cluster_n":int(usable.PSTLCODE.nunique()),
        "original_state_ecoregion_matching_strata_n":int(usable.STATE_ECO.nunique()),
        "number_NWCA_plant_near10_sites":int((usable.NEAR==1).sum()),
        "number_NWCA_plant_far20_sites":int((usable.NEAR==0).sum()),
        "R2":float(m.rsquared),
        "model_type":"original_NWCA_design_weighted_WLS" if design_weight else "site_equal_observational_OLS",
        "source_outfall_near10_beta_and_state_cluster95":None,
        "posthoc_2026_permit_location_not_confirmed_2021_effluent_exposure":True,
    }
    if "NEAR" in m.params:
        x=float(m.params["NEAR"])
        results["source_outfall_near10_beta_and_state_cluster95"]={
            "beta_near10_minus_far20_in_reported_response_units":x,
            "original_US_state_cluster95":list(map(float,m.conf_int().loc["NEAR"])),
            "nominal_p_not_confirmatory":float(m.pvalues["NEAR"])
        }
    for par in ("log10_TN","log10_TP"):
        if par in m.params:
            results[par]={
                "beta_on_log10_nutrient":float(m.params[par]),
                "state_cluster95":list(map(float,m.conf_int().loc[par]))
            }
    return results

def main_data(root,source_dir):
    fx=json.loads((root/EXPOSURE).read_text(encoding="utf-8"))
    if fx["schema"]!="gesui_US_EPA_NWCA2021_568source_NPDES2026_reported_sewage_outfall_GEO_2021permitDates_sourceONLY_v01":
        raise ValueError("Strict municipal WWTP near10 source not frozen")
    gate=fx["source_precommitted_GEO_eligibility_gate"]
    if not gate["permitted_2021_outfall_distance_association_source_support_PASS"]:
        return fx,None,"ORIGINAL_568_SOURCE_OUTFALL_GEOMETRY_GATE_FAILED"
    if not fx["no_original_NWCA_H_ALL_FQAI_or_TNTP_outcomes_accessed_in_source_stage"]:
        raise ValueError("Original WWTP source stage was not blind to new 2021 plant response")
    original=json.loads((root/PREV).read_text(encoding="utf-8"))
    if original["source_568_veg_and_waterchem_joined_site_rows"]!=568 or original["source_joint_TNTP_and_Shannon_nonmissing_quality_sites"]!=499:
        raise ValueError("Original standardized 2021 NWCA main source/output changed")
    files={k:source_dir/p for k,p in FILES.items()}
    paired,_,source_audit=select_and_join(root,files)
    if len(paired)!=568:raise ValueError("Original official NWCA site code join changed")
    pairfx={row["NWCA2021_source_site_visit_SHA256"]:row
            for row in fx["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]}
    if len(pairfx)!=568:raise ValueError("NPDES geometrical stage changed exact 568 source wetland sites")
    x=paired.copy()
    x["SITE_SHA"]=[site_hash(a,b) for a,b in zip(x.SITEID,x.VISITNO)]
    if set(x.SITE_SHA)!=set(pairfx):raise ValueError("Original source-only registered-outfall locations changed the biological site population")
    x["NEAR_SOURCE"]=x.SITE_SHA.map(lambda s:int(pairfx[s]["outfall_rows_within_10km"]>0))
    x["FAR_SOURCE"]=x.SITE_SHA.map(lambda s:int(
        pairfx[s]["nearest_outfall_geodesic_km_within25km_else_null"] is None or
        pairfx[s]["nearest_outfall_geodesic_km_within25km_else_null"]>=20.))
    x["NEAR5_SOURCE"]=x.SITE_SHA.map(lambda s:int(pairfx[s]["outfall_rows_within_5km"]>0))
    if int(x.NEAR_SOURCE.sum())!=128:
        raise ValueError("Pre-outcome source group 128 near10 source eligibility changed")
    if int(x.FAR_SOURCE.sum())!=311:
        raise ValueError("Pre-outcome source group 311 far20 source eligibility changed")
    if ((x.NEAR_SOURCE==1)&(x.FAR_SOURCE==1)).any():
        raise ValueError("Same original wetland source site classified both near and far")
    x["STATE_ECO"]=x.PSTLCODE.astype("string").str.strip()+"|"+x.AGECO9.astype("string").str.strip()
    strata=set(x.loc[x.NEAR_SOURCE==1,"STATE_ECO"])&set(x.loc[x.FAR_SOURCE==1,"STATE_ECO"])
    matched=x.loc[x.STATE_ECO.isin(strata)&((x.NEAR_SOURCE==1)|(x.FAR_SOURCE==1))].copy()
    matched["NEAR"]=matched.NEAR_SOURCE.astype(int)
    if matched.SITE_SHA.duplicated().any():raise ValueError("Some original probability site duplicated in NPDES analysis")
    sizes=Counter(matched.NEAR)
    regions=matched.PSTLCODE.nunique()
    support=bool(sizes[1]>=MIN_NEAR_SOURCE and sizes[0]>=MIN_FAR_SOURCE and
                 regions>=MIN_STATES_SOURCE and len(strata)>=MIN_CELLS_SOURCE)
    return fx,{"original_matched_sites":matched,"all_site_source":x,
               "strata":sorted(strata),"n_source":dict(sizes),"n_regions":int(regions)},None if support else "BIOLOGICAL_SAME_REGION_SAMPLE_SUPPORT_FAIL"

def q(v):
    a=pd.to_numeric(v,errors="coerce").dropna().to_numpy(dtype=float)
    if not len(a):return None
    return {"count":len(a),"mean":float(np.mean(a)),"median":float(np.median(a)),
            "q25":float(np.quantile(a,.25)),"q75":float(np.quantile(a,.75))}

def run(root,source_dir):
    fx,eligible,reason=main_data(root,source_dir)
    out={
        "schema":"gesui_US_NWCA2021_568_probability_plants_NPDES_strict_sewage2026_outfall_near10_far20_match_2021_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_568_sites_ECHO_2026_source_only_geography_SHA256":sha(root/EXPOSURE),
        "original_2021_NWCA_Shannon_actual_TNTP_old_model_source_SHA256":sha(root/PREV),
        "original_EPA_NPDES_outfall_2026_snapshot_SHA256":fx["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"],
        "original_2021_municipal_sewage_source_only_geo_gate":fx["source_precommitted_GEO_eligibility_gate"],
        "geographic_NPDES_source_support_valid":reason is None,
        "stop_reason":reason,
        "2026_original_feature_registered_with_2021_permit_date_NOT_verified_2021_plant_nutrient_or_wetland_hydrological_exposure":True,
        "posthoc_after_original_NWCA_TN_Shannon_negative_result":True,
        "not_randomized_facility_location_or_pollution_exposure":True
    }
    if reason:return out
    d=eligible["original_matched_sites"].copy()
    for k in (PRIMARY_Y,*SECONDARY_Y):
        d[k]=pd.to_numeric(d[k],errors="coerce")
    d.LATANALYS=pd.to_numeric(d.LATANALYS,errors="coerce")
    d.WGTTPCORE=pd.to_numeric(d.WGTTPCORE,errors="coerce")
    if len(d)<100 or d[PRIMARY_Y].notna().sum()<100:
        raise ValueError("Original NWCA H_ALL vegetation not present after source frame")
    out["original_STATE_ECO_source_only_comparison_population"]={
        "all568_source_sites":568,
        "source_geographic_near10":128,
        "source_geographic_far20":311,
        "source_matched_state_ecoregion_near":int((d.NEAR==1).sum()),
        "source_matched_state_ecoregion_far":int((d.NEAR==0).sum()),
        "original_US_states_after_state_ecoregion_matching":int(d.PSTLCODE.nunique()),
        "source_original_state_ecoregion_strata_count":int(d.STATE_ECO.nunique()),
        "original_state_ecoregion_names":eligible["strata"],
        "original_near5_subset_in_same_NWCA_region_source_rows":int((d.NEAR5_SOURCE==1).sum()),
        "source_original_2021_site_status_frozen_before_HALL_WWTP_outcome":True
    }
    desc={}
    for y in (PRIMARY_Y,*SECONDARY_Y):
        desc[y]={"near10":q(d.loc[d.NEAR==1,y]),
                 "far20":q(d.loc[d.NEAR==0,y])}
    out["original_standardized_plants_2021_descriptive_wetland_by_WWTP_source_proximity"]=desc
    models={}
    for y in (PRIMARY_Y,*SECONDARY_Y):
        models[y]=fit(d,y,"NEAR")
    models["HALL_design_weighted_sensitivity"]=fit(d,PRIMARY_Y,"NEAR",design_weight=True)
    out["original_wetland_standardized_plant_WWTP_geodesic_near10_far20_country_model"]=models
    # Biological alternative radius was source fixed as 5km; no optimizing sign.
    near5=d.loc[(d.NEAR5_SOURCE==1)|(d.NEAR==0)].copy()
    near5["NEAR"]=near5.NEAR5_SOURCE.astype(int)
    group=Counter(near5.NEAR)
    if group[1]>=20 and near5.loc[near5.NEAR==1,"PSTLCODE"].nunique()>=6:
        out["source_fixed_strict_5km_sensitivity"]={
            "near":group[1],"far":group[0],
            "states_with_source_near5":int(near5.loc[near5.NEAR==1,"PSTLCODE"].nunique()),
            "HALL_original_sites":fit(near5,"HALL","NEAR")
        }
    else:
        out["source_fixed_strict_5km_sensitivity"]={
            "near":group[1],"far":group[0],"STOP_original_precommitted_20_near_6states":True}
    # Real measured nutrient subset, NOT treated as effluent chemistry.
    actual=valid_rows(eligible["all_site_source"])
    if len(actual)!=499:raise ValueError("Original independently verified 499 actual TN/TP NWCA sites changed")
    actual_matched=actual.loc[actual.STATE_ECO.isin(eligible["strata"]) &
          ((actual.NEAR_SOURCE==1)|(actual.FAR_SOURCE==1))].copy()
    actual_matched["NEAR"]=actual_matched.NEAR_SOURCE.astype(int)
    nnear=int((actual_matched.NEAR==1).sum())
    nfar=int((actual_matched.NEAR==0).sum())
    waters={"n_original_499_quality_plant_and_actual_nutrients":499,
            "n_original_499_quality_joint_source_near":nnear,
            "n_original_499_quality_joint_source_far":nfar,
            "states":int(actual_matched.PSTLCODE.nunique()),
            "measured_water_NP_not_identified_as_NPDES_effluent":True}
    if nnear>=MIN_WATER_NEAR and nfar>=MIN_WATER_FAR and actual_matched.PSTLCODE.nunique()>=MIN_WATER_STATES:
        actual_matched["log10_TN"]=np.log10(pd.to_numeric(actual_matched.NTLRESULT))
        actual_matched["log10_TP"]=np.log10(pd.to_numeric(actual_matched.PTLRESULT)/1000)
        waters["actual_TN_mgN_L"]=fit(actual_matched,"log10_TN","NEAR")
        waters["actual_TP_mgP_L"]=fit(actual_matched,"log10_TP","NEAR")
        waters["HALL_with_TNTP_conditional_not_causal_mediation"]=fit(
            actual_matched,"HALL","NEAR+log10_TN+log10_TP")
        waters["raw_measured_NP_by_candidate_source_proximity"]={
            "near10": {"TN_mgN_L":q(actual_matched.loc[actual_matched.NEAR==1,"NTLRESULT"]),
                        "TP_mgP_L":q(actual_matched.loc[actual_matched.NEAR==1,"PTLRESULT"]/1000.)},
            "far20":{"TN_mgN_L":q(actual_matched.loc[actual_matched.NEAR==0,"NTLRESULT"]),
                       "TP_mgP_L":q(actual_matched.loc[actual_matched.NEAR==0,"PTLRESULT"]/1000.)}
        }
    else:waters["STOP_insufficient_actual_NP_WWTP_proximity_strata"]=True
    out["original_2021_NWCA_499_prior_quality_measured_NP_WWTP_source_proximity_posthoc"]=waters
    out["original_NPDES_near10_admin_proximity_never_direct_receiving_water_plant_nutrient_exposure"]=True
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--nwca-source-dir",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    a=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root,a.nwca_source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":x.get("original_STATE_ECO_source_only_comparison_population"),
        "HALL":x.get("original_wetland_standardized_plant_WWTP_geodesic_near10_far20_country_model",{}).get("HALL"),
        "FQAI":x.get("original_wetland_standardized_plant_WWTP_geodesic_near10_far20_country_model",{}).get("FQAIALL"),
        "5km":x.get("source_fixed_strict_5km_sensitivity"),
        "actualNandP":x.get("original_2021_NWCA_499_prior_quality_measured_NP_WWTP_source_proximity_posthoc")
    },indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
