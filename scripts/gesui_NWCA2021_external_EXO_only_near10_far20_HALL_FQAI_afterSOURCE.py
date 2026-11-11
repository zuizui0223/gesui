#!/usr/bin/env python3
"""Observed standardized US wetland plants vs exactly classified EXTERNAL outfalls.

SOURCE FIRST: the independently acquired 568-site/49,549 External Outfall
subtype classifier was SHA-frozen BEFORE this sensitivity touches H_ALL/FQAI.
2026 registration and 2021 plausible permit dates NOT proven 2021 discharge
or receiving-wetland exposure. Preserves prior 2021 sign-unfavorable results.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

from gesui_us_NWCA2021_actual_site_plant_Shannon_FQAI_TP_TN_after_source_gate import (
    FILES,select_and_join,valid_rows)
from gesui_NWCA2021_frozen_568_wetland_POTW_geodesic_near10_far20_plant_TNTP import (
    fit,site_hash,q)

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
EXT=D/"EPA_NWCA2021_NPDES_EXTERNAL_OUTFALL_ONLY_568SITE_SOURCE_GATE.json"
OLDG=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
OLDM=D/"EPA_NWCA2021_NPDES_WETLAND_HALL_FQAI_AND_ACTUAL_NP_STRICT_NEAR10_FAR20_RESULTS.json"
MIN_NEAR=60
MIN_FAR=60
MIN_STATES=10
MIN_CELLS=15

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def prepare(root,dir):
    source=json.loads((root/EXT).read_text())
    oldg=json.loads((root/OLDG).read_text())
    prev=json.loads((root/OLDM).read_text())
    if source["schema"]!="gesui_NWCA2021_568_source_2026EPA_NPDES_exact_EXTERNAL_OUTFALL_only_pre_plant_response_gate_v01":
        raise ValueError("True EXO source-only audit not available")
    if not source["frozen_before_corrected_source_type_biological_comparison_GATES"]["may_open_new_externalOutfall_source_sensitive_2021_biodiversity_effect"]:
        return source,oldg,prev,None,"PRECOMMITTED_EXTERNAL_OUTFALL_SOURCE_BIOLOGICAL_SUPPORT_FAIL"
    if source["original_EPA_NWCA568_SOURCE_only_SHA256"]!=sha(root/OLDG):
        raise ValueError("Original 568 wetland source-only EPA ECHO data digest altered")
    if source["original_NPDES_external_source_ZIP_SHA256"]!=oldg["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]:
        raise ValueError("NPDES 2026-11-08 official source archive differs")
    if not source["all_568_frozen_nwca_plant_water_locations_source_complete"]:
        raise ValueError("Original 2021 NWCA source site frame not frozen")
    paired,_,original=select_and_join(root,{k:dir/name for k,name in FILES.items()})
    if len(paired)!=568:
        raise ValueError("Original 568 observed species and water source visits changed")
    categories={x["NWCA2021_source_site_visit_SHA256"]:x
                for x in source["all_original_568_site_hashes_and_external_OUTFALL_source_receipts_no_plant_or_nutrient_outcomes"]}
    if len(categories)!=568:raise ValueError("Duplicate 568 source site source-type rows")
    frame=paired.copy()
    frame["SITE_SHA"]=[site_hash(a,b) for a,b in zip(frame.SITEID,frame.VISITNO)]
    if set(frame.SITE_SHA)!=set(categories):
        raise ValueError("Original 2021 NWCA wetland site IDs are not the frozen EXO group")
    frame["NEAR_SOURCE"]=frame.SITE_SHA.map(lambda k:int(categories[k]["outfall_rows_within_10km"]>0))
    frame["FAR_SOURCE"]=frame.SITE_SHA.map(lambda k:int(
        categories[k]["nearest_outfall_geodesic_km_within25km_else_null"] is None or
        categories[k]["nearest_outfall_geodesic_km_within25km_else_null"]>=20.))
    frame["NEAR5_SOURCE"]=frame.SITE_SHA.map(lambda k:int(categories[k]["outfall_rows_within_5km"]>0))
    if frame.NEAR_SOURCE.sum()!=124 or frame.FAR_SOURCE.sum()!=314:
        raise ValueError("Original external-only 124 near /314 far counts changed")
    frame["STATE_ECO"]=frame.PSTLCODE.astype("string").str.strip()+"|"+frame.AGECO9.astype("string").str.strip()
    strata=sorted(set(frame.loc[frame.NEAR_SOURCE==1,"STATE_ECO"])&
                  set(frame.loc[frame.FAR_SOURCE==1,"STATE_ECO"]))
    matched=frame.loc[frame.STATE_ECO.isin(strata)&
                      ((frame.NEAR_SOURCE==1)|(frame.FAR_SOURCE==1))].copy()
    matched["NEAR"]=matched.NEAR_SOURCE.astype(int)
    if (matched.NEAR==1).sum()!=99 or (matched.NEAR==0).sum()!=140 or matched.PSTLCODE.nunique()!=27 or len(strata)!=32:
        raise ValueError("Source-only original same state/ecoregion EXO comparison population changed")
    if matched.SITE_SHA.duplicated().any():raise ValueError("Repeated original source site in EXO response")
    return source,oldg,prev,{"all":frame,"matched":matched,"strata":strata},None

def study(root,dir):
    source,oldgeo,prev,samples,stop=prepare(root,dir)
    out={
        "schema":"gesui_US_NWCA2021_ExternalOutfall2026_real2021_wetland_Shannon_FQAI_near10_far20_AFTER_SOURCE_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "external_source_only_SHA256":sha(root/EXT),
        "old_original_EPA_source_type_unrestricted_GEO_SHA256":sha(root/OLDG),
        "old_any_permitted_feature_NWCA2021_biological_result_SHA256":sha(root/OLDM),
        "original_registered_permit_date_plausible_2021_but_not_actual_historical_2021_effluent":True,
        "original_2021_actual_nutrient_site_concentrations_NOT_proven_derived_from_NPDES_external_outfall":True,
        "frozen_SOURCE_biological_analysis_permitted":stop is None,
        "STOP_reason":stop,
        "different_source_type_posthoc_after_original_weighted_positive_known":True,
    }
    if stop is not None:return out
    d=samples["matched"].copy()
    d["WGTTPCORE"]=pd.to_numeric(d.WGTTPCORE,errors="coerce")
    d["LATANALYS"]=pd.to_numeric(d.LATANALYS,errors="coerce")
    for name in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        d[name]=pd.to_numeric(d[name],errors="coerce")
    if d.HALL.isna().any():
        raise ValueError("Source-matched original NWCA recorded Shannon H_ALL unexpectedly missing")
    out["original_frozen_external_outfall_only_source_sample"]={
        "original_568_same_site_NWCA2021_population":568,
        "actual_ExternalOutfall_registered_near10_source_N":124,
        "actual_ExternalOutfall_registered_far20_source_N":314,
        "same_state_and_eco_matched_source_near10":99,
        "same_state_and_eco_matched_source_far20":140,
        "original_matched_US_states":int(d.PSTLCODE.nunique()),
        "original_state_ecoregion_strata":int(d.STATE_ECO.nunique()),
        "source_2021_community_plots_have_original_probability_design_weight":True,
    }
    models={}
    descript={}
    for y in ("HALL","HNAT","HALIEN","FQAIALL","FQAINAT"):
        models[y+"_site_equal"]=fit(d,y,"NEAR")
        descript[y]={
            "near10":q(d.loc[d.NEAR==1,y]),
            "far20":q(d.loc[d.NEAR==0,y])
        }
    models["HALL_2021_original_design_weighted"]=fit(d,"HALL","NEAR",design_weight=True)
    out["same_original_standardized_wetlands_2021_plant_models"]=models
    out["source_standardized_plant_distribution_by_externalOutfall_proximity"]=descript
    orig=prev["original_wetland_standardized_plant_WWTP_geodesic_near10_far20_country_model"]
    out["original_old_all_feature_same_region_240plot_site_vs_weighted"]={
        "old_equal_240_site_HALL_beta":orig["HALL"]["source_outfall_near10_beta_and_state_cluster95"]["beta_near10_minus_far20_in_reported_response_units"],
        "old_weighted_240_site_HALL_beta":orig["HALL_design_weighted_sensitivity"]["source_outfall_near10_beta_and_state_cluster95"]["beta_near10_minus_far20_in_reported_response_units"],
        "old_near":103,"old_far":137,
        "source_type_corrected_new_near":99,"source_type_corrected_new_far":140
    }
    near5=d.loc[(d.NEAR5_SOURCE==1)|(d.NEAR==0)].copy()
    near5["NEAR"]=near5.NEAR5_SOURCE.astype(int)
    if int((near5.NEAR==1).sum())>=20 and near5.loc[near5.NEAR==1,"PSTLCODE"].nunique()>=6:
        out["source_fixed_5km_externalOutfall_sensitivity"]={
            "near5":int((near5.NEAR==1).sum()),
            "far20":int((near5.NEAR==0).sum()),
            "site_equal":fit(near5,"HALL","NEAR"),
            "design_weighted":fit(near5,"HALL","NEAR",design_weight=True)
        }
    else:out["source_fixed_5km_externalOutfall_sensitivity"]={"STOP_precommitted_min20near_6states":True}
    # Reconciliation only, not a causal mediation model.
    measured=valid_rows(samples["all"])
    if len(measured)!=499:raise ValueError("Original 499 measured nutrient/plant quality sites changed")
    matched=measured.loc[measured.STATE_ECO.isin(samples["strata"])&
                         ((measured.NEAR_SOURCE==1)|(measured.FAR_SOURCE==1))].copy()
    matched["NEAR"]=matched.NEAR_SOURCE.astype(int)
    nn=int((matched.NEAR==1).sum());nf=int((matched.NEAR==0).sum())
    out["original_measured_2021_wetland_N_and_P_sample_by_externalOutfall_proximity"]={
        "old_actual_TNTP_2021_valid_N":499,
        "new_external_source_matched_quality_N":len(matched),
        "near10":nn,"far20":nf,
        "not_measured_effluent_load_or_actual_downstream_receptor":True
    }
    if nn>=50 and nf>=50 and matched.PSTLCODE.nunique()>=10:
        matched["log10_TN"]=np.log10(matched.NTLRESULT)
        matched["log10_TP"]=np.log10(matched.PTLRESULT/1000.)
        m=out["original_measured_2021_wetland_N_and_P_sample_by_externalOutfall_proximity"]
        m["TN_geodesic_near10_vs_far20"]=fit(matched,"log10_TN","NEAR")
        m["TP_geodesic_near10_vs_far20"]=fit(matched,"log10_TP","NEAR")
        m["HALL_with_same_site_TNTP_conditioning_NONCAUSAL"]=fit(matched,"HALL","NEAR+log10_TN+log10_TP")
    else:
        out["original_measured_2021_wetland_N_and_P_sample_by_externalOutfall_proximity"]["STOP_actual_TNTP_comparison_source_n"]=True
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=study(root,a.source_dir)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":result.get("original_frozen_external_outfall_only_source_sample"),
        "models":result.get("same_original_standardized_wetlands_2021_plant_models"),
        "prior":result.get("original_old_all_feature_same_region_240plot_site_vs_weighted"),
        "5km":result.get("source_fixed_5km_externalOutfall_sensitivity"),
        "actualNP":result.get("original_measured_2021_wetland_N_and_P_sample_by_externalOutfall_proximity")
    },indent=2,ensure_ascii=False))
if __name__=="__main__":
    main()
