#!/usr/bin/env python3
"""Adversarial source-report binary presence versus P-mass and P/Q amplitude.

Same 67 source-original 2022 AT sites with true TP, benthic EQR and national
landcover 2021; all source-driven field geographies fixed. Original signed
WorldCover full67 source QA #43888461430 reproduced from identical real public
COG. Separate posthoc reporting-vs-dose test; not causal wastewater claim.
"""
from __future__ import annotations
import argparse,collections,csv,json,math,sys,hashlib
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parents[0]/"20261116_austria_wisa_official_source"))
from WorldCover2021_67_WWTP_PQ_confounded_streamTP_benthic import (
    source_site_landcover,EXCLUDE_SITE,ORIG_REGION,
    MIN_VALID_STATIONS,SCHEMA as WC_SCHEMA
)
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import cohort,P_CONVERSION
from AT68_equal_Q_updown_mass_adversarial_qa import model_cv,comparison

SCHEMA="GESUI_20261116_AUSTRIA67_LANDCOVER_AND_2022_TP_SOURCE_WWTP_REPORT_PRESENCE_VS_P_MASS_Q_V1"
REFERENCE={
    "water_upQ_after_worldcover":4.58218306901158,
    "EQR_upQ_after_worldcover":15.030681825978565,
    "EQR_upQ_after_worldcover_and_real_streamTP":8.388318614326595
}
BACKGROUND=("logQ","logArea","lat","lon")
LANDCLASS=("wc1km_40_fraction","wc1km_50_fraction",
           "wc1km_10_fraction","wc1km_30_fraction")

def original_and_covariates():
    allrows=cohort()
    cases=[dict(r) for r in allrows if r["station"]!=EXCLUDE_SITE]
    if len(cases)!=67 or collections.Counter(z["block"] for z in cases)!=ORIG_REGION:
        raise RuntimeError("Source fixed observed 2022 Austria 67 river stations changed")
    valid=[]
    failures=[]
    for pos,row in enumerate(sorted(cases,key=lambda x:x["station"]),1):
        cover,error=source_site_landcover(row["lat"],row["lon"])
        if cover is None:
            failures.append({"station":"AT"+row["station"],"block":row["block"],"error":error})
            continue
        row.update(cover)
        river_Q=math.expm1(row["logQ"])
        if river_Q<=0:raise RuntimeError("Source original annual receiving river discharge Q invalid")
        row["up_reported"]=float(row["Pup_t_yr"]>0)
        row["down_reported"]=float(row["Pdown_t_yr"]>0)
        row["logPdown_flow_proxy"]=math.log1p(P_CONVERSION*row["Pdown_t_yr"]/river_Q)
        valid.append(row)
        print("GESUI_SOURCE_GENUINE_AT67_ORIGINAL_WORLDCOVER_PKINT_PROGRESS="+str(pos)+"/67",flush=True)
    count=collections.Counter(x["block"] for x in valid)
    if (len(valid)<MIN_VALID_STATIONS or count.keys()!=ORIG_REGION.keys()
        or any(v<5 for v in count.values())):
        return None,{"reason":"fixed_original_WorldCover_support_failed",
            "source_valid_n":len(valid),"source_blocks":dict(count),"failures":failures}
    if len(valid)!=67 or failures:
        # Previously verified complete source support 67/67. A later tile
        # outage/changed source cannot silently change the model population.
        raise RuntimeError("Original source WorldCover 67/67 site frame changed, no 60/67 opportunistic rerun")
    return valid,{"source_valid_n":len(valid),"source_blocks":dict(count),
                  "source_upstream_report_presence_n":sum(x["up_reported"] for x in valid),
                  "source_downstream_report_presence_n":sum(x["down_reported"] for x in valid)}

def design():
    B=list(BACKGROUND)+list(LANDCLASS)
    return B,{
        "M0_landcover":B,
        "M1_up_presence":B+["up_reported"],
        "M2_up_mass":B+["logPup"],
        "M3_up_P_over_Q":B+["logPup_flow_proxy"],
        "M4_up_presence_and_P_over_Q":B+["up_reported","logPup_flow_proxy"],
        "M5_down_presence":B+["down_reported"],
        "M6_down_P_over_Q":B+["logPdown_flow_proxy"],
        "M7_up_mass_and_presence":B+["up_reported","logPup"]},{
        "E0_landcover_plus_real_TP":B+["logTP"],
        "E1_up_presence":B+["logTP","up_reported"],
        "E2_up_mass":B+["logTP","logPup"],
        "E3_up_P_over_Q":B+["logTP","logPup_flow_proxy"],
        "E4_up_presence_and_P_over_Q":B+["logTP","up_reported","logPup_flow_proxy"],
        "E5_down_presence":B+["logTP","down_reported"],
        "E6_down_P_over_Q":B+["logTP","logPdown_flow_proxy"],
        "E7_up_mass_and_presence":B+["logTP","up_reported","logPup"]}

def gain(heldout,base,add):
    return comparison(heldout,base,add)

def all_comp(data,baseline,models):
    full=model_cv(data,baseline,models)
    comparisons={}
    pairs=[("report_presence","M0_landcover","M1_up_presence"),
        ("original_reported_P_mass","M0_landcover","M2_up_mass"),
        ("reported_P_div_river_flow","M0_landcover","M3_up_P_over_Q"),
        ("reported_P_div_flow_AFTER_presence","M1_up_presence","M4_up_presence_and_P_over_Q"),
        ("reported_presence_AFTER_P_div_flow","M3_up_P_over_Q","M4_up_presence_and_P_over_Q"),
        ("reported_P_mass_AFTER_presence","M1_up_presence","M7_up_mass_and_presence"),
        ("downstream_report_presence_placebo","M0_landcover","M5_down_presence"),
        ("downstream_report_P_div_flow_placebo","M0_landcover","M6_down_P_over_Q")]
    return full,comparisons

def compare_fixed(cvdata,baseline,variants,prefix):
    return {name:gain(cvdata,base,expanded)
        for name,base,expanded in variants}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    data,gate=original_and_covariates()
    if data is None:
        result={"schema":SCHEMA,"status":"STOP_2021_WORLD_COVER_SOURCE_RECOVERY_CHANGED",**gate}
        (a.output/"AUSTRIA67_SOURCE_PRESENCE_VS_P_MASS_WITH_WC2021.json").write_text(
            json.dumps(result,indent=2,ensure_ascii=False)+"\n")
        print("GESUI_AT67_SOURCE_PRESENCE_VS_P_MASS_Q="+json.dumps(result,ensure_ascii=False,sort_keys=True),flush=True)
        return
    B,W,E=design()
    water=model_cv(data,"logTP",W)
    eqr=model_cv(data,"eqr",E)
    watervariants=[
        ("up_P_presence_vs_landcover","M0_landcover","M1_up_presence"),
        ("up_P_mass_vs_landcover","M0_landcover","M2_up_mass"),
        ("up_P_div_flow_vs_landcover","M0_landcover","M3_up_P_over_Q"),
        ("up_P_div_flow_AFTER_source_presence","M1_up_presence","M4_up_presence_and_P_over_Q"),
        ("source_presence_AFTER_up_P_div_flow","M3_up_P_over_Q","M4_up_presence_and_P_over_Q"),
        ("up_P_mass_AFTER_source_presence","M1_up_presence","M7_up_mass_and_presence"),
        ("downstream_source_presence_vs_landcover","M0_landcover","M5_down_presence"),
        ("downstream_P_div_flow_vs_landcover","M0_landcover","M6_down_P_over_Q")]
    eqrvariants=[
        ("up_P_presence_AFTER_landcover_and_actual_TP","E0_landcover_plus_real_TP","E1_up_presence"),
        ("up_P_mass_AFTER_landcover_and_actual_TP","E0_landcover_plus_real_TP","E2_up_mass"),
        ("up_P_div_flow_AFTER_landcover_and_actual_TP","E0_landcover_plus_real_TP","E3_up_P_over_Q"),
        ("up_P_div_flow_AFTER_landcover_actual_TP_and_report_presence","E1_up_presence","E4_up_presence_and_P_over_Q"),
        ("report_presence_AFTER_landcover_actual_TP_and_up_P_div_flow","E3_up_P_over_Q","E4_up_presence_and_P_over_Q"),
        ("up_P_mass_AFTER_landcover_actual_TP_and_report_presence","E1_up_presence","E7_up_mass_and_presence"),
        ("downstream_source_presence_AFTER_landcover_and_actual_TP","E0_landcover_plus_real_TP","E5_down_presence"),
        ("downstream_P_div_flow_AFTER_landcover_and_actual_TP","E0_landcover_plus_real_TP","E6_down_P_over_Q")]
    wcomp=compare_fixed(water,"M0_landcover",watervariants,"water")
    ecomp=compare_fixed(eqr,"E0_landcover_plus_real_TP",eqrvariants,"eqr")
    for key,observed in [
      ("water_upQ_after_worldcover",wcomp["up_P_div_flow_vs_landcover"]["source_5block_equal_MSE_added_gain_pct"]),
      ("EQR_upQ_after_worldcover",ecomp["up_P_div_flow_AFTER_landcover_and_actual_TP"]["source_5block_equal_MSE_added_gain_pct"])]:
        expected=REFERENCE["water_upQ_after_worldcover"] if key.startswith("water") else REFERENCE["EQR_upQ_after_worldcover_and_real_streamTP"]
        if abs(observed-expected)>1e-8:
            raise RuntimeError("Original independent 2022 landcover-controlled WWTP-PQ model cannot reproduce exact source gain "+str((key,observed,expected)))
    outfile=a.output/"AUSTRIA67_SOURCE_PRESENCE_VS_P_MASS_WITH_WC2021.json"
    result={"schema":SCHEMA,
       "status":"PASS_AT67_SOURCE_PRESENCE_VS_P_DILUTION_POSTHOC_NONCAUSAL",
       "original_source_2022_true_TP_BQE_2021_WorldCover":gate,
       "source_WWTP_reporting_presence_is_not_true_unexposed_river_control":True,
       "original_2022_source_report_positive_up":int(gate["source_upstream_report_presence_n"]),
       "original_2022_source_report_positive_down":int(gate["source_downstream_report_presence_n"]),
       "original_region_disjoint_Ridge_alpha":10,
       "original_n_source_regions":5,
       "water_TP_source_report_presence_vs_P_mass_in_2021_landcover":wcomp,
       "benthic_EQR_source_report_presence_vs_P_mass_after_actual_TP_and_2021_landcover":ecomp,
       "water_TP_five_region_MSE_per_source_model":water,
       "benthic_EQR_five_region_MSE_per_source_model":eqr,
       "all_results_descriptive_5_region_errors_not_confirmatory":True,
       "source_2022_real_TP_annual_mean_NOT_sewage_fraction":True,
       "source_2021_WorldCover_1km_landcover_NOT_full_50km_upstream_landuse":True,
       "source_2022_local_vegetation_NDVI_smoke5of10_STOP_UNCHANGED":True,
       "no_wastewater_to_riparian_biomass_butterflies_or_ecosystem_service_effect":True}
    outfile.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    small={"n":len(data),"regions":gate["source_blocks"],
      "water_TP":{k:round(v["source_5block_equal_MSE_added_gain_pct"],5) for k,v in wcomp.items()},
      "EQ_after_real_TP":{k:round(v["source_5block_equal_MSE_added_gain_pct"],5) for k,v in ecomp.items()},
      "causal_WWTP_2022_effect_not_identified":True}
    print("GESUI_AT67_SOURCE_PRESENCE_VS_P_MASS_Q="+json.dumps(small,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
