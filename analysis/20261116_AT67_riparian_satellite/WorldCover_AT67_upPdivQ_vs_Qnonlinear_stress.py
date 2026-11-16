#!/usr/bin/env python3
"""Posthoc adversarial exact Austria67 actual 2022 river water/biology source:
test whether WWTP annual reported P / original river Q heldout prediction
is fully explained by NONLINEAR Q alone after original 2021 landcover,
real measured river TP, season and geographic 5-region holdout.

This is not a causal WWTP source or riparian plant biomass estimate.
"""
from __future__ import annotations
import argparse,collections,json,math,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parents[0]/"20261116_austria_wisa_official_source"))
from WorldCover_AT67_WWTP_report_presence_vs_Pdose import original_and_covariates
from AT68_equal_Q_updown_mass_adversarial_qa import model_cv,comparison
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import P_CONVERSION

SCHEMA="GESUI_20261116_AUSTRIA67_WORLD_COVER_Q_NONLINEARITY_VS_UPSTREAM_SOURCE_REPORT_P_Q_V1"
REGIONS={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
ORIGINAL_GAIN_AFTER_WORLD_COVER_AND_REAL_STREAM_TP=8.388318614326595
ORIGINAL_GAIN_AFTER_WORLD_COVER_FOR_RIVER_LOGTP=4.58218306901158
B=["logQ","logArea","lat","lon",
   "wc1km_40_fraction","wc1km_50_fraction",
   "wc1km_10_fraction","wc1km_30_fraction"]
Q_CASES=[
   ("Q2",["source_logQ_squared"]),
   ("Q_inverse",["source_log1p_inverse_Q"]),
   ("Q2_and_inverse",["source_logQ_squared","source_log1p_inverse_Q"]),
   ("Q_unit_annual_P",["source_log_unitmass_P_div_Q"]),
   ("Q_inverse_winsor",["source_log_inverse_Q_winsor_99pct"])
]

def add_source_only_nonlinearity(rows):
    flow=np.array([math.expm1(v["logQ"]) for v in rows],float)
    if not np.isfinite(flow).all() or (flow<=0).any():
        raise RuntimeError("Original source annual water Q nonfinite, zero or negative")
    inv=1./flow
    cap=float(np.quantile(inv,.99))
    if cap<=0:raise RuntimeError("Original Q inverse source pre-model cap invalid")
    for x,q,i in zip(rows,flow,inv):
        x["source_logQ_squared"]=float(x["logQ"]**2)
        x["source_log1p_inverse_Q"]=float(math.log1p(i))
        x["source_log_unitmass_P_div_Q"]=float(math.log1p(P_CONVERSION*i))
        x["source_log_inverse_Q_winsor_99pct"]=float(math.log1p(min(i,cap)))
    return {"source_nonlinear_Q_inverse_winsor_99pct_cap_per_Q_units":cap,
        "source_original_Q_min":float(flow.min()),"source_original_Q_max":float(flow.max())}

def specs(y):
    # Each basis is genuine original independent Q-only source field,
    # matched against SAME added observed 2022 wastewater P/Q.
    base=B+([] if y=="logTP" else ["logTP"])
    models={"baseline":base,
            "baseline_plus_up_reportP_div_Q":base+["logPup_flow_proxy"],
            "baseline_plus_down_reportP_div_Q":base+["logPdown_flow_proxy"],
            "baseline_plus_raw_reportP_mass":base+["logPup"],
            "baseline_plus_source_report_presence":base+["up_reported"]}
    for label,terms in Q_CASES:
        models["flow_"+label]=base+terms
        models["flow_"+label+"_plus_up_P_div_Q"]=base+terms+["logPup_flow_proxy"]
        models["flow_"+label+"_plus_down_P_div_Q"]=base+terms+["logPdown_flow_proxy"]
        models["flow_"+label+"_plus_raw_up_P_mass"]=base+terms+["logPup"]
        models["flow_"+label+"_plus_source_presence"]=base+terms+["up_reported"]
    return models

def evaluate(data,y):
    model=specs(y)
    cv=model_cv(data,y,model)
    cmp={}
    for name in ("baseline_plus_up_reportP_div_Q","baseline_plus_down_reportP_div_Q",
                 "baseline_plus_raw_reportP_mass","baseline_plus_source_report_presence"):
        cmp[name+"_versus_base"]=comparison(cv,"baseline",name)
    for title,terms in Q_CASES:
        base="flow_"+title
        cmp["Qonly_"+title+"_versus_baseline"]=comparison(cv,"baseline",base)
        for name in ("plus_up_P_div_Q","plus_down_P_div_Q","plus_raw_up_P_mass","plus_source_presence"):
            cmp["source_"+name+"_AFTER_Qonly_"+title]=comparison(cv,base,base+"_"+name)
    return {"response":y,"source_5region_out_of_region_metrics":cv,
        "source_Qonly_vs_report_P_Q_matched_comparisons":cmp}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    rows,source=original_and_covariates()
    if rows is None:
        raise RuntimeError("Original WorldCover genuine 67/67 source support no longer passes")
    if len(rows)!=67 or collections.Counter(z["block"] for z in rows)!=REGIONS:
        raise RuntimeError("Original exact AT67 true 2022 aquatic national source geography changed")
    q=add_source_only_nonlinearity(rows)
    water=evaluate(rows,"logTP")
    benthic=evaluate(rows,"eqr")
    wbase=water["source_Qonly_vs_report_P_Q_matched_comparisons"][
          "baseline_plus_up_reportP_div_Q_versus_base"]["source_5block_equal_MSE_added_gain_pct"]
    bbase=benthic["source_Qonly_vs_report_P_Q_matched_comparisons"][
          "baseline_plus_up_reportP_div_Q_versus_base"]["source_5block_equal_MSE_added_gain_pct"]
    if abs(wbase-ORIGINAL_GAIN_AFTER_WORLD_COVER_FOR_RIVER_LOGTP)>1e-8:
        raise RuntimeError("Original source water lnTP P/Q effect not independently reproduced")
    if abs(bbase-ORIGINAL_GAIN_AFTER_WORLD_COVER_AND_REAL_STREAM_TP)>1e-8:
        raise RuntimeError("Original source EQR extra P/Q effect not independently reproduced")
    result={"schema":SCHEMA,"status":"PASS_POSTHOC_FIVE_REGION_Q_NONLINEARITY_WWTP_PQ_ADVERSARIAL",
      "original_same_2022_Austria_field_sites":len(rows),
      "original_source_5_region_sizes":REGIONS,
      "source_worldcover2021_original_67_site_QA":source,
      "source_only_flow_nonlinear_basis":q,
      "same_source_real_2022_stream_TP_prediction":water,
      "same_source_real_2022_benthic_EQR_prediction_adjusting_TP":benthic,
      "original_numeric_WorldCover_2022_water_PQ_gain_reproduced":wbase,
      "original_numeric_WorldCover_2022_benthic_after_TP_PQ_gain_reproduced":bbase,
      "fully_5region_block_disjoint_CV_and_original_Ridge_alpha10":True,
      "original_2022_UP_source_report_P_div_flow_model_not_measured_delivered_WWTP_P":True,
      "source_WorldCover_2021_local_landuse_NOT_upstream_50km_landuse":True,
      "source_2022_Sentinel_NDVI_near_outer_SMOKE_5of10_STILL_STOP":True,
      "posthoc_falsification_same_original_67_sites_not_independent_field_confirmation":True,
      "no_invertebrate_biomass_plant_butterfly_or_caused_wastewater_change_estimate":True}
    (a.output/"AUSTRIA67_SOURCE_WWTP_REPORTED_PQ_VS_Q_NONLINEARITY_AFTER_LOCAL_WORLDCOVER.json").write_text(
      json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    def get(x,k):
        return x["source_Qonly_vs_report_P_Q_matched_comparisons"][k]["source_5block_equal_MSE_added_gain_pct"]
    s={"status":result["status"],
      "n":len(rows),
      "worldcover_source_QA_67":source["source_valid_n"],
      "original_2022_river_TP_PQ_gain_after_landcover":wbase,
      "original_2022_benthic_EQR_PQ_gain_after_landcover_and_realTP":bbase,
      "observed_WQ_Q_only":[{"Qbasis":name,
        "Qalone_gain":get(water,"Qonly_"+name+"_versus_baseline"),
        "reportedP_Q_AFTER_Qbasis":get(water,"source_plus_up_P_div_Q_AFTER_Qonly_"+name),
        "downstream_placebo_AFTER_Qbasis":get(water,"source_plus_down_P_div_Q_AFTER_Qonly_"+name)}
        for name,_ in Q_CASES],
      "observed_BQE_Q_only":[{"Qbasis":name,
        "Qalone_gain":get(benthic,"Qonly_"+name+"_versus_baseline"),
        "reportedP_Q_AFTER_Qbasis":get(benthic,"source_plus_up_P_div_Q_AFTER_Qonly_"+name),
        "downstream_placebo_AFTER_Qbasis":get(benthic,"source_plus_down_P_div_Q_AFTER_Qonly_"+name),
        "Pmass_AFTER_Qbasis":get(benthic,"source_plus_raw_up_P_mass_AFTER_Qonly_"+name)}
        for name,_ in Q_CASES],
      "no_causal_mediation_or_WWTP_terrestrial_effect":True}
    print("GESUI_AUSTRIA67_WC_Q_NONLINEARITY_SOURCE_ADVERSARIAL="+json.dumps(s,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
