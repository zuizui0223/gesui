#!/usr/bin/env python3
"""Austria67 same-2022 observed TP + benthic EQR, 2021 local landscape:
test original HydroRIVERS reach morphology and upstream network complexity
against EEA UWWTD reported N/P/Q source spatial prediction.

Not RiverATLAS full upstream watershed land cover. Not a WWTP causal effect.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/"analysis/20261116_AT67_riparian_satellite"))
sys.path.insert(0,str(ROOT/"analysis/20261116_austria_wisa_official_source"))
sys.path.insert(0,str(HERE))
from WorldCover_AT67_WWTP_report_presence_vs_Pdose import original_and_covariates
from AT68_equal_Q_updown_mass_adversarial_qa import model_cv,comparison
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import P_CONVERSION
from HydroRIVERS_50km_upstream_branch_corridor_source_gate import (
    original_net,upstream,BQE,BQE_BLOB,ORIGINAL_HYDRO_SHA,REGIONS
)

SCHEMA="GESUI_20261117_AT67_NATIVE_HYDRORIVERS_MORPHOLOGY_INDEPENDENT_2022_WWTP_PN_TP_EQR_V1"
REF_WATER_PQ_GAIN=4.58218306901158
REF_EQR_PQ_GAIN=8.388318614326595
BG=("logQ","logArea","lat","lon",
    "wc1km_40_fraction","wc1km_50_fraction",
    "wc1km_10_fraction","wc1km_30_fraction")
MORPH=("source_log_local_catchment","source_strahler_order",
       "source_log_dist_up_headwater","source_log_dist_down_outlet")
NETWORK=("source_log_upstream_50km_branch_n","source_log_upstream_50km_network_total_length")

def blob(p):
    data=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(data))+"\0").encode()+data).hexdigest()

def source_true_reach_station_registry():
    if blob(BQE)!=BQE_BLOB:raise RuntimeError("Original 2022 field BQE source checksum changed")
    with BQE.open(encoding="utf-8-sig",newline="") as f:
        allrows=list(csv.DictReader(f))
    if len(allrows)!=2475:raise RuntimeError("Official EEA original 2022 benthic n2475 drift")
    at={q["site_key"].split("|")[1]:q for q in allrows if q["countryCode"]=="AT"}
    if len(at)!=68:raise RuntimeError("Actual official source Austrian PKINT 68 station identity changed")
    return at

def get_original_morphology(rows,sourcepath):
    g,down,up,length=original_net(sourcepath)
    for k in ("ORD_STRA","CATCH_SKM","DIST_UP_KM","DIST_DN_KM"):
        if k not in g.columns:raise RuntimeError("Actual original HydroRIVERS source missing morphology "+k)
    ids=g["HYRIV_ID"].astype(np.int64).to_numpy()
    unique={int(x):i for i,x in enumerate(ids)}
    if len(unique)!=len(ids):raise RuntimeError("Original river native reach ids not unique")
    actual=source_true_reach_station_registry()
    errors=[];records=[];true_ids=collections.Counter()
    for row in rows:
        station="AT"+row["station"]
        q=actual.get(station)
        if q is None:
            errors.append({"station":station,"reason":"source_original_benthic_station_ID_missing"});continue
        orig_id=int(q["HYRIV_ID"]);true_ids[orig_id]+=1
        ix=unique.get(orig_id)
        if ix is None:
            errors.append({"station":station,"reason":"source_original_2022_HYRIV_ID_not_in_2026_HydroRIVERS"});continue
        src=g.iloc[ix]
        vals={}
        for key in ("ORD_STRA","CATCH_SKM","DIST_UP_KM","DIST_DN_KM"):
            try:v=float(src[key])
            except (TypeError,ValueError):v=math.nan
            vals[key]=v
        if (not all(math.isfinite(v) for v in vals.values()) or
            vals["ORD_STRA"]<1 or vals["ORD_STRA"]>20 or
            min(vals["CATCH_SKM"],vals["DIST_UP_KM"],vals["DIST_DN_KM"])<0):
            errors.append({"station":station,"reason":"invalid_original_native_morphology_values"});continue
        network=upstream(orig_id,up,length)
        if network is None or network.get("source_STOP"):
            errors.append({"station":station,"reason":network.get("source_reason") if network else "source_upstream_morphology_missing"});continue
        row.update({
           "source_log_local_catchment":math.log1p(vals["CATCH_SKM"]),
           "source_strahler_order":vals["ORD_STRA"],
           "source_log_dist_up_headwater":math.log1p(vals["DIST_UP_KM"]),
           "source_log_dist_down_outlet":math.log1p(vals["DIST_DN_KM"]),
           "source_log_upstream_50km_branch_n":math.log1p(network["source_strict_upstream_branching_reach_count"]),
           "source_log_upstream_50km_network_total_length":math.log1p(network["source_network_upstream_km_total_summed_all_branches"]),
           "source_network_upstream_reach_count":network["source_strict_upstream_branching_reach_count"],
           "source_network_upstream_sum_km":network["source_network_upstream_km_total_summed_all_branches"]})
        for name,field in [
            ("Pup_t_yr","upstream_2022_EEA_P_measured_t_per_yr"),
            ("Pdown_t_yr","downstream_2022_EEA_P_t_per_yr_placebo")]:
            if abs(float(q[field])-row[name])>1e-9:
                raise RuntimeError("Same original EEA annual P source field value changed during join")
        Q=math.expm1(row["logQ"])
        if Q<=0:raise RuntimeError("Original source river Q nonpositive")
        Nup=float(q["upstream_2022_EEA_N_measured_t_per_yr"])
        Ndown=float(q["downstream_2022_EEA_N_t_per_yr_placebo"])
        if min(Nup,Ndown)<0 or not math.isfinite(Nup+Ndown):
            raise RuntimeError("Invalid actual source annual N mass")
        row["source_N_upP_Q"]=math.log1p(P_CONVERSION*Nup/Q)
        row["source_N_downP_Q"]=math.log1p(P_CONVERSION*Ndown/Q)
        row["source_P_downP_Q"]=math.log1p(P_CONVERSION*row["Pdown_t_yr"]/Q)
        records.append(row)
    if len(records)!=67 or errors or collections.Counter(z["block"] for z in records)!=REGIONS:
        raise RuntimeError("No source morphology model unless all 67 stations verified; "+str(errors[:5]))
    return records,{"original_reach_IDs_distinct":len(true_ids),
       "source_reaches_with_more_than1_field_station":sum(v>1 for v in true_ids.values()),
       "source_upstream_network_50km_less_than10km_n":sum(x["source_network_upstream_sum_km"]<10 for x in records),
       "source_local_catchment_km2_minmax":[math.expm1(min(x["source_log_local_catchment"] for x in records)),
               math.expm1(max(x["source_log_local_catchment"] for x in records))],
       "source_RiverATLAS_full_upstream_crop_urban_COVARIATES_NOT_ACQUIRED":True}

def model_spec(response,morphology,network):
    base=list(BG)+([] if response=="logTP" else ["logTP"])
    extra=(list(MORPH) if morphology else [])+(list(NETWORK) if network else [])
    B=base+extra
    return {"background":B,"plus_upstream_P_report_per_flow":B+["logPup_flow_proxy"],
       "plus_upstream_N_report_per_flow":B+["source_N_upP_Q"],
       "plus_downstream_P_report_per_flow_PLACEBO":B+["source_P_downP_Q"],
       "plus_downstream_N_report_per_flow_PLACEBO":B+["source_N_downP_Q"],
       "plus_upstream_P_and_N_report_per_flow":B+["logPup_flow_proxy","source_N_upP_Q"]}

def evaluate(rows,y,morph,net):
    models=model_spec(y,morph,net)
    cv=model_cv(rows,y,models)
    comps={}
    for name in models:
        if name=="background":continue
        comps[name]=comparison(cv,"background",name)
    return {"background_features":models["background"],
       "original_2022_5region_ridge_all_loss":cv,
       "original_2022_source_reported_NP_Q_added_predictability":comps}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    cohort67,qa=original_and_covariates()
    if cohort67 is None or len(cohort67)!=67:
        raise RuntimeError("Original same67 source 2021 local cover materialized failed; no hidden field selection")
    raw,integrity=get_original_morphology(cohort67,a.rivers)
    ret={}
    for response in ("logTP","eqr"):
        ret[response]={
           "B_only_worldcover":evaluate(raw,response,False,False),
           "B_plus_original_reach_morphology4":evaluate(raw,response,True,False),
           "B_plus_reach_morphology4_and_network_complexity2":evaluate(raw,response,True,True)}
    first=ret["logTP"]["B_only_worldcover"]["original_2022_source_reported_NP_Q_added_predictability"][
        "plus_upstream_P_report_per_flow"]["source_5block_equal_MSE_added_gain_pct"]
    second=ret["eqr"]["B_only_worldcover"]["original_2022_source_reported_NP_Q_added_predictability"][
        "plus_upstream_P_report_per_flow"]["source_5block_equal_MSE_added_gain_pct"]
    if abs(first-REF_WATER_PQ_GAIN)>1e-8 or abs(second-REF_EQR_PQ_GAIN)>1e-8:
        raise RuntimeError("Original WorldCover/river source numerical baseline not reproduced")
    output={"schema":SCHEMA,
      "status":"PASS_ORIGINAL_AT67_NATIVE_RIVER_STRUCTURE_2022_WWTP_N_P_WATER_BENTHIC_PREDICTION_NONCAUSAL",
      "source_original_N_WWTP_2022_report_P_WWTP_2022_not_same_measured_instream_load":True,
      "source_HydroRIVERS_Europe_v1_orig_SHA256":ORIGINAL_HYDRO_SHA,
      "true_original_2022_water_station_and_2021_landcover_all67":True,
      "source_original_morphology_integrity":integrity,
      "same_original_five_regions":REGIONS,
      "pre_model_Q_EQR_source_identical_baseline_PQ_gains":{"river_lnTP":first,"benthic_EQR_after_real_TP":second},
      "original_2022_TP_and_EQR_morphological_confounds":ret,
      "source_structural_covariates_are_reach_morphology_not_full_upstream_landuse":True,
      "original_RiverATLAS_public_cropped_source_Austria_0_of67_UNCHANGED_STOP":True,
      "original_upstream_10km_landcover_corridor_59of67_below60_UNCHANGED_STOP":True,
      "source_2022_S2_riparian_naturalNDVI_smoke_5of10_below6_STILL_STOP":True,
      "nontemporal_ecology_wastewater_causal_mediation_not_identified":True,
      "no_benthic_biomass_plant_butterfly_ecosystem_service_claim":True}
    (a.output/"AUSTRIA67_HYDROMORPHOLOGY_UPSTREAM_REPORT_NP_Q_TP_BENTHIC_EQR_REGIONCV.json").write_text(
        json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    small={"status":output["status"],"n":len(raw),
           "source67_real_WorldCover2021_QA":qa,
           "reach_integrity":integrity,
           "results":{y:{name:{
               k:round(v["source_5block_equal_MSE_added_gain_pct"],6) for k,v in val[
                 "original_2022_source_reported_NP_Q_added_predictability"].items()}
               for name,val in cases.items()} for y,cases in ret.items()},
           "source_not_true_upstream_landcover_or_causal":True}
    print("GESUI_AUSTRIA67_ORIGINAL_REACH_MORPHOLOGY_WWTP_NP_EQR_STRESS="+json.dumps(
        small,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--rivers",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
