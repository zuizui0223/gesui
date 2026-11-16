#!/usr/bin/env python3
"""FROZEN IF-GATE original EU H7 WQ TNTP to land EIV plant & butterfly HNS.

SOURCE-ONLY Part1 station/year eligibility MUST pass before any biological
response is modeled; independent science support requires >=100 H7×period,
>=40 genuine H6 and >=5 original countries. No WWTP dose, causal water→
plant→butterfly mediation, organism abundance or ecosystem service claims.
"""
from __future__ import annotations
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import (
    ORIGINAL_PANEL_SHA,FACETS_SHA,EIV_SHA,
    source_cluster_regression,cv_h6
)

SOURCE_SCHEMA="GESUI_20261116_CHEMISTRY_FIRST_PART1_WQ_EU232_STRICT_SITE_YEAR_TNTP_V1"
EXACT_COHORT={"source_min_H7period":100,"source_min_H6":40,"source_min_country":5}
PGEOG=["W","crop","urban","riverflow","naturalout","area",
       "nPlant","nButterfly","nAnimal",
       "NUTRI_coverage_all_Plantae_speciesKeys",
       "log1p_eiv_scored_taxa","log1p_qualified_WQ_stations"]
B_EXTRA=["primary_record_coverage_raw","primary_host_completeness"]
RESP={"plant":"Plant_N_incidence","butterfly":"primary_incidence_mean"}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def source_df(a):
    state=json.loads((a.chemistry/"CHEMISTRY_FIRST_EU232_PART1_SOURCE_GATE.json").read_text())
    path=a.chemistry/"CHEMISTRY_FIRST_EU232_PART1_H7_PERIOD_AMBIENT_TNTP.csv"
    if state.get("schema")!=SOURCE_SCHEMA:
        raise RuntimeError("Only original Part1 source-authentic chemistry first cohort is allowed")
    if (state.get("official_original_EEA2026_Part1_ZIP_SHA256") !=
        "ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082"):
        raise RuntimeError("Official original EEA2026 chemistry source SHA invalid")
    result={"source_panel_geometry_sha256":state.get("frozen_H7_WQ_station_year_CSV_SHA256"),
            "source_validated_WQ_stationyears":state.get("source_siteyear_trait_quality_counts",{}),
            "original_primary_chemistry_source_support":state.get("primary_original_GBIF_basin_country_matched_qualified"),
            "source_gate_status":state.get("status"),
            "source_gate_pass":bool(state.get("source_support_gate_pass"))}
    if not result["source_gate_pass"]:
        result["status"]="STOP_SOURCE_CHEMISTRY_SUPPORT"
        return result,None
    if sha(path)!=state.get("part1_primary_H7_period_CSV_SHA256"):
        raise RuntimeError("Chemistry-first original Part1 H7-period source SHA mismatch")
    wq=pd.read_csv(path,dtype={"h":str,"h6":str})
    if wq.duplicated(["h","period"]).any():
        raise RuntimeError("Original WQ H7×period duplicated")
    assert (wq.n_source_valid_WQ_stations.astype(float)>=1).all()
    raw=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if sha(raw)!=ORIGINAL_PANEL_SHA:raise RuntimeError("GBIF/source original 454 original SHA changed")
    main=pd.read_csv(raw,dtype={"h":str,"h6":str})
    if (len(main),main.h.nunique(),main.h6.nunique())!=(454,232,90):
        raise RuntimeError("Historical main panel geo cardinality changed")
    ev=json.loads((a.eiv/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (ev.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not ev.get("hard_source_support_pass") or ev.get("source_facets_SHA256")!=FACETS_SHA
        or ev.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA):
        raise RuntimeError("Original strict 4011 EIV trait source QA not met")
    eligible=pd.read_csv(a.eiv/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",
                          dtype={"h":str,"h6":str})
    traits=pd.read_csv(a.eiv/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",dtype={"h":str})
    if len(eligible)!=454 or len(traits)!=451:
        raise RuntimeError("Original GBIF plant EIV score support changed")
    elig=eligible.loc[eligible.full_source_trait_eligible.astype(str)=="True",
        ["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
         "exact_NUTRI_scored_vascular_accepted_species"]]
    if len(elig)!=451:raise RuntimeError("Original EIV eligible source 451 change")
    d=main.merge(elig,on=["h","period"],how="inner",validate="one_to_one")
    d=d.merge(traits[["h","period","Plant_N_incidence"]],on=["h","period"],
              how="inner",validate="one_to_one")
    d=d.merge(wq,on=["h","period","h6","country"],how="inner",validate="one_to_one")
    if d.duplicated(["h","period"]).any() or len(d)>451:
        raise RuntimeError("Original same-H7 source gain or duplication")
    d["log1p_eiv_scored_taxa"]=np.log1p(
        d.exact_NUTRI_scored_vascular_accepted_species.astype(float))
    d["log1p_qualified_WQ_stations"]=np.log1p(
        d.n_source_valid_WQ_stations.astype(float))
    d["log1p_ambient_TN"]=np.log1p(
        d.original_part1_strict_TN_mgN_L_site_median.astype(float))
    d["log1p_ambient_TP"]=np.log1p(
        d.original_part1_strict_TP_mgP_L_site_median.astype(float))
    # For the original source-regression helper CATS, encode COUNTRY only:
    # majority EEA EQR biological assessment method is inapplicable to WQ
    # station-first sourcing and is NOT used.
    d["country_method"]=d["country"].fillna("__UNKNOWN__").astype(str)
    for key in [*PGEOG,*B_EXTRA,*RESP.values(),"log1p_ambient_TN","log1p_ambient_TP"]:
        if key not in d or not np.isfinite(d[key].astype(float).to_numpy()).all():
            raise RuntimeError("Missing original/source-locked water or land numeric predictor "+key)
    if (d.NUTRI_coverage_all_Plantae_speciesKeys<.40).any():
        raise RuntimeError("Numeric accepted-species source eligibility lost")
    if (d.exact_NUTRI_scored_vascular_accepted_species<10).any():
        raise RuntimeError("Original EIV minimum species below 10")
    result["joined_H7period"]=len(d)
    result["joined_native_H7"]=int(d.h.nunique())
    result["joined_true_H6"]=int(d.h6.nunique())
    result["joined_majority_country"]=int(d.country.nunique())
    result["source_chemistry_periods"]=sorted(d.period.unique().tolist())
    if (len(d)<EXACT_COHORT["source_min_H7period"] or
        d.h6.nunique()<EXACT_COHORT["source_min_H6"] or
        d.country.nunique()<EXACT_COHORT["source_min_country"]):
        result["status"]="STOP_BIOLOGICAL_SAMPLE_SUPPORT"
        return result,None
    if d["log1p_ambient_TN"].std(ddof=0)<1e-10 or d["log1p_ambient_TP"].std(ddof=0)<1e-10:
        result["status"]="STOP_NUMERICAL_NUTRIENT_SOURCE_VARIATION"
        return result,None
    result["status"]="PASS_SOURCE_AND_BIOLOGICAL_OPPORTUNITY_GATE"
    return result,d

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    status,d=source_df(a)
    if d is None:
        (a.output/"CHEMISTRY_FIRST_H7_NUTRI_HNS_IF_GATE_RESULTS.json").write_text(
            json.dumps(status,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
        print("GESUI_CHEMISTRY_FIRST_WATER_LAND_MODEL="+json.dumps(status,ensure_ascii=False,sort_keys=True),flush=True)
        return
    status["original_source_plant_butterfly_scored_eiv_plant_species"]=4011
    status["original_source_butterfly_hns_scored_species"]=186
    status["source_chemical_TN_TP_log_correlation"]=float(np.corrcoef(d.log1p_ambient_TN,d.log1p_ambient_TP)[0,1])
    status["categorical_country_plus_period_not_EEA_method"]=True
    status["true_H6_5fold_Ridge_alpha10"]=True
    status["slope_units"]="original endpoint per 1 SD of log1p(original ambient nutrient mg/L)"
    status["regressions"]={}
    for nutrient in ("TN","TP"):
        k="log1p_ambient_"+nutrient
        if d[k].std(ddof=0)<1e-10:raise RuntimeError("Constant ambient nutrient source predictor")
        for name,y,base in [
            ("plant_NUTRI","Plant_N_incidence",PGEOG),
            ("butterfly_HNS_after_plantN","primary_incidence_mean",PGEOG+B_EXTRA+["Plant_N_incidence"]),
            ("butterfly_HNS_unadjusted_plantN","primary_incidence_mean",PGEOG+B_EXTRA)
        ]:
            vals=source_cluster_regression(d,y,base+[k],k)
            test=cv_h6(d,y,base,k)
            status["regressions"][name+"_plus_"+nutrient]={
                "beta":vals["adjusted_beta_per_1SD_focal"],
                "H6_cluster_beta95":vals["H6_cluster_descriptive95"],
                "H6_out_added_MSE_gain_pct":test["H6_equal_added_MSE_gain_pct"],
                "H6_out_fixed_error95":test["fixed_out_of_fold_error_H6_resampling95_pct"],
                "n_H7period":len(d),"n_original_H6":int(d.h6.nunique()),
                "not_wwtp_source_identified":True}
    status["any_positive_conditional_coeff_not_sewage_or_ecological_mediation"]=True
    (a.output/"CHEMISTRY_FIRST_H7_NUTRI_HNS_IF_GATE_RESULTS.json").write_text(
       json.dumps(status,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_CHEMISTRY_FIRST_WATER_LAND_MODEL="+json.dumps({
        "status":status["status"],"H7period":status["joined_H7period"],
        "H6":status["joined_true_H6"],"countries":status["joined_majority_country"],
        "source_TN_TP_log_pearson":status["source_chemical_TN_TP_log_correlation"],
        "models":status["regressions"],"not_causal":True},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--chemistry",type=Path,required=True)
    p.add_argument("--eiv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
