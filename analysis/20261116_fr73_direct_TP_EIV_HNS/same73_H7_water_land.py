#!/usr/bin/env python3
"""Source-verified 2018/2022 genuine French river TP–benthic EQR–plant EIV N–
butterfly HNS at exactly the same 73 H7 watershed×2 period geographic units.

Previously source-locked physical WQ↔EEA I2M2 identical site river reaches:
263 sites,73 H7,24 genuine H6. GBIF plant/butterfly records are whole-H7,
NOT those physical stations, and do NOT identify local sewage nutrient flux
or biological interaction. All exploratory; no causal SEM or PP selection.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import (
    ORIGINAL_PANEL_SHA,FACETS_SHA,EIV_SHA,
    source_cluster_regression,cv_h6
)

ORIGINAL_FR73_SHA="f2ff4e56f245e1aa8f84f5ac099d68ec1c7e27a786c1067439db388b186dd90b"
FR73_NAME="FR73_HYBAS7_NATIVE_2018_2022_REAL_TP_I2M2_EQR_AND_GBIF_INDEX_DELTA.csv"
SOURCE_H7_NAME="results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
EARLY="2017-2019";LATE="2020-2022"
SOURCE_SCHEMA="GESUI_20261116_FR73_TP_EQR_PLANT_N_HNS_SAME_WATERSHED_V1"
QUALITY_H7_MIN=65
QUALITY_H6_MIN=20
QUALITY_PERIOD_ROWS_MIN=130

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_join(frdir,eivdir):
    file=frdir/FR73_NAME
    if sha(file)!=ORIGINAL_FR73_SHA:
        raise RuntimeError("Genuine original 2018/2022 measured TP/EQR FR73 raw derived archive SHA changed")
    fr=pd.read_csv(file,dtype={"HYBAS_ID_L07":str,"HYBAS_ID_L06":str})
    if len(fr)!=73 or fr.HYBAS_ID_L07.nunique()!=73 or fr.HYBAS_ID_L06.nunique()!=24:
        raise RuntimeError("Original same physical field WQ+EEA source 73/24 topology changed")
    if (fr.original_same_physical_EEA_I2M2_and_WQ_TP_station_n.astype(float)<1).any():
        raise RuntimeError("Original matched station source count invalid")
    raw=ROOT/SOURCE_H7_NAME
    if sha(raw)!=ORIGINAL_PANEL_SHA:
        raise RuntimeError("Source original same-H7/period butterfly GBIF panel SHA mismatch")
    original=pd.read_csv(raw,dtype={"h":str,"h6":str})
    if len(original)!=454 or original.h.nunique()!=232 or original.h6.nunique()!=90:
        raise RuntimeError("Original EU454/232/90 unit cardinality changed")
    gate=json.loads((eivdir/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or gate.get("source_facets_SHA256")!=FACETS_SHA or
        gate.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA or
        not gate.get("hard_source_support_pass") or
        not gate.get("original_EIV_NUTRI_scored_keys_unchanged")):
        raise RuntimeError("Original 4,011 numeric EIV accepted vascular-plant source missing")
    eiv=pd.read_csv(eivdir/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",dtype={"h":str})
    available=pd.read_csv(eivdir/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",dtype={"h":str})
    if len(eiv)!=451 or len(available)!=454:
        raise RuntimeError("Original 451/454 EIV score quality source changed")
    allowed=available.loc[available.full_source_trait_eligible.astype(str)=="True",
             ["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
              "exact_NUTRI_scored_vascular_accepted_species"]]
    if len(allowed)!=451:raise RuntimeError("Source eligibility true rows not 451")
    idx=original.merge(fr[["HYBAS_ID_L07","HYBAS_ID_L06"]].rename(
        columns={"HYBAS_ID_L07":"h","HYBAS_ID_L06":"FR_true_H6"}),on="h",
        how="inner",validate="many_to_one")
    if (len(idx)!=146 or idx.h.nunique()!=73 or idx.h6.nunique()!=24
        or set(idx.period)!={EARLY,LATE} or not (idx.h6==idx.FR_true_H6).all()
        or not (idx.country=="FR").all()):
        raise RuntimeError("Historical source French 73 H7 original two-period H6 country geography mismatch")
    df=idx.merge(allowed,on=["h","period"],how="inner",validate="one_to_one")
    df=df.merge(eiv[["h","period","Plant_N_incidence","Plant_N_record"]],
        on=["h","period"],how="inner",validate="one_to_one")
    both=df.groupby("h").period.nunique()
    hvalid=set(both[both==2].index)
    df=df[df.h.isin(hvalid)].copy()
    df=df.merge(fr,how="left",left_on="h",right_on="HYBAS_ID_L07",validate="many_to_one")
    if df.duplicated(["h","period"]).any() or len(df)!=2*df.h.nunique():
        raise RuntimeError("Original paired French source H7 period not preserved after trait gate")
    status={"source_original_FR73_2018_2022_TP_EQR_H7_csv_sha":ORIGINAL_FR73_SHA,
       "original_full_H7_source_SHA256":ORIGINAL_PANEL_SHA,
       "original_source_H7":73,"original_source_true_H6":24,
       "original_source_real_WQ_EEA_physical_site_pairs":263,
       "source_original_plant_numeric_EIV_species":4011,
       "source_original_butterfly_reference_numeric_species":186,
       "source_trait_complete_two_period_distinct_H7":int(df.h.nunique()),
       "source_trait_complete_true_H6":int(df.h6.nunique()),
       "source_trait_complete_H7period_rows":len(df),
       "precommitted_minimum_paired_distinct_H7":QUALITY_H7_MIN,
       "precommitted_minimum_genuine_H6":QUALITY_H6_MIN,
       "precommitted_minimum_H7period_rows":QUALITY_PERIOD_ROWS_MIN,
       "old_annual_2022_confirmed_upstream_report_wwtp_site_pairs_only":4,
       "old_genuine_upstream_source_positive_H7_only":2,
       "same_water_physical_station_not_same_local_GBIF_observation":True,
       "ambient_stream_TP_not_identified_as_wwtp_P":True,
       "no_original_ecology_recomputed_before_quality_gate":True}
    passgate=(df.h.nunique()>=QUALITY_H7_MIN
        and df.h6.nunique()>=QUALITY_H6_MIN
        and len(df)>=QUALITY_PERIOD_ROWS_MIN)
    status["source_gate_pass"]=bool(passgate)
    if not passgate:
        status["status"]="STOP_ORIGINAL_FR73_TRUE_TWO_PERIOD_EIV_SOURCE_SUPPORT"
        return status,None
    df["log1p_eiv_scored_taxa"]=np.log1p(df.exact_NUTRI_scored_vascular_accepted_species.astype(float))
    df["log1p_FR_observed_2018_2022_TP_benthic_sites"]=np.log1p(
        df.original_same_physical_EEA_I2M2_and_WQ_TP_station_n.astype(float))
    df["source_log_measured_TP"]=np.where(df.period==EARLY,
        df.mean_original_2018_logTP.astype(float),df.mean_original_2022_logTP.astype(float))
    df["actual_French_I2M2_EQR"]=np.where(df.period==EARLY,
        df.mean_original_2018_EQR.astype(float),df.mean_original_2022_EQR.astype(float))
    df["country_method"]="FR::I2M2"
    if ((df.source_log_measured_TP.std(ddof=0)<1e-10)
        or (df.Plant_N_incidence.std(ddof=0)<1e-10)
        or not np.isfinite(df[["source_log_measured_TP","actual_French_I2M2_EQR","Plant_N_incidence","primary_incidence_mean"]].to_numpy(float)).all()):
        raise RuntimeError("Valid original French TP EQ plant HNS sources not finite/variable")
    status["status"]="PASS_ORIGINAL_SOURCE_GEOGRAPHY_TRAIT_READY"
    return status,df

BASE_PLANT=["crop","urban","riverflow","naturalout","area","nEEA","nPlant",
            "NUTRI_coverage_all_Plantae_speciesKeys","log1p_eiv_scored_taxa",
            "W","log1p_FR_observed_2018_2022_TP_benthic_sites"]
BASE_BUTTERFLY=["crop","urban","riverflow","naturalout","area","nEEA",
  "nPlant","nAnimal","nButterfly","primary_record_coverage_raw",
  "primary_host_completeness","NUTRI_coverage_all_Plantae_speciesKeys",
  "log1p_eiv_scored_taxa","W","log1p_FR_observed_2018_2022_TP_benthic_sites"]
BASE_EQR=["crop","urban","riverflow","naturalout","area","nEEA",
          "nPlant","W","log1p_FR_observed_2018_2022_TP_benthic_sites"]

def compute(df,target,base,predictor):
    for key in [target,predictor]+base:
        if key not in df or not np.isfinite(df[key].astype(float)).all():
            raise RuntimeError("Source missing outcome or covariate "+key)
    if df[predictor].std(ddof=0)<1e-10:
        raise RuntimeError("Source predictor constant "+predictor)
    coeff=source_cluster_regression(df,target,base+[predictor],predictor)
    cv=cv_h6(df,target,base,predictor)
    return {"focal":predictor,"response":target,"adjusted_beta_per_source_SD":coeff["adjusted_beta_per_1SD_focal"],
       "true_H6_cluster_descriptive95":coeff["H6_cluster_descriptive95"],
       "H6_out_added_MSE_gain_pct":cv["H6_equal_added_MSE_gain_pct"],
       "H6_out_fixed_error_resampling95":cv["fixed_out_of_fold_error_H6_resampling95_pct"],
       "source_n_H7period":len(df),"source_n_H6":df.h6.nunique(),
       "observational_not_sewage_specific_causal":True}

def paired_deltas(df):
    rows=[]
    for h,grp in df.groupby("h",sort=True):
        if len(grp)!=2 or set(grp.period)!={EARLY,LATE}:raise RuntimeError("Paired source H7 attrition")
        a=grp.loc[grp.period==EARLY].iloc[0];b=grp.loc[grp.period==LATE].iloc[0]
        if a.h6!=b.h6 or abs(a.W-b.W)>1e-10:raise RuntimeError("Same French H7 source changes true H6/static W")
        def delta(k):return float(b[k]-a[k])
        rows.append({"h":h,"h6":a.h6,"country_method":"FR::I2M2","period":"change_2018_2022",
            "W":float(a.W),"Plant_N_incidence":delta("Plant_N_incidence"),
            "d_NUTRI":delta("Plant_N_incidence"),"d_HNS":delta("primary_incidence_mean"),
            "d_TP":delta("source_log_measured_TP"),"d_EQR":delta("actual_French_I2M2_EQR"),
            "d_nPlant":delta("nPlant"),"d_nEEA":delta("nEEA"),
            "d_nButterfly":delta("nButterfly"),
            "d_eivcoverage":delta("NUTRI_coverage_all_Plantae_speciesKeys"),
            "d_scored_taxa":delta("log1p_eiv_scored_taxa"),
            "d_primary_record_coverage_raw":delta("primary_record_coverage_raw"),
            "d_primary_host_completeness":delta("primary_host_completeness")})
    q=pd.DataFrame(rows)
    if len(q)!=df.h.nunique():raise RuntimeError("Source paired observational delta cardinality invalid")
    return q

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    gate,d=source_join(a.fr,a.eiv)
    status={"schema":SOURCE_SCHEMA,**gate,
            "all_original_source_equal_region_weights":True,
            "source_common_H7_period_for_aquatic_plant_butterfly":True,
            "original_WWTP_discharge_not_measured_at_GBIF_patch":True}
    if d is not None:
        q=paired_deltas(d)
        cases=[
         ("same_source_water_EQR_on_real_TP",d,"actual_French_I2M2_EQR",BASE_EQR,"source_log_measured_TP"),
         ("primary_plant_EIV_NUTRI_on_real_TP",d,"Plant_N_incidence",BASE_PLANT,"source_log_measured_TP"),
         ("butterfly_equal_HNS_on_real_TP_without_plant",d,"primary_incidence_mean",BASE_BUTTERFLY,"source_log_measured_TP"),
         ("butterfly_equal_HNS_on_real_TP_after_plant_N",d,"primary_incidence_mean",BASE_BUTTERFLY+["Plant_N_incidence"],"source_log_measured_TP"),
         ("plant_N_on_butterfly_equal_HNS_after_measured_TP",d,"primary_incidence_mean",BASE_BUTTERFLY+["source_log_measured_TP"],"Plant_N_incidence"),
         ("delta_plant_N_on_delta_real_TP",q,"d_NUTRI",["d_nPlant","d_nEEA","d_eivcoverage","d_scored_taxa"],"d_TP"),
         ("delta_butterfly_HNS_on_delta_real_TP_after_delta_plantN",q,"d_HNS",["d_nPlant","d_nEEA","d_nButterfly","d_eivcoverage","d_scored_taxa","d_primary_record_coverage_raw","d_primary_host_completeness","d_NUTRI"],"d_TP"),
         ("delta_butterfly_HNS_on_delta_plantN",q,"d_HNS",["d_nPlant","d_nEEA","d_nButterfly","d_eivcoverage","d_scored_taxa","d_primary_record_coverage_raw","d_primary_host_completeness"],"d_NUTRI")
        ]
        models={}
        for label,data,y,base,x in cases:
            models[label]=compute(data,y,base,x)
        status["models"]=models
        status["n_paired_H7_delta"]=len(q)
        status["true_H6"]=int(d.h6.nunique())
    (a.output/"FR73_MEASURED_TP_EQR_EIV_NUTRI_BUTTERFLY_HNS_ORIGINAL_SAME_GRAIN.json").write_text(
        json.dumps(status,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_FR73_DIRECT_WATER_LAND_TRAIT_TEST="+json.dumps({
       "status":status["status"],"n_scored_H7period":gate["source_trait_complete_H7period_rows"],
       "n_distinct_H7":gate["source_trait_complete_two_period_distinct_H7"],
       "n_true_H6":gate["source_trait_complete_true_H6"],
       "models":status.get("models",{}),
       "still_not_sewage_causal":True},sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--fr",type=Path,required=True)
    p.add_argument("--eiv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
