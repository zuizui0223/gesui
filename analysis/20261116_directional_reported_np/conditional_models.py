#!/usr/bin/env python3
"""Source-audited 2018/2022 reported WWTP directional placebo on plant N/HNS.

Follow an outcome-blind 454 true-H7 source join, 451 EIV-scored H7×periods,
56 original upstream-positive and 61 downstream-positive annual reports.
These are *reported mass/flow proxies*, NOT receiving-water chemistry or
assigned exposed versus untreated field sites. Do NOT claim mediation.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np, pandas as pd

from sys import path as sys_path
sys_path.insert(0,str(Path("analysis/20261115_eu454_plant_nutri").resolve()))
from h6_model_if_gate import (
    ORIGINAL_PANEL_SHA,P_GEOG,B_GEOG,source_cluster_regression,
    cv_h6,design_full,H6_weights
)

DIRECTION_BLOB_SHA="702131afc7d1a7bfb15e80357a506754cae7a70a"
SCHEMA="GESUI_20261116_REPORT_POSITIVE_UPDOWN_EIV_PLANT_BUTTERFLY_SOURCE_GATE_V1"
DIR_FILE=Path("results/20261115_np_direction/H7_454_DIRECTIONAL_SOURCE_AUDIT_PANEL.csv")
PANEL=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
OUTCOMES={
   "plant_species_equal_NUTRI":("Plant_N_incidence",P_GEOG+["W"]),
   "plant_record_weighted_NUTRI":("Plant_N_record",P_GEOG+["W"]),
   "butterfly_species_equal_HNS_adjusted_for_plantN":(
       "primary_incidence_mean",B_GEOG+["Plant_N_incidence"]),
   "butterfly_species_equal_HNS_without_plantN":(
       "primary_incidence_mean",B_GEOG),
   "butterfly_record_weighted_HNS_adjusted_for_plantN":(
       "primary_record_mean",B_GEOG+["Plant_N_incidence"])
}
TYPES={
    "report_indicator":("upstream_report_positive","downstream_report_positive"),
    "reported_N_full_delivery_mgL_log1p":("log1p_upstream_N_equiv","log1p_downstream_N_equiv"),
    "reported_P_full_delivery_mgL_log1p":("log1p_upstream_P_equiv","log1p_downstream_P_equiv")
}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()

def load(a):
    gate=json.loads((a.gate/"UPDOWN_EIV_PLANT_BUTTERFLY_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!=SCHEMA or not gate.get("source_gate_support_pass")
        or gate.get("original_source_rows")!=454
        or gate.get("eligible_source_H7periods")!=451
        or gate.get("eligible_upstream_positive_H6")!=31
        or gate.get("original_H7panel_SHA256")!=ORIGINAL_PANEL_SHA
        or gate.get("original_direction_panel_git_blob")!=DIRECTION_BLOB_SHA):
        raise RuntimeError("No original-outcome-blind 451 site source gate success")
    if sha(PANEL)!=ORIGINAL_PANEL_SHA or blob(DIR_FILE)!=DIRECTION_BLOB_SHA:
        raise RuntimeError("Original two-source H7 panels immutable identifiers changed")
    original=pd.read_csv(PANEL,dtype={"h":str,"h6":str})
    original_direction=pd.read_csv(DIR_FILE,dtype={"h":str,"h6":str})
    if len(original)!=454 or len(original_direction)!=454:
        raise RuntimeError("Genuine original full 454 H7 period sources absent")
    if original.duplicated(["h","period"]).any() or original_direction.duplicated(["h","period"]).any():
        raise RuntimeError("Original H7×period duplicate key")
    for x in ("h6","country","W","primary_record_mean"):
        other=original_direction[["h","period",x]].rename(columns={x:"direction_original_"+x})
        original=original.merge(other,how="left",on=["h","period"],validate="one_to_one")
    if not (original.h6==original.direction_original_h6).all() or not (original.country==original.direction_original_country).all():
        raise RuntimeError("Original reported source H6 or source country identity differs")
    if max(abs(original.W-original.direction_original_W))>1e-10 or max(abs(original.primary_record_mean-original.direction_original_primary_record_mean))>1e-10:
        raise RuntimeError("Original report linkage W/HNS identity differs")
    report_cols=["upstream_report_positive","downstream_report_positive",
       "upstream_N_equiv_mgL","downstream_N_equiv_mgL",
       "upstream_P_equiv_mgL","downstream_P_equiv_mgL"]
    main=original.drop(columns=["direction_original_"+x for x in ("h6","country","W","primary_record_mean")])
    main=main.merge(original_direction[["h","period"]+report_cols],
                    on=["h","period"],validate="one_to_one",how="left")
    source_gate_path=a.eiv/"EU454_FULL_EIV_SOURCE_GATE.json"
    qa=json.loads(source_gate_path.read_text())
    if (qa.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not qa.get("hard_source_support_pass")
        or not qa.get("original_EIV_NUTRI_scored_keys_unchanged")
        or qa.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175
        or qa.get("source_hybrid_score_assigned")!=0):
        raise RuntimeError("Original 4,011 exact EIV source trait corpus not verified")
    elig=pd.read_csv(a.eiv/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv",
                     dtype={"h":str,"h6":str})
    trait=pd.read_csv(a.eiv/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv",
                      dtype={"h":str})
    if len(elig)!=454 or len(trait)!=451 or elig.duplicated(["h","period"]).any() or trait.duplicated(["h","period"]).any():
        raise RuntimeError("EIV prior eligibility/source trait row grain lost")
    keep=elig.loc[elig.full_source_trait_eligible.astype(str)=="True",
         ["h","period","h6","NUTRI_coverage_all_Plantae_speciesKeys",
          "exact_NUTRI_scored_vascular_accepted_species"]]
    d=main.merge(keep,on=["h","period","h6"],how="inner",validate="one_to_one")
    d=d.merge(trait[["h","period","Plant_N_incidence","Plant_N_record","n_scored_species"]],
              on=["h","period"],how="inner",validate="one_to_one")
    if len(d)!=451 or d.h6.nunique()!=90 or d.h.nunique()!=232:
        raise RuntimeError("EIV same-source 451 native H7 periods lost")
    if (d.n_scored_species<10).any() or (d.NUTRI_coverage_all_Plantae_speciesKeys<.4).any():
        raise RuntimeError("Original scored source H7 quality gate bypassed")
    d["log1p_eiv_scored_taxa"]=np.log1p(d.exact_NUTRI_scored_vascular_accepted_species.astype(float))
    d["country_method"]=d.country.fillna("__UNKNOWN__").astype(str)+"::"+d.method.fillna("__UNKNOWN__").astype(str)
    for kind in ("N","P"):
        for direction in ("upstream","downstream"):
            k=f"{direction}_{kind}_equiv_mgL"
            v=d[k].to_numpy(float)
            if not np.isfinite(v).all() or (v<0).any():raise RuntimeError("Original source mass negative/nonfinite")
            d[f"log1p_{direction}_{kind}_equiv"]=np.log1p(v)
    n_up=int(d.upstream_report_positive.sum())
    n_down=int(d.downstream_report_positive.sum())
    if n_up!=56 or n_down!=61:
        raise RuntimeError("Frozen EEA annual WWTP report support n56/n61 changed")
    return d

def compare(d,target,base,addition):
    # Equal true-H6 assessment, fixed 5fold original HUC parents, train-only
    # scaling, H6 2000 fixed-error resampling. Not a treatment-effect CI.
    try:
        full,columns=design_full(d,base+[addition])
        designcond=float(np.linalg.cond(full*np.sqrt(H6_weights(d))[:,None]))
        fit=source_cluster_regression(d,target,base+[addition],addition)
        cv=cv_h6(d,target,base,addition)
        return {"status":"FIT",
           "nominal_H6_WLS_beta":fit["adjusted_beta_per_1SD_focal"],
           "H6_cluster95":fit["H6_cluster_descriptive95"],
           "nominal_p":fit["nominal_p_descriptive"],
           "added_H6_out_of_sample_MSE_gain_percent":cv["H6_equal_added_MSE_gain_pct"],
           "fixed_heldout_errors_H6_bootstrap95":cv["fixed_out_of_fold_error_H6_resampling95_pct"],
           "n_H7period":len(d),"n_true_H6":d.h6.nunique(),
           "n_design_columns":len(columns),"source_design_condition_number":designcond,
           "condition_gt100_warning":designcond>100,
           "not_causal_or_truly_measured_receiving_water":True}
    except (ValueError,RuntimeError,np.linalg.LinAlgError) as ex:
        return {"status":"STOP_MODEL_MATRIX_OR_SUPPORT","reason":type(ex).__name__+": "+str(ex)[:230],
                "no_arbitrary_subgroup_or_covariate_replacement":True}

def main(a):
    d=load(a)
    original_plant=source_cluster_regression(d,"Plant_N_incidence",P_GEOG+["W"],"W")
    original_butterfly=source_cluster_regression(d,"primary_incidence_mean",
                                                  B_GEOG+["Plant_N_incidence"],"Plant_N_incidence")
    if abs(original_plant["adjusted_beta_per_1SD_focal"]-0.21524238037690294)>1e-9:
        raise RuntimeError("Original EU451 W→NUTRI positive source effect fails exact reproduction")
    if abs(original_butterfly["adjusted_beta_per_1SD_focal"]-0.19555118559308948)>1e-9:
        raise RuntimeError("Original EU451 NUTRI→butterfly source effect fails exact reproduction")
    out={"schema":"GESUI_20261116_DIRECTIONAL_REPORTED_NP_NUTRI_BUTTERFLY_HNS_NESTED_H6_V1",
         "source_original_h7_sha256":ORIGINAL_PANEL_SHA,
         "source_original_direction_git_blob":DIRECTION_BLOB_SHA,
         "n_exact_H7period":len(d),"n_true_H6":d.h6.nunique(),"n_H7":d.h.nunique(),
         "n_upstream_positive_report":int(d.upstream_report_positive.sum()),
         "n_downstream_positive_report":int(d.downstream_report_positive.sum()),
         "original_source_baselines_unchanged":{
            "W_to_plant_incidence_NUTRI_beta":original_plant["adjusted_beta_per_1SD_focal"],
            "plantN_to_butterfly_incidence_HNS_beta":original_butterfly["adjusted_beta_per_1SD_focal"]},
         "original_EIV_4011_scored_nutrient_affinity_NOT_leaf_or_river_N":True,
         "NOT_TRUE_UNTREATED_CONTROL_OR_ACTUAL_DELIVERED_EFFLUENT":True,
         "NOT_A_CAUSAL_MEDIATION_OR_ECOSYSTEM_SERVICE_TEST":True,
         "original_sources_and_models":{}}
    summary={}
    for label,(target,bg) in OUTCOMES.items():
        tmp={}
        main_summary={}
        for key,(up,down) in TYPES.items():
            # Reporting opportunities always remain in nutrient intensity
            # baseline; a source-unreported cell is never considered clean.
            common=bg if key=="report_indicator" else bg+[
                "upstream_report_positive","downstream_report_positive"]
            after_down=common+[down]
            after_up=common+[up]
            # Avoid a duplicate source indicator when the exposure itself IS
            # the reporting indicator; both opposite indicators are included
            # exactly once.
            if key=="report_indicator":
                after_down=bg+[down]
                after_up=bg+[up]
            u=compare(d,target,after_down,up)
            z=compare(d,target,after_up,down)
            tmp[key]={"upstream_added_after_downstream":u,
                      "downstream_added_after_upstream":z,
                      "same_451_EIV_source_rows":True}
            main_summary[key]={"up_gain":u.get("added_H6_out_of_sample_MSE_gain_percent"),
                               "up_gain95":u.get("fixed_heldout_errors_H6_bootstrap95"),
                               "up_beta":u.get("nominal_H6_WLS_beta"),
                               "up_beta95":u.get("H6_cluster95"),
                               "down_gain":z.get("added_H6_out_of_sample_MSE_gain_percent"),
                               "down_gain95":z.get("fixed_heldout_errors_H6_bootstrap95"),
                               "down_beta":z.get("nominal_H6_WLS_beta"),
                               "down_beta95":z.get("H6_cluster95"),
                               "up_status":u["status"],"down_status":z["status"]}
        out["original_sources_and_models"][label]=tmp
        summary[label]=main_summary
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EU451_UPDOWN_REPORTED_NP_PLANT_BUTTERFLY_NESTED_MODELS.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    compact={"n_original_H7period":len(d),"n_true_H6":d.h6.nunique(),
        "original_W_to_plant_N_beta_reproduced":out["original_source_baselines_unchanged"]["W_to_plant_incidence_NUTRI_beta"],
        "original_plantN_to_HNS_beta_reproduced":out["original_source_baselines_unchanged"]["plantN_to_butterfly_incidence_HNS_beta"],
        "outcomes":summary,
        "reported_source_not_actual_river_NP_or_clean_controls":True}
    print("GESUI_EU451_REPORTED_UPDOWN_NP_PLANT_HNS_RESULT="+json.dumps(
        compact,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--eiv",type=Path,required=True)
    p.add_argument("--gate",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
