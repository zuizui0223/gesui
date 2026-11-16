#!/usr/bin/env python3
"""Same original EU451 H7×period water–plant NUTRI bridge, posthoc observation.

No source data requery, no wastewater-at-site exposure, no true cross-boundary
flux or causal path. Main question: does actual aquatic benthic EQR covary with
the NEW plant nutrient habitat-affinity score after common geography W,
original EEA/GBIF recording efforts and method/period controls, and does
aquatic EQR improve butterfly HNS prediction beyond plant nutrient affinity?
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis/20261115_eu454_plant_nutri"))
from h6_model_if_gate import (ORIGINAL_PANEL_SHA, FACETS_SHA, EIV_SHA, B_GEOG,
    P_GEOG, cv_h6, source_cluster_regression, design_full, H6_weights)

EQR = "EQR_original_site_equal_basin_period_mean"
P_N = "Plant_N_incidence"
HNS = "primary_incidence_mean"

# Background includes the original resource and hydrological geography, EEA
# and GBIF recording effort; EIV scored fraction/number. No HNS outcome.
EQR_BG = ["crop","urban","riverflow","naturalout","area",
          "nEEA","nAnimal","nPlant","nButterfly",
          "NUTRI_coverage_all_Plantae_speciesKeys",
          "log1p_eiv_scored_taxa","W"]
SEED = 20261116


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_cohort(source):
    gate = json.loads((source/"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (gate.get("schema") != "GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not gate.get("hard_source_support_pass")
        or gate.get("source_facets_SHA256") != FACETS_SHA
        or gate.get("source_EIV",{}).get("original_source_SHA256") != EIV_SHA):
        raise RuntimeError("Source-authentic previous EIV trait gate is not valid")
    original = ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    if sha(original) != ORIGINAL_PANEL_SHA:
        raise RuntimeError("Original fixed source panel hash changed")
    df = pd.read_csv(original, dtype={"h":str,"h6":str})
    if (len(df),df.h.nunique(),df.h6.nunique()) != (454,232,90):
        raise RuntimeError("Original original sample units changed")
    score = pd.read_csv(source/"EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv", dtype={"h":str})
    eligible = pd.read_csv(source/"EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv", dtype={"h":str})
    eligible=eligible.loc[eligible.full_source_trait_eligible.astype(str)=="True",
       ["h","period","NUTRI_coverage_all_Plantae_speciesKeys","exact_NUTRI_scored_vascular_accepted_species"]]
    if len(score)!=451 or len(eligible)!=451:
        raise RuntimeError("Source quality gate no longer selects exactly 451")
    if score.duplicated(["h","period"]).any() or eligible.duplicated(["h","period"]).any():
        raise RuntimeError("Duplicate original source row")
    d=df.merge(eligible,on=["h","period"],how="inner",validate="one_to_one")
    d=d.merge(score[["h","period","Plant_N_incidence","Plant_N_record"]],
              on=["h","period"],how="inner",validate="one_to_one")
    if len(d)!=451 or d.h6.nunique()!=90 or d.h.nunique()!=232:
        raise RuntimeError("Source-grain changed after merge")
    d["log1p_eiv_scored_taxa"] = np.log1p(
        d.exact_NUTRI_scored_vascular_accepted_species.astype(float))
    d["country_method"] = (
        d.country.fillna("__UNKNOWN__").astype(str) + "::" +
        d.method.fillna("__UNKNOWN__").astype(str))
    if not np.isfinite(d[[EQR,P_N,HNS,"W"]+EQR_BG].to_numpy(dtype=float)).all():
        raise RuntimeError("Actual original EQR, plant NUTRI, HNS, background incomplete")
    return d


def residual_H6_correlation(df):
    """Source-adjusted EQR-PlantN residual correlation over TRUE H6 means.

    Same original 451×period, no pseudoindependent pooling. Different
    endpoint-specific control sets, so supplementary descriptive only.
    """
    def resid(response,control):
        X,labs=design_full(df,control)
        w=H6_weights(df)
        m=sm.WLS(df[response].to_numpy(float),X,weights=w).fit()
        return df[response].to_numpy(float)-m.fittedvalues
    a=resid(EQR,EQR_BG)
    b=resid(P_N,P_GEOG+["W"])
    d=pd.DataFrame({"h6":df.h6.to_numpy(),"eqr":a,"plantN":b})
    d=d.groupby("h6",sort=True,as_index=False)[["eqr","plantN"]].mean()
    x=d.eqr.to_numpy(float);y=d.plantN.to_numpy(float)
    observed=float(np.corrcoef(x,y)[0,1])
    if not np.isfinite(observed):
        raise RuntimeError("No source-residual variation")
    rng=np.random.default_rng(SEED)
    n=len(d)
    ij=rng.integers(0,n,size=(2000,n))
    bs=np.asarray([np.corrcoef(x[z],y[z])[0,1] for z in ij],float)
    bs=bs[np.isfinite(bs)]
    return {"n_source_H6":n,"r_H6_mean_residual":observed,
            "fixed_residual_H6_bootstrap_descriptive95":np.quantile(bs,[.025,.975]).tolist(),
            "no_environmental_exposure_or_biological_mediation_proved":True}


def main(a):
    d=source_cohort(a.source)
    eqr=source_cluster_regression(d,EQR,EQR_BG+[P_N],P_N)
    cv=cv_h6(d,EQR,EQR_BG,P_N)
    # Conditional terrestrial test: does aquatic status add anything when
    # original EU plant-N and original geography/effort are already known?
    hns=source_cluster_regression(d,HNS,B_GEOG+[P_N,EQR],EQR)
    hns_cv=cv_h6(d,HNS,B_GEOG+[P_N],EQR)
    # Reverse co-response is an alternate predictive endpoint, not causal
    # reversal; use exact unchanged 451 cohort.
    plant=source_cluster_regression(d,P_N,P_GEOG+["W",EQR],EQR)
    plant_cv=cv_h6(d,P_N,P_GEOG+["W"],EQR)
    residual=residual_H6_correlation(d)
    out={"schema":"GESUI_ORIGINAL_EU451_H7_EQR_PLANTN_HNS_NONCAUSAL_BRIDGE_V1",
         "status":"POSTHOC_SAME_SOURCE_WATERSHED_DIAGNOSTIC",
         "fixed_original_H7_panel_sha":ORIGINAL_PANEL_SHA,
         "eligible_original_H7_periods":len(d),
         "eligible_original_H7":int(d.h.nunique()),
         "eligible_true_H6":int(d.h6.nunique()),
         "matched_periods":sorted(d.period.unique().tolist()),
         "EQR_actual_EEA_source_plantN_habitat_indicator_not_riparian_flux":True,
         "source_controls_EQR":EQR_BG,
         "source_controls_PlantN":P_GEOG+["W"],
         "source_controls_HNS":B_GEOG+[P_N],
         "primary_EQR_outcome_plantN_predictor":{
            "H6_cluster":eqr,"out_of_H6_predictive_added_value":cv},
         "secondary_HNS_outcome_EQR_after_plantN":{
            "H6_cluster":hns,"out_of_H6_predictive_added_value":hns_cv},
         "secondary_plantN_outcome_EQR_added_after_W":{
            "H6_cluster":plant,"out_of_H6_predictive_added_value":plant_cv},
         "source_adjusted_H6_mean_EQR_plantN_residuals":residual,
         "noncausal_no_direct_measurement_of_wastewater_N_reaching_plants":True,
         "exploratory_multiple_tests_no_confirmatory_probability":True}
    a.output.mkdir(parents=True,exist_ok=True)
    p=a.output/"GESUI_EU451_EQR_PLANTN_HNS_OBSERVATIONAL_BRIDGE.json"
    p.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    brief={"status":out["status"],"n":len(d),"H6":d.h6.nunique(),
           "EQR_on_PlantN":{"beta":eqr["adjusted_beta_per_1SD_focal"],
                  "cluster95":eqr["H6_cluster_descriptive95"],
                  "H6_out_gain":cv["H6_equal_added_MSE_gain_pct"],
                  "H6_out95":cv["fixed_out_of_fold_error_H6_resampling95_pct"]},
           "HNS_on_EQR_after_PlantN":{"beta":hns["adjusted_beta_per_1SD_focal"],
                  "cluster95":hns["H6_cluster_descriptive95"],
                  "H6_out_gain":hns_cv["H6_equal_added_MSE_gain_pct"],
                  "H6_out95":hns_cv["fixed_out_of_fold_error_H6_resampling95_pct"]},
           "PlantN_on_EQR_after_W":{"beta":plant["adjusted_beta_per_1SD_focal"],
                  "cluster95":plant["H6_cluster_descriptive95"],
                  "H6_out_gain":plant_cv["H6_equal_added_MSE_gain_pct"],
                  "H6_out95":plant_cv["fixed_out_of_fold_error_H6_resampling95_pct"]},
           "H6_residuals":residual}
    print("GESUI_EU451_EQR_PLANT_N_BRIDGE="+json.dumps(brief,ensure_ascii=False,sort_keys=True),flush=True)


if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
