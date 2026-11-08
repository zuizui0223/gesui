#!/usr/bin/env python3
"""POSTHOC P/N annual reported source-mass source-label permutations and outlier QC.

Uses only SOURCE-LOCKED earlier 2022 DE/AT butterfly *aggregate photo* counts
and original 2021 Plantae-only effort. No 2022 butterfly API/photos acquired.
Preserves initial predeclared outfall-count primary and all negative results.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import gesui_de_at_2022_posthoc_source_block_cv as audit

COUNTRY_DRAWS=399
WITHIN_TWO_DEG_DRAWS=199
SEED=20261108
SOURCE_NAMES={"reported_P_t_per_year":"public_EEA_total_reported_P_t_year",
              "reported_N_t_per_year":"public_EEA_total_reported_N_t_year"}

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda:f.read(1048576),b""):h.update(c)
    return h.hexdigest()

def permute(rows,field,mode,draws,seed,observed):
    rng=np.random.default_rng(seed)
    vals=[]
    reasons=Counter()
    for i in range(draws):
        o=audit.get_override(rows,field,rng,blockwithin=(mode=="within_2degree"))
        try:
            result=audit.cv(rows,field,override=o)
            val=result["heldout_relative_MSE_reduction_percent"]
            if val is not None and math.isfinite(val):vals.append(float(val))
            else:reasons["nonfinite_cv_gain"]+=1
        except (ValueError,np.linalg.LinAlgError) as e:
            reasons[type(e).__name__]+=1
    return {
      "shuffle_definition": "country-preserving P/N source mass randomization" if mode=="within_country" else "preserve country×2deg block while randomly reassigning source mass among cells inside it",
      "fixed_draws_requested":draws,"valid":len(vals),"failed_draws":dict(reasons),
      "seed":seed,"observed_relative_MSE_gain_pct":observed,
      "null_median_pct":float(np.median(vals)) if vals else None,
      "null_95pct_percentile_interval":[float(x) for x in np.quantile(vals,[0.025,0.975])] if len(vals)>=50 else None,
      "fraction_random_label_gains_ge_observed_including_plus1":float((1+sum(v>=observed for v in vals))/(1+len(vals))) if vals else None,
      "not_randomized_WWTP_assignment_or_butterfly_experiment":True,
    }

def drop_top_per_country(rows,field):
    removed=[]
    for c in ("DE","AT"):
        r=[x for x in rows if x["countryCode"]==c]
        if not r:raise ValueError("No country strata for maximum source-mass deletion")
        # Use ORIGINAL SOURCE ONLY field, never identify influential butterfly response cells.
        largest=max(r,key=lambda x:(float(x[field]),float(x["lat_floor"]),float(x["lon_floor"])))
        removed.append(largest)
    selected=[r for r in rows if r not in removed]
    if len(selected)!=len(rows)-2:raise ValueError("Must remove one source-heavy tile in each country")
    res=audit.cv(selected,field,return_records=True)
    res["source_only_removed_tiles"]=[{
        "country":x["countryCode"],
        "original_EEA_source_1degree_lat_floor":x["lat_floor"],
        "original_EEA_source_1degree_lon_floor":x["lon_floor"],
        "exposure_source_mass_t_per_year":float(x[field]),
        "butterfly_outcomes_NOT_used_in_selection":True}
        for x in removed]
    return res

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-result",required=True,type=Path)
    p.add_argument("--plant-gate",required=True,type=Path)
    p.add_argument("--first-block-cv",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    first=json.loads(a.first_block_cv.read_text())
    original=json.loads(a.source_result.read_text())
    plant=json.loads(a.plant_gate.read_text())
    if not first.get("posthoc_NOT_preregistered_before_2022_photo_results"):
        raise ValueError("Previous source model not labeled POSTHOC")
    if first["original_2022_photo_source_result_SHA256"]!=sha(a.source_result):
        raise ValueError("2022 source-original photo data changed after first block validation")
    if first["original_2021_nonbutterfly_plant_effort_gate_SHA256"]!=sha(a.plant_gate):
        raise ValueError("Source-only prior-year plant effort source changed after prior validation")
    if not plant["gate_2021_effort_eligibility_before_2022_butterfly_outcome"] or not original["2022_photo_share_model_gate"]:
        raise ValueError("Original support gate not passed")
    rows=audit.getrows(original)
    if len(rows)!=44 or Counter(r["countryCode"] for r in rows)!=Counter({"DE":29,"AT":15}):
        raise ValueError("Original sample changed; no outcome-based selection allowed")
    tests={}
    for name,field in SOURCE_NAMES.items():
        previous=first["CV_results"][name]
        observed=audit.cv(rows,field,return_records=False)
        gain=observed["heldout_relative_MSE_reduction_percent"]
        if abs(gain-previous["heldout_relative_MSE_reduction_percent"])>1e-10:
            raise ValueError("Source outcome regression replay differs from first block CV")
        tests[name]={
            "same_original_source_panel_out_of_region_base_and_augmented_MSE_replayed":observed,
            "country_permutation_negative_control":permute(rows,field,"within_country",COUNTRY_DRAWS,SEED+10 if name=="reported_P_t_per_year" else SEED+20,gain),
            "within_original_2degree_block_permutation":permute(rows,field,"within_2degree",WITHIN_TWO_DEG_DRAWS,SEED+30 if name=="reported_P_t_per_year" else SEED+40,gain),
            "fixed_remove_highest_2022_source_mass_country_only":drop_top_per_country(rows,field),
        }
        print("GESUI_LYCAENIDAE_2022_POSTHOC_SOURCE_MASS_P_N_STRESS",
          json.dumps({"source":name,"heldout_MSE_improvement_percent":gain,
            "country_shuffle_fraction_ge_observed":tests[name]["country_permutation_negative_control"]["fraction_random_label_gains_ge_observed_including_plus1"],
            "within2deg_shuffle_fraction_ge_observed":tests[name]["within_original_2degree_block_permutation"]["fraction_random_label_gains_ge_observed_including_plus1"],
            "removal_two_biggest_mass_tiles_gain_pct":tests[name]["fixed_remove_highest_2022_source_mass_country_only"]["heldout_relative_MSE_reduction_percent"]},
          sort_keys=True),flush=True)
    out={
      "schema":"gesui_2022_Lycaenidae_posthoc_P_and_N_reported_source_mass_shuffle_influence_v01",
      "date_utc":datetime.now(timezone.utc).isoformat(),
      "POSTHOC_after_positive_P_N_CV_being_viewed":True,
      "first_per_source_2022_geoblock_CV_result_SHA256":sha(a.first_block_cv),
      "source_2022_family_photo_original_result_SHA256":sha(a.source_result),
      "source_2021_nonfocal_plant_effort_original_gate_SHA256":sha(a.plant_gate),
      "source_P_and_N_reported_annual_mass_collinearity_rho":original["exploratory_wastewater_source_vs_public_photo_reporting_descriptive"]["separate_P_N_2022_annual_WWTP_reported_mass_log1p_rank_correlation"],
      "original_2022_source_cells":len(rows),
      "original_source_by_country":{"DE":29,"AT":15},
      "source_mass_P_and_N_original_heldout_tests_and_falsification":tests,
      "source_count_predeclared_primary_remains_separate":True,
      "no_butterfly_observations_or_photos_added":True,
      "not_identifying_chemical_P_versus_N_ecological_effect":True,
      "not_identifying_wastewater_causing_butterfly_abundance_or_puddling":True,
      "remaining_confounds":[
          "source intensity correlates with 2021 plant observation effort especially in Austria",
          "one degree cross border land fraction and built/cropland land cover not independently controlled",
          "species-level host plants and butterfly photography observer effort not modeled",
          "WWTP reported source mass is not delivered freshwater/soil nutrient dose",
          "nonrandom human settlement/WWTP reporting intensity versus insect photo coverage",
      ],
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_LYCAENIDAE_P_N_WWTP_MASS_POSTHOC_NULL_AND_INFLUENCE_VERIFICATION_SUCCESS",flush=True)

if __name__=="__main__":main()
