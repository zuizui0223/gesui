#!/usr/bin/env python3
"""2022 German WISE6 real river TP source-vs-sourcezero UPSTREAM/DOWNSTREAM
matched quartets: matched spatial difference-in-differences at 8 source+
8 zero monitoring paths (32 distinct water stations).

Source quartets frozen before opening original WISE6 2022 TP numeric values;
quality/overlapping-month gate committed separately before any deltas.
The result is NON-CAUSAL and POSTHOC to earlier wastewater/photo signals:
outfall source not randomized, tributary/mixing and station seasons differ.
"""
from __future__ import annotations
import argparse,hashlib,itertools,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np

EXPECTED_ORIGINAL_QUARTET_BLOB="8bd67be25968311ca16be4026888daa783fbe69f"
MIN_QUALIFYING=6
MIN_GEO_BLOCKS=4
PSEUDO_MG_P_L=0.001
BOOT_N=9999
BOOT_SEED=20261109

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def blob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def med(d):
    return float(np.median(np.asarray(d,float)))

def annual(site):
    if not site["original_site_numeric_4day_3month_LOQ_gate"]:
        raise ValueError("Read TP value only from original quantified 2022 sites passing predeclared LOQ/date QC")
    x=site["original_2022_annual_siteTP_median_of_daily_medians_mgP_L"]
    if not isinstance(x,(int,float)) or not math.isfinite(x) or x<0:
        raise ValueError("Missing original nonnegative quantified TP median")
    return float(x)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-source-only-quartets",required=True,type=Path)
    ap.add_argument("--original-numeric-QC",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    if blob(a.original_source_only_quartets)!=EXPECTED_ORIGINAL_QUARTET_BLOB:
        raise ValueError("Original non-numeric 8 source/8 registered-zero quartet source selection SHA drift")
    src=json.loads(a.original_source_only_quartets.read_text())
    numbers=json.loads(a.original_numeric_QC.read_text())
    if not src.get("predeclared_8_quartet_sourceonly_quality_GATE") or src["original_WISE6_2022_total_TP_numeric_sample_values_OPENED"]:
        raise ValueError("Original source sample did not precede TP numerical outcomes")
    if numbers.get("schema")!="gesui_DE_2022_WISE6_original_8WWTPsource_vs8zero_matched_river_monitor_quartet_32site_numeric_QC_v01":
        raise ValueError("Wrong original 2022 water station numeric source")
    if numbers["original_source_frozen_8X8_quartets_Git_blob_SHA1"]!=EXPECTED_ORIGINAL_QUARTET_BLOB:
        raise ValueError("Original TP numeric source not linked to pre-numeric selected quartets")
    prior=numbers["source_original_32_monitoring_sites_2022_numeric_TP_QC"]
    if len(prior)!=32:raise ValueError("Wrong 32 original source chosen stations")
    station={q["original_country_scheme_site_SHA256"]:q for q in prior}
    if len(station)!=32:raise ValueError("Duplicate original WISE6 TP station ID")
    complete=numbers["original_source_8matched_quartet_numeric_completion"]
    if len(complete)!=8:raise ValueError("Original 8 physical-source quartets changed")
    gate=bool(numbers["original_2022_quartet_numeric_data_eligibility_GATE"])
    qctotal=0;observed=[]
    if gate:
        for group in complete:
            if not group["all_original_four_station_TP_numeric_source_eligible"]:
                continue
            original=group["original_source_and_control_2station_arms"]
            x={}
            for label in ("source","registered_source_zero"):
                arm=original[label]
                uu=station[arm["source_original_upstream_station_hash"]]
                dd=station[arm["source_original_downstream_station_hash"]]
                if not arm["original_source_2station_pair_numeric_QC_and_2month_overlap"]:
                    raise ValueError("Incomplete original monthly-site pair admitted despite numeric gate")
                if arm["common_original_quantified_months_n"]<2:
                    raise ValueError("Missing >=2 original common qualified months between upstream and downstream")
                up=annual(uu);down=annual(dd)
                x[label]={
                    "upstream_TP_mgP_L":up,
                    "downstream_TP_mgP_L":down,
                    "delta_TP_mgP_L":down-up,
                    "logTP_delta_down_minus_up":math.log(down+PSEUDO_MG_P_L)-math.log(up+PSEUDO_MG_P_L),
                    "registered_strict_source_between_P_t_year":float(arm["source_only_original_intervening_WWTP_P_t_year"]),
                    "source_route_path_km":float(arm["source_full_reach_path_km"]),
                    "common_original_months":arm["common_original_2022_quantified_months"],
                }
            if x["source"]["registered_strict_source_between_P_t_year"]<=0:
                raise ValueError("Original source quartet has no strictly intervening registered P source")
            if x["registered_source_zero"]["registered_strict_source_between_P_t_year"]!=0:
                raise ValueError("Original clean registered-source-zero quartet has reported P source between")
            d=x["source"]["logTP_delta_down_minus_up"]-x["registered_source_zero"]["logTP_delta_down_minus_up"]
            observed.append({
              "original_preoutcome_quartet_index":group["original_source_quartet_index"],
              "source_only_original_country2degree_block":group["original_source_UPSTREAM_country2deg_block"],
              "original_source_only_matching_cost":group["source_only_quartet_geometry_matching_cost"],
              "source_pair":x["source"],
              "registered_source_zero_pair":x["registered_source_zero"],
              "matched_difference_of_logTP_river_gradients_source_minus_control":d,
              "source_P_load_reported_tonnes_not_actual_downstream_P_dose":True
            })
        if len(observed)!=numbers["source_quartets_valid_4station_4days_3months_2commonmonths_LOQ"]:
            raise ValueError("Number of admitted 4station quartets not equal to original quality-only frozen gate")
        if len(observed)<MIN_QUALIFYING or len(set(x["source_only_original_country2degree_block"] for x in observed))<MIN_GEO_BLOCKS:
            raise ValueError("Eligible source-site observed magnitudes mismatch 6quartet /4block gate")
    out={
      "schema":"gesui_2022_DE_WISE6_realTP_8source8sourcezero_matched_same_river_quartet_original32sites_observational_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "source_original_8source8registeredZero_quartet_Gitblob_SHA1":EXPECTED_ORIGINAL_QUARTET_BLOB,
      "source_numeric_32site_2022_WISE6_QC_SHA256":sha(a.original_numeric_QC),
      "source_preopening_quartet_selection_SHA256":sha(a.original_source_only_quartets),
      "POSTHOC_overall_hypothesis_after_original_Lycaenidae_source_positive":True,
      "original_frozen_8_source_vs_8_control_quartets":8,
      "source_numeric_4site_2commonmonth_quartets_qualified":numbers["source_quartets_valid_4station_4days_3months_2commonmonths_LOQ"],
      "original_numeric_gate_permits_matched_gradient_model":gate,
      "original_numeric_quality_gate_requires_ge6_quartets_across_ge4_blocks":True,
      "only_nonnegative_actual_2022_riverTP_mgP_L_original":True,
      "source_unexposed_means_no_QUALIFYING_REGISTERED_WWTP_between_sites_not_sewage_free":True,
      "matched_log_gradient_fixed_additive_constant_mgP_L":PSEUDO_MG_P_L,
      "no_Lycaenidae_photo_population_behavior_measured":True,
      "no_hydrologic_outlet_source_randomization_or_causal_claim":True,
    }
    if not gate:
        out["stage"]="STOP_BEFORE_SOURCE_VS_REGISTERED_ZERO_RIVER_TP_GRADIENT_COMPARISON"
        out["reason"]="Original 32site numerical/LOQ or >=2 overlapping months support lacks >=6 complete quartets across >=4 source2deg blocks"
        out["per_quartet_outcomes_not_computed"]=[]
    else:
        vals=np.asarray([x["matched_difference_of_logTP_river_gradients_source_minus_control"] for x in observed],float)
        n=len(vals)
        signed_means=[]
        for signs in itertools.product((-1,1),repeat=n):
            signed_means.append(float(np.mean(vals*np.asarray(signs,float))))
        actual=float(np.mean(vals))
        exact=sum(abs(v)>=abs(actual)-1e-12 for v in signed_means)/len(signed_means)
        rng=np.random.default_rng(BOOT_SEED)
        draws=np.asarray([float(np.mean(rng.choice(vals,size=n,replace=True))) for _ in range(BOOT_N)])
        out["stage"]="COMPLETED_SOURCE_VS_SOURCEZERO_MATCHED_ORIGINAL_RIVER_GRADIENT_OBSERVATIONAL_RESULT"
        out["original_2022_TP_gradient_per_eligible_quartet"]=observed
        out["matched_source_vs_registered_zero_summary"]={
           "n_source_control_quartets":n,
           "n_disjoint_original_monitoring_stations":4*n,
           "independent_country2degree_source_blocks":len(set(x["source_only_original_country2degree_block"] for x in observed)),
           "mean_difference_of_logTP_gradients":actual,
           "median_difference_of_logTP_gradients":float(np.median(vals)),
           "sign_positive_quartets":int((vals>0).sum()),
           "sign_negative_quartets":int((vals<0).sum()),
           "source_upstream_to_downstream_logTP_gradient_median":
               med([x["source_pair"]["logTP_delta_down_minus_up"] for x in observed]),
           "registered_zero_upstream_to_downstream_logTP_gradient_median":
               med([x["registered_source_zero_pair"]["logTP_delta_down_minus_up"] for x in observed]),
           "descriptive_9999_quartet_resample_mean95pct_interval":
               [float(z) for z in np.quantile(draws,[.025,.975])],
           "exact_two_sided_sign_flip_reference_fraction":
               float(exact),
           "exact_sign_flip_permutation_n":len(signed_means),
           "sign_flip_is_NOT_randomized_WWTP_association_or_causal_p_value":True,
        }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_DE_2022_8SOURCE_VS_8SOURCEZERO_TRUE_TP_SAME_RIVER_MATCHED_GRADIENT",
      json.dumps({"stage":out["stage"],"eligible_quartets":out["source_numeric_4site_2commonmonth_quartets_qualified"],
        "summary":out.get("matched_source_vs_registered_zero_summary"),"no_causal_claim":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
