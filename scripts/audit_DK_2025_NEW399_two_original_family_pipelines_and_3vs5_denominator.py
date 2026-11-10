#!/usr/bin/env python3
"""Independent offline cross-replay of two 2025 new399 licensed eBMS family counts.

Source-locked 399 site hashes, original 2025 survey effort, 6 accepted butterfly
families, 2022 registered EEA WWTP source IDs are verified row by row.
Separately verify the opposite signs of 3-family vs 5-family pooled fraction,
the 2025 Lycaenidae count PER documented date, and six within-1° area contrasts.
Not an intervention, causal wastewater inference, or preregistered test.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import defaultdict, Counter
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
PRIMARY=ROOT+"DENMARK_2025_NEW399_SOURCE_LOCKED_ACCEPTED_BUTTERFLY_FAMILY_COUNTS.json"
INDEPENDENT=ROOT+"DK_EBMS_2025_NEW399_ORIGINAL_FAMILY_COUNT_BEFORE_DISTANCE_CONTRAST.json"
GEO=ROOT+"DENMARK_2025_NEW399_LYCAENIDAE_2022_WWTP_AGGREGATE_REGION_COMPONENT_SENSITIVITY.json"
FIVE=ROOT+"DENMARK_2025_NEW399_THREE_VS_FIVE_FAMILY_SIGN_REVERSAL_INDEPENDENT_QC.json"
FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae","Riodinidae")
THREE=("Lycaenidae","Pieridae","Papilionidae")
FIVE_FAMILIES=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def odds(counts, selected):
    ly=counts["Lycaenidae"]
    others=sum(counts[f] for f in selected if f!="Lycaenidae")
    return math.log((ly+.5)/(others+.5)) if ly+others>0 else None

def tally(rows):
    counts={fam:sum(r["original_2025_family_individualCount"][fam] for r in rows) for fam in FAMILIES}
    days=sum(r["original_2025_distinct_dates_with_any_published_Occurrence"] for r in rows)
    return {"site_count":len(rows),"any_occurrence_linked_dates":days,"family_counts":counts,
            "Lycaenidae_per_ANY_Occurrence_linked_date":counts["Lycaenidae"]/days if days else None,
            "Lycaenidae_share_3family":counts["Lycaenidae"]/sum(counts[f] for f in THREE),
            "Lycaenidae_share_5family":counts["Lycaenidae"]/sum(counts[f] for f in FIVE_FAMILIES),
            "logodds_3family":odds(counts,THREE),"logodds_5family":odds(counts,FIVE_FAMILIES)}

def verify(root):
    files={k:root/p for k,p in (("prior_frozen_399_family",PRIMARY),("new_independent_original_399_family",INDEPENDENT),
                              ("previous_region_component_analysis",GEO),("previous_5family_change",FIVE))}
    records={k:json.loads(path.read_text(encoding="utf-8")) for k,path in files.items()}
    original=records["prior_frozen_399_family"]
    newer=records["new_independent_original_399_family"]
    geo=records["previous_region_component_analysis"]
    denom=records["previous_5family_change"]
    x=original["all399_source_site_original_family_counts_private_hash_aggregate_only"]
    y=newer["original_2025_NEW399_individualCount_source_site_records"]
    assert len(x)==len(y)==399
    aa={r["site_id_SHA256"]:r for r in x}
    bb={r["site_SHA256"]:r for r in y}
    assert len(aa)==len(bb)==399 and set(aa)==set(bb)
    mismatches=[]
    for key in sorted(aa):
        a,b=aa[key],bb[key]
        checks={
          "source_band":(a["source_band"],b["source_original_2022_registered_outfall_band"]),
          "original_1degree_region":(a["original_1degree_region"],b["original_country_1deg_region"]),
          "source_outfall_id":(a["2022_nearest_registered_EEA_FW_outfall_SHA256"],b["nearest_original_2022_WWTP_outfall_SHA256"]),
          "registered_outfall_distance_km":(a["2022_nearest_qualifying_WWTP_source_distance_km"],b["original_2025_source_distance_km"]),
          "source_5km_component":(a["original_WWTP_spatial_component_5km_SHA256"],b["near_outfall_5km_spatial_cluster_SHA256"]),
          "source_10km_component":(a["original_WWTP_spatial_component_10km_SHA256"],b["near_outfall_10km_spatial_cluster_SHA256"]),
          "all_Event_dates":(a["original_all_Event_distinct_date_count"],b["original_2025_distinct_Event_dates"]),
          "any_Occurrence_linked_dates":(a["original_any_Occurrence_linked_distinct_date_count"],b["original_2025_distinct_dates_with_any_published_Occurrence"]),
          "preoutcome_gate":(a["pre_2025_family_any_occurrence_linked_site_source_gate_pass"],b["original_2025_2plus_any_Occurrence_date_source_gate"])
        }
        for fam in FAMILIES:checks[fam]=(a["family_individual_count_original"][fam],b["original_2025_family_individualCount"][fam])
        for name,(m,n) in checks.items():
            if m!=n:mismatches.append({"sitehash":key,"field":name})
    if mismatches:raise ValueError("Independent source row mismatch: "+str(mismatches[:10]))
    assert newer["original_source_SHA256"]["frozen_preoutcome_NEW399_site_source"]==original["original_source_399_2025_event_any_occurrence_gate_SHA256"]
    qualified=[r for r in y if r["original_2025_2plus_any_Occurrence_date_source_gate"]]
    assert len(qualified)==286
    bygroup={k:[r for r in qualified if r["source_original_2022_registered_outfall_band"]==k]
             for k in ("near5","remote10","intermediate")}
    assert [len(bygroup[k]) for k in ("near5","remote10","intermediate")]==[70,92,124]
    near=tally(bygroup["near5"])
    far=tally(bygroup["remote10"])
    compare={
        "original_new399_cohort_near":near,"original_new399_cohort_remote":far,
        "Lycaenidae_recording_per_published_date_near_over_remote":
           near["Lycaenidae_per_ANY_Occurrence_linked_date"]/far["Lycaenidae_per_ANY_Occurrence_linked_date"],
        "pooled_three_family_share_oddsratio_near_over_remote":
           math.exp(near["logodds_3family"]-far["logodds_3family"]),
        "pooled_five_family_share_oddsratio_near_over_remote":
           math.exp(near["logodds_5family"]-far["logodds_5family"])
    }
    allregions=defaultdict(lambda:defaultdict(list))
    for r in qualified:
        allregions[r["original_country_1deg_region"]][r["source_original_2022_registered_outfall_band"]].append(r)
    matched=sorted(k for k,x in allregions.items() if x["near5"] and x["remote10"])
    assert len(matched)==6
    perregion=[]
    for key in matched:
        near=tally(allregions[key]["near5"]);far=tally(allregions[key]["remote10"])
        perregion.append({
            "region_1degree":key,
            "source_site_n":near["site_count"],"remote_site_n":far["site_count"],
            "Lycaenidae_3family_composition_logodds_near_minus_remote":near["logodds_3family"]-far["logodds_3family"],
            "Lycaenidae_5family_composition_logodds_near_minus_remote":near["logodds_5family"]-far["logodds_5family"],
            "Lycaenidae_recording_rate_near_over_remote":
                 near["Lycaenidae_per_ANY_Occurrence_linked_date"]/far["Lycaenidae_per_ANY_Occurrence_linked_date"] if far["Lycaenidae_per_ANY_Occurrence_linked_date"] else None,
            "near_3family_share":near["Lycaenidae_share_3family"],"remote_3family_share":far["Lycaenidae_share_3family"],
            "near_5family_share":near["Lycaenidae_share_5family"],"remote_5family_share":far["Lycaenidae_share_5family"],
        })
    equal3=sum(r["Lycaenidae_3family_composition_logodds_near_minus_remote"] for r in perregion)/6
    equal5=sum(r["Lycaenidae_5family_composition_logodds_near_minus_remote"] for r in perregion)/6
    assert abs(compare["Lycaenidae_recording_per_published_date_near_over_remote"]-geo["pooled_new399_2025_butterfly_family_count_vs_daily_recorded_rate"]["Lycaenidae_per_linked_date_rate_ratio"])<1e-12
    assert abs(compare["pooled_three_family_share_oddsratio_near_over_remote"]-denom["three_family_vs_five_family_denominator_response"]["three_family"]["pooled_Lycaenidae_share_OR_near_divided_remote"])<1e-12
    assert abs(compare["pooled_five_family_share_oddsratio_near_over_remote"]-denom["three_family_vs_five_family_denominator_response"]["five_family"]["pooled_Lycaenidae_share_OR_near_divided_remote"])<1e-12
    assert abs(equal3-geo["equal_1degree_region_unweighted_mean_pooled_family_log_odds_near_minus_remote"])<1e-12
    assert abs(equal5-denom["three_family_vs_five_family_denominator_response"]["five_family"]["equal_1degree_region_mean_log_OR"])<1e-12
    return {
        "schema":"gesui_DK2025_NEW399_2independent_raw_source_aggregate_exact_replay_and_denominator_sign_v01",
        "archived_original_input_sha256":{k:sha(path) for k,path in files.items()},
        "original_new399_2025_raw_licensed_Event_Occurrence_tables_NOT_accessed_for_this_offline_audit":True,
        "original_previous_107_outcome_sites_excluded":True,
        "original_source_399_independently_reconstructed_site_keys":len(aa),
        "six_accepted_family_exact_site_count_and_source_QA_mismatches":len(mismatches),
        "source_frozen_before_2025_family_counts":True,
        "preselected_original_reporting_2plus_any_occurrence_dates_n":len(qualified),
        "independent_NEW399_near_remote_denominator_effects":compare,
        "region_equal_weight_families_logodds": {
            "n_country_1degree_near_and_remote_regions":len(matched),
            "3family_equal_region_mean_logOR":equal3,
            "3family_equal_region_or":math.exp(equal3),
            "5family_equal_region_mean_logOR":equal5,
            "5family_equal_region_or":math.exp(equal5),
            "regions_positive_3family":sum(r["Lycaenidae_3family_composition_logodds_near_minus_remote"]>0 for r in perregion),
            "regions_positive_5family":sum(r["Lycaenidae_5family_composition_logodds_near_minus_remote"]>0 for r in perregion),
            "original_six_regions":perregion,
        },
        "do_not_infer_2025_sewage_positive_or_negative_population_effect":True,
        "five_family_vs_three_family_is_an_exploratory_posthoc_denominator_stress_test":True,
        "registered_source_WWTP_5km_proximity_not_actual_hydrological_exposure":True,
        "full_Denmark_original_eBMS_license":"CC-BY-NC 4.0",
    }

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out",required=True,type=Path)
    x=a.parse_args()
    r=verify(Path(__file__).resolve().parents[1])
    x.out.parent.mkdir(parents=True,exist_ok=True)
    x.out.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "row_mismatches":r["six_accepted_family_exact_site_count_and_source_QA_mismatches"],
        "cohort":r["independent_NEW399_near_remote_denominator_effects"],
        "regions":{k:v for k,v in r["region_equal_weight_families_logodds"].items() if k!="original_six_regions"}
    },indent=2))
