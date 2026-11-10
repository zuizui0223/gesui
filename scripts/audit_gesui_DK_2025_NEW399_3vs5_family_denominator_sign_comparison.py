#!/usr/bin/env python3
"""Denmark eBMS 2025 NEW399: fixed outcome denominator sensitivity.

POSTHOC: original 2025 published/licensed individualCount of accepted
Lycaenidae, Pieridae, Papilionidae, Nymphalidae and Hesperiidae.
A positive 3-family SHARE cannot establish increased Lycaenidae abundance.
No source phosphorus dose or physically directed river exposure is measured.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random,statistics
from collections import defaultdict
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
FAMILY=ROOT+"DENMARK_2025_NEW399_SOURCE_LOCKED_ACCEPTED_BUTTERFLY_FAMILY_COUNTS.json"
SPATIAL=ROOT+"DENMARK_2025_NEW399_LYCAENIDAE_2022_WWTP_AGGREGATE_REGION_COMPONENT_SENSITIVITY.json"
OLD=ROOT+"DENMARK_2025_ORIGINAL107_LYCAENIDAE_OUTOFTIME_FAMILY_AND_MONTH_RESULTS.json"
GROUPS={
    "three_family":["Lycaenidae","Pieridae","Papilionidae"],
    "five_family":["Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae"],
}
SEED=20261110


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def counts(items,families):
    l=sum(x["family_individual_count_original"]["Lycaenidae"] for x in items)
    total=sum(sum(x["family_individual_count_original"][f] for f in families) for x in items)
    day=sum(x["original_any_Occurrence_linked_distinct_date_count"] for x in items)
    return {"sites":len(items),"original_published_survey_dates":day,"Lycaenidae_recorded_individualCount":l,
            "original_selected_family_individualCount_total":total,"proportion":l/total if total else None,
            "rate_per_published_survey_date":l/day if day else None,
            "logit":math.log((l+.5)/(total-l+.5)) if total else None}


def old_counts(data,source,families):
    a=data["original_107site_full_year_three_year_source_band_counts"][source]["2025"]
    c=a["species_counts_by_accepted_family"]
    den=sum(c[f] for f in families)
    l=c["Lycaenidae"]
    return {"Lycaenidae":l,"total_families":den,"proportion":l/den if den else None,
            "original_repeated_old_source_site_rows":a["site_locations_n"]}


def main_analysis(root):
    p={key:root/name for key,name in
       (("family_source",FAMILY),("spatial_effect",SPATIAL),("old107_benchmark",OLD))}
    d={k:json.loads(f.read_text()) for k,f in p.items()}
    family=d["family_source"];spatial=d["spatial_effect"];old=d["old107_benchmark"]
    assert family["focal_family_outcome_support_GATE"]["original_new399_2025_family_effect_descriptive_contrasts_permitted"]
    assert spatial["source_original_sha256"]["original_licensed_family_count_for_NEW399"]==sha(p["family_source"])
    all399=family["all399_source_site_original_family_counts_private_hash_aggregate_only"]
    sample=[x for x in all399 if x["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]]
    near=[x for x in sample if x["source_band"]=="near5"]
    far=[x for x in sample if x["source_band"]=="remote10"]
    assert (len(sample),len(near),len(far))==(286,70,92)
    regions=defaultdict(lambda:defaultdict(list))
    for s in sample:regions[s["original_1degree_region"]][s["source_band"]].append(s)
    shared=sorted(r for r,b in regions.items() if b["near5"] and b["remote10"])
    assert len(shared)==6

    response={}
    for key,fams in GROUPS.items():
        n=counts(near,fams);f=counts(far,fams)
        assert n["Lycaenidae_recorded_individualCount"]==382
        assert f["Lycaenidae_recorded_individualCount"]==644
        rec=[]
        for region in shared:
            aa=counts(regions[region]["near5"],fams)
            bb=counts(regions[region]["remote10"],fams)
            rec.append({
                "country_1degree_region":region,"near":aa,"remote":bb,
                "original_near_minus_remote_fraction":aa["proportion"]-bb["proportion"],
                "original_near_minus_remote_logodds":aa["logit"]-bb["logit"],
                "source_near_family_fraction_is_higher":aa["proportion"]>bb["proportion"],
            })
        dif=[r["original_near_minus_remote_logodds"] for r in rec]
        mean=statistics.mean(dif)
        random_instance=random.Random(SEED)
        draws=[statistics.mean(random_instance.choice(dif) for __ in dif) for _ in range(10000)]
        draws.sort()
        loo=[{"removed_original_region":r["country_1degree_region"],
              "geometric_mean_near_vs_remote_fractional_family_OR_remaining_five_regions":
                  math.exp(statistics.mean(x for j,x in enumerate(dif) if rec[j]["country_1degree_region"]!=r["country_1degree_region"]))}
             for r in rec]
        response[key]={
            "defined_original_families":fams,
            "new_2025_n399_event_supported_near":n,
            "new_2025_n399_event_supported_remote":f,
            "pooled_Lycaenidae_proportion_delta_near_minus_remote_pp":100*(n["proportion"]-f["proportion"]),
            "pooled_Lycaenidae_share_OR_near_divided_remote":math.exp(n["logit"]-f["logit"]),
            "actual_original_Lycaenidae_per_2025_ANY_published_date_ratio":n["rate_per_published_survey_date"]/f["rate_per_published_survey_date"],
            "original_matched_1degree_region_relative_directions":rec,
            "regions_showing_positive_near_Lycaenidae_fraction":sum(x["source_near_family_fraction_is_higher"] for x in rec),
            "equal_1degree_region_mean_log_OR":mean,
            "equal_1degree_region_geometric_mean_OR":math.exp(mean),
            "posthoc_10000_region_bootstrap_descriptive95_log_OR":[draws[int(.025*(len(draws)-1))],draws[int(.975*(len(draws)-1))]],
            "6_region_leave_one_region_out":loo,
            "old107_original_2025_near":old_counts(old,"near5",fams),
            "old107_original_2025_remote":old_counts(old,"remote10",fams),
            "selecting_three_vs_five_family_denominator_is_POSTHOC_after_new399_2025_counts":True,
        }
    if not (response["three_family"]["pooled_Lycaenidae_share_OR_near_divided_remote"]>1 and
            response["five_family"]["pooled_Lycaenidae_share_OR_near_divided_remote"]<1):
        raise ValueError("Expected originally observed three-versus-five family sign change not replayed")
    assert response["three_family"]["regions_showing_positive_near_Lycaenidae_fraction"]==2
    assert response["five_family"]["regions_showing_positive_near_Lycaenidae_fraction"]==1
    assert abs(response["three_family"]["equal_1degree_region_geometric_mean_OR"]-0.9079846220349959)<1e-10
    assert abs(response["five_family"]["equal_1degree_region_geometric_mean_OR"]-0.4044555416232778)<1e-10
    families={}
    for f in GROUPS["five_family"]:
        n=family["source_band_2025_family_original_counts_and_documented_visit_days"]["near5"]["accepted_families_per_any_published_Occurrence_linked_date"][f]
        r=family["source_band_2025_family_original_counts_and_documented_visit_days"]["remote10"]["accepted_families_per_any_published_Occurrence_linked_date"][f]
        families[f]={"original_near_individualCount_per_any_occ_linked_date":n,
                     "original_remote_individualCount_per_any_occ_linked_date":r,
                     "relative_near_vs_remote_recording_rate":n/r if r else None}
    return{
        "schema":"gesui_DK_2025_NEW399_original_3vs5_family_denominator_sign_contrast_v01",
        "original_frozen_result_sha256":{k:sha(v) for k,v in p.items()},
        "posthoc_denominator_choice_interpretation_no_confirmatory_P_value":True,
        "source_2022_EEA_outfall_geography_is_not_actual_waterborne_butterfly_exposure":True,
        "three_family_vs_five_family_denominator_response":response,
        "complete_five_family_rate_per_publicly_documented_2025_survey_day":families,
        "registered_WWTP_outfall_2025_geodesic_source_connected_components":17,
        "qualifying_new399_original_2025_source_near_sites":70,
        "qualifying_new399_original_2025_source_remote_sites":92,
        "no_evidence_of_direct_WWTP_2025_butterfly_population_change":True,
        "original_species_name_taxonomy_frozen_before_NEW399_family_count":True,
        "original_public_licensed_aggregates_only_CC_BY_NC_4_0":True,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    d=main_analysis(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:{"pooled_OR":v["pooled_Lycaenidae_share_OR_near_divided_remote"],
                          "six_region_mean_OR":v["equal_1degree_region_geometric_mean_OR"],
                          "source_regions_positive":v["regions_showing_positive_near_Lycaenidae_fraction"],
                          "near_2025_share":v["new_2025_n399_event_supported_near"]["proportion"],
                          "remote_2025_share":v["new_2025_n399_event_supported_remote"]["proportion"]}
                      for k,v in d["three_family_vs_five_family_denominator_response"].items()},indent=2))


if __name__=="__main__":
    main()
