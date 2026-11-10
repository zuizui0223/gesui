#!/usr/bin/env python3
"""Independent offline audit: 2025 near vs remote Lycaenidae month/year QC.

Do not re-open licensed raw observations here; trust only SHA-tracked aggregate
source files and independently recompute group-family per-documented-date
ratios and original four-source exposure effective sample.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
ORIGINAL=ROOT+"DENMARK_2025_ORIGINAL107_LYCAENIDAE_OUTOFTIME_FAMILY_AND_MONTH_RESULTS.json"
LINK=ROOT+"DK_EBMS_2025_SAME_SITE_ANY_OCCURRENCE_LINKED_DAY_QUALITY_SENSITIVITY.json"
MONTH=ROOT+"DENMARK_2025_SAME_MONTH_OCCURRENCE_LINKED_FAMILY_SENSITIVITY.json"
QA=ROOT+"DK_EBMS_2025_ORIGINAL_EVENT_OCCURRENCE_PUBLICATION_MISSINGNESS_QA.json"
FAMILIES=("Lycaenidae","Pieridae","Nymphalidae","Hesperiidae")

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def annual_rate(panel,fam):
    v=panel["recorded_individualCount_family"][fam]
    d=panel["published_Occurrence_linked_distinct_survey_days"]
    return v/d if d else None

def monthly_rate(panel,fam):
    v=panel["family_individualCount_original"][fam]
    d=panel["any_Occurrence_linked_distinct_survey_days_in_matched_months"]
    return v/d if d else None

def reconstruct(root):
    files={name:root/f for name,f in (
       ("original_2025",ORIGINAL),("2024_25_linked_same_site",LINK),
       ("month_original_replay",MONTH),("Event_Occurrence_link_quality",QA))}
    a={k:json.loads(p.read_text()) for k,p in files.items()}
    orig=a["original_2025"];link=a["2024_25_linked_same_site"]
    month=a["month_original_replay"];q=a["Event_Occurrence_link_quality"]
    if not orig["all_prior_214_2023_2024_family_individual_sums_and_days_reproduced"]:
        raise ValueError("Original 214 year panel was not reproduced")
    if not month["original_2025_family_individualCount_reproduced_exactly"]:
        raise ValueError("Original 2025 family counts changed")
    if q["source_original_row_qc"]["full_original_Event_core_rows"]!=15606 or q["source_original_row_qc"]["full_original_Occurrence_extension_rows"]!=21617:
        raise ValueError("Original eBMS Event or Occurrence source count mismatch")
    normal=link["quality_matched_2024_2025_linked_2plus_same_sites"]
    assert normal["eligible_same_source_site_n"]["near5"]==14
    assert normal["eligible_same_source_site_n"]["remote10"]==27
    annual={}
    for fam in FAMILIES:
        n=normal["group_year"]["near5"]
        r=normal["group_year"]["remote10"]
        ra=annual_rate(n["2025"],fam)/annual_rate(n["2024"],fam) if annual_rate(n["2024"],fam) else None
        rb=annual_rate(r["2025"],fam)/annual_rate(r["2024"],fam) if annual_rate(r["2024"],fam) else None
        old=normal["relative_Lycaenidae_change_2024_to_2025"][fam]
        if ra is None or rb is None or abs((ra/rb)-old["rate_ratio_of_ratios"])>1e-12:
            raise ValueError("Published same-source 2024/25 change fails independent replay")
        annual[fam]=ra/rb
    specific={}
    for label,v in month["configurations"].items():
        counts=v["all_years"]
        near=counts["near5"];far=counts["remote10"]
        for year in v["matched_years"]:
            assert counts["near5"][str(year)]["qualified_site_rows"]==v["matched_site_counts_by_source_band"]["near5"]
            assert counts["remote10"][str(year)]["qualified_site_rows"]==v["matched_site_counts_by_source_band"]["remote10"]
        ratios={}
        if 2024 in v["matched_years"] and 2025 in v["matched_years"]:
            for fam in FAMILIES:
                nr=monthly_rate(near["2025"],fam)/monthly_rate(near["2024"],fam) if monthly_rate(near["2024"],fam) else None
                fr=monthly_rate(far["2025"],fam)/monthly_rate(far["2024"],fam) if monthly_rate(far["2024"],fam) else None
                x=v["relative_2024_2025_change"][fam]
                if nr is None or fr is None or abs(x["near_over_remote_rate_ratio_ratio"]-nr/fr)>1e-12:
                    raise ValueError("Original source linked-month rate ratio failed")
                ratios[fam]=nr/fr
        outfall=v["outfall_specific_near_source_2025"]
        assert len(outfall)<=4 and len(outfall)>=2
        for group in outfall:
            assert group["2025_group"]["any_Occurrence_linked_distinct_survey_days_in_matched_months"]>0
        specific[label]={
            "qualified_public_site_rows":v["matched_site_counts_by_source_band"],
            "independent_near_registered_source_groups_with_2025_recorded_Lycaenidae":v["WWTP_near_cluster_count_with_reported_individuals"],
            "source_groups_with_any_published_Occurrence_in_matched_months":len(outfall),
            "relative_2024_to_2025_Lycaenidae_annual_rate_ratio_near_over_remote":ratios.get("Lycaenidae"),
            "all_family_relative_2024_to_2025_recording_rate_changes":ratios,
            "2025_near_Lycaenidae_per_published_occurrence_date":monthly_rate(near["2025"],"Lycaenidae"),
            "2025_remote_Lycaenidae_per_published_occurrence_date":monthly_rate(far["2025"],"Lycaenidae"),
            "2025_near_vs_remote_Lycaenidae_relative_recording_rate":monthly_rate(near["2025"],"Lycaenidae")/monthly_rate(far["2025"],"Lycaenidae"),
            "source_cluster_leave_one_out_relative_change":v["four_outfall_leave_one_out_2024_2025_ratio_sensitivity"],
        }
    assert abs(specific["2024_2025_all_common_months"]["relative_2024_to_2025_Lycaenidae_annual_rate_ratio_near_over_remote"]-2.2748932393553276)<1e-9
    assert abs(annual["Lycaenidae"]-2.2965577342047934)<1e-9
    assert specific["2024_2025_all_common_months"]["qualified_public_site_rows"]["near5"]==15
    assert specific["2024_2025_all_common_months"]["qualified_public_site_rows"]["remote10"]==28
    return{
        "schema":"gesui_DK_2025_sameMonth_occ_link_original_family_audit_independent_QA_v01",
        "original_sources_SHA256":{k:sha(p) for k,p in files.items()},
        "source_year_WWTP_2022_static_neither_measured_2025_dose_nor_randomized_treatment":True,
        "original_2025_outcome_known_before_this_secondary_quality_audit":True,
        "accepted_families_replayed_without_re-reading_original_individual_identifiers":True,
        "Event_with_no_published_Occurrence_is_UNCERTAIN_not_confirmed_zero_butterflies":True,
        "independent_original_2024_2025_2plus_linked_days_source_band_ratio_of_ratios":annual,
        "2023_2024_2025_threeyear_matched_2plus_linked_same_site_count":link["quality_matched_2023_2024_2025_linked_2plus_same_sites"]["eligible_same_source_site_n"],
        "original_published_Occurrence_linked_same_month_stratum_results":specific,
        "study_replicated_source_clusters_too_few_for_confident_population_sewage_claim":True,
    }

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=reconstruct(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"twoYearLinked":x["independent_original_2024_2025_2plus_linked_days_source_band_ratio_of_ratios"],
                      "seasonMatched":{k:v["relative_2024_to_2025_Lycaenidae_annual_rate_ratio_near_over_remote"] for k,v in x["original_published_Occurrence_linked_same_month_stratum_results"].items()}},indent=2))
