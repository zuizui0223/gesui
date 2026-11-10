#!/usr/bin/env python3
"""Freeze original 2025 DK eBMS survey effort on SAME 107 2023/24 site locations.

No 2025 butterfly occurrence data, individuals, taxa, or effect are read here.
Original 2022 EEA registered outfall proximity stays unchanged (4 source clusters).
"""
from __future__ import annotations
import argparse, hashlib, json, math
from collections import Counter, defaultdict
from datetime import datetime,timezone
from pathlib import Path

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
SITE=ROOT+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
TAX=ROOT+"DENMARK_EBMS_V14_ORIGINAL_SCIENTIFICNAME_GBIF_FAMILY_BLIND_TAXONOMY.json"
PANEL=ROOT+"DENMARK_EBMS_2023_2024_FAMILY_SITEYEAR_COUNTS_QC_BEFORE_WWTP_CONTRAST.json"
PREVIOUS=ROOT+"DENMARK_2023_2024_SAME_SITE_MONTH_LYCAENIDAE_POSTHOC_FAMILY_SOURCE_CLUSTER_QC.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def hav(a,b):
    r=math.pi/180
    dlat=(a["lat"]-b["lat"])*r
    dlon=(a["lon"]-b["lon"])*r
    z=math.sin(dlat/2)**2+math.cos(a["lat"]*r)*math.cos(b["lat"]*r)*math.sin(dlon/2)**2
    return 12742*math.asin(min(1.,math.sqrt(z)))

def freeze(root):
    paths={k:root/v for k,v in {"original_site_and_source":SITE,"frozen_taxonomy":TAX,"original_2023_24_siteyear_counts":PANEL,"earlier_season_analysis":PREVIOUS}.items()}
    obj={k:json.loads(p.read_text(encoding="utf-8")) for k,p in paths.items()}
    src=obj["original_site_and_source"];tax=obj["frozen_taxonomy"]
    if src["butterfly_species_2023_2024_occurrence_extension_READ"] is not False:
        raise ValueError("Source-only original data was contaminated with butterfly response")
    if not tax["original_individualCount_values_NOT_READ"]:
        raise ValueError("Taxonomy not originally fixed independently of butterfly counts")
    if src["original_eBMS_2026_v1_4_source"]["source_archive_SHA256"]!="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b":
        raise ValueError("Original Denmark archive SHA differs")
    if len(src["original_2022_EEA_DK_qualifying_Fw_outfalls"])!=203:
        raise ValueError("Registered original DK source cohort unexpectedly changed")
    prev=obj["original_2023_24_siteyear_counts"]
    if len(prev["site_year_family_counts_original_public_CC_BY_NC_aggregate"])!=214:
        raise ValueError("Original 107 site two-year count panel not preserved")
    cohort=[s for s in src["original_2023_2024_repeated_event_site_source_metadata"]
            if s["repeat_2023_2024_dates_gate"]
            and s["precise_500m_explicit_uncertainty_gate"]
            and s["rounded5_site_coord_scatter_le2m"]]
    if len(cohort)!=107 or len({s["site_id_SHA256"] for s in cohort})!=107:
        raise ValueError("Original 107-site cohort must not be reselected for 2025 response")
    near=[s for s in cohort if s["nearest_qualified_EEA2022_DK_FW_outfall_km"]<=5]
    far=[s for s in cohort if s["nearest_qualified_EEA2022_DK_FW_outfall_km"]>=10]
    mid=[s for s in cohort if 5<s["nearest_qualified_EEA2022_DK_FW_outfall_km"]<10]
    if (len(near),len(far),len(mid))!=(23,38,46):raise ValueError("Same original exposure bands changed")
    sources=src["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    nearest={}
    for s in near:
        ordered=sorted((hav(s,x),x["outfall_id_SHA256"]) for x in sources)
        if abs(ordered[0][0]-s["nearest_qualified_EEA2022_DK_FW_outfall_km"])>1e-8:
            raise ValueError("Original source-proximity geometry changed")
        nearest[s["site_id_SHA256"]]=ordered[0][1]
    if len(set(nearest.values()))!=4:
        raise ValueError("Original four registered near-source systems are not preserved")
    all_valid=all(s["visit_days_2025_descriptive"]>=2 for s in cohort)
    grouped={}
    for name,arr in (("near5",near),("remote10",far),("intermediate",mid)):
        grouped[name]={
            "original_2023_24_eligible_sites":len(arr),
            "2025_sites_with_at_least_two_original_distinct_survey_days":sum(s["visit_days_2025_descriptive"]>=2 for s in arr),
            "2025_original_distinct_survey_days":sum(s["visit_days_2025_descriptive"] for s in arr),
            "2025_days_min_max":[min(s["visit_days_2025_descriptive"] for s in arr),max(s["visit_days_2025_descriptive"] for s in arr)],
            "2025_existing_DWCA_Event_metadata_only_not_butterfly_counts":True,
        }
    shared23={s["country_1deg_lat_lon_cell"] for s in near}&{s["country_1deg_lat_lon_cell"] for s in far}
    do_gate=bool(all_valid and len(near)>=18 and len(far)>=30 and len(shared23)>=3 and len(set(nearest.values()))==4)
    # The existence of site-year opportunity is source stage; not a species-detection gate.
    if not do_gate: raise ValueError("Original 2025 source-only gate fails; DO NOT open species outcomes")
    out={
        "schema":"gesui_denmark_2025_same107_original_source_ONLY_time_holdout_before_individualCount_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "is_prior_to_2025_family_individualCount_read_in_this_study":True,
        "original_2022_EEA_source_fixed_not_a_2025_exposure_change":True,
        "reuse_original_107_2023_2024_sites_without_response_selection":True,
        "all_original_2023_2024_family_results_known_before_2025_design":True,
        "original_source_file_SHA256":{k:sha(p) for k,p in paths.items()},
        "original_DK_2022_qualifying_registered_FW_outfall_n":203,
        "original_near_registered_independent_outfall_systems":4,
        "source_system_near_site_rows":dict(Counter(nearest.values())),
        "shared_1degree_near_remote_cells":sorted(shared23),
        "repeated_2025_source_only_group_sites_and_days":grouped,
        "source_sample_gate_prespecified": {
            "min_days_2025_per_each_frozen_site":2,
            "min_2025_near_sites":18,"min_2025_remote_sites":30,
            "min_1degree_cells_containing_near_and_remote":3,
            "min_distinct_registered_WWTP_source_clusters":4,
            "2025_family_outcome_access_permitted":do_gate
        },
        "original_2025_repeated_site_year_metadata_fixed":[
            {"site_id_SHA256":s["site_id_SHA256"],"region_1degree":s["country_1deg_lat_lon_cell"],
             "2025_original_visit_days":s["visit_days_2025_descriptive"],
             "2023_original_visit_days":s["visit_days_2023"],
             "2024_original_visit_days":s["visit_days_2024"],
             "2022_source_distance_km":s["nearest_qualified_EEA2022_DK_FW_outfall_km"],
             "source_band":"near5" if s in near else "remote10" if s in far else "intermediate",
             "nearest_original_WWTP_source_id_SHA256":nearest.get(s["site_id_SHA256"])}
            for s in sorted(cohort,key=lambda s:s["site_id_SHA256"])
        ],
        "2025_butterfly_occurrence_data_rows_opened":0,
        "2025_Lycaenidae_individualCount_read":False,
        "no_new_registered_WWTP_or_treatment_status_checked_after_outcomes":True,
    }
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=freeze(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"groups":out["repeated_2025_source_only_group_sites_and_days"],
                      "shared_1degree_cells":out["shared_1degree_near_remote_cells"],
                      "source_systems":out["original_near_registered_independent_outfall_systems"],
                      "2025_butterfly_data_access":out["source_sample_gate_prespecified"]["2025_family_outcome_access_permitted"]},indent=2))

if __name__=="__main__":
    main()
