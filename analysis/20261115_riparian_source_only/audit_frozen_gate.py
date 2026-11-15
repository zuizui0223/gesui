#!/usr/bin/env python3
"""STRICT source-only replay of the frozen 2026-11-15 riparian HNS gate.

No species HNS is opened, calculated, regressed or imputed in this script.
Do not use older 2/5 km photo gates as interchangeable with frozen 0-200 m /
500-1500 m arm opportunity counts. Unavailable is NOT an ecological zero.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path

FROZEN_PROTOCOL_GIT_BLOB = "b9ce415c36d398eee77d8318a20aef4873e0bf44"
H7_AGGREGATE_SHA256 = "1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
SOURCE_BLOBS = {
    "dk2025": "16a89e29f6bfe4aa6c39542b28f06e85c3880823",
    "de2022": "1c1af1b5be5d061e8d7d2ca901f09156f7cc40e1",
    "de2022_stop": "ec76f9a33f4021650c0a1af4cb0b8ddf971b8acf",
}
MIN_PAIRS = 24
MIN_H6 = 12
MIN_COUNTRIES = 3
MIN_CELLS_PER_PAIR = 8
MIN_FIVEFAM_RECORDS_PER_CELL = 20
MIN_SCORED_LYCAENID_SPECIES_PER_CELL = 2
MIN_OBSERVER_DATES_PER_CELL = 5
MIN_LYCAENID_SCORE_COVERAGE = 0.70
DATA_DATE = "2022-05-01/2022-09-30"
EXPECTED_BANDS_M = ("0-200", "500-1500")
FROZEN_KEYWORDS = [
    "at least 24 nonoverlapping matched",
    "at least 12 genuine H6",
    "3 countries",
    "0–200m",
    "500–1500m",
    "2022-05-01 to 2022-09-30",
    ">=20 unique source-georeferenced five-family",
    ">=2 accepted Lycaenidae",
    ">=70% Lycaenidae",
    ">=5 distinct observer-dates",
    "retrospective preregistration",
]


def git_blob(raw: bytes) -> str:
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()


def load_verified_blob(file: Path, expected_git_sha: str):
    raw=file.read_bytes()
    got=git_blob(raw)
    if got!=expected_git_sha:
        raise ValueError(f"Source Git blob SHA mismatch: {file} {got} != {expected_git_sha}")
    return json.loads(raw)


def count_relation_subset(rows, label):
    r=[x for x in rows if x.get("physical_source_2022_to_site2025_class")==label]
    return {"source_site_relations":len(r),
            "distinct_site_SHA256":len(set(x["original_public_site_SHA256"] for x in r)),
            "distinct_reported_outfall_SHA256":len(set(x["original_registered_outfall_SHA256"] for x in r)),
            "distinct_monitor_reaches":len(set(x["original_butterfly_monitor_centroid_HYRIV_ID"] for x in r)),
            "site_uncertainty_classification_not_standalone_occurrences":True}


def replay(args):
    fraw=args.freeze.read_bytes()
    sha=git_blob(fraw)
    if sha!=FROZEN_PROTOCOL_GIT_BLOB:
        raise ValueError("Frozen one-shot protocol changed or was replaced")
    frozen=fraw.decode("utf-8")
    for item in FROZEN_KEYWORDS:
        if item not in frozen: raise ValueError("Locked protocol anchor missing: "+item)
    dk=load_verified_blob(args.dk,SOURCE_BLOBS["dk2025"])
    de=load_verified_blob(args.de,SOURCE_BLOBS["de2022"])
    stop=load_verified_blob(args.de_stop,SOURCE_BLOBS["de2022_stop"])
    if hashlib.sha256(args.h7.read_bytes()).hexdigest()!=H7_AGGREGATE_SHA256:
        raise ValueError("Original H7 source SHA mismatch")
    # Only schema and H7/H6 grouping keys are accessed; HNS values are NOT read by index.
    with args.h7.open(encoding="utf-8",newline="") as f:
        reader=csv.reader(f)
        columns=next(reader)
        hidx=columns.index("h")
        h6idx=columns.index("h6")
        nrows=0
        h7ids=set()
        h6ids=set()
        for record in reader:
            nrows+=1
            h7ids.add(record[hidx])
            h6ids.add(record[h6idx])
    if nrows!=454 or len(h7ids)!=232 or len(h6ids)!=90:
        raise ValueError("Original H7 source-grain drift")
    forbidden=[
      "decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
      "occurrenceID","eventID","samplingProtocol","recordedBy",
      "individualCount","latitude","longitude"
    ]
    if any(c in columns for c in forbidden):
        raise ValueError("H7 panel schema altered: point-grain fields may exist; re-review feasibility")
    if dk["original_Danish_2022_EEA_source_points_n"]!=203:
        raise ValueError("Original DK 2022 source denominator changed")
    if dk["site_HydroRIVERS_reach_snap_le500m_n"]!=60:
        raise ValueError("DK site snap denominator drift")
    if dk["site_500m_nearest_river_is_NOT_necessarily_riparian_or_butterfly_puddling_habitat"] is not True:
        raise ValueError("Original DK source safety flag absent")
    rel=dk["all_outfall_site_directed_relation_records_without_butterfly_responses"]
    if len(rel)!=92: raise ValueError("DK original relation denominator drift")
    dkcounts={
      "relation_rows":len(rel),
      "survey_year":2025,
      "reported_source_year":2022,
      "original_public_survey_site_centroids":286,
      "source_reported_outfalls":203,
      "site_snap_le500m":int(dk["site_HydroRIVERS_reach_snap_le500m_n"]),
      "site_reach_not_snapped":sum(x["physical_source_2022_to_site2025_class"]=="SITE_REACH_NOT_SNAPPED" for x in rel),
      "same_reach_order_ambiguous":sum(x["physical_source_2022_to_site2025_class"]=="SAME_REACH_ORDER_UNKNOWN" for x in rel),
      "strict_downstream":count_relation_subset(rel,"DOWNSTREAM_CONNECTED_LE10KM"),
      "strict_upstream":count_relation_subset(rel,"UPSTREAM_CONNECTED_LE10KM"),
      "compatible_with_frozen_primary_2022_year":False,
      "original_centroid_500m_snap_NOT_0_200m_occurrence_coordinate":True,
    }
    if (dkcounts["site_reach_not_snapped"]!=69 or
        dkcounts["same_reach_order_ambiguous"]!=4 or
        dkcounts["strict_downstream"]["source_site_relations"]!=10 or
        dkcounts["strict_downstream"]["distinct_reported_outfall_SHA256"]!=2 or
        dkcounts["strict_upstream"]["source_site_relations"]!=4 or
        dkcounts["strict_upstream"]["distinct_reported_outfall_SHA256"]!=2):
        raise ValueError("Original DK directed relation values drifted")
    if (de["original_source_frame_SHA256"] if "original_source_frame_SHA256" in de else de["source_frame_SHA256"])!="fc9aa81b28aabb338a57c533eb81e5366cfc5432f7edfc5577838655b94678f2":
        raise ValueError("DE source frame mismatch")
    if de["original_total_station_radius_requests"]!=200 or de["new_2026_11_10_plant_requests_validated"]!=136:
        raise ValueError("Original DE effort source queries changed")
    d5=de["by_radius"]["5km"]
    if d5["original_source_pairs_n"]!=73 or d5["original_WISE6_independent_pair_geo_quartets"]!=5:
        raise ValueError("DE original geometry/proxy opportunity changed")
    if not de["all_predeclared_2km_and_5km_source_only_gates_fail"]:
        raise ValueError("Original DE source-only STOP unexpectedly reversed")
    if stop["actual_public_2022_butterfly_photo_API_queries"]!=0 or stop["by_radius"]["5km"]["status"]!="STOP_PREBUTTERFLY_GATES_FAILED":
        raise ValueError("Original DE butterfly biological gate had unexpected outcome work")
    dec={
      "schema":"GESUI_FROZEN_WITHIN_H7_RIPARIAN_HNS_SOURCE_ONLY_AUDIT_V1",
      "frozen_protocol_git_blob":FROZEN_PROTOCOL_GIT_BLOB,
      "date_window":DATA_DATE,
      "locked_bands_m":list(EXPECTED_BANDS_M),
      "locked_primary":"Lycaenidae incidence-equal larval-host HNS, 47 scored species",
      "minimum_complete_8cell_matched_source_reference_pairs":MIN_PAIRS,
      "minimum_distinct_H6":MIN_H6,
      "minimum_countries":MIN_COUNTRIES,
      "required_local_cells_min":MIN_PAIRS*MIN_CELLS_PER_PAIR,
      "necessary_total_5family_occurrence_opportunities_min":MIN_PAIRS*MIN_CELLS_PER_PAIR*MIN_FIVEFAM_RECORDS_PER_CELL,
      "necessary_local_cell_observer_date_opportunities_min":MIN_PAIRS*MIN_CELLS_PER_PAIR*MIN_OBSERVER_DATES_PER_CELL,
      "min_scored_lycaenid_species_each_cell":MIN_SCORED_LYCAENID_SPECIES_PER_CELL,
      "min_scored_lycaenid_record_fraction_each_cell":MIN_LYCAENID_SCORE_COVERAGE,
      "raw_source_sha_reverified":{
         "dk_git_blob":SOURCE_BLOBS["dk2025"],
         "de_git_blob":SOURCE_BLOBS["de2022"],
         "de_stop_git_blob":SOURCE_BLOBS["de2022_stop"],
         "H7_aggregate_SHA256":H7_AGGREGATE_SHA256},
      "replayed_preexisting_source_support":{
          "denmark_2025_survey_year_not_primary":dkcounts,
          "germany_2022_legacy_2km_5km_proxy_not_primary_0_200m_500_1500m":{
               "orig_pairs":d5["original_source_pairs_n"],
               "orig_site_points":100,
               "5km_legacy_photo_supported_independent_quartets":d5["original_WISE6_independent_pair_geo_quartets"],
               "5km_legacy_photo_opportunity_supported_pairs":d5["photo_opportunity_eligible_pairs"],
               "5km_legacy_outcome_api_queries":stop["actual_public_2022_butterfly_photo_API_queries"],
               "original_2km_and_5km_preoutcome_gates_failed":True,
               "not_verified_as_2022_8cell_design":True},
          "H7_454_aggregate":{"rows":454,"H7":232,"H6":90,
                "n_columns":len(columns),
                "no_occurrence_coordinates_or_observer_events":True,
                "cannot_assign_near_inland_or_2022_summer":True}},
      "qualified_8cell_source_reference_pairs_currently_verified":0,
      "zero_means_no_source_qualified_complete_pair_documented_NOT_ecological_zero":True,
      "qualifying_2022_8cell_point_cohort_assembled":False,
      "no_local_HNS_response_computed_or_fit":True,
      "biological_outcome_test_executed":False,
      "decision":"STOP_DATA_OPPORTUNITY_NOT_YET_DEMONSTRATED",
      "reason_codes":["NO_VERIFIED_2022_EIGHT_CELL_SOURCE_REFERENCE_PAIR_COHORT",
        "DK_2025_BUTTERFLY_YEAR_MISMATCH_AND_TWO_INDEPENDENT_OUTFALLS",
        "DE_2022_LEGACY_GEOMETRY_EFFORT_GATE_FAILED_NOT_SAME_NEW_DISTANCE_BANDS",
        "H7_PANEL_SPECIESKEY_AGGREGATES_NO_POINT_COORDINATES",
        "REQUIRED_OBSERVER_DAYS_AND_WITHIN_H7_MATCHED_SOURCE_REFERENCE_UNAVAILABLE"],
      "no_inference_about_presence_absence_of_wastewater_to_land_impact":True,
    }
    if dec["qualified_8cell_source_reference_pairs_currently_verified"]<MIN_PAIRS and dec["biological_outcome_test_executed"]:
        raise AssertionError("BIOLOGICAL ANALYSIS IS FORBIDDEN WHEN SOURCE GATE FAILS")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dec,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("GESUI_SOURCE_GATE_ONLY="+json.dumps({
       "status":dec["decision"],
       "original_dk_relations":len(rel),
       "dk_downstream_sites":dkcounts["strict_downstream"]["distinct_site_SHA256"],
       "dk_independent_downstream_outfalls":dkcounts["strict_downstream"]["distinct_reported_outfall_SHA256"],
       "de_old_5km_quartets_not_primary":d5["original_WISE6_independent_pair_geo_quartets"],
       "H7_nonspatial_original_rows":454,
       "verified_complete_2022_8cell_pairs":0,
       "biological_test_ran":False},sort_keys=True),flush=True)
    print("GESUI_FROZEN_SOURCE_GATE_FULL_JSON="+json.dumps(dec,ensure_ascii=False,sort_keys=True,separators=(",",":")),flush=True)
    return dec


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--dk",type=Path,required=True)
    ap.add_argument("--de",type=Path,required=True)
    ap.add_argument("--de-stop",type=Path,required=True)
    ap.add_argument("--h7",type=Path,required=True)
    ap.add_argument("--freeze",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    replay(ap.parse_args())


if __name__=="__main__":
    main()
