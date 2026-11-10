#!/usr/bin/env python3
"""Independent OFFLINE arithmetic/source-cohort QA of DK2025 river-direction gate.

Checks recorded source-only reach pairs and joins, no original eBMS individual
or family counts, no raw raster or licensed observation. This is NOT an
independent re-read of the third-party HydroRIVERS raw graph, whose original
source SHA was checked by the separate source-only workflow.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from datetime import datetime,timezone
from collections import Counter,defaultdict
from pathlib import Path

D="results/ebms_denmark_2023_2024_WWTP_source_only/"
SOURCE=D+"DK2025_NEW399_286SITE_203WWTP_SHA_LOCKED_HYDRORIVERS_DIRECTED_SOURCE_FEASIBILITY.json"
PRE=D+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
EVENT=D+"DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
RELATIONS={
    "DOWNSTREAM_CONNECTED_LE10KM","DOWNSTREAM_CONNECTED_10TO50KM",
    "UPSTREAM_CONNECTED_LE10KM","UPSTREAM_CONNECTED_10TO50KM",
    "SITE_REACH_NOT_SNAPPED","SOURCE_REACH_NOT_SNAPPED",
    "SAME_REACH_ORDER_UNKNOWN",
    "BOTH_SNAPPED_NO_VERIFIED_DIRECTED_ROUTE_LE50KM"
}
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def unique_count(items,k):
    return len({x[k] for x in items})

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    fp={"directional":root/SOURCE,"original399_source":root/PRE,"originalEventEEA":root/EVENT}
    files={k:json.loads(v.read_text()) for k,v in fp.items()}
    d=files["directional"]
    prev=files["original399_source"]
    e=files["originalEventEEA"]
    assert d["source_only_no_2025_butterfly_taxon_or_individualCount_read"] is True
    assert d["third_party_HydroRIVERS_not_official_verified"] is True
    assert d["original_source_sha256"]["original_2025_new399_source_only_site_cohort"]==sha(fp["original399_source"])
    assert d["original_source_sha256"]["original_eBMS_Event_only_EEA2022_outfalls"]==sha(fp["originalEventEEA"])
    assert d["original_source_sha256"]["third_party_Europe_HydroRIVERS_v1_geoparquet"]=="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
    assert prev["no_2025_butterfly_family_response_accessed_for_new_sites"] is True
    assert prev["source_eligible_new_sites_with_2plus_any_published_occurrence_dates_by_band"]=={
        "intermediate":124,"near5":70,"remote10":92
    }
    orig_sites={x["site_id_SHA256"]:x for x in prev["all_original_399_source_only_new_site_2025_documented_visit_metadata"]}
    good={k:v for k,v in orig_sites.items() if v["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]}
    assert len(good)==286
    source203={x["outfall_id_SHA256"] for x in e["original_2022_EEA_DK_qualifying_Fw_outfalls"]}
    assert len(source203)==203
    snaps=d["source_only_original_286_butterfly_monitor_centroids_nearest_HYRIV"]
    outlets=d["source_only_original_203_EEA_outfall_nearest_HYRIV"]
    assert len(snaps)==286 and len(outlets)==203
    near_sites={k for k,v in good.items() if v["source_band"]=="near5"}
    assert len(near_sites)==70
    assert {r["site_id_SHA256"] for r in snaps}==set(good)
    assert {r["outfall_SHA256"] for r in outlets}==source203
    for s in snaps:
        if s["original_source_band"]!=good[s["site_id_SHA256"]]["source_band"]:
            raise ValueError("Site source-geography label changed")
        if s["original_linked_2025_monitor_dates"]!=good[s["site_id_SHA256"]]["original_linked_occurrence_dates"]:
            raise ValueError("Source-only original Event linkage count changed")
        if (s["river_reach_ID"] is None)!=(s["site_500m_river_snap_m"] is None):
            raise ValueError("Source river reach ID and nearest distance inconsistent")
        if s["site_500m_river_snap_m"] is not None and not 0<=s["site_500m_river_snap_m"]<=500:
            raise ValueError("Site far from predeclared 500m river snap")
    for s in outlets:
        if (s["river_reach_ID"] is None)!=(s["outfall_500m_river_snap_m"] is None):
            raise ValueError("Original registered outlet snap incomplete")
        if s["outfall_500m_river_snap_m"] is not None and not 0<=s["outfall_500m_river_snap_m"]<=500:
            raise ValueError("Outfall beyond predeclared 500m river snap")
    site_lookup={s["site_id_SHA256"]:s for s in snaps}
    source_lookup={s["outfall_SHA256"]:s for s in outlets}
    pair=d["all_outfall_site_directed_relation_records_without_butterfly_responses"]
    if len(pair)!=92 or unique_count(pair,"original_public_site_SHA256")!=70:
        raise ValueError("Source-distance<=5km original pair count changed")
    if len({(x["original_public_site_SHA256"],x["original_registered_outfall_SHA256"]) for x in pair})!=len(pair):
        raise ValueError("Duplicate original WWTP×monitor source pair")
    for r in pair:
        label=r["physical_source_2022_to_site2025_class"]
        if label not in RELATIONS:raise ValueError("Unknown original direction category")
        if r["original_public_site_SHA256"] not in near_sites:
            raise ValueError("Near outfall pair not from original fixed near5 sites")
        if r["original_registered_outfall_SHA256"] not in source203:
            raise ValueError("Outfall ID not in original official registered 203 sources")
        s=site_lookup[r["original_public_site_SHA256"]]
        w=source_lookup[r["original_registered_outfall_SHA256"]]
        if r["original_butterfly_monitor_centroid_HYRIV_ID"]!=s["river_reach_ID"]:
            raise ValueError("Source site snap reach not identical across records")
        if r["original_source_outfall_HYRIV_ID"]!=w["river_reach_ID"]:
            raise ValueError("Registered outfall snap reach mismatch")
        if r["original_geodesic_WWTP_to_site_km"]>5.00000001:
            raise ValueError("Original 5km source candidate selection changed")
        dr=r["directed_source_to_site_reach_approx_route_km"]
        up=r["directed_site_to_source_reach_approx_route_km"]
        if label=="DOWNSTREAM_CONNECTED_LE10KM":
            if not (dr is not None and 0<=dr<=10 and up is None and s["river_reach_ID"]!=w["river_reach_ID"]):
                raise ValueError("Invalid original STRICT downstream direction")
        if label=="UPSTREAM_CONNECTED_LE10KM":
            if not (up is not None and 0<=up<=10 and dr is None and s["river_reach_ID"]!=w["river_reach_ID"]):
                raise ValueError("Invalid original STRICT upstream direction")
        if label=="SAME_REACH_ORDER_UNKNOWN":
            if not(s["river_reach_ID"] is not None and s["river_reach_ID"]==w["river_reach_ID"]):
                raise ValueError("Same reach mislabeled a strict direction")
        if label=="SITE_REACH_NOT_SNAPPED" and s["river_reach_ID"] is not None:
            raise ValueError("Previously unsnapped site switched source status")
        if label=="SOURCE_REACH_NOT_SNAPPED" and w["river_reach_ID"] is not None:
            raise ValueError("Previously unsnapped WWTP switched source status")
    stats=Counter(r["physical_source_2022_to_site2025_class"] for r in pair)
    assert dict(stats)==d["site_outfall_geodesic_within5km_pair_relation_counts"]
    assert stats==Counter({
        "DOWNSTREAM_CONNECTED_LE10KM":10,"UPSTREAM_CONNECTED_LE10KM":4,
        "SAME_REACH_ORDER_UNKNOWN":4,"BOTH_SNAPPED_NO_VERIFIED_DIRECTED_ROUTE_LE50KM":3,
        "SITE_REACH_NOT_SNAPPED":69,"SOURCE_REACH_NOT_SNAPPED":2
    })
    down=[r for r in pair if r["physical_source_2022_to_site2025_class"]=="DOWNSTREAM_CONNECTED_LE10KM"]
    upstream=[r for r in pair if r["physical_source_2022_to_site2025_class"]=="UPSTREAM_CONNECTED_LE10KM"]
    dc=d["strict_site10km_DOWNSTREAM_original_WWTP_groups"]
    uc=d["strict_site10km_UPSTREAM_direction_control"]
    for rows,expected in ((down,dc),(upstream,uc)):
        assert len({x["original_public_site_SHA256"] for x in rows})==expected["distinct_near_source_sites"]
        assert len({x["original_registered_outfall_SHA256"] for x in rows})==expected["distinct_reported_FW_outfall_IDs"]
        assert len({x["source_original_geographic_1degree_region"] for x in rows})==expected["independent_original_1degree_regions"]
    assert dc["distinct_near_source_sites"]==10 and dc["distinct_reported_FW_outfall_IDs"]==2
    assert uc["distinct_near_source_sites"]==4 and uc["distinct_reported_FW_outfall_IDs"]==2
    assert len({x["original_butterfly_monitor_centroid_HYRIV_ID"] for x in down})==2
    assert len({x["original_butterfly_monitor_centroid_HYRIV_ID"] for x in upstream})==2
    assert d["fixed_source_status_pre_butterfly_FEASIBILITY_GATE"]["new_hydrologically_directed_2025_Lycaenidae_family_comparison_Permitted"] is False
    out={
       "schema":"gesui_DK_2025_286site_registeredWWTP203_HydroRIVERS_direction_source_only_offline_QA_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "source_files_sha256":{k:sha(v) for k,v in fp.items()},
       "original_river_source_digest_verified_by_original_action_not_reopened_in_this_offline_audit":True,
       "original_286_event_only_2025_monitor_site_keys_reproduced":True,
       "original_203_registered_outfalls_reproduced":True,
       "original_near5_source_monitor_site_count":len(near_sites),
       "original_site_source_pairs_le5km":len(pair),
       "original_direction_classification_counts":dict(stats),
       "independent_near5_observer_points_that_snap_within_500m":sum(s["river_reach_ID"] is not None for s in snaps if s["original_source_band"]=="near5"),
       "strict_DOWNSTREAM_reach_path":{
           "public_site_distinct_n":unique_count(down,"original_public_site_SHA256"),
           "actual_distinct_butterfly_monitor_river_reach_ids":unique_count(down,"original_butterfly_monitor_centroid_HYRIV_ID"),
           "actual_distinct_registered_outfall_IDs":unique_count(down,"original_registered_outfall_SHA256"),
           "actual_distinct_2022_outfall_river_reach_ids":unique_count(down,"original_source_outfall_HYRIV_ID"),
           "original_1degree_geography_units":unique_count(down,"source_original_geographic_1degree_region"),
           "site_count_per_original_registered_WWTP_ID":dict(Counter(x["original_registered_outfall_SHA256"] for x in down)),
           "note":"These are network-related 2025 eBMS site CENTROIDS only, not sampled wastewater-receiving butterfly habitat; ten sites collapse to two downstream reach IDs."
       },
       "strict_UPSTREAM_reach_path":{
           "distinct_eBMS_site_centroids":unique_count(upstream,"original_public_site_SHA256"),
           "distinct_registered_2022_EEA_outfalls":unique_count(upstream,"original_registered_outfall_SHA256"),
           "distinct_original_1degree_regions":unique_count(upstream,"source_original_geographic_1degree_region"),
       },
       "original_prespecified_downstream_and_upstream_site_source_support_PASS":False,
       "new_2025_Lycaenidae_outcomes_not_read_to_estimate_directed_WWTP_effect":True,
       "explicit_limit_third_party_HydroRIVERS_source_not_validated_against_official":True,
       "registered_2022_outfall_location_is_not_measured_2025_water_chemistry_or_puddling_site":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in out.items() if k in (
        "original_direction_classification_counts","original_near5_source_monitor_site_count",
        "strict_DOWNSTREAM_reach_path","strict_UPSTREAM_reach_path",
        "original_prespecified_downstream_and_upstream_site_source_support_PASS")},indent=2))

if __name__=="__main__":
    main()
