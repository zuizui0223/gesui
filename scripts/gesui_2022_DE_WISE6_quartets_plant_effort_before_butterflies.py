#!/usr/bin/env python3
"""2026-11-09 ORIGINAL-COHORT source-only iNat 2019-2021 Plantae photo gate.

The original 8 WISE6 2022 German source-vs-registered-zero river
quartets are fixed independently of ANY butterfly photo or family record.
This stage only looks up pre-2022 Plantae photo opportunity, not butterflies.
It cannot establish wastewater exposure, riverbank habitat, or causality.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ORIGINAL_QUARTET_PATH = "results/wise6_2022_national_DEAT_TP_sourceonly/GERMANY_2022_WISE6_TP_8SOURCE_8ZERO_SAME_RIVER_MONITOR_QUARTETS_SOURCE_ONLY_LOCK.json"
ORIGINAL_STATIONS_PATH = "results/wise6_2022_national_DEAT_TP_sourceonly/NATIONAL_DE_AT_2022_WISE6_RIVER_TP_PUBLIC_STATION_SOURCE_ONLY_QC.json"
API = "https://api.inaturalist.org/v1/observations"
STUDY_DATE_START, STUDY_DATE_END = "2019-01-01", "2021-12-31"
TAXON_PLANTAE = 47126
DE_PLACE_ID = 7207
RADII = (2, 5)
PAUSE_SECS = 1.15
# Pre-outcome effort gates, not adjusted to achieve a butterfly-positive result.
PREDECLARED_THRESHOLDS = {
    "2km": {"min_plant_photos_per_arm": 5, "min_four_arm_minmax_ratio": 0.15, "min_eligible_quartets": 5, "min_geographic_2deg_blocks": 4},
    "5km": {"min_plant_photos_per_arm": 5, "min_four_arm_minmax_ratio": 0.15, "min_eligible_quartets": 4, "min_geographic_2deg_blocks": 3},
}
REQUIRED_SUCCESSFUL_REQUESTS = 32 * len(RADII)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def get_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def haversine(lat1, lon1, lat2, lon2):
    r = math.pi / 180
    dlat = (lat2 - lat1) * r
    dlon = (lon2 - lon1) * r
    a = math.sin(dlat/2)**2 + math.cos(lat1*r)*math.cos(lat2*r)*math.sin(dlon/2)**2
    return 12742 * math.asin(min(1., math.sqrt(a)))


def load_fixed(root):
    qfile = root / ORIGINAL_QUARTET_PATH
    sfile = root / ORIGINAL_STATIONS_PATH
    q = get_json(qfile)["selected_original_8_source_vs_clean_zero_monitor_quartets"]
    stations = get_json(sfile)["nationwide_original_source_2022_station_registry"]
    byid = {s["original_country_scheme_site_SHA256"]: s for s in stations}
    if len(q) != 8 or len(byid) != len(stations):
        raise ValueError("Expected fixed original 8 complete quartets and unique station hashes")
    rows = []
    for i, item in enumerate(q, 1):
        quartet_block = item["source"]["upstream_country2degree_block"]
        for arm in ("source", "registered_source_zero"):
            for side in ("upstream", "downstream"):
                source_id = item[arm][side + "_station_SHA256"]
                if source_id not in byid:
                    raise ValueError("Frozen quartet monitoring station hash absent from original national registry")
                s = byid[source_id]
                if s["countryCode"] != "DE":
                    raise ValueError("Original German quartet contains non-DE station")
                rows.append({
                    "original_quartet_number": i, "original_2022_source_class": arm,
                    "side_of_river_monitor_pair": side, "original_monitor_hash": source_id,
                    "countryCode": "DE",
                    "original_quartet_2deg_block": quartet_block,
                    "latitude": float(s["original_WISE6_WGS84_lat"]),
                    "longitude": float(s["original_WISE6_WGS84_lon"]),
                    "original_2022_TP_sampling_dates": int(s["original_2022_TP_independent_sample_dates"]),
                    "source_pair_original_2022_reported_P_t_y": float(item["source"]["real_exact_source_between"]["strictly_between_original_P_reported_t_year"]),
                    "original_within_pair_reach_length_km": float(item[arm]["full_reach_path_to_downstream_reach_start_km"]),
                })
    assert len(rows) == 32 and len({x["original_monitor_hash"] for x in rows}) == 32
    assert all(x["original_2022_TP_sampling_dates"] >= 4 for x in rows)
    for i, x in enumerate(rows):
        x["min_distance_to_other_original_monitor_km"] = min(
            haversine(x["latitude"], x["longitude"], y["latitude"], y["longitude"])
            for j, y in enumerate(rows) if i != j
        )
    return rows, {"original_quartets_sha256": digest(qfile), "national_sites_sha256": digest(sfile)}


def query_plant(row, radius):
    args = {
        "taxon_id": TAXON_PLANTAE, "place_id": DE_PLACE_ID,
        "lat": round(row["latitude"], 5),
        "lng": round(row["longitude"], 5),
        "radius": radius,
        "d1": STUDY_DATE_START, "d2": STUDY_DATE_END,
        "photos": "true", "geo": "true",
        "geoprivacy": "open", "taxon_geoprivacy": "open",
        "quality_grade": "research", "captive": "false",
        "acc_below": 500, "per_page": 1
    }
    url = API + "?" + urlencode(args)
    last_error = None
    for attempt in range(4):
        try:
            req = Request(url, headers={
                "User-Agent": "gesui-lycaenidae-2022-directed-river-effort-audit/0.1 (GitHub zuizui0223/gesui)",
                "Accept": "application/json",
            })
            with urlopen(req, timeout=40) as response:
                raw = response.read()
                status = response.status
            payload = json.loads(raw)
            n = payload["total_results"]
            if status != 200 or not isinstance(n, int) or n < 0:
                raise ValueError("Original iNaturalist response lacked valid total_results")
            return {
                "status": status,
                "original_API_url": url,
                "original_JSON_response_SHA256": hashlib.sha256(raw).hexdigest(),
                "count_public_2019_2021_research_grade_georeferenced_plantae_photos": n,
                "no_butterfly_outcome_was_requested": True,
            }
        except (HTTPError, URLError, ValueError, TimeoutError) as e:
            last_error = str(e)
            if attempt < 3:
                time.sleep(5 * (attempt + 1))
    return {
        "status": "FAILED",
        "original_API_url": url, "error": last_error,
        "count_public_2019_2021_research_grade_georeferenced_plantae_photos": None,
        "no_butterfly_outcome_was_requested": True,
    }


def evaluate(rows):
    report = {}
    for rad in RADII:
        name = f"{rad}km"
        grouped = {}
        for row in rows:
            grouped.setdefault(row["original_quartet_number"], []).append(row)
        gates = []
        for idx, arms in sorted(grouped.items()):
            dists = [
                haversine(a["latitude"], a["longitude"], b["latitude"], b["longitude"])
                for j, a in enumerate(arms) for b in arms[j+1:]
            ]
            # Explicit buffer disjointness within the quartet is a nonnegotiable gate;
            # 5km overlapping circles may NOT be independent contrasts.
            disjoint = min(dists) > 2*rad
            val = [a["independent_previous_year_plantae_photo_counts"][name]["count_public_2019_2021_research_grade_georeferenced_plantae_photos"] for a in arms]
            no_error = all(isinstance(v, int) for v in val)
            ratio = min(val)/max(val) if no_error and max(val) > 0 else 0.0
            t = PREDECLARED_THRESHOLDS[name]
            passed = bool(no_error and disjoint and min(val) >= t["min_plant_photos_per_arm"] and ratio >= t["min_four_arm_minmax_ratio"])
            gates.append({
                "original_quartet_number": idx,
                "geographic_source_2degree_block": arms[0]["original_quartet_2deg_block"],
                "all_4_sampling_arm_source_photo_sufficient": passed,
                "four_circles_nonoverlap_within_quartet": disjoint,
                "min_pairwise_km": min(dists),
                "four_site_plant_photo_counts_by_original_arm": {
                    a["original_2022_source_class"] + "_" + a["side_of_river_monitor_pair"]: a["independent_previous_year_plantae_photo_counts"][name]["count_public_2019_2021_research_grade_georeferenced_plantae_photos"] for a in arms
                },
                "minimum_to_maximum_plant_photo_effort_ratio": ratio,
                "all_queries_valid": no_error,
            })
        valid = [a for a in gates if a["all_4_sampling_arm_source_photo_sufficient"]]
        distinct = len({a["geographic_source_2degree_block"] for a in valid})
        t = PREDECLARED_THRESHOLDS[name]
        report[name] = {
            "predeclared_gate": t,
            "original_quartets": len(gates),
            "eligible_quartets": len(valid),
            "eligible_geographic_blocks": distinct,
            "photo_outcome_query_allowed": bool(
                len(valid) >= t["min_eligible_quartets"]
                and distinct >= t["min_geographic_2deg_blocks"]
            ),
            "disjoint_source_only_quartet_details": gates,
        }
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    rows, sha = load_fixed(root)
    for n, row in enumerate(rows):
        row["independent_previous_year_plantae_photo_counts"] = {}
        for radius in RADII:
            q = query_plant(row, radius)
            row["independent_previous_year_plantae_photo_counts"][f"{radius}km"] = q
            print(f"plant-only source QC {n+1}/32 {radius}km {row['original_quartet_number']} "
                  f"{row['original_2022_source_class']}/{row['side_of_river_monitor_pair']}: "
                  f"{q['count_public_2019_2021_research_grade_georeferenced_plantae_photos']}", flush=True)
            time.sleep(PAUSE_SECS)
    gate = evaluate(rows)
    all_queries = sum(1 for row in rows for v in row["independent_previous_year_plantae_photo_counts"].values()
                      if v["status"] == 200)
    out = {
        "schema": "gesui_original32_WISE6_quartet_2019_2021_plant_photo_prebutterfly_outcome_gate_v01",
        "UTC": datetime.now(timezone.utc).isoformat(),
        "original_frozen_sourceSHA256": sha,
        "date_min": STUDY_DATE_START, "date_max": STUDY_DATE_END,
        "no_butterfly_species_or_photos_requested": True,
        "original_monitoring_sites_n": 32, "original_quartets_n": 8,
        "independent_original_2degree_blocks_n": len({r["original_quartet_2deg_block"] for r in rows}),
        "radii_km": list(RADII), "plant_taxon_id": TAXON_PLANTAE,
        "successful_queries": all_queries, "planned_queries": REQUIRED_SUCCESSFUL_REQUESTS,
        "any_query_failed_not_biological_zero": all_queries < REQUIRED_SUCCESSFUL_REQUESTS,
        "prebutterfly_effort_gate": gate,
        "original_station_source_without_butterfly_outcomes": rows,
        "interpretation": "Photo effort gate, not nutrient or butterfly benefit. WWTP reported 2022 source versus clean registered-zero river paths from original source-only selection; clean does not mean unpolluted. 2019-2021 plant photos are a rough independent observer opportunity proxy. Even a passed 2022 butterfly photo gate cannot establish adult substrate feeding or causal wastewater effects."
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "all_source_requests_valid": all_queries == REQUIRED_SUCCESSFUL_REQUESTS,
        "source_effort_gate": {k: {
            "eligible_quartets": v["eligible_quartets"], "eligible_regions": v["eligible_geographic_blocks"],
            "go_to_butterfly_outcome": v["photo_outcome_query_allowed"]
        } for k, v in gate.items()}
    }, indent=2), flush=True)
    # Partial result is preserved by the workflow even if a public API request failed.


if __name__ == "__main__":
    main()
