#!/usr/bin/env python3
"""Validate an all-plant site x period x taxon table; no empirical model is run."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

REQUIRED = (
    "site_id", "period", "taxon_id", "catchment_id", "habitat_type",
    "origin_status", "occurrence_records", "site_total_records",
    "downstream_status", "outfall_id", "network_distance_km",
    "ww_volume_m3_yr", "effluent_N_kg_yr", "effluent_N_source",
    "trait_N", "trait_source",
)
HABITATS = {"aquatic", "riparian", "terrestrial", "unknown"}
ORIGINS = {"native", "alien", "unknown"}
DIRECTIONS = {"downstream", "upstream", "unconnected", "unknown"}
N_SOURCES = {"measured_concentration_x_flow", "reported_facility_mass_load", "modeled_mass_load"}

def _nonnegative(value: str, integer: bool = False) -> float:
    x = float(value)
    if x < 0 or x != x or x in (float("inf"), float("-inf")):
        raise ValueError("negative or non-finite quantity")
    if integer and not x.is_integer():
        raise ValueError("count is not integer")
    return x

def check(path: Path) -> list[str]:
    errors: list[str] = []
    keys: set[tuple[str, str, str]] = set()
    catchments: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        cols = set(reader.fieldnames or [])
        if not set(REQUIRED) <= cols:
            return ["missing columns: " + ", ".join(sorted(set(REQUIRED) - cols))]
        for i, row in enumerate(reader, 2):
            issues = []
            def get(k: str) -> str:
                return (row.get(k) or "").strip()
            key = (get("site_id"), get("period"), get("taxon_id"))
            if not all(key) or not get("catchment_id"):
                issues.append("missing site/period/taxon/catchment key")
            if key in keys:
                issues.append("duplicate site-period-taxon")
            keys.add(key)
            if key[0] in catchments and catchments[key[0]] != get("catchment_id"):
                issues.append("site changes catchment")
            catchments[key[0]] = get("catchment_id")
            if get("habitat_type") not in HABITATS:
                issues.append("invalid habitat_type")
            if get("origin_status") not in ORIGINS:
                issues.append("invalid origin_status")
            if get("downstream_status") not in DIRECTIONS:
                issues.append("invalid downstream_status")
            try:
                n = _nonnegative(get("occurrence_records"), integer=True)
                total = _nonnegative(get("site_total_records"), integer=True)
                if total == 0 or n > total:
                    issues.append("site record denominator invalid")
            except ValueError:
                issues.append("invalid occurrence or total record count")
            if get("downstream_status") == "downstream":
                if not get("outfall_id") or not get("network_distance_km"):
                    issues.append("downstream record needs outfall and along-network distance")
            for field in ("network_distance_km", "ww_volume_m3_yr", "effluent_N_kg_yr"):
                if get(field):
                    try:
                        _nonnegative(get(field))
                    except ValueError:
                        issues.append("invalid " + field)
            if get("effluent_N_kg_yr") and get("effluent_N_source") not in N_SOURCES:
                issues.append("N mass must have a measured/reported/modeled mass source, not WW volume alone")
            if get("effluent_N_source") and not get("effluent_N_kg_yr"):
                issues.append("N source specified without N mass")
            if get("trait_N") and not get("trait_source"):
                issues.append("trait_N needs a traceable trait source")
            if issues:
                errors.append("row " + str(i) + ": " + "; ".join(issues))
    return errors

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("csv", type=Path)
    args = p.parse_args()
    problems = check(args.csv)
    if problems:
        for problem in problems:
            print(problem)
        raise SystemExit(1)
    print("PASS: site-period-taxon data contract (not an ecological analysis)")

if __name__ == "__main__":
    main()
