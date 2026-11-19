#!/usr/bin/env python3
"""SOURCE-only audit: externally verified actual sludge use vs ROSEAU candidate coverage.

An absent ROSEAU candidate is never classified as unexposed. All source records
are site/commune evidence, never an independent basin-level biotic effect.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from datetime import datetime, timezone

YEARS = {str(y) for y in range(2017, 2023)}
SOURCE_TAG = "listed_in_STEU_spreading_communes_with_positive_STEU_agricultural_tMS_NOT_field_exposure"
EVIDENCE_TAG = "PERFORMED_SLUDGE_SPREADING_IN_COMMUNE"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_candidates(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("Candidate source CSV empty")
    codes = {}
    for row in rows:
        code = (row.get("INSEE_code") or "").strip()
        if not re.fullmatch(r"\d{5}", code) or code in codes:
            raise ValueError("Invalid or duplicate INSEE source candidate code")
        years = set((row.get("candidate_source_years") or "").split(";"))
        if not years or not years <= YEARS:
            raise ValueError("Invalid candidate year codes")
        if row.get("source_interpretation") != SOURCE_TAG:
            raise ValueError("Candidate file semantics changed; cannot audit")
        codes[code] = years
    return codes


def read_external_positive(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("No externally verified source events")
    out = []
    seen = set()
    for row in rows:
        code = (row.get("INSEE_code") or "").strip()
        year = (row.get("application_year") or "").strip()
        key = (code, year, row.get("source_station_SANDRE", ""))
        if not re.fullmatch(r"\d{5}", code) or year not in YEARS or key in seen:
            raise ValueError("Invalid or duplicate external positive")
        if row.get("evidence_class") != EVIDENCE_TAG or not row.get("primary_report_url", "").startswith("https://"):
            raise ValueError("External source must document performed spreading with provenance")
        seen.add(key)
        out.append(row)
    return out


def audit(candidate_path, external_path):
    candidates = read_candidates(candidate_path)
    external = read_external_positive(external_path)
    cases = []
    for row in external:
        code, year = row["INSEE_code"], row["application_year"]
        years = candidates.get(code, set())
        cases.append({
            "INSEE_code": code,
            "application_year": year,
            "source_station_SANDRE": row["source_station_SANDRE"],
            "verified_source_evidence_class": row["evidence_class"],
            "in_ROSEAU_positive_station_commune_candidate_same_year": year in years,
            "in_ROSEAU_candidate_union_2017to22": code in candidates,
            "primary_report_url": row["primary_report_url"],
            "source_locator": row["source_locator"],
            "known_application_event_not_plant_or_butterfly_exposure": True,
        })
    missing_year = sum(not x["in_ROSEAU_positive_station_commune_candidate_same_year"] for x in cases)
    unique_source_station_years = len({(x["source_station_SANDRE"],x["application_year"]) for x in cases})
    unique_confirmed_commune_years = len({(x["INSEE_code"],x["application_year"]) for x in cases})
    return {
        "schema": "gesui_FR_ROSEAU_independent_confirmed_actual_spreading_coverage_v1",
        "run_UTC": datetime.now(timezone.utc).isoformat(),
        "source_files_SHA256": {"ROSEAU_candidate_csv": sha256(candidate_path),
                                "external_verified_csv": sha256(external_path)},
        "candidate_union_2017to22": len(candidates),
        "candidate_communes_in_2018": sum("2018" in years for years in candidates.values()),
        "external_verified_commune_station_year_cases": len(cases),
        "external_unique_commune_years": unique_confirmed_commune_years,
        "external_unique_source_station_years": unique_source_station_years,
        "these_communes_are_not_independent_plant_or_H6_replicates": True,
        "external_cases_not_in_2018_ROSEAU_candidate_list": missing_year,
        "external_cases_not_in_2017to22_union": sum(not x["in_ROSEAU_candidate_union_2017to22"] for x in cases),
        "cases": cases,
        "case_study_not_sampling_representative": True,
        "absence_from_list_is_not_absence_of_field_application": True,
        "cannot_estimate_candidate_capture_probability_from_this_case": True,
        "independent_H6_n_for_new_biotic_model": 0,
        "biological_model_source_gate": "STOP",
        "no_NUTRI_HNS_EQR_values_read": True,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--external", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = audit(args.candidates, args.external)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n")
    print("GESUI_FR_EXTERNAL_APPLIED_SLUDGE_SOURCE_COVERAGE", json.dumps({k:v for k,v in result.items() if k not in ("cases", "source_files_SHA256")}, sort_keys=True))


if __name__ == "__main__":
    main()
