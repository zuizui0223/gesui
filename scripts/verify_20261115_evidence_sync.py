#!/usr/bin/env python3
"""Verify exact source-aligned gesui 2026-11-15 derived CSV result tables.
This CI does NOT refit local models, reproduce all 3508 files, or identify causal SEM.
"""
import csv
import hashlib
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "results" / "20261115_sync"
TABLES = {
    "HNS_MODEL_SUMMARY.csv": ("c2a9d8be7501026663c2744ffd672c0d70b9d17680c90c9dc2c66918ac0d852a", 15),
    "HNS_GEOGRAPHY_AND_TEMPORAL.csv": ("49233615c5813251277fe32263e2d963498ac61db4cc84dee61b637e23fd8bd0", 5),
    "BUTTERFLY_FAMILY_NETWORK_PREDICTION.csv": ("874c83dcbd9d52938cdddc9f70c96e93c71c0681dd1166592da79e69ab0be1f6", 12),
    "PIERIDAE_SOURCE_AND_SPECIES_SENSITIVITY.csv": ("13baf6ec66a891a76c80c9a5de4216fb6a078a25f8e6e9e3f364fbd574fab92d", 8),
}

def near(value, expected):
    assert math.isfinite(float(value)) and abs(float(value) - expected) <= 1e-9, (value, expected)

def run():
    data = {}
    for name, (expected_sha, expected_rows) in TABLES.items():
        source = ROOT / name
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        assert digest == expected_sha, (name, "source SHA drift", digest)
        with source.open(encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == expected_rows, (name, len(rows), expected_rows)
        data[name] = rows
    main = {r["analysis"]:r for r in data["HNS_MODEL_SUMMARY.csv"]}
    diag = {r["analysis"]:r for r in data["HNS_GEOGRAPHY_AND_TEMPORAL.csv"]}
    assert len(main)==15 and len(diag)==5, "Duplicate HNS model name"
    near(main["primary_record"]["W_beta"],0.1300624446205031)
    near(diag["geography_record"]["W_beta"],0.10802204501635151)
    near(diag["geography_incidence"]["W_beta"],0.03417579058128888)
    near(diag["geography_record"]["H6_gain_pct"],5.374035658152443)
    assert float(diag["geography_record"]["H6_low"]) < 0 < float(diag["geography_record"]["H6_high"])
    assert float(diag["geography_incidence"]["beta_low"]) < 0 < float(diag["geography_incidence"]["beta_high"])
    assert (int(main["primary_record"]["n"]),int(main["primary_record"]["H7"]),int(main["primary_record"]["H6"]))==(454,232,90)
    assert (int(diag["temporal_record"]["n"]),int(diag["temporal_record"]["H7"]),int(diag["temporal_record"]["H6"]))==(222,222,89)
    families = {(r["cohort"],r["outcome"]):r for r in data["BUTTERFLY_FAMILY_NETWORK_PREDICTION.csv"]}
    assert len(families)==12, "Duplicate butterfly family outcome"
    pier = families[("coverage80","Pieridae_host_support")]
    sens = data["PIERIDAE_SOURCE_AND_SPECIES_SENSITIVITY.csv"]
    assert len({(r["test"],r["cohort"]) for r in sens})==len(sens)
    original = next(r for r in sens if r["test"]=="original" and r["cohort"]=="coverage80")
    near(pier["H6_gain_pct"],float(original["H6_gain_pct"]))
    assert (int(pier["n_rows"]),int(pier["H7"]),int(pier["H6"]))==(390,217,83)
    assert float(families[("coverage80","geographic_null_link_excess")]["H6_gain_pct"]) < 0
    country=next(r for r in sens if r["test"]=="species_composition_adjusted" and r["cohort"]=="coverage80")
    assert float(country["country_CI_low"]) < 0 < float(country["country_CI_high"])
    audit=json.loads((ROOT/"VERIFICATION_20261115.json").read_text())
    assert audit["module_manifest_entries"]["total"]==3508
    assert audit["unit_tests"]["total_pass"]==42 and audit["unit_tests"]["total_fail"]==0
    assert audit["sha_mismatches"]==0
    print(json.dumps({
      "status":"PASS", "git_csv_sha256_verified":4,
      "model_rows":{k:len(v) for k,v in data.items()},
      "local_tests_previously_run_not_reexecuted_by_this_workflow":42,
      "local_sha_entries_previously_checked_not_reread_by_workflow":3508,
      "causal_path_or_BACI_inference":False},sort_keys=True))
if __name__=="__main__":run()
