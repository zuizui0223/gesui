#!/usr/bin/env python3
"""Repair ONLY original-source GBIF speciesKey HTTP429 failures in the frozen
234-H7 × two periods × six taxon-group cohort. Never zero-impute missing taxa.
Archive all first-run successes unchanged. Source repair is not ecology."""
import collections
import hashlib
import json
import os
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from gesui_H7_original_fivefamily_plants_speciesKey_retention import (
    SOURCE_SHA, PERIODS, FAMILIES, read_frozen, task,
)

EXPECTED_FAILED = {
    ("2070020830","2020-2022","Papilionidae"),
    ("2070026940","2017-2019","Plantae"),
    ("2070290810","2020-2022","Plantae"),
    ("2070493320","2020-2022","Papilionidae"),
    ("2070499660","2020-2022","Hesperiidae"),
    ("2070502920","2020-2022","Plantae"),
    ("2070016510","2020-2022","Plantae"),
    ("2070548700","2020-2022","Papilionidae"),
}
FIELDS = ("Plantae",) + FAMILIES
def identity(r):
    return (str(r["HYBAS_ID_L07"]), r["period"], r["group"])

def main():
    basins, families = read_frozen()
    basin_by_id = {str(z["HYBAS_ID_L07"]): z for z in basins}
    assert len(basin_by_id) == 234
    original = {}
    shard_audit = []
    for i in range(4):
        folder = pathlib.Path(os.environ["SPECIESKEY_SHARD_DIR_"+str(i)])
        rows_path = folder / ("H7_speciesKey_full_shard"+str(i)+".json")
        audit_path = folder / ("H7_speciesKey_full_shard"+str(i)+"_audit.json")
        audit = json.loads(audit_path.read_text())
        if audit["source_Frozen_234_source_SHA"] != SOURCE_SHA or audit["shard"] != i or audit["mode"] != "full":
            raise RuntimeError("Wrong authentic original source shard or SHA")
        rows = json.loads(rows_path.read_text())
        if len(rows) != audit["n_target_source_facets"] or len(rows) != audit["n_original_H7"]*12:
            raise RuntimeError("Incomplete or row-inconsistent original source artifact")
        shard_audit.append(audit)
        for r in rows:
            k = identity(r)
            if k in original or k[0] not in basin_by_id or k[1] not in PERIODS or k[2] not in FIELDS:
                raise RuntimeError("Duplicate/invalid H7×period×group source identity")
            original[k] = r
    all_required = {(h,p,g) for h in basin_by_id for p in PERIODS for g in FIELDS}
    if set(original) != all_required or len(original) != 2808:
        raise RuntimeError("Source original 2808 opportunities not complete")
    missing = {k for k,v in original.items() if v["source_status"] != "ok"}
    if missing != EXPECTED_FAILED:
        raise RuntimeError("Observed first-run source failure set does not match signed 8-case 429 source gate: "+repr(sorted(missing)))
    for k in missing:
        q = original[k]
        errors=[e for receipt in q["source_request_receipts"] for e in receipt.get("errors",[])]
        if not errors or not all(str(e).startswith("429:") for e in errors):
            raise RuntimeError("Found a non-HTTP429 source failure, do not silently retry")
    out = pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_original_source_retry429"
    out.mkdir(parents=True,exist_ok=True)
    attempts=[]
    # Serial source requests, cool-down after the 4-shard 12-worker initial run.
    time.sleep(45)
    for j,k in enumerate(sorted(missing)):
        if j: time.sleep(25)
        h,p,g=k
        new_row=task(basin_by_id[h],p,g,families)
        ok=new_row["source_status"]=="ok"
        attempts.append({"identity":list(k),"status":"recovered" if ok else "still_missing",
                         "source_query_count":new_row.get("source_query_occurrence_count"),
                         "speciesKey_count":new_row.get("source_speciesKey_count")})
        if ok:
            original[k]=new_row
        print(json.dumps(attempts[-1],sort_keys=True),flush=True)
    allrows=[original[k] for k in sorted(original)]
    if len(allrows)!=2808: raise RuntimeError("Merged original source row drift")
    remain=[identity(r) for r in allrows if r["source_status"]!="ok"]
    keybins=collections.defaultdict(dict)
    for r in allrows:
        keybins[(str(r["HYBAS_ID_L07"]),r["period"])][r["group"]]=r
    diag={}
    for threshold in (0.8,0.9,0.95):
        eligible=0
        for d in keybins.values():
            if len(d)!=6:raise RuntimeError("Original group identity drift")
            if all(z["source_status"]=="ok" and
                   z["source_facet_coverage"] is not None and
                   z["source_facet_coverage"]>=threshold for z in d.values()):
                eligible+=1
        diag[str(threshold)]=eligible
    audit={
        "origin_source_SHA":SOURCE_SHA,
        "initial_run":43592211212,
        "expected_unique_H7":234,
        "expected_H7_period":468,
        "expected_taxon_group_facets":2808,
        "initial_success":2800,
        "initial_HTTP429_missing":8,
        "retry_success":8-len(remain),
        "final_missing":len(remain),
        "remaining_original_H7_period_group":[list(k) for k in sorted(remain)],
        "all_six_species_facetted_coverage_threshold_counts":diag,
        "source_request_sha_retained_per_query":True,
        "new_source_requery_differs_in_acquisition_time":True,
        "no_zero_imputation_or_ecological_inference":True,
        "shard_first_run_audits":shard_audit,
        "attempts":attempts,
    }
    (out/"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json").write_text(json.dumps(allrows,ensure_ascii=False,sort_keys=True))
    (out/"GESUI_H7_234x2x6_SPECIESKEY_429_RECOVERY_AUDIT.json").write_text(json.dumps(audit,indent=2,sort_keys=True))
    print(json.dumps({k:v for k,v in audit.items() if k not in ("shard_first_run_audits","attempts")},sort_keys=True),flush=True)
    if remain:
        raise SystemExit("Original source HTTP429 unresolved; these values stay missing, not biological zero")

if __name__=="__main__":
    main()
