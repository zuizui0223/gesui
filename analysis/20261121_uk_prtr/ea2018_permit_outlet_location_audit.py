"""Audit recorded permit outlet locations on the unchanged historical cohort.

Private source rows never enter output. Version-linked current-publication
geometry does not prove a 2018 survey or a formal EEA/native outlet identity.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import urllib.error
import zipfile
from pathlib import Path

import ofwat2018_permit_sample_identity_audit as identity
from permit_outlet_location_support import FIELDS, location_index, summarize_locations
from permit_sample_identity_support import KEY_FIELDS, summarize_candidates

SUCCESS = "HISTORICAL_PERMIT_OUTLET_LOCATION_ATTRIBUTES_AUDITED_EEA_NATIVE_ID_UNRESOLVED"
SCHEMA = "GESUI_UK2018_PERMIT_OUTLET_LOCATION_AUDIT_V1"
BASELINE_PATH = identity.ROOT / "results/uk2018_permit_identity_audit_20261121/run_44907990658.json"
BASELINE_SHA = "40e318cffc0d2abb507df73f2933ce4a239dd71d56bc200a2f209e7cbba2bc48"
# Reconstructed from the accepted identity implementation and frozen source bytes.
# The earlier result did not record this additional station/five-key fingerprint.
CORROBORATED_PAIR_SHA = "1ce2b5be966498d282a5f3ae23bff896fa27c9b63712d67baba3ea7881f308ae"
STAGES = {"initialization", "ofwat_workbook", "permit_database", "historical_EEA",
          "current_WQE_registry", "fixed_candidate_cohort", "outlet_locations", "serialization"}
BASELINE = {
    "candidate_points": 930, "with_historical_ofwat_pair": 346,
    "with_exactly_one_ofwat_key": 346, "with_multiple_ofwat_keys": 0,
    "with_ea_all_corroborated_pair": 295, "with_ea_active_corroborated_pair": 16,
    "with_corroboration_only_in_all_not_active": 279, "without_ofwat_pair": 584,
}


def verify_registry_pages(registry_qa, baseline_path=BASELINE_PATH):
    """Pin coordinate-bearing source bytes, not just station membership."""
    identity.verify_file(baseline_path, BASELINE_SHA)
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))["current_WQE_registry"]

    def complete_pages(qa):
        if qa.get("points") != 5722 or qa.get("station_set_sha256") != identity.REGISTRY_SHA:
            raise ValueError("Accepted registry identity changed")
        receipts = qa.get("page_receipts", [])
        if len(receipts) != 23:
            raise ValueError("Complete accepted registry page coverage missing")
        pages = {}
        for receipt in receipts:
            skip = receipt["skip"]
            if skip in pages or skip not in range(0, 5722, 250):
                raise ValueError("Registry page duplicated or outside accepted coverage")
            if (receipt["totalItems"] != 5722
                    or receipt["n_source_members"] != min(250, 5722 - skip)):
                raise ValueError("Registry page size changed")
            pages[skip] = tuple(receipt[field] for field in
                                ("source_sha256", "n_bytes", "n_source_members", "totalItems"))
        if set(pages) != set(range(0, 5722, 250)):
            raise ValueError("Accepted registry page missing")
        return pages

    if complete_pages(registry_qa) != complete_pages(baseline):
        raise ValueError("Accepted coordinate-bearing registry source bytes changed")
    return {"accepted_identity_run": 44907990658, "accepted_result_sha256": BASELINE_SHA,
            "verified_pages": 23, "all_original_coordinate_bearing_page_bytes_match": True}


def audit_fixed_candidates(candidates, pairs, corroborated, active, eligible_keys,
                           location_rows, matching, registry, *,
                           expected_pair_sha256=CORROBORATED_PAIR_SHA):
    support = summarize_candidates(candidates, pairs, corroborated, active)
    if any(support.get(key) != value for key, value in BASELINE.items()):
        raise ValueError("Accepted historical identity support changed")
    if set(matching) != set(candidates):
        raise ValueError("Complete original candidate mapping missing")
    subset = {(point, key) for point, key in corroborated if point in candidates}
    if len(subset) != 295:
        raise ValueError("Accepted candidate pair support changed")
    canonical = sorted([[point, list(key)] for point, key in subset])
    pair_sha = hashlib.sha256(json.dumps(canonical, ensure_ascii=False,
                                        separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()
    if pair_sha != expected_pair_sha256:
        raise ValueError("Accepted station and full permit key identities changed")
    index, qa = location_index(location_rows, eligible_keys)
    if any(statuses.get("absent_key", 0) for statuses in qa["field_status_counts"].values()):
        raise ValueError("Location export omitted a required historical key")
    return {
        "prior_identity_partition": support,
        "fixed_295_native_permit_pair_sha256": pair_sha,
        "national_location_QA": qa,
        "fixed_295_location_support": summarize_locations(subset, index, matching, registry),
    }


def location_csv_rows(stream):
    """Reject malformed exports; explicit empty fields remain source missingness."""
    fields = KEY_FIELDS["consents"] + FIELDS
    reader = csv.reader(stream, strict=True)
    header = next(reader, [])
    if any(header.count(field) != 1 for field in fields):
        raise ValueError("Required original location column missing or duplicated")
    positions = {field: header.index(field) for field in fields}
    for number, row in enumerate(reader):
        if number >= 1000000:
            raise ValueError("Location source exceeds row bound")
        if len(row) != len(header):
            raise ValueError("Location export record width differs from header")
        yield {field: row[position] for field, position in positions.items()}


def mdb_location_rows(path, executable):
    with subprocess.Popen([executable, "-T", "%Y-%m-%d", str(path), "consents_all"],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as process:
        try:
            yield from location_csv_rows(io.TextIOWrapper(process.stdout, encoding="utf-8", errors="strict"))
            if process.wait(timeout=30) != 0:
                raise ValueError("Location table export failed")
        except BaseException:
            process.kill()
            process.wait(timeout=10)
            raise


def produce(args, out):
    import eea_UWWTD2018_UK_outlet_source_gate as eea
    import eea2018_outlet_EA_TP_2018_ChemPop_river_gate as spatial
    with tempfile.TemporaryDirectory(prefix="gesui_permit_outlet_private_") as folder:
        temp = Path(folder)
        out["stage"] = "ofwat_workbook"
        workbook = args.ofwat or temp / "source.xlsx"
        receipts = []
        if args.ofwat is None:
            for url in (identity.OFWAT_URL, identity.OFWAT_ARCHIVE_URL):
                try:
                    receipts.append(identity.download(url, workbook, 30 * 1024 * 1024))
                    break
                except urllib.error.HTTPError as exc:
                    receipts.append({"url": url, "http": exc.code})
                    if url == identity.OFWAT_ARCHIVE_URL:
                        raise
        out["ofwat_source"] = {**identity.verify_file(workbook, identity.OFWAT_SHA),
                                "retrieval_receipts": receipts,
                                "verified_private_cache_used": args.ofwat is not None}
        pairs, raw_pairs, out["ofwat_identity_QA"] = identity.read_ofwat(workbook)
        out["stage"] = "permit_database"
        database = args.permit_db or temp / "source.accdb"
        if args.permit_db is None:
            archive = temp / "permits.zip"
            out["permit_archive_receipt"] = identity.download(identity.PERMIT_URL, archive, 80 * 1024 * 1024)
            with zipfile.ZipFile(archive) as z, database.open("wb") as target:
                if z.getinfo(identity.PERMIT_MEMBER).file_size > 600 * 1024 * 1024:
                    raise ValueError("Permit database exceeds declared byte cap")
                with z.open(identity.PERMIT_MEMBER) as source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(chunk)
        out["permit_database_source"] = {**identity.verify_file(database, identity.PERMIT_SHA),
                                         "verified_private_cache_used": args.permit_db is not None}
        all_index = identity.consent_index(identity.mdb_rows(database, "consents_all", args.mdb_export))
        active_index = identity.consent_index(identity.mdb_rows(database, "consents_active", args.mdb_export))
        out["consents_all_QA"], out["consents_active_QA"] = all_index[3], active_index[3]
        corroborated, active, out["EA_corroboration_QA"] = identity.corroborate(
            pairs, raw_pairs, all_index, active_index, identity.mdb_rows(database, "Effluents", args.mdb_export))
        out["stage"] = "historical_EEA"
        eea_archive = args.eea_zip or temp / "eea.zip"
        if args.eea_zip is None:
            out["EEA_archive_receipt"] = identity.download(eea.URL, eea_archive, eea.MAX_BYTES)
        coordinates, counts, out["EEA_2018_source"] = identity.original_tables(eea_archive)
        out["stage"] = "current_WQE_registry"
        registry, out["WQE_registry"] = identity.current_registry(args.registry_dir)
        out["WQE_coordinate_source_freeze"] = verify_registry_pages(out["WQE_registry"])
        out["stage"] = "fixed_candidate_cohort"
        candidates, out["fixed_cohort"] = identity.fixed_candidates(coordinates, counts, registry)
        matching, _ = spatial.select_current_SA_2018_UWWTD_historic_outlets(coordinates, counts, registry)
        out["stage"] = "outlet_locations"
        out.update(audit_fixed_candidates(candidates, pairs, corroborated, active, all_index[0],
                                         mdb_location_rows(database, args.mdb_export), matching, registry))
        if out["national_location_QA"]["source_rows"] != out["consents_all_QA"]["source_rows"]:
            raise ValueError("Location extraction did not reproduce full source row count")
        out["status"] = SUCCESS
        out["stage"] = "serialization"


def guarded_run(path, producer):
    out = {"schema": SCHEMA, "status": "PERMIT_OUTLET_LOCATION_AUDIT_STOP",
           "stage": "initialization", "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "no_TP_flow_BIO_values_or_river_rerouting_used": True,
           "no_raw_identifiers_coordinates_names_or_individual_distances_published": True,
           "formal_EEA_outlet_to_permit_identity_still_unverified": True,
           "permit_dates_are_not_coordinate_survey_dates": True}
    exit_code = 0
    try:
        producer(out)
        if out["status"] != SUCCESS:
            raise ValueError("Source producer did not complete")
        encoded = json.dumps(out, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n"
    except Exception as exc:
        exit_code = 1
        out = {"schema": SCHEMA, "status": "PERMIT_OUTLET_LOCATION_AUDIT_STOP",
               "error_type": type(exc).__name__,
               "error_stage": out.get("stage") if out.get("stage") in STAGES else "unknown_stage",
               "raw_exception_message_and_partial_results_suppressed": True}
        encoded = json.dumps(out, indent=2, sort_keys=True) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(encoded, encoding="utf-8")
    return exit_code


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--ofwat", type=Path)
    parser.add_argument("--permit-db", type=Path)
    parser.add_argument("--eea-zip", type=Path)
    parser.add_argument("--registry-dir", type=Path)
    parser.add_argument("--mdb-export", default="mdb-export")
    args = parser.parse_args()
    try:
        code = guarded_run(args.out, lambda out: produce(args, out))
        result = json.loads(args.out.read_text(encoding="utf-8"))
        print(json.dumps({key: result[key] for key in ("status", "error_type", "error_stage", "fixed_295_location_support")
                          if key in result}, sort_keys=True), flush=True)
        return code
    except Exception as exc:
        print(json.dumps({"status": "PERMIT_OUTLET_LOCATION_OUTPUT_STOP", "error_type": type(exc).__name__}), flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
