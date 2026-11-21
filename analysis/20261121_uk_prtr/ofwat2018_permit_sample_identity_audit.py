#!/usr/bin/env python3
"""Historical Ofwat sample/permit identity on the unchanged 930 source candidates.

Downloads and private source caches never enter the output directory. No TP,
flow values, biological records, source IDs, coordinates or names are emitted.
Ofwat's year summary is not a per-sample version assignment or an EEA crosswalk.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

from permit_sample_identity_support import (
    EA_TO_WQE, KEY_FIELDS, date_overlap_status, ofwat_pairs, permit_key,
    summarize_candidates, text,
)

ROOT = Path(__file__).resolve().parents[2]
SUCCESS = "HISTORICAL_OFWAT_SAMPLE_PERMIT_IDENTITIES_AUDITED_UWWTD_OUTLET_ID_UNRESOLVED"
OFWAT_URL = "https://www.ofwat.gov.uk/wp-content/uploads/2023/06/Water-company-phosphorus-discharge-data-2011-2021.xlsx"
OFWAT_ARCHIVE_URL = "https://webarchive.nationalarchives.gov.uk/ukgwa/20260203161134id_/" + OFWAT_URL
OFWAT_SHA = "b7699aa3202366b15baab4c415e0702a5fb8f0401249b415903b8e763939b3c4"
PERMIT_URL = "https://environment.data.gov.uk/api/file/download?fileDataSetId=a54fdea1-7769-4b22-a518-10d51fed6f33&fileName=Consented+Discharges+to+Controlled+Waters+with+Conditions.zip"
PERMIT_MEMBER = "Consented Discharges to Controlled Waters with Conditions.accdb"
PERMIT_SHA = "d664ea13da990c064006c6bfe32b3f76bc682c8b188bca55be33a3c8077136a7"
REGISTRY_SHA = "bf4cb7027196d8de54c68ea6ef1956f44d17e9909dd5f4b080de750a74ff328e"
COHORT_SHA = "d55a7773bd85c59bf137dda4dac835cd0a1f77d8912a1164a422ac48ee7cee07"
PAIR_SHA = "b7e73bf03a76a376a44245d6a4835f0cdcd6bb84142c7625d925f89aa3de9da1"
OFWAT_FIELDS = ("year", "sample_samplingpoint_notation", "EA_REGION", "permit_number",
                "version", "outlet_number", "effluent_number", "EFFECTIVE_DATE")
SHEETS = ("Sites with P permits", "Sites without P Permits")
FINAL_TYPE = "SEWAGE DISCHARGES - FINAL/TREATED EFFLUENT - WATER COMPANY"
STAGES = {"initialization", "ofwat_workbook", "permit_database", "historical_EEA",
          "current_WQE_registry", "fixed_candidate_cohort", "identity_aggregation", "serialization"}


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_file(path, expected):
    actual = file_sha(path)
    if actual != expected:
        raise ValueError("Frozen official source bytes changed")
    return {"sha256": actual, "bytes": Path(path).stat().st_size}


def download(url, path, cap):
    req = urllib.request.Request(url, headers={"User-Agent": "gesui-historical-permit-sampling-source-audit"})
    total = 0
    with urllib.request.urlopen(req, timeout=60) as response, path.open("wb") as target:
        if response.status != 200:
            raise ValueError("Source HTTP status not 200")
        for chunk in iter(lambda: response.read(131072), b""):
            total += len(chunk)
            if total > cap:
                raise ValueError("Source exceeds declared byte cap")
            target.write(chunk)
    return {"url": url, "http": 200, "bytes": total, "sha256": file_sha(path)}


def project_sheet(source_rows):
    """Refresh column positions at repeated headers; retain identifiers only."""
    positions = None
    result = []
    header_rows = []
    for row_number, row in enumerate(source_rows, 1):
        labels = [text(v) for v in row]
        if "sample_samplingpoint_notation" in labels:
            if any(labels.count(k) != 1 for k in OFWAT_FIELDS):
                raise ValueError("Required source header missing or duplicated")
            positions = {k: labels.index(k) for k in OFWAT_FIELDS}
            header_rows.append(row_number)
            continue
        if positions is None:
            continue
        selected = {k: row[i] if i < len(row) else None for k, i in positions.items()}
        if any(v is not None for v in selected.values()):
            result.append(selected)
    if positions is None:
        raise ValueError("Source identity header not found")
    return result, {"header_count": len(header_rows), "header_row_numbers": header_rows,
                    "selected_fields": list(OFWAT_FIELDS), "projected_rows": len(result)}


def read_ofwat(path):
    from openpyxl import load_workbook
    book = load_workbook(path, read_only=True, data_only=True)
    rows, qa, sheet_sets = [], {}, {}
    try:
        for name in SHEETS:
            selected, schema = project_sheet(book[name].iter_rows(values_only=True))
            pairs, counts = ofwat_pairs(selected)
            qa[name] = {"schema": schema, "identity_counts": counts,
                        "year_row_counts": dict(sorted(collections.Counter(
                            text(r["year"]) for r in selected if text(r["year"]).isdigit()).items()))}
            sheet_sets[name] = pairs
            rows.extend(selected)
    finally:
        book.close()
    pairs, totals = ofwat_pairs(rows)
    raw_pairs = {(text(r["sample_samplingpoint_notation"]), permit_key(r, "ofwat", normalize_numbers=False))
                 for r in rows if text(r["year"]) in {"2018", "2018.0"}
                 and (text(r["sample_samplingpoint_notation"]), permit_key(r, "ofwat")) in pairs}
    qa["union"] = totals
    qa["identity_pairs_shared_by_both_sheets"] = len(sheet_sets[SHEETS[0]] & sheet_sets[SHEETS[1]])
    qa["sheets_are_not_independent_or_assumed_disjoint_P_permit_groups"] = True
    return pairs, raw_pairs, qa


def mdb_rows(path, table, executable):
    kind = "consents" if table.startswith("consents_") else "effluents"
    extra = (("EFFECTIVE_DATE", "REVOCATION_DATE", "EFF_TYPE_DESCRIPTION") if kind == "consents"
             else ("EFF_SAMPLE_POINT", "SPT_DESC"))
    fields = KEY_FIELDS[kind] + extra
    with subprocess.Popen([executable, "-T", "%Y-%m-%d", str(path), table],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as process:
        try:
            reader = csv.DictReader(io.TextIOWrapper(process.stdout, encoding="utf-8", errors="strict"))
            if not set(fields) <= set(reader.fieldnames or []):
                raise ValueError("Required Access source fields missing")
            for i, row in enumerate(reader):
                if i >= 1000000:
                    raise ValueError("Source table exceeds row bound")
                yield {k: row[k] for k in fields}
            if process.wait(timeout=30) != 0:
                raise ValueError("Access source export failed")
        except BaseException:
            process.kill()
            process.wait(timeout=10)
            raise


def consent_index(rows):
    groups = collections.defaultdict(set)
    raw_keys = collections.defaultdict(set)
    multiplicity = collections.Counter()
    count = collections.Counter()
    for row in rows:
        count["source_rows"] += 1
        key = permit_key(row, "consents")
        if key is None:
            count["invalid_key_rows"] += 1
            continue
        groups[key].add(tuple(text(row.get(f)) for f in
                              ("EFFECTIVE_DATE", "REVOCATION_DATE", "EFF_TYPE_DESCRIPTION")))
        raw_keys[key].add(permit_key(row, "consents", normalize_numbers=False))
        multiplicity[key] += 1
    statuses = {}
    for key, values in groups.items():
        if len(values) != 1:
            statuses[key] = "conflicting_required_source_fields"
            continue
        effective, revoked, kind = next(iter(values))
        if kind != FINAL_TYPE:
            statuses[key] = "not_water_company_final_treated"
        else:
            statuses[key] = date_overlap_status({"EFFECTIVE_DATE": effective, "REVOCATION_DATE": revoked})
    eligible = {k for k, status in statuses.items() if status == "overlaps_by_recorded_dates"}
    qa = {"source_rows": count["source_rows"], "invalid_key_rows": count["invalid_key_rows"],
          "distinct_canonical_keys": len(groups),
          "duplicate_key_groups": sum(n > 1 for n in multiplicity.values()),
          "duplicate_rows_beyond_first_key": sum(n - 1 for n in multiplicity.values()),
          "key_status_counts": dict(sorted(collections.Counter(statuses.values()).items())),
          "eligible_distinct_keys": len(eligible),
          "date_screen_only_not_complete_2018_register": True}
    return eligible, statuses, raw_keys, qa


def corroborate(pairs, raw_pairs, all_index, active_index, effluent_rows):
    eligible, statuses, raw_keys, _ = all_index
    active_eligible = active_index[0]
    groups = collections.defaultdict(set)
    effluent_raw = collections.defaultdict(set)
    total = 0
    for row in effluent_rows:
        total += 1
        key = permit_key(row, "effluents")
        if key not in eligible:
            continue
        groups[key].add((text(row.get("EFF_SAMPLE_POINT")), text(row.get("SPT_DESC"))))
        effluent_raw[key].add(permit_key(row, "effluents", normalize_numbers=False))
    corroborated, active, strict = set(), set(), set()
    fates = collections.Counter()
    for pair in pairs:
        point, key = pair
        if key not in eligible:
            fates[statuses.get(key, "no_exact_key_in_consents_all")] += 1
            continue
        values = groups.get(key, set())
        if not values:
            fates["no_exact_key_in_Effluents"] += 1
            continue
        if len(values) != 1:
            fates["conflicting_Effluents_point_or_type"] += 1
            continue
        local, kind = next(iter(values))
        if kind != FINAL_TYPE:
            fates["Effluents_not_water_company_final_treated"] += 1
            continue
        if not local:
            fates["Effluents_sampling_point_missing"] += 1
            continue
        if EA_TO_WQE[key[0]] + "-" + local != point:
            fates["Effluents_sampling_point_disagrees"] += 1
            continue
        corroborated.add(pair)
        if key in active_eligible:
            active.add(pair)
        if any((point, raw) in raw_pairs for raw in raw_keys[key] & effluent_raw[key]):
            strict.add(pair)
        fates["corroborated"] += 1
    qa = {"Effluents_source_rows": total, "eligible_keys_with_Effluents": len(groups),
          "eligible_keys_with_conflicting_effluent_fields": sum(len(v) != 1 for v in groups.values()),
          "ofwat_pair_fates": dict(sorted(fates.items())),
          "corroborated_pairs_without_numeric_key_normalization": len(strict),
          "additional_pairs_from_integer_padding_normalization": len(corroborated - strict),
          "same_original_5_keys_and_cross_publication_sample_notation_required": True,
          "publications_share_EA_WIMS_lineage_not_independent_measurements": True}
    return corroborated, active, qa


def original_tables(path):
    import eea_UWWTD2018_UK_outlet_source_gate as frozen
    tables, manifest = {}, {}
    with zipfile.ZipFile(path) as archive:
        for name, expected in frozen.EXPECT.items():
            data = archive.read(frozen.BASE + name)
            actual = hashlib.sha256(data).hexdigest()
            if actual != expected:
                raise ValueError("Frozen EEA scientific table changed")
            tables[name] = list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
            manifest[name] = {"sha256": actual, "bytes": len(data)}
    reporters = [r for r in tables["ReportPeriod.csv"] if r["repCode"].strip() == "UK5"
                 and r["repReportedPeriod"].strip() == "2018" and "2018" in r["repSituationAt"]]
    plants = [r for r in tables["UWWTPS.csv"] if r["repCode"].strip() == "UK5"]
    outlets = [r for r in tables["DischargePoints.csv"] if r["repCode"].strip() == "UK5"]
    if len(reporters) != 1 or len(plants) != 1799 or len(outlets) != 1809:
        raise ValueError("Frozen EEA reporting population changed")
    coordinates = {(r["repCode"].strip(), r["uwwCode"].strip(), r["dcpCode"].strip()):
                   (float(r["dcpLongitude"]), float(r["dcpLatitude"]))
                   for r in outlets if r["dcpWaterBodyType"].strip() == "FW"}
    if len(coordinates) != 1498:
        raise ValueError("Frozen EEA freshwater population changed")
    counts = collections.Counter(r["uwwCode"].strip() for r in outlets)
    qa = {"scientific_members": manifest, "harmonised_UK5_plants": len(plants),
          "harmonised_UK5_discharge_rows": len(outlets), "freshwater_outlets": len(coordinates),
          "repCode_UK5_not_all_rows_in_original_UK_submission": True,
          "annual_treated_volume_nonblank": sum(bool(r.get("uwwWasteWaterTreated", "").strip()) for r in plants),
          "annual_volume_method_nonblank": sum(bool(r.get("uwwMethodWasteWaterTreated", "").strip()) for r in plants)}
    return coordinates, counts, qa


def current_registry(cache):
    sys.path.insert(0, str(ROOT / "analysis" / "20261120_wessex_wwtp"))
    import ea_national_SA_registry_to_ChemPop_source_ceiling as registry
    points, receipts = {}, []
    region_counts = collections.Counter()
    for page in range(40):
        skip = page * 250
        if cache is None:
            members, receipt = registry.get_page(skip)
        else:
            data = (cache / ("page_" + str(skip) + ".json")).read_bytes()
            doc = json.loads(data)
            members = doc["member"]
            receipt = {"skip": skip, "totalItems": doc["totalItems"], "n_source_members": len(members),
                       "source_sha256": hashlib.sha256(data).hexdigest(), "n_bytes": len(data),
                       "input_mode": "private_cached_official_response"}
        if receipt["totalItems"] != 5722 or len(members) > 250:
            raise ValueError("WQE SA registry count or page size changed")
        receipts.append(receipt)
        for item in members:
            name = text(item.get("notation"))
            region = text((item.get("region") or {}).get("notation"))
            if (item.get("samplingPointType") or {}).get("notation") != "SA":
                raise ValueError("Non-SA row in source-filtered registry")
            if not name or name in points or not name.startswith(region + "-") or region not in EA_TO_WQE.values():
                raise ValueError("Registry point identity or region invalid")
            coord, _ = registry.as_bng((item.get("geometry") or {}).get("asWKT"))
            if coord is None:
                raise ValueError("Original SA geometry missing or invalid")
            points[name] = coord
            region_counts[region] += 1
        if len(members) < 250:
            break
    digest = hashlib.sha256("\n".join(sorted(points)).encode()).hexdigest()
    if len(points) != 5722 or digest != REGISTRY_SHA:
        raise ValueError("Complete frozen WQE SA station set changed or incomplete")
    return points, {"points": len(points), "station_set_sha256": digest,
                    "page_receipts": receipts, "region_counts": dict(sorted(region_counts.items()))}


def fixed_candidates(coordinates, counts, registry):
    import platform
    import pyproj
    import eea2018_outlet_EA_TP_2018_ChemPop_river_gate as spatial
    matching, qa = spatial.select_current_SA_2018_UWWTD_historic_outlets(coordinates, counts, registry)
    station_sha = hashlib.sha256("\n".join(sorted(matching)).encode()).hexdigest()
    canonical = [[s, list(matching[s]["outlet_id"]), list(matching[s]["outlet_coord"])]
                 for s in sorted(matching)]
    pair_sha = hashlib.sha256(json.dumps(canonical, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    native_pair_sha = hashlib.sha256(json.dumps(
        [[s, list(matching[s]["outlet_id"])] for s in sorted(matching)],
        ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
    if len(matching) != 930 or station_sha != COHORT_SHA or pair_sha != PAIR_SHA:
        raise ValueError("Existing 930 station or complete pairing fingerprint changed")
    return set(matching), {"candidate_points": len(matching), "station_set_sha256": station_sha,
                           "complete_pair_mapping_sha256": pair_sha, "both_existing_fingerprints_verified": True,
                           "native_pair_mapping_sha256": native_pair_sha,
                           "native_pair_hash_is_additional_not_a_replacement_gate": True,
                           "spatial_runtime": {"python": platform.python_version(),
                                               "pyproj": pyproj.__version__, "PROJ": pyproj.proj_version_str},
                           "geographical_match_remains_distinct_from_formal_UWWTD_ID_match": True}


def produce(args, out):
    import eea_UWWTD2018_UK_outlet_source_gate as frozen
    out["source_urls"] = {"ofwat_original": OFWAT_URL, "ofwat_UK_National_Archives": OFWAT_ARCHIVE_URL,
                          "ea_permits": PERMIT_URL, "eea_original": frozen.URL,
                          "wqe_registry": "https://environment.data.gov.uk/water-quality/data/sampling-point"}
    with tempfile.TemporaryDirectory(prefix="gesui_permit_identity_private_") as temp:
        folder = Path(temp)
        out["stage"] = "ofwat_workbook"
        workbook = args.ofwat or folder / "source.xlsx"
        receipts = []
        if args.ofwat is None:
            for url in (OFWAT_URL, OFWAT_ARCHIVE_URL):
                try:
                    receipts.append(download(url, workbook, 30 * 1024 * 1024))
                    break
                except urllib.error.HTTPError as exc:
                    receipts.append({"url": url, "http": exc.code})
                    if url == OFWAT_ARCHIVE_URL:
                        raise
        out["ofwat_source"] = {**verify_file(workbook, OFWAT_SHA), "retrieval_receipts": receipts,
                                "verified_private_cache_used": args.ofwat is not None}
        pairs, raw_pairs, out["national_ofwat_2018_identity"] = read_ofwat(workbook)
        out["stage"] = "permit_database"
        database = args.permit_db or folder / "source.accdb"
        if args.permit_db is None:
            archive = folder / "permits.zip"
            out["permit_archive_receipt"] = download(PERMIT_URL, archive, 80 * 1024 * 1024)
            with zipfile.ZipFile(archive) as z, database.open("wb") as target:
                if z.getinfo(PERMIT_MEMBER).file_size > 600 * 1024 * 1024:
                    raise ValueError("Permit database exceeds declared byte cap")
                with z.open(PERMIT_MEMBER) as source:
                    for chunk in iter(lambda: source.read(1024 * 1024), b""):
                        target.write(chunk)
        out["permit_database_source"] = {**verify_file(database, PERMIT_SHA),
                                         "verified_private_cache_used": args.permit_db is not None}
        all_index = consent_index(mdb_rows(database, "consents_all", args.mdb_export))
        active_index = consent_index(mdb_rows(database, "consents_active", args.mdb_export))
        out["consents_all"] = all_index[3]
        out["consents_active"] = active_index[3]
        corroborated, active, out["EA_corroboration"] = corroborate(
            pairs, raw_pairs, all_index, active_index, mdb_rows(database, "Effluents", args.mdb_export))
        out["national_support"] = summarize_candidates({p for p, _ in pairs}, pairs, corroborated, active)
        out["stage"] = "historical_EEA"
        eea_zip = args.eea_zip or folder / "eea.zip"
        if args.eea_zip is None:
            out["eea_archive_receipt"] = download(frozen.URL, eea_zip, frozen.MAX_BYTES)
        coordinates, counts, out["EEA_2018_source"] = original_tables(eea_zip)
        out["stage"] = "current_WQE_registry"
        registry, out["current_WQE_registry"] = current_registry(args.registry_dir)
        out["stage"] = "fixed_candidate_cohort"
        candidates, out["fixed_cohort"] = fixed_candidates(coordinates, counts, registry)
        out["stage"] = "identity_aggregation"
        out["fixed_930_support"] = summarize_candidates(candidates, pairs, corroborated, active)
        out["ofwat_2018_points_in_current_SA_registry"] = len({p for p, _ in pairs} & set(registry))
        out["status"] = SUCCESS
        out["stage"] = "serialization"


def guarded_run(path, producer):
    out = {"schema": "GESUI_UK2018_HISTORICAL_OFWAT_PERMIT_SAMPLE_ID_AUDIT_V1",
           "status": "SOURCE_IDENTITY_AUDIT_STOP", "stage": "initialization",
           "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "no_TP_flow_or_biology_values_used_in_identity_analysis": True,
           "no_raw_identifiers_or_coordinates_in_public_output": True,
           "UWWTD_native_outlet_to_permit_outlet_identity_still_unverified": True,
           "annual_summary_not_sample_date_specific_permit_assignment": True}
    code = 0
    try:
        producer(out)
        if out["status"] != SUCCESS:
            raise ValueError("Source producer did not complete")
        encoded = json.dumps(out, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    except Exception as exc:
        code = 1
        out = {"schema": "GESUI_UK2018_HISTORICAL_OFWAT_PERMIT_SAMPLE_ID_AUDIT_V1",
               "status": "SOURCE_IDENTITY_AUDIT_STOP", "error_type": type(exc).__name__,
               "error_stage": out.get("stage") if out.get("stage") in STAGES else "unknown_stage",
               "raw_exception_message_and_source_rows_suppressed": True}
        encoded = json.dumps(out, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(encoded, encoding="utf-8")
    return code


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ofwat", type=Path)
    p.add_argument("--permit-db", type=Path)
    p.add_argument("--eea-zip", type=Path)
    p.add_argument("--registry-dir", type=Path)
    p.add_argument("--mdb-export", default="mdb-export")
    args = p.parse_args()
    try:
        code = guarded_run(args.out, lambda out: produce(args, out))
        result = json.loads(args.out.read_text(encoding="utf-8"))
        print(json.dumps({k: result[k] for k in ("status", "error_type", "error_stage", "fixed_930_support")
                          if k in result}, sort_keys=True), flush=True)
        return code
    except Exception as exc:
        print(json.dumps({"status": "SOURCE_IDENTITY_OUTPUT_STOP", "error_type": type(exc).__name__}), flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
