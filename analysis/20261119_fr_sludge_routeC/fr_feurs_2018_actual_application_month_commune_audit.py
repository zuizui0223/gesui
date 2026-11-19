#!/usr/bin/env python3
"""Audit 43 anonymized *realized* 2018 Feurs sludge field-application rows.

The 2018 Feurs official agronomic report supplies independent validation totals:
Table 7: 43 records, 161.53 applied ha, 1,918 wet tonnes;
Table 6: 761t March, 50t April, 613t July, 214t September,
         208t October, 72t November.
Rows are recorded anonymously; NO farmer/holding/parcel identifiers or
biological outcomes are read. This is ONE station-year source cluster.
"""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

STEU = "0442094S0002"
SOURCE_PDF = "https://www.auvergne-rhone-alpes.developpement-durable.gouv.fr/IMG/pdf/10_step_feurs_cas_par-cas_annexe9_epandage_2018bilan.pdf"
EXPECTED_COLUMNS = (
    "INSEE_code", "year", "month", "applied_area_ha", "wet_sludge_t",
    "source_printed_page",
)
EXPECTED_CODES = {"42066", "42174", "42251", "42155", "42319", "42094", "42143", "42135"}
EXPECTED_MONTHLY_WET_T = {3: 761, 4: 50, 7: 613, 9: 214, 10: 208, 11: 72}
EXPECTED_APPLICATION_ROWS = 43
EXPECTED_TOTAL_AREA_HA = Decimal("161.53")
EXPECTED_TOTAL_WET_T = Decimal("1918")


def summary(rows):
    if len(rows) != EXPECTED_APPLICATION_ROWS:
        raise ValueError("Table 7 count no longer equals 43")
    by_month, by_commune, by_commune_month = defaultdict(list), defaultdict(list), defaultdict(list)
    accepted = []
    for r in rows:
        if set(r) != set(EXPECTED_COLUMNS):
            raise ValueError("Unexpected columns (possibly private information)")
        code = r["INSEE_code"]
        if code not in EXPECTED_CODES or r["year"] != "2018":
            raise ValueError("Invalid 2018 source municipality/year")
        try:
            month = int(r["month"])
            area = Decimal(r["applied_area_ha"])
            wet = Decimal(r["wet_sludge_t"])
            page = int(r["source_printed_page"])
        except (TypeError, ValueError, InvalidOperation) as exc:
            raise ValueError("Malformed application record") from exc
        if month not in EXPECTED_MONTHLY_WET_T or page not in (16, 17):
            raise ValueError("Invalid original Table 7 month/page")
        if not area.is_finite() or not wet.is_finite() or area <= 0 or wet <= 0:
            raise ValueError("Missing, negative or nonfinite applied quantity")
        # Identical area/tonnes/month values may be legitimate distinct log entries.
        accepted.append((code, month, area, wet))
        by_month[month].append((area, wet))
        by_commune[code].append((area, wet))
        by_commune_month[(code, month)].append((area, wet))
    if set(by_commune) != EXPECTED_CODES:
        raise ValueError("Missing/unexpected actual application commune")
    area_total = sum((x[2] for x in accepted), Decimal("0"))
    wet_total = sum((x[3] for x in accepted), Decimal("0"))
    if area_total != EXPECTED_TOTAL_AREA_HA or wet_total != EXPECTED_TOTAL_WET_T:
        raise ValueError("Original Table 7 reported annual total mismatch")
    for month, total in EXPECTED_MONTHLY_WET_T.items():
        if sum((q[1] for q in by_month[month]), Decimal("0")) != total:
            raise ValueError("Original Table 6 independent monthly total mismatch")
    def aggregate(group):
        return {
            "application_log_rows": len(group),
            "actual_applied_area_ha": float(sum((x[0] for x in group), Decimal("0"))),
            "actual_wet_sludge_t": float(sum((x[1] for x in group), Decimal("0"))),
        }
    month_summary = {str(k): aggregate(v) for k, v in sorted(by_month.items())}
    commune_summary = {str(k): aggregate(v) for k, v in sorted(by_commune.items())}
    commune_month_summary = {f"{code}|{month:02d}": aggregate(v) for (code, month), v in sorted(by_commune_month.items())}
    march_july = month_summary["3"]["actual_wet_sludge_t"] + month_summary["7"]["actual_wet_sludge_t"]
    return {
        "schema": "gesui_FR_Fe u rs_2018_month_commune_real_application_source_v1".replace(" ", ""),
        "station_SANDRE": STEU,
        "application_year": 2018,
        "original_source_pdf": SOURCE_PDF,
        "source_locator": "2018 Feurs bilan agronomique Tables 6–7, printed pp.15–17",
        "primary_pdf_binary_hash_not_available": True,
        "source_table7_actual_application_rows": len(rows),
        "municipalities_with_actual_application": len(by_commune),
        "original_annual_applied_area_ha": float(area_total),
        "original_annual_applied_wet_t": float(wet_total),
        "month": month_summary,
        "commune": commune_summary,
        "commune_month_source_opportunity_cells": len(by_commune_month),
        "commune_month": commune_month_summary,
        "commune_month_not_confirmed_specific_H7_field_exposure": True,
        "march_july_wet_t": march_july,
        "march_july_share_of_2018_wet_t": round(march_july / float(wet_total), 6),
        "independent_realized_stations_years_in_current_case": 1,
        "log_rows_and_communes_not_independent_stations_or_H6": True,
        "actual_parcel_coordinates_not_present_in_anonymized_extract": True,
        "actual_sludge_N_P_assays_not_assigned_to_field_or_month": True,
        "not_evidence_of_plant_or_butterfly_N_uptake": True,
        "same_2018_aqueous_measured_TN_TP_export_not_identified": True,
        "biological_analysis_gate": "STOP",
        "no_NUTRI_HNS_EQR_values_read": True,
    }


def audit_csv(path):
    p = Path(path)
    with p.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != EXPECTED_COLUMNS:
            raise ValueError("CSV has unexpected column order or privacy fields")
        rows = list(reader)
    result = summary(rows)
    result["anonymized_source_csv_sha256"] = hashlib.sha256(p.read_bytes()).hexdigest()
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    result = audit_csv(args.input)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("FEURS_2018_ANON_REAL_SPREADING_MONTH_COMMUNE", json.dumps(result, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
