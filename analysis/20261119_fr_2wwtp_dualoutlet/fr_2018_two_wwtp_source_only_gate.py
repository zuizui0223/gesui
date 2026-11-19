#!/usr/bin/env python3
"""SOURCE-only 2018 dual-outlet evidence from two independent French WWTPs.

Use original document's hydrologic vs agronomic measurements as separate
variables. No water annual N/P masses are shown in source documents, and
one treatment-water outflow value is NOT receiving-river total discharge if
bypass telemetry is unreliable. No plant, butterfly, or aquatic response data.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

MAIN_COLS=(
    "station_SANDRE","year","station_name","receiving_river",
    "land_realized_bilan_proven","land_application_communes",
    "land_application_log_rows","land_actual_area_ha","land_actual_wet_t",
    "land_actual_dry_t","land_dry_matter_lime_status",
    "water_A4_treated_outlet_m3","water_A4_flow_source_status",
    "water_bypass_volume_reliability","annual_aqueous_total_N_measured_t",
    "annual_aqueous_total_P_measured_t","annual_land_total_N_measured_t",
    "annual_land_total_elemental_P_measured_t","actual_field_geometry_ready",
    "land_application_months","water_source_uri","land_source_uri",
    "source_pdf_binary_sha_verified"
)
EVENT_COLS=(
    "station_SANDRE","year","month","commune",
    "register_event_rows","realized_mass_by_group_unidentified"
)
EXPECTED_STATIONS={
    "0442094S0002":{"name":"Feurs","water":"Loire","communes":8,
      "events":43,"area":Decimal("161.53"),"wet":Decimal("1918"),
      "dry":Decimal("481.31"),"lime":"lime_included",
      "months":("03","04","07","09","10","11"),
      "water_A4":None,"water_flow_status":"permit_outlet_only",
      "bypass":"not_audited"},
    "0517347V004":{"name":"Saint-Jean-d'Angely-Moulinveau-La-Vergne",
      "water":"Boutonne","communes":2,"events":8,
      "area":Decimal("74.20"),"wet":Decimal("845"),
      "dry":Decimal("162.45"),"lime":"unlimed",
      "months":("03","09"),"water_A4":343221,
      "water_flow_status":"reported_A4_treated_outlet",
      "bypass":"unreliable_A5_telemetry"},
}
SJ_EVENTS={("03","Ternant"):2,("03","Mazeray"):3,
           ("09","Ternant"):2,("09","Mazeray"):1}


def decimal_exact(x):
    try:
        r=Decimal(x)
    except (InvalidOperation,TypeError,ValueError) as exc:
        raise ValueError("Missing or malformed positive original mass") from exc
    if not r.is_finite() or r<=0:
        raise ValueError("Negative, zero or nonfinite realized applications")
    return r


def getcsv(path, header):
    with Path(path).open(newline="",encoding="utf-8-sig") as h:
        reader=csv.DictReader(h)
        if tuple(reader.fieldnames or ())!=header:
            raise ValueError("Unexpected SOURCE columns (could leak PII or confuse units)")
        return list(reader)


def audit(source_rows,event_rows):
    if len(source_rows)!=2:raise ValueError("Two independent 2018 WWTP records required")
    seen=set()
    out=[]
    for r in source_rows:
        if tuple(r)!=MAIN_COLS:
            raise ValueError("Unexpected SOURCE row schema")
        s=r["station_SANDRE"]
        if s not in EXPECTED_STATIONS or s in seen or r["year"]!="2018":
            raise ValueError("Invalid or duplicated ST EU-year pairing")
        seen.add(s)
        exp=EXPECTED_STATIONS[s]
        if r["station_name"]!=exp["name"] or r["receiving_river"]!=exp["water"]:
            raise ValueError("Original source ID or river changed")
        if r["land_realized_bilan_proven"]!="true" or r["actual_field_geometry_ready"]!="false":
            raise ValueError("Cannot infer actual field polygons from agronomy table")
        if r["source_pdf_binary_sha_verified"]!="false":
            raise ValueError("No independently hashed primary binary for this source")
        if any(r[x] for x in ("annual_aqueous_total_N_measured_t",
                "annual_aqueous_total_P_measured_t","annual_land_total_N_measured_t",
                "annual_land_total_elemental_P_measured_t")):
            raise ValueError("Do not manufacture measured water or land N/P from source totals")
        if any(not r[k].startswith("https://") for k in ("water_source_uri","land_source_uri")):
            raise ValueError("Untraced source URL")
        if (r["land_dry_matter_lime_status"]!=exp["lime"] or
            r["water_A4_flow_source_status"]!=exp["water_flow_status"] or
            r["water_bypass_volume_reliability"]!=exp["bypass"]):
            raise ValueError("Lime or water bypass classification changed")
        area=decimal_exact(r["land_actual_area_ha"])
        wet=decimal_exact(r["land_actual_wet_t"])
        dry=decimal_exact(r["land_actual_dry_t"])
        if (area,wet,dry)!=(exp["area"],exp["wet"],exp["dry"]):
            raise ValueError("Original actual dry/wet land amount changed")
        if (int(r["land_application_communes"])!=exp["communes"] or
            int(r["land_application_log_rows"])!=exp["events"]):
            raise ValueError("Original field application events or commune count changed")
        months=tuple(r["land_application_months"].split(";"))
        if months!=exp["months"]:
            raise ValueError("Do not substitute 2019/2021 or planned sludge years")
        a4=r["water_A4_treated_outlet_m3"]
        if (None if a4=="" else int(a4))!=exp["water_A4"]:
            raise ValueError("Outlet A4 must be observed as reported, not guessed")
        out.append({
            "station_SANDRE":s,"year":2018,"receiving_river":exp["water"],
            "realized_sludge_area_ha":float(area),
            "realized_wet_sludge_t":float(wet),
            "realized_dry_sludge_t_NOT_CROSS_PLANT_COMPARABLE":float(dry),
            "sludge_dry_matter_lime_status":exp["lime"],
            "actual_application_log_rows":exp["events"],
            "actual_application_communes":exp["communes"],
            "realized_application_months":list(months),
            "water_treated_outlet_A4_m3":exp["water_A4"],
            "water_A4_volume_is_not_confirmed_total_received_with_unreliable_A5":s=="0517347V004",
            "original_water_report_has_outlet_NGL_PT_monitoring_not_annual_kg":s=="0517347V004",
            "same_year_measured_annual_aqueous_N_and_P_export_available":False,
            "same_year_field_specific_measured_N_P2O5_mass_available":False,
            "receiving_field_geometries_not_identified":True,
            "no_biological_values_read":True
        })
    if seen!=set(EXPECTED_STATIONS):
        raise ValueError("Expected original two stations not both identified")
    if len(event_rows)!=4:
        raise ValueError("Saint-Jean anonymized application ledger groups expected 4")
    pairs=set()
    monthcounts=Counter()
    communes=set()
    for row in event_rows:
        if tuple(row)!=EVENT_COLS or row["station_SANDRE"]!="0517347V004" or row["year"]!="2018":
            raise ValueError("Wrong station, year or possible PII in ledger")
        key=(row["month"],row["commune"])
        if key in pairs or key not in SJ_EVENTS:raise ValueError("Duplicate or unexpected ledger group")
        pairs.add(key)
        if int(row["register_event_rows"])!=SJ_EVENTS[key]:
            raise ValueError("Original 2018 ledger event count no longer reconciles")
        if row["realized_mass_by_group_unidentified"]!="true":
            raise ValueError("Event-level sludge masses must remain unspecified")
        monthcounts[row["month"]]+=int(row["register_event_rows"])
        communes.add(row["commune"])
    if pairs!=set(SJ_EVENTS) or monthcounts!=Counter({"03":5,"09":3}) or len(communes)!=2:
        raise ValueError("2018 agronomy ledger month/commune mismatch")
    # No independent H6, land nutrient mass or water nutrient source in this cohort.
    return {
        "schema":"gesui_FR2018_2independent_WWTP_realized_spreading_water_outlets_SOURCE_only_v1",
        "source_year":2018,
        "independent_identified_STEU_year_units":2,
        "distinct_receiving_rivers_from_original_sources":2,
        "independent_H6_assigned_and_audited":0,
        "2018_realized_plant_application_ledger_rows":sum(x["events"] for x in EXPECTED_STATIONS.values()),
        "land_actual_area_ha_sum_TWO_SOURCES_DESCRIPTIVE_ONLY":float(sum((x["area"] for x in EXPECTED_STATIONS.values()),Decimal(0))),
        "land_actual_wet_t_sum_TWO_SOURCES_DESCRIPTIVE_ONLY":float(sum((x["wet"] for x in EXPECTED_STATIONS.values()),Decimal(0))),
        "sludge_total_dry_mass_across_STEUs_liming_mismatch_not_pooled":True,
        "one_plant_same_year_A4_treated_outlet_m3_reported":1,
        "Saint_Jean_2018_A4_m3":343221,
        "Saint_Jean_2018_by_pass_A5_instrumentation_unreliable":True,
        "Saint_Jean_2018_influent_NGL_and_PT_loads_must_not_be_treated_as_outflow":True,
        "Saint_Jean_2018_aqueous_A4_NGL_PT_monitoring_not_annual_export":True,
        "Saint_Jean_month_register_rows":dict(sorted(monthcounts.items())),
        "actual_field_geometries_verified_and_linked":0,
        "2018_plant_years_with_both_measured_aqueous_TN_TP_and_realized_land_N_P":0,
        "NUTRI_HNS_EQR_values_read":False,
        "realized_field_and_water_source_pairing":"2018 SAME_ORIGIN_OUTLETS; NO NUTRIENT_FLUX_PAIR",
        "coupled_N_P_mass_balance_gate":"STOP",
        "biotic_response_and_causal_path_gate":"STOP",
        "source_rows":out,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--plants",type=Path,required=True)
    p.add_argument("--events",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    source=getcsv(a.plants,MAIN_COLS)
    events=getcsv(a.events,EVENT_COLS)
    r=audit(source,events)
    r["source_receipts_SHA256"]={
        "2018_plant_source_rows":hashlib.sha256(a.plants.read_bytes()).hexdigest(),
        "2018_SJ_month_commune_counts":hashlib.sha256(a.events.read_bytes()).hexdigest(),
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("GESUI_FR_2018_TWO_WWTP_DUAL_OUTLET_SOURCE",json.dumps({
        k:v for k,v in r.items() if k not in ("source_rows","source_receipts_SHA256")
    },sort_keys=True))


if __name__=="__main__":
    main()
