#!/usr/bin/env python3
"""One French station in 2018: 12 treated-water A4 months vs realized sludge event months.

Measured outlet water volume is NOT a measured TN/PT export, river total water,
or independent monthly causal replicates. Month-specific sludge N/P and wet
amounts are unknown. Bypass telemetry A5 unreliable in original report.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

COLUMNS=("station_SANDRE","year","month","reported_treated_A4_volume_m3",
         "actual_sludge_application_log_events_in_month",
         "actual_sludge_wet_mass_by_month_t",
         "monthly_outlet_TN_mass_measured_kg","monthly_outlet_TP_mass_measured_kg",
         "A5_bypass_volume_reliable")
A4=(36168,28509,32800,26858,27330,24918,24108,25408,25399,26716,29400,35607)
EVENTS=(0,0,5,0,0,0,0,0,3,0,0,0)


def audit(rows):
    if len(rows)!=12:raise ValueError("Require all 12 months of 2018")
    monthly=[]
    for i,row in enumerate(rows):
        if tuple(row)!=COLUMNS or row["station_SANDRE"]!="0517347V004" or row["year"]!="2018":
            raise ValueError("Invalid station/year/field contract")
        if row["month"]!=f"{i+1:02d}":
            raise ValueError("Month order/coverage changed")
        if int(row["reported_treated_A4_volume_m3"])!=A4[i]:
            raise ValueError("Original annual A4 table monthly volume changed")
        if int(row["actual_sludge_application_log_events_in_month"])!=EVENTS[i]:
            raise ValueError("Original agronomy register month-event count changed")
        if any(row[field] for field in ("actual_sludge_wet_mass_by_month_t",
                          "monthly_outlet_TN_mass_measured_kg",
                          "monthly_outlet_TP_mass_measured_kg")):
            raise ValueError("Unmeasured monthly mass N/P must stay missing")
        if row["A5_bypass_volume_reliable"]!="false":
            raise ValueError("Failed bypass meter cannot become confirmed zero")
        monthly.append({"month":f"{i+1:02d}","A4_treated_water_volume_m3":A4[i],
                        "actual_sludge_log_events":EVENTS[i],
                        "actual_sludge_wet_mass_by_month_t":None,
                        "aqueous_export_TN_kg":None,"aqueous_export_TP_kg":None})
    annual=sum(A4)
    if annual!=343221 or sum(EVENTS)!=8:
        raise ValueError("Published annual A4 total or 2018 ledger failed closure")
    pulse_months=[r for r in monthly if r["actual_sludge_log_events"]>0]
    if [r["month"] for r in pulse_months]!=["03","09"]:
        raise ValueError("Sludge pulses shifted")
    pulse_A4_m3=sum(x["A4_treated_water_volume_m3"] for x in pulse_months)
    return {
        "schema":"gesui_FR_SaintJean2018_one_STEU_water_month_land_pulse_SOURCE_only_v1",
        "STEU":"0517347V004","year":2018,
        "annual_A4_treated_effluent_m3_reported":annual,
        "annual_A4_is_not_total_receiving_river_wastewater_when_A5_unreliable":True,
        "positive_A4_treated_outlet_flow_months":sum(x>0 for x in A4),
        "months_with_realized_sludge_application_logs":[x["month"] for x in pulse_months],
        "realized_agricultural_sludge_ledger_events":sum(EVENTS),
        "water_A4_volume_m3_in_months_with_sludge_events":pulse_A4_m3,
        "fraction_annual_A4_treated_water_in_sludge_event_months_DESCRIPTIVE_ONLY":round(pulse_A4_m3/annual,6),
        "monthly_actual_sludge_N_P_mass_unidentified":True,
        "monthly_aqueous_NGL_PT_export_mass_unidentified":True,
        "A5_bypass_meter_malfunction_preserved":True,
        "NGL_PT_2018_influent_annual_load_not_effluent":True,
        "independent_station_years":1,
        "12_months_are_not_independent_stations":True,
        "no_EQR_NUTRI_HNS_values_read":True,
        "monthly_causal_coupled_nutrient_transfer_gate":"STOP",
        "monthly_source_opportunity_timing":"GO",
        "months":monthly
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    with args.input.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=COLUMNS:
            raise ValueError("Unexpected monthly source schema, including private columns")
        rows=list(reader)
    r=audit(rows)
    r["input_SHA256"]=hashlib.sha256(args.input.read_bytes()).hexdigest()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(r,sort_keys=True,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GESUI_SJ_2018_WATER_OUTFLOW_AND_REAL_SLUDGE_MONTH",json.dumps({
        k:v for k,v in r.items() if k not in ("months","input_SHA256")
    },sort_keys=True))


if __name__=="__main__":
    main()
