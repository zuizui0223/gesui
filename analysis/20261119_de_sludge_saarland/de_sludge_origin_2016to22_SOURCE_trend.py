#!/usr/bin/env python3
"""German 2016–22 sewage-sludge source-routing time series, NOT field exposure.

Originals: Destatis GENESIS 32214-0001 national annually; Saarland 2023
Statistisches Jahrbuch table 13 Umwelt no.8, printed page 143.
All quantities are t dry matter, attributed by WWTP origin, not land receiving site.
Only descriptive origin-flow comparison; no N/P, river, or biotic outcomes.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

NATIONAL_SOURCE="https://genesis.destatis.de/datenbank/online/statistic/32214/table/32214-0001"
SAARLAND_SOURCE="https://www.statistischebibliothek.de/mir/servlets/MCRFileNodeServlet/SLHeft_derivate_00012170/SL%20Jahrbuch%202023.pdf"
HEADER=("year","DE_total_dry_t","DE_agricultural_dry_t","DE_landscape_dry_t",
        "DE_other_material_dry_t","DE_incineration_mono_dry_t",
        "DE_incineration_co_dry_t","DE_incineration_unknown_dry_t",
        "DE_other_disposal_dry_t","SL_total_dry_t","SL_agricultural_dry_t",
        "SL_landscape_dry_t","SL_thermal_dry_t","SL_other_dry_t")
YEARS=list(range(2016,2023))

def origin_routing(rows):
    if len(rows)!=len(YEARS):raise ValueError("Source must contain seven published years")
    out=[]
    for i,row in enumerate(rows):
        if tuple(row)!=HEADER:
            raise ValueError("Unexpected schema: cannot silently reinterpret source year/types")
        try:
            year=int(row["year"])
            v={k:int(row[k]) for k in HEADER if k!="year"}
        except (ValueError,KeyError,TypeError) as exc:
            raise ValueError("NA/negative/noninteger source units prohibited") from exc
        if year!=YEARS[i] or any(n<0 for n in v.values()):
            raise ValueError("Invalid time order or source mass")
        de_mass=(v["DE_agricultural_dry_t"]+v["DE_landscape_dry_t"]+
                 v["DE_other_material_dry_t"]+v["DE_incineration_mono_dry_t"]+
                 v["DE_incineration_co_dry_t"]+v["DE_incineration_unknown_dry_t"]+
                 v["DE_other_disposal_dry_t"])
        sl_mass=(v["SL_agricultural_dry_t"]+v["SL_landscape_dry_t"]+
                 v["SL_thermal_dry_t"]+v["SL_other_dry_t"])
        if de_mass!=v["DE_total_dry_t"] or sl_mass!=v["SL_total_dry_t"]:
            raise ValueError("Published total does not close by disposal routes")
        if v["SL_agricultural_dry_t"]>v["DE_agricultural_dry_t"] or v["SL_total_dry_t"]>v["DE_total_dry_t"]:
            raise ValueError("State values exceed their national subset")
        rest_total=v["DE_total_dry_t"]-v["SL_total_dry_t"]
        rest_agric=v["DE_agricultural_dry_t"]-v["SL_agricultural_dry_t"]
        if rest_total<=0 or rest_agric<0:
            raise ValueError("National less Saarland must remain positive")
        de_thermal=(v["DE_incineration_mono_dry_t"]+
                    v["DE_incineration_co_dry_t"]+
                    v["DE_incineration_unknown_dry_t"])
        out.append({
            "year":year,
            "DE_total_dry_t":v["DE_total_dry_t"],
            "DE_agricultural_dry_t":v["DE_agricultural_dry_t"],
            "DE_agriculture_share":v["DE_agricultural_dry_t"]/v["DE_total_dry_t"],
            "DE_thermal_dry_t":de_thermal,
            "DE_thermal_share":de_thermal/v["DE_total_dry_t"],
            "SL_total_dry_t":v["SL_total_dry_t"],
            "SL_agricultural_dry_t":v["SL_agricultural_dry_t"],
            "SL_agriculture_share":v["SL_agricultural_dry_t"]/v["SL_total_dry_t"],
            "RestDE_total_dry_t":rest_total,
            "RestDE_agricultural_dry_t":rest_agric,
            "RestDE_agriculture_share":rest_agric/rest_total
        })
    by_year={x["year"]:x for x in out}
    first,last=by_year[2017],by_year[2022]
    change=lambda key: (last[key]/first[key]-1)
    diff_of_share_changes=(last["SL_agriculture_share"]-first["SL_agriculture_share"])-(
                           last["RestDE_agriculture_share"]-first["RestDE_agriculture_share"])
    return {
       "schema":"gesui_DE_sludge_2016to22_stated_origin_routes_source_only_v1",
       "primary_source_National_GENESIS_table":NATIONAL_SOURCE,
       "primary_source_Saarland_2023_yearbook_page_143":SAARLAND_SOURCE,
       "unit":"tonnes_dry_matter_per_year",
       "annual_rows":out,
       "2017to2022":{
          "DE_agri_percent_change":round(100*change("DE_agricultural_dry_t"),3),
          "SL_agri_percent_change":round(100*change("SL_agricultural_dry_t"),3),
          "RestDE_agri_percent_change":round(100*change("RestDE_agricultural_dry_t"),3),
          "SL_minus_RestDE_change_in_agri_share_percentage_points_DESCRIPTIVE_ONLY":round(100*diff_of_share_changes,3),
          "DE_agri_share_change_percentage_points":round(100*(last["DE_agriculture_share"]-first["DE_agriculture_share"]),3),
          "DE_thermal_share_change_percentage_points":round(100*(last["DE_thermal_share"]-first["DE_thermal_share"]),3),
       },
       "national_includes_Saarland_avoid_pseudoreplication":True,
       "non_overlapping_comparison_uses_Saarland_vs_rest_DE":True,
       "data_release_dates_and_revisions_differ":True,
       "regional_allocation_by_WWTP_origin_not_actually_receiving_farmland":True,
       "routing_share_is_not_plant_or_butterfly_exposure":True,
       "sludge_dry_mass_is_not_N_or_P_mass":True,
       "thermal_share_not_evidence_of_more_aqueous_N_P_discharge":True,
       "calendar_years_not_independent_experiments":True,
       "no_hypothesis_test_or_p_values":True,
       "paired_water_and_land_2018_N_P_export_mass_complete":False,
       "source_GO_for_biological_model":False,
       "no_EQR_NUTRI_HNS_used":True,
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--input",required=True,type=Path)
    a.add_argument("--out",required=True,type=Path)
    args=a.parse_args()
    raw=args.input.read_bytes()
    with args.input.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=HEADER:
            raise ValueError("Source field schema differs")
        rows=list(reader)
    result=origin_routing(rows)
    result["original_source_transcription_SHA256"]=hashlib.sha256(raw).hexdigest()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GESUI_DE_SLUDGE_ORIGIN_ONLY_SOURCE_TREND",json.dumps(result["2017to2022"],sort_keys=True))

if __name__=="__main__": main()
