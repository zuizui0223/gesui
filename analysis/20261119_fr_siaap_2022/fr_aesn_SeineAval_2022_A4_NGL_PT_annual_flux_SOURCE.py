#!/usr/bin/env python3
"""2022 AESN Seine Aval same-year monthly A4 effluent NGL/PT source flux.

Source discovery and original SHA were independently fixed BEFORE numeric
values were opened in Github Actions #44536921058 (file 55,938,440 bytes,
SHA256 5fbed...). Read only exact Seine Aval SANDRE & original source 2022,
A4 outlet and A3 inlet/A5 bypass comparison, codes 1551=NGL 1350=PT;
monthly indicator units exactly kg, Etat_calcul exactly VALIDE. A4 volume
1552 is m3; no mixing input and output or A4 and bypass. Substantive results
are official *validated monthly calculated autosurveillance indicators*,
not individual gravimetric direct assays or proven full-river load.
"""
import argparse,csv,hashlib,json,sys
from collections import defaultdict
from decimal import Decimal,InvalidOperation
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fr_aesn_2022_performance_exact_sandre_source_coverage import (
    fetch_catalog,original_resource,download_ranges,EXPECTED_COLS
)

PIN_SHA="5fbed27273b0c2c1d4590e9286fa52e724578b5cb19e68369fd1fb8e77347af7"
PIN_BYTES=55938440
STEU="037800501000"
POINTS=("A3","A4","A5")
PARAMS=("1551","1350")
VOLUME="1552"
YEAR="2022"
ORIGINAL_GHA_DISCOVERY=44536921058

def decimal_value(s):
    raw=str(s or "").strip().replace(" ","").replace("\u00a0","")
    if not raw:raise ValueError("Validated monthly indicator has no numeric result")
    try:z=Decimal(raw.replace(",","."))
    except InvalidOperation as ex:raise ValueError("Invalid original monthly mass") from ex
    if not z.is_finite() or z<0 or z>Decimal("1e12"):
        raise ValueError("Source mass negative, nonfinite or out-of-range")
    return z

def scan(path):
    records={}
    source_station_rows=0
    with Path(path).open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f,delimiter=";")
        if tuple(reader.fieldnames or ())!=EXPECTED_COLS:
            raise ValueError("Official source CSV columns changed since source-only audit")
        for row in reader:
            code=str(row.get("Station_Code SANDRE") or "").strip()
            if code not in (STEU,"FR"+STEU):
                continue
            if str(row.get("Station_Nom") or "").strip().upper()!="SEINE AVAL":
                raise ValueError("Original SANDRE ID and station name conflict")
            if str(row.get("Annee") or "").strip()!=YEAR:
                raise ValueError("Wrong original source year for same-source exposure")
            source_station_rows+=1
            point=str(row.get("Point_Prelevement_Code") or "").strip().upper()
            param=str(row.get("Parametre_Code SANDRE") or "").strip()
            if (point not in POINTS or param not in PARAMS) and not (point=="A4" and param==VOLUME):
                continue
            month=int(str(row.get("Mois") or "").strip())
            if not 1<=month<=12:raise ValueError("Noncalendar reporting month")
            key=(point,param,month)
            if key in records:raise ValueError("Duplicate ST EU outlet/parameter/month")
            state=str(row.get("Etat_calcul") or "").strip().upper()
            unit=str(row.get("Resultat_Unite_Symbole") or "").strip().lower().replace("³","3")
            unit_code=str(row.get("Resultat_Unite_Code") or "").strip()
            if state!="VALIDE":raise ValueError("Unvalidated or estimated result status")
            if param==VOLUME:
                if (unit_code,unit)!=("115","m3"):raise ValueError("A4 water volume unit changed")
            else:
                if (unit_code,unit)!=("67","kg"):raise ValueError("N/P source mass not in kg")
            records[key]=decimal_value(row.get("Resultat_Valeur"))
    if source_station_rows!=528:raise ValueError("2022 station source count changed from frozen 528")
    required={(point,param,month) for point in POINTS for param in PARAMS for month in range(1,13)}
    required|={("A4",VOLUME,month) for month in range(1,13)}
    missing=required-set(records)
    if missing:
        raise ValueError("Required 2022 monthly data incomplete: "+str(sorted(missing)[:10]))
    if set(records)!=required:
        raise ValueError("Unexpected indicator values not in source predeclared cohort")
    monthly=[]
    for month in range(1,13):
        row={"month":f"{month:02d}"}
        for point in POINTS:
            row[f"{point}_NGL_kg"]=float(records[(point,"1551",month)])
            row[f"{point}_PT_kg"]=float(records[(point,"1350",month)])
        row["A4_treated_water_m3"]=float(records[("A4",VOLUME,month)])
        monthly.append(row)
    annual={}
    for point in POINTS:
        for code,label in (("1551","NGL"),("1350","PT")):
            total=sum((records[(point,code,m)] for m in range(1,13)),Decimal(0))
            annual[f"{point}_{label}_kg_2022"]=float(total)
            annual[f"{point}_{label}_t_2022"]=float(total/Decimal(1000))
    volume=sum((records[("A4",VOLUME,m)] for m in range(1,13)),Decimal(0))
    annual["A4_treated_water_m3_2022"]=float(volume)
    # Useful comparison but not all flow pathways. Only actual treated A4.
    return {
       "schema":"gesui_AESN2022_SeineAval_monthly_validated_calculated_A4_NGL_PT_source_mass_v1",
       "official_station_SANDRE":STEU,
       "original_reporting_year":2022,
       "original_source_frozen_ahead_of_values_sha256":PIN_SHA,
       "original_source_discovery_run":ORIGINAL_GHA_DISCOVERY,
       "original_STEU_rows":source_station_rows,
       "monthly_2022_source_indicator_count":len(records),
       "monthly":monthly,
       "annual_source_flux_indicators":annual,
       "indicator_aggregation":"Sum of 12 monthly VALIDE results, original SANDRE code 1551/1350 unit kg; A4 treated-outlet only",
       "A3_influent_not_effluent":True,
       "A5_bypass_separate_not_summed_into_A4":True,
       "A4_treated_effluent_is_not_total_plant_or_river_flux":True,
       "source_monthly_indicators_not_direct_daily_measurements":True,
       "not_a_closed_WWTP_N_P_mass_balance":True,
       "sludge_field_specific_N_P_mass_and_plant_uptake_not_inferred":True,
       "no_EQR_NUTRI_HNS_used":True,
       "biological_effect_causal_path_GATE":"STOP",
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    a.work.mkdir(parents=True,exist_ok=True)
    dest=a.work/"stepcol-indicateurs-mensuels-2022.csv"
    result={"schema":"gesui_AESN2022_SeineAval_A4_NGL_PT_sameyear_source_only_execution_v1",
            "status":"SOURCE_ACCESS_STOP","expected_SHA256":PIN_SHA,
            "STEU":STEU,"year":2022,"biological_outcomes_read":False}
    try:
        cat,catsha=fetch_catalog()
        url,size=original_resource(cat)
        if size!=PIN_BYTES:raise ValueError("Source file size changed since original frozen discovery")
        receipt=download_ranges(url,size,dest)
        if receipt["original_SHA256"]!=PIN_SHA:
            raise ValueError("Source original SHA differs from pre-numeric-value freeze")
        result["original_source_catalog_SHA256"]=catsha
        result["original_source_receipt"]={k:v for k,v in receipt.items() if k!="source_range_receipts"}
        result["actual_source"]=scan(dest)
        result["status"]="OFFICIAL_2022_VALIDATED_MONTHLY_A4_NGL_PT_MASS_RECONSTRUCTED"
    except Exception as e:
        result["source_error"]=type(e).__name__+": "+str(e)[:450]
    finally:
        dest.unlink(missing_ok=True)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    s=result.get("actual_source",{}).get("annual_source_flux_indicators",{})
    print("GESUI_SEINEAVAL_2022_A4_NGL_PT_ANNUAL_OFFICIAL_INDICATORS",json.dumps({
        "status":result["status"],"error":result.get("source_error"),
        "annual_official_indicator_values":s,"biology_used":False},sort_keys=True))
if __name__=="__main__":main()
