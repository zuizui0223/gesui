#!/usr/bin/env python3
"""Official EA new Water Quality Explorer: bounded 2017–18 treated-STW TP source.

Source-only census of actual 0348 Total Phosphorus as P (mg/L) sample
observations, limited to native Water Company FINAL/TREATED EFFLUENT (SA)
near Bath within 50 km. This is not a nationwide count, nor an outfall
network join nor proof of 2018 N/P annual kilograms or causality.

The 2025-12 EA archive migration documents POST /data/observation with
API-Version:1, limit and skip. Refuse to extrapolate a truncated page search.
Do not publish site ID, street, easting/northing, or raw concentration rows;
only dataset-level aggregate source coverage and source SHA receipts.
"""
from __future__ import annotations
import argparse,collections,csv,datetime,hashlib,io,json,math,time
import urllib.error,urllib.parse,urllib.request
from pathlib import Path

HOST="https://environment.data.gov.uk/water-quality"
PATH="/data/observation"
ANALYTE_TOTAL_P_AS_P="0348"
BOUND={"latitude":"51.383","longitude":"-2.363","radius":"50"}
YEARS=(2017,2018)
LIMIT=250
MAX_PAGES_PER_YEAR=12
MAX_BYTES_PER_PAGE=210000
UA="gesui-EA-WQE-native-STW-final-effluent-2017-18-original-census/20261120"
HEADERS={"User-Agent":UA,"API-Version":"1","Accept":"text/csv",
         "CSV-Header":"present","Content-Type":"application/json"}
SA="SEWAGE DISCHARGES - FINAL/TREATED EFFLUENT - WATER COMPANY"
SB="SEWAGE DISCHARGES - STW STORM OVERFLOW/STORM TANK - WATER COMPANY"
UA_TYPE="SEWAGE DISCHARGES - FINAL/TREATED EFFLUENT - NOT WATER COMPANY"
REQUIRED={
  "samplingPoint.notation","samplingPoint.samplingPointType",
  "sampleMaterialType","samplingPurpose","phenomenonTime",
  "determinand.notation","determinand.prefLabel","result","unit"
}

def fetch(year,skip):
    params={
        "dateFrom":str(year)+"-01-01","dateTo":str(year)+"-12-31",
        "determinand":ANALYTE_TOTAL_P_AS_P,"limit":LIMIT,"skip":skip,
        **BOUND}
    url=HOST+PATH+"?"+urllib.parse.urlencode(params)
    request=urllib.request.Request(url,method="POST",data=b"",headers=HEADERS)
    receipt={"year":year,"skip":skip,"status":None,"content_type":None,
             "nbytes":None,"sha256":None,"error":None}
    try:
        with urllib.request.urlopen(request,timeout=50) as r:
            receipt["status"]=int(r.status)
            receipt["content_type"]=r.headers.get("Content-Type","")
            raw=r.read(MAX_BYTES_PER_PAGE+1)
        receipt["nbytes"]=len(raw)
        if len(raw)>MAX_BYTES_PER_PAGE:
            raise RuntimeError("Original WQE response exceeded bounded source page")
        receipt["sha256"]=hashlib.sha256(raw).hexdigest()
        if receipt["status"]==204:return [],receipt
        if receipt["status"]!=200:
            raise RuntimeError("Unexpected official HTTP source response status")
        text=raw.decode("utf-8-sig","replace")
        parser=csv.DictReader(io.StringIO(text))
        head=parser.fieldnames or []
        if not REQUIRED.issubset(set(head)):
            raise RuntimeError("Source schema lacks required native WQE observation keys")
        rows=list(parser)
        if len(rows)>LIMIT:
            raise RuntimeError("Official page exceeded fixed 250 records")
        receipt["n_source_rows"]=len(rows)
        return rows,receipt
    except urllib.error.HTTPError as err:
        receipt["status"]=err.code
        receipt["error"]="HTTPError "+str(err)[:220]
        return None,receipt
    except Exception as err:
        receipt["error"]=type(err).__name__+": "+str(err)[:320]
        return None,receipt

def parse(year):
    rows=[]
    receipts=[]
    truncated=False
    complete=False
    for i in range(MAX_PAGES_PER_YEAR):
        page,source=fetch(year,i*LIMIT)
        receipts.append(source)
        if page is None:break
        if len(page)<LIMIT:
            rows.extend(page);complete=True;break
        rows.extend(page)
        if i==MAX_PAGES_PER_YEAR-1:truncated=True
    alltypes=collections.Counter()
    sampletypes=collections.Counter()
    purposes=collections.Counter()
    units=collections.Counter()
    months=collections.Counter()
    n_p_numeric=0
    n_censored=0
    n_invalid_date=0
    n_mislabelled=0
    treated_codes=set()
    treated_sample_dates=set()
    n_real_SA=0
    n_real_SA_UWWTD=0
    n_real_SA_mgL=0
    n_non_water_company_final=0
    n_stormoverflow=0
    n_river_lake=0
    sa_site_counts=collections.Counter()
    for r in rows:
        recorded_date=r["phenomenonTime"].strip()
        if not recorded_date.startswith(str(year)):
            n_invalid_date+=1
        name=r["samplingPoint.samplingPointType"].strip()
        material=r["sampleMaterialType"].strip()
        purpose=r["samplingPurpose"].strip()
        unit=r["unit"].strip()
        alltypes[name]+=1
        sampletypes[material]+=1
        purposes[purpose]+=1
        units[unit]+=1
        if r["determinand.notation"].strip()!=ANALYTE_TOTAL_P_AS_P:
            n_mislabelled+=1
            continue
        val=r["result"].strip()
        try:
            value=float(val)
            if math.isfinite(value) and value>=0:n_p_numeric+=1
            else:n_censored+=1
        except (ValueError,TypeError):n_censored+=1
        if name==SA:
            treated_codes.add(r["samplingPoint.notation"].strip())
            sa_site_counts[r["samplingPoint.notation"].strip()]+=1
            treated_sample_dates.add((r["samplingPoint.notation"].strip(),recorded_date[:10]))
            n_real_SA+=1
            if "UWWTD" in purpose.upper():n_real_SA_UWWTD+=1
            if unit.upper()=="MILLIGRAM PER LITRE":n_real_SA_mgL+=1
            if len(recorded_date)>=7:months[recorded_date[:7]]+=1
        elif name==UA_TYPE:
            n_non_water_company_final+=1
        elif name==SB:
            n_stormoverflow+=1
        else:
            n_river_lake+=1
    return {
      "year":year,
      "requested_radius":"50 km from Bath SW England (test geography only)",
      "source_pagination_complete":bool(complete),
      "source_pagination_limit_hit":truncated,
      "source_any_HTTP_error":any(x["error"] is not None for x in receipts),
      "source_pages_requested":len(receipts),
      "source_page_receipts":receipts,
      "source_observation_rows_seen":len(rows),
      "original_determinand_mismatched_rows":n_mislabelled,
      "source_year_mismatch_rows":n_invalid_date,
      "observed_source_sample_point_type_counts":dict(alltypes),
      "observed_source_material_type_counts":dict(sampletypes),
      "observed_source_sampling_purpose_counts":dict(purposes),
      "observed_source_unit_counts":dict(units),
      "source_nonnegative_numeric_results_in_full_page":n_p_numeric,
      "source_censored_or_invalid_results":n_censored,
      "source_valid_SA_water_company_treated_effluent_rows":n_real_SA,
      "source_valid_SA_unique_native_sampling_points":len(treated_codes),
      "source_valid_SA_native_sampling_point_date_pairs":len(treated_sample_dates),
      "source_valid_SA_records_UWWTD_monitoring_purpose":n_real_SA_UWWTD,
      "source_valid_SA_records_mg_per_L":n_real_SA_mgL,
      "source_other_non_WaterCompany_UA_rows":n_non_water_company_final,
      "source_storm_overflow_SB_rows":n_stormoverflow,
      "source_other_river_lake_groundwater_rows":n_river_lake,
      "source_SA_month_count":dict(sorted(months.items())),
      "source_SA_site_coding_SHA256":hashlib.sha256(("\n".join(sorted(treated_codes))).encode()).hexdigest(),
      "_SA_station_set":treated_codes,
      "_SA_station_counts":sa_site_counts,
    }

def main(a):
    yearly=[parse(yr) for yr in YEARS]
    common=set.intersection(*(z["_SA_station_set"] for z in yearly))
    data={
      "schema":"GESUI_20261120_EA_WQE_2017_18_STW_TP_SA_SAMPLE_SOURCE_CENSUS_V1",
      "status":"SOURCE_CENSUS_BOUNDED_COMPLETE",
      "source_only_no_biological_outcome_or_flow_computed":True,
      "official_WQE_original_determinand_code":ANALYTE_TOTAL_P_AS_P,
      "official_WQE_observation_units":"mg/L only when actual original row labelled MILLIGRAM PER LITRE",
      "original_years":list(YEARS),
      "original_search_geography":BOUND,
      "not_a_nationwide_source_population":True,
      "not_original_WWTP_discharge_mass_kg_per_year":True,
      "not_downstream_connection_to_ChemPop":True,
      "original_STW_final_effluent_type_required":SA,
      "exclude_water_company_storm_overflow_type":SB,
      "exclude_non_water_company_treated_type":UA_TYPE,
      "yearly":[],
      "cross_year":{
         "n_SA_stations_seen_in_both_years_in_same_test_geography":len(common),
         "source_SA_shared_stations_SHA256":hashlib.sha256(("\n".join(sorted(common))).encode()).hexdigest(),
      }
    }
    for v in yearly:
        out={k:value for k,value in v.items() if not k.startswith("_")}
        data["yearly"].append(out)
        if not (out["source_pagination_complete"] and
                out["original_determinand_mismatched_rows"]==0 and
                out["source_year_mismatch_rows"]==0):
            data["status"]="SOURCE_CENSUS_INCOMPLETE_OR_SOURCE_GATE_STOP"
    data["source_gates"]={
        "has_verified_WaterCompany_SA_final_TP_2017_18_record":any(y["source_valid_SA_records_mg_per_L"]>0 for y in data["yearly"]),
        "has_repeated_native_SA_station_two_years":len(common)>0,
        "at_least_30_official_STW_SA_points_reliably_across_two_years":len(set.union(*(x["_SA_station_set"] for x in yearly)))>=30,
        "original_SA_sites_to_ChemPop_river_sites_route_verified":False,
        "original_measured_WWTP_discharge_total_P_or_N_kg_yr_verified":False,
        "geographic_pilot_only_not_source_model_authorized":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    summary={
        "status":data["status"],"source_gates":data["source_gates"],
        "cross_year":data["cross_year"],
        "yearly":[{
           k:z[k] for k in (
             "year","source_pagination_complete","source_pagination_limit_hit",
             "source_pages_requested","source_observation_rows_seen",
             "source_valid_SA_water_company_treated_effluent_rows",
             "source_valid_SA_unique_native_sampling_points",
             "source_valid_SA_native_sampling_point_date_pairs",
             "source_valid_SA_records_UWWTD_monitoring_purpose",
             "source_valid_SA_records_mg_per_L",
             "source_other_non_WaterCompany_UA_rows","source_storm_overflow_SB_rows",
             "source_other_river_lake_groundwater_rows",
             "source_sample_point_type_error_rows" if "source_sample_point_type_error_rows" in z else "source_censored_or_invalid_results",
           )
        } for z in data["yearly"]]}
    print("GESUI_EA_WQE_2017_18_SA_TP_CENSUS="+json.dumps(summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
