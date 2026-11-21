#!/usr/bin/env python3
"""Authentic EEA 2017-18 UWWTD v07 r01 UK discharge point source-only audit.

A zero-UK count in the STRICT WWTP both-P-and-N-measured-positive 2018
subset does NOT logically establish that England's original historic
DischargePoints.csv had zero official UK outlet geometries. Verify distinct
2018 ReportPeriod reporter keys, UK/GB plant rows, native outfall coordinates,
source-type FW and observed plant states by EXACT official original member
SHA256 before any release-to-river joining. Source metadata only, no row
geography, native facility IDs or WWTP source discharge masses rehosted.
Original annual P/N missing ≠ zero emissions. Physical 2018 original
reported UWWTD discharge locations may differ from EA WQE SA sample chamber
or modern permit outlets; join evidence remains necessary.

Historical source: archived EEA public Waterbase UWWTD reporting 2017-18,
which was previously confirmed for the gesui 2018 EQR source-only work.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,tempfile,urllib.request,zipfile
from pathlib import Path

URL="https://sdi.eea.europa.eu/datashare/s/mGmXg3zcEqKQLTG/download"
BASE="eea_t_uwwtd-waterbase_p_2017-2018_v07_r01/Waterbase_UWWTD_v8_csv/"
EXPECT={
 "DischargePoints.csv":"fe96e70df6b169222b4aaab78201f687ecfadfe8bec2bfe8554a93784ee112da",
 "ReportPeriod.csv":"30bed8a962218c1a2c772ba1bd151fc16d058874f77a4aef263ba2d83dfeb9fd",
 "UWWTPS.csv":"bfd60e8ecc1539d4e4c31ea1182eafef2653e8ac46680e09f268f1a889d7fa48"
}
MAX_BYTES=100*1024*1024
UA="gesui-official-EEA-2018-UK-physical-DischargePoints-source-audit/20261121"
def finite(x):
    try:
        f=float(str(x).strip())
        return f if abs(f)<10000 and f==f else None
    except (ValueError,TypeError):return None

def main(a):
    result={"schema":"GESUI_20261121_UWWTD2018_UK_DISCHARGEPOINTS_SOURCE_ONLY_V1",
        "source_only":True,"original_official_url":URL,
        "historical_2018_source_plant_data_not_2022":True,
        "strict_measured_both_PN_zero_UK_does_not_imply_no_UK_outlets":True,
        "not_any_new_BMWP_or_plant_ecological_response":True,
        "not_proof_EA_SA_monitoring_point_is_licensed_UWWTD_outfall":True,
        "status":"SOURCE_STOP"}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_EEA_2018_UK_UWWTD_") as tmp:
            f=Path(tmp)/"EEA_original.zip"
            req=urllib.request.Request(URL,headers={"User-Agent":UA})
            h=hashlib.sha256();n=0
            with urllib.request.urlopen(req,timeout=100) as r,f.open("wb") as dest:
                if r.status!=200:raise RuntimeError("Official EEA archive HTTP!="+"200")
                while True:
                    chunk=r.read(131072)
                    if not chunk:break
                    n+=len(chunk)
                    if n>MAX_BYTES:raise RuntimeError("Official EEA original archive unexpectedly large")
                    h.update(chunk);dest.write(chunk)
            result["official_zip_bytes"]=n
            result["official_archive_SHA256"]=h.hexdigest()
            with zipfile.ZipFile(f) as z:
                filenames=z.namelist()
                manifest={}
                tables={}
                for name,digest in EXPECT.items():
                    full=BASE+name
                    data=z.read(full)
                    sha=hashlib.sha256(data).hexdigest()
                    if sha!=digest:raise RuntimeError("Original EEA 2018 source member SHA changed: "+name)
                    manifest[name]={"bytes":len(data),"SHA256":sha,
                                    "columns":None}
                    with io.TextIOWrapper(io.BytesIO(data),encoding="utf-8-sig",errors="replace",newline="") as wrapper:
                        reader=csv.DictReader(wrapper)
                        manifest[name]["columns"]=reader.fieldnames
                        tables[name]=list(reader)
            if len(filenames)!=23:raise RuntimeError("Historical original EEA 23 files source membership changed")
        periods=tables["ReportPeriod.csv"]
        plants=tables["UWWTPS.csv"]
        outfalls=tables["DischargePoints.csv"]
        counts=collections.Counter()
        reporterperiod=collections.defaultdict(list)
        code_to_country={}
        for row in periods:
            code=str(row.get("repCode","")).strip()
            country=str(row.get("rptMStateKey","")).strip()
            year=str(row.get("repReportedPeriod","")).strip()
            situation=str(row.get("repSituationAt","")).strip()
            if code:
                reporterperiod[code].append((country,year,situation))
                code_to_country[code]=country
                counts["period_by_"+country+"_"+year]+=1
        uk_reps=sorted([code for code,rows in reporterperiod.items()
                    if any(c in {"UK","GB","GBR","United Kingdom"} and y=="2018" and "2018" in sit
                           for c,y,sit in rows)])
        result["source_report_period_summary"]={
           "original_report_period_n":len(periods),
           "report_period_country_and_year_counts":dict(sorted(counts.items())),
           "all_reporter_code_count":len(reporterperiod),
           "UK_or_GB_repcodes_with_explicit_2018_report_year":uk_reps,
           "original_2018_ref_reporters":sorted(
               c for c,rows in reporterperiod.items() if any(y=="2018" and "2018" in si for _,y,si in rows))
        }
        reported_all_uk=[(k,r) for k in reporterperiod for r in reporterperiod[k]
                         if r[0] in {"UK","GB","GBR","United Kingdom"}]
        result["original_all_period_UK_or_GB_period_rows"]=len(reported_all_uk)
        subset={}
        for label,rows in (("plants",plants),("outfalls",outfalls)):
            count=collections.Counter(str(r.get("repCode","")).strip() for r in rows)
            matching=[r for r in rows if str(r.get("repCode","")).strip() in uk_reps]
            result["source_"+label]={
               "original_all_reporters_n":len(rows),
               "original_unique_repCode_count":len(count),
               "original_num_rows_with_UK_GB_2018_reportCode":len(matching),
               "source_2018_UK_candidate_repCode_counts":{c:count[c] for c in uk_reps},
            }
            subset[label]=matching
        plant_rows=subset["plants"]
        out_rows=subset["outfalls"]
        plantkeys={(r.get("repCode",""),r.get("uwwCode","")) for r in plant_rows}
        count_keys=collections.Counter((r.get("repCode",""),r.get("uwwCode","")) for r in out_rows)
        coordstatus=collections.Counter()
        wb=collections.Counter()
        surface=collections.Counter()
        matched_keys=set()
        unique_reported_outlet_n=0
        unique_coord_outlet_n=0
        measured_n=0;measured_p=0;measured_both=0
        for r in plant_rows:
            try:
                n_val=float(r.get("uwwNDischargeMeasured",""))
                n_ok=n_val>0
            except (ValueError,TypeError):n_ok=False
            try:
                p_val=float(r.get("uwwPDischargeMeasured",""))
                p_ok=p_val>0
            except (ValueError,TypeError):p_ok=False
            measured_n+=int(n_ok);measured_p+=int(p_ok);measured_both+=int(n_ok and p_ok)
        for row in out_rows:
            key=(row.get("repCode",""),row.get("uwwCode",""))
            if key in plantkeys:matched_keys.add(key)
            lat,lon=finite(row.get("dcpLatitude")),finite(row.get("dcpLongitude"))
            if lat is not None and lon is not None and 49<=lat<=61 and -12<=lon<=4:
                unique_coord_outlet_n+=1
                coordstatus["original_uk_geo_plausible_outfall"]+=1
            else:coordstatus["original_missing_or_outside_UK_bbox"]+=1
            wb[str(row.get("dcpWaterBodyType","")).strip()]+=1
            surface[str(row.get("dcpSurfaceWaters","")).strip()]+=1
            unique_reported_outlet_n+=1
        result.update({
           "status":"OFFICIAL_EEA_UWWTD_2018_UK_OUTLET_SOURCE_CENSUS_COMPLETED",
           "official_member_SHA256":manifest,
           "original_UK_or_GB_2018_outlet_QA":{
               "UK_2018_reporters":len(uk_reps),
               "original_2018_UK_plant_rows":len(plant_rows),
               "original_2018_UK_outlet_rows":len(out_rows),
               "outlet_rows_with_valid_plausible_UK_coordinate":unique_coord_outlet_n,
               "outlet_geo_status":dict(coordstatus),
               "outlet_waterbody_type_codes":dict(wb),
               "outlet_surfacewater_codes":dict(surface),
               "original_distinct_UK_plant_keys_with_any_matched_outlet":len(matched_keys),
               "original_2018_UK_plant_keys_with_exactly_one_outlet":
                   sum(count_keys[k]==1 for k in plantkeys),
               "original_2018_UK_plant_keys_with_2plus_outlets":
                   sum(count_keys[k]>=2 for k in plantkeys),
               "original_2018_UK_plant_keys_with_no_outlet":
                   sum(count_keys[k]==0 for k in plantkeys),
               "original_2018_UK_plant_rows_with_measured_N_positive":measured_n,
               "original_2018_UK_plant_rows_with_measured_P_positive":measured_p,
               "original_2018_UK_plant_rows_with_both_measured_positive":measured_both,
               "original_2018_UK_outlet_may_not_be_same_as_2018_EA_SA_samplepoint_or_PRTR_facility":True,
               "physical_outlet_2018_geography_if_exists_does_NOT_confirm_WWTP_source_kgP_or_BMWP_effect":True
           },
           "source_gates":{
               "UK_2018_report_explicitly_confirmed":bool(uk_reps),
               "historical_UK_outlet_records_exist":unique_reported_outlet_n>0,
               "historical_UK_outlet_geography_available":unique_coord_outlet_n>0,
               "EA_SA_2018_monitoring_point_to_UWWTD_native_2018_plant_outlet_crosswalk":False,
               "measured_UK_2018_P_and_N_both_positive_UWWTD_reports_sufficient":measured_both>=30,
               "true_2017_outlet_and_bio_2017_unconditionally_usable":False,
               "source_BMWP_causal_model_not_authorized":True
           }
        })
    except Exception as e:
        result.update({"status":"SOURCE_OR_ARCHIVE_STOP","error":type(e).__name__+": "+str(e)[:700]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EEA_2018_UK_OUTLET_SOURCE_GATE="+json.dumps({
       "status":result["status"],"error":result.get("error"),
       "UK_period_repcodes":result.get("source_report_period_summary",{}).get("UK_or_GB_repcodes_with_explicit_2018_report_year"),
       "plants":result.get("source_plants"),
       "outlets":result.get("source_outfalls"),
       "UK_source_QA":result.get("original_UK_or_GB_2018_outlet_QA"),
       "gates":result.get("source_gates")
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
