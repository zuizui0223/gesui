#!/usr/bin/env python3
"""Probe EEA WISE6 individual-sample DisaggregatedData at original EQR station IDs.

Complements source-only annual AggregatedData probe. No terrestrial or aquatic
biological outcome values opened; 2018/2022 true TN/TP sample availability only.
"""
import collections
import hashlib
import json
import os
import pathlib
import time
import urllib.parse
import urllib.request
import urllib.error

ENDPOINT="https://discodata.eea.europa.eu/sql"
TABLES=[
 "[WISE_SOE].[latest].[Waterbase_T_WISE6_DisaggregatedData]",
 "[WISE_SOE].[v1r1].[Waterbase_T_WISE6_DisaggregatedData]"
]
SITES=[
 ("HU|HU101178391|euMonitoringSiteCode",2018),
 ("HR|HR29010|euMonitoringSiteCode",2022),
 ("IT|IT051096|euMonitoringSiteCode",2018),
 ("IT|IT03N008001091012LO1|euMonitoringSiteCode",2022),
]
TOTAL={"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"}
def query(sql,page=1,hits=150):
    url=ENDPOINT+"?"+urllib.parse.urlencode({"query":sql,"p":page,"nrOfHits":hits})
    try:
        req=urllib.request.Request(url,headers={"Accept":"application/json",
                 "User-Agent":"gesui-original-EEA-disaggregated-samples-site-year-precheck/20261115"})
        with urllib.request.urlopen(req,timeout=45) as r: raw=r.read(750000)
        j=json.loads(raw)
        data=j.get("results") if isinstance(j,dict) else None
        return {"ok":isinstance(data,list),"source_rows":data if isinstance(data,list) else [],
                "response_sha":hashlib.sha256(raw).hexdigest(),"url":url}
    except urllib.error.HTTPError as e:
        raw=e.read(512)
        return {"ok":False,"status":e.code,"error_digest":hashlib.sha256(raw).hexdigest(),"url":url,
                "source_rows":[]}
    except Exception as e:
        return {"ok":False,"error":type(e).__name__+":"+str(e)[:180],"url":url,
                "source_rows":[]}

def sampled_time(row):
    year=row.get("phenomenonTimeReferenceYear",row.get("phenomenonTimeSamplingDate_year"))
    if year is not None and str(year).strip():
        try: return int(year)
        except (ValueError,TypeError):pass
    date=row.get("phenomenonTimeSamplingDate","")
    if len(str(date))>=4 and str(date)[:4].isdigit():return int(str(date)[:4])
    return None

def diagnostic(site,year,rows):
    c,s,scheme=site.split("|")
    native=[r for r in rows if r.get("monitoringSiteIdentifier")==s
            and r.get("monitoringSiteIdentifierScheme")==scheme]
    byyear=collections.Counter(str(sampled_time(r)) for r in native)
    matched=[r for r in native if sampled_time(r)==year and r.get("parameterWaterBodyCategory")=="RW"
             and r.get("procedureAnalysedMatrix")=="W"]
    nutrient={}
    for label,code in TOTAL.items():
        original=[r for r in matched if r.get("observedPropertyDeterminandCode")==code]
        qual=[]
        for r in original:
            try:
                value=float(r.get("resultObservedValue"))
                unit=r.get("resultUom")
                badloq=str(r.get("resultQualityObservedValueBelowLOQ")).lower() in ("true","1")
                if value>=0 and not badloq and unit==("mg{N}/L" if label=="TN" else "mg{P}/L"):
                    qual.append(r)
            except (ValueError,TypeError):pass
        nutrient[label]={"source_raw_samples":len(original),"strict_numeric_samples":len(qual),
                          "native_source_units":dict(collections.Counter(str(r.get("resultUom")) for r in original))}
    return {"original_site":site,"original_target_year":year,
       "raw_source_rows_for_exact_native_site":len(native),
       "raw_year_counts":dict(byyear),
       "target_year_river_water_samples":len(matched),
       "original_NP_source":nutrient,
       "chemistry_source_eligibility_any_nutrient":any(nutrient[k]["strict_numeric_samples"]>=3 for k in TOTAL),
       "no_biological_outcome_examined":True}

def main():
    report={"schema":"GESUI_20261115_WISE6_DISAGGREGATED_EEA_SOURCE_SITE_PILOT_V1",
            "source_population":"4 original biological monitoring site IDs from original n3525",
            "studied_years":[2018,2022],
            "no_comparison_to_untreated_control":True,"no_ecological_outcomes_examined":True,
            "probes":[]}
    valid_table=None
    for table in TABLES:
        result=query("SELECT TOP 1 * FROM "+table)
        report["probes"].append({"type":"schema","table":table,
             "source_status":result["ok"],"response_sha":result.get("response_sha"),
             "fields":list(result["source_rows"][0]) if result["source_rows"] else [],
             "error_code":result.get("status"),"error":result.get("error")})
        print("DISAGG_SCHEMA",json.dumps(report["probes"][-1]),flush=True)
        if result["ok"] and len(result["source_rows"])==1:
            valid_table=table;break
    report["verified_table"]=valid_table
    if valid_table:
        for site,year in SITES:
            code=site.split("|")[1]
            if "'" in code: raise RuntimeError("Unsafe native ID")
            sql="SELECT * FROM "+valid_table+" WHERE monitoringSiteIdentifier = '"+code+"'"
            source=[];receipts=[];complete=False
            for p in range(1,5):
                r=query(sql,page=p,hits=150)
                receipts.append({"page":p,"ok":r["ok"],"source_sha":r.get("response_sha"),
                    "source_count":len(r["source_rows"]),"error_code":r.get("status")})
                if not r["ok"]:break
                source.extend(r["source_rows"])
                if len(r["source_rows"])<150:
                    complete=True;break
            data=diagnostic(site,year,source) if complete else None
            record={"type":"exact_original_site","station":site,"year":year,
                    "source_complete":complete,"n_source_rows":len(source),
                    "page_source_receipts":receipts,"source_diagnostic":data,
                    "source_truncated_must_not_be_zero":not complete}
            report["probes"].append(record)
            print("DISAGG_SITE",json.dumps({
              "id":site,"year":year,"complete":complete,"rows":len(source),
              "target_year_rows":data.get("target_year_river_water_samples") if data else None,
              "TN":data.get("original_NP_source",{}).get("TN") if data else None,
              "TP":data.get("original_NP_source",{}).get("TP") if data else None},sort_keys=True),flush=True)
            time.sleep(1)
    report["year_match_count"]=sum(bool(x.get("source_diagnostic") and
        x["source_diagnostic"].get("target_year_river_water_samples",0)>0)
        for x in report["probes"] if x["type"]=="exact_original_site")
    root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_disagg_pilot"
    root.mkdir(parents=True,exist_ok=True)
    (root/"DISAGGREGATED_EXACT_SITE_SOURCE_GATE.json").write_text(json.dumps(report,indent=2,ensure_ascii=False))
    print("DISAGG_SOURCE_END",json.dumps({"table":valid_table,"target_year_sites_with_samples":report["year_match_count"],
              "original_sites_tested":len(SITES)}),flush=True)
if __name__=="__main__":main()
