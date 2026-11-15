#!/usr/bin/env python3
"""Source-only official EEA WISE SoE Water Quality API preflight.
Do not open ecological outcome fields or infer zero sewage from absent reports.
"""
import csv, hashlib, json, os, pathlib, urllib.parse, urllib.request

DATASETS=(
    "[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]",
    "[WISE_SOE].[v1].[Waterbase_T_WISE6_AggregatedData]",
    "[WISE_SOE].[v2].[Waterbase_T_WISE6_AggregatedData]",
    "[WISE_SOE].[v1r1].[Waterbase_T_WISE6_AggregatedData]"
)

def source_ids():
    base=pathlib.Path(os.environ["ORIGINAL_H7_454_DIR"])
    with (base/"EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv").open() as f:
        cohort=list(csv.DictReader(f))
    keys={(r["HYBAS_ID_L07"],r["PERIOD"]) for r in cohort}
    assert len(cohort)==454 and len(keys)==454 and len({x[0] for x in keys})==232
    path=pathlib.Path(os.environ["ORIGINAL_EEA_SITE_H7_MAP_DIR"])/"EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
    source=json.loads(path.read_bytes())
    rows=[r for r in source["candidate_source_site_period_rows"]
          if (str(r["HYBAS_ID_L07"]),r["period"]) in keys
          and r["panel"] in ("EEA18","EEA22")]
    if len(rows)!=3525 or len({r["physical_id"] for r in rows})!=2485:
        raise RuntimeError("Source-original 3525 station-years / 2485 stations changed")
    for r in rows:
        code=r["native_id"].split("|")
        if len(code)!=3 or code[0]!=r["country"] or not all(code):
            raise RuntimeError("Invalid original WISE monitoring identifier")
    return rows

def probe(query):
    url="https://discodata.eea.europa.eu/sql?"+urllib.parse.urlencode(
        {"query":query,"p":1,"nrOfHits":3})
    try:
        req=urllib.request.Request(url,headers={
          "Accept":"application/json","User-Agent":"gesui-source-only-EEA-WISE6/20261115"})
        with urllib.request.urlopen(req,timeout=32) as resp:
            raw=resp.read(65536)
            status=resp.status
        j=json.loads(raw)
        values=j.get("results") if isinstance(j,dict) else None
        return {"ok":isinstance(values,list),"status":status,"url":url,
            "source_sha256":hashlib.sha256(raw).hexdigest(),
            "raw_prefix":raw[:500].decode("utf-8","replace"),
            "keys":list(j) if isinstance(j,dict) else [],
            "rows":len(values) if isinstance(values,list) else None,
            "column_names":list(values[0]) if values and isinstance(values[0],dict) else [],
            "first_row_example":values[0] if values and isinstance(values[0],dict) else None}
    except Exception as e:
        return {"ok":False,"error":type(e).__name__+":"+str(e)[:250],"url":url}

def main():
    sites=source_ids()
    output={"schema":"GESUI_EEA_WISE6_2026_SAME_SITE_SOURCE_ONLY_PRECHECK_V1",
            "original_study_rows":454,"original_site_year_rows":3525,
            "unique_original_sites":2485,"no_ecological_values_opened":True,
            "no_missing_chemistry_imputed":True,"probes":[]}
    chosen=None
    for table in DATASETS:
        x=probe("SELECT TOP 1 * FROM "+table)
        x["table"]=table
        output["probes"].append(x)
        print("TABLE",json.dumps({"table":table,"ok":x["ok"],
              "http":x.get("status"),"keys":x.get("keys"),
              "rows":x.get("rows"),"columns":x.get("column_names"),"error":x.get("error")}),flush=True)
        if x["ok"] and x.get("rows"):
            chosen=table
            break
    output["verified_table"]=chosen
    if chosen:
        ids=sorted({r["physical_id"] for r in sites})
        for native in (ids[0],ids[len(ids)//2],ids[-1]):
            country,station,scheme=native.split("|")
            if any("'" in x for x in (station,scheme)):
                raise RuntimeError("Unsafe native station code")
            sql=("SELECT TOP 3 monitoringSiteIdentifier, monitoringSiteIdentifierScheme,"
                 " parameterWaterBodyCategory, observedPropertyDeterminandCode,"
                 " phenomenonTimeReferenceYear, resultMeanValue, resultUom,"
                 " resultNumberOfSamples, resultQualityMeanBelowLOQ FROM "+chosen+
                 " WHERE monitoringSiteIdentifier = '"+station+"'"
                 " AND monitoringSiteIdentifierScheme = '"+scheme+"'"
                 " AND phenomenonTimeReferenceYear IN (2018,2022)")
            x=probe(sql)
            x["source_native_id"]=native
            output["probes"].append(x)
            print("STATION",json.dumps({"id":native,"ok":x["ok"],"rows":x.get("rows"),
                   "error":x.get("error")}),flush=True)
    out=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_source_probe"
    out.mkdir(parents=True,exist_ok=True)
    (out/"WISE6_2026_SOURCE_PREFLIGHT.json").write_text(json.dumps(output,indent=2,sort_keys=True))
    print("STUDY_SOURCE_GATE",json.dumps({"study":454,"sites":2485,"station_year":3525,
                                         "table":chosen,"source_only":True}),flush=True)
    if not chosen:raise SystemExit("WISE SoE API source access unavailable: no chemistry join/fit.")
if __name__=="__main__":
    main()
