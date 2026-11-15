#!/usr/bin/env python3
"""Official EEA WISE SoE site-year and TOTAL N/P exact-code pilot.
This is source-quality work only, not an ecological outcome fit.
"""
import hashlib,json,os,pathlib,time,urllib.parse,urllib.request,urllib.error

TABLE="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
API="https://discodata.eea.europa.eu/sql"
QUERY_COLUMNS="monitoringSiteIdentifier,monitoringSiteIdentifierScheme,parameterWaterBodyCategory,observedPropertyDeterminandCode,observedPropertyDeterminandLabel,phenomenonTimeReferenceYear,procedureAnalysedMatrix,resultUom,resultNumberOfSamples,resultQualityMeanBelowLOQ,resultMeanValue"
PAIRS=[
 ("site_fr2018","monitoringSiteIdentifier = 'FRGR04096745'"),
 ("site_fr2018_year","monitoringSiteIdentifier = 'FRGR04096745' AND phenomenonTimeReferenceYear = 2018"),
 ("site_at2022","monitoringSiteIdentifier = 'ATFW10000177'"),
 ("site_at2022_year","monitoringSiteIdentifier = 'ATFW10000177' AND phenomenonTimeReferenceYear = 2022"),
 ("code_tn","observedPropertyDeterminandCode = 'EEA_31615-01-7'"),
 ("code_tp","observedPropertyDeterminandCode = 'CAS_7723-14-0'"),
 ("tn_year","observedPropertyDeterminandCode = 'EEA_31615-01-7' AND phenomenonTimeReferenceYear = 2022"),
]
def query(sql):
    u=API+"?"+urllib.parse.urlencode({"query":sql,"p":1,"nrOfHits":100})
    try:
        req=urllib.request.Request(u,headers={"User-Agent":"gesui-20261115-EEA-WISE6-site-source-pilot/1","Accept":"application/json"})
        with urllib.request.urlopen(req,timeout=30) as response: raw=response.read(700000)
        obj=json.loads(raw)
        rows=obj.get("results") if isinstance(obj,dict) else None
        if not isinstance(rows,list):return {"ok":False,"error":"no_results_array","raw_preview":str(obj)[:300],"url":u}
        # Only original EEA measured source fields, never biodiversity endpoints.
        return {"ok":True,"url":u,"response_sha256":hashlib.sha256(raw).hexdigest(),
                "source_row_count":len(rows),"all_rows":rows}
    except urllib.error.HTTPError as e:
        b=e.read(500)
        return {"ok":False,"status":e.code,"error_body_sha256":hashlib.sha256(b).hexdigest(),
                "error":"HTTPError","url":u}
    except Exception as e:
        return {"ok":False,"error":type(e).__name__+":"+str(e)[:210],"url":u}
def main():
    probes=[]
    for label,pred in PAIRS:
        sql="SELECT TOP 100 "+QUERY_COLUMNS+" FROM "+TABLE+" WHERE "+pred
        item=query(sql)
        item["probe_name"]=label
        if item["ok"]:
            rows=item.pop("all_rows")
            item["all_rows_for_audit"]=rows
            exact=[r for r in rows if r.get("parameterWaterBodyCategory")=="RW"
                  and r.get("procedureAnalysedMatrix") in ("W","W-DIS","W-TOT")
                  and r.get("phenomenonTimeReferenceYear") in (2018,2022)
                  and r.get("observedPropertyDeterminandCode") in ("EEA_31615-01-7","CAS_7723-14-0")]
            item["actual_TNTP_RW_2018_2022_in_first100"]=len(exact)
            item["example_source_ids"]=[str(r.get("monitoringSiteIdentifier")) for r in exact[:5]]
        probes.append(item)
        print("MEASURED_SOURCE",json.dumps({k:v for k,v in item.items() if k not in ("url","all_rows_for_audit")},sort_keys=True),flush=True)
        time.sleep(2)
    out=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_exact_source_pilot"
    out.mkdir(parents=True,exist_ok=True)
    (out/"SOURCE_PILOT_QA.json").write_text(json.dumps({
        "schema":"GESUI_20261115_WISE6_EXACT_CODE_SITE_2022_SOURCE_PILOT_V1",
        "nutrient_code_map":{"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"},
        "not_nitrate_or_orthophosphate":True,"no_ecological_outcome_examined":True,
        "no_unreported_as_zero":True,"query_receipts":probes},indent=2,ensure_ascii=False))
    print("END_SOURCE_PILOT",len(probes),flush=True)
if __name__=="__main__":main()
