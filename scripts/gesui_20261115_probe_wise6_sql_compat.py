#!/usr/bin/env python3
"""Official EEA WISE SoE Discodata SQL compatibility probe (source-only).
Do not open original biological endpoints or infer zero chemistry on HTTP403.
"""
import hashlib,json,os,pathlib,time,urllib.parse,urllib.request,urllib.error

TABLE="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
BASE="https://discodata.eea.europa.eu/sql"
QUERIES=(
 ("repeat_baseline","SELECT TOP 1 * FROM "+TABLE),
 ("projection","SELECT TOP 1 monitoringSiteIdentifier FROM "+TABLE),
 ("year_numeric","SELECT TOP 1 * FROM "+TABLE+" WHERE phenomenonTimeReferenceYear = 2022"),
 ("country_literal","SELECT TOP 1 * FROM "+TABLE+" WHERE countryCode = 'FR'"),
 ("existing_native_site","SELECT TOP 1 * FROM "+TABLE+" WHERE monitoringSiteIdentifier = 'FRGR04067200'"),
 ("distinct_year_recent","SELECT TOP 1 phenomenonTimeReferenceYear FROM "+TABLE+" ORDER BY phenomenonTimeReferenceYear DESC"),
)
def request(q):
    u=BASE+"?"+urllib.parse.urlencode({"query":q,"p":1,"nrOfHits":1})
    headers={"User-Agent":"gesui-EEA-WISE6-SQL-source-gate-20261115/1","Accept":"application/json"}
    req=urllib.request.Request(u,headers=headers)
    try:
        with urllib.request.urlopen(req,timeout=29) as resp:
            raw=resp.read(60000)
            code=resp.status
        payload=json.loads(raw)
        vals=payload.get("results") if isinstance(payload,dict) else None
        return {"url":u,"status":code,"is_record_response":isinstance(vals,list),
                "row_count":len(vals) if isinstance(vals,list) else None,
                "result_sha256":hashlib.sha256(raw).hexdigest(),
                "result_preview":vals[:1] if isinstance(vals,list) else str(payload)[:240]}
    except urllib.error.HTTPError as e:
        raw=e.read(512)
        return {"url":u,"status":e.code,"error":"HTTPError",
                "error_sha256":hashlib.sha256(raw).hexdigest(),"error_preview":raw.decode("utf-8","replace")[:300]}
    except Exception as e:
        return {"url":u,"error":type(e).__name__+":"+str(e)[:250]}
def main():
    out=[]
    for name,q in QUERIES:
        time.sleep(3)
        item=request(q);item["probe_name"]=name
        out.append(item)
        print("WISE6_QUERY",json.dumps({k:v for k,v in item.items() if k not in ("url","result_preview","error_preview")},sort_keys=True),flush=True)
    dest=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_query_syntax"
    dest.mkdir(parents=True,exist_ok=True)
    obj={"schema":"GESUI_WISE_SOE_DISCODATA_QUERY_COMPATIBILITY_V1",
         "data_queried":"Water Quality ICM AggregatedData only",
         "pre_existing_source_exact_site_support":3525,
         "no_ecological_outcomes_opened":True,
         "no_water_quality_nonmatches_equated_to_zero":True,"query_receipts":out}
    (dest/"QUERY_COMPATIBILITY_AUDIT.json").write_text(json.dumps(obj,sort_keys=True,indent=2))
    print("WISE6_QA",json.dumps({"status_counts":{str(s):sum(a.get("status")==s for a in out) for s in set(a.get("status") for a in out)},"queries":len(out)}),flush=True)
if __name__=="__main__":main()
