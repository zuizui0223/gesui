#!/usr/bin/env python3
"""Source-only check of whether SQL length or TOP/page syntax triggers HTTP403."""
import hashlib,json,os,pathlib,time,urllib.parse,urllib.request,urllib.error
T="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
P=[
 ("top1_site","SELECT TOP 1 * FROM "+T+" WHERE monitoringSiteIdentifier = 'FRGR04096745'"),
 ("top10_site","SELECT TOP 10 * FROM "+T+" WHERE monitoringSiteIdentifier = 'FRGR04096745'"),
 ("all_site","SELECT * FROM "+T+" WHERE monitoringSiteIdentifier = 'FRGR04096745'"),
 ("site_year","SELECT TOP 1 * FROM "+T+" WHERE monitoringSiteIdentifier = 'FRGR04096745' AND phenomenonTimeReferenceYear = 2018"),
 ("top1_tn","SELECT TOP 1 * FROM "+T+" WHERE observedPropertyDeterminandCode = 'EEA_31615-01-7'"),
 ("top1_tp","SELECT TOP 1 * FROM "+T+" WHERE observedPropertyDeterminandCode = 'CAS_7723-14-0'"),
 ("all_year","SELECT * FROM "+T+" WHERE phenomenonTimeReferenceYear = 2022"),
]
def get(sql):
    url="https://discodata.eea.europa.eu/sql?"+urllib.parse.urlencode({"query":sql,"p":1,"nrOfHits":10})
    try:
        q=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui-20261115-source-shortsql/1"})
        with urllib.request.urlopen(q,timeout=25) as r:raw=r.read(100000)
        data=json.loads(raw)
        rows=data.get("results") if isinstance(data,dict) else None
        return {"status":200,"n":len(rows) if isinstance(rows,list) else None,
            "first":rows[0] if rows else None,"raw_sha":hashlib.sha256(raw).hexdigest(),"url":url}
    except urllib.error.HTTPError as e:
        raw=e.read(500)
        return {"status":e.code,"error_sha":hashlib.sha256(raw).hexdigest(),"url":url}
    except Exception as e:return {"error":str(e)[:250],"url":url}
def main():
    out=[]
    for label,sql in P:
        x=get(sql);x["probe"]=label
        out.append(x)
        print("SHORTSQL",json.dumps({"probe":label,"status":x.get("status"),"n":x.get("n"),
                                     "error":x.get("error"),"first_id":(x.get("first") or {}).get("monitoringSiteIdentifier")}),flush=True)
        time.sleep(2)
    root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_shortsql"
    root.mkdir(parents=True,exist_ok=True)
    (root/"SHORT_SQL_ACCESS_AUDIT.json").write_text(json.dumps({
     "schema":"GESUI_WISE6_SHORT_SQL_ACCESS_V1","no_biotic_data":True,
     "no_missing_values_imputed":True,"requests":out},indent=2))
if __name__=="__main__":main()
