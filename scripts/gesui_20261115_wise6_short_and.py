#!/usr/bin/env python3
"""Short same-site AND nutrient exact-code API gates; source-only."""
import hashlib,json,os,pathlib,urllib.parse,urllib.request,urllib.error,time
T="[WISE_SOE].[latest].[Waterbase_T_WISE6_AggregatedData]"
Q=[
 ("known_tn_site","SELECT * FROM "+T+" WHERE monitoringSiteIdentifier = 'BERW_04800'"),
 ("known_tn_site_code","SELECT * FROM "+T+" WHERE monitoringSiteIdentifier = 'BERW_04800' AND observedPropertyDeterminandCode = 'EEA_31615-01-7'"),
 ("known_tp_site","SELECT * FROM "+T+" WHERE monitoringSiteIdentifier = 'FISW_1542'"),
 ("known_tp_site_code","SELECT * FROM "+T+" WHERE monitoringSiteIdentifier = 'FISW_1542' AND observedPropertyDeterminandCode = 'CAS_7723-14-0'"),
 ("year2022_tn","SELECT * FROM "+T+" WHERE phenomenonTimeReferenceYear = 2022 AND observedPropertyDeterminandCode = 'EEA_31615-01-7'"),
]
def go(sql):
    u="https://discodata.eea.europa.eu/sql?"+urllib.parse.urlencode({"query":sql,"p":1,"nrOfHits":100})
    try:
        req=urllib.request.Request(u,headers={"Accept":"application/json","User-Agent":"gesui-WISE6-short-AND-20261115/1"})
        with urllib.request.urlopen(req,timeout=32) as r:b=r.read(300000)
        obj=json.loads(b);rows=obj.get("results",[])
        return {"ok":isinstance(rows,list),"n":len(rows),"rows":rows[:100],
            "response_sha":hashlib.sha256(b).hexdigest(),"url":u}
    except urllib.error.HTTPError as e:
        b=e.read(512);return {"ok":False,"http_status":e.code,"error_digest":hashlib.sha256(b).hexdigest(),"url":u}
    except Exception as e:return {"ok":False,"error":str(e)[:250],"url":u}
def main():
    ret=[]
    for name,q in Q:
        x=go(q);x["name"]=name;ret.append(x)
        print("AND_PROBE",json.dumps({"name":name,"ok":x["ok"],"n":x.get("n"),
                                      "status":x.get("http_status"),"error":x.get("error")}),flush=True)
        time.sleep(2)
    p=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_wise6_and_probe"
    p.mkdir(parents=True,exist_ok=True)
    (p/"SQL_AND_SOURCE_QA.json").write_text(json.dumps({
        "schema":"GESUI_20261115_WISE6_SHORT_AND_SOURCE_V1",
        "no_ecology":True,"probes":ret},indent=2))
if __name__=="__main__":main()
