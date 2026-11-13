#!/usr/bin/env python3
"""Official Water Quality Portal NRSA 2023-24 biological source ACCESS probe.
Not a macro benthic dataset or negative ecological effect. Select four fixed
national NARS source sites before any WQP response, test exact identified
monitoring location / activity / project/year; never substitute nearby stations.

Use bounded WQP3 and legacy biological queries, query only a small
source-blind 2023+2024 site sample, and preserve explicit HTTP/status evidence.
"""
from __future__ import annotations
import collections,csv,hashlib,io,json,os,pathlib,time,urllib.error,urllib.parse,urllib.request
R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui2324_wqp_official"
R.mkdir(parents=True,exist_ok=True)
SITE="https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv"
SITE_SHA="f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013"
HOST="https://www.waterqualitydata.us"
# No mass station query or broad unfiltered biological result dump.
MIME="text/csv"
def sha(b):return hashlib.sha256(b).hexdigest()
def get(url,max_bytes=2*1024*1024):
 try:
  req=urllib.request.Request(url,headers={"User-Agent":"gesui-public-NRSA-WQP-biology-source-only/20261114",
                                         "Accept":"text/csv,application/json,*/*"})
  with urllib.request.urlopen(req,timeout=35) as r:
   typ=r.headers.get("Content-Type","")
   raw=r.read(max_bytes+1)
   if len(raw)>max_bytes:
    return {"http_status":r.status,"content_type":typ,
            "status":"TOO_LARGE_DO_NOT_INFER_OR_PARSE","max_bytes":max_bytes}
   hdr=raw[:256].decode("utf-8-sig",errors="replace")
   if "<html" in hdr.lower() or "<!doctype html" in hdr.lower():
    return {"http_status":r.status,"status":"HTML_NOT_DATA","first_text":hdr[:160]}
   rows=raw.decode("utf-8-sig",errors="replace").splitlines()
   rhead=(rows[0] if rows else "")
   return {"http_status":r.status,"content_type":typ,
           "sha256":sha(raw),"response_bytes":len(raw),
           "CSV_header":rhead[:800],"CSV_records_after_header":max(0,len(rows)-1),
           "sample_public_column_names":next(csv.reader([rhead]),[])[:30] if rhead else []}
 except Exception as e:
  return {"http_status":getattr(e,"code",None),"status":"SOURCE_UNAVAILABLE_OR_QUERY_UNSUPPORTED",
          "error":type(e).__name__+":"+str(e)[:240]}
def sites():
 p=R/"nrsa2324_sites.csv"
 req=urllib.request.Request(SITE,headers={"User-Agent":"gesui-EPA-WQP-research-source-check/20261114"})
 with urllib.request.urlopen(req,timeout=105) as r,p.open("wb") as w:
  while b:=r.read(1048576):w.write(b)
 if sha(p.read_bytes())!=SITE_SHA:raise RuntimeError("Native original EPA 2324 SiteInfo SHA changed")
 frozen=json.loads(pathlib.Path(os.environ["ORIGINAL_PRE_BIOFILM_UID_2129"]).read_text())
 z=frozen["2023_24"]["uid_bridge_manifest"]
 valid=set(z["water"].values())&set(z["biomass"].values())&set(z["chla"].values())
 if len(valid)!=2129:raise RuntimeError("Original 2023/24 source cohort changed")
 rows=[]
 with p.open(encoding="utf-8-sig",newline="") as f:
  for r in csv.DictReader(f):
   if r.get("UID") in valid and r.get("SITE_ID") and r.get("DATE_COL"):
    date=str(r["DATE_COL"])
    year=str(r.get("YEAR","")).strip()
    if year in ("2023","2024"):
     rows.append({"UID":r["UID"],"SITE_ID":r["SITE_ID"],"YEAR":year,"date":date})
 p.unlink()
 for y in ("2023","2024"):
  if sum(x["YEAR"]==y for x in rows)<2:raise RuntimeError("Source year missing in original EPA SiteInfo")
 out=[]
 for year in ("2023","2024"):
  entries=sorted((r for r in rows if r["YEAR"]==year),key=lambda x:sha(x["UID"].encode()))
  out.extend(entries[:2])
 if len({x["SITE_ID"] for x in out})!=4:
  raise RuntimeError("Original source fixed WQP test sites are not distinct")
 return out
def main():
 chosen=sites()
 report={"schema":"gesui_NR​​SA2023_24_GENUINE_WQP_SOURCE_ONLY_BIOLOGICAL_RETRIEVAL_TEST_v1",
         "source_EPA_2324_siteinfo_sha256":SITE_SHA,"native_frozen_source_2129_before_WQP_result_sign":True,
         "WQP_legacy_and_3point0_official_publisher_API_not_3rd_party_scrape":True,
         "genuine_sites_tested_n":len(chosen),"tested_queries":[],
         "sample_method":"2 native UID each 2023 and 2024; lowest SHA256 UID, predetermined before portal result",
         "no_inference_if_zero_or_nonmatching":True,"no_biological_numeric_effects_or_absence_claims":True}
 for s in chosen:
  site=s["SITE_ID"]
  for stid in (site, "EPA_NRSA-"+site):
   for path,other in [
     ("/data/Station/search",{"mimeType":"csv","providers":"STORET"}),
     ("/data/Result/search",{"mimeType":"csv","providers":"STORET","dataProfile":"biological",
                            "startDateLo":"2023-01-01","startDateHi":"2024-12-31"}),
     ("/wqx3/Station/search",{"mimeType":"csv"})]:
    args={"siteid":stid,**other}
    url=HOST+path+"?"+urllib.parse.urlencode(args)
    value=get(url)
    info={"native_original_SITE_ID":site,"native_UID_sha256":sha(s["UID"].encode()),
          "native_date":s["date"],"query_url":url,"source":value}
    report["tested_queries"].append(info)
    print("WQP_SOURCE_TEST",s["YEAR"],path,("original" if stid==site else "WQXprefixed"),
          "HTTP",value.get("http_status"),"count",value.get("CSV_records_after_header"),
          "status",value.get("status","DATA_OR_EMPTY"),flush=True)
 out=R/"NRSA2324_WQP_BIOLOGICAL_SOURCE_PROBE.json"
 out.write_text(json.dumps(report,sort_keys=True,indent=2))
 print("WQP_SOURCE_PROBE_COMPLETE",json.dumps({"n_sites":4,"queries":len(report["tested_queries"])}))
if __name__=="__main__":main()
