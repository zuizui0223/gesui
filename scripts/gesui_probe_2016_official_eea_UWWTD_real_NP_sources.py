#!/usr/bin/env python3
"""Find 2016 historical UWWTD actual *measured* discharge public source.

Probe only official EEA/europa catalogue URLs and public share metadata.
Do not read biological EQR numeric outcomes, do not replace missing 2016
mass by 2018/2020/2022 source-year reports, no pseudo-replication.
"""
from __future__ import annotations
import argparse,hashlib,json,re,urllib.request,urllib.parse,urllib.error
from datetime import datetime,timezone
from pathlib import Path
from html.parser import HTMLParser

YEAR=2016
SOURCE_ID="eea_t_uwwtd-waterbase_p_2015-2016_v06_r00"
PAGE="https://data.europa.eu/data/datasets/36115534-1c21-4c5e-9ec8-2ad22e5ffbb7"
URLS={
 "eea_2016_original_legacy_page":"https://www.eea.europa.eu/data-and-maps/data/waterbase-uwwtd-urban-waste-water-treatment-directive-6",
 "eea_2016_original_legacy_download":"https://www.eea.europa.eu/data-and-maps/data/waterbase-uwwtd-urban-waste-water-treatment-directive-6/at_download/file",
 "eea_2016_original_catalogue_hub":"https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715",
 "europa_2016_open_dataset":PAGE,
 "europa_2016_linked_turtle":PAGE+".ttl",
 "eea_2016_original_catalogue_path":"https://sdi.eea.europa.eu/data/"+SOURCE_ID,
}
UA="gesui-EEA-2016-third-independent-macro-original-source-availability/1.0"
MAX=1024*1024
class Links(HTMLParser):
 def __init__(self):
  super().__init__();self.found=[]
 def handle_starttag(self,tag,attrs):
  if tag in ("a","link","meta"):
   a=dict(attrs);u=a.get("href") or a.get("content")
   if u and any(z in str(u).lower() for z in ("download","2015-2016","2016","uwwtd","waterbase","tabular","dataset")):
    self.found.append(u)
def fetch(url,method):
 headers={"User-Agent":UA,"Accept":"text/html,application/json,text/turtle,application/zip,*/*"}
 if method=="GET":headers["Range"]="bytes=0-1048575"
 req=urllib.request.Request(url,headers=headers,method=method)
 try:
  with urllib.request.urlopen(req,timeout=48) as r:
   body=r.read(MAX if method=="GET" else 0)
   return {"http_status":r.status,"final_url":r.url,"content_type":r.headers.get("Content-Type"),
     "content_length_header":r.headers.get("Content-Length"),"first_1MB_SHA256":hashlib.sha256(body).hexdigest(),
     "content_starts_with_PK_ZIP":body.startswith(b"PK\x03\x04"),
     "retrieved_bytes":len(body),"body":body}
 except urllib.error.HTTPError as e:
  return {"http_status":e.code,"final_url":e.url,"error":str(e),"body":b""}
 except Exception as e:
  return {"final_url":url,"error":type(e).__name__+": "+str(e),"body":b""}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 results={}
 for key,url in URLS.items():
  rec=fetch(url,"GET")
  raw=rec.pop("body")
  if not rec.get("content_starts_with_PK_ZIP",False) and raw:
   decoded=raw.decode("utf-8","replace")
   l=Links();l.feed(decoded)
   rec["candidate_public_href_evidence"]=list(dict.fromkeys(l.found))[:45]
   rec["source_id_in_preview"]=SOURCE_ID in decoded
   rec["year2016_in_preview"]="2016" in decoded
   rec["preview_excerpt"]=decoded[:850] if len(decoded)<=2000 else ""
  results[key]=rec
  print("EEA2016_SOURCE_DISCOVERY",key,json.dumps({
    k:v for k,v in rec.items() if k not in ("candidate_public_href_evidence",)},ensure_ascii=False),flush=True)
 manifest={
  "schema":"gesui_EEA2016_3rd_benthic_year_OFFICIAL_source_access_and_period_definition_v01",
  "utc":datetime.now(timezone.utc).isoformat(),
  "2016_original_cycle_selected_BEFORE_2016_BQE_EQR_scores":True,
  "official_EEA_2015_2016_datahub_dataset_id":SOURCE_ID,
  "official_EU_dataset_listing":PAGE,
  "first_biology_score_2016_read":False,
  "any_reported_annual_mass_for_2016_accepted":False,
  "official_public_reference_access_probes":results,
  "older_public_2016_source_may_be_replaced_by_newer_uwwtd_reporting":True,
  "different_2018_2022_positive_benthic_prediction_results_UNCHANGED":True,
  "cannot_assume_2016_reporting_code_ensures_individual_sample_year":True
 }
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(manifest,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
 print("EEA2016_CYCLE_SOURCE_ONLY",json.dumps({"n_probed":len(results),
    "n_http200":sum(v.get("http_status")==200 for v in results.values()),
    "n_zip_preview":sum(v.get("content_starts_with_PK_ZIP") for v in results.values()),
    "NO_2016_EQR":True},indent=2),flush=True)
if __name__=="__main__":main()
