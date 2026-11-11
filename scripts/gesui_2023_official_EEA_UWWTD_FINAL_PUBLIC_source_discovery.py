#!/usr/bin/env python3
"""EEA official 2022-23 FINAL PUBLIC wastewater source metadata/source URL gate.

Year 2023 EEA Biology EQR site geometry is already frozen separately; never
read numeric EQR scores while finding plant-year measured TN/TP. Probe official
2026 published public metadata (not restricted internal), extract exact EEA
download/geo/ESRI service URL if available, do not assume 2022=2023 source.
"""
from __future__ import annotations
import argparse,hashlib,json,re,urllib.request,urllib.error,html
from datetime import datetime,timezone
from pathlib import Path
from html.parser import HTMLParser

TP_PUBLIC="88b3b08d-3b1b-4a65-9406-365b218e09da"
EEA="https://sdi.eea.europa.eu"
PUBLIC_ROOT="https://www.eea.europa.eu/en/datahub/datahubitem-view"
URLS={
 "public_treatment_inspire_json":f"{EEA}/catalogue/inspire-pds/api/records/{TP_PUBLIC}?language=all",
 "public_treatment_datahub_json":f"{EEA}/catalogue/datahub/api/records/{TP_PUBLIC}?language=all",
 "public_treatment_catalogue_api":f"{EEA}/catalogue/srv/api/records/{TP_PUBLIC}",
 "public_treatment_catalogue_xml":f"{EEA}/catalogue/srv/api/records/{TP_PUBLIC}/formatters/xml?approved=true",
 "public_treatment_eea_datahub":f"{PUBLIC_ROOT}/6ffd0af1-64eb-4cce-9970-0028e01ee64e",
 "public_discharge_eea_datahub":f"{PUBLIC_ROOT}/21874828-fa7a-4e7e-8a0a-52ec7d92f99f",
 "public_tabular_eea_datahub":f"{PUBLIC_ROOT}/6244937d-1c2c-47f5-bdf1-33ca01ff1715",
}
UA={"Accept":"application/json,text/html,application/xml,*/*","User-Agent":"gesui-2023-benthic-2023-sewage-public-final-source-only/1.0"}
MAX=600_000
PAT=re.compile(r'https?://[^\s<>"\'\\]+')
class Parser(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  vals=dict(attrs)
  if tag in ("a","link","meta"):
   u=vals.get("href") or vals.get("content")
   if u and ("http" in u or "download" in u):self.links.append(u)
def probe(url):
 try:
  req=urllib.request.Request(url,headers=UA)
  with urllib.request.urlopen(req,timeout=55) as r:
   raw=r.read(MAX)
   source=raw.decode("utf-8","replace")
   out={"HTTP":r.status,"final_url":r.url,"ContentType":r.headers.get("Content-Type"),
        "ContentLength":r.headers.get("Content-Length"),
        "first600kB_sha256":hashlib.sha256(raw).hexdigest(),
        "retrieved_initial_bytes":len(raw),
        "2023_and_FINAL_in_initial_source":("2023" in source and "FINAL" in source.upper()),
        "2022_and_FINAL_in_initial_source":("2022" in source and "FINAL" in source.upper())}
   candidate=PAT.findall(source)
   x=Parser();x.feed(source)
   candidate.extend(x.links)
   cleaned=[]
   for x in candidate:
    x=html.unescape(x.strip("\\u0027\\u0022)\\u003e;,."))
    if any(t in x.lower() for t in ("datashare","download","2022-2023","uwwtd","rest/services","datastore","geopackage",".gpkg")):
     cleaned.append(x[:650])
   out["official_public_source_URL_candidates_unverified"]=list(dict.fromkeys(cleaned))[:65]
   try:
    parsed=json.loads(raw)
    if isinstance(parsed,dict):
     out["json_top_keys"]=sorted(parsed.keys())[:100]
     out["catalogue_record_ID_matches_public_UUID"]=TP_PUBLIC in json.dumps(parsed)[:20000]
     out["public_catalogue_distribution_link_count"]=len(cleaned)
   except Exception:pass
   return out
 except Exception as e:
  return {"error":type(e).__name__+": "+str(e)[:200],"URL":url}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 out={
  "schema":"gesui_2023_BQE_2022_23_EEA_UWWTD_FINAL_PUBLIC_measured_NP_OFFICIAL_ACCESS_ONLY_v01",
  "UTC":datetime.now(timezone.utc).isoformat(),
  "final_EEA_UWWTD_2023_public_TP_metadata_UUID":TP_PUBLIC,
  "2022_23_public_source_UNLIKE_internal_restricted_release":True,
  "2023_source_year_only_NOT_assumed_from_2022_2023_catalogue_period":True,
  "source_2023_FRESHWATER_outfall_plus_plant_must_both_support_measured_annual_NP":True,
  "public_catalogue_probes":{k:probe(url) for k,url in URLS.items()},
  "original_EEA_BQE2023_numeric_EQR_outcomes_not_read":True,
  "2018_2022_successful_actual_NP_benthic_results_not_used_to_select_2023_treatment_source":True
 }
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 for k,v in out["public_catalogue_probes"].items():
  print("PUBLIC_EEA_UWWTD_FINAL2023_SOURCE",k,json.dumps({q:z for q,z in v.items() if q not in ("official_public_source_URL_candidates_unverified",)},sort_keys=True),flush=True)
  if v.get("official_public_source_URL_candidates_unverified"):
   print("ORIGINAL_2023_DOWNLOAD_CANDIDATES",k,json.dumps(v["official_public_source_URL_candidates_unverified"][:20]),flush=True)
if __name__=="__main__":main()
