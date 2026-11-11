#!/usr/bin/env python3
"""NEW EEA WISE-1 WATERBASE EMISSIONS 2026 official source, year2023 metadata.

Separate emission database from UWWTD 2023 reporting exercise (which mainly
uses 2022 source-reference year). Assess whether *actual* reported 2023 N/P
emissions uniquely link to municipal WWTP FW outlet and monitored river.
A year included in overall 1977–2025 coverage is NOT source/site completeness.
No numerical benthic invertebrate EQR, any taxonomic outcomes or model read.
"""
from __future__ import annotations
import argparse,hashlib,json,re,urllib.request,urllib.error,html
from datetime import datetime,timezone
from pathlib import Path
from html.parser import HTMLParser

DOI="115940ff-51d1-42ac-ae14-699043b953ca"
ROOT="https://sdi.eea.europa.eu"
URLS={
 "official_EEA_Waterbase_emissions2026_catalogue_svr":ROOT+"/catalogue/srv/api/records/"+DOI,
 "official_EEA_Waterbase_emissions2026_datahub":ROOT+"/catalogue/datahub/api/records/"+DOI+"?language=all",
 "official_EEA_Waterbase_emissions2026_INPSIRE":ROOT+"/catalogue/inspire-pds/api/records/"+DOI+"?language=all",
 "official_public_datahub_page":"https://www.eea.europa.eu/en/datahub/datahubitem-view/5bd950b1-52b9-42f3-b0d7-bf51ff3056b7",
 "official_europa_open_dataset":"https://data.europa.eu/data/datasets/"+DOI+"?locale=en"
}
SIZE=850_000
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,t,a):
  if t in ("link","a","meta"):
   o=dict(a)
   u=o.get("href") or o.get("content")
   if u and any(z in str(u).lower() for z in ("download","data","webdav","datashare","emissions","dictionary")):
    self.links.append(str(u))
def fetch(url):
 req=urllib.request.Request(url,headers={"User-Agent":"gesui-WISE1-2026-emissions-2023-source-only/1.0","Accept":"text/html,application/json,application/xml,*/*"})
 try:
  with urllib.request.urlopen(req,timeout=70) as r:
   raw=r.read(SIZE)
   source=raw.decode("utf-8","replace")
   links=Links();links.feed(source)
   possible=re.findall(r"https?://[^\s\"'<>\\]+",source)
   urls=[]
   for v in links.links+possible:
    v=html.unescape(v).strip(");,.")
    if any(x in v.lower() for x in ("datashare","webdav","emissions","download","spatialite","sdi.eea")):
     urls.append(v[:550])
   return {"http":r.status,"final_url":r.url,"content_type":r.headers.get("Content-Type"),
      "content_length":r.headers.get("Content-Length"),"first850k_sha256":hashlib.sha256(raw).hexdigest(),
      "first_read_bytes":len(raw),"source_candidate_official_URLs":list(dict.fromkeys(urls))[:55],
      "value_2023_is_mentioned":("2023" in source),"source_2026_1977_2025_is_mentioned":("2025" in source)}
 except Exception as e:
  return {"url":url,"error":type(e).__name__+": "+str(e)[:260]}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 checks={k:fetch(v) for k,v in URLS.items()}
 out={
   "schema":"gesui_2023_original_WISE1_Emissions2026_public_1977to2025_SOURCE_ONLY_inventory_v01",
   "UTC":datetime.now(timezone.utc).isoformat(),
   "catalogue_EEA_2026_emissions_datahub":"https://www.eea.europa.eu/en/datahub/datahubitem-view/5bd950b1-52b9-42f3-b0d7-bf51ff3056b7",
   "original_EEA_2026_Emissions_1977to2025_dataset_identifier":"eea_t_waterbase-emissions_p_1977-2025_v01_r00",
   "source_2023_annual_observed_POINT_source_MUNICIPAL_NP_site_support_UNVERIFIED":True,
   "point_emission_country_only_are_not_2023_upstream_wetland_or_river_plant_discharge":True,
   "source_website_and_metadata_access_checks":checks,
   "no_EEA2023_BQE_numeric_outcomes_accessed":True,
   "positive_original_2018_and_2022_model_results_UNCHANGED":True
 }
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 for k,v in checks.items():
  print("EEA2026_EMISSIONS2023_SOURCE_META",k,json.dumps({z:y for z,y in v.items() if z!="source_candidate_official_URLs"}),flush=True)
  if "source_candidate_official_URLs" in v:
   print("EEA2026_EMISSIONS2023_PUBLIC_HREFS",k,json.dumps(v["source_candidate_official_URLs"][:24]),flush=True)
if __name__=="__main__":main()
