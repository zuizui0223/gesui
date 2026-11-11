#!/usr/bin/env python3
"""Official EEA Waterbase Emissions 2026 WISE-1 source/grain field dictionary.

Official downloadable CSV dataset definition for emission values and 
geographic monitoring object. Does NOT inspect numeric emission data, does 
not assume a facility-level point source or measured WWTP N/P site exposure.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,urllib.request
from pathlib import Path
from datetime import datetime,timezone

ID="115940ff-51d1-42ac-ae14-699043b953ca"
B="https://sdi.eea.europa.eu/catalogue/api/records/"+ID+"/attachments/"
FILES={
 "Emissions":"Waterbase_v2022_1_dataset_definition_T_WISE1_Emissions.csv",
 "SpatialObject":"Waterbase_v2022_1_dataset_definition_S_WISE1_SpatialObject_DerivedData.csv",
 "RiverineInputLoads":"Waterbase_v2022_1_dataset_definition_T_WISE1_RiverineInputLoads.csv"
}
def fetch(url):
 req=urllib.request.Request(url,headers={"User-Agent":"gesui-2023-WISE1-2026-source-fields-only","Accept":"text/csv,*/*"})
 with urllib.request.urlopen(req,timeout=48) as r:
  raw=r.read(1_500_000)
  if len(raw)>=1_500_000:raise ValueError("Original official EEA field definition bigger than 1.5MB")
  if not raw:raise ValueError("Original official EEA field definition source empty")
  ct=r.headers.get("Content-Type","")
  if "html" in ct.lower() or raw[:20].lower().startswith(b"<html"):raise ValueError("EEA source is HTML not real CSV")
  decoded=raw.decode("utf-8-sig","replace")
  head=decoded.splitlines()[0]
  sep=max([",",";","\t"],key=head.count)
  rows=list(csv.DictReader(io.StringIO(decoded),delimiter=sep))
  if not rows:raise ValueError("EEA field definition parsed zero metadata rows")
  return {"source_URL":r.url,"status":r.status,"content_type":ct,
          "original_bytes":len(raw),"source_SHA256":hashlib.sha256(raw).hexdigest(),
          "delimiter":sep,"metadata_columns":list(rows[0].keys()),
          "all_field_definition_source_records":rows[:180],
          "record_count":len(rows)}
def main():
 a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True);opts=a.parse_args()
 out={"schema":"gesui_EEA_WaterbaseEmissions2026_2023_ORIGINAL_WISE1_schema_field_dictionary_SOURCE_ONLY_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "reference_year_2023_NOT_INFERRED_from_1977_2025_dataset_time_range":True,
      "actual_WWTP_P_N_load_station_site_coverage_NOT_ESTABLISHED":True,
      "EEA_2026_original_WISE1_field_definitions":{},
      "all_2023_BQE_biological_scores_NOT_OPENED":True}
 for key,p in FILES.items():
  v=fetch(B+p)
  out["EEA_2026_original_WISE1_field_definitions"][key]=v
  print("ORIGINAL_EEA_EMISSIONS_2026_SCHEMAS",key,"source",v["source_SHA256"],"rows",v["record_count"],"header",v["metadata_columns"],"preview",json.dumps(v["all_field_definition_source_records"][:6]),flush=True)
 opts.out.parent.mkdir(parents=True,exist_ok=True)
 opts.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()
