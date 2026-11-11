#!/usr/bin/env python3
"""Official JRC GREEN 1990-2021 source inventory, no 2016 aquatic EQR.

JRC 2016 domestic/industrial nutrient point-source SHARE at European sea
river-basin outlets is a MODELLED BASIN OUTLET apportionment, not observed
2016 individual WWTP effluent at an upstream BQE monitoring station.
Collect official item metadata, optionally bounded source archive schema
when the original file is modest in size; do not inspect BQE outcomes.
"""
from __future__ import annotations
import argparse,hashlib,json,urllib.request,urllib.error
from datetime import datetime,timezone
from pathlib import Path

DATASET="https://data.jrc.ec.europa.eu/dataset/fc72adb8-e70e-4ad9-b1f5-1a6c712394a6"
ITEM="https://arcgis-maps.jrc.ec.europa.eu/portal/sharing/rest/content/items/8b105c72109243798913853ad9d05072"
DATA=ITEM+"/data"
HEADER={"Accept":"application/json,application/zip,*/*","User-Agent":"gesui-JRC-GREEN-2016-source-only-1990-2021/1.0"}
def ask(url,method,preview):
    req=urllib.request.Request(url,headers=HEADER,method=method)
    try:
        with urllib.request.urlopen(req,timeout=55) as r:
            raw=r.read(preview if method=="GET" else 0)
            return {"http_status":r.status,"final_url":r.url,"content_type":r.headers.get("Content-Type"),
                    "content_length":r.headers.get("Content-Length"),
                    "content_disposition":r.headers.get("Content-Disposition"),
                    "preview_raw_sha256":hashlib.sha256(raw).hexdigest(),
                    "preview_is_ZIP":raw.startswith(b"PK\x03\x04"),
                    "preview_first_bytes":raw[:60].hex(),
                    "preview_json":json.loads(raw) if raw.startswith(b"{") else None}
    except Exception as e:
        return {"error":type(e).__name__+": "+str(e)[:350],"requested_url":url}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    source={
      "schema":"gesui_2016_independentJRC_GREEN1990to2021_pointsource_basinOutflow_SOURCE_METADATA_ONLY_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "JRC_official_catalogue":DATASET,
      "JRC_official_original_research_data_DOI":"10.2905/JRC.GG4HJK8",
      "JRC_historical_2016_is_modelled_nutrient_apportionment_at_sea_basin_outlets":True,
      "JRC_modelled_point_sources_combine_domestic_and_industrial_not_POTW_effluent_only":True,
      "official_file_GDB_ARCGIS_sharing_source_URL":DATA,
      "source_access_schema_first":{
        "original_item_json":ask(ITEM+"?f=json","GET",500_000),
        "original_item_metadata":ask(ITEM+"/info/metadata/metadata.xml","HEAD",0),
        "original_gdb_data_HEAD":ask(DATA,"HEAD",0),
        "original_gdb_data_first_128KB":ask(DATA,"GET",131072),
      },
      "original_2016_EEA_BQE_values_read":False,
      "NO_2016_real_WWTP_NP_tonnes_per_yr_inferred":True,
      "cannot_interpret_basin_outlet_modeled_point_share_as_nutrient_delivered_to_each_upstream_BQE_station":True
    }
    # Avoid publishing huge licence/content preview; metadata may contain
    # contact persons, only retain safe item fields.
    o=source["source_access_schema_first"]["original_item_json"].get("preview_json")
    if isinstance(o,dict):
        safe={x:o.get(x) for x in ("id","title","size","type","typeKeywords","access","created","modified",
                                  "name","spatialReference","numViews","url","licenseInfo") if x in o}
        source["source_access_schema_first"]["original_item_json"]["preview_json"]=safe
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(source,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"JRC_item":source["source_access_schema_first"]["original_item_json"],
        "official_file":source["source_access_schema_first"]["original_gdb_data_HEAD"],
        "preview":source["source_access_schema_first"]["original_gdb_data_first_128KB"]},indent=2))
if __name__=="__main__":main()
