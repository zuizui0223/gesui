#!/usr/bin/env python3
"""Read JRC native river basin source schema, inventory and attribute-only pilot.

Strictly SOURCE-ONLY. Importantly GREEN point BasinID != soil basins HydroID
unless a later independent polygon overlay proves equivalence.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from datetime import datetime,timezone

SOIL="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer"
WATER="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
CONFIG={
 "soil_atmosN_2015":(SOIL,5,"Basin_3035_HydroID","MEAN_mgsqmy","Basin_3035_AreaSqKm"),
 "soil_mineralN_2010_19":(SOIL,10,"Basin_3035_HydroID","MEAN_kgha","Basin_3035_AreaSqKm"),
 "soil_organicN_2010_19":(SOIL,15,"Basin_3035_HydroID","MEAN_kgha","Basin_3035_AreaSqKm"),
 "soil_fixationN_2010_19":(SOIL,20,"Basin_3035_HydroID","MEAN_kgha","Basin_3035_AreaSqKm"),
 "water_pointN_2014_18":(WATER,5,"BasinID","PsN","AreaSqKm"),
}
HEADERS={"User-Agent":"gesui-native-JRC-riverbasin-point-vs-soilN-source-only/20261119",
         "Accept":"application/json"}
MAX_LAYER_OBJECTS=200000
SAMPLE_OBJECTS=25

def fetch(url, params):
    uri=url+"?"+urllib.parse.urlencode({"f":"json",**params})
    last=None
    for attempt in range(4):
        try:
            req=urllib.request.Request(uri,headers=HEADERS)
            with urllib.request.urlopen(req,timeout=70) as f:
                raw=f.read(15_000_000)
            result=json.loads(raw)
            if not isinstance(result,dict) or "error" in result:
                raise ValueError(f"ArcGIS original JSON error for {url}: {result}")
            return result,hashlib.sha256(raw).hexdigest()
        except (urllib.error.URLError,urllib.error.HTTPError,TimeoutError,ValueError) as e:
            last=str(e)
            if attempt<3:time.sleep(2*(attempt+1))
    raise RuntimeError(f"Official source access failed; not a zero source: {url}; {last}")

def validate_meta(meta,name,key,val,area):
    names={r.get("name") for r in meta.get("fields",[])}
    expected={"OBJECTID",key,val,area}
    if not expected.issubset(names) or meta.get("type")!="Feature Layer":
        raise ValueError(f"Original native JRC schema changed {name}: missing {expected-names}")
    if meta.get("geometryType")!="esriGeometryPolygon":
        raise ValueError(f"Original {name} layer is not polygon geometry")
    return sorted(names)

def summary_layer(name):
    base,lid,key,value,area=CONFIG[name]
    endpoint=f"{base}/{lid}"
    meta,hmeta=fetch(endpoint,{})
    fields=validate_meta(meta,name,key,value,area)
    count,hcount=fetch(endpoint+"/query",{"where":"1=1","returnCountOnly":"true"})
    ids,hids=fetch(endpoint+"/query",{"where":"1=1","returnIdsOnly":"true"})
    n=count.get("count")
    seq=ids.get("objectIds")
    if type(n) is not int or not(100<=n<=MAX_LAYER_OBJECTS):
        raise ValueError(f"Original {name} count missing or invalid: {n}")
    if not isinstance(seq,list) or len(seq)!=n or len(set(seq))!=n:
        raise ValueError(f"Original {name} objectIDs incomplete {n} vs {len(seq) if isinstance(seq,list) else type(seq)}")
    chosen=sorted(seq)[:SAMPLE_OBJECTS]
    data,hdata=fetch(endpoint+"/query",{
        "objectIds":",".join(map(str,chosen)),
        "outFields":",".join(["OBJECTID",key,value,area]),
        "returnGeometry":"false"
    })
    rows=[z.get("attributes",{}) for z in data.get("features",[])]
    if len(rows)!=len(chosen) or {r.get("OBJECTID") for r in rows}!=set(chosen):
        raise ValueError(f"Original {name} sample IDs incomplete")
    keys=[r.get(key) for r in rows]
    if any(type(z) is not int for z in keys) or len(set(keys))!=len(keys):
        raise ValueError(f"Original {name} sample native IDs missing/duplicated")
    sample=[]
    for r in rows:
        sample.append({key:r.get(key),value:r.get(value),area:r.get(area)})
    report={
      "name":name,"URL":endpoint,"native_geo_key_field":key,
      "value_field":value,"model_unit":"source-defined as in official map service",
      "official_original_feature_count":n,
      "sampled_original_native_distinct_key_n":len(set(keys)),
      "objectIDs_all_unique_verified":True,
      "original_feature_layer_fields":fields,
      "meta_sha256":hmeta,"count_sha256":hcount,"IDs_sha256":hids,
      "sample_attributes_sha256":hdata,
      "source_sample_values":sample,
      "full_source_attributes_NOT_YET_DOWNLOADED":True,
      "new_ecological_outcomes_NOT_OPENED":True
    }
    print("JRC_NATIVE_RIVERBASIN_SOURCE_Layer",name,
          json.dumps({"feature_count":n,"first5_sample":sample[:5],"native_key":key},ensure_ascii=False),flush=True)
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    result={}
    for name in CONFIG:
        result[name]=summary_layer(name)
    compatible={
      "source_native_schema_jrc_soil_all4_HydroID":all(
        result[n]["native_geo_key_field"]=="Basin_3035_HydroID" for n in CONFIG if n.startswith("soil_")),
      "source_native_schema_jrc_point_BasinID":result["water_pointN_2014_18"]["native_geo_key_field"]=="BasinID",
      "minimum_1000_native_soil_feature_rows_all_layers":all(
        result[n]["official_original_feature_count"]>=1000 for n in CONFIG if n.startswith("soil_")),
      "minimum_100_water_pointN_native_feature_rows":result["water_pointN_2014_18"]["official_original_feature_count"]>=100,
      "real_HydroID_cross_layer_full_overlap_NOT_YET_TESTED":True,
      "pointBasinID_vs_soilHydroID_spatial_overlap_NOT_TESTED":True,
      "original_2018_water_x_2017to19_land_joint_SOURCE_STOP_remains":True,
    }
    report={
      "schema":"gesui_JRC_20261119_native_riverbasin_soilN_water_pointN_schema_source_first_v1",
      "date_UTC":datetime.now(timezone.utc).isoformat(),
      "native_layers":result,
      "schema_and_inventory_support":compatible,
      "source_only_step1_inventory_GO_to_full_native_HydroID_census":all([
        compatible["source_native_schema_jrc_soil_all4_HydroID"],
        compatible["source_native_schema_jrc_point_BasinID"],
        compatible["minimum_1000_native_soil_feature_rows_all_layers"],
        compatible["minimum_100_water_pointN_native_feature_rows"]
      ]),
      "no_numeric_BQE_EQR_EIV_NUTRI_butterfly_HNS_opened":True,
      "do_not_assume_same_source_key":True,
    }
    (a.out/"JRC_NATIVE_RIVERBASIN_POINT_VS_SOILN_SOURCE_INVENTORY.json").write_text(
      json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("JRC_NATIVE_INVENTORY_GATE",json.dumps(compatible,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
