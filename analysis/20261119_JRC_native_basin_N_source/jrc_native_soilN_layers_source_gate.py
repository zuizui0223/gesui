#!/usr/bin/env python3
"""Read JRC native river-basin soil N source layers, without biota or EQR values.

GREEN point source basin-ID namespace is deliberately NOT crosswalked to
soil HydroID; the two are separate GIS products until polygons are aligned.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE_SOIL="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer"
BASE_WATER="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
SOURCES={
 "soil_atmos_N_mg_m2_2015":(BASE_SOIL,5,"MEAN_mgsqmy","soil"),
 "soil_mineral_N_kg_ha_2010_19":(BASE_SOIL,10,"MEAN_kgha","soil"),
 "soil_organic_N_kg_ha_2010_19":(BASE_SOIL,15,"MEAN_kgha","soil"),
 "soil_fixation_N_kg_ha_2010_19":(BASE_SOIL,20,"MEAN_kgha","soil"),
 "water_model_point_N_t_2014_18":(BASE_WATER,5,"PsN","point"),
}
PAGE=400
MIN_COMPLETE_SOIL_BASINS=100
MAX_SOURCE_OBJECTS=50000

def request(url,params):
    query=urllib.parse.urlencode({"f":"json",**params})
    final=url+"?"+query
    last=None
    for retry in range(4):
        try:
            r=urllib.request.Request(final,headers={
              "User-Agent":"GESUI-20161119-academic-basin-source-only/1.0",
              "Accept":"application/json"
            })
            with urllib.request.urlopen(r,timeout=90) as inp:
                payload=inp.read(12_000_000)
            obj=json.loads(payload)
            if not isinstance(obj,dict) or obj.get("error"):
                raise ValueError("Official source returned ArcGIS JSON error: "+str(obj.get("error")))
            return obj,hashlib.sha256(payload).hexdigest()
        except (OSError,ValueError,TimeoutError) as e:
            last=str(e)
            if retry!=3:time.sleep(retry+1)
    raise RuntimeError("Unavailable original JRC layer "+url+": "+str(last))

def validate_id(value):
    if isinstance(value,bool) or not isinstance(value,int):
        raise ValueError("Original spatial basin ID must be an integer")
    return value

def nonnegative(value):
    try:v=float(value)
    except (TypeError,ValueError):return None
    return v if math.isfinite(v) and v>=0 else None

def audit_layer(label,base,lid,value_field,grain):
    url=f"{base}/{lid}"
    meta,mhash=request(url,{})
    fields={str(x.get("name")) for x in meta.get("fields",[])}
    if meta.get("type")!="Feature Layer" or meta.get("geometryType")!="esriGeometryPolygon":
        raise ValueError("Not a genuine polygon feature layer: "+label)
    geomid="Basin_3035_HydroID" if "Basin_3035_HydroID" in fields else ("HydroID" if "HydroID" in fields else None)
    if grain=="point":geomid="BasinID" if "BasinID" in fields else None
    oid=str(meta.get("objectIdField") or "OBJECTID")
    if geomid is None or oid not in fields or value_field not in fields:
        raise ValueError("Original schema lacks basin ID/value: "+label+" "+str(sorted(fields)[:50]))
    countobj,chash=request(url+"/query",{"where":"1=1","returnCountOnly":"true"})
    idobj,ihash=request(url+"/query",{"where":"1=1","returnIdsOnly":"true"})
    count=countobj.get("count")
    ids=idobj.get("objectIds")
    if type(count) is not int or not (1<=count<=MAX_SOURCE_OBJECTS):
        raise ValueError("Original source count absent / too large "+label+" "+str(count))
    if not isinstance(ids,list) or len(ids)!=count or len(set(ids))!=count:
        raise ValueError("Original IDs and count do not match "+label)
    ids=sorted(validate_id(v) for v in ids)
    ids_covered=set()
    records={}
    page_hashes=[]
    output_fields=sorted({oid,geomid,value_field})
    for k in range(0,len(ids),PAGE):
        batch=ids[k:k+PAGE]
        obj,hashpage=request(url+"/query",{
            "objectIds":",".join(map(str,batch)),
            "outFields":",".join(output_fields),
            "returnGeometry":"false"
        })
        feature=obj.get("features",[])
        if len(feature)!=len(batch):
            raise ValueError("ArcGIS truncated original query page "+label+" at "+str(k))
        got={a.get("attributes",{}).get(oid) for a in feature}
        if got!=set(batch):
            raise ValueError("Original JRC OBJECTID page IDs drift "+label)
        page_hashes.append(hashpage)
        for ent in feature:
            attrs=ent["attributes"]
            oidv=validate_id(attrs.get(oid))
            bid=attrs.get(geomid)
            if bid is None or isinstance(bid,bool) or isinstance(bid,(dict,list)):
                raise ValueError("Missing/invalid original native basin ID "+label)
            key=str(bid)
            if key in records:
                raise ValueError("Original JRC basin ID is duplicated "+label+" key "+key)
            if oidv in ids_covered:raise ValueError("OBJECTID repeated in pages")
            ids_covered.add(oidv)
            records[key]=nonnegative(attrs.get(value_field))
        time.sleep(.025)
    if len(ids_covered)!=count:
        raise ValueError("Not all original ArcGIS feature OIDs delivered "+label)
    values=list(records.values())
    report={
      "layer":lid,"url":url,"original_name":meta.get("name"),
      "original_feature_geometry_type":meta.get("geometryType"),
      "source_id_field":geomid,"source_value_field":value_field,
      "source_n_feature":count,"source_n_spatial_id":len(records),
      "number_missing_or_invalid_numeric":sum(x is None for x in values),
      "number_true_zero":sum(x==0 for x in values),
      "number_positive":sum(x is not None and x>0 for x in values),
      "meta_response_sha256":mhash,"count_response_sha256":chash,
      "id_inventory_response_sha256":ihash,
      "attribute_page_count":len(page_hashes),
      "attribute_page_sha256_of_hashlist":hashlib.sha256("\n".join(page_hashes).encode()).hexdigest(),
      "all_source_field_names":sorted(fields),
      "original_metadata_type":meta.get("type"),
      "spatial_ID_is_not_assumed_to_match_HydroBASINS_H7":True
    }
    return records,report

def source_gate(records):
    preferred=["soil_atmos_N_mg_m2_2015","soil_mineral_N_kg_ha_2010_19","soil_organic_N_kg_ha_2010_19"]
    full=sorted(set.intersection(*(set(records[k]) for k in preferred)))
    all4=sorted(set.intersection(*(set(records[k]) for k in list(SOURCES)[:4])))
    common=[k for k in full if all(records[v][k] is not None for v in preferred)]
    nonzero={v:sum(records[v][k]>0 for k in common) for v in preferred}
    varied={v:len({records[v][k] for k in common})>1 for v in preferred}
    checks={
        "source_only_complete_same_native_soil_basin_at_least100":len(common)>=MIN_COMPLETE_SOIL_BASINS,
        "source_only_nonzero_soil_basins_at_least30_each":all(z>=30 for z in nonzero.values()),
        "all_three_soil_model_variables_vary":all(varied.values())
    }
    return {
       "native_soil_basins_3layers_key_intersection":len(full),
       "native_soil_basins_all4_key_intersection":len(all4),
       "native_soil_basins_3layers_complete_nonmissing":len(common),
       "positive_native_basin_per_variable":nonzero,
       "native_value_varies":varied,"hard_gates":checks,
       "soil_source_first_native_basin_pass":all(checks.values()),
       "water_point_N_basins_unique":len(records["water_model_point_N_t_2014_18"]),
       "soil_source_hydroID_not_equated_to_point_source_BasinID":True,
       "soil_source_hydroID_not_equated_to_HydroBASINS_H7":True,
       "biotic_2x2_source_STOP_remains":True,
       "no_EEA_EOL_or_EIV_HNS_numerical_outcomes_read":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    values={}
    receipts={}
    for key,(base,lid,field,kind) in SOURCES.items():
        print("INSPECT_JRC_NATIVE_SOURCE",key,flush=True)
        rows,report=audit_layer(key,base,lid,field,kind)
        values[key]=rows
        receipts[key]=report
        print("SOURCE_SUMMARY",key,report["source_n_spatial_id"],report["number_missing_or_invalid_numeric"],flush=True)
    gate=source_gate(values)
    summary={
      "schema":"gesui_20261119_JRC_NATIVE_SOIL_BASINS_0p1deg_1km_model_inputs_SOURCE_FIRST_V1",
      "utc":datetime.now(timezone.utc).isoformat(),
      "jrc_official_soil_layers":receipts,
      "source_only_gate":gate,
      "soil_model_period":"EMEP2015 for deposition; DAYCENT2010-2019 for fertilizer and fixation",
      "water_model_period":"JRC GREEN 2014-2018 mean point N domestic plus industrial",
      "raster_original_numeric_pixels_NOT_DOWNLOADED":True,
      "original_HydroBASINS_H7_crosswalk_NOT_DONE":True,
      "do_not_make_biotic_prediction":True
    }
    (args.out/"JRC_NATIVE_BASIN_SOURCE_ONLY_QA.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    allids=sorted(set.union(*(set(z) for z in values.values())))
    with (args.out/"JRC_NATIVE_BASIN_SOURCE_VALUES_UNMATCHED_GEOS.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f)
        writer.writerow(["native_numeric_ID_UNMATCHED","pointN_MODEL_BasinID_value_not_same_as_soil_HydroID",*list(SOURCES)[:4]])
        for k in allids:
            writer.writerow([k,values["water_model_point_N_t_2014_18"].get(k),
                *[values[v].get(k) for v in list(SOURCES)[:4]]])
    print("GESUI_JRC_NATIVE_RIVERBASIN_SOURCE_ONLY",json.dumps(gate,sort_keys=True),flush=True)

if __name__=="__main__":main()
