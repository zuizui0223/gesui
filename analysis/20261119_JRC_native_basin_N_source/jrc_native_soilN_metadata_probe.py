#!/usr/bin/env python3
"""Metadata and deterministic first-record probe for very large JRC native layers.

No biological values, no full 74k catchment extraction and no H7 joins.
"""
import argparse,hashlib,json,time,urllib.request,urllib.parse
from pathlib import Path
from datetime import datetime,timezone

BASE="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer"
POINT="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
LAYERS={
 "atmos2015_soil_native":(BASE,5,"MEAN_mgsqmy","Basin_3035_HydroID"),
 "mineralN_soil_native":(BASE,10,"MEAN_kgha","Basin_3035_HydroID"),
 "organicN_soil_native":(BASE,15,"MEAN_kgha","Basin_3035_HydroID"),
 "fixationN_soil_native":(BASE,20,"MEAN_kgha","Basin_3035_HydroID"),
 "pointN_water_native":(POINT,5,"PsN","BasinID")
}
def get(url,parameters):
    url=url+"?"+urllib.parse.urlencode({"f":"json",**parameters})
    last=None
    for tries in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-native-JRC-source-geometry-probe-20261119"})
            with urllib.request.urlopen(req,timeout=60) as f:raw=f.read(5000000)
            obj=json.loads(raw)
            if not isinstance(obj,dict) or obj.get("error"):raise ValueError(str(obj.get("error")))
            return obj,hashlib.sha256(raw).hexdigest()
        except (OSError,ValueError,TimeoutError) as e:
            last=e
            if tries<3:time.sleep(tries+1)
    raise RuntimeError("Original JRC service request unavailable: "+str(last))
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    result={}
    for key,(base,lid,val,kfield) in LAYERS.items():
        url=f"{base}/{lid}"
        meta,metasha=get(url,{})
        fields={x.get("name") for x in meta.get("fields",[])}
        count,csha=get(url+"/query",{"where":"1=1","returnCountOnly":"true"})
        n=count.get("count")
        if not isinstance(n,int) or n<=0:raise ValueError("No original source count "+key)
        orig=meta.get("objectIdField") or "OBJECTID"
        if (not orig in fields or val not in fields or kfield not in fields or
            meta.get("type")!="Feature Layer" or meta.get("geometryType")!="esriGeometryPolygon"):
            raise ValueError("JRC official polygon source schema mismatch: "+key+
                             " missing "+str(set((orig,val,kfield))-fields))
        sample,sig=get(url+"/query",{
            "where":"1=1","outFields":",".join((orig,kfield,val)),"returnGeometry":"false",
            "orderByFields":orig+" ASC","resultRecordCount":3
        })
        if len(sample.get("features",[]))!=min(3,n):
            raise ValueError("JRC first-source-IDs page truncated "+key)
        attrs=[x["attributes"] for x in sample["features"]]
        if len({x[orig] for x in attrs})!=len(attrs):raise ValueError("Duplicated native OID sample")
        result[key]={
          "original_layer":url,"count":n,"idfield":kfield,"valuefield":val,
          "objectIdField":orig,"source_fields":sorted(fields),"metadata_sha":metasha,
          "source_count_sha":csha,"first3_source_attributes_sha":sig,
          "first3_source_records":attrs,"original_srs":meta.get("extent",{}).get("spatialReference"),
          "full_source_not_downloaded_by_design":True
        }
        print("OFFICIAL_NATIVE_SOURCE_METADATA",key,"count",n,"key",kfield,"sample",attrs[:1],flush=True)
    rasters={}
    for lid in [21,22,23,24]:
        meta,sh=get(f"{BASE}/{lid}",{})
        rasters[str(lid)]={"type":meta.get("type"),"name":meta.get("name"),"metadata_sha":sh,
                           "no_guaranteed_numeric_download":True}
    out={
      "schema":"gesui_JRC_soil_EMEP_2015_DAYCENT_2010to19_native_basin_raster_METADATA_PREOUTCOME_v1",
      "utc":datetime.now(timezone.utc).isoformat(),
      "layers":result,
      "original_raster_layers":rasters,
      "native_SOIL_HydroID_and_WATER_BasinID_must_not_be_numeric_joined":True,
      "previous_water_land_NUTS2_2x2_SOURCE_STOP_preserved":True,
      "no_species_EQR_NUTRI_HNS_values_read":True
    }
    p=args.out/"JRC_HYDROID_74K_NATIVE_FEATURES_METADATA_ONLY.json"
    p.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print("JRC_NATIVE_SOIL_META_PROBE_COMPLETE",json.dumps({k:{a:result[k][a] for a in ("count","idfield","valuefield")} for k in result}),flush=True)
if __name__=="__main__":main()
