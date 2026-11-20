#!/usr/bin/env python3
"""HydroRIVERS v1 Europe+Middle East original network source-only availability.

Download official HydroSHEDS ZIP (NOT republished in Git), verify archive
shape, original .shp/.shx/.dbf/.prj sidecars, the native NEXT_DOWN directed
river reach identifier, and original coordinate reference. No station
snapping, WWTP/BMWP values, river-routed exposure, effect models or maps
are created in this phase.

HydroRIVERS v1 includes only reaches with area >=10km2 or flow >=0.1 cms;
omitted headwaters, canals and drainage ditches cannot be treated as no
downstream water connectivity.
"""
from __future__ import annotations

import argparse,csv,hashlib,io,json,urllib.request,zipfile,tempfile
from pathlib import Path

URL="https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip"
MAX_BYTES=130*1024*1024
UA="gesui-official-HydroRIVERS-v1-Europe-net-source-gate-20261120"
REQUIRED={"HYRIV_ID","NEXT_DOWN","LENGTH_KM"}
def main(a):
    out={"schema":"GESUI_20261120_HYDRORIVERS_EU_V1_ORIGINAL_SOURCE_GATE_V1",
         "source_only":True,"provider":"HydroSHEDS HydroRIVERS v1 Europe+Middle East",
         "official_url":URL,"status":"SOURCE_STOP",
         "does_not_measure_2017_18_source_discharge_or_topological_chem_delivery":True}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_hyriv_v1_") as root:
            f=Path(root)/"original.zip"
            sha=hashlib.sha256();count=0
            req=urllib.request.Request(URL,headers={"User-Agent":UA})
            with urllib.request.urlopen(req,timeout=120) as resp,f.open("wb") as dest:
                if resp.status!=200:raise RuntimeError("Original HydroRIVERS source HTTP != 200")
                while True:
                    block=resp.read(2**20)
                    if not block:break
                    count+=len(block)
                    if count>MAX_BYTES:raise RuntimeError("HydroRIVERS EU source above safety limit")
                    sha.update(block);dest.write(block)
            if count<10*1024*1024 or not zipfile.is_zipfile(f):
                raise RuntimeError("Unexpected small/nonZIP official Europe source")
            with zipfile.ZipFile(f) as z:
                files=[x for x in z.infolist() if not x.is_dir()]
                extensions={Path(i.filename).suffix.lower() for i in files}
                if not {".shp",".shx",".dbf",".prj"}.issubset(extensions):
                    raise RuntimeError("Missing native source shapefile sidecars")
                names=[i.filename for i in files]
                shp=[x for x in files if x.filename.lower().endswith(".shp")]
                if len(shp)!=1:raise RuntimeError("Unexpected multiple HydroRIVERS river layer shapefiles")
                # PyShp can extract source attribute definition without loading all lines.
                import shapefile
                p=Path(root)/"official"
                p.mkdir()
                for i in files:
                    if i.filename.lower().endswith((".dbf",".shx",".shp",".prj",".cpg")):
                        dest=p/Path(i.filename).name
                        with z.open(i) as src,dest.open("wb") as dst:
                            import shutil
                            shutil.copyfileobj(src,dst)
                source_file=p/Path(shp[0].filename).name
                with shapefile.Reader(str(source_file)) as sf:
                    fields=[f[0] for f in sf.fields if f[0]!="DeletionFlag"]
                    first=next(sf.iterRecords(),None)
                    n_records=len(sf)
                    bbox=list(sf.bbox)
                if not REQUIRED.issubset(set(fields)):
                    raise RuntimeError("HydroRIVERS source lacks NEXT_DOWN/HYRIV_ID/length")
                # No original reach geometries or records republished.
                out.update({
                    "status":"OFFICIAL_HYDRORIVERS_V1_EU_SOURCE_SCHEMA_AND_TOPOLOGY_FIELDS_VERIFIED",
                    "source_zip_SHA256":sha.hexdigest(),
                    "source_zip_bytes":count,
                    "source_n_members":len(files),
                    "shp_member_SHA256":hashlib.sha256(source_file.read_bytes()).hexdigest(),
                    "original_source_layer_name":Path(shp[0].filename).name,
                    "native_attribute_field_names":fields,
                    "native_reach_record_count":n_records,
                    "original_degree_bbox":bbox,
                    "original_source_projection_member":{
                        "name":next((x for x in names if x.endswith(".prj")),None),
                        "schema_member_present":".prj" in extensions
                    },
                    "source_HAS_directed_NEXT_DOWN":True,
                    "original_high_resolution_small_channels_missing_possible":True
                })
    except Exception as ex:
        out.update({"status":"SOURCE_STOP",
                    "source_error":type(ex).__name__+": "+str(ex)[:550]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("GESUI_OFFICIAL_EU_HYDRORIVERS_SOURCE_GATE="+json.dumps({
        "status":out["status"],
        "error":out.get("source_error"),
        "n_zip_bytes":out.get("source_zip_bytes"),
        "source_zip_SHA256":out.get("source_zip_SHA256"),
        "n_original_river_reaches":out.get("native_reach_record_count"),
        "native_fields":out.get("native_attribute_field_names"),
        "source_bbox":out.get("original_degree_bbox"),
        "HAS_NEXT_DOWN":out.get("source_HAS_directed_NEXT_DOWN"),
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    main(p.parse_args())
