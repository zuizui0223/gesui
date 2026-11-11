#!/usr/bin/env python3
"""Official JRC GREEN 1990-2021 European point-source basin-load GDB SOURCE.

Verify original JRC item, 47MB zip, extract safely into CI scratch, read GDB
schema only: layers and source grain/year/country fields. Do NOT load EEA
biological EQR numeric values or treat GREEN modeled outlet nutrients as
observed N/P from an upstream sewage plant at an EQR monitoring site.
"""
from __future__ import annotations
import argparse,hashlib,json,zipfile,sys,os
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter
import pyogrio
import numpy as np
import pandas as pd

SOURCE_URL="https://arcgis-maps.jrc.ec.europa.eu/portal/sharing/rest/content/items/8b105c72109243798913853ad9d05072/data"
SOURCE_BYTES=47068808
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for x in iter(lambda:f.read(2*1024*1024),b""):h.update(x)
    return h.hexdigest()
def result(archive,tmp):
    if archive.stat().st_size!=SOURCE_BYTES or not zipfile.is_zipfile(archive):
        raise ValueError("Original official JRC 1990–2021 basin load archive missing or size changed")
    tmp.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        names=[i.filename for i in z.infolist()]
        if len(names)>4000:raise ValueError("JRC GDB unexpected directory size")
        for name in names:
            if Path(name).is_absolute() or ".." in Path(name).parts:
                raise ValueError("Original public GDB zip path is unsafe")
        gdbs=sorted({name.split("/",1)[0] for name in names if ".gdb/" in name.lower() and "/" in name})
        if len(gdbs)!=1:raise ValueError("Original JRC archive needs one .gdb folder")
        total=sum(i.file_size for i in z.infolist())
        if total>1500_000_000:raise ValueError("Original model geodatabase above source-only budget")
        z.extractall(tmp)
    gdb=tmp/gdbs[0]
    layers=pyogrio.list_layers(gdb)
    layerout=[]
    for name,geom in layers:
        inf=pyogrio.read_info(gdb,layer=name)
        fields=[str(x) for x in inf.get("fields",[])]
        candidates=[field for field in fields if any(q in field.lower() for q in (
            "year","time","1990","2016","2021","waste","psource","point","domestic","industry",
            "industrial","discharge","tn","tp","nitrogen","phosphorus","load","basin",
            "catch","river","src","n_point","p_point"))]
        layerout.append({
            "layer":str(name),"geometry_type":str(geom),
            "original_feature_count":int(inf.get("features") or 0),
            "original_GDB_field_count":len(fields),
            "original_GDB_fields":fields,
            "2016_and_point_source_model_fields":candidates,
            "source_geometry_CRS":str(inf.get("crs"))[:220],
            "layer_total_indexed_official_source":True
        })
    if not layerout:raise ValueError("Original JRC model GDB has no readable layers")
    return {
       "schema":"gesui_original_JRC_GREEN_1990to2021_basin_pointsource_2016_GDB_schema_SOURCE_ONLY_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "official_original_2026_JRC_BasinModel_GDB_url":SOURCE_URL,
       "official_archive_SHA256":sha(archive),
       "official_archive_bytes":archive.stat().st_size,
       "ZIP_GDB_top_folder":gdbs[0],
       "ZIP_file_count":len(names),
       "ZIP_uncompressed_bytes":total,
       "official_OPEN_public_ArcGIS_GDB_layers":layerout,
       "2016_modelled_point_source_is_NOT_measured_2016_POTW_annual_NP_outfall_mass":True,
       "JRC_GREEN_NP_values_are_at_river_basin_SEA_OUTLETS_not_EEA_upstream_site":True,
       "2022_2018_positive_BQE_EEA_observation_values_NOT_read":True,
       "2016_NP_BQE_model_P_value_not_fit":True,
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--sourcezip",type=Path,required=True)
    p.add_argument("--tmp",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    out=result(args.sourcezip,args.tmp)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"sha":out["official_archive_SHA256"],
        "source_GDB_layers":[{"layer":z["layer"],"feature_count":z["original_feature_count"],
        "fields":z["2016_and_point_source_model_fields"],"total_fields":z["original_GDB_field_count"]}
        for z in out["official_OPEN_public_ArcGIS_GDB_layers"]]},indent=2))
if __name__=="__main__":
    main()
