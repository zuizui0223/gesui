#!/usr/bin/env python3
"""Genuine official Europe HydroRIVERS shapefile (v1) bounded source audit.

Read only published official HydroRIVERS v1 Europe and Middle East shapefile
ZIP, verify archive and river reach source keys; do not compute biodiversity
outcomes or map any individual's GBIF position. License forbids unapproved
redistribution of original HydroRIVERS data; upload only source metadata.
"""
from __future__ import annotations
import argparse,hashlib,io,json,zipfile
from pathlib import Path
from urllib.parse import urlparse
import requests,shapefile

URL="https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip"
MAX_BYTES=90*1024*1024
UA="gesui-20261116-HydroRIVERS-official-Europe-source-metadata-only/1"

def download(path):
    h=hashlib.sha256();size=0
    with requests.get(URL,headers={"User-Agent":UA,"Accept":"application/zip"},stream=True,timeout=(15,90)) as r:
        if r.status_code!=200:raise RuntimeError("Official HydroRIVERS source HTTP "+str(r.status_code))
        if urlparse(r.url).hostname!="data.hydrosheds.org":raise RuntimeError("Redirected off hydrosheds source")
        with path.open("wb") as f:
            for block in r.iter_content(1048576):
                if not block:continue
                size+=len(block)
                if size>MAX_BYTES:raise RuntimeError("Unexpected huge source zip")
                h.update(block);f.write(block)
    if not zipfile.is_zipfile(path):raise RuntimeError("HydroRIVERS original archive is not genuine ZIP")
    return {"source_url":URL,"source_file_SHA256":h.hexdigest(),"original_bytes":size}

def inspect(p):
    with zipfile.ZipFile(p) as z:
        entries=z.infolist()
        by_base={}
        for ent in entries:
            if ent.is_dir():continue
            n=ent.filename
            ext=n.rsplit(".",1)[-1].lower()
            if ext in ("shp","shx","dbf","prj"):
                b=n.rsplit(".",1)[0]
                by_base.setdefault(b,{})[ext]=n
        group=[g for g in by_base.values() if {"shp","shx","dbf"}.issubset(g)]
        if len(group)!=1:raise RuntimeError("Unexpected official HydroRIVERS shapefile count "+str(len(group)))
        group=group[0]
        shp=shapefile.Reader(shp=io.BytesIO(z.read(group["shp"])),
                    shx=io.BytesIO(z.read(group["shx"])),
                    dbf=io.BytesIO(z.read(group["dbf"])))
        names=[x[0] for x in shp.fields[1:]]
        original_count=len(shp)
        if "HYRIV_ID" not in names:raise RuntimeError("Original river reach ID not found")
        if original_count<1000:raise RuntimeError("Original EU HydroRIVERS number of reaches unreasonably low")
        types={str(shp.shape(i).shapeTypeName) for i in range(min(100,original_count))}
        inspect_names=[x for x in ["HYRIV_ID","NEXT_DOWN","HYBAS_L12","LENGTH_KM","DIS_AV_CMS",
                "ORD_STRA","ORD_FLOW","DIST_DN_KM","DIST_UP_KM"] if x in names]
        sample=[]
        for i in range(min(original_count,100)):
            row=shp.record(i)
            sample.append({k:row[names.index(k)] for k in inspect_names if k!="HYRIV_ID"})
        prj=z.read(group["prj"]).decode("utf-8",errors="replace") if "prj" in group else ""
        if not prj:raise RuntimeError("Original projection file absent")
        out={"original_river_reach_count":original_count,
             "native_HYRIV_ID_exists":True,
             "native_downstream_edge_exists":"NEXT_DOWN" in names,
             "native_HydroBASINS_lower_level_overlap_field_exists":"HYBAS_L12" in names,
             "field_names":names,
             "first100_geometry_type_names":sorted(types),
             "source_projection_header_excerpt":prj[:240],
             "original_zip_members":[{"filename":x.filename,"original_bytes":x.file_size,
                   "compression":x.compress_type} for x in entries if not x.is_dir()][:20],
             "first100_downstream_no_coordinates_or_ID_values_exported":True,
             "no_original_hydrography_geometry_data_rehosted":True}
        return out

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    f=a.output/"TEMP_HydroRIVERS_v10_eu_shp.zip"
    src=download(f)
    inspect_out=inspect(f)
    f.unlink()
    out={"schema":"GESUI_20261116_FIG4_OFFICIAL_HYDRORIVERS_EU_ROAD_DISTANCE_SOURCE_QA_V1",
         "official_HydroRIVERS_v1_EU_zip":src,
         "hydro_source_geometry_schema":inspect_out,
         "published_HydroRIVERS_v1_threshold_min_upstream_km2":10,
         "published_HydroRIVERS_v1_threshold_min_natural_discharge_m3_s":0.1,
         "small_tributaries_below_threshold_may_be_missing":True,
         "not_a_H7_river_proximity_observation_or_biological_EIV_score":True,
         "source_only_no_EQR_HNS_plantN_or_W_outcome_opened":True}
    (a.output/"HYDRORIVERS_EU_FIG4_ORIGINAL_RIVER_NETWORK_SOURCE_QA.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_FIG4_HYDRORIVERS_EU_SOURCE="+json.dumps({
        "source":src,"n_reaches":inspect_out["original_river_reach_count"],
        "fields":inspect_out["field_names"],
        "shape_types_first100":inspect_out["first100_geometry_type_names"],
        "native_downstream":inspect_out["native_downstream_edge_exists"],
        "original_lower_HYBAS_link":inspect_out["native_HydroBASINS_lower_level_overlap_field_exists"],
        "source_ready_no_distance_band_model_run":True},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path);main(p.parse_args())
