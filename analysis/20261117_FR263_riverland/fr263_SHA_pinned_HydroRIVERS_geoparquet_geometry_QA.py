#!/usr/bin/env python3
"""Only verifies genuine FR263 EEA stations against original HydroRIVERS EU v1 linework.

No TP/EQR/vegetation numerical outcomes are read. No rubber-band snapping:
every station must retain its original HYRIV_ID and the true official geometry.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
import pathlib
import statistics
import tempfile
import zipfile
from collections import Counter
import geopandas as gpd
import pandas as pd
from pyproj import Transformer
from shapely.geometry import Point
from shapely.ops import substring

HYDRO_SOURCE_SHA = "5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"

SOURCE_JSON_SHA = "2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264"
N_SITES = 263
N_REACH = 262
N_H7 = 73
N_H6 = 24

def file_sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        while b:=f.read(1<<20):h.update(b)
    return h.hexdigest()

def read_original(geopath, sourcepath):
    original_bytes=sourcepath.read_bytes()
    if hashlib.sha256(original_bytes).hexdigest()!=SOURCE_JSON_SHA:
        raise RuntimeError("Original locked date-only 263-station list SHA changed")
    dat=json.loads(original_bytes)
    ss=dat["original_source_only_exact_same_EEA2018_2022_biological_site_and_annual_measured_TP_support"]
    if len(ss)!=N_SITES:raise RuntimeError("Original source station n differs")
    if not all(s["source_2018_eligible"] and s["source_2022_eligible"] for s in ss):
        raise RuntimeError("Original 2018/2022 date source eligibility is not complete")
    s=pd.DataFrame([{
        "site_key":x["EEA2022_site_key"],"H7":str(x["HYBAS_ID_L07"]),
        "H6":str(x["HYBAS_ID_L06_parent"]), "WQ":x["original_WQ_ID"]}
        for x in ss])
    p=pd.read_csv(geopath,dtype={"site_key":str,"HYRIV_ID":str})
    if len(p)!=N_SITES or p.site_key.nunique()!=N_SITES or p.HYRIV_ID.nunique()!=N_REACH:
        raise RuntimeError("Original 263 coordinates and 262 reach IDs changed")
    if p["source_EEA2022_Git_blob_sha1"].nunique()!=1 or \
       p["source_EEA2022_Git_blob_sha1"].iloc[0]!="325f199e985d9862ffee5b0edb205c61a7be3833":
        raise RuntimeError("Source EEA2022 spatial provenance changed")
    d=s.merge(p,on="site_key",validate="1:1")
    if len(d)!=N_SITES or d.H7.nunique()!=N_H7 or d.H6.nunique()!=N_H6:
        raise RuntimeError("Source station ID join or parent groups changed")
    if not d[["lat","lon"]].notna().all().all():
        raise RuntimeError("Nonfinite source coordinates")
    return d

def read_original_hydrorivers(hydrozip,wanted):
    # This is a strict SHA-locked third-party GeoParquet conversion of HydroRIVERS
    # v1 Europe, NOT the official ZIP. Previous official-source curl HTTP403
    # remains a standalone failed source attempt; never silently recast.
    if file_sha(hydrozip)!=HYDRO_SOURCE_SHA:
        raise RuntimeError("Pinned HydroRIVERS v1 thirdparty GeoParquet SHA256 mismatch")
    g=gpd.read_parquet(hydrozip)
    if "HYRIV_ID" not in g.columns or g.crs is None:
        raise RuntimeError("Source converted GeoParquet missing original ID or CRS")
    g["HYRIV_ID"]=g.HYRIV_ID.astype(str)
    if g.HYRIV_ID.duplicated().any():
        raise RuntimeError("Duplicate source native HYRIV ID")
    chosen=g[g.HYRIV_ID.isin(wanted)].copy()
    if len(chosen)!=len(wanted):
        raise RuntimeError(f"True HYRIV ID source mismatch: found {len(chosen)}/{len(wanted)}")
    if chosen.geometry.isna().any() or chosen.geometry.is_empty.any():
        raise RuntimeError("Original river geometry missing")
    from shapely.geometry import LineString
    from shapely.ops import linemerge
    lines={}
    for r in chosen.itertuples(index=False):
        geom=getattr(r,g.geometry.name)
        if geom.geom_type=="MultiLineString":
            geom=linemerge(geom)
        if geom.geom_type!="LineString" or not geom.is_valid:
            raise RuntimeError("One original river reach cannot be resolved into one valid connected LineString")
        lines[str(r.HYRIV_ID)]=geom
    crs=chosen.crs
    metr=gpd.GeoSeries(lines.values(),index=list(lines),crs=crs).to_crs("EPSG:3035")
    return dict(zip(metr.index,metr.values)),len(g)


def audit(args):
    args.out.mkdir(parents=True,exist_ok=True)
    d=read_original(args.geo,args.datesource)
    lines,n_geo=read_original_hydrorivers(args.hydrorivers,set(d.HYRIV_ID))
    xy=Transformer.from_crs("EPSG:4326","EPSG:3035",always_xy=True)
    rows=[]
    for r in d.itertuples(index=False):
        pt=Point(*xy.transform(float(r.lon),float(r.lat)))
        line=lines[str(r.HYRIV_ID)]
        offset=pt.distance(line)
        foot=line.project(pt)
        seg=substring(line,max(0,foot-500),min(line.length,foot+500))
        seglen=seg.length if seg.geom_type=="LineString" else 0.
        snap_250=offset<=250
        segment_900=seglen>=900
        rows.append(dict(
          site_key=r.site_key,HYRIV_ID=str(r.HYRIV_ID),H7=str(r.H7),H6=str(r.H6),
          eaa_lat=float(r.lat),eea_lon=float(r.lon),
          source_reach_length_m=round(line.length,2),
          original_station_to_original_reach_m=round(offset,2),
          foot_to_original_reach_start_m=round(foot,2),
          foot_to_original_reach_end_m=round(line.length-foot,2),
          original_along_reach_window_length_m=round(seglen,2),
          same_original_reach_within_250m=bool(snap_250),
          segment_at_least_900m=bool(segment_900),
          geometry_pass=bool(snap_250 and segment_900)))
    s=pd.DataFrame(rows).sort_values("site_key")
    ok=s[s.geometry_pass]
    r={
      "schema":"GESUI_20261117_FR263_OFFICIAL_HYDRORIVERS_TRUE_REACH_GEOMETRY_GATE_V1",
      "status":"SOURCE_THIRDPARTY_GEOPARQUET_GEOMETRY_SUCCESS_NO_LAND_PIXELS_OR_ECOLOGY_FIT",
      "source_original_station_n":N_SITES,
      "source_original_H7":N_H7,
      "source_original_H6":N_H6,
      "original_required_reach_ids":N_REACH,
      "original_HydroRIVERS_bbox_features_scanned":n_geo,
      "original_reach_ids_found":len(lines),\n      "thirdparty_conversion_NOT_official_original_zip_byte_verification":True,
      "n_station_with_reach_distance_le_250m":int(s.same_original_reach_within_250m.sum()),
      "n_station_with_reach_window_ge_900m":int(s.segment_at_least_900m.sum()),
      "n_station_joint_geometry_pass":len(ok),
      "joint_geometry_pass_unique_H7":ok.H7.nunique(),
      "joint_geometry_pass_unique_H6":ok.H6.nunique(),
      "geometry_pass_suffices_new_predeclared_minima":bool(len(ok)>=120 and ok.H7.nunique()>=40 and ok.H6.nunique()>=15),
      "station_to_original_reach_offset_m":{
          "min":float(s.original_station_to_original_reach_m.min()),
          "median":float(s.original_station_to_original_reach_m.median()),
          "p90":float(s.original_station_to_original_reach_m.quantile(.9)),
          "max":float(s.original_station_to_original_reach_m.max())},
      "source_sha256":{
         "original_source_station_geo_csv":file_sha(args.geo),
         "original_date_only_source_json":file_sha(args.datesource),
         "SHA_pinned_HydroRIVERS_v1_thirdparty_GeoParquet":file_sha(args.hydrorivers)},
      "ecology_not_read":True,
      "no_HRVPP_pixels_or_riparian_response_measured":True,
      "next_source_stage":"Actual 10m 2018/2022 HR-VPP TPROD+QFLAG pixels, original upstream/matched cover strata, or an independently frozen genuine public optical alternative"}
    s.to_csv(args.out/"FR263_ORIGINAL_HYRIV_SOURCE_GEOMETRY_QA.csv",index=False)
    (args.out/"FR263_ORIGINAL_HYRIV_SOURCE_GEOMETRY_REPORT.json").write_text(
        json.dumps(r,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(r,ensure_ascii=False,indent=2),flush=True)
    return r

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geo",type=pathlib.Path,required=True)
    p.add_argument("--datesource",type=pathlib.Path,required=True)
    p.add_argument("--hydrorivers",type=pathlib.Path,required=True)
    p.add_argument("--out",type=pathlib.Path,required=True)
    audit(p.parse_args())

if __name__=="__main__":main()
