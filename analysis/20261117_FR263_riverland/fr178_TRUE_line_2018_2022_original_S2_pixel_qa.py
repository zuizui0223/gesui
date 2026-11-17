#!/usr/bin/env python3
"""FR178 source QA: genuine 2018/2022 10m S2 NDVI at true riverline distance.

B1=20-100m, B2=100-250m, B3=250-500m, B4=500-1000m from exact
HYRIV_ID polyline (station +/-500m along the true projected riverline).
Original source 2021 WorldCover natural land classes are applied equally in
2018 and 2022. TP, benthic EQR and wastewater values are never loaded.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
from contextlib import ExitStack
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.windows import from_bounds,Window
from rasterio.warp import reproject
from rasterio.enums import Resampling
from pyproj import Transformer
from shapely import points as shp_points, distance as shp_distance
from shapely.geometry import Point
from shapely.ops import transform as shp_transform, substring, linemerge
from fr263_2018_2022_S2_COG_full24H6_pixel_QA import asset,worldcover,summarize_pair

SCHEMA="GESUI_FR178_2018_2022_ORIGINAL_10M_TRUE_RIVERLINE_S2_PIXEL_SHARD_V1"
IN_SCHEMA="GESUI_FR178_2018_2022_S2_SCENE_METADATA_SOURCE_SHARD_V1"
GEOSHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
INPUT_SHA="9a5d11a6c50ba49160da46097873c3d5eedda881e0a173f1aa90968c94950b11"
SOURCE_PRESPEC_N=178
NATURAL=(10,20,30,90)
DIST_BANDS=((20,100),(100,250),(250,500),(500,1000))
BAND_MIN_NATURAL=100
COVER_CLASS_MIN=10
MONTHS=(6,7,8)
MAX_OFF_M=250
MIN_ALONG_M=900

def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        while b:=f.read(1<<20):h.update(b)
    return h.hexdigest()

def source(sourcefile,riverpath):
    j=json.loads(sourcefile.read_text())
    if j.get("schema")!=IN_SCHEMA or j.get("source_FR178_frozen_csv_SHA256")!=INPUT_SHA or j.get("original_178_site_population")!=178:
        raise RuntimeError("Original fixed STAC shard metadata/signature changed")
    rows=j["source_sites"]
    if len(rows)!=j["source_shard_count_sites"] or not rows:
        raise RuntimeError("Original source shard not exhaustively enumerated")
    if sha(riverpath)!=GEOSHA:raise RuntimeError("Original source v1 HydroRIVERS conversion SHA changed")
    gd=gpd.read_parquet(riverpath)
    if gd.crs is None or "HYRIV_ID" not in gd:raise RuntimeError("River source CRS/ID absent")
    gd["HYRIV_ID"]=gd["HYRIV_ID"].astype(str)
    ids=set(x["HYRIV_ID"] for x in rows)
    if gd.HYRIV_ID.duplicated().any():raise RuntimeError("Original river IDs duplicated")
    selected=gd[gd.HYRIV_ID.isin(ids)].copy()
    if len(selected)!=len(ids):raise RuntimeError("Original exact water station HYRIV_ID not in true river")
    geo={}
    to=Transformer.from_crs("EPSG:4326","EPSG:3035",always_xy=True)
    for v in selected.to_crs("EPSG:3035").itertuples(index=False):
        geom=getattr(v,gd.geometry.name)
        if geom.geom_type=="MultiLineString":geom=linemerge(geom)
        if geom.geom_type!="LineString" or not geom.is_valid:raise RuntimeError("Real river geometry invalid")
        geo[str(v.HYRIV_ID)]=geom
    segments={}
    for x in rows:
        pt=Point(*to.transform(x["original_lon"],x["original_lat"]))
        river=geo[x["HYRIV_ID"]]
        off=pt.distance(river)
        pos=river.project(pt)
        segment=substring(river,max(0,pos-500),min(river.length,pos+500))
        if off>MAX_OFF_M or segment.geom_type!="LineString" or segment.length<MIN_ALONG_M:
            raise RuntimeError("Source-frozen true reach eligibility changed before optical values")
        segments[x["site_key"]]=segment
    return j,segments

def one_scene(scene,site,segment):
    vals={k:asset(scene,k) for k in ("red","nir","scl")}
    if abs(vals["red"][1]-vals["nir"][1])>1e-12 or abs(vals["red"][2]-vals["nir"][2])>1e-12:
        raise RuntimeError("Original real red/nir scaling mismatch")
    lat=float(site["original_lat"]);lon=float(site["original_lon"])
    WC,tile=worldcover(lat,lon)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_HTTP_MAX_RETRY="2",GDAL_HTTP_RETRY_DELAY="1",
      GDAL_CACHEMAX=64,VSI_CACHE="TRUE",VSI_CACHE_SIZE=1000000):
      with ExitStack() as st:
        src={q:st.enter_context(rasterio.open("/vsicurl/"+vals[q][0]))
             for q in ("red","nir","scl")}
        world=st.enter_context(rasterio.open("/vsicurl/"+WC))
        red=src["red"];nir=src["nir"];scl=src["scl"]
        if red.crs is None or nir.crs!=red.crs or scl.crs!=red.crs or world.crs is None or red.transform!=nir.transform:
            raise RuntimeError("Source original red NIR SCL/land raster projection inconsistent")
        if not (8<=abs(red.transform.a)<=12 and 8<=abs(red.transform.e)<=12):
            raise RuntimeError("Original Sentinel2 red native resolution not 10m")
        to_raster=Transformer.from_crs("EPSG:3035",red.crs,always_xy=True)
        segment_raster=shp_transform(to_raster.transform,segment)
        bounds=segment_raster.buffer(1020).bounds
        w=from_bounds(*bounds,transform=red.transform)
        w=Window(int(math.floor(w.col_off)),int(math.floor(w.row_off)),
                 int(math.ceil(w.width))+2,int(math.ceil(w.height))+2)
        rraw=red.read(1,window=w,boundless=True,fill_value=0)
        nraw=nir.read(1,window=w,boundless=True,fill_value=0)
        t=red.window_transform(w)
        sc=np.zeros(rraw.shape,dtype=np.uint8)
        wc=np.zeros(rraw.shape,dtype=np.uint8)
        reproject(source=rasterio.band(scl,1),destination=sc,
          dst_transform=t,dst_crs=red.crs,resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
        reproject(source=rasterio.band(world,1),destination=wc,
          dst_transform=t,dst_crs=red.crs,resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
        rval=rraw.astype(np.float32)*vals["red"][1]+vals["red"][2]
        nval=nraw.astype(np.float32)*vals["nir"][1]+vals["nir"][2]
        good=(rval>.005)&(nval>.005)&(rval<1)&(nval<1)&np.isin(sc,[4,5])&np.isin(wc,NATURAL)
        raw_eligible=int(good.sum())
        if raw_eligible<2*BAND_MIN_NATURAL:
            return {"pixel_source_pass":False,"reason":"too_few_clear_valid_natural_source_pixels","source_valid_n":raw_eligible}
        yi,xi=np.nonzero(good)
        xx=t.c+(xi+.5)*t.a+(yi+.5)*t.b
        yy=t.f+(xi+.5)*t.d+(yi+.5)*t.e
        pp=shp_points(xx,yy)
        d=shp_distance(pp,segment_raster)
        bands={}
        for lower,upper in DIST_BANDS:
            key=f"{lower}_{upper}m"
            idx=(d>lower)&(d<=upper)
            bands[key]={"n_valid_source_pixels":int(idx.sum()),"class_pixel_counts":{
               str(k):int(((wc[yi,xi]==k)&idx).sum()) for k in NATURAL}}
        # Require original prespecified source support in ALL FOUR river-distance
        # bands, not only the near/far contrast. This gate uses pixel support,
        # never TP/EQR or directional NDVI effects.
        if any(v["n_valid_source_pixels"]<BAND_MIN_NATURAL for v in bands.values()):
            return {"pixel_source_pass":False,"reason":"original_any_of_four_river_distance_bands_under_100_valid_pixels",
                    "distance_band_source_counts":bands}
        near=(d>20)&(d<=100)
        far=(d>500)&(d<=1000)
        if min(int(near.sum()),int(far.sum()))<BAND_MIN_NATURAL:
            return {"pixel_source_pass":False,"reason":"insufficient_near_or_far_actual_riverline_pixels",
                    "distance_band_source_counts":bands}
        classes=sorted(set(wc[yi,xi][near])&set(wc[yi,xi][far]))
        NDVI=(nval-rval)/(nval+rval+1e-12)
        nd=NDVI[yi,xi]
        matched=[]
        for k in classes:
            a=near&(wc[yi,xi]==k);b=far&(wc[yi,xi]==k)
            if int(a.sum())<COVER_CLASS_MIN or int(b.sum())<COVER_CLASS_MIN:continue
            matched.append({"class":int(k),"n_near":int(a.sum()),"n_far":int(b.sum()),
                "median_near":float(np.median(nd[a])),"median_far":float(np.median(nd[b]))})
        if not matched:
            return {"pixel_source_pass":False,"reason":"no_both_distance_bands_common_source_cover_class",
                    "distance_band_source_counts":bands}
        return {"pixel_source_pass":True,"reason":"original_real_pixel_qa_pass",
           "distance_band_source_counts":bands,
           "cover_matched_original_real_pixel_NDVI":matched,
           "source_scene_scale":vals["red"][1],"source_scene_offset":vals["red"][2],
           "static_WorldCover2021_tile":tile,
           "original_actual_hydrorivers_true_segment_2018_22_same":True}

def run(a):
    source_data,segments=source(a.stac,a.rivers)
    outcome=[];err=Counter()
    for idx,x in enumerate(source_data["source_sites"],1):
        passed_months=[];attempts=[]
        if x.get("STAC_scene_metadata_two_years_gate"):
          for m in x["original_JJA_same_calendar_months"]:
            year={}
            for yr in (2018,2022):
                scenes=sorted((sc for sc in x["years"][str(yr)]["source_selected_scenes"] if sc["month"]==m),
                      key=lambda y:(y["cloud_pct"],y["datetime"],y["scene_id"]))
                if len(scenes)>2:raise RuntimeError("Original source STAC >two scenes/month")
                for scene in scenes:
                    try:res=one_scene(scene,x,segments[x["site_key"]])
                    except Exception as e:
                        res={"pixel_source_pass":False,"reason":"original_COG_access_or_raster_error",
                             "error_type":type(e).__name__,"error":str(e)[:250]}
                    attempts.append({"year":yr,"month":m,"scene":scene["scene_id"],"reason":res.get("reason"),
                                     "pixel_source_pass":res.get("pixel_source_pass")})
                    if res["pixel_source_pass"]:
                        year[yr]={"month":m,**res};break
                    err[res.get("reason","unknown")]+=1
            if 2018 in year and 2022 in year:
                p=summarize_pair(year[2018],year[2022])
                if p is not None:passed_months.append(p)
        qualify=len(passed_months)>=2
        outcome.append({"site_key":x["site_key"],"H7":x["H7"],"H6":x["H6"],
            "HYRIV_ID":x["HYRIV_ID"],"origin_STAC_source_bothyear":x["STAC_scene_metadata_two_years_gate"],
            "original_TP_JJA_same_months":x["original_JJA_same_calendar_months"],
            "source_two_years_TP_compatible_months_with_real_NDVI_n":len(passed_months),
            "source_two_years_real_satellite_qa_PASS":qualify,
            "source_original_near100_far1000_delta_NDVI_by_month":passed_months,
            "source_real_COG_attempts":attempts})
        print("GESUI_FR178_TRUE_RIVER_PIXEL_SHARD="+str(source_data["source_shard_index"])+
              " SITE="+str(idx)+"/"+str(len(source_data["source_sites"]))+
              " qualified="+str(qualify)+" pairedmonths="+str(len(passed_months)),flush=True)
    eligible=[x for x in outcome if x["source_two_years_real_satellite_qa_PASS"]]
    j={"schema":SCHEMA,"status":"COMPLETE_SOURCE_QA_SHARD_NO_AQUATIC_ECOLOGY_MODEL",
      "source_original_shard":source_data["source_shard_index"],
      "original_source_S2_STAC_json_SHA256":sha(a.stac),
      "source_HydroRIVERS_converted_geo_SHA256":GEOSHA,
      "source_total_178_population":SOURCE_PRESPEC_N,"source_shard_sites":len(outcome),
      "same_calendar_month_two_years_real_true_riverline_optical_sites":len(eligible),
      "predeclared_original_true_river_distance_bands_m":DIST_BANDS,
      "source_original_segment_along_channel_m":"station projected to same HYRIV +/-500m",
      "near_clear_natural_min100_far_clear_natural_min100":True,
      "source_2021_WorldCover_original_land_class_mask":NATURAL,
      "source_failed_scene_counts":dict(err),
      "original_EEA_river_TP_EQR_butterflies_sewage_values_not_loaded":True,
      "source_remote_sensed_NDVI_not_actual_biomass_or_nitrogen_flux":True,
      "records":outcome}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/f"FR178_TRUE_RIVERLINE_S2_PIXEL_SHARD{source_data['source_shard_index']}.json").write_text(
        json.dumps(j,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_FR178_PIXEL_COMPLETE="+json.dumps({
        "shard":source_data["source_shard_index"],"n":len(outcome),"source_real_pixel_passing_n":len(eligible)}),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stac",type=Path,required=True)
    p.add_argument("--rivers",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
if __name__=="__main__":main()
