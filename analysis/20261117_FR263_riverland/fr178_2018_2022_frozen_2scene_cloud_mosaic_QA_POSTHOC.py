#!/usr/bin/env python3
"""POSTHOC source-only technical sensitivity: combine TWO ALREADY SOURCE-FROZEN
S2 scenes/month to mitigate complementary clouds for the exact FR178 cohort.

This DOES NOT revise the original independent all-four-bands 116/178 STOP,
does NOT choose dates/classes/radii on TP/EQR or NDVI signs, and does NOT
establish wastewater-derived terrestrial nutrient flux. It creates a different,
explicitly exploratory paired-month 2018/22 optical measurement.
"""
from __future__ import annotations
import argparse
import collections
import json
import math
from pathlib import Path
from contextlib import ExitStack
import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.warp import reproject
from rasterio.enums import Resampling
from shapely import points as shp_points, distance as shp_distance
from fr178_TRUE_line_2018_2022_original_S2_pixel_qa import (
    source, sha, DIST_BANDS, BAND_MIN_NATURAL, NATURAL, COVER_CLASS_MIN
)
from fr263_2018_2022_S2_COG_full24H6_pixel_QA import asset,worldcover,summarize_pair

SCHEMA="GESUI_FR178_S2_FROZEN_2SCENE_SAME_MONTH_CLOUD_MOSAIC_SOURCE_QA_POSTHOC_V1"
N_EXPECTED=178
N_H7_REQUIRED=40
N_H6_REQUIRED=15
N_SITES_REQUIRED=120
MONTHS=(6,7,8)
COG_HOST="https://sentinel-cogs.s3.us-west-2.amazonaws.com/"
TARGET_CRS="EPSG:3035"
RESOLUTION_M=10
MAX_SOURCE_SCENES=2
GATE="original_file_SCL4_or5_radiometrically_valid_and_WC2021_natural_mask"

def source_scene(scene):
    vals={q:asset(scene,q) for q in ("red","nir","scl")}
    if abs(vals["red"][1]-vals["nir"][1])>1e-12 or abs(vals["red"][2]-vals["nir"][2])>1e-12:
        raise RuntimeError("Signed scene red/NIR physical reflectance scale+offset changed")
    return vals

def true_fixed_grid(segment):
    # A 10m grid anchored to EPSG3035 multiples of 10m ensures identical
    # target pixels across both years, scenes, and months for every station.
    xmin,ymin,xmax,ymax=segment.buffer(1010).bounds
    xmin=math.floor(xmin/10)*10
    ymax=math.ceil(ymax/10)*10
    width=math.ceil((xmax-xmin)/10)
    height=math.ceil((ymax-ymin)/10)
    if not (100<=width<=650 and 100<=height<=650):
        raise RuntimeError("Original metric bank grid is oversized/undersized")
    tr=from_origin(xmin,ymax,10,10)
    return tr,int(height),int(width)

def good_pixel_ndvi(scene,target_transform,shape):
    vals=source_scene(scene)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_HTTP_MAX_RETRY="2",
                      GDAL_HTTP_RETRY_DELAY="1",GDAL_CACHEMAX=64,VSI_CACHE="TRUE",
                      VSI_CACHE_SIZE=1000000):
      with ExitStack() as st:
        src={k:st.enter_context(rasterio.open("/vsicurl/"+vals[k][0])) for k in ("red","nir","scl")}
        red,nir,scl=src["red"],src["nir"],src["scl"]
        if red.crs is None or nir.crs!=red.crs or scl.crs!=red.crs:
            raise RuntimeError("Original source scene band CRS mismatch")
        if red.transform!=nir.transform or red.width!=nir.width or red.height!=nir.height:
            raise RuntimeError("Original S2 red/nir source grid mismatch")
        if not 8<=abs(red.transform.a)<=12:
            raise RuntimeError("Original red/NIR native resolution not 10m")
        rr=np.zeros(shape,dtype=np.float32)
        nn=np.zeros(shape,dtype=np.float32)
        sc=np.zeros(shape,dtype=np.uint8)
        for name,dst in (("red",rr),("nir",nn),("scl",sc)):
            reproject(
                source=rasterio.band(src[name],1),destination=dst,
                dst_transform=target_transform,dst_crs=TARGET_CRS,
                src_nodata=0,dst_nodata=0,resampling=Resampling.nearest)
        r=rr*vals["red"][1]+vals["red"][2]
        n=nn*vals["nir"][1]+vals["nir"][2]
        valid=(np.isin(sc,[4,5])&(r>0.005)&(n>0.005)&(r<1)&(n<1))
        ndvi=(n-r)/(n+r+1e-12)
        valid&=np.isfinite(ndvi)&(ndvi>=-1)&(ndvi<=1)
        return np.where(valid,ndvi,0).astype("float32"),valid

def get_original_wc(site,tr,shape):
    url,_=worldcover(float(site["original_lat"]),float(site["original_lon"]))
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
                      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_HTTP_MAX_RETRY="2",
                      GDAL_HTTP_RETRY_DELAY="1",GDAL_CACHEMAX=64):
      with rasterio.open("/vsicurl/"+url) as wc:
        if wc.crs is None:raise RuntimeError("Original WorldCover raster source lacks CRS")
        classes=np.zeros(shape,dtype=np.uint8)
        reproject(source=rasterio.band(wc,1),destination=classes,
                  dst_transform=tr,dst_crs=TARGET_CRS,
                  src_nodata=0,dst_nodata=0,resampling=Resampling.nearest)
    return classes

def band_masks(tr,shape,segment):
    yi,xi=np.indices(shape,dtype=np.int32)
    coordsx=tr.c+(xi+.5)*tr.a
    coordsy=tr.f+(yi+.5)*tr.e
    pp=shp_points(coordsx,coordsy)
    dd=shp_distance(pp,segment)
    return {f"{lo}_{hi}m":((dd>lo)&(dd<=hi)) for lo,hi in DIST_BANDS}

def source_fixed_month(site,segment,year,month,tr,shape,wc,bands):
    scenes=sorted(
        (s for s in site["years"][str(year)]["source_selected_scenes"] if s["month"]==month),
        key=lambda x:(x["cloud_pct"],x["datetime"],x["scene_id"])
    )
    if len(scenes)>MAX_SOURCE_SCENES:
        raise RuntimeError("Original source preselected more than two scenes per month")
    sums=np.zeros(shape,dtype=np.float32)
    counts=np.zeros(shape,dtype=np.uint8)
    attempts=[]
    for scene in scenes:
        try:
            nd,valid=good_pixel_ndvi(scene,tr,shape)
            sums+=nd
            counts+=valid.astype("uint8")
            attempts.append({"source_scene":scene["scene_id"],
               "source_scene_UTC":scene["datetime"],
               "scene_tile_cloud_pct":scene["cloud_pct"],
               "original_good_pixel_n":int(valid.sum()),"error":None})
        except (OSError,RuntimeError,ValueError,rasterio.errors.RasterioIOError) as ex:
            attempts.append({"source_scene":scene["scene_id"],
                 "source_scene_UTC":scene["datetime"],
                 "source_scene_UTC_source_year":year,
                 "error":type(ex).__name__+":"+str(ex)[:150]})
    natural=(counts>0)&np.isin(wc,NATURAL)
    pixcounts={k:int(np.sum(natural&mask)) for k,mask in bands.items()}
    common={"year":year,"month":month,"source_scene_n_attempted":len(scenes),
            "source_scene_attempts":attempts,
            "original_pixel_mosaic_method":"arithmetic mean valid NDVI across up to 2 frozen 2018/22 same-month scenes, equal 10m LAEA grid, never max NDVI",
            "band_counts":pixcounts}
    if any(pixcounts[k]<BAND_MIN_NATURAL for k in bands):
        return {**common,"pixel_source_pass":False,
                "reason":"at_least_one_of_original_four_bands_under_100_composite_natural_pixels"}
    near=bands["20_100m"]
    far=bands["500_1000m"]
    merged=np.zeros(shape,dtype=np.float32)
    np.divide(sums,counts,out=merged,where=counts>0)
    classmatched=[]
    for cl in NATURAL:
        i=(near&natural&(wc==cl))
        o=(far&natural&(wc==cl))
        ni=int(i.sum());no=int(o.sum())
        if ni<COVER_CLASS_MIN or no<COVER_CLASS_MIN:continue
        classmatched.append({
            "class":int(cl),"n_near":ni,"n_far":no,
            "median_near":float(np.median(merged[i])),
            "median_far":float(np.median(merged[o]))
        })
    if not classmatched:
        return {**common,"pixel_source_pass":False,
                "reason":"no_class_with_ten_pixels_near_and_far_in_same_composite"}
    return {**common,"pixel_source_pass":True,"reason":"source_2scene_composite_pass",
            "cover_matched_original_real_pixel_NDVI":classmatched}

def synthetic_smoke():
    # Numeric test: cloud holes in two scenes fill each other.
    a=np.array([[.1,0],[0,.5]],dtype="float32")
    b=np.array([[0,.2],[.4,0]],dtype="float32")
    va=np.array([[True,False],[False,True]])
    vb=np.array([[False,True],[True,False]])
    s=a+b;c=va.astype("uint8")+vb.astype("uint8")
    out=np.zeros_like(s);np.divide(s,c,out=out,where=c>0)
    assert np.allclose(out,np.array([[.1,.2],[.4,.5]]))
    assert int((c>0).sum())==4
    return True

def run(args):
    synthetic_smoke()
    j,segments=source(args.stac,args.rivers)
    results=[]
    failures=collections.Counter()
    for k,site in enumerate(j["source_sites"],1):
        tr,h,w=true_fixed_grid(segments[site["site_key"]])
        shape=(h,w)
        wc=get_original_wc(site,tr,shape)
        bands=band_masks(tr,shape,segments[site["site_key"]])
        allmonths=[];summary=[];scene_errors=[]
        if site["STAC_scene_metadata_two_years_gate"]:
            for month in site["original_JJA_same_calendar_months"]:
                d={}
                for year in (2018,2022):
                    d[year]=source_fixed_month(site,segments[site["site_key"]],year,month,tr,shape,wc,bands)
                    summary.append({**{q:d[year][q] for q in ("year","month","pixel_source_pass","reason",
                         "band_counts","source_scene_n_attempted")}})
                    if not d[year]["pixel_source_pass"]:
                        failures[(year,d[year]["reason"])]+=1
                if d[2018]["pixel_source_pass"] and d[2022]["pixel_source_pass"]:
                    matched=summarize_pair(d[2018],d[2022])
                    if matched is not None:allmonths.append(matched)
        passed=len(allmonths)>=2
        results.append({"site_key":site["site_key"],"H7":site["H7"],"H6":site["H6"],
                        "HYRIV_ID":site["HYRIV_ID"],
                        "source_original_TP_JJA_same_months":site["original_JJA_same_calendar_months"],
                        "source_two_years_real_mosaic_pixel_eligible":passed,
                        "source_2018_2022_real_mosaic_months_n":len(allmonths),
                        "posthoc_composite_delta_NDVI_by_month":allmonths,
                        "original_monthwise_source_QA":summary})
        print("GESUI_FR178_MOSAIC_SOURCE_SHARD="+str(j["source_shard_index"])+
              " SITE="+str(k)+"/"+str(len(j["source_sites"]))+
              " goodmonths="+str(len(allmonths))+" eligible="+str(passed),flush=True)
    eligible=[r for r in results if r["source_two_years_real_mosaic_pixel_eligible"]]
    report={"schema":SCHEMA,"status":"POSTHOC_COMPOSITE_SHARD_COMPLETED_NO_AQUATIC_MODEL",
      "original_source_shard":j["source_shard_index"],
      "original_source_178_pop":N_EXPECTED,"original_shard_site_count":len(results),
      "original_source_STAC_json_SHA256":sha(args.stac),
      "genuine_pinned_HydroRIVERS_converted_source_SHA256":"5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f",
      "posthoc_TP_compatible_two_year_2scene_mosaic_eligible_sites":len(eligible),
      "fixed_original_4_river_bands_m":[[int(a),int(b)] for a,b in DIST_BANDS],
      "fixed_min_100_natural_clear_10m_pixels_per_band_both_years":True,
      "original_max_two_STAC_scenes_per_month":True,
      "original_static_2021_WorldCover_landclasses":list(NATURAL),
      "posthoc_source_improvement_does_not_revive_original_116_below120_STOP":True,
      "the_study_reuses_original_178_cohort_and_is_NOT_independent_replication":True,
      "source_only_errors_by_year_and_type":[{"year":k[0],"reason":k[1],"n":v}
                                             for k,v in sorted(failures.items())],
      "no_river_TP_EQR_or_wwtp_values_loaded":True,
      "observations":results}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/f"FR178_2SCENE_MOSAIC_SOURCE_SHARD{j['source_shard_index']}.json").write_text(
          json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_FR178_MOSAIC_SOURCE="+json.dumps(
        {"shard":j["source_shard_index"],"n_source":len(results),"n_source_eligible":len(eligible)},
         ensure_ascii=False),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stac",type=Path,required=True)
    p.add_argument("--rivers",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:
        print({"status":"PASS_SYNTHETIC_MOSAIC_CLOUD_GAP_FILL"} if synthetic_smoke() else {})
        return
    run(args)

if __name__=="__main__":
    main()
