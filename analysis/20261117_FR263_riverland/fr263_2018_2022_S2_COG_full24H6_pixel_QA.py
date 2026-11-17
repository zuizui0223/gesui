#!/usr/bin/env python3
"""Original France 24-H6 2018/2022 real S2 COG clear-land pixel pilot.

A source-quality study, NO river TP/EQR values opened or modeled. This pilot
uses original field-station-centered 0-250m and 250-500m radial rings; it is
NOT a georeferenced riverline buffer nor evidence of wastewater nitrogen flux.
Only 2021 static WorldCover cover classes, not 2018 vegetation type stability.
"""
from __future__ import annotations
import argparse, collections, hashlib, json, math
from pathlib import Path
from contextlib import ExitStack
import numpy as np
import rasterio
from rasterio.warp import reproject
from rasterio.windows import from_bounds, Window
from rasterio.enums import Resampling
from pyproj import Transformer

INPUT_SCHEMA="GESUI_20261117_FR263_S2_2018_2022_24H6_PILOT_SOURCE_OPPORTUNITY_V1"
OUTPUT_SCHEMA="GESUI_20261117_FR263_S2_2018_2022_TRUE_10M_PIXEL_FULL24H6_SOURCE_QA_V1"
INPUT_SHA="22414dec230bc9bf9c308d823f166303cc843143c6ade922af6cb44b1f2c3f3c"
COG_HOST="https://sentinel-cogs.s3.us-west-2.amazonaws.com/"
WC_HOST="https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
NATURAL=(10,20,30,90)
MINS={"near_SCL":100,"far_SCL":200,"near_natural":40,"far_natural":100,"per_class":10}
MONTHS=(6,7,8)
YEAR1=2018
YEAR2=2022
NEAR=250.
FAR=500.
SMOKE_H6=24
SMOKE_PASS=18

def sha(b):return hashlib.sha256(b).hexdigest()

def load(path):
    raw=path.read_bytes()
    if sha(raw)!=INPUT_SHA:
        raise RuntimeError("Previously signed FR263 24H6 STAC metadata source SHA mismatched")
    data=json.loads(raw)
    if data.get("schema")!=INPUT_SCHEMA or data.get("source_supported_bothyear_pilot_n")!=24 or len(data.get("pilot_sites",[]))!=24:
        raise RuntimeError("Original full24 same-station 2018/2022 STAC source gate no longer passes")
    r=data["pilot_sites"]
    if len({x["H6"] for x in r})!=24:
        raise RuntimeError("Original one-per-H6 source identity changed")
    if any(not x["source_supported_both_2018_and2022_two_summer_months"] for x in r):
        raise RuntimeError("Old signed original source said 24/24 but site entries changed")
    # Hash H6 identities only, not image greenness, river TP, EQR or WWTP
    selected=sorted(r,key=lambda x:(sha(x["H6"].encode()),x["H6"]))[:SMOKE_H6]
    if len(selected)!=SMOKE_H6:raise RuntimeError("Preselected 24 H6 not available")
    return selected

def worldcover(lat,lon):
    la=math.floor(lat/3)*3;lo=math.floor(lon/3)*3
    tile=("N" if la>=0 else "S")+f"{abs(la):02d}"+("E" if lo>=0 else "W")+f"{abs(lo):03d}"
    return WC_HOST+f"ESA_WorldCover_10m_2021_v200_{tile}_Map.tif",tile

def asset(scene,k):
    v=scene["bands"][k]
    url=v["href"]
    if not isinstance(url,str) or not url.startswith(COG_HOST) or not url.endswith(".tif"):
        raise RuntimeError("Original Sentinel2 COG href changed from signed source")
    if k=="scl":return url,None,None
    rb=v["raster_bands"]
    if len(rb)!=1 or "scale" not in rb[0]:raise RuntimeError("Actual 2018/22 original radiometry unspecified")
    sc=float(rb[0]["scale"]);offset=float(rb[0].get("offset",0))
    if not (0<sc<=1 and -1<=offset<=1):raise RuntimeError("Original reflectance scale/offset invalid")
    return url,sc,offset

def scene_pixels(scene,site):
    vals={k:asset(scene,k) for k in ("red","nir","scl")}
    if (abs(vals["red"][1]-vals["nir"][1])>1e-12 or
        abs(vals["red"][2]-vals["nir"][2])>1e-12):
        raise RuntimeError("Original red/nir source scale/offset disagree")
    lat=float(site["original_lat"]);lon=float(site["original_lon"])
    wc_url,wctile=worldcover(lat,lon)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
      CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_HTTP_MAX_RETRY="2",
      GDAL_HTTP_RETRY_DELAY="1",GDAL_CACHEMAX=64,VSI_CACHE="TRUE",VSI_CACHE_SIZE=1000000):
      with ExitStack() as stack:
        img={k:stack.enter_context(rasterio.open("/vsicurl/"+vals[k][0]))
               for k in ("red","nir","scl")}
        wc=stack.enter_context(rasterio.open("/vsicurl/"+wc_url))
        red=img["red"];nir=img["nir"];scl=img["scl"]
        if (red.crs is None or nir.crs!=red.crs or scl.crs!=red.crs or wc.crs is None or
          red.transform!=nir.transform or red.width!=nir.width or red.height!=nir.height):
            raise RuntimeError("Original 10m reflectance assets no longer co-registered")
        if not (8<=abs(red.transform.a)<=12 and 8<=abs(red.transform.e)<=12):
            raise RuntimeError("Source red reflectance native grid not 10m")
        p=Transformer.from_crs("EPSG:4326",red.crs,always_xy=True)
        x,y=p.transform(lon,lat)
        win=from_bounds(x-515,y-515,x+515,y+515,red.transform)
        win=Window(int(math.floor(win.col_off)),int(math.floor(win.row_off)),
                   int(math.ceil(win.width))+2,int(math.ceil(win.height))+2)
        rr=red.read(1,window=win,boundless=True,fill_value=0)
        nn=nir.read(1,window=win,boundless=True,fill_value=0)
        t=red.window_transform(win)
        ss=np.zeros(rr.shape,dtype=np.uint8)
        cc=np.zeros(rr.shape,dtype=np.uint8)
        reproject(source=rasterio.band(scl,1),destination=ss,
                  dst_transform=t,dst_crs=red.crs,
                  src_nodata=0,dst_nodata=0,resampling=Resampling.nearest)
        reproject(source=rasterio.band(wc,1),destination=cc,
                  dst_transform=t,dst_crs=red.crs,
                  src_nodata=0,dst_nodata=0,resampling=Resampling.nearest)
        yi,xi=np.indices(rr.shape)
        xx=t.c+(xi+.5)*t.a+(yi+.5)*t.b
        yy=t.f+(xi+.5)*t.d+(yi+.5)*t.e
        dist2=(xx-x)**2+(yy-y)**2
        near=dist2<=NEAR**2
        far=(dist2>NEAR**2)&(dist2<=FAR**2)
        r=rr.astype(np.float32)*vals["red"][1]+vals["red"][2]
        n=nn.astype(np.float32)*vals["nir"][1]+vals["nir"][2]
        physically_valid=(r>.005)&(n>.005)&(r<1)&(n<1)
        clear=(np.isin(ss,[4,5]) & physically_valid)
        natural=(np.isin(cc,NATURAL) & clear)
        counts={"near_total":int(near.sum()),"far_total":int(far.sum()),
           "near_clear":int((near&clear).sum()),"far_clear":int((far&clear).sum()),
           "near_natural":int((near&natural).sum()),"far_natural":int((far&natural).sum())}
        if (counts["near_clear"]<MINS["near_SCL"] or counts["far_clear"]<MINS["far_SCL"] or
            counts["near_natural"]<MINS["near_natural"] or counts["far_natural"]<MINS["far_natural"]):
            return {"pixel_source_pass":False,"failure":"insufficient_real_cloudfree_land_or_natural_pixels",
                    "counts":counts,"worldcover_tile":wctile}
        ndvi=(n-r)/(n+r+1e-12)
        valid=(ndvi>=-1)&(ndvi<=1)&np.isfinite(ndvi)
        matches=[]
        for k in NATURAL:
            mask=(cc==k)&natural&valid
            i=near&mask;o=far&mask
            ni=int(i.sum());no=int(o.sum())
            if ni<MINS["per_class"] or no<MINS["per_class"]:continue
            matches.append({"class":k,"n_near":ni,"n_far":no,
                "median_near":float(np.median(ndvi[i])),
                "median_far":float(np.median(ndvi[o]))})
        return {"pixel_source_pass":bool(matches),
          "failure":None if matches else "no_same_natural_cover_class_both_original_rings",
          "counts":counts,"original_scene_reflectance_scale":vals["red"][1],
          "original_scene_reflectance_offset":vals["red"][2],
          "worldcover_tile":wctile,
          "cover_matched_original_real_pixel_NDVI":matches}

def try_month(site,year,month):
    scenes=sorted((v for v in site["years"][str(year)]["source_selected_scenes"]
          if v["month"]==month),
          key=lambda x:(x["cloud_pct"],x["datetime"],x["scene_id"]))
    if len(scenes)>2:raise RuntimeError("Original source scene attempted >2 per month")
    attempts=[]
    for x in scenes:
        try:r=scene_pixels(x,site)
        except Exception as e:
            r={"pixel_source_pass":False,"failure":"original_public_COG_access_or_raster_error",
               "error_type":type(e).__name__,"error":str(e)[:250]}
        rec={"year":year,"month":month,"scene":x["scene_id"],"datetime":x["datetime"],
             "scene_tile_cloud_pct":x["cloud_pct"],**r}
        attempts.append(rec)
        if r.get("pixel_source_pass"):
            return rec,attempts
    return None,attempts

def summarize_pair(yr1,yr2):
    a={c["class"]:c for c in yr1["cover_matched_original_real_pixel_NDVI"]}
    b={c["class"]:c for c in yr2["cover_matched_original_real_pixel_NDVI"]}
    common=sorted(set(a)&set(b))
    if not common:return None
    weights=np.array([min(a[k]["n_near"],b[k]["n_near"]) for k in common],float)
    weights=weights/weights.sum()
    p=lambda z:float(np.dot(weights,np.array([z[k]["median_near"]-z[k]["median_far"] for k in common])))
    return {"same_calendar_month":yr1["month"],"common_worldcover_classes":common,
      "ndvi_near_minus_far_2018":p(a),"ndvi_near_minus_far_2022":p(b),
      "delta_near_minus_far_2022_minus_2018":p(b)-p(a),
      "source_2021_common_cover_weight_from_min_2018_2022_clear_near_pixel_count":weights.tolist()}

def run(args):
    selected=load(args.stac)
    full=[]
    for i,site in enumerate(selected,1):
        months=[]
        attempted=[]
        for month in MONTHS:
            a,qa=try_month(site,YEAR1,month)
            b,qb=try_month(site,YEAR2,month)
            attempted.extend(qa+qb)
            if a is not None and b is not None:
                pair=summarize_pair(a,b)
                if pair:months.append(pair)
        ok=len(months)>=2
        row={"H6":site["H6"],"H7":site["H7"],"site_key":site["site_key"],
          "original_HYRIV_ID":site["HYRIV_ID"],
          "n_matched_2018_2022_calendar_months":len(months),
          "source_valid_2_month_bothyears":ok,
          "matched_months":months,
          "year1_year2_scene_pixel_all_attempts":attempted,
          "real_biotic_TP_EQR_not_loaded":True}
        if months:
            row["descriptive_delta_NDVI_near_minus_far_2022_minus_2018"]=float(
              np.median([z["delta_near_minus_far_2022_minus_2018"] for z in months]))
        full.append(row)
        print(f"GESUI_FR24H6_S2_TRUE_PIXEL_FIXED8_SITE={i}/{len(selected)} H6={site['H6']} pairedmonths={len(months)} gate={ok}",flush=True)
    pass_sites=[x for x in full if x["source_valid_2_month_bothyears"]]
    report={"schema":OUTPUT_SCHEMA,
      "status":"PASS_REAL_S2_2018_2022_PIXEL_FULL24_NO_ECOLOGY" if len(pass_sites)>=SMOKE_PASS else
               "STOP_REAL_S2_2018_2022_PIXEL_FULL24_COVERAGE_NOT_BIOLOGICAL_NULL",
      "input_metadata_sha256":INPUT_SHA,"original_FR263_H6_source_universe":24,
      "original_H6_fixed_24_by_SHA256_id":[x["H6"] for x in selected],
      "n_fixed_original_H6_tested":24,"n_two_years_real_cloudfree_matched_months_eligible":len(pass_sites),
      "predeclared_next_full24_min_pass":SMOKE_PASS,
      "optical_2018_2022_bands_with_each_scene_original_scale_and_offset":True,
      "source_static_WorldCover2021_10m_natural_classes":list(NATURAL),
      "source_SCL_clear_classes":[4,5],
      "predeclared_pixel_thresholds":MINS,
      "fixed_rings_field_station_center_m":[0,250,500],
      "fixed_months":[6,7,8],"both_years_same_calendar_month_comparisons":True,
      "full24_source_still_not_109_riparian_cover_matched_station_cohort":True,
      "original_riverline_distance_not_measured":True,
      "this_is_remote_NDVI_not_genuine_biomass_or_plant_species_diversity":True,
      "2021_landclass_used_for_2018_2022_mask_not_true_two_year_landuse_stability":True,
      "source_no_river_TP_EQR_or_wastewater_effect_estimated":True,
      "source_fixed_24_site_results":full}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"FR263_S2_2018_2022_24H6_REAL_COG_PIXEL_OPPORTUNITY_QA.json").write_text(
         json.dumps(report,indent=2,allow_nan=False,ensure_ascii=False)+"\n")
    print("GESUI_FR263_REAL_PIXEL_FULL24_RESULT="+json.dumps({
      "status":report["status"],"n_stations":SMOKE_H6,"n_two_years_matched":len(pass_sites),
      "n_eligible_for_next_stage":int(len(pass_sites)>=SMOKE_PASS)}),flush=True)
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stac",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    run(p.parse_args())

if __name__=="__main__":
    main()
