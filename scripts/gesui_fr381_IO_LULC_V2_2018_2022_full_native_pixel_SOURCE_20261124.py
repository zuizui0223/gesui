#!/usr/bin/env python3
"""Actual ORIGINAL 2018 and 2022 IO Annual LULC V2 paired 10m raster class source.

Across the SAME 381 source-locked French 2018/2022 EEA sites, calculate only
source-land-cover class fraction deltas (Trees/Rangeland/Crops/Built) with exact
paired cloud/nodata/land masks. BIOLOGICAL EQR NUMERIC IS NEVER READ here.
A hard SOURCE GO/STOP is emitted before any future ecology model.
"""
from __future__ import annotations
import argparse,collections,contextlib,hashlib,json,pathlib,sys,urllib.request
import numpy as np,rasterio
from rasterio.windows import from_bounds
from pyproj import Transformer
sys.path.insert(0,"scripts")
import gesui_FR381_TCD2018_2022_quality_tile_matched_source_20261124 as old
import gesui_fr381_IO_LULC2018_2022_real_site_UTM_source_20261124 as source
import gesui_fr381_IO_LULC_2018_2022_5fold_real_pixel_pilot_20261124 as pilot

COLL="io-lulc-annual-v02"
CLASSES={"trees":2,"flooded_vegetation":4,"crops":5,"built":7,"bare":8,"rangeland":11}
CLOUD=10;NODATA=0;WATER=1;SNOW=9
RINGS=(("near",0,250,200),("far",250,500,400))
SOURCE_MIN_VALID_SITES=300
MIN_FOLDS=30
SOURCE_CLOUD_LIMIT=.05
def original_data(args):
 prev=json.loads(args.site_qa.read_text())
 oldpilot=json.loads(args.pilot.read_text())
 if (prev.get("schema")!="GESUI_FR381_IO_LULC_V2_2018_22_EEA_RZ_ORIGINAL_SOURCE_SITE_MGRS_MAP_v1" or
     prev.get("source_geographical_GO") is not True or
     oldpilot.get("schema")!="GESUI_FR381_5FOLD_IO_LULC_V2_ORIGINAL_2018_2022_NATIVE_PIXEL_SOURCE_QA_v1"):
  raise RuntimeError("Original official EEA source-pinned site or source pixel pilot missing")
 fr,rz=old.sites(args)
 tile_source,meta=source.source_items()
 if len(fr)!=381 or sorted(tile_source)!=sorted(prev["original_UTM_supercell_STAC_source_items"]):
  raise RuntimeError("Source exact 381 French sites or source STAC MGRS original tile keys drift")
 matched=pilot.site_candidates(fr,tile_source)
 if len(matched)!=381:raise RuntimeError("Source paired all-381 original sites changed")
 oldcount=prev["source_france_sites_by_original_MGRS_tiles"]["per_supercell_number_of_sites"]
 now=collections.Counter(x["prefix"] for x in matched)
 if dict(sorted(now.items()))!=oldcount:
  raise RuntimeError("Exact source 381 MGRS native V2 site distribution changed")
 return fr,matched,oldpilot
def source_urls(prefixes,token):
 out={}
 for prefix in prefixes:
  out[prefix]={}
  for yr in ("2018","2022"):
   url,metadata=pilot.native_COG(prefix,yr,token)
   out[prefix][yr]=(url,metadata)
 return out
def native_read(ds2018,ds2022,site):
 x1,y1=Transformer.from_crs(4326,ds2018.crs,always_xy=True).transform(site["lon"],site["lat"])
 if not (ds2018.crs==ds2022.crs and ds2018.transform==ds2022.transform and
         ds2018.width==ds2022.width and ds2018.height==ds2022.height and
         ds2018.count==ds2022.count==1):
  raise RuntimeError("Native original 2018/2022 source COG class grids misaligned")
 if not (ds2018.res[0]>9.9 and ds2018.res[0]<10.1 and ds2018.res[1]>9.9 and ds2018.res[1]<10.1):
  raise RuntimeError("Source original V2 10m pixel resolution not nominal 10m")
 if not(ds2018.bounds.left+505<x1<ds2018.bounds.right-505 and ds2018.bounds.bottom+505<y1<ds2018.bounds.top-505):
  return {"status":"FAIL_BUFFER_CROSSES_NATIVE_TILE_EDGE","eligible":False,"quality_reason":"full500m not within chosen native tile"}
 window=from_bounds(x1-505,y1-505,x1+505,y1+505,ds2018.transform).round_offsets().round_lengths()
 try:
  d18=ds2018.read(1,window=window)
  d22=ds2022.read(1,window=window)
 except Exception as e:
  return {"status":"FAIL_HTTP_NATIVE_RASTER_READ","eligible":False,"quality_reason":type(e).__name__}
 if d18.shape!=d22.shape or not (99<=d18.shape[0]<=105 and 99<=d18.shape[1]<=105):
  raise RuntimeError("Source 2018/22 class pixel window rows/cols not spatially matching")
 if not set(np.unique(d18)).issubset({0,1,2,4,5,7,8,9,10,11}) or not set(np.unique(d22)).issubset({0,1,2,4,5,7,8,9,10,11}):
  return {"status":"FAIL_ORIGINAL_UNDEFINED_LAND_CLASS","eligible":False,"quality_reason":"invalid native class code"}
 tx=ds2018.window_transform(window)
 rr,cc=np.indices(d18.shape)
 XX=tx.c+tx.a*(cc+.5)+tx.b*(rr+.5)
 YY=tx.f+tx.d*(cc+.5)+tx.e*(rr+.5)
 dist=(XX-x1)**2+(YY-y1)**2
 out={"eligible":True,"native_2018_2022_raster_geospatial_grid_equal":True,
      "status":"PASS_ORIGINAL_PAIRED_SOURCE_CLASS_GRID",
      "original_2018_2022_pair_raster_code_support_equal_UTM":True}
 for band,mi,ma,needed in RINGS:
  mask=(dist>mi*mi)&(dist<=ma*ma)
  terra=(~np.isin(d18,[WATER,SNOW]))&(~np.isin(d22,[WATER,SNOW]))&mask
  q=terra&(~np.isin(d18,[CLOUD,NODATA]))&(~np.isin(d22,[CLOUD,NODATA]))
  n=int(q.sum());land=int(terra.sum())
  loss=1-n/max(1,land)
  cls={}
  for name,cat in CLASSES.items():
   cls[name+"_frac_2018"]=float(np.mean(d18[q]==cat)) if n else None
   cls[name+"_frac_2022"]=float(np.mean(d22[q]==cat)) if n else None
   cls[name+"_delta_frac"]=float(np.mean(d22[q]==cat)-np.mean(d18[q]==cat)) if n else None
  passes=n>=needed and land>0 and loss<=SOURCE_CLOUD_LIMIT
  out[band]={"source_pair_annulus_m":[mi,ma],"source_original_ring_pixels":int(mask.sum()),
    "source_original_2year_land_candidate_pixel_count":land,
    "source_paired_year_valid_pixels":n,"source_cloud_nodata_loss_fraction":float(loss),
    "source_expected_min_valid_pixels":needed,"source_max_missing_fraction":SOURCE_CLOUD_LIMIT,
    "source_2018_2022_native_categorical_pixel_quality_pass":bool(passes),
    "class":cls}
  out["eligible"]=out["eligible"] and passes
 if not out["eligible"]:out["status"]="SOURCE_PIXEL_QA_FAIL_STRICT_ORIGINAL_GATE"
 return out
def run(args):
 frame,matched,previous_pilot=original_data(args)
 token=pilot.GET("https://planetarycomputer.microsoft.com/api/sas/v1/token/"+COLL).get("token")
 if not token:raise RuntimeError("Original open PlanetaryComputer public SAS token not returned")
 selected=sorted(set(x["prefix"] for x in matched))
 refs=source_urls(selected,token)
 features=[]
 problems=collections.Counter()
 per_fold=collections.Counter()
 by_reg=collections.Counter()
 tile_stats={}
 with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
       CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",GDAL_HTTP_MAX_RETRY="2",
       GDAL_HTTP_RETRY_DELAY="1",GDAL_CACHEMAX=256):
  with contextlib.ExitStack() as stack:
   raster={}
   for code in selected:
    try:
     r18=stack.enter_context(rasterio.open(refs[code]["2018"][0]))
     r22=stack.enter_context(rasterio.open(refs[code]["2022"][0]))
    except Exception as e:
     # DO NOT leak signed public SAS query in exception messages.
     raise RuntimeError("Original V2 native public class raster opening failed: "+type(e).__name__) from None
    if r18.crs!=r22.crs or r18.transform!=r22.transform or r18.width!=r22.width or r18.height!=r22.height:
     raise RuntimeError("Original 2018 and 2022 IO V2 native COG source grids not identical for "+code)
    raster[code]=(r18,r22)
    tile_stats[code]={"native_original_CRS":str(r18.crs),
      "source_dimensions_wh":[r18.width,r18.height],
      "year2018_source_original_STAC_ID":refs[code]["2018"][1]["source_item_ID"],
      "year2022_source_original_STAC_ID":refs[code]["2022"][1]["source_item_ID"],
      "source_native_pair_2018_2022_geotransform_exact":True}
   sourceframe=frame.set_index("site_key")
   for i,s in enumerate(matched):
    site_id=s["site_key"]
    where=sourceframe.loc[site_id]
    fold=int(s["fold"])
    region=str(where["region_3deg_block"])
    a,b=raster[s["prefix"]]
    r=native_read(a,b,s)
    problems[r["status"]]+=1
    if r["eligible"]:
     per_fold[fold]+=1
     by_reg[region]+=1
    item={"site_hash_sha256":hashlib.sha256(site_id.encode()).hexdigest(),
      "original_2018_fixed_source_spatial_fold":fold,
      "original_2018_fixed_source_3deg_region":region,
      "source_original_MGRS_supercell":s["prefix"],
      "source_quality_pass":bool(r["eligible"]),
      "source_pixel_status":r["status"],
      "real_2018_2022_category_fraction_source":r}
    features.append(item)
    if (i+1)%75==0:
     print("FR381_IO_LULC_SOURCE_QUALITY_PROGRESS",json.dumps({
      "source_total_sites_visited":i+1,
      "source_eligible_so_far":sum(per_fold.values()),
      "source_fail_types_so_far":dict(problems)},sort_keys=True),flush=True)
 if len(features)!=381 or len({x["site_hash_sha256"] for x in features})!=381:
  raise RuntimeError("Original France381 source sites no longer unique")
 eligible=sum(per_fold.values())
 geo_source_GO=eligible>=SOURCE_MIN_VALID_SITES and min(per_fold.values() or [0])>=MIN_FOLDS and len(by_reg)>=10 and all(per_fold[g]>=MIN_FOLDS for g in range(1,6))
 outputs={
  "schema":"GESUI_REAL_FR381_2018_2022_IO_10M_CLASS_SOURCE_YEARPAIRED_PIXELS_V1",
  "original_data_provider":"Impact Observatory Esri Microsoft 10m annual 9-class V2 (CC-BY4.0)",
  "provider_original_STAC_collection":COLL,
  "native_genuine_category_codes":CLASSES,
  "source_original_EEA_FR381_sites_before_outcome":381,
  "source_original_years":[2018,2022],
  "original_EEA_site_2018_and_2022_numeric_EQR_NOT_OPENED":True,
  "source_native_2018_2022_same_grid_for_each_UTM_tile":True,
  "source_mapped_mgsr_tile_source":tile_stats,
  "source_france_original_site_quality_eligible":eligible,
  "source_original_site_pixel_quality_failures_by_reason":dict(problems),
  "source_original_2018_fixed_fold_valid_site_counts":dict(sorted(per_fold.items())),
  "source_valid_original_2018_3deg_region_count":len(by_reg),
  "source_full_census_original_site_quality_records":features,
  "source_predeclared_min_sites":SOURCE_MIN_VALID_SITES,
  "source_predeclared_min_each_5fold":MIN_FOLDS,
  "source_predeclared_min_regions":10,
  "source_predeclared_cloud_nodata_loss_limit":SOURCE_CLOUD_LIMIT,
  "source_strict_high_quality_go":bool(geo_source_GO),
  "source_status":"SOURCE_GO_ORIGINAL_LULC2018_2022_SOURCE_QUALITY_PIXEL_PASS_ECOLOGY_NEXT" if geo_source_GO else "SOURCE_STOP_ORIGINAL_LULC2018_2022_PIXEL_QUALITY_SUPPORT",
  "no_POTW_biology_or_association_fitted":True}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(outputs,sort_keys=True,indent=2)+"\n")
 print("FR381_REAL_IO_LULC_2018_2022_SOURCE_FULL_CENSUS",json.dumps({
  "status":outputs["source_status"],"original_source_sites":381,
  "valid_source_sites":eligible,
  "folds":outputs["source_original_2018_fixed_fold_valid_site_counts"],
  "regions":outputs["source_valid_original_2018_3deg_region_count"],
  "failures":outputs["source_original_site_pixel_quality_failures_by_reason"],
  "original_source_tiles":tile_stats},sort_keys=True),flush=True)
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq18","eq22","site_qa","pilot","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 run(p.parse_args())
if __name__=="__main__":main()
