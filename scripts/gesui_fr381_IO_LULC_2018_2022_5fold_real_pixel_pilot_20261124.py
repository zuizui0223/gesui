#!/usr/bin/env python3
"""Original public signed Planetary Computer V2 annual 10m LULC numeric pixel QA.

Exactly one source-hash-picked France original EEA monitoring site from each
of five frozen 2018 geographical folds, 2018+2022 pair. Original EQR numeric
response never read. Native categorical class values, not styled PNG.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib,sys,collections,urllib.parse,urllib.request
import numpy as np
import rasterio
from rasterio.windows import from_bounds
from pyproj import Transformer
sys.path.insert(0,"scripts")
import gesui_FR381_TCD2018_2022_quality_tile_matched_source_20261124 as old
import gesui_fr381_IO_LULC2018_2022_real_site_UTM_source_20261124 as source

COLL="io-lulc-annual-v02"
STAC="https://planetarycomputer.microsoft.com/api/stac/v1"
SAS="https://planetarycomputer.microsoft.com/api/sas/v1/token/"+COLL
YEARS=("2018","2022")
CODES={0,1,2,4,5,7,8,9,10,11}
MIN_NEAR=200
MIN_FAR=400
def GET(url):
 req=urllib.request.Request(url,headers={"User-Agent":"gesui impact observatory native V2 France original 10m quality gate/20261124","Accept":"application/json"})
 with urllib.request.urlopen(req,timeout=35) as r:
  x=json.loads(r.read(120000))
  if r.status!=200:raise RuntimeError("Original publisher STAC source status not200")
  return x
def site_candidates(fr,source_items):
 out=[]
 for _,r in fr.iterrows():
  lon=float(r.lon_2018);lat=float(r.lat_2018);fold=int(r.pre_2018_source_only_spatial_fold)
  zone=int((lon+180)//6)+1
  prefix=[]
  for item,yearmap in source_items.items():
   A=yearmap.get("2018",{});B=yearmap.get("2022",{})
   if (A.get("source_utm_epsg")!=32600+zone or B.get("source_utm_epsg")!=32600+zone or
       A.get("status")!="ORIGINAL_SOURCE_STAC_EXACT_ITEM" or
       B.get("status")!="ORIGINAL_SOURCE_STAC_EXACT_ITEM"):continue
   if source.covers(A["source_publisher_bbox"],lon,lat,True):prefix.append(item)
  if len(prefix)!=1:raise RuntimeError("Original real France site source UTM MGRS class pairing ambiguous")
  out.append({"site_key":r.site_key,"lon":lon,"lat":lat,"fold":fold,"prefix":prefix[0]})
 return out
def verify_source(args):
 source_qa=json.loads(args.site_qa.read_text())
 if source_qa.get("schema")!="GESUI_FR381_IO_LULC_V2_2018_22_EEA_RZ_ORIGINAL_SOURCE_SITE_MGRS_MAP_v1" or not source_qa["source_geographical_GO"]:
  raise RuntimeError("Real source France381 2018/22 UTM geogate not signed")
 fr,rz=old.sites(args)
 original,src=source.source_items()
 if len(fr)!=381 or len(original)!=len(source_qa["original_UTM_supercell_STAC_source_items"]):
  raise RuntimeError("Original true EEA source or public STAC changed")
 if sorted(original)!=sorted(source_qa["original_UTM_supercell_STAC_source_items"]):
  raise RuntimeError("Source original published STAC item key list changed")
 candidate=site_candidates(fr,original)
 if len(candidate)!=381:raise RuntimeError("Original site source census 381 changed")
 selected=[]
 for fold in range(1,6):
  options=[x for x in candidate if x["fold"]==fold]
  options.sort(key=lambda q:hashlib.sha256(q["site_key"].encode()).hexdigest())
  chosen=options[0]
  if chosen["prefix"]!=source_qa["source_france_sites_by_original_MGRS_tiles"]["source_hash_selected_native_cell_per_fivefold"][str(fold)]["original_source_supercell"]:
   raise RuntimeError("Original fixed 5fold site original source tile selection drift")
  selected.append(chosen)
 return selected
def native_COG(prefix,year,token):
 item=GET(STAC+"/collections/"+COLL+"/items/"+prefix+"-"+year)
 if item.get("id")!=prefix+"-"+year or item.get("collection")!=COLL:raise RuntimeError("Original STAC true-year native item changed")
 dat=(item.get("assets") or {}).get("data") or {}
 if "tiff" not in str(dat.get("type","")).lower():raise RuntimeError("Original native asset no longer raw categorical TIFF")
 href=dat.get("href","")
 if not href.startswith("https://") or "io-lulc" not in href:raise RuntimeError("Original source COG asset missing")
 sep="&" if "?" in href else "?"
 return href+sep+token,{"source_item_ID":item.get("id"),"source_UTM_epsg":item.get("properties",{}).get("proj:epsg"),
                         "source_asset_host":urllib.parse.urlparse(href).netloc,
                         "source_asset_path":urllib.parse.urlparse(href).path}
def read_tiny(site,prefix,token):
 rasters={}
 meta={}
 for year in YEARS:
  url,detail=native_COG(prefix,year,token)
  # Rasterio only reads source HTTP byte-ranges for local 500m disk.
  # Never print signed URL/exception message (contains signed token).
  with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
       GDAL_HTTP_MAX_RETRY="2",GDAL_HTTP_RETRY_DELAY="1",GDAL_CACHEMAX=64):
   try:
    with rasterio.open(url) as ds:
     if ds.count!=1 or ds.crs is None or abs(ds.res[0]-10)>0.1 or abs(ds.res[1]-10)>0.1:
      raise RuntimeError("Original native raster layer not single original 10m categorical class")
     transform=Transformer.from_crs(4326,ds.crs,always_xy=True)
     x,y=transform.transform(site["lon"],site["lat"])
     if not(ds.bounds.left+510<x<ds.bounds.right-510 and ds.bounds.bottom+510<y<ds.bounds.top-510):
      raise RuntimeError("Original same-UTM native source raster grid cannot support 500m buffer")
     w=from_bounds(x-505,y-505,x+505,y+505,ds.transform).round_offsets().round_lengths()
     if w.width<99 or w.width>105 or w.height<99 or w.height>105:
      raise RuntimeError("Original native source raster window pixel size not 10m")
     d=ds.read(1,window=w)
     ww=ds.window_transform(w)
     rows,cols=np.indices(d.shape)
     XX=ww.c+ww.a*(cols+.5)+ww.b*(rows+.5)
     YY=ww.f+ww.d*(cols+.5)+ww.e*(rows+.5)
     circle=(XX-x)**2+(YY-y)**2
     if not set(np.unique(d)).issubset(CODES):
      raise RuntimeError("Original IO categorical raster has undefined LULC code")
     rasters[year]={"image":d,"distance2":circle,"crs":str(ds.crs),
       "source_geotransform":list(ds.transform.to_gdal()),
       "source_raster_width":ds.width,"source_raster_height":ds.height,
       "source_native_nodata":ds.nodata,
       "center_UTM_epsg":ds.crs.to_epsg()}
     meta[year]={**detail,"native_CRS":str(ds.crs),
       "native_resolution_m":list(ds.res),"native_tile_size":[ds.width,ds.height],
       "native_nodata":ds.nodata,
       "source_unique_category_codes_for_tiny_disk":sorted(int(v) for v in np.unique(d))}
   except Exception as e:
    # No secrets in logs from rasterio exceptions.
    raise RuntimeError("Original public signed raster tiny window failure "+type(e).__name__) from None
 a=rasters["2018"];b=rasters["2022"]
 if (a["crs"]!=b["crs"] or a["source_geotransform"]!=b["source_geotransform"] or
     a["image"].shape!=b["image"].shape or a["source_raster_width"]!=b["source_raster_width"]):
  raise RuntimeError("2018/22 original V2 UTM raster grid not truly aligned!")
 outputs={}
 for name,minr,maxr in (("near_0to250",0,250),("far_250to500",250,500)):
  within=(a["distance2"]>minr*minr)&(a["distance2"]<=maxr*maxr)
  v1=a["image"];v2=b["image"]
  terrestrial=(~np.isin(v1,[1,9]))&(~np.isin(v2,[1,9]))&within
  valid=terrestrial&(~np.isin(v1,[0,10]))&(~np.isin(v2,[0,10]))
  total=int(terrestrial.sum())
  n=int(valid.sum())
  if n < (MIN_NEAR if name.startswith("near") else MIN_FAR):
   raise RuntimeError("Source fixed 500m ring lacks valid original 10m land class pixels")
  cloud_loss=1-(n/max(1,total))
  if cloud_loss>.05:raise RuntimeError("Original paired cloud+nodata mask removes >5 pct")
  outputs[name]={"original_pixel_circles_source":int(within.sum()),
    "two_year_terrestrial_candidate_pixels":total,
    "native_both_year_accepted_land_pixels":n,
    "source_annual_cloud_nodata_loss_pct":100*cloud_loss,
    "original_tree_percent_2018_on_fixed_mask":float(100*np.mean(v1[valid]==2)),
    "original_tree_percent_2022_on_fixed_mask":float(100*np.mean(v2[valid]==2)),
    "original_tree_class_delta_percent_points":float(100*(np.mean(v2[valid]==2)-np.mean(v1[valid]==2)))}
 return meta,outputs
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq18","eq22","site_qa","out"):p.add_argument("--"+k,type=pathlib.Path,required=True)
 args=p.parse_args()
 fixed=verify_source(args)
 signed=GET("https://planetarycomputer.microsoft.com/api/sas/v1/token/"+COLL)
 token=signed.get("token")
 if not isinstance(token,str) or len(token)<20:raise RuntimeError("Original public SAS API no valid short-lived token")
 output={}
 for site in fixed:
  fold=str(site["fold"])
  meta,frac=read_tiny(site,site["prefix"],token)
  output[fold]={"source_original_UTM_supercell":site["prefix"],
                "native_original_annual_raster_meta":meta,
                "source_true_paired_year_land_category_tiny_disk":frac}
 report={"schema":"GESUI_FR381_5FOLD_IO_LULC_V2_ORIGINAL_2018_2022_NATIVE_PIXEL_SOURCE_QA_v1",
  "original_FRANCE_same_site_2018_22_source_sites":381,
  "preoutcome_source_fixed_five_fold_geo_sites":5,
  "original_source_2018_2022_STAC_collection":COLL,
  "original_public_SAS_evidence_only_token_not_logged":True,
  "actual_source_native_raw_categorical_10m_class_values_read":True,
  "original_fixed_five_fold_site_raster_QA":output,
  "source_both_year_native_UTM_grid_identical":True,
  "source_source_only_one_HUC_EEA_site_per_fixed_fold":True,
  "source_still_MISSING_all_381_sites_quality_and_class_counts":True,
  "original_EEA_2018_2022_EQR_numeric_READ":False,
  "source_status":"SOURCE_ORIGINAL_FIVE_FOLD_RAW_PIXEL_QA_GO_NOT_YET_381_SITE_SOURCE"}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FR381_IO_LULC_2018_2022_PILOT_NATIVE_PIXEL_QA",json.dumps({
  "status":report["source_status"],
  "pilots":{k:{"supercell":v["source_original_UTM_supercell"],
    "near":v["source_true_paired_year_land_category_tiny_disk"]["near_0to250"],
    "far":v["source_true_paired_year_land_category_tiny_disk"]["far_250to500"],
    "rasters":{yr:{"epsg":m["source_UTM_epsg"],"resolution":m["native_resolution_m"],
                    "classes":m["source_unique_category_codes_for_tiny_disk"]} for yr,m in v["native_original_annual_raster_meta"].items()}}
    for k,v in output.items()}},sort_keys=True),flush=True)
if __name__=="__main__":main()
