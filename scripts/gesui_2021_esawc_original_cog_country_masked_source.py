#!/usr/bin/env python3
"""Original ESA WorldCover 2021 categorical COG land cover for fixed DE/AT WWTP grids.

Reads *only* source-frozen EEA 1degree cells and original public WorldCover
COG class pixels plus Natural Earth country polygon. Never loads 2022 butterfly
photo outcomes or local BMS sensitive sites. Source-only posthoc confounder QC.
"""
from __future__ import annotations
import argparse,contextlib,csv,hashlib,json,math,os,time,urllib.request
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_bounds
from rasterio.windows import from_bounds as window_from_bounds

EA2022_SHA="926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce"
NATURAL_EARTH="https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson"
AWS_PREFIX="https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
LAND_COVER_CLASSES={0:"NoData",10:"tree",20:"shrubland",30:"grassland",40:"cropland",50:"built_up",
                    60:"bare_sparse",70:"snow_ice",80:"water",90:"wetlands",95:"mangroves",100:"moss_lichen"}
VALID_CLASS_CODES=set(LAND_COVER_CLASSES)
COUNTRY_NAMES={"DE":"Germany","AT":"Austria"}
RES_MAIN=512
RES_SENS=256
SOURCE_COUNTS={"DE":32,"AT":21}
MAX_CELLS=53

def sha_bytes(b):return hashlib.sha256(b).hexdigest()
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for c in iter(lambda:f.read(1024*1024),b""):h.update(c)
    return h.hexdigest()

def get_natural_earth():
    req=urllib.request.Request(NATURAL_EARTH,headers={
        "User-Agent":"gesui-2021-original-ESA-worldcover-categorical-landcover-pilot/0.1",
        "Accept":"application/geo+json,application/json,*/*"})
    with urllib.request.urlopen(req,timeout=90) as r:raw=r.read(6000000)
    if len(raw)>5000000:raise ValueError("Unexpected 50m Natural Earth source size")
    obj=json.loads(raw)
    features=obj.get("features",[])
    country={}
    for c,name in COUNTRY_NAMES.items():
        hits=[x for x in features if (str(x.get("properties",{}).get("ISO_A2","")).strip()==c
                    or str(x.get("properties",{}).get("ADM0_A3","")).strip()==("DEU" if c=="DE" else "AUT"))
                    and x.get("geometry")]
        if len(hits)!=1:
            raise ValueError("Natural Earth admin0 country unique geometry missing: "+c+" hits "+str(len(hits)))
        country[c]=hits[0]["geometry"]
    return country,{"original_Natural_Earth_admin0_50m_geojson_URL":NATURAL_EARTH,
        "original_geojson_bytes":len(raw),"original_geojson_SHA256":sha_bytes(raw),
        "original_country_features":{c:COUNTRY_NAMES[c] for c in country},
        "original_country_boundary_is_not_2021_dynamic_border_or_exact_coastline":True}

def tile_id(lat,lon):
    tile_lat=math.floor(lat/3)*3
    tile_lon=math.floor(lon/3)*3
    if tile_lat<0 or tile_lon<0:raise ValueError("Expected northern/eastern European tile")
    return f"N{tile_lat:02d}E{tile_lon:03d}"

def class_stats(arr,geom,latlo,side):
    countrymask=rasterize([(geom,1)],out_shape=(side,side),
            transform=from_bounds(latlo[1],latlo[0],latlo[1]+1,latlo[0]+1,side,side),
            fill=0,all_touched=False,dtype="uint8").astype(bool)
    frac_country_mask=float(np.mean(countrymask))
    latvalues=latlo[0]+(1.-(np.arange(side,dtype=float)+.5)/side)
    rowcos=np.cos(np.radians(latvalues))
    weights=np.broadcast_to(rowcos[:,None],(side,side))
    masked=arr[countrymask]
    if not set(int(x) for x in np.unique(masked)).issubset(VALID_CLASS_CODES):
        raise ValueError("Original raw COG contains unknown WorldCover class IDs")
    present=countrymask&(arr!=0)
    terrestrial=present&(arr!=80)
    def area(mask):return float(np.sum(weights[mask]))
    fullweight=float(np.sum(weights))
    land_w=area(terrestrial)
    if land_w<=0:land_w=0
    frac_country_terrestrial_cell=land_w/fullweight
    km2=land_w*(111.32**2)/(side**2)
    class_fraction={}
    for classnum in (10,30,40,50,90):
        class_fraction[LAND_COVER_CLASSES[classnum]+"_share_of_country_terrestrial"]=(area(terrestrial&(arr==classnum))/land_w if land_w>0 else None)
    return {
      "source_tile_sampling_side":side,
      "country_border_mask_sample_points":int(np.sum(countrymask)),
      "country_terrestrial_land_sample_points":int(np.sum(terrestrial)),
      "country_border_cover_fraction_cell":frac_country_mask,
      "country_terrestrial_land_fraction_cell":frac_country_terrestrial_cell,
      "original_source_NODATA_fraction_within_country":float(np.mean(masked==0)) if len(masked)>0 else None,
      "original_source_water_fraction_within_country":float(np.mean(masked==80)) if len(masked)>0 else None,
      "country_terrestrial_land_area_km2_approx":km2,
      **class_fraction,
      "country_class_code_sample_counts":{str(int(k)):int(v) for k,v in zip(*np.unique(masked,return_counts=True))},
      "source_nodata_is_not_built_up_zero":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-manifest",required=True,type=Path)
    p.add_argument("--source-tile-csv",required=True,type=Path)
    p.add_argument("--source-tile-qc",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--max-cells",type=int,default=MAX_CELLS)
    args=p.parse_args()
    if sha(args.source_manifest)!=EA2022_SHA:raise ValueError("Original EEA 2022 source full data digest changed")
    qc=json.loads(args.source_tile_qc.read_text())
    if qc["source_original_sha256"]!=EA2022_SHA or qc["source_tile_CSV_SHA256"]!=sha(args.source_tile_csv):
        raise ValueError("Original EEA source-only raster grid not SHA locked")
    cells=list(csv.DictReader(args.source_tile_csv.open()))
    if len(cells)!=53 or Counter(c["countryCode"] for c in cells)!=Counter(SOURCE_COUNTS):
        raise ValueError("Original 2022 WWTP 53 source grid sample changed")
    if not(1<=args.max_cells<=MAX_CELLS):raise ValueError("Invalid fixed max cells")
    geo,geoQC=get_natural_earth()
    cached={}
    rows=[]
    try:
        with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
          CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
          GDAL_HTTP_MAX_RETRY="3",GDAL_HTTP_RETRY_DELAY="1",
          GDAL_HTTP_TIMEOUT="55",
          AWS_NO_SIGN_REQUEST="YES",
          VSI_CACHE="TRUE",VSI_CACHE_SIZE="33554432"):
          for i,c in enumerate(cells[:args.max_cells]):
            lat=float(c["lat_floor"]);lon=float(c["lon_floor"])
            country=c["countryCode"]
            tile=tile_id(lat,lon)
            srcurl=AWS_PREFIX+f"ESA_WorldCover_10m_2021_v200_{tile}_Map.tif"
            if tile not in cached:
                ds=rasterio.open(srcurl)
                if ds.crs.to_epsg()!=4326 or ds.count!=1 or ds.width<10000 or ds.height<10000 or (ds.nodata not in (0,None)):
                    raise ValueError("WorldCover original categorical COG has unexpected CRS, 10m shape or bands")
                if ds.bounds.left>lon or ds.bounds.bottom>lat or ds.bounds.right<lon+1 or ds.bounds.top<lat+1:
                    raise ValueError("Original 3 degree COG does not cover fixed 1 degree cell")
                cached[tile]=ds
                print("GESUI_WORLDCOVER2021_OFFICIAL_NEW_ORIGINAL_COG",json.dumps({
                   "original_url":srcurl,"CRS":str(ds.crs),"shape":[ds.height,ds.width],
                   "original_overview_decimations":ds.overviews(1),"nodata":ds.nodata}),flush=True)
            src=cached[tile]
            window=window_from_bounds(lon,lat,lon+1,lat+1,transform=src.transform)
            arr_main=src.read(1,window=window,out_shape=(RES_MAIN,RES_MAIN),
                    resampling=Resampling.nearest)
            arr_sens=src.read(1,window=window,out_shape=(RES_SENS,RES_SENS),
                    resampling=Resampling.nearest)
            if arr_main.shape!=(RES_MAIN,RES_MAIN) or arr_sens.shape!=(RES_SENS,RES_SENS):
                raise ValueError("Incomplete original COG window read")
            main=class_stats(arr_main,geo[country],(lat,lon),RES_MAIN)
            sens=class_stats(arr_sens,geo[country],(lat,lon),RES_SENS)
            deviations={}
            for k in ("built_up_share_of_country_terrestrial",
                  "cropland_share_of_country_terrestrial","tree_share_of_country_terrestrial",
                  "grassland_share_of_country_terrestrial"):
                v=main.get(k);z=sens.get(k)
                deviations[k]=abs(v-z) if v is not None and z is not None else None
            good=(main["country_terrestrial_land_fraction_cell"]>=0.05 and
                  main["country_terrestrial_land_sample_points"]>=512 and
                  main["original_source_NODATA_fraction_within_country"] is not None and
                  main["original_source_NODATA_fraction_within_country"]<=.05)
            stable=(all(deviations[k] is not None and deviations[k]<=.10
                    for k in ("built_up_share_of_country_terrestrial","cropland_share_of_country_terrestrial")))
            row={
             "countryCode":country,"lat_floor":lat,"lon_floor":lon,
             "original_EEA_2022_qualified_WWTP_count":int(c["public_EEA_both_measured_positive_unique_FW_outfalls"]),
             "original_ESAWC_2021_3deg_source_tile":tile,
             "original_ESAWC_2021_COG_public_URL":srcurl,
             "orig_WorldCover_categorical_sample_main":main,
             "orig_WorldCover_categorical_sample_secondary_256":sens,
             "original_512_vs256_absolute_class_fraction_differences":deviations,
             "original_landcountry_geometry_and_NODATA_gate":bool(good),
             "WorldCover_resolution_reliability_gate":bool(stable),
             "no_2022_butterfly_photo_or_species_information_accessed":True
            }
            rows.append(row)
            print("GESUI_DE_AT_ESAWC2021_2022_SOURCE_REGION_INDEPENDENT_LANDCLASS",json.dumps({
               "progress":i+1,"expected":args.max_cells,"country":country,"lat":lat,"lon":lon,
               "sourceTile":tile,
               "country_land_fraction":main["country_terrestrial_land_fraction_cell"],
               "source_WC2021_built_frac":main["built_up_share_of_country_terrestrial"],
               "source_WC2021_crop_frac":main["cropland_share_of_country_terrestrial"],
               "source_quality_pass":good,"stability":stable},sort_keys=True),flush=True)
    finally:
        for ds in cached.values():
            try:ds.close()
            except Exception:pass
    # Every original WWTP source cell is kept in source-only QC regardless masks.
    qa={
      "schema":"gesui_worldcover2021_original_categorical_COG_country_masked_EEA_DE_AT_1deg_source_only_v01",
      "timestamp":datetime.now(timezone.utc).isoformat(),
      "original_EEA_source_manifest_sha256":EA2022_SHA,
      "original_EEA_2022_source_tile_csv_SHA256":sha(args.source_tile_csv),
      "original_natural_earth_admin0_source":geoQC,
      "ESAWC_2021_original_public_COG_v200_S3_base":AWS_PREFIX,
      "original_WorldCover_product_DOI":"10.5281/zenodo.7254221",
      "original_WorldCover_categorical_raster_codes":LAND_COVER_CLASSES,
      "note_original_COG_10m_class_map_resampled_nearest_512_approx_not_10m_landcensus":True,
      "original_country_mask_only_from_NaturalEarth50m_polygon":True,
      "original_53_source_cells_total":len(cells),
      "source_only_processed_cells":len(rows),
      "source_only_requested_max_cells":args.max_cells,
      "unique_public_COG_3degree_tiles_accessed":sorted(cached),
      "source_land_quality_pass_count":sum(r["original_landcountry_geometry_and_NODATA_gate"] for r in rows),
      "source_resolution_stability_pass_count":sum(r["WorldCover_resolution_reliability_gate"] for r in rows),
      "landcover_resampling_main_side":RES_MAIN,
      "landcover_resampling_secondary_side":RES_SENS,
      "landcover_proxy_is_not_city_urbanization_or_exact_2022_landcover":True,
      "NO_2022_BUTTERFLY_PHOTO_RESULTS_OPENED":True,
      "NO_ACTUAL_P_OR_N_CHEMISTRY_AT_BUTTERFLY_SUBSTRATE":True,
      "source_cells":rows,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(qa,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_WORLD_COVER_2021_ORIGINAL_SOURCE_CATEGORICAL_LAND_QUALITY_GATE",
      json.dumps({k:qa[k] for k in ("original_53_source_cells_total","source_only_processed_cells","source_land_quality_pass_count","source_resolution_stability_pass_count")}),flush=True)

if __name__=="__main__":main()
