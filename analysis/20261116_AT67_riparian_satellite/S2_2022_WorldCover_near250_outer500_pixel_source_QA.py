#!/usr/bin/env python3
"""2022 real Austrian river-field-site near 250m and outer 250-500m S2 optical
greenness technical quality, class-matched 2021 WorldCover stable landcover.

Source-only SATELLITE pixel availability test, not an ecology/sewage model.
Only genuine original 2022 TP+benthic EQR STATION IDs chosen pre-pixel; DO NOT
use their TP, EQR, WWTP P to select scenes or sites. Optical NDVI != biomass.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,sys
from contextlib import ExitStack
from pathlib import Path
import numpy as np
import rasterio
from rasterio.errors import RasterioIOError
from rasterio.windows import Window,from_bounds
from rasterio.warp import reproject
from rasterio.enums import Resampling
from pyproj import Transformer

SCHEMA="GESUI_20261116_AT67_RIVER_FIELD_S2_2022_WORLD_COVER_CLASSMATCHED_PIXEL_SOURCE_V1"
STAC_SCHEMA="GESUI_20261116_AT67_REAL_MONITOR_2022_S2_B04_B08_SCL_SCENE_METADATA_V1"
ORIG_BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
REGIONS={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
S2_HOST="https://sentinel-cogs.s3.us-west-2.amazonaws.com/"
WC_HOST="https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
NATURAL={10,20,30,90}
SOURCE_MONTHS=(6,7,8)
MIN_SCL_INNER=100
MIN_SCL_OUTER=200
MIN_NATURAL_INNER=40
MIN_NATURAL_OUTER=100
CLASS_COMMON_MIN=10
INNER_M=250.
OUTER_M=500.

def sha(b):return hashlib.sha256(b).hexdigest()

def load_source(directory):
    file=directory/"AUSTRIA67_EEA2022_SAME_STATION_S2_2022_SCENE_SOURCE_GATE.json"
    original=json.loads(file.read_text())
    if (original.get("schema")!=STAC_SCHEMA or
        not original.get("source_gate_pass") or
        original.get("source_scene_supported_sites_by_geographic_region")!=REGIONS or
        original.get("original_benthic_site_source_gitblob")!=ORIG_BQE_BLOB):
        raise RuntimeError("Original 67 Austrian 2022 same-actual-water-station source gate failed")
    rows=original.get("source_site_months",[])
    if len(rows)!=67 or len({q["original_BMLUK_PKINT_2022_true_water_benthos"] for q in rows})!=67:
        raise RuntimeError("Original 67 authenticated source station coordinates changed")
    images=json.loads((directory/"AUSTRIA67_S2_2022_JUN_AUG_SELECTED_SCENE_RED_NIR_SCL.json").read_text())
    if len(images)!=394 or len({x["original_national_PKINT"] for x in images})!=67:
        raise RuntimeError("Original 394 original 2022 red/nir/scl scene item source set changed")
    grouped=collections.defaultdict(list)
    for s in images:
        if int(s["month"]) not in SOURCE_MONTHS:
            raise RuntimeError("Source scene month is not original June-August")
        grouped[(s["original_national_PKINT"],int(s["month"]))].append(s)
    for k,lis in grouped.items():
        if len(lis)>2:
            raise RuntimeError("Original STAC source preselected >2 images/month/site")
        lis.sort(key=lambda x:(x["tile_cloud_fraction_pct"],x["original_scene_UTC"],x["original_stac_scene_id"]))
    return original,rows,grouped

def select(rows,subset):
    if subset=="all":
        return sorted(rows,key=lambda r:r["original_BMLUK_PKINT_2022_true_water_benthos"])
    groups=collections.defaultdict(list)
    for row in rows:groups[row["fixed_original_benthic_geo_block"]].append(row)
    chosen=[]
    for b in sorted(REGIONS):
        q=sorted(groups[b],key=lambda z:(
            sha(z["original_BMLUK_PKINT_2022_true_water_benthos"].encode()),
            z["original_BMLUK_PKINT_2022_true_water_benthos"]))
        if len(q)<2:raise RuntimeError("Not enough original stations for deterministic smoke subset")
        chosen.extend(q[:2])
    if len(chosen)!=10:raise RuntimeError("Expected 2 original sites per 5 true regions")
    return chosen

def worldcover_url(lat,lon):
    la=int(math.floor(lat/3)*3)
    lo=int(math.floor(lon/3)*3)
    tile=("N" if la>=0 else "S")+f"{abs(la):02d}"+("E" if lo>=0 else "W")+f"{abs(lo):03d}"
    return WC_HOST+f"ESA_WorldCover_10m_2021_v200_{tile}_Map.tif",tile

def asset(item,name):
    meta=item["original_bands"][name]
    url=meta["href"]
    if not isinstance(url,str) or not url.startswith(S2_HOST) or not url.endswith(".tif"):
        raise RuntimeError("Original 2022 Sentinel asset is not pinned original COG source")
    if name=="scl":return url,None,None
    b=meta.get("source_raster_band_schema") or []
    if len(b)!=1 or "scale" not in b[0]:
        raise RuntimeError("Original 2022 red/NIR radiometric scale and offset absent")
    sc=float(b[0]["scale"]);off=float(b[0].get("offset",0))
    if not (0<sc<=1 and -1<=off<=1):
        raise RuntimeError("Original 2022 reflectance band source scale or offset invalid")
    return url,sc,off

def pixel_scene(scene,station):
    vals={k:asset(scene,k) for k in ("red","nir","scl")}
    if (abs(vals["red"][1]-vals["nir"][1])>1e-12 or
        abs(vals["red"][2]-vals["nir"][2])>1e-12):
        raise RuntimeError("Original 2022 red and nir asset radiometry differs unexpectedly")
    lat=float(station["original_EEA2022_lat"])
    lon=float(station["original_EEA2022_lon"])
    WC,tile=worldcover_url(lat,lon)
    with rasterio.Env(
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
        GDAL_HTTP_MAX_RETRY="2",
        GDAL_HTTP_RETRY_DELAY="1",
        GDAL_CACHEMAX=64,
        VSI_CACHE="TRUE",VSI_CACHE_SIZE=1000000):
      with ExitStack() as st:
        source={q:st.enter_context(rasterio.open("/vsicurl/"+vals[q][0]))
                for q in ("red","nir","scl")}
        world=st.enter_context(rasterio.open("/vsicurl/"+WC))
        red=source["red"];nir=source["nir"];scl=source["scl"]
        if (red.crs is None or nir.crs!=red.crs or
            scl.crs!=red.crs or world.crs is None or
            red.transform!=nir.transform or red.width!=nir.width or red.height!=nir.height):
            raise RuntimeError("True 2022 S2 original red/nir SCL source coordinate/resolution mismatch")
        if not (8<=abs(red.transform.a)<=12 and 8<=abs(red.transform.e)<=12):
            raise RuntimeError("Original 2022 band red native cell is not ~10m")
        converter=Transformer.from_crs("EPSG:4326",red.crs,always_xy=True)
        cx,cy=converter.transform(lon,lat)
        window=from_bounds(cx-520,cy-520,cx+520,cy+520,red.transform)
        window=Window(int(math.floor(window.col_off)),int(math.floor(window.row_off)),
                      int(math.ceil(window.width))+2,int(math.ceil(window.height))+2)
        red_raw=red.read(1,window=window,boundless=True,fill_value=0)
        nir_raw=nir.read(1,window=window,boundless=True,fill_value=0)
        target_tr=red.window_transform(window)
        sc=np.zeros(red_raw.shape,dtype=np.uint8)
        wc=np.zeros(red_raw.shape,dtype=np.uint8)
        reproject(source=rasterio.band(scl,1),destination=sc,dst_transform=target_tr,
                  dst_crs=red.crs,resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
        reproject(source=rasterio.band(world,1),destination=wc,dst_transform=target_tr,
                  dst_crs=red.crs,resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
        rr,cc=np.indices(red_raw.shape)
        xx=target_tr.c+(cc+.5)*target_tr.a+(rr+.5)*target_tr.b
        yy=target_tr.f+(cc+.5)*target_tr.d+(rr+.5)*target_tr.e
        dist2=(xx-cx)**2+(yy-cy)**2
        inner=dist2<=INNER_M**2
        outer=(dist2>INNER_M**2)&(dist2<=OUTER_M**2)
        reflect_r=red_raw.astype(np.float32)*vals["red"][1]+vals["red"][2]
        reflect_n=nir_raw.astype(np.float32)*vals["nir"][1]+vals["nir"][2]
        phys=(reflect_r>.005)&(reflect_n>.005)&(reflect_r<1)&(reflect_n<1)
        clear=np.isin(sc,[4,5])&phys
        natural=np.isin(wc,list(NATURAL))&clear
        counts={
          "inner10m_total_pix":int(inner.sum()),"outer10m_total_pix":int(outer.sum()),
          "inner_SCL4_5_clearland_pix":int((inner&clear).sum()),
          "outer_SCL4_5_clearland_pix":int((outer&clear).sum()),
          "inner_2021_worldcover_natural_clear_pix":int((inner&natural).sum()),
          "outer_2021_worldcover_natural_clear_pix":int((outer&natural).sum()),
          "worldcover_2021_classes_inner":{
             str(k):int((inner&(wc==k)).sum()) for k in (10,20,30,40,50,80,90)},
          "worldcover_2021_tile_ID":tile,
          "original_2022_redNIR_reflectance_offset_applied":vals["red"][2]}
        enough=(
          counts["inner_SCL4_5_clearland_pix"]>=MIN_SCL_INNER and
          counts["outer_SCL4_5_clearland_pix"]>=MIN_SCL_OUTER and
          counts["inner_2021_worldcover_natural_clear_pix"]>=MIN_NATURAL_INNER and
          counts["outer_2021_worldcover_natural_clear_pix"]>=MIN_NATURAL_OUTER)
        if not enough:
            return {"pixel_qa_pass":False,"reason":"insufficient_clear_land_or_2021_natural_pixels",
                    **counts}
        vi=(reflect_n-reflect_r)/(reflect_n+reflect_r+1e-12)
        vi=np.where(np.isfinite(vi)&(vi>=-1)&(vi<=1),vi,np.nan)
        class_items=[]
        for k in sorted(NATURAL):
            mask=(wc==k)&natural&np.isfinite(vi)
            n=int((inner&mask).sum())
            m=int((outer&mask).sum())
            if n<CLASS_COMMON_MIN or m<CLASS_COMMON_MIN:continue
            class_items.append({
                "source_2021_WorldCover_class":k,
                "source_original_clear_inner_pixels":n,
                "source_original_clear_outer_pixels":m,
                "source_inner_2022_S2_NDVI_median":float(np.median(vi[inner&mask])),
                "source_outer_2022_S2_NDVI_median":float(np.median(vi[outer&mask]))})
        if not class_items:
            return {"pixel_qa_pass":False,"reason":"no_common_natural_landcover_class_in_both_rings",
                    **counts}
        weight=np.asarray([x["source_original_clear_inner_pixels"] for x in class_items],float)
        weight=weight/weight.sum()
        classdiff=np.asarray([
          x["source_inner_2022_S2_NDVI_median"]-x["source_outer_2022_S2_NDVI_median"] for x in class_items])
        paired=float(np.dot(weight,classdiff))
        return {"pixel_qa_pass":True,"reason":"ok",
            "source_2022_same_scene_classmatched_near_minus_outer_NDVI":paired,
            "source_common_natural_cover_classes":class_items,
            "source_2022_inner_NDVI_all_natural_land":float(np.nanmedian(vi[inner&natural])),
            "source_2022_outer_NDVI_all_natural_land":float(np.nanmedian(vi[outer&natural])),
            **counts}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    source,sites,source_scenes=load_source(a.stac)
    chosen=select(sites,a.subset)
    results=[];trials=[]
    for i,site in enumerate(chosen,1):
        pk=site["original_BMLUK_PKINT_2022_true_water_benthos"]
        monthvals={}
        for month in SOURCE_MONTHS:
            cand=source_scenes.get((pk,month),[])
            attempts=[];accepted=None
            for item in cand:
                try:
                    r=pixel_scene(item,site)
                except (RasterioIOError,RuntimeError,OSError,ValueError) as ex:
                    r={"pixel_qa_pass":False,"reason":"source_COG_pixel_read_error",
                       "error_type":type(ex).__name__,"message":str(ex)[:240]}
                record={"original_PKINT":pk,"source_region":site["fixed_original_benthic_geo_block"],
                   "month":month,"original_scene_id":item["original_stac_scene_id"],
                   "original_scene_utc":item["original_scene_UTC"],
                   "source_tile_cloud_pct":item["tile_cloud_fraction_pct"],**r}
                attempts.append(record)
                if r.get("pixel_qa_pass"):
                    accepted=record;break
            trials.extend(attempts)
            monthvals[str(month)]={"technical_pixel_pass":bool(accepted),
               "scene_attempts":len(attempts),
               "adopted":accepted if accepted else None}
        good=[v["adopted"] for v in monthvals.values() if v["technical_pixel_pass"]]
        accepted_year=len(good)>=2
        paired=float(np.median([
          v["source_2022_same_scene_classmatched_near_minus_outer_NDVI"] for v in good])) if good else None
        results.append({"original_2022_genuine_aquatic_station_PKINT":pk,
           "original_source_geo_region":site["fixed_original_benthic_geo_block"],
           "source_one_year_2022_valid_months":sorted(int(k) for k,v in monthvals.items() if v["technical_pixel_pass"]),
           "source_site_has_2_months_natural_landcover_and_clear_pixel":accepted_year,
           "site_2022_l2a_s2_near_minus_outer_land_classmatched_NDVI_descriptive_only":paired,
           "source_monthwise_pixel_qa":monthvals,
           "not_original_local_wastewater_P_uptake_or_plant_biomass":True})
        print("GESUI_AT67_2022_S2_WC_"+a.subset.upper()+"_SOURCE_SITE="+str(i)+"/"+str(len(chosen))+
              " "+pk+" validmonths="+str(len(good)),flush=True)
    good=[q for q in results if q["source_site_has_2_months_natural_landcover_and_clear_pixel"]]
    groups=collections.Counter(q["original_source_geo_region"] for q in good)
    if a.subset=="smoke":
        eligible=len(good)>=6 and sum(v>0 for v in groups.values())>=4
    else:
        eligible=len(good)>=40 and sum(v>=5 for v in groups.values())>=4
    result={"schema":SCHEMA,
      "status":"PASS_PIXEL_TECHNICAL_SOURCE_GATE_NO_ECOLOGY_FIT" if eligible else
               "STOP_SOURCE_NATURAL_LAND_RIPARIAN_LOCAL_PIXEL_COVERAGE",
      "fixed_subset":"original_hashed_two_PKINT_per_geographical_region" if a.subset=="smoke" else "all_original67",
      "original_same_site_2022_STAC_metadata_SOURCE_GATE":{"schema":STAC_SCHEMA,
          "original_source_sites":67,"original_scene_preselected":394},
      "original_source_2022_station_count_tested":len(chosen),
      "source_precommitted_pixel_months":list(SOURCE_MONTHS),
      "source_near_river_circle_radius_m":INNER_M,
      "source_outer_same_station_annulus_m":[INNER_M,OUTER_M],
      "source_ESA_WorldCover2021_v200_natural_classes":sorted(NATURAL),
      "source_ESA_WorldCover2021_v200_cropland_builtwater_excluded":True,
      "source_min_SCL_clear_pixel_near_outer":[MIN_SCL_INNER,MIN_SCL_OUTER],
      "source_min_WorldCover_natural_clear_pixel_near_outer":[MIN_NATURAL_INNER,MIN_NATURAL_OUTER],
      "source_min_matched_worldcover_class_pixels_per_zone":CLASS_COMMON_MIN,
      "original_2022_n_stations_2_month_natural_land_pixel_QA":len(good),
      "original_2022_geographic_blocks_with_valid_pixels":dict(groups),
      "predeclared_source_pixel_gate_pass":bool(eligible),
      "source_site_pixel_qa":results,
      "source_remote_COG_error_counts":dict(collections.Counter(
          z.get("reason") for z in trials if z.get("reason")!="ok")),
      "real_outfall_load_and_water_TP_or_EQR_not_opened_in_this_sensor_QA":True,
      "worldcover_2021_landcover_not_verified_real_riverbank_biological_vegetation":True,
      "original_source_2022_NDVI_is_not_biomass_species_richness_or_EcosystemService":True,
      "no_aquatic_terrestrial_causal_or_butterfly_model_fitted":True}
    out=a.output/("AUSTRIA2022_"+a.subset.upper()+"_S2_WORLDCOVER_RIPARIAN_PIXELS_SOURCE_QA.json")
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_AT67_2022_S2_WORLDCOVER_"+a.subset.upper()+"_SOURCE_GATE="+json.dumps({
       "status":result["status"],"source_sites_tested":len(chosen),
       "source_valid_n":len(good),"regions":dict(groups),
       "pixel_source_failures":result["source_remote_COG_error_counts"],
       "no_vegetation_vs_EQR_or_WWTP_effect_fitted":True
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--stac",type=Path,required=True)
    p.add_argument("--subset",choices=("smoke","all"),required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
