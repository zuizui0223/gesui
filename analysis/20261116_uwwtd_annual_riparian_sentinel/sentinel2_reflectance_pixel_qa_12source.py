#!/usr/bin/env python3
"""Radiometry, cloud, geometrical coverage only for 12 EEA measured-2022-NP
outfall system × genuine Sentinel-2 L2A 2018+2022 comparable local pixels.

THIS IS A SENSOR DATA QUALITY PILOT, NOT A WASTEWATER EFFECT OR BIOMASS
ESTIMATE. Original STAC source metadata frozen in successful run 43851791936.
Original 2018 L2A band scale/offset 0.0001/0 vs original 2022 0.0001/−0.1
must be applied correctly. True GBIF/organism outcomes are NOT accessed.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math
from contextlib import ExitStack
from pathlib import Path
import numpy as np
import rasterio
from rasterio.errors import RasterioIOError
from rasterio.windows import from_bounds,Window
from rasterio.enums import Resampling
from rasterio.warp import reproject
from pyproj import Transformer

SOURCE_SCHEMA="GESUI_20261116_EEA2022_MEASURED_NP_ONE_FW_OUTFALL_DIRECTED_HYDRORIVERS_SOURCE_GEO_V1"
STAC_SCHEMA="GESUI_20261116_EEA2022_TRUE_NP_SENTINEL2_LEGACY_L2A_2018_2022_SOURCE_SCENE_GATES_V1"
SCHEMA="GESUI_20261116_EEA_TRUE_NP_RIPARIAN_S2_12SOURCE_L2A_REFL_SCL_PIXEL_QA_V1"
OUTFALLS_PER_COUNTRY=6
MONTHS=(6,7,8)
YEARS=(2018,2022)
SITE_RADIUS_M=250
MIN_CLOUD_CLEAR_LAND_PIXELS=100
MIN_COMMON_YEARLY_MONTHS=2
MIN_TECHNICAL_SOURCE_TOTAL=8
MIN_COUNTRY_TECHNICAL_SOURCE=3
EXPECTED_URL_PREFIX="https://sentinel-cogs.s3.us-west-2.amazonaws.com/"

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def sources(directory):
    gate=json.loads((directory/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").read_text())
    fp=directory/"EEA2022_SATELLITE_SOURCE_BLIND_INDEPENDENT_OUTFALL_PILOT.csv"
    if (gate.get("schema")!=SOURCE_SCHEMA or not gate.get("original_source_pilot_pass")
        or gate.get("pilot_source_CSV_SHA256")!=sha(fp)):
        raise RuntimeError("Original published outfall local-source pilot provenance invalid")
    with fp.open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=48 or len({q["original_EEA2022_plant_id"] for q in rows})!=48:
        raise RuntimeError("Source genuine 48 country-balanced NP outfall pilot changed")
    samples=[]
    for cc in ("AT","DE"):
        group=[q for q in rows if q["EEA_original_countryCode"]==cc]
        if len(group)!=24:raise RuntimeError("Original source EEA 2022 pilot country denominator changed")
        group.sort(key=lambda q:(hashlib.sha256(q["original_EEA2022_plant_id"].encode()).hexdigest(),q["original_EEA2022_plant_id"]))
        samples.extend(group[:OUTFALLS_PER_COUNTRY])
    if len(samples)!=12:raise RuntimeError("Fixed source-smoke case cohort must have exactly 12 independent outlets")
    return samples,gate

def scenes(directory,src_gate):
    gate=json.loads((directory/"EEA2022_SENTINEL2_2018_2022_RIPARIAN_STAC_SCENE_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!=STAC_SCHEMA or not gate.get("source_STAC_gate_pass") or
        gate.get("origin_genuine_EEA2022_source_SHA")!=src_gate.get("source_EEA2022_original_manifest_Gitblob") or
        gate.get("source_geo_pilot_manifest_SHA")!=src_gate.get("pilot_source_CSV_SHA256") or
        gate.get("source_outfalls_with_both_years_SCENE_NOT_local_cloud_pass")!=48 or
        gate.get("original_STAC_collection")!="sentinel-2-l2a"):
        raise RuntimeError("Source original 2018+2022 SAME legacy COG collection L2A metadata not validated")
    s=json.loads((directory/"EEA2022_SOURCE_SCENE_RED_NIR_SCL_METADATA.json").read_text())
    if len(s)<48*2*3:raise RuntimeError("Original source STAC monthly scene set implausibly sparse")
    grouped=collections.defaultdict(list)
    for item in s:
        key=(item["original_EEA2022_plant_id"],int(item["year"]),int(item["calendar_month"]))
        grouped[key].append(item)
    for k in grouped:
        grouped[k].sort(key=lambda v:(v["original_tile_cloud_percent"],v["scene_datetime_UTC"],v["stac_scene_id"]))
    return grouped,gate

def radiometry(scene,band):
    asset=scene["asset_original_source"][band]
    href=asset["original_public_COG_href"]
    if not isinstance(href,str) or not href.startswith(EXPECTED_URL_PREFIX) or not href.endswith(".tif"):
        raise RuntimeError("STAC source public 2026 COG href does not match pinned AWS Sentinel COG domain")
    info=asset.get("original_raster_bands") or []
    if band=="scl":
        return href,None,None
    if len(info)!=1 or "scale" not in info[0]:
        raise RuntimeError("Missing actual source reflectance raster:bands scale for "+band)
    scale=float(info[0]["scale"])
    offset=float(info[0].get("offset",0.))
    if not (scale>0 and math.isfinite(scale) and math.isfinite(offset)):
        raise RuntimeError("Source numeric band radiometric mapping invalid")
    return href,scale,offset

def scene_pair(scene,source):
    bands={}
    for k in ("red","nir","scl"):
        bands[k]=radiometry(scene,k)
    if (abs(bands["red"][1]-bands["nir"][1])>1e-10 or
        abs(bands["red"][2]-bands["nir"][2])>1e-10):
        raise RuntimeError("Source original red vs NIR reflectance calibration differs")
    result={}
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
              CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
              GDAL_HTTP_MAX_RETRY="2",GDAL_HTTP_RETRY_DELAY="1",
              GDAL_CACHEMAX=64,VSI_CACHE="TRUE",VSI_CACHE_SIZE=1000000):
      with ExitStack() as st:
        # Original provider STAC references COG asset; no external download or source rehost.
        a={k:st.enter_context(rasterio.open("/vsicurl/"+bands[k][0])) for k in bands}
        rd=a["red"];nr=a["nir"];scl=a["scl"]
        if (rd.crs is None or nr.crs!=rd.crs or scl.crs!=rd.crs or
            rd.count!=1 or nr.count!=1 or scl.count!=1):
            raise RuntimeError("Original S2 GeoTIFF CRS/band count mismatch")
        if (rd.transform!=nr.transform or rd.width!=nr.width or rd.height!=nr.height):
            raise RuntimeError("S2 original red NIR original resolutions not congruent")
        if not (8<=abs(rd.transform.a)<=12 and 8<=abs(rd.transform.e)<=12):
            raise RuntimeError("S2 original RED source pixel resolution not ~10m")
        converter=Transformer.from_crs("EPSG:4326",rd.crs,always_xy=True)
        for role in ("UPSTREAM","DOWNSTREAM"):
            lat=float(source[f"new_{role}_1500m_lat"])
            lon=float(source[f"new_{role}_1500m_lon"])
            cx,cy=converter.transform(lon,lat)
            win=from_bounds(cx-270,cy-270,cx+270,cy+270,rd.transform)
            col_off=int(math.floor(win.col_off));row_off=int(math.floor(win.row_off))
            width=int(math.ceil(win.width))+2;height=int(math.ceil(win.height))+2
            window=Window(col_off,row_off,width,height)
            red_raw=rd.read(1,window=window,boundless=True,fill_value=0)
            nir_raw=nr.read(1,window=window,boundless=True,fill_value=0)
            tr=rd.window_transform(window)
            scl_re=np.zeros(red_raw.shape,dtype=np.uint8)
            reproject(source=rasterio.band(scl,1),destination=scl_re,
                dst_transform=tr,dst_crs=rd.crs,
                resampling=Resampling.nearest,src_nodata=0,dst_nodata=0)
            # 250m ACTUAL metric disc inside original source projected field location.
            rr,cc=np.indices(red_raw.shape)
            px=tr.c+(cc+.5)*tr.a+(rr+.5)*tr.b
            py=tr.f+(cc+.5)*tr.d+(rr+.5)*tr.e
            inside=((px-cx)**2+(py-cy)**2)<=SITE_RADIUS_M**2
            land=inside & np.isin(scl_re,[4,5])
            rs,ro=bands["red"][1:]
            ns,no=bands["nir"][1:]
            red=red_raw.astype(np.float32)*rs+ro
            nir=nir_raw.astype(np.float32)*ns+no
            physical=(red>0)&(nir>0)&(red<1.0)&(nir<1.0)&((red+nir)>.01)
            good=land & physical
            npix=int(good.sum());total=int(inside.sum())
            if npix:
                vi=(nir[good]-red[good])/(nir[good]+red[good])
                vi=vi[np.isfinite(vi)&(vi>=-1)&(vi<=1)]
            else:vi=np.array([],dtype=float)
            result[role]={"source_circle_native_10m_total_pixels":total,
                 "scene_pixel_local_SCL_clear_land_count":npix,
                 "source_pixel_clear_land_fraction":npix/max(1,total),
                 "source_valid_NDVI_count":len(vi),
                 "source_land_clear_NDVI_median":float(np.median(vi)) if len(vi) else None,
                 "source_SCL_land4_pixels":int((inside&(scl_re==4)).sum()),
                 "source_SCL_bare5_pixels":int((inside&(scl_re==5)).sum()),
                 "source_SCL_cloud8_9_10_pixels":int((inside&np.isin(scl_re,[8,9,10])).sum()),
                 "source_SCL_water6_pixels":int((inside&(scl_re==6)).sum()),
                 "source_redNIR_reflectance_offset_applied":ro,
                 "source_pixel_technical_gate_pass":len(vi)>=MIN_CLOUD_CLEAR_LAND_PIXELS}
    return result

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    chosen,sourcegate=sources(a.geography)
    group,scenegate=scenes(a.stac,sourcegate)
    quality=[];stats=collections.Counter();report=[]
    for study in chosen:
        fac=study["original_EEA2022_plant_id"]
        yearmonths={str(y):{} for y in YEARS}
        for year in YEARS:
            for month in MONTHS:
                scenechoices=group.get((fac,year,month),[])[:2]
                attempts=[]
                for item in scenechoices:
                    try:
                        ret=scene_pair(item,study)
                        passed=all(ret[role]["source_pixel_technical_gate_pass"] for role in ("UPSTREAM","DOWNSTREAM"))
                        info={"original_stac_scene_id":item["stac_scene_id"],
                          "original_scene_datetime_UTC":item["scene_datetime_UTC"],
                          "original_tile_cloud_pct":item["original_tile_cloud_percent"],
                          "source_year":year,"source_month":month,
                          "source_EEA_outfall_id":fac,"source_EEA_country":study["EEA_original_countryCode"],
                          "source_reflectance_rednir_asset_offset":radiometry(item,"red")[2],
                          "source_two_local_circles_pixel_qa":ret,
                          "local_pixel_qa_pass":bool(passed)}
                        attempts.append(info)
                        if passed:break
                    except (RasterioIOError,OSError,ValueError,RuntimeError) as ex:
                        attempts.append({"original_stac_scene_id":item["stac_scene_id"],
                            "source_COG_read_error_type":type(ex).__name__,
                            "source_COG_read_error_message":str(ex)[:250],
                            "local_pixel_qa_pass":False})
                good=next((v for v in attempts if v["local_pixel_qa_pass"]),None)
                if good is not None:
                    yearmonths[str(year)][str(month)]={"technical_pass":True,
                             "original_stac_scene_id":good["original_stac_scene_id"],
                             "original_scene_datetime_UTC":good["original_scene_datetime_UTC"],
                             "local_upstream_clear_land_NDVI_median":good["source_two_local_circles_pixel_qa"]["UPSTREAM"]["source_land_clear_NDVI_median"],
                             "local_downstream_clear_land_NDVI_median":good["source_two_local_circles_pixel_qa"]["DOWNSTREAM"]["source_land_clear_NDVI_median"],
                             "source_technical_scene_payload":good}
                else:
                    yearmonths[str(year)][str(month)]={"technical_pass":False,
                             "original_scene_count_tried":len(attempts),
                             "attempted_scenes":attempts}
                quality.extend(attempts)
                stats[(year,month,"passed" if good else "failed")]+=1
        common=[m for m in MONTHS if yearmonths["2018"][str(m)]["technical_pass"] and yearmonths["2022"][str(m)]["technical_pass"]]
        source_pass=len(common)>=MIN_COMMON_YEARLY_MONTHS
        report.append({"original_EEA2022_plant_id":fac,
            "original_EEA_country":study["EEA_original_countryCode"],
            "original_2022_source_reported_positive_measured_NP":True,
            "source_EEA_original_annual_N_t_year":float(study["source_N_t_year_reported_measured_2022"]),
            "source_EEA_original_annual_P_t_year":float(study["source_P_t_year_reported_measured_2022"]),
            "source_valid_pair_year_common_June_July_Aug_months":common,
            "source_minimum_two_common_months":source_pass,
            "source_monthwise_technical_metrics":yearmonths,
            "NDVI_is_sensor_landgreen_not_wastewater_causal_or_vegetation_biomass":True})
        print("GESUI_NDVI_SOURCE_RASTER_SMOKE_COMPLETE="+fac+" year_common="+str(len(common)),flush=True)
    eligible=[x for x in report if x["source_minimum_two_common_months"]]
    cc=collections.Counter(x["original_EEA_country"] for x in eligible)
    supported=(len(eligible)>=MIN_TECHNICAL_SOURCE_TOTAL and all(cc.get(k,0)>=MIN_COUNTRY_TECHNICAL_SOURCE for k in ("AT","DE")))
    out={"schema":SCHEMA,
       "status":"PASS_SENTINEL_COG_12SOURCE_LOCAL_PIXEL_QUALITY_FOR_NEXT_MASKED_STUDY" if supported else
                "STOP_SENTINEL_COG_12SOURCE_LOCAL_PIXEL_QUALITY",
       "original_source_true_2022_annual_measured_NP_outfalls":1156,
       "original_source_geom_source_pilot":48,
       "original_same_collection_S2_STAC_2018_2022_source_gate_pass":True,
       "source_selected_original_outfall_cases":len(chosen),
       "source_country_6_AT_6_DE_responseblind":True,
       "original_same_calendar_months_tested":list(MONTHS),
       "original_geometrical_site_radius_m":SITE_RADIUS_M,
       "source_mask_SCL_accept_codes":[4,5],
       "source_mask_SCL_excludes_cloud_shadow_water_snow":True,
       "source_clear_land_minimum_pixels_per_scene_per_original_bank":MIN_CLOUD_CLEAR_LAND_PIXELS,
       "source_annual_2018_vs2022_raster_original_offsets_honoured":True,
       "source_technical_outfalls_with_two_calendar_months_both_years":len(eligible),
       "source_technical_pass_by_country":dict(cc),
       "source_technical_acceptance_gate_pass":bool(supported),
       "original_S2_COG_pixel_attempt_type_counts":dict(collections.Counter(
          (v.get("source_COG_read_error_type") or ("passed" if v.get("local_pixel_qa_pass") else "quality_failed"))
           for v in quality)),
       "source_monthly_qa_counts":{f"{k[0]}_{k[1]}_{k[2]}":v for k,v in sorted(stats.items())},
       "cases":report,
       "source_ESA_WorldCover_mask_and_real_reach_controls_not_applied":True,
       "satellite_L2A_NDVI_not_plant_biomass_or_plant_species_richness":True,
       "actual_same_source_freshwater_TNTP_not_measured_at_bank":True,
       "no_direct_EQR_or_BMS_butterfly_source_overlap_assumed":True,
       "no_wastewater_or_aquatic_terrestrial_treatment_effect_fitted":True}
    (a.output/"EEA2022_12SOURCE_L2A_LOCAL_PIXEL_RADIOMETRIC_SCL_TECHNICAL_QA.json").write_text(
      json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2022_S2_L2A_PIXEL_QA="+json.dumps({
       "status":out["status"],
       "source_12_case_pilot":len(chosen),
       "source_case_year_pixel_pass":len(eligible),
       "by_country":dict(cc),
       "error_types":out["original_S2_COG_pixel_attempt_type_counts"],
       "source_raster_reflectance_offsets_applied":True,
       "no_ecological_effect_inferred":True
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geography",required=True,type=Path)
    p.add_argument("--stac",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
