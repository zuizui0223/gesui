#!/usr/bin/env python3
"""Source-only original Sentinel-2 2018/2022 scenes at actually measured 2022
EEA-UWWTD nitrogen/phosphorus freshwater outfalls, true upstream/downstream
directed 1.5km local sites. Original pixels NDVI/GBIF ecology NOT opened.

Earth Search STAC v1, public Sentinel-2 C1 Level2A 10m surface reflectance.
Actual red / nir / scl COGs checked in genuine scene metadata.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,time
from pathlib import Path
import requests
from shapely.geometry import shape,box

STAC_URL="https://earth-search.aws.element84.com/v1/search"
STAC_COLLECTION="sentinel-2-c1-l2a"
SOURCE_SCHEMA="GESUI_20261116_EEA2022_MEASURED_NP_ONE_FW_OUTFALL_DIRECTED_HYDRORIVERS_SOURCE_GEO_V1"
OUT_SCHEMA="GESUI_20261116_EEA2022_TRUE_NP_SENTINEL2_L2A_2018_2022_SOURCE_SCENE_GATES_V1"
YEARS=(2018,2022)
MONTHS=(5,6,7,8,9)
CLOUD_LIMIT=40
MIN_SEASON_SCENE_DATES=3
MIN_SUPPORTED_OUTFALLS=20
MIN_BY_COUNTRY=8

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_source(folder):
    q=json.loads((folder/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").read_text())
    p=folder/"EEA2022_SATELLITE_SOURCE_BLIND_INDEPENDENT_OUTFALL_PILOT.csv"
    if (q.get("schema")!=SOURCE_SCHEMA or
        q.get("source_status")!="PASS_SOURCE_GEOMETRY_FOR_STAC_METADATA_ONLY" or
        not q.get("original_source_pilot_pass") or
        not p.exists() or sha(p)!=q.get("pilot_source_CSV_SHA256")):
        raise RuntimeError("Original genuine reported facility/riparian source pilot not authenticated")
    with p.open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=48 or collections.Counter(r["EEA_original_countryCode"] for r in rows)!={"AT":24,"DE":24}:
        raise RuntimeError("Exactly 48 response-blind original spatial source outfalls not preserved")
    return rows,q

def original_scene_search(session,bbox,year):
    search={"collections":[STAC_COLLECTION],
      "bbox":bbox,
      "datetime":f"{year}-05-01T00:00:00Z/{year}-09-30T23:59:59Z",
      "query":{"eo:cloud_cover":{"lte":CLOUD_LIMIT}},
      "limit":100}
    seen=set();features=[]
    url=STAC_URL;method="POST";payload=search
    for page in range(7):
        error=None
        for retry in range(4):
            try:
                if method=="POST":
                    r=session.post(url,json=payload,timeout=(15,90))
                else:
                    r=session.get(url,timeout=(15,90))
                if r.status_code in (429,500,502,503,504):
                    raise requests.RequestException("STAC retryable "+str(r.status_code))
                r.raise_for_status()
                obj=r.json()
                if not isinstance(obj,dict) or "features" not in obj:
                    raise RuntimeError("STAC does not return original scene features")
                break
            except (requests.RequestException,ValueError,RuntimeError) as ex:
                error=ex
                if retry==3:raise RuntimeError("Official STAC source unavailable: "+str(error))
                time.sleep(3+4*retry)
        for item in obj["features"]:
            if item.get("id") not in seen:
                features.append(item)
                seen.add(item.get("id"))
        more=next((t for t in obj.get("links",[]) if t.get("rel")=="next"),None)
        if not more:break
        if len(features)>550:
            raise RuntimeError("More than 550 original scene items for exact source outfall-year, need stronger metadata paging")
        url=more.get("href","")
        if not url.startswith("https://earth-search.aws.element84.com/"):
            raise RuntimeError("STAC pagination off expected source domain")
        method=more.get("method","GET").upper()
        payload=more.get("body",None)
        if method not in ("POST","GET"):
            raise RuntimeError("Unsupported STAC source pagination method")
    else:
        raise RuntimeError("STAC scene original source pages not exhausted")
    return features

def local_geo_boxes(row):
    rectangles=[]
    for role in ("UPSTREAM","DOWNSTREAM"):
        lat=float(row[f"new_{role}_1500m_lat"])
        lon=float(row[f"new_{role}_1500m_lon"])
        if not (45<=lat<=56 and 5<=lon<=18):
            raise RuntimeError("Original DE AT source official outlet coordinates outside reasonable study envelope")
        half_lat=250/111100
        half_lon=250/(111100*math.cos(math.radians(lat)))
        rectangles.append(box(lon-half_lon,lat-half_lat,lon+half_lon,lat+half_lat))
    x0=min(x.bounds[0] for x in rectangles);y0=min(x.bounds[1] for x in rectangles)
    x1=max(x.bounds[2] for x in rectangles);y1=max(x.bounds[3] for x in rectangles)
    return rectangles,[x0,y0,x1,y1]

def good_scene(item,regions,year):
    if item.get("collection")!=STAC_COLLECTION:
        return None
    prop=item.get("properties") or {}
    dt=str(prop.get("datetime") or "")
    if not dt.startswith(str(year)):
        return None
    try:
        month=int(dt[5:7]);cloud=float(prop["eo:cloud_cover"])
    except (TypeError,ValueError,KeyError):
        return None
    if month not in MONTHS or not math.isfinite(cloud) or not(0<=cloud<=CLOUD_LIMIT):
        return None
    if item.get("geometry") is None:
        return None
    scene=shape(item["geometry"])
    if not all(scene.covers(x) for x in regions):
        return None
    assets=item.get("assets") or {}
    if not all(key in assets for key in ("red","nir","scl")):
        return None
    refs={}
    for name in ("red","nir","scl"):
        a=assets[name]
        href=a.get("href","")
        if not (isinstance(href,str) and href.startswith("https://")):
            return None
        bandlist=a.get("raster:bands") or []
        ref={"original_public_COG_href":href,
             "type":a.get("type"),
             "original_raster_bands":bandlist[:1]}
        refs[name]=ref
    return {"stac_scene_id":item["id"],"scene_datetime_UTC":dt,
            "calendar_month":month,"original_tile_cloud_percent":cloud,
            "asset_original_source":refs,
            "original_scene_footprint_contains_both_250m_radius_site_boxes":True}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    source,gate=read_source(a.source)
    session=requests.Session()
    session.headers.update({"User-Agent":"gesui-source-original-Sentinel2-EEA2022/20261116","Accept":"application/json"})
    scenes=[];groups=[];not_qualified=[]
    for j,row in enumerate(sorted(source,key=lambda x:x["original_EEA2022_plant_id"]),1):
        rect,bbox=local_geo_boxes(row)
        years={}
        for year in YEARS:
            raw=original_scene_search(session,bbox,year)
            relevant=[]
            for item in raw:
                x=good_scene(item,rect,year)
                if x is not None:relevant.append(x)
            unique={}
            for x in relevant:unique[x["stac_scene_id"]]=x
            good=list(unique.values())
            count_days=len({x["scene_datetime_UTC"][:10] for x in good})
            available_months=sorted({x["calendar_month"] for x in good})
            monthly=[]
            for m in MONTHS:
                options=sorted((x for x in good if x["calendar_month"]==m),
                        key=lambda v:(v["original_tile_cloud_percent"],v["scene_datetime_UTC"],v["stac_scene_id"]))
                monthly.extend(options[:2])
            yrs={"source_original_STAC_scene_features_returned":len(raw),
                 "source_cloud_geo_and_COG_asset_valid_scenes":len(good),
                 "source_valid_independent_acquisition_days":count_days,
                 "source_valid_months":available_months,
                 "original_source_selected_max_two_lowcloud_scenes_per_month":len(monthly),
                 "source_sufficient_scene_dates":count_days>=MIN_SEASON_SCENE_DATES}
            years[str(year)]=yrs
            for scene in monthly:
                scenes.append({"original_EEA2022_plant_id":row["original_EEA2022_plant_id"],
                  "original_country":row["EEA_original_countryCode"],
                  "year":year,**scene})
        both=years["2018"]["source_sufficient_scene_dates"] and years["2022"]["source_sufficient_scene_dates"]
        months_common=set(years["2018"]["source_valid_months"])&set(years["2022"]["source_valid_months"])
        group={"original_EEA2022_plant_id":row["original_EEA2022_plant_id"],
            "country":row["EEA_original_countryCode"],
            "source_report_2022_N_t_year":float(row["source_N_t_year_reported_measured_2022"]),
            "source_report_2022_P_t_year":float(row["source_P_t_year_reported_measured_2022"]),
            "original_surface_buffer_radius_m":250,
            "actual_river_1500m_each_side":True,
            "source_comparable_2018_2022_scenes_both_years":both,
            "source_common_valid_calendar_months":sorted(months_common),
            "source_2018":years["2018"],"source_2022":years["2022"],
            "original_Sentinel_pixels_not_read":True,
            "source_year2022_annual_NP_report_not_Sentinel_scene_level_discharge":True}
        groups.append(group)
        print("GESUI_ORIGINAL_S2_STAC_SOURCE_OUTFALL_COMPLETE="+row["original_EEA2022_plant_id"]+" "+str(j)+"/"+str(len(source)),flush=True)
        time.sleep(.16)
    if len(groups)!=48 or not scenes:
        raise RuntimeError("Source authentic Sentinel2 original scene feasibility list invalid")
    passed=[v for v in groups if v["source_comparable_2018_2022_scenes_both_years"]]
    counts=collections.Counter(x["country"] for x in passed)
    supported=len(passed)>=MIN_SUPPORTED_OUTFALLS and counts["AT"]>=MIN_BY_COUNTRY and counts["DE"]>=MIN_BY_COUNTRY
    report={"schema":OUT_SCHEMA,
        "status":"PASS_ORIGINAL_S2_STAC_SCENE_SOURCE_FOR_NDVI_PROTOCOL" if supported else "STOP_S2_SAME_YEAR_SCENE_SOURCE_OPPORTUNITY",
        "origin_genuine_EEA2022_source_SHA":gate["source_EEA2022_original_manifest_Gitblob"],
        "source_geo_pilot_manifest_SHA":gate["pilot_source_CSV_SHA256"],
        "original_source_pilot_n":len(source),
        "original_STAC_endpoint":STAC_URL,
        "original_STAC_collection":STAC_COLLECTION,
        "source_cloud_threshold_pct":CLOUD_LIMIT,
        "source_query_growing_season_years":YEARS,
        "source_query_growing_season_months":MONTHS,
        "minimum_source_original_scene_acquisition_days_per_year":MIN_SEASON_SCENE_DATES,
        "source_outfalls_with_both_years_SCENE_NOT_local_cloud_pass":len(passed),
        "source_passing_country_counts":dict(counts),
        "source_min_outfalls":MIN_SUPPORTED_OUTFALLS,
        "source_min_each_DE_AT":MIN_BY_COUNTRY,
        "source_monthly_original_scene_assets":len(scenes),
        "source_STAC_gate_pass":supported,
        "source_pilot_outfall_seasonal_eligibility":groups,
        "source_item_genuine_original_footprint_and_COG_red_nir_scl_checked":True,
        "actual_circles_local_scl_cloud_mask_NOT_YET_TESTED":True,
        "2018_and_2022_annual_NP_outfall_reporting_comparable":False,
        "actual_local_NDVI_landcover_or_plant_biomass_NOT_CALCULATED":True,
        "unmeasured_other_sources_and_reference_rivers":True,
        "no_causal_sewage_or_benthic_butterfly_effect":True}
    (a.output/"EEA2022_SENTINEL2_2018_2022_RIPARIAN_STAC_SCENE_SOURCE_GATE.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    (a.output/"EEA2022_SOURCE_SCENE_RED_NIR_SCL_METADATA.json").write_text(
       json.dumps(scenes,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2022_S2_STAC_METADATA_SOURCE_GATE="+json.dumps({
       "status":report["status"],"source_pilot_outfalls":len(source),
       "source_both_year_scenes":len(passed),"by_country":dict(counts),
       "seasonal_scene_assets":len(scenes),"source_gate_pass":supported,
       "NDVI_plant_EQR_butterfly_not_analyzed":True},sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
