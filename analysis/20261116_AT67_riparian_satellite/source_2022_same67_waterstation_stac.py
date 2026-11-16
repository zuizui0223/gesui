#!/usr/bin/env python3
"""SOURCE-ONLY Sentinel2 2022 scene opportunity around exactly 67 original
Austrian river field stations with same-2022 measured TP & benthic EQR and
source-model 2022 upstream reported WWTP P. Exclude one prior mismatched
BMLUK vs EEA geo point BEFORE satellite imagery sampling.

Does not open NDVI/reflectance pixel, ecology effects, actual leaf P uptake.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math,time
from pathlib import Path
import requests
from shapely.geometry import shape,box

ROOT=Path(__file__).resolve().parents[2]
EEA_BQE=ROOT/"results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
SOURCE="https://earth-search.aws.element84.com/v1/search"
COLLECTION="sentinel-2-l2a"
SCHEMA="GESUI_20261116_AT67_REAL_MONITOR_2022_S2_B04_B08_SCL_SCENE_METADATA_V1"
EXCLUDED="ATFW54110117"
MONTHS=(6,7,8)
REGIONS={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
MIN_SITES=40
MIN_GROUPS=4
MIN_PER_SUPPORTED_GROUP=5

def blob(p):
    s=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(s))+"\0").encode()+s).hexdigest()
def sha(data):return hashlib.sha256(data).hexdigest()

def original_sites():
    if blob(EEA_BQE)!=BQE_BLOB:
        raise RuntimeError("Original 2022 Austrian benthic 68 source file Git blob changed")
    with EEA_BQE.open(encoding="utf-8-sig",newline="") as f:allrows=list(csv.DictReader(f))
    if len(allrows)!=2475:raise RuntimeError("Source official EEA2022 2475 benthic sites changed")
    at=[q for q in allrows if q["countryCode"]=="AT"]
    if len(at)!=68 or len({q["site_key"] for q in at})!=68:
        raise RuntimeError("Original EEA observed station-year 68 identities no longer unique")
    x=[]
    for q in at:
        identifier=q["site_key"].split("|")
        if len(identifier)!=3 or identifier[0]!="AT" or identifier[2]!="euMonitoringSiteCode":
            raise RuntimeError("Original observed 2022 BQE PKINT identity changed")
        if identifier[1]==EXCLUDED:continue
        lat=float(q["lat"]);lon=float(q["lon"])
        if not(45<=lat<=51 and 8<=lon<=19):raise RuntimeError("Source Austrian observed station coordinate invalid")
        x.append({"PKINT":identifier[1],"lat":lat,"lon":lon,
           "block":q["region_3deg_block"],"BQE_method":q["procedureClassificationSystem"],
           "source_2022_EEA_benthic_observed_and_national_FW1_TP_site":True,
           "original_river_native_HYRIV_ID":q["HYRIV_ID"]})
    if len(x)!=67 or collections.Counter(z["block"] for z in x)!=REGIONS:
        raise RuntimeError("Original 67 site independent coordinates and region-groups changed")
    return x

def bounds(q,metres=500):
    lat=q["lat"];lon=q["lon"]
    dlat=metres/111100
    dlon=metres/(111100*math.cos(math.radians(lat)))
    return [lon-dlon,lat-dlat,lon+dlon,lat+dlat]

def get_scenes(session,bbox):
    payload={
       "collections":[COLLECTION],
       "bbox":bbox,
       "datetime":"2022-05-01T00:00:00Z/2022-09-30T23:59:59Z",
       "query":{"eo:cloud_cover":{"lte":40}},
       "limit":100}
    method="POST";url=SOURCE;body=payload
    seen=set();results=[]
    for page in range(9):
        for attempt in range(4):
            try:
                if method=="POST":r=session.post(url,json=body,timeout=(15,100))
                else:r=session.get(url,timeout=(15,100))
                if r.status_code in (429,500,502,503,504):
                    raise requests.RequestException("Temporary source STAC HTTP "+str(r.status_code))
                r.raise_for_status()
                raw=r.json()
                if not isinstance(raw,dict) or not isinstance(raw.get("features"),list):
                    raise RuntimeError("Official STAC item structure not recognized")
                break
            except (requests.RequestException,ValueError,RuntimeError) as ex:
                if attempt==3:raise RuntimeError("Public STAC source inaccessible "+str(ex))
                time.sleep(2+attempt*3)
        for item in raw["features"]:
            if item.get("id") not in seen:
                results.append(item);seen.add(item.get("id"))
        cursor=next((p for p in raw.get("links",[]) if p.get("rel")=="next"),None)
        if not cursor:break
        url=cursor.get("href","")
        if not url.startswith("https://earth-search.aws.element84.com/"):
            raise RuntimeError("Original STAC page next URL source changed")
        method=cursor.get("method","GET").upper()
        body=cursor.get("body",{})
        if method not in ("GET","POST"):
            raise RuntimeError("Original STAC pagination incompatible")
    else:raise RuntimeError("Original STAC source pages not exhausted")
    return results

def valid_scene(item,rectangle):
    if item.get("collection")!=COLLECTION:return None
    pr=item.get("properties",{})
    dt=str(pr.get("datetime",""))
    if not dt.startswith("2022-"):return None
    try:month=int(dt[5:7]);cloud=float(pr.get("eo:cloud_cover"))
    except (ValueError,TypeError):return None
    if month not in MONTHS or not (0<=cloud<=40):return None
    if item.get("geometry") is None:return None
    if not shape(item["geometry"]).covers(rectangle):return None
    assets=item.get("assets",{})
    info={}
    for band in ("red","nir","scl"):
        asset=assets.get(band,{})
        url=asset.get("href","")
        if not isinstance(url,str) or not url.startswith("https://") or not url.endswith(".tif"):
            return None
        orig=asset.get("raster:bands",[])
        if band!="scl" and (len(orig)!=1 or "scale" not in orig[0]):
            return None
        info[band]={"href":url,"source_raster_band_schema":orig[:1],"type":asset.get("type")}
    return {"original_stac_scene_id":item["id"],"original_scene_UTC":dt,
       "month":month,"tile_cloud_fraction_pct":cloud,"original_bands":info,
       "source_STAC_footprint_covers_full_500m_station_buffer":True}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    places=original_sites()
    session=requests.Session()
    session.headers.update({"User-Agent":"gesui-Austrian67-riparian-stac-metadata-source/20261116","Accept":"application/json"})
    counts=[];scenes=[]
    for i,site in enumerate(sorted(places,key=lambda q:q["PKINT"]),1):
        bbox=bounds(site)
        src=get_scenes(session,bbox)
        rectangle=box(*bbox)
        valid=[v for q in src if (v:=valid_scene(q,rectangle)) is not None]
        selected=[]
        for m in MONTHS:
            pick=sorted((q for q in valid if q["month"]==m),
                 key=lambda q:(q["tile_cloud_fraction_pct"],q["original_scene_UTC"],q["original_stac_scene_id"]))[:2]
            for q in pick:
                scenes.append({"original_national_PKINT":site["PKINT"],
                    "original_benthic_geo_block":site["block"],**q})
            selected.extend(pick)
        supported=len({q["month"] for q in valid})>=2 and len({q["original_scene_UTC"][:10] for q in selected})>=2
        counts.append({
          "original_BMLUK_PKINT_2022_true_water_benthos":site["PKINT"],
          "original_EEA2022_lat":site["lat"],"original_EEA2022_lon":site["lon"],
          "fixed_original_benthic_geo_block":site["block"],
          "original_Sentinel2_2022_source_images_found":len(src),
          "original_2022_valid_land_bands_scene_ids":len(valid),
          "2022_valid_source_months":sorted({q["month"] for q in valid}),
          "source_original_selected_monthly_images":len(selected),
          "source_Jun_to_Aug_geographically_source_supported":supported,
          "actual_clear_riverbank_pixels_NOT_YET_VALIDATED":True,
          "no_original_2022_NDVI_or_biotic_effect_read":True})
        print("GESUI_AT67_RIVERBANK_STAC_SOURCE_SITE="+str(i)+"/67 "+site["PKINT"]+
              " months="+str(len(counts[-1]["2022_valid_source_months"])),flush=True)
        time.sleep(.1)
    n=sum(x["source_Jun_to_Aug_geographically_source_supported"] for x in counts)
    per=collections.Counter(q["fixed_original_benthic_geo_block"] for q in counts
                            if q["source_Jun_to_Aug_geographically_source_supported"])
    supported=[k for k,v in per.items() if v>=MIN_PER_SUPPORTED_GROUP]
    gate=n>=MIN_SITES and len(supported)>=MIN_GROUPS
    res={"schema":SCHEMA,
       "status":"PASS_2022_SAME_AQUATIC_STATION_YEAR_S2_SCENE_SOURCE_METADATA" if gate else "STOP_S2_2022_GEO_SCENE_OPPORTUNITY",
       "original_benthic_site_source_gitblob":BQE_BLOB,
       "original_AT_EEA2022_benthic_full_station_count":68,
       "prior_GZUEV_government_vs_EEA2022_coordinate_anomaly_station_excluded_before_S2":EXCLUDED,
       "original_source_2022_same_site_stream_TP_and_EQR_after_geometry_gate":len(places),
       "original_3degree_geo_regions_after_predeclared_coordinate_gate":REGIONS,
       "source_STAC_2022_collection":COLLECTION,
       "source_may_to_september_2022_search":True,
       "fixed_source_Jun_Jul_Aug_months":list(MONTHS),
       "max_2_S2_images_preselected_by_TILE_CLOUD_percent_per_month":True,
       "site_circular_buffer_radius_m":250,
       "site_outer_background_annulus_m":[250,500],
       "scene_tile_could_cover_all_500m_source_count":n,
       "source_scene_supported_sites_by_geographic_region":dict(per),
       "source_regions_with_at_least5_source_sites":len(supported),
       "predeclared_scene_gate_min_sites":MIN_SITES,
       "predeclared_scene_gate_min_regions":MIN_GROUPS,
       "source_gate_pass":gate,
       "source_site_months":counts,
       "source_scenes_preselected_original_metadata_n":len(scenes),
       "optical_actual_2022_SCL_cloud_pixel_greenness_NOT_YET_SAMPLED":True,
       "WorldCover2021_fixed_vegetation_landclass_mask_NOT_YET_SAMPLED":True,
       "not_same_2022_outfall_and_vegetation_causal_NP_delivery":True,
       "no_riparian_green_NDVI_or_insect_or_sewage_EQR_effect_estimated":True}
    ajson=a.output/"AUSTRIA67_EEA2022_SAME_STATION_S2_2022_SCENE_SOURCE_GATE.json"
    ajson.write_text(json.dumps(res,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    (a.output/"AUSTRIA67_S2_2022_JUN_AUG_SELECTED_SCENE_RED_NIR_SCL.json").write_text(
       json.dumps(scenes,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_AUSTRIA67_SAME_STATION_2022_SENTINEL_SCENE_SOURCE="+json.dumps({
        "status":res["status"],"eligible_sites":n,"n_source_regions":len(supported),
        "region_counts":dict(per),"scenes":len(scenes),
        "no_ndvi_or_ecological_effect":True
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
