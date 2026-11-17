#!/usr/bin/env python3
"""Authentic 2018+2022 Sentinel-2 STAC source opportunity for fixed FR178.

The cohort is source-frozen before seeing TP/EQR/NDVI. Four stable H6 shards,
no source substitutions and no invented cloud-free land pixel availability.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time
from pathlib import Path
from types import SimpleNamespace
import requests
from shapely.geometry import box,shape

SCHEMA="GESUI_FR178_2018_2022_S2_SCENE_METADATA_SOURCE_SHARD_V1"
INPUT_SHA="9a5d11a6c50ba49160da46097873c3d5eedda881e0a173f1aa90968c94950b11"
SOURCE_ID_SHA="76ba45169e58901a2cd82a18244d022773074963a6ac983b527f06302a4aedb8"
API="https://earth-search.aws.element84.com/v1/search"
COLLECTION="sentinel-2-l2a"
MONTHS=(6,7,8)
N_SHARDS=4

def sha(b):return hashlib.sha256(b).hexdigest()

def bbox(r,metres=1800):
    lat=float(r.lat);lon=float(r.lon)
    dlat=metres/111132.0
    dlon=metres/(111320.0*math.cos(math.radians(lat)))
    return [lon-dlon,lat-dlat,lon+dlon,lat+dlat]

def search_scenes(session,r,year):
    geom=bbox(r)
    query={
      "collections":[COLLECTION],
      "bbox":geom,
      "datetime":f"{year}-05-01T00:00:00Z/{year}-09-30T23:59:59Z",
      "query":{"eo:cloud_cover":{"lte":40}},
      "limit":100}
    all_items=[];seen=set();url=API;method="POST";body=query
    for page in range(7):
        for attempt in range(4):
            try:
                x=session.post(url,json=body,timeout=(15,80)) if method=="POST" else \
                  session.get(url,timeout=(15,80))
                if x.status_code in (429,500,502,503,504):
                    raise requests.RequestException("Transient STAC status="+str(x.status_code))
                x.raise_for_status()
                parsed=x.json()
                if not isinstance(parsed.get("features"),list):
                    raise RuntimeError("Not STAC feature list")
                break
            except (requests.RequestException,ValueError,RuntimeError) as e:
                if attempt==3:raise RuntimeError("Source STAC read unavailable: "+str(e))
                time.sleep(2+attempt*2)
        for item in parsed["features"]:
            if item.get("id") not in seen:
                all_items.append(item);seen.add(item.get("id"))
        nextp=next((a for a in parsed.get("links",[]) if a.get("rel")=="next"),None)
        if not nextp:break
        url=nextp.get("href","")
        if not url.startswith("https://earth-search.aws.element84.com/"):
            raise RuntimeError("Unexpected STAC next-page URL")
        method=nextp.get("method","GET").upper()
        body=nextp.get("body",{})
        if method not in ("GET","POST"):raise RuntimeError("STAC next pagination method unsupported")
    else:raise RuntimeError("STAC pagination >7 pages; stop source")
    out=[]
    for it in all_items:
        props=it.get("properties",{})
        dt=str(props.get("datetime",""))
        if not dt.startswith(str(year)+"-"):continue
        try:month=int(dt[5:7]);cloud=float(props["eo:cloud_cover"])
        except (ValueError,KeyError,TypeError):continue
        if month not in MONTHS or cloud<0 or cloud>40:continue
        if it.get("collection")!=COLLECTION or not it.get("geometry"):continue
        try:
            if not shape(it["geometry"]).covers(box(*geom)):continue
        except ValueError:continue
        assets=it.get("assets",{});item_bands={}
        for k in ("red","nir","scl"):
            a=assets.get(k,{})
            href=a.get("href","")
            if not isinstance(href,str) or not href.startswith("https://sentinel-cogs.s3.us-west-2.amazonaws.com/") or not href.endswith(".tif"):
                break
            rb=a.get("raster:bands") or []
            if k!="scl" and (len(rb)!=1 or "scale" not in rb[0]):
                break
            item_bands[k]={"href":href,"raster_bands":rb[:1]}
        if len(item_bands)!=3:continue
        out.append({"scene_id":it.get("id"),"datetime":dt,"month":month,
             "cloud_pct":cloud,"bands":item_bands})
    picks=[]
    for m in MONTHS:
        group=sorted((a for a in out if a["month"]==m),
                     key=lambda a:(a["cloud_pct"],a["datetime"],a["scene_id"]))[:2]
        picks.extend(group)
    return {"stac_items_received":len(all_items),
            "valid_original_items":len(out),
            "selected_months":sorted({x["month"] for x in picks}),
            "source_selected_scenes":picks}


def run(a):
    raw=a.source.read_bytes()
    if sha(raw)!=INPUT_SHA:raise RuntimeError("FR178 original fixed source roster SHA changed")
    rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    if len(rows)!=178 or len({x["site_key"] for x in rows})!=178 or len({x["H6"] for x in rows})!=23:
        raise RuntimeError("Original fixed 178 same physical site frame invalid")
    ids=sorted(x["site_key"] for x in rows)
    if sha("\n".join(ids).encode())!=SOURCE_ID_SHA:
        raise RuntimeError("Original source site IDs or sorting changed")
    if a.shard<0 or a.shard>=N_SHARDS:
        raise ValueError("Unsupported original H6 shard")
    cohort=sorted((x for x in rows if int(sha(x["H6"].encode())[:8],16)%N_SHARDS==a.shard),
                  key=lambda x:(x["H6"],x["site_key"]))
    if not cohort:raise RuntimeError("No source H6 for fixed shard")
    ses=requests.Session()
    ses.headers.update({"User-Agent":"gesui-original-FR178-river-bank-2018-2022-STAC-scene-only/20261117","Accept":"application/geo+json,application/json"})
    result=[]
    for i,original in enumerate(cohort,1):
        position=SimpleNamespace(lat=float(original["eaa_lat"]),lon=float(original["eea_lon"]))
        yrs={}
        for year in (2018,2022):
            try:
                yrs[str(year)]=search_scenes(ses,position,year)
                yrs[str(year)]["source_request_success"]=True
            except (requests.RequestException,RuntimeError,ValueError) as e:
                yrs[str(year)]={"source_request_success":False,"source_error_type":type(e).__name__,
                                  "source_error":str(e)[:300],
                                  "selected_months":[],"source_selected_scenes":[]}
        joined=all(yrs[str(y)]["source_request_success"] and len(yrs[str(y)]["selected_months"])>=2
                   for y in (2018,2022))
        result.append({
           "site_key":original["site_key"],"H6":original["H6"],"H7":original["H7"],
           "HYRIV_ID":original["HYRIV_ID"],"original_lat":position.lat,"original_lon":position.lon,
           "original_JJA_same_calendar_months":[int(z) for z in original["original_JJA_same_calendar_months"].split(",")],
           "original_reach_offset_m":float(original["original_station_to_original_reach_m"]),
           "STAC_scene_metadata_two_years_gate":joined,"years":yrs})
        print("GESUI_FR178_STAC_SOURCE_SHARD="+str(a.shard)+" SITE="+str(i)+"/"+str(len(cohort))+
              " SUPPORT="+str(joined),flush=True)
    eligible=[x for x in result if x["STAC_scene_metadata_two_years_gate"]]
    report={"schema":SCHEMA,"status":"SOURCE_STAC_SCENE_SHARD_COMPLETED_NO_PIXELS",
      "source_FR178_frozen_csv_SHA256":INPUT_SHA,"original_178_site_population":178,
      "original_58_H7_23_H6":True,"source_shard_index":a.shard,"source_shard_count":N_SHARDS,
      "source_shard_count_sites":len(cohort),"source_eligible_bothyear_2month_sites":len(eligible),
      "source_scene_search_bbox_halfwidth_m":1800,"max_two_original_scenes_per_JJA_month":True,
      "max_original_STAC_scene_cloud_fraction":40,
      "cloudfree_pixel_eligibility_not_yet_measured":True,
      "no_river_TP_EQR_butterfly_or_NDVI_outcomes_read":True,
      "source_sites":result}
    a.out.mkdir(parents=True,exist_ok=True)
    path=a.out/f"FR178_ORIGINAL_STAC_2018_2022_SHARD{a.shard}.json"
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_FR178_STAC_FINISHED="+json.dumps({"shard":a.shard,"n":len(cohort),"two_years_source":len(eligible)}),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
if __name__=="__main__":main()
