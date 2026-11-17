#!/usr/bin/env python3
"""Outcome-blind 2018/2022 public Sentinel-2 source-availability screen for FR263.

Exact original EEA station coordinates, H6-group balanced hash sample fixed before
opening imagery, and no TP, benthic, vegetation or NDVI magnitude in this stage.
Source availability alone does NOT pass the genuine cloud-free riverbank pixel gate.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
import math
import pathlib
import time
import pandas as pd
import requests
from shapely.geometry import box,shape

SCHEMA="GESUI_20261117_FR263_S2_2018_2022_24H6_PILOT_SOURCE_OPPORTUNITY_V1"
API="https://earth-search.aws.element84.com/v1/search"
COLLECTION="sentinel-2-l2a"
MONTHS=(6,7,8)
MIN_SUPPORTED_H6=18
ORIGINAL_GITBLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
SOURCE_DATE_SHA="2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264"

def sha(b):return hashlib.sha256(b).hexdigest()

def source_sites(geo,datesource):
    raw=datesource.read_bytes()
    if sha(raw)!=SOURCE_DATE_SHA:raise RuntimeError("Original 263 source dates SHA mismatch")
    d=json.loads(raw)["original_source_only_exact_same_EEA2018_2022_biological_site_and_annual_measured_TP_support"]
    if len(d)!=263:raise RuntimeError("Original 263 source date list changed")
    original=pd.DataFrame([{
       "site_key":z["EEA2022_site_key"],"H7":str(z["HYBAS_ID_L07"]),
       "H6":str(z["HYBAS_ID_L06_parent"])} for z in d])
    coords=pd.read_csv(geo,dtype={"site_key":str,"HYRIV_ID":str})
    if len(coords)!=263 or coords.site_key.nunique()!=263 or coords.HYRIV_ID.nunique()!=262:
        raise RuntimeError("Original pre-outcome 263 geography changed")
    if coords.source_EEA2022_Git_blob_sha1.nunique()!=1 or coords.source_EEA2022_Git_blob_sha1.iloc[0]!=ORIGINAL_GITBLOB:
        raise RuntimeError("Original EEA geographic source blob no longer matches")
    joined=original.merge(coords,on="site_key",validate="1:1")
    if len(joined)!=263 or joined.H6.nunique()!=24 or joined.H7.nunique()!=73:
        raise RuntimeError("Original 73 H7 / 24 H6 source lineage differs")
    return joined

def bbox(r,metres=1000):
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

def run(args):
    stations=source_sites(args.geo,args.datesource)
    chosen=(stations.assign(source_hash=stations.site_key.map(lambda x:sha(x.encode())))
                   .sort_values(["H6","source_hash","site_key"]).drop_duplicates("H6").sort_values("H6"))
    if len(chosen)!=24:raise RuntimeError("Original 24 H6 sample not exactly one/H6")
    sess=requests.Session()
    sess.headers.update({"User-Agent":"gesui-fr263-source-only-s2-stac-20261117/1",
                         "Accept":"application/geo+json,application/json"})
    records=[]
    for i,site in enumerate(chosen.itertuples(index=False),1):
        attempts={};counts={}
        for year in (2018,2022):
            q=search_scenes(sess,site,year)
            counts[str(year)]=q
        passes=all(len(counts[str(year)]["selected_months"])>=2 for year in (2018,2022))
        records.append({
          "site_key":site.site_key,"H6":site.H6,"H7":site.H7,
          "HYRIV_ID":site.HYRIV_ID,
          "original_lat":float(site.lat),"original_lon":float(site.lon),
          "source_supported_both_2018_and2022_two_summer_months":passes,
          "years":counts})
        print(f"FR263_S2_24H6_SOURCE={i}/24 {site.H6} eligible_both={passes} months2018={counts['2018']['selected_months']} months2022={counts['2022']['selected_months']}",flush=True)
    eligible=sum(z["source_supported_both_2018_and2022_two_summer_months"] for z in records)
    out={
      "schema":SCHEMA,
      "status":"PASS_SCENE_METADATA_ONLY_NO_PIXELS" if eligible>=MIN_SUPPORTED_H6 else "STOP_SCENE_METADATA_FEASIBILITY_NO_PIXELS",
      "predeclared_original_sample":"One original FR263 physical site selected by SHA256(site_key) per each of 24 native H6 parents before reading S2 scenes, EQR, TP or NDVI",
      "original_FR263_site_universe":263,"original_H7":73,"original_H6":24,
      "n_H6_fixed_pilot":24,"min_bothyear_supported_H6_for_followup":MIN_SUPPORTED_H6,
      "source_supported_bothyear_pilot_n":eligible,
      "original_years":[2018,2022],"fixed_months":[6,7,8],
      "original_STAC_collection":COLLECTION,
      "fixed_max_scene_tile_cloud":40,
      "max_original_scenes_per_month_year_site":2,
      "optical_pilot_month_gate_at_least2_each_year":True,
      "valid_actual_pixel_coverage_unknown":True,
      "original_riverline_bank_distance_NOT_YET_EXTRACTED":True,
      "any_2018_2022_plant_productivity_or_NDVI_NOT_COMPUTED":True,
      "full_original_source_sample_after_metadata_pilot_is_STILL_263":True,
      "source_station_selection_sha256_geo":sha(args.geo.read_bytes()),
      "source_station_selection_sha256_date":SOURCE_DATE_SHA,
      "pilot_sites":records}
    args.out.mkdir(parents=True,exist_ok=True)
    path=args.out/"FR263_24H6_S2_2018_2022_ORIGINAL_STAC_METADATA_SOURCE_QA.json"
    path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    print("FR263_S2_PUBLIC_METADATA_RESULT "+json.dumps({
        k:out[k] for k in ["status","n_H6_fixed_pilot","source_supported_bothyear_pilot_n","min_bothyear_supported_H6_for_followup"]}),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--geo",required=True,type=pathlib.Path)
    p.add_argument("--datesource",required=True,type=pathlib.Path)
    p.add_argument("--out",required=True,type=pathlib.Path)
    run(p.parse_args())

if __name__=="__main__":main()
