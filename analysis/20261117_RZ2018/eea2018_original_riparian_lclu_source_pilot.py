#!/usr/bin/env python3
"""Original source-only 2018 EEA WWTP/benthic stations × 2018 RZ LCLU polygon source.

New observational common-site terrestrial habitat SOURCE QA, NOT a new
WWTP->terrestrial causal ecological model, biodiversity or plant biomass.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,collections,time
from pathlib import Path
from urllib.parse import urlparse
import requests
from shapely.geometry import shape,Point
from shapely.ops import transform,unary_union
from pyproj import Transformer

SOURCE=Path("results/eea2018_independent_source_frozen_20261108/2018_SOURCE_ONLY_BQE_WWTP_NP_UPSTREAM_AND_NEGATIVE_CONTROL.csv")
SOURCE_SHA="5a074f5182711c4ab74ef9bdacfa3918e9fcb434dd0b75c67072a0f06459c3e1"
URL="https://copernicus.discomap.eea.europa.eu/arcgis/rest/services/RiparianZones/RZ_2018/MapServer/0"
LAYER_ID=0
EXPECTED_SOURCE_ROWS=6364
BENTHIC="EEA_13-01-4"
BUFFER_NEAR=250.0
BUFFER_FAR=500.0
COVER_MIN=0.70
EXPECTED_FIELD_NAMES={"OBJECTID","CODE_1_18","CODE_4_18","NODATA_18"}
SCHEMA="GESUI_2018_ORIGINAL_EEA_UWWTD_BENTHIC_RZ2018_250M_500M_SOURCE_ONLY_V1"
PROJ=Transformer.from_crs("EPSG:4326","EPSG:3035",always_xy=True)
UA="gesui2018-official-riparian-sourceonly-20261117-no-ecology"
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def hash_site(s):return hashlib.sha256(s.encode()).hexdigest()
def source(full_census=False):
    if sha(SOURCE)!=SOURCE_SHA:raise RuntimeError("2018 true source-frozen CSV SHA256 changed")
    with SOURCE.open(encoding="utf-8-sig",newline="") as f:
        frame=list(csv.DictReader(f))
    if len(frame)!=EXPECTED_SOURCE_ROWS:raise RuntimeError("Original 2018 EEA source record count drift")
    if len({(r["site_key"],r["observedPropertyDeterminandBiologyEQRCode"]) for r in frame})!=len(frame):
        raise RuntimeError("Original 2018 site×BQE duplicate")
    if any(r["phenomenonTimeReferenceYear"]!="2018" for r in frame):
        raise RuntimeError("Original 2018 reference year mismatch")
    if any(col in frame[0] for col in ("resultNormalisedEQRValue","numeric_EQR_2018","resultEQRValue","normalised_EQR")):
        raise RuntimeError("Source-only fixture must not include numeric 2018 EQR outcomes")
    b=[r for r in frame if r["observedPropertyDeterminandBiologyEQRCode"]==BENTHIC]
    if len(b)<2600 or len({r["site_key"] for r in b})!=len(b):
        raise RuntimeError("Original ~2639 2018 benthic source-only population altered")
    for r in b:
        try:
            r["_lat"]=float(r["lat"]);r["_lon"]=float(r["lon"])
            r["_P"]=float(r["2018_measured_upstream_P_t_per_year"])
            r["_N"]=float(r["2018_measured_upstream_N_t_per_year"])
        except (ValueError,TypeError) as exc:raise RuntimeError("Source missing valid river-site coords or reporter") from exc
        if not (math.isfinite(r["_lat"]) and math.isfinite(r["_lon"]) and r["_P"]>=0 and r["_N"]>=0):
            raise RuntimeError("Invalid 2018 source coordinates/reporting")
    groups=collections.defaultdict(list)
    for r in b:groups[r["countryCode"]].append(r)
    pick=[]
    for country in sorted(groups):
        g=sorted(groups[country],key=lambda r:hash_site(r["site_key"]))
        if len(g)<2: continue
        yes=next((r for r in g if r["_P"]>0),None)
        no=next((r for r in g if r["_P"]==0),None)
        if yes and no:
            pick.extend([yes,no])
        else: pick.extend(g[:2])
    if full_census: pick=sorted(b,key=lambda r:r["site_key"])
    sample={ "original_benthic_all_2018_source_n":len(b),
      "source_original_countries":len(groups),
      "original_source_positive_P_2018_n":sum(r["_P"]>0 for r in b),
      "selected_site_n":len(pick),
      "selected_countries_n":len({r["countryCode"] for r in pick}),
      "prespecified_two_per_country_by_original_P_eligible_and_SHA":not full_census,
      "full_census_all_original_2018_benthic_source_sites":bool(full_census)}
    return pick,sample

def service(session):
    r=session.get(URL,params={"f":"json"},timeout=(12,45));r.raise_for_status()
    d=r.json()
    if d.get("error"): raise RuntimeError("RZ source service error: "+str(d["error"])[:150])
    fields={x.get("name") for x in d.get("fields",[])}
    if (d.get("id")!=LAYER_ID or d.get("geometryType")!="esriGeometryPolygon" or
         not EXPECTED_FIELD_NAMES.issubset(fields)):
        raise RuntimeError("Original Copernicus 2018 LCLU polygon classes missing or wrong layer")
    return {"source_layer_id":d["id"],"source_geometry_type":d["geometryType"],
      "source_field_names_required":sorted(EXPECTED_FIELD_NAMES),
      "source_max_record_count":d.get("maxRecordCount"),
      "source_first_year":2018,
      "source_56_thematic_lclu_classes_from_official_EEA":True}

def source_site(session,r):
    lat,lon=r["_lat"],r["_lon"]
    # Conservative bounding rectangle >500m in latitude/longitude, then exactly
    # clip within 250/500 metre original source station in projected EPSG3035.
    dy=0.006;dx=0.006/max(math.cos(math.radians(lat)),0.25)
    bbox=json.dumps({"xmin":lon-dx,"ymin":lat-dy,"xmax":lon+dx,"ymax":lat+dy,
                     "spatialReference":{"wkid":4326}},separators=(",",":"))
    params={"f":"geojson","where":"1=1","geometry":bbox,
            "geometryType":"esriGeometryEnvelope","inSR":4326,
            "spatialRel":"esriSpatialRelIntersects","outFields":
              "OBJECTID,CODE_1_18,CODE_4_18,NODATA_18",
            "outSR":4326,"returnGeometry":"true","resultRecordCount":2000}
    url=URL+"/query"
    j={"site_key":r["site_key"],"country":r["countryCode"],
       "source_original_region_3deg":r["region_3deg_block"],
       "source_2018_upstream_P_positive_report":bool(r["_P"]>0),
       "source_river_station_geom_valid":True}
    try:
        resp=session.get(url,params=params,timeout=(12,80))
        j["source_http_status"]=resp.status_code
        resp.raise_for_status()
        if not resp.headers.get("Content-Type","").startswith(("application/json","application/geo+json","text/plain")):
            raise RuntimeError("RZ source non JSON/GeoJSON content-type")
        data=resp.json()
        if data.get("error") or data.get("type")!="FeatureCollection":
            raise RuntimeError("Official REST RZ layer query error / not GeoJSON: "+str(data.get("error",data.get("type")))[:130])
        feats=data.get("features",[])
        if data.get("exceededTransferLimit") or len(feats)>=2000:
            raise RuntimeError("Official RZ query paginated/truncated; no approximate area allowed")
        x,y=PROJ.transform(lon,lat)
        pt=Point(x,y)
        near=pt.buffer(BUFFER_NEAR,resolution=32)
        outer=pt.buffer(BUFFER_FAR,resolution=32).difference(near)
        area_polys=[]
        classes=set()
        landcover_group=collections.defaultdict(list)
        for f in feats:
            props=f.get("properties") or {}
            if any(k not in props for k in ("CODE_1_18","CODE_4_18","NODATA_18")):
                raise RuntimeError("LCLU2018 classification fields not returned")
            if props["NODATA_18"] not in (0,"0",None):continue
            code=props["CODE_4_18"]
            if code is None:continue
            g=transform(PROJ.transform,shape(f["geometry"]))
            if not g.is_valid:g=g.buffer(0)
            if g.is_empty:continue
            area_polys.append(g);classes.add(int(code))
            landcover_group[int(props["CODE_1_18"])].append(g)
        if area_polys:
            multi=unary_union(area_polys)
            nc=float(multi.intersection(near).area/near.area)
            oc=float(multi.intersection(outer).area/outer.area)
        else:nc=0.;oc=0.
        if not (-1e-6<=nc<=1.00001 and -1e-6<=oc<=1.00001):
            raise RuntimeError("Original LCLU RZ polygons overlap or wrong area")
        # All Level1 terrestrial habitat classes are original publisher labels,
        # NOT an optimized biological response or measured plant biomass.
        by_near={str(k):float(unary_union(gg).intersection(near).area/near.area)
                 for k,gg in sorted(landcover_group.items())}
        by_outer={str(k):float(unary_union(gg).intersection(outer).area/outer.area)
                 for k,gg in sorted(landcover_group.items())}
        if sum(by_near.values())>1.0001 or sum(by_outer.values())>1.0001:
            raise RuntimeError("Published polygon Level1 categories overlap")
        j.update({"original_polygon_features_n":len(feats),
          "source_Level1_class_area_frac_of_entire_near250m_circle":by_near,
          "source_Level1_class_area_frac_of_entire_outer250to500m_ring":by_outer,
          "original_classified_features_n":len(area_polys),
          "original_mapped_unique_classes_n":len(classes),
          "original_2018_near250m_LCLU_coverage":nc,
          "original_2018_outer250to500m_LCLU_coverage":oc,
          "source_near_outer_0p70_eligible":bool(nc>=COVER_MIN and oc>=COVER_MIN),
          "original_EEA_RZ2018_not_yearly_plant_production":True})
    except Exception as exc:
        j["source_error"]=type(exc).__name__+": "+str(exc)[:220]
        j["source_near_outer_0p70_eligible"]=False
    return j

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    result={"schema":SCHEMA,"status":"STARTED_SOURCE_ONLY_NO_ECOLOGY",
      "source_original_Github_2018_fixed_data_sha256":SOURCE_SHA,
      "source_original_RZ2018_provider_url":URL,
      "source_Frozen_radii_m":[BUFFER_NEAR,BUFFER_FAR],
      "source_frozen_polygon_coverage_min":COVER_MIN,
      "prelock_GO_min_eligible_sites":10,
      "prelock_GO_min_countries":5,
      "prelock_GO_min_independent_spatial3deg_blocks":5,
      "prelock_GO_min_source_report_positive_sites":3,
      "no_EQR_or_WWTP_biological_effect_model_fitted":True,
      "not_a_terrestrial_biomass_or_species_richness_measure":True,
      "prior_AT67_S2_riparian_source_STOP_unchanged":True}
    try:
        chosen,qa=source(a.full_census)
    except Exception as exc:
        result.update({"status":"STOP_ORIGINAL_EEA2018_SOURCE",
          "source_error":type(exc).__name__+": "+str(exc)[:240]})
        save(args.output,result);return
    result["source_original_EEA2018_2018_benthic_population"]=qa
    session=requests.Session()
    session.headers.update({"User-Agent":UA,"Accept":"application/json"})
    try:result["source_RZ2018_official_layer_metadata"]=service(session)
    except Exception as exc:
        result.update({"status":"STOP_OFFICIAL_RZ2018_PUBLIC_REST_METADATA",
          "source_error":type(exc).__name__+": "+str(exc)[:240]})
        save(args.output,result);return
    xx=[]
    for i,z in enumerate(chosen,1):
        row=source_site(session,z);xx.append(row)
        print("GESUI_RZ2018_PILOT_"+str(i)+"_of_"+str(len(chosen))+
          " country="+str(row["country"])+" eligible="+str(row["source_near_outer_0p70_eligible"])+
          " polygons="+str(row.get("original_polygon_features_n",0)),flush=True)
    result["original_source_pilot_true_station_classification"]=xx
    good=[r for r in xx if r["source_near_outer_0p70_eligible"]]
    result["source_pilot_eligible_n"]=len(good)
    result["source_pilot_eligible_countries"]=len({x["country"] for x in good})
    result["source_pilot_eligible_3degree_regions"]=len({x["source_original_region_3deg"] for x in good})
    result["source_pilot_eligible_upstream_P_positive_sites"]=sum(x["source_2018_upstream_P_positive_report"] for x in good)
    if a.full_census:
        passed=(len(good)>=200 and result["source_pilot_eligible_countries"]>=6
             and result["source_pilot_eligible_3degree_regions"]>=10
             and result["source_pilot_eligible_upstream_P_positive_sites"]>=20)
        result["status"]=("PASS_2018_RZ2018_ALL_ORIGINAL_BENTHIC_SITES_SOURCE_ONLY_NO_ECOLOGY"
                          if passed else "STOP_PRELOCKED_FULL_2018_RZ_SOURCE_SUPPORT")
        result["source_frozen_full_census_thresholds"]={
            "min_original_2018_benthic_sites_with_both_rings":200,
            "min_countries":6,"min_spatial3degree_regions":10,
            "min_original_report_positive_upstream_sites":20}
    else:
        passed=(len(good)>=10 and result["source_pilot_eligible_countries"]>=5
            and result["source_pilot_eligible_3degree_regions"]>=5
            and result["source_pilot_eligible_upstream_P_positive_sites"]>=3)
        result["status"]=("PASS_2018_PUBLIC_RIPARIAN_LCLU_SOURCE_PILOT_ONLY_NO_ECOLOGY"
            if passed else "STOP_PRELOCKED_2018_COMMON_SITE_RZ_PILOT_SOURCE_SUPPORT")
    save(args.output,result)

def save(folder,result):
    (folder/"EEA2018_2018_RIPARIAN_ZONES_REST_OFFICIAL_SOURCE_ONLY_QA.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2018_RZ2018_COMMON_SOURCE_GATE="+json.dumps({
       k:result.get(k) for k in ("status","source_original_EEA2018_2018_benthic_population",
           "source_pilot_eligible_n","source_pilot_eligible_countries",
           "source_pilot_eligible_3degree_regions","source_pilot_eligible_upstream_P_positive_sites",
           "source_error")},
       ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--full-census",action="store_true",help="New source-first 2639 original station census")
    a=p.parse_args()
    main(a)
