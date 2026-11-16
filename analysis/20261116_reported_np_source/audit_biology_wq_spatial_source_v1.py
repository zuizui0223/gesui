#!/usr/bin/env python3
"""Original 2026 EEA WQ spatial vs true 2018/22 EEA Biology station source.

SOURCE-ONLY, NO TNTP, EQR, BUTTERFLY OR NUTRI OUTCOMES.
Exact native country/site/scheme joins first. When original Biology site
coordinates exist, use them; otherwise exact original ID in official WQ
spatial to place a subset of Biology stations (report denominator separately).
Distances are PRESELECTED 250/1000/5000 m ONLY. They are Euclidean earth
great-circle proximity, not river-network connectivity, sewage exposure or
same original station chemistry. Never treat near sites as exact IDs.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
import requests
from scipy.spatial import cKDTree

STATION_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
STATION_NAME="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
H7_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
WQ_URL="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_SpatialObjects_DerivedData-csv.zip"
WQ_ZIP_SHA="aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d"
WQ_ZIP_SIZE=2797036
WQ_CSV_NAME="Waterbase_v2025_1_S_WISE6_SpatialObject_DerivedData.csv"
R=6371008.8
CUTOFFS=[250.,1000.,5000.]
SOURCE_PROTOCOL="GESUI_EEA2026_BIOLOGY_WQ_SPATIAL_MATCH_SOURCE_ONLY_V1"
YEARS={"EEA18":(2018,"2017-2019"),"EEA22":(2022,"2020-2022")}
REQUIRED_WQ_FIELDS=["countryCode","monitoringSiteIdentifier",
      "monitoringSiteIdentifierScheme","lat","lon"]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def only(root,name):
    items=list(root.rglob(name))
    if len(items)!=1:raise RuntimeError("Expected exactly one original source "+name)
    return items[0]
def coord(lat,lon):
    try:
        a=float(lat);b=float(lon)
        if not math.isfinite(a) or not math.isfinite(b) or not (-90<=a<=90 and -180<=b<=180):
            return None
        # Treat 0/0 as invalid placeholder for Europe.
        if abs(a)<1e-10 and abs(b)<1e-10:return None
        return (a,b)
    except (TypeError,ValueError,OverflowError):return None
def sphere(c):
    lat,lon=c
    a,b=np.radians([lat,lon])
    return np.array([math.cos(a)*math.cos(b),math.cos(a)*math.sin(b),math.sin(a)])
def haversine(a,b):
    x=np.radians(a);y=np.radians(b)
    da=x[0]-y[0];db=x[1]-y[1]
    h=math.sin(da/2)**2+math.cos(x[0])*math.cos(y[0])*math.sin(db/2)**2
    return 2*R*math.asin(min(1.,math.sqrt(max(0.,h))))
def source_sites(args):
    original=only(args.station_dir,STATION_NAME)
    if digest(original)!=STATION_SHA:raise RuntimeError("Frozen station provenance SHA changed")
    obj=json.loads(original.read_text())
    records=obj.get("candidate_source_site_period_rows")
    if len(records)!=6330:raise RuntimeError("Original 6330 station source records changed")
    panel=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
    if digest(panel)!=H7_SHA:raise RuntimeError("Main original H7 geography source SHA changed")
    with panel.open(encoding="utf-8-sig",newline="") as f:
        p=list(csv.DictReader(f))
    selected={(x["h"],x["period"]):x["h6"] for x in p}
    if len(p)!=454 or len(selected)!=454:raise RuntimeError("Original H7/period panel changed")
    records2=[]
    key_status=collections.Counter()
    src_fields=sorted(set(k for x in records[:200] for k in x))
    coords_names=[
      ("lat","lon"),("latitude","longitude"),
      ("source_lat","source_lon"),("source_latitude","source_longitude"),
      ("station_lat","station_lon"),("station_latitude","station_longitude"),
      ("site_lat","site_lon")
    ]
    candidate=[pair for pair in coords_names if all(n in src_fields for n in pair)]
    for rec in records:
        ty=YEARS.get(rec.get("panel"))
        if ty is None:continue
        year,period=ty
        h=rec.get("HYBAS_ID_L07")
        if h is None or (str(h),period) not in selected:continue
        key=rec.get("native_id")
        if key!=rec.get("physical_id") or not isinstance(key,str):
            raise RuntimeError("Original EEA physical key mismatch")
        split=key.split("|")
        if len(split)!=3 or not all(split):raise RuntimeError("Original EEA station id not CC|ID|scheme")
        if rec.get("country")!=split[0]:raise RuntimeError("Original station country mismatch")
        xy=None;origin=None
        for name in candidate:
            if coord(rec.get(name[0]),rec.get(name[1])) is not None:
                xy=coord(rec.get(name[0]),rec.get(name[1]));origin="original_6330_"+name[0]+name[1]
                break
        records2.append({"key":tuple(split),"year":year,"period":period,"h":str(h),
                          "h6":str(selected[(str(h),period)]),"coords":xy,"coord_source":origin})
    if len(records2)!=3525 or len(set((x["key"],x["year"]) for x in records2))!=3525:
        raise RuntimeError("Frozen 3525 original native Biology site-years not reproduced")
    site=collections.defaultdict(list)
    for r in records2:site[r["key"]].append(r)
    if len(site)!=2485:raise RuntimeError("Frozen 2485 distinct Biology station keys changed")
    return site,records2,{"station_source_fields":src_fields,
      "eligible_original_coord_field_pairs":candidate,
      "source_6330_geometry_SHA256":STATION_SHA,
      "original_biology_native_site_n":len(site),
      "original_biology_native_site_year_n":len(records2),
      "source_biology_sites_with_latlon":sum(
          any(q["coords"] is not None for q in rows) for rows in site.values())}

def wq_sites(outfolder):
    p=outfolder/"EEA_TEMP_2026_WQ_SpatialObjects.zip"
    size=0;dig=hashlib.sha256()
    with requests.get(WQ_URL,headers={
           "User-Agent":"gesui-20261116-original-2026-wq-spatial-site-only/0.1",
           "Accept":"application/zip"},timeout=(12,90),stream=True) as r:
        if r.status_code!=200 or r.url.split("/")[2]!="sdi.eea.europa.eu":
            raise RuntimeError("Official WQ original ZIP unavailable")
        with p.open("wb") as f:
            for block in r.iter_content(65536):
                if not block:continue
                size+=len(block)
                if size>WQ_ZIP_SIZE:raise RuntimeError("Publisher original ZIP too large")
                dig.update(block);f.write(block)
    if size!=WQ_ZIP_SIZE or dig.hexdigest()!=WQ_ZIP_SHA or not zipfile.is_zipfile(p):
        raise RuntimeError("2026 publisher spatial source ZIP SHA/size mismatch")
    d=collections.defaultdict(list)
    raw_coord_nonempty=0
    coordinate_failure=collections.Counter()
    with zipfile.ZipFile(p) as z:
        paths=[q.filename for q in z.infolist() if q.filename.endswith(".csv")]
        if paths!=[WQ_CSV_NAME]:raise RuntimeError("Unexpected published WQ spatial CSV member")
        with z.open(paths[0]) as f:
            reader=csv.DictReader(io.TextIOWrapper(f,encoding="utf-8-sig",newline=""))
            if any(k not in reader.fieldnames for k in REQUIRED_WQ_FIELDS):
                raise RuntimeError("Original WQ scheme/coordinate source fields changed")
            n=0;goodcoord=0;keynonnull=0
            for row in reader:
                n+=1
                key=(row["countryCode"],row["monitoringSiteIdentifier"],
                     row["monitoringSiteIdentifierScheme"])
                if not all(key):continue
                keynonnull+=1
                raw_coord=bool(str(row["lat"]).strip()) and bool(str(row["lon"]).strip())
                if raw_coord:raw_coord_nonempty+=1
                location=coord(row["lat"],row["lon"])
                if location is not None:goodcoord+=1
                elif raw_coord:
                    try:
                        lat=float(row["lat"]);lon=float(row["lon"])
                        if abs(lat)<1e-10 and abs(lon)<1e-10:
                            coordinate_failure["source_zero_zero_placeholder"]+=1
                        elif not math.isfinite(lat) or not math.isfinite(lon):
                            coordinate_failure["source_nonfinite"]+=1
                        elif not (-90<=lat<=90 and -180<=lon<=180):
                            coordinate_failure["source_outside_WGS84"]+=1
                        else:coordinate_failure["source_other_invalid"]+=1
                    except (TypeError,ValueError):
                        coordinate_failure["source_unparseable"]+=1
                d[key].append(location)
    p.unlink()
    if n!=90859 or keynonnull!=88914 or raw_coord_nonempty!=71592:
        raise RuntimeError(f"Original spatial source row or original nonempty coordinate count drift: {n}/{keynonnull}/{raw_coord_nonempty}")
    if goodcoord>raw_coord_nonempty or sum(coordinate_failure.values())!=raw_coord_nonempty-goodcoord:
        raise RuntimeError("Malformed source coordinate validity or discrepancy audit")
    coords={}
    ambiguous=0
    for k,seq in d.items():
        known=[q for q in seq if q is not None]
        if not known:continue
        center=known[0]
        if any(haversine(center,q)>50 for q in known[1:]):
            ambiguous+=1;continue
        coords[k]=center
    return coords,{"EEA_official_2026_spatial_source_SHA256":dig.hexdigest(),
                   "EEA_WQ_spatial_raw_rows":n,
                   "original_usable_site_identifiers":keynonnull,
                   "site_identifier_rows_with_original_nonempty_latlon":raw_coord_nonempty,
                   "site_identifiers_with_valid_latlon_rows":goodcoord,
                   "source_nonempty_but_not_valid_coordinates_n":raw_coord_nonempty-goodcoord,
                   "source_nonempty_but_not_valid_reason_counts":dict(coordinate_failure),
                   "unique_complete_keys_with_unambiguous_coords":len(coords),
                   "keys_with_disagreeing_coords_excluded":ambiguous,
                   "site_original_source_coords_not_republished":True}

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    bio,sy,source=source_sites(args)
    water,wqsource=wq_sites(args.output)
    exact_keys=set(bio)&set(water)
    observations=[]
    for key,rows in bio.items():
        original_co=[r["coords"] for r in rows if r["coords"] is not None]
        if original_co and max(haversine(original_co[0],xy) for xy in original_co)>50:
            raise RuntimeError("Original EEA Biology site moves >50m over original observation years")
        if original_co:
            c=original_co[0];origin="original_EEA_Biology_coordinates"
        elif key in water:
            c=water[key];origin="EEA_WQ_2026_same_exact_Biology_ID_geo"
        else:c=None;origin="unknown_no_coordinates"
        observations.append({"key":key,"coords":c,"source":origin,"h6":{r["h6"] for r in rows},
                             "h7":{r["h"] for r in rows},
                             "years":{r["year"] for r in rows},
                             "siteyears":len(rows)})
    located=[x for x in observations if x["coords"] is not None]
    water_keys=list(water)
    xyz=np.array([sphere(water[k]) for k in water_keys])
    tree=cKDTree(xyz)
    mincount=collections.Counter()
    sites_by_dist={d:set() for d in CUTOFFS}
    site_years_by_dist={d:0 for d in CUTOFFS}
    h7_by_dist={d:set() for d in CUTOFFS}
    h6_by_dist={d:set() for d in CUTOFFS}
    country_by_dist={d:set() for d in CUTOFFS}
    distance_source=collections.Counter()
    exact_with_geo=set(exact_keys)
    nearest=[]
    for q in located:
        distance_source[q["source"]]+=1
        radius=2*math.sin(CUTOFFS[-1]/(2*R))
        ids=tree.query_ball_point(sphere(q["coords"]),radius)
        nearby=[]
        for i in ids:
            key=water_keys[i]
            if key==q["key"] or key[0]!=q["key"][0]:continue
            d=haversine(q["coords"],water[key])
            if d<=CUTOFFS[-1]+1e-7:nearby.append((float(d),key))
        nearby.sort(key=lambda x:x[0])
        # Sites may have thousands of nearby WQ locations. Store only counts.
        best=nearby[0][0] if nearby else None
        nearest.append(best)
        for cutoff in CUTOFFS:
            n=sum(d<=cutoff for d,_ in nearby)
            if n:
                sites_by_dist[cutoff].add(q["key"])
                site_years_by_dist[cutoff]+=q["siteyears"]
                h7_by_dist[cutoff].update(q["h7"])
                h6_by_dist[cutoff].update(q["h6"])
                country_by_dist[cutoff].add(q["key"][0])
                mincount[(str(int(cutoff)),"at_least_one_other_WQ_site")]+=1
    raw_source_count=len(set(k for k in bio))
    if raw_source_count!=2485:raise RuntimeError("Number of genuine Biology sites drifted")
    eligible={
        str(int(d))+"m":{
            "n_original_Biology_site_with_at_least_one_other_WQ_station":len(sites_by_dist[d]),
            "n_original_Biology_site_years_may_have_other_WQ_station":site_years_by_dist[d],
            "n_original_H7_with_nearby_other_WQ_station":len(h7_by_dist[d]),
            "n_true_H6_with_nearby_other_WQ_station":len(h6_by_dist[d]),
            "n_country_groups":len(country_by_dist[d]),
            "same_official_H7_polygon_proven":False,
            "same_river_reach_or_updown_proven":False,
            "water_quality_TNTP_in_target_year_at_nearby_station_NOT_VERIFIED":True,
        } for d in CUTOFFS
    }
    out={
        "schema":SOURCE_PROTOCOL,
        "original_EEA_Biology_source":source,
        "official_2026_WQ_SpatialObjects":wqsource,
        "original_EEA_native_site_years":len(sy),
        "original_EEA_distinct_native_Biology_sites":raw_source_count,
        "exact_same_station_identifier_and_scheme_with_verified_WQ_geography":len(exact_keys),
        "original_Biology_sites_with_any_valid_coords":len(located),
        "bio_coordinate_source_counts":dict(distance_source),
        "source_candidate_nearby_distances_metres":CUTOFFS,
        "distance_gate_source_only":eligible,
        "n_Biology_with_no_valid_geographic_coordinate":len(observations)-len(located),
        "strict_original_2026_same_station_year_TNTP_source_gate_still_STOP":True,
        "original_union_part1_part2_report_positive_H7periods":9,
        "original_same_station_year_TNTP_report_positive_minimum_H7periods":30,
        "source_250m_1000m_5000m_are_not_exposure_or_within_H7_proofs":True,
        "no_butterfly_plant_EIV_HNS_W_or_EQR_effect_opened":True,
        "no_nearby_observational_chemistry_fitted":True,
        "no_causal_sewage_measured_nitrogen_effect":True
    }
    (args.output/"EEA2026_BIOLOGY_VS_WQ_SPATIAL_SOURCE_ONLY_GATE.json").write_text(
        json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_BIOLOGY_WQ_SPATIAL_SOURCE_QA="+json.dumps(out,ensure_ascii=False,sort_keys=True),flush=True)
    print("SPATIAL_CANDIDATES_STATUS=ONLY_SOURCE_OPPORTUNITY_NOT_MEASURED_TNTP_AT_SAME_SITE",flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--station-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
