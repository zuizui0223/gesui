#!/usr/bin/env python3
"""Official original HydroBASINS L7 polygon check for nearby EEA Water Quality stations.

Primary: original 3525/2485 genuine EEA benthic Biology source station-years,
with real original WGS84 site coordinates and HydroBASINS native level7 IDs.
WQ spatial site catalogue: official July2026 ZIP SHA-pinned 90859 raw rows.
Do not claim same monitoring station, same HydroRIVERS reach, hydrologic WWTP
delivery, local field nutrient exposure or plant-butterfly causal effect.

Stage only: freeze ALL original ≤1km same-H7 WQ candidate IDs to source-only
Actions artifact before inspecting chemistry. No biology outcomes read.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import Point

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_reported_np_source"))
from audit_biology_wq_spatial_source_v1 import (
    source_sites,wq_sites,sphere,haversine,R,CUTOFFS,STATION_SHA,WQ_ZIP_SHA
)
sys.path.insert(0,str(ROOT/"analysis/20261115_h6_source_expansion"))
from map_official_eu501_to_h6 import fetch as official_fetch,geom_records as official_geom_records

EXPECTED_RAW_NEAR_SITE={"250":86,"1000":327,"5000":1770}
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
SOURCE_SCHEMA="GESUI_20261116_EEA2026_ORIGINAL_NATIVE_H7_GEOMETRIC_WQ_NEAR_SITE_OPPORTUNITY_V1"

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    biological,siteyears,bsrc=source_sites(a)
    water,wsrc=wq_sites(a.output)
    if (len(biological),len(siteyears))!=(2485,3525):
        raise RuntimeError("Original exact EEA Biology source site/years changed")
    h7_keys={str(rec["h"]) for rec in siteyears}
    if len(h7_keys)!=232:raise RuntimeError("Original HydroBASINS analysis H7 count changed")
    blob,geoqa=official_fetch(7)
    h7_ids,h7_polygons,shapeqa=official_geom_records(blob,7,h7_keys)
    if len(h7_ids)!=232 or set(h7_ids)!=h7_keys or shapeqa["n_unusable"]:
        raise RuntimeError("Official Level7 original 232 polygon geometry source not fully available")
    polys=dict(zip(h7_ids,h7_polygons))
    keys=list(water)
    xyz=np.asarray([sphere(water[k]) for k in keys])
    idx=cKDTree(xyz)
    raw_counter={int(c):set() for c in CUTOFFS}
    h7_counter={int(c):set() for c in CUTOFFS}
    wq_ids={int(c):set() for c in CUTOFFS}
    biologic_sites={int(c):set() for c in CUTOFFS}
    bio_siteyears={int(c):set() for c in CUTOFFS}
    source_h6={int(c):set() for c in CUTOFFS}
    source_country={int(c):set() for c in CUTOFFS}
    co_located=Counter()
    bio_not_in_h7=[]
    candidates=[]
    all_exact_geo=0
    for k,orig in biological.items():
        coordlist=[r["coords"] for r in orig if r["coords"] is not None]
        if len(coordlist)!=len(orig):
            raise RuntimeError("Original EEA native station coordinates incomplete")
        center=coordlist[0]
        if any(haversine(center,xy)>50 for xy in coordlist[1:]):
            raise RuntimeError("Original Biology physical station moved >50m between source periods")
        if k in water:all_exact_geo+=1
        indexes=idx.query_ball_point(sphere(center),2*math.sin(5000/(2*R)))
        nearby=[]
        for v in indexes:
            wkey=keys[v]
            if k==wkey or k[0]!=wkey[0]:continue
            meters=haversine(center,water[wkey])
            if meters<=5000+1e-8:nearby.append((meters,wkey))
        nearby.sort(key=lambda d:(d[0],d[1]))
        for cutoff in CUTOFFS:
            if any(d<=cutoff for d,_ in nearby):raw_counter[int(cutoff)].add(k)
        for record in orig:
            h=str(record["h"]);parent=str(record["h6"]);year=int(record["year"])
            g=polys[h]
            if not g.is_valid:raise RuntimeError("Official H7 original geometry invalid: "+h)
            if not g.covers(Point(center[1],center[0])):
                bio_not_in_h7.append({"country":k[0],"h7":h,"h6":parent,"year":year})
                continue
            for dist,wkey in nearby:
                p=Point(water[wkey][1],water[wkey][0])
                if not g.covers(p):continue
                for cutoff in CUTOFFS:
                    if dist<=cutoff:
                        c=int(cutoff)
                        h7_counter[c].add((h,record["period"]))
                        biologic_sites[c].add(k)
                        bio_siteyears[c].add((k,year))
                        source_h6[c].add(parent)
                        source_country[c].add(k[0])
                        wq_ids[c].add(wkey)
                        if dist<=25:co_located[str(c)]+=1
                if dist<=1000:
                    candidates.append({
                       "original_biology_country":k[0],
                       "original_biology_station_identifier":k[1],
                       "original_biology_station_identifier_scheme":k[2],
                       "original_biology_calendar_year":year,
                       "original_biology_period":record["period"],
                       "original_biology_HYBAS_ID_L07":h,
                       "original_true_HYBAS_ID_L06":parent,
                       "official_WQ_country":wkey[0],
                       "official_WQ_station_identifier":wkey[1],
                       "official_WQ_station_identifier_scheme":wkey[2],
                       "straight_geodesic_distance_m":round(dist,3),
                       "possible_WQ_alias_within25m":dist<=25,
                       "both_original_station_coordinates_in_same_native_H7_polygon":True,
                       "year_specific_TNTP_source_presence_NOT_YET_CHECKED":True,
                       "source_water_reach_or_updown_NOT_YET_CHECKED":True,
                       "different_native_station_ID":True
                    })
    # Strict replay of previous unconditional all-original-geometry opportunity
    # before interpreting now verified same native basin subset.
    if {str(c):len(raw_counter[int(c)]) for c in CUTOFFS}!=EXPECTED_RAW_NEAR_SITE:
        raise RuntimeError("Old 2485 original site-only nearest-WQ pre-H7 counts not reproduced")
    if all_exact_geo!=2485:
        raise RuntimeError("Genuine 2485 same-identifiers official WQ spatial source geolocations changed")
    by={}
    for cutoff in CUTOFFS:
        c=int(cutoff)
        by[str(c)+"m"]={
           "original_biology_sites_with_other_WQ_same_country_and_H7":len(biologic_sites[c]),
           "original_biology_station_years_with_other_WQ_same_H7":len(bio_siteyears[c]),
           "original_H7_periods_with_near_other_WQ_station":len(h7_counter[c]),
           "original_H7_unique_with_near_other_WQ_station":len({x[0] for x in h7_counter[c]}),
           "original_true_H6_unique_with_near_other_WQ_station":len(source_h6[c]),
           "original_country_groups":len(source_country[c]),
           "distinct_official_WQ_station_IDs":len(wq_ids[c]),
           "source_candidates_not_same_EEA_Biology_station":True,
           "not_chemistry_observed":True,
           "not_river_reach_hydrologically_aligned":True}
    if len(candidates)>0:
        key=lambda x:(x["original_biology_country"],x["original_biology_station_identifier"],
                      x["original_biology_station_identifier_scheme"],x["original_biology_calendar_year"],
                      x["official_WQ_country"],x["official_WQ_station_identifier"],
                      x["official_WQ_station_identifier_scheme"])
        if len(set(key(x) for x in candidates))!=len(candidates):
            raise RuntimeError("Duplicate same original Biology station-year×WQ station keys")
    if len(biologic_sites[1000])>len(raw_counter[1000]):
        raise RuntimeError("H7 valid candidate exceeds original raw physical proximity")
    report={
       "schema":SOURCE_SCHEMA,
       "source_original_Biology_station_SHA256":STATION_SHA,
       "source_original_2026_WQ_spatial_ZIP_SHA256":WQ_ZIP_SHA,
       "source_official_HydroBASINS_Level7":geoqa,
       "source_official_polygon_qa":shapeqa,
       "n_original_eea_biology_site_year":len(siteyears),
       "n_original_eea_biology_physical_sites":len(biological),
       "n_original_target_H7":len(h7_keys),
       "n_original_bio_sites_outside_original_native_H7":len(set((x["country"],x["h7"],x["year"]) for x in bio_not_in_h7)),
       "n_original_site_years_outside_H7_geometry":len(bio_not_in_h7),
       "original_bio_outside_H7_diagnostic_examples":bio_not_in_h7[:12],
       "raw_unconditional_other_water_station_250_1000_5000":{str(c):len(raw_counter[int(c)]) for c in CUTOFFS},
       "source_verified_same_NATIVE_H7_opportunity_by_distance":by,
       "same_original_2485_station_IDs_in_WQ_spatial":all_exact_geo,
       "source_candidate_original_bio_year_times_near_WQ_ids_1000m":len(candidates),
       "source_candidate_alias_dist_le25m_n_h7_year_pairs":int(co_located["1000"]),
       "official_2026_wq_station_source_metadata_QA":wsrc,
       "original_Biology_site_geometry_source_QA":bsrc,
       "strict_same_Biology_site_and_year_TNTP_support_original_STILL_9_H7period":True,
       "any_near_WQ_ID_with_year_specific_TNTP_UNKNOWN":True,
       "not_source_same_river_reach_or_outfall_exposure":True,
       "not_spatial_reassignment_of_original_H7":True,
       "no_original_W_EQR_NUTRI_HNS_results_accessed":True,
       "no_causal_or_biodiversity_conclusion":True
    }
    a.output.mkdir(parents=True,exist_ok=True)
    with (a.output/"EEA2026_ORIGINAL_H7_GEOGRAPHIC_NEAR_WQ_CANDIDATE_SITE_YEAR_KEYS_SOURCE_ONLY.csv").open(
        "w",encoding="utf-8",newline="") as f:
        rows=sorted(candidates,key=lambda x:(
           x["original_biology_country"],x["original_biology_station_identifier"],x["original_biology_calendar_year"],
           x["straight_geodesic_distance_m"],x["official_WQ_station_identifier"]))
        if rows:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    (a.output/"EEA2026_OFFICIAL_H7_GEOGRAPHIC_BIOLOGY_NEAR_WQ_SOURCE_GATE.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_EEA2026_SAME_H7_WQ_GEOGRAPHIC_GATE="+json.dumps(
        {k:v for k,v in report.items() if k not in (
           "original_Biology_site_geometry_source_QA","official_2026_wq_station_source_metadata_QA")},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
    print("STATION_YEAR_CHEMISTRY_MODEL_STATUS=NOT_YET_CHECKED_NEARBY_OFFICIAL_WQ_SITES",flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--station-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    main(ap.parse_args())
