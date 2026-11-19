#!/usr/bin/env python3
"""Fixed 30-ID polygon geometry audit for JRC soil HydroID vs point BasinID.

All candidates come from a full-data, source-only native value+area match.
The API polygon equality check is required before claiming common native basins.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime,timezone

from shapely.geometry import Polygon
from shapely.ops import unary_union
from shapely.validation import make_valid
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from jrc_soilN_waterpoint_native_inventory import CONFIG,fetch

SAMPLE_N=30
MIN_PASS=27
AREA_MATCH_MIN=.95
POLYGON_NATIVE_AREA_MAX_REL_ERROR=.20

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def signed_area(ring):
    if not ring or len(ring)<4:return 0.
    return sum((float(a[0])*float(b[1])-float(b[0])*float(a[1])) for a,b in zip(ring,ring[1:]))*.5

def esri_polygon(geom):
    if not isinstance(geom,dict) or not isinstance(geom.get("rings"),list):
        raise ValueError("ArcGIS original polygon rings missing")
    ext,holes=[],[]
    for points in geom["rings"]:
        if len(points)<4:raise ValueError("ArcGIS polygon ring too short")
        poly=Polygon(points)
        if poly.is_empty or poly.area<=0:continue
        if signed_area(points)<0:ext.append(poly)
        else:holes.append(poly)
    if not ext:
        raise ValueError("No ESRI clockwise outer ring; ambiguous geometry")
    out=unary_union(ext)
    if holes:out=out.difference(unary_union(holes))
    if not out.is_valid:out=make_valid(out)
    if out.is_empty or out.area<=0:raise ValueError("Original ESRI polygon invalid")
    return out

def original_geom(name,source_id):
    base,layer,key,value,area=CONFIG[name]
    url=f"{base}/{layer}/query"
    result,request_sha=fetch(url,{
      "where":f"{key}={int(source_id)}",
      "outFields":",".join(["OBJECTID",key,value,area]),
      "returnGeometry":"true",
      "outSR":3035,
      "geometryPrecision":2,
      "resultRecordCount":2
    })
    rows=result.get("features",[])
    if len(rows)!=1 or rows[0].get("attributes",{}).get(key)!=source_id:
        raise ValueError("Original JRC polygon query missing or duplicate native ID")
    geom=esri_polygon(rows[0].get("geometry"))
    att=rows[0]["attributes"]
    native_area=float(att[area])
    if not math.isfinite(native_area) or native_area<=0:
        raise ValueError("Original basin area invalid")
    return geom,att,request_sha

def check(source_csv,source_report):
    receipt=json.loads(Path(source_report).read_text())
    if not receipt.get("no_EEA_EQR_EIV_NUTRI_HNS_outcome_opened"):
        raise ValueError("Not source-only full JRC census")
    if receipt.get("candidateCSVsha256")!=sha(source_csv):
        raise ValueError("Original JRC full-source candidate CSV SHA changed")
    with open(source_csv,newline="") as f:
        allrows=list(csv.DictReader(f))
    ids=[int(x["soil_hydroid"]) for x in allrows]
    if len(ids)!=len(set(ids)) or len(ids)<SAMPLE_N:
        raise ValueError("Native area-matching sample size/identity invalid")
    # Deterministic SOURCE-independent sample; pick no IDs based on data values.
    pick=sorted(ids,key=lambda v:(hashlib.sha256(f"gesui-fixed-20261119-{v}".encode()).hexdigest(),v))[:SAMPLE_N]
    return pick

def audit(ids):
    matched=[]
    all_geos=[]
    for sid in ids:
        geo1,att1,h1=original_geom("soil_atmosN_2015",sid)
        geo2,att2,h2=original_geom("water_pointN_2014_18",sid)
        numerator=geo1.intersection(geo2).area
        f1=float(numerator/geo1.area)
        f2=float(numerator/geo2.area)
        a1=float(att1["Basin_3035_AreaSqKm"])
        a2=float(att2["AreaSqKm"])
        measured_ratio1=geo1.area/(a1*1e6)
        measured_ratio2=geo2.area/(a2*1e6)
        geom_ok=f1>=AREA_MATCH_MIN and f2>=AREA_MATCH_MIN
        area_consistent=(abs(1-measured_ratio1)<=POLYGON_NATIVE_AREA_MAX_REL_ERROR
                         and abs(1-measured_ratio2)<=POLYGON_NATIVE_AREA_MAX_REL_ERROR)
        row={
          "native_integer_ID_not_complete_identity_evidence":sid,
          "jrc_soilHydroID_to_JRC_pointBasinID_tested":True,
          "soil_hydro_vs_point_intersection_fraction_of_soil":f1,
          "soil_hydro_vs_point_intersection_fraction_of_point":f2,
          "soil_polygon_projection_area_over_published_area":measured_ratio1,
          "point_polygon_projection_area_over_published_area":measured_ratio2,
          "actual_geometry_bidirectional_area95_PASS":geom_ok,
          "original_native_area_consistent80to120pct":area_consistent,
          "sample_geo_pass":geom_ok and area_consistent,
          "original_soil_query_SHA256":h1,
          "original_point_query_SHA256":h2,
        }
        matched.append(row)
        all_geos.append((sid,geo1))
        print("GESUI_JRC_POLYGON_SOURCE_MATCH",json.dumps({
          "id":sid,"soil_fraction":f1,"point_fraction":f2,
          "soil_area_rel":measured_ratio1,"point_area_rel":measured_ratio2,
          "sample_pass":row["sample_geo_pass"]
        }),flush=True)
    nested_pairs=[]
    for i,(a,pa) in enumerate(all_geos):
        for b,pb in all_geos[i+1:]:
            f=pa.intersection(pb).area/min(pa.area,pb.area)
            if f>.05:
                nested_pairs.append({"id_a":a,"id_b":b,
                     "intersection_over_smaller_area":float(f)})
    good=sum(x["sample_geo_pass"] for x in matched)
    support={
      "original_polygon_sample_tested":len(matched),
      "sample_polygon_geo_pass_n":good,
      "sample_polygon_geo_pass_threshold":MIN_PASS,
      "bidirectional_overlap_min_ratio":AREA_MATCH_MIN,
      "source_polygon_pilot_geospatial_match_PASS":good>=MIN_PASS,
      "sample_pairwise_native_overlapping_or_nested_basin_pairs_n":len(nested_pairs),
      "sample_pairs_that_are_independent_nonoverlap_NOT_GUARANTEED":True,
      "all_source_ids_geometry_equivalence_NOT_PROVED":True,
      "2018water_2017to19plant_common_biotic_2x2_STOP":True,
    }
    return matched,nested_pairs,support

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-csv",type=Path,required=True)
    p.add_argument("--source-report",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    ids=check(a.source_csv,a.source_report)
    rows,nested,support=audit(ids)
    with (a.out/"JRC_NATIVE_HYDROID_POINTBASINID_GEO_PILOT30.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    receipt={
      "schema":"gesui_JRC_native_soil_HydroID_point_BasinID_fixed30_source_geom_audit_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "full_JRC_candidateCSV_SHA256":sha(a.source_csv),
      "original_full_JRC_source_gate_report_SHA256":sha(a.source_report),
      "original_deterministic_sample_id_list":ids,
      "source_polygon_geometry_support":support,
      "native_hydrology_nested_geometry_sample_pairs":nested,
      "no_biological_EQR_EIV_butterfly_values_opened":True,
    }
    (a.out/"JRC_NATIVE_SOIL_POINT_GEOMETRY_SAMPLE30_GATE.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("JRC_NATIVE_SOURCE_30POLYGON_GATE",json.dumps(support,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
