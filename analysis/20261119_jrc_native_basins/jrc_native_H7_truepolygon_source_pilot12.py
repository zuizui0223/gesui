#!/usr/bin/env python3
"""H7 x JRC native polygon SOURCE-ONLY coverage pilot; never read EIV/HNS/EQR.

Input: original 226 species-source-eligible H7 IDs, true original GBIF facet
WKT source receipts, full 9739 native JRC point/soil area candidate IDs.
Pilot preselects 12 original H7 deterministically before geometry requests.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
from urllib.parse import urlparse,parse_qs
import pandas as pd
import numpy as np
from pyproj import Transformer
from shapely import wkt
from shapely.ops import transform,unary_union
from shapely.validation import make_valid

sys.path.insert(0,str(Path(__file__).resolve().parent))
from jrc_native_point_soil_polygon30_gate import esri_polygon
from jrc_soilN_waterpoint_native_inventory import WATER,fetch

PANEL=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv")
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
N_PILOT=12
REQUIRED_PASS=8
H6_MIN=4
COVERAGE_MIN=.80
MAX_DUPLICATION=1.05
MAX_BBOX_NATIVE=3000
PAGE=20
MAX_LOCAL_SOURCE_AREA_MULTIPLIER=4
REQUIRE_LOCAL_H7=8
MIN_DISTINCT_DOMINANT_SOURCE_IDS=8

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()

def original_h7(eligfile,facetfile):
    if digest(PANEL)!=PANEL_SHA or digest(facetfile)!=FACETS_SHA:
        raise ValueError("Frozen H7 454 original panel or 2808-facet original source SHA drift")
    p=pd.read_csv(PANEL,usecols=["h","h6","period"],dtype=str)
    if len(p)!=454 or p.h.nunique()!=232 or p.h6.nunique()!=90:
        raise ValueError("Frozen original 454-panel cohort changed")
    e=pd.read_csv(eligfile,usecols=["h","period","full_source_trait_eligible"],dtype=str)
    if len(e)!=454 or e.duplicated(["h","period"]).any():
        raise ValueError("Original EIV source quality support rows changed")
    q=e.full_source_trait_eligible.str.lower().map({"true":True,"false":False})
    if q.isna().any() or q.sum()!=451:
        raise ValueError("Original 451 quality-eligible H7-periods not reproduced")
    e["good"]=q
    d=p[p.period=="2017-2019"].merge(e[["h","period","good"]],on=["h","period"],validate="one_to_one")
    if len(d)!=227 or sum(d.good)!=226:
        raise ValueError("2017-2019 genuine H7/EIV study population not original")
    eligible=d[d.good].copy()
    sample=eligible.h.tolist()
    sample=sorted(sample,key=lambda x:hashlib.sha256(("gesui-20261119-JRC-H7-pilot:"+x).encode()).hexdigest())[:N_PILOT]
    h6=dict(zip(eligible.h,eligible.h6))
    facets=json.loads(Path(facetfile).read_text())
    raw={}
    for row in facets:
        if row.get("group")!="Plantae" or row.get("period")!="2017-2019":continue
        k=str(row["HYBAS_ID_L07"])
        if k not in sample:continue
        if k in raw:raise ValueError("Original study H7 Plantae facet duplicate")
        polys=[]
        for rec in row["source_request_receipts"]:
            arr=parse_qs(urlparse(rec["url"]).query).get("geometry",[])
            if len(arr)!=1 or hashlib.sha256(arr[0].encode()).hexdigest()!=rec["component_SHA"]:
                raise ValueError("Genuine H7 source request polygon component SHA changed")
            polys.append(wkt.loads(arr[0]))
        if not polys:raise ValueError("No original geometry component for H7 "+k)
        v=unary_union(polys)
        if not v.is_valid:v=make_valid(v)
        raw[k]=v
    if len(raw)!=N_PILOT:
        raise ValueError("Source-only selected H7 WKT exact 12 geometry incomplete")
    proj=Transformer.from_crs(4326,3035,always_xy=True).transform
    geo={h:transform(proj,g) for h,g in raw.items()}
    if any(g.is_empty or g.area<=0 for g in geo.values()):
        raise ValueError("Genuine H7 geometry after projection invalid")
    return sample,h6,geo

def candidates(csv_path,json_path,geom_path):
    receipt=json.loads(Path(json_path).read_text())
    geom=json.loads(Path(geom_path).read_text())
    if digest(csv_path)!=receipt["candidateCSVsha256"] or digest(csv_path)!=geom["full_JRC_candidateCSV_SHA256"]:
        raise ValueError("Native JRC candidate original CRC changed")
    if geom["source_polygon_geometry_support"]["sample_polygon_geo_pass_n"]!=30:
        raise ValueError("Independent fixed30 native geometry source cannot establish common JRC native key")
    with open(csv_path,newline="") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=9739 or len({r["soil_hydroid"] for r in rows})!=len(rows):
        raise ValueError("Original JRC common candidate 9739 native basin IDs changed")
    return {int(r["soil_hydroid"]):r for r in rows}

def spatial_point_basin_ids(h7geom):
    # H7 bounding box is used only to discover candidates; actual
    # native polygon∩true H7 geometry determines real source coverage.
    xy=h7geom.bounds
    url=f"{WATER}/5/query"
    base={
      "where":"1=1","geometry":",".join(str(round(v,2)) for v in xy),
      "geometryType":"esriGeometryEnvelope","inSR":3035,
      "spatialRel":"esriSpatialRelIntersects"
    }
    count,csha=fetch(url,{**base,"returnCountOnly":"true"})
    ids,isha=fetch(url,{**base,"returnIdsOnly":"true"})
    n=count.get("count")
    objs=ids.get("objectIds")
    if type(n)is not int or n<0 or n>MAX_BBOX_NATIVE or not isinstance(objs,list) or len(objs)!=n or len(set(objs))!=n:
        raise ValueError("Original ArcGIS bbox source has invalid objectID coverage")
    return sorted(objs),csha,isha

def source_overlaps(geom, fullsource):
    object_ids,ch,ih=spatial_point_basin_ids(geom)
    uri=f"{WATER}/5/query"
    parts=[]
    source_ids=set()
    source_pieces=[]
    count=0
    for j in range(0,len(object_ids),PAGE):
        ids=object_ids[j:j+PAGE]
        result,h=fetch(uri,{
          "objectIds":",".join(map(str,ids)),
          "outFields":"OBJECTID,BasinID,AreaSqKm,PsN",
          "returnGeometry":"true","outSR":3035,"geometryPrecision":2
        })
        features=result.get("features",[])
        if len(features)!=len(ids) or {x.get("attributes",{}).get("OBJECTID") for x in features}!=set(ids):
            raise ValueError("Original ArcGIS native polygon page incomplete")
        for entry in features:
            a=entry.get("attributes",{})
            basin=int(a["BasinID"])
            if basin not in fullsource:continue
            native=esri_polygon(entry.get("geometry"))
            clip=native.intersection(geom)
            if clip.is_empty or clip.area<1e-4:continue
            if basin in source_ids:raise ValueError("Duplicate JRC GREEN BasinID within H7 overlap")
            source_ids.add(basin)
            parts.append(clip)
            original_native_area=float(fullsource[basin]["soil_area_sqkm"])
            if original_native_area<=0:
                raise ValueError("Native basin candidate original area invalid")
            source_pieces.append((basin,clip,original_native_area))
            count+=1
    combined=unary_union(parts) if parts else None
    coverage=float(combined.area/geom.area) if combined and not combined.is_empty else 0.
    rawarea=sum(x.area for x in parts)
    duplication=float(rawarea/combined.area) if combined and combined.area>0 else None
    h7_area_km2=geom.area/1e6
    local_coverage={}
    for multiplier in (1,2,4,10):
        good=[clip for basin,clip,a in source_pieces if a<=multiplier*h7_area_km2]
        cov=unary_union(good).area/geom.area if good else 0.0
        local_coverage[str(multiplier)]=float(cov)
    dominant=max(source_pieces,key=lambda x:x[1].area) if source_pieces else None
    dominated_ID=int(dominant[0]) if dominant else None
    dominant_native_area_ratio=float(dominant[2]/h7_area_km2) if dominant else None
    return {
      "JRC_point_bbox_original_objectids":len(object_ids),
      "JRC_native_P_D_modelled_supported_positive_overlap_ids":count,
      "coverage_union_of_native_basinPolygons_over_original_H7":coverage,
      "sum_over_union_native_polygon_overlap_duplication":duplication,
      "source_area_PASS_80pct":coverage>=COVERAGE_MIN,
      "source_overlap_no_geographic_pseudoreplication_PASS":duplication is not None and duplication<=MAX_DUPLICATION,
      "bbox_source_count_SHA256":ch,"bbox_source_ids_SHA256":ih,
      "native_basinid_sample_first10":sorted(source_ids)[:10],
      "all_native_basinIDs_actual_geometry_intersection":sorted(source_ids),
      "source_area_ratio_at_dominant_native_basin":dominant_native_area_ratio,
      "dominant_native_pointN_soilN_BasinID":dominated_ID,
      "source_polygon_covered_fraction_by_native_area_at_most_multiplied_H7":local_coverage,
      "local_source_spatial_resolution_4x_area80pct_PASS":local_coverage["4"]>=COVERAGE_MIN,
      "not_BQE_or_NUTRI_observed":True
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--facets",type=Path,required=True)
    ap.add_argument("--eiv-eligibility",type=Path,required=True)
    ap.add_argument("--native-source-csv",type=Path,required=True)
    ap.add_argument("--native-source-report",type=Path,required=True)
    ap.add_argument("--native-geom30-report",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    v=ap.parse_args()
    v.out.mkdir(parents=True,exist_ok=True)
    ids,h6,geoms=original_h7(v.eiv_eligibility,v.facets)
    native=candidates(v.native_source_csv,v.native_source_report,v.native_geom30_report)
    results=[]
    for h in ids:
        try:
            z=source_overlaps(geoms[h],native)
            z["source_status"]="ok"
        except Exception as exc:
            z={"source_status":"SOURCE_QUERY_OR_GEOMETRY_STOP","source_failure":str(exc)}
        z["genuine_original_H7"]=h
        z["parent_H6"]=h6[h]
        z["true_H7_area_km2"]=geoms[h].area/1e6
        results.append(z)
        print("GESUI_JRC_H7_NATIVE_SOURCE_PILOT",json.dumps(z,ensure_ascii=False),flush=True)
    acceptable=[x for x in results if x.get("source_area_PASS_80pct")
                and x.get("source_overlap_no_geographic_pseudoreplication_PASS")]
    unique_dominants={x["dominant_native_pointN_soilN_BasinID"] for x in results
                      if x.get("dominant_native_pointN_soilN_BasinID") is not None}
    local_pass=[x for x in results if x.get("local_source_spatial_resolution_4x_area80pct_PASS")]
    qc={
       "source_resolution_post_pilot_diagnostic_NOT_original_coverage_gate":True,
       "native_area_lte_4x_H7_80pct_local_coverage_H7_n":len(local_pass),
       "distinct_native_dominant_sourceIDs":len(unique_dominants),
       "post_pilot_stricter_source_resolution_adequacy":len(local_pass)>=REQUIRE_LOCAL_H7
                             and len(unique_dominants)>=MIN_DISTINCT_DOMINANT_SOURCE_IDS,
       "original_EIV_2017to19_source_only_eligible_H7":226,
       "fixed_pilot_H7_n":N_PILOT,
       "pilot_original_good_H7_geo_n":len(acceptable),
       "pilot_source_GO_for_future_full_H7_coverage":len(acceptable)>=REQUIRED_PASS
                and len({x["parent_H6"] for x in acceptable})>=H6_MIN,
       "pilot_good_H6":len({x["parent_H6"] for x in acceptable}),
       "pilot_bad_H7_n":N_PILOT-len(acceptable),
       "pilot_good_threshold":REQUIRED_PASS,
       "native_2018water_2017to19land_2x2_original_STOP_remains":True,
       "SOURCE_basinN_model_geography_not_actual_discharge_dose":True,
       "biotic_values_never_loaded":True
    }
    report={"schema":"gesui_20261119_H7_vs_native_JRC_point_soilN_truepolygon_source_pilot12_v1",
            "UTC":datetime.now(timezone.utc).isoformat(),
            "original_fixed_pilot_IDs":ids,
            "original_H7_panel_SHA256":digest(PANEL),
            "original_H7_facet_SHA256":digest(v.facets),
            "original_EIV_quality_source_SHA256":digest(v.eiv_eligibility),
            "original_JRC_9739_candidate_SHA256":digest(v.native_source_csv),
            "original_JRC_geom30_SHA256":digest(v.native_geom30_report),
            "pilot":results,"gate":qc}
    (v.out/"JRC_NATIVE_POINTN_SOILN_H7_GENUINE_12_SOURCE_PILOT.json").write_text(
      json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("JRC_H7_NATIVE_SOURCE_PILOT_GATE",json.dumps(qc,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
