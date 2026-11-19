#!/usr/bin/env python3
"""Strict source-only 12-H7 geometry pilot against JRC native basin features.

NEVER used for ecological inference or cherry-picking source/model regions.
The 12 H7 are chosen deterministically from the SOURCE quality-supported
2017-2019 original cohort, without EIV/HNS/EQR outcome values.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time,urllib.parse,urllib.request
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import geopandas as gpd
from shapely import wkt
from shapely.geometry import shape
from shapely.ops import unary_union,transform
from pyproj import Transformer

BASE="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer"
WATER="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
SOURCES={
 "atmos2015":(BASE,5,"MEAN_mgsqmy","Basin_3035_HydroID"),
 "mineralN":(BASE,10,"MEAN_kgha","Basin_3035_HydroID"),
 "organicN":(BASE,15,"MEAN_kgha","Basin_3035_HydroID"),
 "fixationN":(BASE,20,"MEAN_kgha","Basin_3035_HydroID"),
 "pointN":(WATER,5,"PsN","BasinID"),
}
FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
PAGE=150
N_PILOT=12
MAX_CANDIDATES_PER_H7_LAYER=1600
PROJECT=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform

def hashpath(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def rq(url,params):
    url=url+"?"+urllib.parse.urlencode(params)
    last=None
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"GESUI-JRC-native-H7-geometry-pilot-20261119",
                "Accept":"application/json"
            })
            with urllib.request.urlopen(req,timeout=80) as fp:raw=fp.read(28_000_000)
            data=json.loads(raw)
            if not isinstance(data,dict) or data.get("error"):
                raise ValueError("JRC ArcGIS JSON problem "+str(data.get("error")))
            return data,hashlib.sha256(raw).hexdigest()
        except (ValueError,OSError,TimeoutError) as e:
            last=e
            if attempt<3:time.sleep(1+attempt)
    raise RuntimeError("JRC official native basin spatial query unavailable: "+str(last))

def original_h7(facets,elig,panel):
    if hashpath(facets)!=FACETS_SHA or hashpath(panel)!=PANEL_SHA:
        raise ValueError("Original H7 source facet or panel drift")
    import pandas as pd
    d=pd.read_csv(panel,usecols=["h","h6","country","period"],dtype=str)
    t=pd.read_csv(elig,usecols=["h","period","full_source_trait_eligible"],dtype=str)
    if len(d)!=454 or len(t)!=454 or d.duplicated(["h","period"]).any():
        raise ValueError("Original 454 H7-period cohort drift")
    support=t.full_source_trait_eligible.str.lower().map({"true":True,"false":False})
    if support.isna().any() or int(support.sum())!=451:raise ValueError("Original 451 EIV source support drift")
    t["source_valid"]=support
    sub=d[d.period=="2017-2019"].merge(t[["h","period","source_valid"]],on=["h","period"],validate="one_to_one")
    ids=sorted(sub.loc[sub.source_valid,"h"].tolist())
    if len(ids)!=226:raise ValueError("Old 2017-19 H7 source quality cohort differs")
    # Outcome-blind spread across the original ordered H7 IDs, no replacement.
    ix=sorted({round(i*(len(ids)-1)/(N_PILOT-1)) for i in range(N_PILOT)})
    selected=[ids[i] for i in ix]
    raw=json.loads(Path(facets).read_text())
    if len(raw)!=2808:raise ValueError("Original GBIF facets source total changed")
    lookup={}
    for r in raw:
        if r.get("group")!="Plantae" or r.get("period")!="2017-2019":continue
        h=str(r.get("HYBAS_ID_L07"))
        if h not in selected:continue
        if h in lookup:raise ValueError("Duplicate original selected H7 WKT receipts")
        geoms=[]
        for q in r["source_request_receipts"]:
            values=urllib.parse.parse_qs(urllib.parse.urlsplit(q["url"]).query).get("geometry",[])
            if len(values)!=1:raise ValueError("Original H7 source query geometry missing")
            if hashlib.sha256(values[0].encode()).hexdigest()!=q["component_SHA"]:
                raise ValueError("Original H7 component WKT checksum differs")
            geoms.append(wkt.loads(values[0]))
        geom=unary_union(geoms)
        if geom.is_empty or geom.area<=0:raise ValueError("Invalid H7 source polygon")
        lookup[h]=geom
    if len(lookup)!=len(selected):
        raise ValueError("Original H7 GBIF WKT source not complete")
    info=sub.set_index("h")[["h6","country"]].to_dict("index")
    return [(h,lookup[h],info[h]["h6"],info[h]["country"]) for h in selected]

def geo_source(h7,geom,model,base,lid,valfield,idfield):
    bounds=geom.bounds
    spatial={
      "where":"1=1",
      "geometry":",".join(str(v) for v in bounds),
      "geometryType":"esriGeometryEnvelope",
      "inSR":4326,
      "spatialRel":"esriSpatialRelIntersects",
    }
    endpoint=f"{base}/{lid}/query"
    count,_=rq(endpoint,{"f":"json",**spatial,"returnCountOnly":"true"})
    n=count.get("count")
    if not isinstance(n,int) or n<0:raise ValueError("Invalid basin spatial count")
    if n>MAX_CANDIDATES_PER_H7_LAYER:
        return {"h7":h7,"source":model,"status":"too_many_bbox_candidates",
                "bbox_candidates":n,"coverage_area_ratio":None}
    features=[]
    receipts=[]
    for pos in range(0,n,PAGE):
        data,stamp=rq(endpoint,{"f":"geojson",**spatial,
            "outFields":",".join(("OBJECTID",idfield,valfield)),
            "returnGeometry":"true","outSR":4326,
            "orderByFields":"OBJECTID ASC","resultOffset":pos,
            "resultRecordCount":min(PAGE,n-pos)})
        part=data.get("features",[])
        if len(part)!=min(PAGE,n-pos):
            raise ValueError("Incomplete ArcGIS GeoJSON bbox page for "+model)
        receipts.append(stamp)
        features.extend(part)
    if len(features)!=n:raise ValueError("ArcGIS spatial query incomplete")
    original_coverage=[]
    numeric_cover=[]
    mass=[]
    found=set()
    pgeom=transform(PROJECT,geom)
    if not pgeom.is_valid:raise ValueError("Original H7 WKT invalid")
    for feature in features:
        attrs=feature.get("properties") or {}
        oid=attrs.get("OBJECTID")
        native=attrs.get(idfield)
        if oid is None or native is None or oid in found:
            raise ValueError("Native JRC feature OID absent or repeated")
        found.add(oid)
        g=shape(feature["geometry"])
        if g.is_empty or not g.is_valid:raise ValueError("Original JRC native polygon invalid")
        ge=transform(PROJECT,g)
        overlap=pgeom.intersection(ge)
        if overlap.is_empty:continue
        area=overlap.area
        if area<=0:continue
        original_coverage.append(overlap)
        try:
            z=float(attrs.get(valfield))
            if not math.isfinite(z) or z<0:z=None
        except (TypeError,ValueError):z=None
        if z is not None:
            numeric_cover.append(overlap)
            if model!="pointN":mass.append((area,z))
    union=unary_union(original_coverage) if original_coverage else None
    numunion=unary_union(numeric_cover) if numeric_cover else None
    cov=union.area/pgeom.area if union else 0
    numeric=numunion.area/pgeom.area if numunion else 0
    if cov>1.02 or numeric>1.02:
        raise ValueError("JRC geographic overlay larger than original H7. Check projections")
    weighted=sum(a*v for a,v in mass)/sum(a for a,_ in mass) if mass else None
    return {"h7":h7,"source":model,"status":"complete_source_bbox_query",
            "bbox_candidates":n,"true_intersecting_native_catchment_features":len(original_coverage),
            "native_H7_overlap_area_pct":100*cov,
            "native_H7_valid_model_value_overlap_pct":100*numeric,
            "soil_H7_area_weighted_native_catchment_model_mean_NOT_original_raster":weighted,
            "original_query_page_SHA256_of_pages":hashlib.sha256("\n".join(receipts).encode()).hexdigest()}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--facets",type=Path,required=True)
    p.add_argument("--eligibility",type=Path,required=True)
    p.add_argument("--panel",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    source=original_h7(args.facets,args.eligibility,args.panel)
    result=[]
    for h,geom,h6,country in source:
        for model,(base,lid,field,namespace) in SOURCES.items():
            record=geo_source(h,geom,model,base,lid,field,namespace)
            record["original_true_H6"]=h6
            record["original_source_country_label"]=country
            result.append(record)
            print("NATIVE_H7_SOURCE_PILOT",h,model,record.get("bbox_candidates"),record.get("native_H7_valid_model_value_overlap_pct"),flush=True)
    status=Counter(x["status"] for x in result)
    data={
        "schema":"gesui_20161119_JRC_native_basin_vs_H7_original_WKT_outcome_blind_12pilot_v1",
        "utc":datetime.now(timezone.utc).isoformat(),
        "selected_H7_n":len(source),"selected_H7":sorted(x[0] for x in source),
        "real_native_H7_polygon_from_original_GBIF_query":True,
        "source_H7_facets_original_SHA":hashpath(args.facets),
        "source_EIV_eligibility_SHA":hashpath(args.eligibility),
        "original_454_H7_panel_SHA":hashpath(args.panel),
        "status_counts":dict(status),
        "source_only_results":result,
        "DO_NOT_USE_PILOT_FOR_BIOTIC_REGRESSION_OR_NEW_THRESHOLD":True,
        "no_EEQ_R_EIV_NUTRI_HNS_OUTCOME_VALUES_OPENED":True,
        "native_soil_HydroID_and_point_BasinID_not_equated":True,
        "source_soil_model_basin_means_NOT_original_EMEP_raster_pixels":True,
        "original_2018_water_land_2x2_source_STOP_still_binding":True
    }
    (args.out/"JRC_NATIVE_BASIN_12_H7_GEOMETRY_PILOT_SOURCE_ONLY.json").write_text(
        json.dumps(data,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    with (args.out/"JRC_NATIVE_BASIN_12_H7_SOURCE_ONLY_COVERAGE.csv").open("w",newline="") as f:
        cols=list(result[0])
        w=csv.DictWriter(f,fieldnames=cols)
        w.writeheader();w.writerows(result)
    print("GESUI_JRC_H7_SOURCE_GEOMETRY_PILOT_DONE",json.dumps({
        "n_H7":len(source),"status_counts":dict(status),
        "atmos_covered_gt90":sum(x["source"]=="atmos2015" and
            x.get("native_H7_valid_model_value_overlap_pct") is not None and
            x["native_H7_valid_model_value_overlap_pct"]>=90 for x in result),
        "all3_soil_source_cover_gt90":sum(all(any(x["h7"]==h and x["source"]==v and
            x.get("native_H7_valid_model_value_overlap_pct") is not None and
            x["native_H7_valid_model_value_overlap_pct"]>=90 for x in result)
            for v in ("atmos2015","mineralN","organicN")) for h,_,_,_ in source),
        "NOT_A_FINAL_SOURCE_GATE":True
    }),flush=True)

if __name__=="__main__":main()
