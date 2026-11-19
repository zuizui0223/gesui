#!/usr/bin/env python3
"""Original-cohort 226-H7 SOURCE-only JRC native-basin geometry census, sharded.

The 2018 water x 2017-19 land 2x2 remains SOURCE-STOP. No EQR, plant NUTRI
or butterfly HNS outcome values are loaded. Official H7 query WKT and
quality-eligibility metadata are SHA-locked. Namespaced native basin IDs must
never be compared by literal equality between independent JRC services.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, time, urllib.parse, urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from shapely import wkt, make_valid
from shapely.geometry import Polygon
from shapely.ops import unary_union, transform
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
EXPECTED_H7=226
SHARDS=12
PAGE=130
MAX_BBOX_CANDIDATES=1600
P=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform

def filehash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for data in iter(lambda:f.read(1024*1024),b""):h.update(data)
    return h.hexdigest()

def req(endpoint,params):
    url=endpoint+"?"+urllib.parse.urlencode(params)
    last=None
    for attempt in range(4):
        try:
            q=urllib.request.Request(url,headers={
              "User-Agent":"gesui-original-226H7-native-JRC-SOURCE-only-20261119",
              "Accept":"application/json",
            })
            with urllib.request.urlopen(q,timeout=65) as f:raw=f.read(35_000_000)
            data=json.loads(raw)
            if not isinstance(data,dict) or "error" in data:
                raise ValueError("Official JRC ArcGIS error "+str(data.get("error")))
            return data,hashlib.sha256(raw).hexdigest()
        except (OSError,ValueError,TimeoutError) as exc:
            last=exc
            if attempt<3:time.sleep(attempt+1)
    raise RuntimeError("Original JRC source API unavailable: "+str(last))

def load_roster(elig,panel):
    if filehash(panel)!=PANEL_SHA:raise ValueError("Immutable 454 H7 biology-containing panel source changed (values never read)")
    p=pd.read_csv(panel,usecols=["h","h6","country","period"],dtype=str)
    q=pd.read_csv(elig,usecols=["h","period","full_source_trait_eligible"],dtype=str)
    if len(p)!=454 or len(q)!=454 or p.duplicated(["h","period"]).any() or q.duplicated(["h","period"]).any():
        raise ValueError("Old 454 H7-period source inventory changed")
    q["quality"]=q.full_source_trait_eligible.str.lower().map({"true":True,"false":False})
    if q.quality.isna().any() or int(q.quality.sum())!=451:
        raise ValueError("Original 451/454 EIV SOURCE quality labels differ")
    first=p[p.period.eq("2017-2019")].merge(q[["h","period","quality"]],on=["h","period"],validate="one_to_one")
    if len(first)!=227:raise ValueError("Original first-period 227 H7 changed")
    good=first.loc[first.quality].copy()
    if len(good)!=EXPECTED_H7 or good.h.duplicated().any():
        raise ValueError("Original quality-supported first-period 226 H7 changed")
    return good.set_index("h").loc[sorted(good.h)].reset_index()

def load_shard_h7(facets,roster,shard):
    if filehash(facets)!=FACETS_SHA:raise ValueError("Original 2808-facet WKT SHA changed")
    selected={h for i,h in enumerate(roster.h) if i%SHARDS==shard}
    data=json.loads(Path(facets).read_text())
    if len(data)!=2808:raise ValueError("Original source GBIF species-facet roster changed")
    found={}
    for r in data:
        if r.get("group")!="Plantae" or r.get("period")!="2017-2019":continue
        h=str(r.get("HYBAS_ID_L07"))
        if h not in selected:continue
        if h in found:raise ValueError("Duplicate H7 source polygon query facets")
        fragments=[]
        for z in r["source_request_receipts"]:
            geom=urllib.parse.parse_qs(urllib.parse.urlsplit(z["url"]).query).get("geometry",[])
            if len(geom)!=1 or hashlib.sha256(geom[0].encode()).hexdigest()!=z["component_SHA"]:
                raise ValueError("Original GBIF native H7 WKT SHA mismatch")
            fragments.append(wkt.loads(geom[0]))
        g=unary_union(fragments)
        if g.is_empty or not g.is_valid or g.area<=0:raise ValueError("Invalid original H7 geometry "+h)
        found[h]=g
    if set(found)!=selected:raise ValueError("One or more original H7 source polygon receipts missing")
    return found

def signed(coords):
    return .5*sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(coords,coords[1:]+coords[:1]))

def esri_shape(rings):
    outer=[]
    holes=[]
    invalid_ring_count=0
    for coords in rings:
        if len(coords)<4:continue
        geom=Polygon(coords)
        if geom.is_empty or geom.area<=0:continue
        if not geom.is_valid:
            invalid_ring_count+=1
            geom=make_valid(geom)
        if geom.is_empty:raise ValueError("Cannot repair native original Esri polygon ring")
        if signed(coords)<0:outer.append(geom)
        else:holes.append(geom)
    if not outer:raise ValueError("Native JRC EsriJSON geometry did not contain a clockwise shell")
    g=unary_union(outer)
    if holes:g=g.difference(unary_union(holes))
    if not g.is_valid:g=make_valid(g)
    if g.is_empty or not g.is_valid or g.area<=0:
        raise ValueError("Native JRC polygon invalid after original Esri ring decoding")
    return g,invalid_ring_count

def numeric(x):
    try:v=float(x)
    except (TypeError,ValueError):return None
    return v if math.isfinite(v) and v>=0 else None

def h7_layer(h,original,projected,model,base,lid,field,idfield):
    box=original.bounds
    spatial={
      "where":"1=1","geometry":",".join(str(x) for x in box),
      "geometryType":"esriGeometryEnvelope","inSR":4326,
      "spatialRel":"esriSpatialRelIntersects"
    }
    endpoint=f"{base}/{lid}/query"
    obj,count_sha=req(endpoint,{"f":"json",**spatial,"returnCountOnly":"true"})
    count=obj.get("count")
    if type(count) is not int or count<0:raise ValueError("Invalid native JRC regional spatial count")
    empty={
      "h7":h,"source":model,"original_bbox_count":count,
      "native_geom_coverage_pct":None,"numeric_geom_coverage_pct":None,
      "max_single_native_id":None,"max_single_native_coverage_pct":None,
      "max_numeric_single_native_coverage_pct":None,
      "soil_model_basin_weighted_mean":None,
      "model_pointN_t_per_native_area_km2_proxy":None,
      "num_orig_Esri_ring_repairs":None,
      "number_geom_features_intersecting_H7":None,
      "native_model_source_ids_intersecting_H7":[],
      "query_count_SHA256":count_sha,
      "query_page_hashes_SHA256":None,
      "geometry_area_overlap_fraction_gt1p02":False,
      "material_original_jrc_ring_repair":False,
    }
    if count>MAX_BBOX_CANDIDATES:
        return {**empty,"status":"SOURCE_STOP_BBOX_TOO_MANY_FEATURES"}
    records={}
    hashes=[]
    base_fields=["OBJECTID",idfield,field]
    # Point source PsN is a total per basin. Convert only if official basin
    # AreaSqKm is provided; never area-average tonnes as local H7 totals.
    if model=="pointN":base_fields+=["AreaSqKm"]
    for start in range(0,count,PAGE):
        reqparam={"f":"json",**spatial,
            "outFields":",".join(base_fields),
            "returnGeometry":"true","outSR":4326,
            "orderByFields":"OBJECTID ASC",
            "resultOffset":start,"resultRecordCount":min(PAGE,count-start)}
        data,stamp=req(endpoint,reqparam)
        features=data.get("features",[])
        if len(features)!=min(PAGE,count-start):
            raise ValueError("Truncated native EsriJSON page "+model+" for "+h)
        hashes.append(stamp)
        for r in features:
            a=r.get("attributes",{})
            oid=a.get("OBJECTID")
            native=a.get(idfield)
            if oid is None or native is None or oid in records:
                raise ValueError("Duplicated/invalid original native source objectID")
            geom,fixes=esri_shape(r["geometry"]["rings"])
            shape_proj=transform(P,geom)
            intersection=shape_proj.intersection(projected)
            if intersection.is_empty or intersection.area<=0:
                records[oid]=None
                continue
            area=float(intersection.area)
            source_val=numeric(a.get(field))
            native_area=numeric(a.get("AreaSqKm")) if model=="pointN" else None
            density=source_val/native_area if (source_val is not None and native_area and native_area>0) else None
            if model=="pointN" and source_val is not None and density is None:
                raise ValueError("Point source total t lacks original positive basin area")
            records[oid]={
                "source_id":str(native),"area":area,"intersection":intersection,
                "model_value":source_val,"point_intensity":density,"repaired_original_rings":fixes
            }
        time.sleep(.025)
    if len(records)!=count:raise ValueError("Original ArcGIS OIDs incomplete")
    active=[x for x in records.values() if x is not None]
    geom_union=unary_union([x["intersection"] for x in active]) if active else None
    valid=[x for x in active if x["model_value"] is not None]
    valid_union=unary_union([x["intersection"] for x in valid]) if valid else None
    parea=projected.area
    geocover=100*(geom_union.area/parea if geom_union else 0)
    numcover=100*(valid_union.area/parea if valid_union else 0)
    if geocover>102 or numcover>102:raise ValueError("GIS overlap area exceeded H7 by >2%")
    # A 'single native basin' criterion is different from the union coverage.
    ordered=sorted(active,key=lambda x:(-x["area"],x["source_id"]))
    largest=ordered[0] if ordered else None
    uniqueIds=len({x["source_id"] for x in active})==len(active)
    if not uniqueIds:raise ValueError("Same JRC native source ID repeated inside one H7; do not pseudo-sum")
    numeric_largest=max((x["area"] for x in valid),default=0)
    numeric_measure=[(x["area"],x["model_value"]) for x in valid] if model!="pointN" else [
       (x["area"],x["point_intensity"]) for x in valid if x["point_intensity"] is not None]
    measure=sum(a*v for a,v in numeric_measure)/sum(a for a,_ in numeric_measure) if numeric_measure else None
    return {
        **empty,"status":"COMPLETE_OFFICIAL_JRC_SOURCE_ONLY",
        "native_geom_coverage_pct":geocover,"numeric_geom_coverage_pct":numcover,
        "max_single_native_id":largest["source_id"] if largest else None,
        "max_single_native_coverage_pct":100*largest["area"]/parea if largest else 0,
        "max_numeric_single_native_coverage_pct":100*numeric_largest/parea,
        "soil_model_basin_weighted_mean":measure if model!="pointN" else None,
        "model_pointN_t_per_native_area_km2_proxy":measure if model=="pointN" else None,
        "num_orig_Esri_ring_repairs":sum(x["repaired_original_rings"] for x in active),
        "number_geom_features_intersecting_H7":len(active),
        "native_model_source_ids_intersecting_H7":[x["source_id"] for x in ordered],
        "query_page_hashes_SHA256":hashlib.sha256("\n".join(hashes).encode()).hexdigest(),
        "geometry_area_overlap_fraction_gt1p02":geocover>102,
        "material_original_jrc_ring_repair":sum(x["area"] for x in active if x["repaired_original_rings"])>0.05*parea
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--facets",type=Path,required=True)
    parser.add_argument("--eligibility",type=Path,required=True)
    parser.add_argument("--panel",type=Path,required=True)
    parser.add_argument("--shard",type=int,required=True)
    parser.add_argument("--n-shards",type=int,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    if args.n_shards!=SHARDS or not 0<=args.shard<SHARDS:
        raise ValueError("Fixed total shard count and deterministic index cannot be changed")
    roster=load_roster(args.eligibility,args.panel)
    allgeom=load_shard_h7(args.facets,roster,args.shard)
    args.out.mkdir(parents=True,exist_ok=True)
    rows=[]
    for h,info in roster.iterrows():
        ident=str(info["h"])
        if ident not in allgeom:continue
        geom=allgeom[ident]
        proj=transform(P,geom)
        for model,(base,lid,field,namespace) in SOURCES.items():
            try:
                record=h7_layer(ident,geom,proj,model,base,lid,field,namespace)
            except Exception as exc:
                record={"h7":ident,"source":model,"status":"TECHNICAL_SOURCE_RETRIEVAL_STOP",
                        "technical_error":str(exc)[:800]}
            record["true_h6"]=str(info["h6"])
            record["original_source_country"]=str(info["country"])
            rows.append(record)
            print("GESUI_NATIVE_226_SOURCE",args.shard,ident,model,
                  record["status"],round(record.get("numeric_geom_coverage_pct") or 0,3),
                  round(record.get("max_numeric_single_native_coverage_pct") or 0,3),
                  flush=True)
    if len(rows)!=len(allgeom)*len(SOURCES) or len({(x["h7"],x["source"]) for x in rows})!=len(rows):
        raise ValueError("Fixed H7 x source inventory incomplete or duplicated")
    count=Counter(x["status"] for x in rows)
    result={
        "schema":"gesui_JRC_originalH7_2017to19_native_source_226_full_sharded_v1",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "shard":args.shard,"n_shards":SHARDS,
        "original_total_H7_EIV_eligible_first_period":EXPECTED_H7,
        "original_H7_original_panel_SHA256":filehash(args.panel),
        "original_2808_facets_SHA256":filehash(args.facets),
        "EIV_451_trait_quality_eligibility_SHA256":filehash(args.eligibility),
        "H7_processed":len(allgeom),"number_model_cells":len(rows),
        "source_status_counts":dict(count),
        "source_rows":rows,
        "JRC_GREEN_pointN_total_not_a_site_received_mass":True,
        "pointN_density_proxy_derived_using_original_JRC_native_basin_area":True,
        "native_GIS_basis":"JRC EsriJSON original polygon rings EPSG4326 over original source GBIF H7 WKT EPSG4326",
        "no_EQR_NUTRI_HNS_biological_outcome_numeric_read":True,
        "old_2018_2017to19_NUTS2_biotic_2x2_SOURCE_STOP_still_binding":True
    }
    filepath=args.out/f"JRC_NATIVE_SOURCE_H7_SHARD_{args.shard:02d}_OF_12.json"
    filepath.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_NATIVE_FULL_SOURCE_SHARD_RESULT",json.dumps({
        "shard":args.shard,"H7":len(allgeom),
        "source_status_counts":dict(count),
        "strict_coeligible_h7_in_shard":sum(
            all(any(z["h7"]==h and z["source"]==m and
                    z["status"]=="COMPLETE_OFFICIAL_JRC_SOURCE_ONLY" and
                    (z.get("max_numeric_single_native_coverage_pct") or 0)>=90
                    for z in rows)
                for m in ("atmos2015","mineralN","organicN","pointN"))
            for h in allgeom),
        "no_biotic_values_read":True
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()
