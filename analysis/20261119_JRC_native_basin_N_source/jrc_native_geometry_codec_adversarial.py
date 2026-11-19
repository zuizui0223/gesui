#!/usr/bin/env python3
"""Independent GIS QA of original ArcGIS EsriJSON vs GeoJSON polygon interpretation.

Three fixed source-only H7: one unflagged and two requiring >5% repair in
the 12-H7 native JRC pilot. This diagnoses geometry codecs, NOT biological
phenotypes, and cannot qualify all 226 H7 for ecological modelling.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time,urllib.parse,urllib.request
from pathlib import Path
from datetime import datetime,timezone
from shapely.geometry import Polygon,shape
from shapely.ops import transform,unary_union
from shapely import make_valid
from pyproj import Transformer

BASIN_SOURCE="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/MIS_Inputs_to_soil/MapServer"
POINT_SOURCE="https://water.jrc.ec.europa.eu/arcgis/rest/services/inmap/Losses_to_water_MIS/MapServer"
FOCUS=["2070012580","2070021030","2070498990"]
LAYERS={
  "atmos":(BASIN_SOURCE,5,"Basin_3035_HydroID","MEAN_mgsqmy","Basin_3035_AreaSqKm"),
  "mineral":(BASIN_SOURCE,10,"Basin_3035_HydroID","MEAN_kgha","Basin_3035_AreaSqKm"),
  "pointN":(POINT_SOURCE,5,"BasinID","PsN","AreaSqKm")
}
PANELSHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PROJECT=Transformer.from_crs("EPSG:4326","EPSG:6933",always_xy=True).transform

def sha(b):return hashlib.sha256(b).hexdigest()
def http(url,args):
    last=None
    for j in range(4):
        try:
            req=urllib.request.Request(url+"?"+urllib.parse.urlencode(args),headers={"User-Agent":"gesui-native-esri-json-geometry-source-audit"})
            with urllib.request.urlopen(req,timeout=70) as f:raw=f.read(15000000)
            result=json.loads(raw)
            if "error" in result:raise ValueError(str(result["error"]))
            return result,sha(raw)
        except (TimeoutError,OSError,ValueError) as e:
            last=e
            if j<3:time.sleep(j+1)
    raise RuntimeError("JRC raw source unavailable: "+str(last))

def signed(r):
    return sum((a[0]*b[1]-b[0]*a[1]) for a,b in zip(r,r[1:]+r[:1]))/2

def from_esri(rings):
    exterior=[]
    interior=[]
    for coordinates in rings:
        if len(coordinates)<4:continue
        poly=Polygon(coordinates)
        if poly.is_empty or poly.area<=0:continue
        s=signed(coordinates)
        if s<0:exterior.append(make_valid(poly))
        else:interior.append(make_valid(poly))
    if not exterior:raise ValueError("No Esri clockwise original polygon shell; original native geometry ambiguous")
    g=unary_union(exterior)
    if interior:g=g.difference(unary_union(interior))
    if not g.is_valid:g=make_valid(g)
    if g.is_empty:raise ValueError("Empty original Esri rings geometry")
    return g

def pilot_geom(original):
    if sha(Path(original).read_bytes())!=PANELSHA:raise ValueError("Original 2808-facet SHA drift")
    raw=json.loads(Path(original).read_text())
    geom={}
    for z in raw:
        if z.get("period")!="2017-2019" or z.get("group")!="Plantae" or str(z.get("HYBAS_ID_L07")) not in FOCUS:continue
        h=str(z["HYBAS_ID_L07"])
        from shapely import wkt
        items=[]
        for r in z["source_request_receipts"]:
            query=urllib.parse.parse_qs(urllib.parse.urlsplit(r["url"]).query).get("geometry",[])
            if len(query)!=1 or sha(query[0].encode())!=r["component_SHA"]:raise ValueError("Original H7 source WKT receipt drift")
            items.append(wkt.loads(query[0]))
        geom[h]=unary_union(items)
    if set(geom)!=set(FOCUS):raise ValueError("Original focus H7 source geometries incomplete")
    return geom

def collect(h,geom,name,base,layer,idfield,valuefield,areafield):
    bbox=geom.bounds
    args={"where":"1=1","geometry":",".join(map(str,bbox)),
          "geometryType":"esriGeometryEnvelope","inSR":4326,
          "spatialRel":"esriSpatialRelIntersects",
          "returnGeometry":"true","outSR":4326,
          "outFields":",".join(["OBJECTID",idfield,valuefield,areafield]),
          "orderByFields":"OBJECTID ASC","resultRecordCount":500}
    endpoint=f"{base}/{layer}/query"
    count,_=http(endpoint,{"f":"json",**args,"returnCountOnly":"true"})
    n=count.get("count")
    if type(n) is not int or n>500 or n<0:raise ValueError("Fixed QA source bbox too many native polygons")
    sj,js=http(endpoint,{"f":"json",**args})
    gj,gs=http(endpoint,{"f":"geojson",**args})
    x=sj.get("features",[])
    y=gj.get("features",[])
    if len(x)!=n or len(y)!=n:raise ValueError("Native EsriJSON/GeoJSON incomplete")
    orig={z["attributes"]["OBJECTID"]:z for z in x}
    altered={z["properties"]["OBJECTID"]:z for z in y}
    if set(orig)!=set(altered):raise ValueError("Independent JSON geometry codecs did not return same native OIDs")
    pH7=transform(PROJECT,geom)
    areas={"make_valid":[],"buffer_zero":[],"esri_rings":[]}
    native_area_deltas=[]
    invalid=0
    for oid,z in orig.items():
        g=shape(altered[oid]["geometry"])
        invalid+=int(not g.is_valid)
        variants={
          "make_valid":make_valid(g) if not g.is_valid else g,
          "buffer_zero":g.buffer(0) if not g.is_valid else g,
          "esri_rings":from_esri(z["geometry"]["rings"])
        }
        published=z["attributes"].get(areafield)
        for typ,item in variants.items():
            if item.is_empty or not item.is_valid:raise ValueError("Failed polygon decoder "+typ)
            pa=transform(PROJECT,item)
            inter=pa.intersection(pH7)
            if inter.area>0:areas[typ].append(inter)
            if published is not None and float(published)>0 and typ=="make_valid":
                native_area_deltas.append((pa.area/1000000-float(published))/float(published))
    covered={}
    for typ,shapes in areas.items():
        union=unary_union(shapes) if shapes else None
        covered[typ]=100*(union.area/pH7.area if union else 0.)
    differences={k:covered[k]-covered["make_valid"] for k in covered}
    return {"h7":h,"source":name,"bbox_source_n":n,
            "native_invalid_GeoJSON_feature_n":invalid,
            "spatial_H7_coverage_percent_by_decoder":covered,
            "difference_from_makevalid_percentage_points":differences,
            "native_published_area_relative_differences_for_make_valid":native_area_deltas,
            "esri_original_query_SHA":js,"geojson_query_SHA":gs,
            "codec_agreement_within_1pct":all(abs(x)<1 for x in differences.values())}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-h7-facets",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    polygons=pilot_geom(a.original_h7_facets)
    rows=[]
    for h in FOCUS:
        for name,(url,layer,key,val,area) in LAYERS.items():
            z=collect(h,polygons[h],name,url,layer,key,val,area)
            rows.append(z)
            print("SOURCE_GEOM_CODEC_COMPARE",h,name,z["native_invalid_GeoJSON_feature_n"],z["spatial_H7_coverage_percent_by_decoder"],flush=True)
    out={"schema":"gesui_20261119_JRC_native_source_esri_vs_geojson_geom_quality_source_only_3H7",
         "UTC":datetime.now(timezone.utc).isoformat(),
         "fixed_3H7":[{"id":x,"meaning":"original SOURCE-only pilot affected/unaffected"} for x in FOCUS],
         "records":rows,
         "all_three_source_geometry_decoders_agree_within_1pct":all(z["codec_agreement_within_1pct"] for z in rows),
         "original_h7_facets_SHA":sha(a.original_h7_facets.read_bytes()),
         "no_EQR_EIV_NUTRI_HNS_VALUES_OPENED":True,
         "NOT_A_FINAL_JRC_SOURCE_OR_BIOTIC_GATE":True,
         "original_2018_water_land_2x2_SOURCE_STOP_remains":True}
    (a.out/"JRC_ESRIJSON_VS_GEOJSON_NATIVE_GEOMETRY_QA.json").write_text(
       json.dumps(out,indent=2,sort_keys=True)+"\n")
    print("GESUI_JRC_NATIVE_GEO_CODEC_PILOT",json.dumps({
       "native_qa_H7":len(FOCUS),"total_source_comparisons":len(rows),
       "all_decoders_area_agree_within_1pct":out["all_three_source_geometry_decoders_agree_within_1pct"]
    }),flush=True)
if __name__=="__main__":main()
