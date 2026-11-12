#!/usr/bin/env python3
"""Stage 0, BEFORE NEW GBIF count calls: source-pinned original HydroBASINS
Level7 1219 native river polygon WKT manifest.

GBIF accepted polygons require WKT POLYGON, not MULTIPOLYGON. Preserve every
Polygon component and original holes. Any unqueryable polygon is SOURCE MISSING.
Never use harbor-centered point radius as a subbasin substitute.
"""
import collections,hashlib,io,json,os,pathlib,time,urllib.request,zipfile
import shapefile
from shapely.geometry import shape
from shapely.geometry.polygon import orient
from shapely import wkt

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_hybas7_original_polygons"
ROOT.mkdir(parents=True,exist_ok=True)
NATIVE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
sourcepath=pathlib.Path(os.environ["HYBAS7_FROZEN_NATIVE_SOURCE_MANIFEST"])
buf=sourcepath.read_bytes()
if hashlib.sha256(buf).hexdigest()!=NATIVE_SHA:
    raise RuntimeError("Original 2018/22 EPA/EEA basin source sample changed")
original=json.loads(buf)
samples=original["candidate_source_site_period_rows"]
allids={str(s["HYBAS_ID_L07"]) for s in samples if s["HYBAS_ID_L07"] is not None}
if len(allids)!=1219:raise RuntimeError("Original river monitored-basin Level7 frame not 1219")
original_source_sha={"eu":"b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44",
                     "na":"4e5150f0fe3a2c3d17b21938ec5c6ca1959d136497f572a0203af6b1842e501d"}
SHA=hashlib.sha256
def fetch(region):
    url=f"https://data.hydrosheds.org/file/hydrobasins/standard/hybas_{region}_lev07_v1c.zip"
    err=None
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-original-L7-GBIF-geometry-precounts-source/1","Accept":"application/zip,*/*"})
            with urllib.request.urlopen(req,timeout=80) as response:data=response.read(85000000)
            if SHA(data).hexdigest()!=original_source_sha[region]:
                raise RuntimeError("Official HydroBASINS original ZIP SHA changed")
            return data,url
        except Exception as e:
            err=e;time.sleep(i+1)
    raise RuntimeError("Official original HydroBASINS ZIP unavailable: "+str(err))
# GBIF requires POLYGON WKT; record every component, do not drop holes/islands.
TOLS=[0,.002,.004,.008,.015,.03,.06]
def encode_poly(p):
    if p.is_empty or not p.is_valid or p.geom_type!="Polygon":return None
    org=p.area
    if org<=0:return None
    for tol in TOLS:
        simplified=p if tol==0 else p.simplify(tol,preserve_topology=True)
        if simplified.is_empty or simplified.geom_type!="Polygon" or not simplified.is_valid:
            continue
        ratio=simplified.area/org
        if not .95<=ratio<=1.05:continue
        wk=wkt.dumps(orient(simplified,sign=1),rounding_precision=5,trim=True)
        if len(wk)<=3500:
            return {"wkt":wk,"wkt_SHA256":SHA(wk.encode()).hexdigest(),
              "original_degree_area":org,"simplified_degree_area_ratio":ratio,
              "simplification_tolerance_degrees":tol,"wkt_chars":len(wk),
              "holes_retained":len(simplified.interiors),
              "n_exterior_vertices":len(simplified.exterior.coords)}
    return None
basins={}
source_metrics={}
for region in ("eu","na"):
    data,url=fetch(region)
    z=zipfile.ZipFile(io.BytesIO(data))
    files=z.namelist()
    candidates={}
    for path in files:
        if "lev07" in path.lower() and path.rsplit(".",1)[-1].lower() in {"shp","shx","dbf"}:
            base,ext=path.rsplit(".",1)
            candidates.setdefault(base,{})[ext.lower()]=path
    found=[v for v in candidates.values() if set(v)=={"shp","shx","dbf"}]
    if len(found)!=1:raise RuntimeError("Original level-7 official shape zip source ambiguous")
    base=found[0]
    reader=shapefile.Reader(
        shp=io.BytesIO(z.read(base["shp"])),
        shx=io.BytesIO(z.read(base["shx"])),
        dbf=io.BytesIO(z.read(base["dbf"])))
    cols=[f[0] for f in reader.fields[1:]]
    if "HYBAS_ID" not in cols or "SUB_AREA" not in cols:raise RuntimeError("HydroBASINS original columns changed")
    n=0;hits=0;invalid=0
    for sr in reader.iterShapeRecords():
        n+=1
        fields=dict(zip(cols,sr.record))
        hid=str(fields["HYBAS_ID"])
        if hid not in allids:continue
        hits+=1
        g=shape(sr.shape.__geo_interface__)
        if not g.is_valid or g.is_empty:
            invalid+=1
            basins[hid]={"HYBAS_ID_L07":hid,"region":region,
                 "status":"invalid_original_geometry_not_repaired"}
            continue
        if g.geom_type=="Polygon":components=[g]
        elif g.geom_type=="MultiPolygon":components=list(g.geoms)
        else:
            basins[hid]={"HYBAS_ID_L07":hid,"region":region,"status":"unsupported_geometry_type_"+g.geom_type}
            continue
        polygons=[encode_poly(c) for c in components]
        if any(x is None for x in polygons):
            basins[hid]={"HYBAS_ID_L07":hid,"region":region,
                 "status":"some_original_components_not_queryable","components":len(components)}
        else:
            basins[hid]={"HYBAS_ID_L07":hid,"region":region,
                 "status":"all_original_polygons_WKT_queryable",
                 "sub_area_km2":float(fields["SUB_AREA"]),
                 "up_area_km2":float(fields["UP_AREA"]),
                 "components":len(components),
                 "polygon_parts":polygons,
                 "original_GEOM_area_square_degrees":g.area}
    source_metrics[region]={"official_zip_url":url,"source_zip_SHA256":original_source_sha[region],
                "source_feature_count":n,"original_river_monitored_HYBAS_ID_rows":hits,
                "original_invalid_selected_geometries":invalid}
if set(basins)!=allids:
    raise RuntimeError("Some original target Level7 basin source IDs missing from official EU/NA files")
eligible=[b for b in basins.values() if b["status"]=="all_original_polygons_WKT_queryable"]
N=len(eligible)
sizes=[p["wkt_chars"] for b in eligible for p in b["polygon_parts"]]
source_gate=(N>=.85*1219)
# Exact deterministic source 4-shard decision, independent of any GBIF result.
shards={str(i):0 for i in range(4)}
for b in eligible:
    shard=int(SHA(b["HYBAS_ID_L07"].encode()).hexdigest()[:8],16)%4
    b["original_source_shard"]=shard
    shards[str(shard)]+=1
out={
  "schema":"HYBAS7_1219_original_WKT_components_PRE_GBIF_OUTCOME_v1",
  "source_geometry_original_2018_2022_EEA_EPA_manifest_SHA256":NATIVE_SHA,
  "official_ZIP_source_receipts":source_metrics,
  "original_complete_river_monitored_L7_basin_count":1219,
  "original_WKT_eligible_basins":N,
  "original_WKT_ineligible_basins":1219-N,
  "original_WKT_full_polygon_components":sum(b["components"] for b in eligible),
  "polygons_holes_retained_count":sum(p["holes_retained"] for b in eligible for p in b["polygon_parts"]),
  "polygon_count_with_simplification":sum(p["simplification_tolerance_degrees"]>0 for b in eligible for p in b["polygon_parts"]),
  "max_WKT_characters":max(sizes) if sizes else 0,
  "n_selected_SOURCE_shards":shards,
  "predeclared_min_basin_source_geometry_coverage":.85,
  "source_gate_pass":source_gate,
  "GBIF_biological_response_counts_NOT_YET_QUERIED":True,
  "no_H3_or_buffer_substitution":True,
  "unqueryable_basin_source_ids":[b["HYBAS_ID_L07"] for b in basins.values() if b["status"]!="all_original_polygons_WKT_queryable"]
}
archive={"schema":out["schema"],
         "source_original_SHA256":NATIVE_SHA,
         "source_identities":original["source_HydroBASINS_region_SHA256"],
         "basins":sorted(basins.values(),key=lambda b:b["HYBAS_ID_L07"]),
         "geography_source_gate_pass":source_gate}
manifest=(json.dumps(archive,sort_keys=True,separators=(",",":"),ensure_ascii=False)+"\n").encode()
out["original_WKT_source_manifest_SHA256"]=SHA(manifest).hexdigest()
(ROOT/"HYBAS7_1219_ORIGINAL_NATIVE_POLYGONS_WKT_BEFORE_NEW_GBIF_COUNTS.json").write_bytes(manifest)
(ROOT/"HYBAS7_1219_POLYGON_NATIVE_SOURCE_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("HYBAS7_ORIGINAL_POLYGON_BEFORE_GBIF_RESULT_START",flush=True)
print(json.dumps(out,sort_keys=True),flush=True)
print("HYBAS7_ORIGINAL_POLYGON_BEFORE_GBIF_RESULT_END",flush=True)
if not source_gate:
    raise RuntimeError("Predeclared 85% original polygon source feasibility gate FAILED; no GBIF new counts")
