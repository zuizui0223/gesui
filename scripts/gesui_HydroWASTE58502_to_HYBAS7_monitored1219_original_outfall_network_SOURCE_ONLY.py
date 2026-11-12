#!/usr/bin/env python3
"""Official source-only HydroWASTE 58502 static discharge crosswalk to original
HydroBASINS L7 EU+North America and 1219 original river-monitoring subbasins.

One row per H7 basin; report local registered outfall and cumulative to L7
basin OUTLET. This is NOT dose at any original EEA or EPA monitoring station,
nor actual WWTP nitrogen/phosphorus 2018/2022, nor actual adult emergence.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,time,urllib.request,zipfile
import shapefile
from shapely.geometry import shape as polyshape,Point
from shapely.strtree import STRtree

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HYDROWASTE_HYBAS7_source_crosswalk"
ROOT.mkdir(parents=True,exist_ok=True)
WASTE_SOURCE_ZIP_SHA="6e19d8f31d95d533b2b29056b823eeb915d933fa22ea090f234dc38f6427d044"
WASTE_RAW_SHA="af040d9b4fdc93bab2eb8d1b52985f5bacc6b0db9f16795c29e31c10fb615dd4"
H7_MONITOR_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
REGION_SHA={
  "eu":"b9b810f4db271f0fa6b39fecf01784041b948944b3c3f296a393766388e86b44",
  "na":"4e5150f0fe3a2c3d17b21938ec5c6ca1959d136497f572a0203af6b1842e501d"
}
h7_path=pathlib.Path(os.environ["HYBAS7_FROZEN_NATIVE_SOURCE_MANIFEST"])
if hashlib.sha256(h7_path.read_bytes()).hexdigest()!=H7_MONITOR_SHA:
    raise RuntimeError("Original 1219 source river-basin candidates changed")
source=json.loads(h7_path.read_text())
wanted={str(r["HYBAS_ID_L07"]) for r in source["candidate_source_site_period_rows"] if r["HYBAS_ID_L07"] is not None}
if len(wanted)!=1219:raise RuntimeError("Original H7 river monitoring basin source sample changed")
wastepath=pathlib.Path(os.environ["HYDROWASTE_ORIGINAL_PUBLISHER_ZIP"])
archive=wastepath.read_bytes()
if hashlib.sha256(archive).hexdigest()!=WASTE_SOURCE_ZIP_SHA:raise RuntimeError("Publisher HydroWASTE ZIP differs from original source")
z=zipfile.ZipFile(io.BytesIO(archive))
data=z.read("HydroWASTE_v10.csv")
if hashlib.sha256(data).hexdigest()!=WASTE_RAW_SHA:raise RuntimeError("Original HydroWASTE CSV sha differs")
readme=z.read("README.txt").decode("iso-8859-1")
p=readme.find("WASTE_DIS")
print("Original publisher WASTE_DIS unit raw README context",readme[max(0,p-200):p+370].replace("\n"," "),flush=True)
print("Original publisher RIVER_DIS unit raw README context",readme[max(0,readme.find("RIVER_DIS")-150):readme.find("RIVER_DIS")+280].replace("\n"," "),flush=True)
reader=csv.DictReader(io.StringIO(data.decode("iso-8859-1"),newline=""))
heads=reader.fieldnames or []
REQUIRED={"WASTE_ID","LAT_OUT","LON_OUT","HYRIV_ID","WASTE_DIS","QUAL_WASTE","STATUS","LEVEL","RIVER_DIS"}
if not REQUIRED.issubset(set(heads)):raise RuntimeError("Original source unit/readme columns changed")
plants=list(reader)
if len(plants)!=58502 or len({r["WASTE_ID"] for r in plants})!=58502:
    raise RuntimeError("Original global HydroWASTE plant ID census changed")
# Read ALL EU/NA L7 subbasins; avoid assigning only target polygons because
# upstream sources may originate in other local subbasins.
allpolys=[]
allid=[]
areas=[]
up_areas=[]
NEXT={}
shape_audit={}
for region in ("eu","na"):
    url=f"https://data.hydrosheds.org/file/hydrobasins/standard/hybas_{region}_lev07_v1c.zip"
    req=urllib.request.Request(url,headers={"Accept":"application/zip,*/*","User-Agent":"gesui-HydroWASTE-true-L7-outfall-route-source/1"})
    with urllib.request.urlopen(req,timeout=85) as r:
        bb=r.read(85000000)
    if hashlib.sha256(bb).hexdigest()!=REGION_SHA[region]:
        raise RuntimeError("Official 2018/2022 HydroBASINS source ZIP unexpectedly modified")
    zs=zipfile.ZipFile(io.BytesIO(bb))
    paths=zs.namelist();groups={}
    for f in paths:
        if "lev07" in f.lower() and f.rsplit(".",1)[-1].lower() in {"shp","shx","dbf"}:
            a,ext=f.rsplit(".",1)
            groups.setdefault(a,{})[ext.lower()]=f
    opts=[x for x in groups.values() if set(x)=={"shp","shx","dbf"}]
    if len(opts)!=1:raise RuntimeError("Official H7 ZIP shapefile not uniquely available")
    q=opts[0]
    shp=shapefile.Reader(shp=io.BytesIO(zs.read(q["shp"])),shx=io.BytesIO(zs.read(q["shx"])),
                           dbf=io.BytesIO(zs.read(q["dbf"])))
    columns=[z[0] for z in shp.fields[1:]]
    if not {"HYBAS_ID","NEXT_DOWN","UP_AREA","SUB_AREA"}.issubset(set(columns)):
        raise RuntimeError("Missing genuine H7 flow-graph topology")
    n_original=0;invalid=0
    for rec in shp.iterShapeRecords():
        n_original+=1
        d=dict(zip(columns,rec.record))
        hid=str(d["HYBAS_ID"])
        if hid in NEXT:raise RuntimeError("Duplicate full source basin ID across regions")
        nxt=str(d["NEXT_DOWN"])
        NEXT[hid]=nxt
        geom=polyshape(rec.shape.__geo_interface__)
        if not geom.is_valid or geom.is_empty:
            invalid+=1
            continue
        allpolys.append(geom)
        allid.append(hid)
        areas.append(float(d["SUB_AREA"]))
        up_areas.append(float(d["UP_AREA"]))
    shape_audit[region]={"n_source_raw_polygons":n_original,"n_invalid_original_geometry_excluded":invalid,
      "n_valid_polygons":n_original-invalid,"official_hybas7_zip_SHA256":REGION_SHA[region]}
if not wanted.issubset(set(allid)):
    raise RuntimeError("Some target 1219 monitored basin sources missing after original shape integrity")
tree=STRtree(allpolys)
print("Original true HydroBASINS polygons source total",len(allid),"and original HydroWASTE WWTP",len(plants),flush=True)
source_counts=collections.Counter()
status_all=collections.Counter()
localcount=collections.Counter();localflow=collections.Counter()
quality_by_basin=collections.defaultdict(collections.Counter)
hits=collections.defaultdict(list)
missing_geom=0;not_in_eu_na=0;ambiguous=0
def number(s):
    try:
        a=float(s)
        return a if math.isfinite(a) else None
    except (ValueError,TypeError):return None
for plant in plants:
    status_all[plant["STATUS"]]+=1
    lat=number(plant["LAT_OUT"]);lon=number(plant["LON_OUT"])
    if lat is None or lon is None or not -90<=lat<=90 or not -180<=lon<=180:
        missing_geom+=1;continue
    pt=Point(lon,lat)
    ix=tree.query(pt,predicate="intersects")
    if len(ix)==0:
        not_in_eu_na+=1;continue
    inside=[int(i) for i in ix if allpolys[int(i)].contains(pt)]
    if len(inside)==1:
        idx=inside[0]
    elif len(ix)==1:
        idx=int(ix[0])
    else:
        ambiguous+=1;continue
    hid=allid[idx]
    v=number(plant["WASTE_DIS"])
    if v is None:
        source_counts["missing_float_waste_dis"]+=1
    elif v<0:
        source_counts["negative_waste_dis"]+=1
    localcount[hid]+=1
    if v is not None and v>=0:
        localflow[hid]+=v
        source_counts["mapped_outfall_numeric_nonnegative_flow"]+=1
    quality_by_basin[hid][plant["QUAL_WASTE"]]+=1
    source_counts["mapped_outfall_to_unique_EU_NA_HYBAS7"]+=1
    if hid in wanted:
        hits[hid].append(plant["WASTE_ID"])
source_counts.update({
 "original_58502_WWTP":len(plants),
 "original_EU_NA_basin_features":len(NEXT),
 "original_EU_NA_valid_polygons":len(allpolys),
 "non_EU_NA_or_unassigned_outfalls":not_in_eu_na,
 "ambiguous_catchment_boundary_outfalls":ambiguous,
 "invalid_outfall_point":missing_geom
})
# Directed drainage graph on entire EU+NA Level7 (not all globe), a
# cumulative load reaching the BASIN OUTLET. NEVER a source attribution to a
# monitoring station inside the basin. Upstream outside included source extent
# or invalid polygons may be missing; record separate source coverage warning.
upstream_nodes=collections.defaultdict(list)
indeg=collections.Counter({k:0 for k in NEXT})
missing_links=0
for hid,next_id in NEXT.items():
    if next_id in NEXT and next_id!=hid and next_id not in {"0","-1"}:
        upstream_nodes[next_id].append(hid)
        indeg[next_id]+=1
    elif next_id not in ("0","-1",hid):
        missing_links+=1
queue=collections.deque([k for k in NEXT if indeg[k]==0])
cumulativecount=collections.Counter(localcount)
cumulativeflow=collections.Counter(localflow)
traversed=0
while queue:
    b=queue.popleft();traversed+=1
    dst=NEXT[b]
    if dst in NEXT and dst!=b:
        cumulativecount[dst]+=cumulativecount[b]
        cumulativeflow[dst]+=cumulativeflow[b]
        indeg[dst]-=1
        if indeg[dst]==0:queue.append(dst)
if traversed!=len(NEXT):
    raise RuntimeError("Source graph next-down topology cyclic: not all original H7 basins processed")
out=[]
pointindex={allid[i]:i for i in range(len(allid))}
for hid in sorted(wanted):
    idx=pointindex[hid]
    out.append({
     "HYBAS_ID_L07":hid,
     "LOCAL_SUB_AREA_KM2":areas[idx],
     "UPSTREAM_AREA_KM2":up_areas[idx],
     "HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT":localcount[hid],
     "HYDROWASTE_STATIC_LOCAL_MODELLED_WASTEWATER_DISCHARGE_M3_PER_DAY":round(localflow[hid],6),
     "HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT":cumulativecount[hid],
     "HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY":round(cumulativeflow[hid],6),
     "ORIGINAL_QUAL_WASTE_LOCAL_LEVELS":dict(quality_by_basin[hid]),
     "UPSTREAM_BASIN_OUTLET_NOT_RECEIVER_STATION_DOSE":True,
     "NO_MEASURED_2018_2022_NP_TPY_OR_SOURCE_YEAR":True
    })
if len(out)!=1219 or any(z["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT"]<
    z["HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT"] for z in out):
    raise RuntimeError("Expected 1219 basins or flow-routing monotonicity failed")
localpos=sum(z["HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT"]>0 for z in out)
up_pos=sum(z["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT"]>0 for z in out)
sourcegate=bool(len(allpolys)>13000 and source_counts["mapped_outfall_to_unique_EU_NA_HYBAS7"]>=4000
  and localpos>50 and up_pos>localpos and traversed==len(NEXT))
summary={
 "schema":"GESUI_TRUE_GLOBAL_HYDROWASTE58K_EU_NA_MONITORED_1219_HYBAS7_STATIC_SOURCE_ONLY_v1",
 "original_HYDROWASTE_ZIP_SHA256":WASTE_SOURCE_ZIP_SHA,
 "original_HYDROWASTE_CSV_SHA256":WASTE_RAW_SHA,
 "original_monitoring_HYBAS7_source_SHA256":H7_MONITOR_SHA,
 "original_hydrological_geometries":shape_audit,
 "original_outfall_all_countries_counts":dict(status_all),
 "original_source_QA":dict(source_counts),
 "original_H7_outlet_topology_all_basins":len(NEXT),
 "original_hydrological_topology_edges_outside_regions":missing_links,
 "topologically_sorted_L7_nodes":traversed,
 "monitored_L7_basins":len(out),
 "monitored_L7_basins_at_least_1_LOCAL_original_WWTP":localpos,
 "monitored_L7_basins_at_least_1_UPSTREAM_TO_BASIN_OUTLET_WWTP":up_pos,
 "n_unique_regional_plants_with_LOCAL_outfall_in_target":sum(z["HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT"] for z in out),
 "source_gate_pass":sourcegate,
 "no_true_station_specific_upstream_delivered_flow":True,
 "no_EEA_2018_2022_true_NP_link_or_biological_effect_estimations":True,
 "static_modelled_estimated_WASTE_DIS_NOT_annual_2018_2022_sewage_load":True
}
csvpath=ROOT/"HYBAS7_1219_OFFICIAL_HYDROWASTE_STATIC_UPSTREAM_OUTFALL_PRESSURE_SOURCE_ONLY.csv"
with csvpath.open("w",newline="") as f:
    fields=[k for k in out[0] if k!="ORIGINAL_QUAL_WASTE_LOCAL_LEVELS"]
    writer=csv.DictWriter(f,fieldnames=fields)
    writer.writeheader()
    writer.writerows([{k:v for k,v in r.items() if k in fields} for r in out])
summary["original_1219_source_csv_SHA256"]=hashlib.sha256(csvpath.read_bytes()).hexdigest()
(ROOT/"HYDROWASTE_EU_NA_1219_HYBAS7_OUTFALL_SOURCE_AUDIT.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
print("HYDROWASTE_1219_WWTP_RIVER_SOURCE_JOIN_RESULT_START",flush=True)
print(json.dumps(summary,sort_keys=True),flush=True)
print("HYDROWASTE_1219_WWTP_RIVER_SOURCE_JOIN_RESULT_END",flush=True)
if not sourcegate:raise RuntimeError("Predefined original WWTP static source crosswalk density gate failed; no effects")
