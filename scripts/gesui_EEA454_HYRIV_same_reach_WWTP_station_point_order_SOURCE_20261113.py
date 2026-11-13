#!/usr/bin/env python3
"""SOURCE ONLY original 2018/2022 EEA site vs original WWTP on SAME original
HYRIV reach (Europe). Uses original HydroRIVERS full SHP vector polyline
geometry and NEXT_DOWN neighbor endpoint to infer local reach orientation,
rather than assuming stored polyline coordinate ordering.

No EEA EQR/GBIF biodiversity magnitude accessed. Original pre-declared
<=500m point-to-line and endpoint-join gates, >=500m source-site along-line
gap and >=2*(station+outfall snaps) are hard input QA. Ambiguous source stays
unclassified, not 0/positive intervention. Native WASTE_DIS m3/day is static.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,statistics,time,urllib.request,zipfile
import shapefile
from pyproj import CRS,Transformer
from shapely.geometry import Point,LineString

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_EU454_HYRIV_same_reach_geometry"
ROOT.mkdir(parents=True,exist_ok=True)
ZIP_URL="https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip"
SOURCE_RIVER_SHA="500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5"
ORIGINAL_WASTE_SHA="6e19d8f31d95d533b2b29056b823eeb915d933fa22ea090f234dc38f6427d044"
ORIGINAL_WASTE_CSV_SHA="af040d9b4fdc93bab2eb8d1b52985f5bacc6b0db9f16795c29e31c10fb615dd4"
ORIGINAL_H7_SITE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
EEA={
 "EEA18":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
def gitsha(raw):return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
def reach(x):
    try:
        n=int(float(x))
        return str(n) if n>0 else None
    except (ValueError,TypeError):return None
def num(x):
    try:
        n=float(x)
        return n if math.isfinite(n) else None
    except (ValueError,TypeError):return None

# Cohort identities only. Responses in source CSV MUST NOT be used for source
# geometry eligibility or point ordering.
with open(os.environ["ORIGINAL_454_H7_STUDY_COHORT"],encoding="utf-8-sig",newline="") as f:
    fixed={(q["HYBAS_ID_L07"],q["PERIOD"]) for q in csv.DictReader(f)}
if len(fixed)!=454 or len({q[0] for q in fixed})!=232:
    raise RuntimeError("Original 454 basinperiod n232 H7 cohort changed")
nativepath=pathlib.Path(os.environ["ORIGINAL_EEA_STATION_H7_SOURCE"])
if hashlib.sha256(nativepath.read_bytes()).hexdigest()!=ORIGINAL_H7_SITE_SHA:
    raise RuntimeError("Original EEA monitoring source H7 native IDs changed")
native=json.loads(nativepath.read_text())["candidate_source_site_period_rows"]
stations={k:[] for k in sorted(fixed)}
for q in native:
    if q["panel"] not in EEA or q["HYBAS_ID_L07"] is None:continue
    key=(q["HYBAS_ID_L07"],"2017-2019" if q["panel"]=="EEA18" else "2020-2022")
    if key in stations:stations[key].append((q["panel"],q["native_id"]))
if any(not v for v in stations.values()):raise RuntimeError("Original EEA monitoring station H7 source missing")
eearows={}
for panel,(path,sha) in EEA.items():
    raw=pathlib.Path(path).read_bytes()
    if gitsha(raw)!=sha:raise RuntimeError("Original EEA published source CSV Git SHA changed")
    records=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
    if len(records)!={"EEA18":2634,"EEA22":2475}[panel]:raise RuntimeError("EEA original source year row count changed")
    eearows[panel]={r["site_key"]:r for r in records}
    if len(eearows[panel])!=len(records):raise RuntimeError("Original EEA site ID repeated in year")
site=[]
for (h,period),v in stations.items():
    for panel,sid in v:
        z=eearows[panel][sid]
        # Source-only lat/lon, geographic reach, monitoring method/country
        # (NO EQR or GBIF numerical endpoints are accessed).
        site.append({"H7":h,"PERIOD":period,"panel":panel,"site_id":sid,
           "reach":reach(z["HYRIV_ID"]),"lat":num(z["lat"]),"lon":num(z["lon"]),
           "country":z["countryCode"],"river_DIS_AV_CMS":num(z["DIS_AV_CMS"])})
if len(site)!=3525:raise RuntimeError("Source 3525 original monitoring year rows not preserved")

# Publisher original WWTP source
p=pathlib.Path(os.environ["ORIGINAL_WASTE_PUBLISHER_ZIP"])
b=p.read_bytes()
if hashlib.sha256(b).hexdigest()!=ORIGINAL_WASTE_SHA:raise RuntimeError("HydroWASTE original ZIP sha changed")
z=zipfile.ZipFile(io.BytesIO(b))
rw=z.read("HydroWASTE_v10.csv")
if hashlib.sha256(rw).hexdigest()!=ORIGINAL_WASTE_CSV_SHA:raise RuntimeError("HydroWASTE original raw CSV sha changed")
W=list(csv.DictReader(io.StringIO(rw.decode("latin-1"),newline="")))
if len(W)!=58502 or len({w["WASTE_ID"] for w in W})!=58502:
    raise RuntimeError("Genuine original 58502 WWTP registry changed")
byw=collections.defaultdict(list)
for w in W:
    id=reach(w["HYRIV_ID"])
    if id is None:continue
    x=num(w["LON_OUT"]);y=num(w["LAT_OUT"]);f=num(w["WASTE_DIS"])
    if x is None or y is None or f is None or f<0 or not (-90<=y<=90 and -180<=x<=180):continue
    byw[id].append({"plant_id":w["WASTE_ID"],"lat":y,"lon":x,"f":f,
       "quality":str(w["QUAL_WASTE"]),"status":str(w["STATUS"])})
site_reaches={x["reach"] for x in site if x["reach"] is not None}
candidate_reaches=sorted(site_reaches&set(byw))
site_with_same=sum(x["reach"] in byw for x in site)
if site_with_same!=471:
    raise RuntimeError("Original independent source study had 471 EEA site-years with same-reach registered WWTP; sample drift "+str(site_with_same))
print("SAME_REACH_SOURCE_CANDIDATES original_station_years",site_with_same,
      "native_reaches",len(candidate_reaches),"original_facs",sum(len(byw[x]) for x in candidate_reaches),flush=True)

# Original full published river shapefile. Iterate DBF records to map exact
# HYRIV_ID -> shapefile index and NEXT_DOWN; use SHX for random-access only
# targeted reach and next-down neighbor geometry. Do NOT publish raw shape.
err=None
for t in range(3):
    try:
        req=urllib.request.Request(ZIP_URL,headers={"Accept":"application/zip,application/octet-stream,*/*",
                   "User-Agent":"gesui-original-HYRIV-same-reach-geometry-source-20261113/1"})
        with urllib.request.urlopen(req,timeout=85) as rr: river_bytes=rr.read(95000000)
        if hashlib.sha256(river_bytes).hexdigest()!=SOURCE_RIVER_SHA:
            raise RuntimeError("Original official HydroRIVERS Europe original ZIP SHA changed")
        break
    except Exception as e:
        err=e;time.sleep(2+2*t)
else:raise RuntimeError("Original HydroRIVERS EU source inaccessible "+str(err))
zz=zipfile.ZipFile(io.BytesIO(river_bytes))
groups={}
for fn in zz.namelist():
    if "hydrorivers_v10_eu" in fn.lower() and fn.lower().endswith((".shp",".shx",".dbf")):
        name,ext=fn.rsplit(".",1)
        groups.setdefault(name,{})[ext.lower()]=fn
m=[x for x in groups.values() if set(x)=={"shp","shx","dbf"}]
if len(m)!=1:raise RuntimeError("Original HydroRIVERS SHP/SHX/DBF not unique")
part=m[0]
reader=shapefile.Reader(shp=io.BytesIO(zz.read(part["shp"])),
    shx=io.BytesIO(zz.read(part["shx"])),dbf=io.BytesIO(zz.read(part["dbf"])),encoding="latin-1")
headers=[q[0] for q in reader.fields[1:]]
if "HYRIV_ID" not in headers or "NEXT_DOWN" not in headers:
    raise RuntimeError("Original HydroRIVERS authoritative river attribute missing")
id_col=headers.index("HYRIV_ID"); next_col=headers.index("NEXT_DOWN")
lookup={};down={}
for i,rec in enumerate(reader.iterRecords()):
    k=reach(rec[id_col]);v=reach(rec[next_col])
    if k is None or k in lookup:raise RuntimeError("Original native river ID null or duplicate")
    lookup[k]=i
    if k in candidate_reaches:down[k]=v
if len(lookup)!=938544:raise RuntimeError("Original official HydroRIVERS 938544 reach network changed")
if not set(candidate_reaches).issubset(lookup):
    raise RuntimeError("Some original candidate same-reach EEA/WWTP IDs outside official EU network")

def as_line(id):
    if id is None or id not in lookup:return None
    shp=reader.shape(lookup[id])
    if not shp.points or len(shp.parts)!=1 or len(shp.points)<2:return None
    pts=shp.points
    if any(not (-180<=x<=180 and -90<=y<=90) for x,y in pts):return None
    if pts[0]==pts[-1]:return None
    return pts
def projected(pts,lon,lat):
    proj=CRS.from_proj4(f"+proj=aeqd +lat_0={lat:.8f} +lon_0={lon:.8f} +datum=WGS84 +units=m +no_defs")
    tf=Transformer.from_crs("EPSG:4326",proj,always_xy=True)
    return [tf.transform(x,y) for x,y in pts],tf
def get_geometry(k):
    pts=as_line(k)
    if pts is None:return {"status":"not_single_line"}
    nxt=down[k]
    np=as_line(nxt)
    if np is None:return {"status":"no_valid_NEXT_DOWN_neighbor"}
    lon0,lat0=pts[0]
    src_xy,tf=projected(pts,lon0,lat0)
    dst_xy=[tf.transform(x,y) for x,y in np]
    line=LineString(src_xy)
    if line.length<=250:return {"status":"reach_short_250m_or_less"}
    dline=LineString(dst_xy)
    # Endpoint to next downstream segment endpoint: no assumption about
    # whether source LINE coordinate 0 is its upstream or downstream end.
    a=Point(src_xy[0]);b=Point(src_xy[-1])
    nextends=[Point(dst_xy[0]),Point(dst_xy[-1])]
    da=min(a.distance(p) for p in nextends)
    db=min(b.distance(p) for p in nextends)
    if min(da,db)>500 or abs(da-db)<150:
        return {"status":"source_line_nextdown_endpoint_ambiguous",
                "closest_nextdown_endpoint_m":float(min(da,db)),
                "endpoints_advantage_m":float(abs(da-db))}
    downstream_which="first" if da<db else "last"
    return {"status":"oriented_from_NEXT_DOWN_neighbor",
       "polyline":line,"transform":tf,
       "river_length_projected_m":float(line.length),
       "NEXT_DOWN_reach":nxt,
       "downstream_endpoint":downstream_which,
       "downstream_connection_gap_m":float(min(da,db)),
       "other_endpoint_connection_gap_m":float(max(da,db))}
geometry={k:get_geometry(k) for k in candidate_reaches}
geometry_status=collections.Counter(q["status"] for q in geometry.values())
print("SAME_REACH_SOURCE_ORIENTATION_QUALITY",dict(geometry_status),flush=True)

# Source-pair ordering: each original station and each same native HYRIV
# outfall point must be <=500 m from original polyline; require 500m
# along-reach separation AND 2*both projected errors. Ambiguous stays NA.
source_rows=[];paircounts=collections.Counter();original_countries=collections.Counter()
for s in site:
    if s["reach"] not in byw:continue
    g=geometry[s["reach"]]
    orig=byw[s["reach"]]
    row={"HYBAS_ID_L07":s["H7"],"PERIOD":s["PERIOD"],
      "original_EEA_station":s["site_id"],"original_EEA_native_HYRIV_ID":s["reach"],
      "country":s["country"],"original_registered_same_reach_facilities_n":len(orig),
      "geometry_source_status":g["status"],"source_confident_upstream_plant_n":0,
      "source_confident_downstream_plant_n":0,"source_uncertain_same_reach_plant_n":0,
      "source_confident_upstream_WASTE_DIS_m3_day":0.,
      "source_confident_downstream_WASTE_DIS_m3_day":0.,
      "source_unresolved_same_reach_WASTE_DIS_m3_day":0.,
      "source_position_order_eligible":False}
    original_countries[s["country"]]+=1
    if g["status"]!="oriented_from_NEXT_DOWN_neighbor" or s["lat"] is None or s["lon"] is None:
        row["source_uncertain_same_reach_plant_n"]=len(orig)
        row["source_unresolved_same_reach_WASTE_DIS_m3_day"]=sum(q["f"] for q in orig)
        paircounts["source_not_geometry_orientation_eligible"]+=len(orig)
        source_rows.append(row);continue
    line=g["polyline"];tf=g["transform"]
    px=Point(tf.transform(s["lon"],s["lat"]))
    gap=line.distance(px)
    row["station_source_original_snap_distance_m"]=float(gap)
    if gap>500:
        row["source_uncertain_same_reach_plant_n"]=len(orig)
        row["source_unresolved_same_reach_WASTE_DIS_m3_day"]=sum(q["f"] for q in orig)
        paircounts["source_EEA_station_snap_distance_over_500m"]+=len(orig)
        source_rows.append(row);continue
    src_proj=float(line.project(px))
    if g["downstream_endpoint"]=="first":src_proj=line.length-src_proj
    row["source_station_linear_km_from_upstream_end"]=src_proj/1000
    row["source_position_order_eligible"]=True
    for item in orig:
        pt=Point(tf.transform(item["lon"],item["lat"]))
        ww_gap=line.distance(pt)
        if ww_gap>500:
            decision="source_outfall_snap_distance_over_500m"
        else:
            wwproj=float(line.project(pt))
            if g["downstream_endpoint"]=="first":wwproj=line.length-wwproj
            dist=src_proj-wwproj
            if abs(dist)<500 or abs(dist)<2*(gap+ww_gap):
                decision="source_station_spacing_too_close_to_infer_direction"
            elif dist>0:
                decision="true_same_reach_source_upstream_of_station"
            else:
                decision="true_same_reach_source_downstream_of_station"
        paircounts[decision]+=1
        if decision=="true_same_reach_source_upstream_of_station":
            row["source_confident_upstream_plant_n"]+=1
            row["source_confident_upstream_WASTE_DIS_m3_day"]+=item["f"]
        elif decision=="true_same_reach_source_downstream_of_station":
            row["source_confident_downstream_plant_n"]+=1
            row["source_confident_downstream_WASTE_DIS_m3_day"]+=item["f"]
        else:
            row["source_uncertain_same_reach_plant_n"]+=1
            row["source_unresolved_same_reach_WASTE_DIS_m3_day"]+=item["f"]
    if row["original_registered_same_reach_facilities_n"]!=(
       row["source_confident_upstream_plant_n"]+row["source_confident_downstream_plant_n"]+
       row["source_uncertain_same_reach_plant_n"]):
        raise RuntimeError("Source conservation: some same-reach facilities classified twice")
    source_rows.append(row)

if len(source_rows)!=471:raise RuntimeError("Exactly 471 original same-reach EEA station-years source entries lost")
pairs_total=sum(r["original_registered_same_reach_facilities_n"] for r in source_rows)
up=sum(r["source_confident_upstream_plant_n"] for r in source_rows)
dn=sum(r["source_confident_downstream_plant_n"] for r in source_rows)
unknown=sum(r["source_uncertain_same_reach_plant_n"] for r in source_rows)
if pairs_total!=up+dn+unknown:
    raise RuntimeError("Source conservation: total paired station×plant count changed")
sufficient={
 "n_stationyears_with_ORIENTED_reach":sum(r["geometry_source_status"]=="oriented_from_NEXT_DOWN_neighbor" for r in source_rows),
 "n_stationyears_with_original_position_order_eligible":sum(r["source_position_order_eligible"] for r in source_rows),
 "n_distinct_original_station_IDs_with_position_order_eligible":len({r["original_EEA_station"] for r in source_rows if r["source_position_order_eligible"]}),
 "n_original_countries_with_position_order_eligible":len({r["country"] for r in source_rows if r["source_position_order_eligible"]}),
 "n_source_station_WWTP_confident_upstream_pairs":up,
 "n_source_station_WWTP_confident_downstream_pairs":dn,
 "n_source_station_WWTP_original_order_UNRESOLVED_pairs":unknown,
 "n_original_2018_22_H7_source_basinperiods_with_at_least1_upstream_same_reach_source":len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in source_rows if r["source_confident_upstream_plant_n"]>0}),
 "n_original_2018_22_H7_source_basinperiods_with_at_least1_downstream_same_reach_source":len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in source_rows if r["source_confident_downstream_plant_n"]>0}),
}
gate=bool(
 sufficient["n_stationyears_with_original_position_order_eligible"]>=50
 and sufficient["n_distinct_original_station_IDs_with_position_order_eligible"]>=25
 and sufficient["n_original_countries_with_position_order_eligible"]>=5
 and up>=20)
report={
 "schema":"GESUI_EEA454_SOURCE_ONLY_TRUE_HYRIV_SAME_REACH_WWTP_POINT_ORDER_QA_20261113_v1",
 "original_HydroRIVERS_publisher_v10_eu_ZIP_SHA256":SOURCE_RIVER_SHA,
 "original_HydroRIVERS_exact_river_reach_count":len(lookup),
 "original_HydroWASTE_source_publisher_ZIP_SHA256":ORIGINAL_WASTE_SHA,
 "original_HydroWASTE_original_CSV_SHA256":ORIGINAL_WASTE_CSV_SHA,
 "original_EEA2018_and_2022_monitoring_GitBlobs":{q:v[1] for q,v in EEA.items()},
 "source_fixed_original_native_H7_monitor_geo_SHA256":ORIGINAL_H7_SITE_SHA,
 "original_study_454_H7_periods":len(fixed),
 "original_study_3525_stationyear_rows":len(site),
 "original_same_reach_471_stationyear_source":len(source_rows),
 "same_reach_native_river_reaches_with_registered_facility":len(candidate_reaches),
 "original_reach_geometry_orientation_rejection_counts":dict(geometry_status),
 "original_source_pair_classification_count":dict(paircounts),
 "original_source_pair_totals":{"all":pairs_total,"confident_UPSTREAM":up,"confident_DOWNSTREAM":dn,
                                  "unresolved":unknown},
 "original_country_counts_in_candidate_same_reach_stationyears":dict(original_countries),
 "predeclared_hydrological_geometry_quality_gates":{
   "endpoint_join_nearest_m_at_most":500,
   "near_vs_far_endpoint_advantage_m_at_least":150,
   "native_reach_min_length_m":250,
   "station_snap_distance_to_original_line_max_m":500,
   "WWTP_outfall_snap_distance_to_original_line_max_m":500,
   "source_station_separation_at_least_m":500,
   "source_station_separation_at_least_twice_total_snap_error":True
 },
 "source_original_geo_eligibility":sufficient,
 "source_extension_ecological_model_predeclared_gate_PASS":gate,
 "no_causal_dose_estimation_or_biotic_outcomes_read":True,
 "no_raw_HydroRIVERS_or_registered_original_WWTP_standalone_data_redistributed":True,
 "limitations":[
   "The original native river vertex orientation is NOT assumed; inferred from NEXT_DOWN neighbor line endpoint connectivity.",
   "A whole river reach could have wrong original geometry or station point snap; conservative endpoint and offset quality gates attempt to prevent incorrect ordering.",
   "Same-reach registered facilities cannot be considered discharged active plants without monitoring; most source status unknown.",
   "Only source-verified confidence subset can add point sources; other unresolved source WWTP are a source uncertainty upper bound, never true 0 recorded emission.",
   "This classification is posthoc relative to older outlet +11.1% and strictly-upstream +1.4% predictive EQR results; biological response cannot determine geometry eligibility."
 ]}
derived={
 "schema":report["schema"],
 "original_geo_source_454_frame":True,
 "same_reach_EEA_station_year_source_decisions":source_rows,
 "original_source_HydroRIVERS_ZIP_SHA256":SOURCE_RIVER_SHA,
 "no_EEA_EQR_plant_or_butterfly_outcomes_opened":True
}
raw=(json.dumps(derived,sort_keys=True,separators=(",",":"))+"\n").encode()
f=ROOT/"EEA454_HYRIV_SAME_REACH_POSITION_ORDER_SOURCE_ONLY.json"
f.write_bytes(raw)
report["source_471_original_stationyear_point_order_derived_manifest_SHA256"]=hashlib.sha256(raw).hexdigest()
(ROOT/"EEA454_ORIGINAL_SAME_REACH_WWTP_POSITION_ORDER_QA.json").write_text(json.dumps(report,indent=2,sort_keys=True))
print("HYRIV_SAME_REACH_SOURCE_ONLY_RESULT_START",flush=True)
print(json.dumps({k:v for k,v in report.items() if k!="limitations"},sort_keys=True),flush=True)
print("HYRIV_SAME_REACH_SOURCE_ONLY_RESULT_END",flush=True)
if not gate:
    raise RuntimeError("Same-reach position ordering hard source gate failed. Do not estimate new ecological model")
