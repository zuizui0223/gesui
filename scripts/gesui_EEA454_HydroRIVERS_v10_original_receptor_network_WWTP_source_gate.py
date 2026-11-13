#!/usr/bin/env python3
"""Original HydroRIVERS EU 68MB SHP ZIP DBF, native EEA 2018/22 observed
HYRIV_ID receiver reach, official HydroWASTE 58502 outfall HYRIV_ID.

SOURCE ONLY, NO EQR/Plant/Butterfly values opened, and no standalone
republication of licensed HydroRIVERS network. Computes true directed,
STRICT UPSTREAM (excluding unresolved same-reach sources) 50 km/full route
and full-upstream modeled flow. Explicit downstream 50 km direction-control.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,statistics,time,urllib.request,zipfile
import shapefile

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HydroRIVERS_EEA_RECEPTOR_SOURCE"
ROOT.mkdir(parents=True,exist_ok=True)
RIVERS_URL="https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip"
WASTE_ZIP_SHA="6e19d8f31d95d533b2b29056b823eeb915d933fa22ea090f234dc38f6427d044"
WASTE_CSV_SHA="af040d9b4fdc93bab2eb8d1b52985f5bacc6b0db9f16795c29e31c10fb615dd4"
H7_SITE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
EEA_BLOBS={
 "EEA18":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
def gitsha(x):return hashlib.sha1(("blob "+str(len(x))+"\0").encode()+x).hexdigest()
def reachkey(x):
    if x is None:return None
    s=str(x).strip()
    if not s or s.lower() in ("none","na","nan","null","0","-1"):return None
    try:
        k=int(float(s))
        if k<=0:return None
        return str(k)
    except ValueError:return None
def num(v):
    try:
        a=float(v)
        return a if math.isfinite(a) else None
    except (TypeError,ValueError):return None

# Validate sample selection using original 454-row RESPONSE-FREE identity. The
# previous CSV includes biological columns but ONLY keys are read in this
# source contract. This script never touches numerical EQR/GBIF responses.
with open(os.environ["ORIGINAL_FIXED_454_H7_COHORT"],encoding="utf-8-sig",newline="") as f:
    fixedkeys={(r["HYBAS_ID_L07"],r["PERIOD"]) for r in csv.DictReader(f)}
if len(fixedkeys)!=454 or len({k[0] for k in fixedkeys})!=232:
    raise RuntimeError("Original 454 source frame or 232 original H7 basins changed")
geo=pathlib.Path(os.environ["ORIGINAL_EEA_TO_H7_SOURCE_MANIFEST"])
if hashlib.sha256(geo.read_bytes()).hexdigest()!=H7_SITE_SHA:
    raise RuntimeError("Original 2018/2022 monitoring station geography SHA mismatch")
old=json.loads(geo.read_text())["candidate_source_site_period_rows"]
if len(old)!=6330:raise RuntimeError("Original native H7 station mapping n6330 changed")
site_refs=collections.defaultdict(list)
for r in old:
    if r["panel"] not in ("EEA18","EEA22") or r["HYBAS_ID_L07"] is None:continue
    period="2017-2019" if r["panel"]=="EEA18" else "2020-2022"
    key=(r["HYBAS_ID_L07"],period)
    if key in fixedkeys:site_refs[key].append(r["native_id"])
if set(site_refs)!=fixedkeys:raise RuntimeError("Missing original 454 native EEA station source keys")
source_stations={}
for kind,(path,sha) in EEA_BLOBS.items():
    raw=pathlib.Path(path).read_bytes()
    if gitsha(raw)!=sha:raise RuntimeError("Original EEA source Git blob differs "+kind)
    data=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
    if len(data)!={"EEA18":2634,"EEA22":2475}[kind]:
        raise RuntimeError("Unexpected true original EEA river station year rows")
    source_stations[kind]={r["site_key"]:{
        "reach":reachkey(r["HYRIV_ID"]),
        "lat":num(r["lat"]),"lon":num(r["lon"]),
        "observed_site_DIS_AV_CMS_original":num(r["DIS_AV_CMS"]),
        "country":r["countryCode"],
        "BQE_method_source":r["procedureClassificationSystem"]
        } for r in data}
    if len(source_stations[kind])!=len(data):
        raise RuntimeError("Original 2018/2022 river station ID duplication")
source_rows=[]
for (hid,period),ids in sorted(site_refs.items()):
    k="EEA18" if period=="2017-2019" else "EEA22"
    for site in ids:
        if site not in source_stations[k]:raise RuntimeError("Original H7 monitoring station identifier not in EEA source")
        source_rows.append({"H7":hid,"PERIOD":period,"EEA_original_native_site":site,**source_stations[k][site]})
if len(source_rows)!=3525:
    raise RuntimeError("Original fixed 454 H7 study had 3525 station-year sources; a source ID changed "+str(len(source_rows)))

# Authenticate the original public HydroWASTE publisher source, including 58k
# original unique registered facility records and source QUAL_WASTE/STATUS.
p=pathlib.Path(os.environ["ORIGINAL_HYDROWASTE_PUBLISHER_ZIP"])
wb=p.read_bytes()
if hashlib.sha256(wb).hexdigest()!=WASTE_ZIP_SHA:
    raise RuntimeError("Original HydroWASTE original publisher ZIP sha changed")
zw=zipfile.ZipFile(io.BytesIO(wb))
csvraw=zw.read("HydroWASTE_v10.csv")
if hashlib.sha256(csvraw).hexdigest()!=WASTE_CSV_SHA:
    raise RuntimeError("Original HydroWASTE publisher source CSV sha changed")
waste=list(csv.DictReader(io.StringIO(csvraw.decode("latin-1"),newline="")))
if len(waste)!=58502 or len({r["WASTE_ID"] for r in waste})!=58502:
    raise RuntimeError("Original HydroWASTE total or unique wastewater facility census changed")
if not all(set(["HYRIV_ID","WASTE_DIS","QUAL_WASTE","STATUS","WASTE_ID"]).issubset(r) for r in waste[:1]):
    raise RuntimeError("HydroWASTE original source fields not accessible")

def fetch(url):
    err=None
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/zip,application/octet-stream,*/*",
                "User-Agent":"gesui-HydroRIVERS-Original-Upstream-Station-SourceAudit-20261113/1.0"})
            with urllib.request.urlopen(req,timeout=85) as resp:
                b=resp.read(95000000)
                resolved=resp.url
            if len(b)>=95000000 or not zipfile.is_zipfile(io.BytesIO(b)):
                raise RuntimeError("Not a genuine original HydroRIVERS Europe ZIP")
            return b,resolved
        except Exception as e:
            err=e
            time.sleep(2+3*attempt)
    raise RuntimeError("Official original HydroRIVERS source fetch unavailable "+str(err))
rbytes,resolved=fetch(RIVERS_URL)
rsha=hashlib.sha256(rbytes).hexdigest()
z=zipfile.ZipFile(io.BytesIO(rbytes))
dbfs=[m for m in z.infolist() if m.filename.lower().endswith(".dbf") and "hydrorivers_v10_eu" in m.filename.lower()]
if len(dbfs)!=1:raise RuntimeError("Original HydroRIVERS native Europe DBF not unique "+str([q.filename for q in dbfs]))
dbf=dbfs[0]
if dbf.file_size>=550_000_000:
    raise RuntimeError("Unexpected large EU river reach attribute source; do not take unbounded memory")
rraw=z.read(dbf)
rdbf_sha=hashlib.sha256(rraw).hexdigest()
reader=shapefile.Reader(dbf=io.BytesIO(rraw),encoding="latin-1")
heads=[field[0] for field in reader.fields[1:]]
REQUIRED={"HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS"}
if not REQUIRED.issubset(set(heads)):
    raise RuntimeError("Genuine original HydroRIVERS v1 EU source schema changed "+str(heads))
# Preserve native graph 1:1 IDs without simplifying or georeferencing to other
# sources. The original reach length will be used only for the 50km path test.
nextdown={};length={};reachflow={};badsource=collections.Counter()
for q in reader.iterRecords():
    data=dict(zip(heads,q))
    id=reachkey(data["HYRIV_ID"])
    if id is None or id in nextdown:
        raise RuntimeError("Original HydroRIVERS EU source has null or repeated reach key")
    dest=reachkey(data["NEXT_DOWN"])
    le=num(data["LENGTH_KM"])
    dis=num(data["DIS_AV_CMS"])
    if le is None or le<0:
        badsource["bad_source_reach_length"]+=1
        le=0
    if dis is None or dis<0:badsource["bad_source_estimated_reach_flow"]+=1
    nextdown[id]=dest
    length[id]=le
    reachflow[id]=dis
n_reach=len(nextdown)
if n_reach<100000:raise RuntimeError("Source Europe HydroRIVERS reaches implausibly few")
missingout=collections.Counter()
for id,dest in nextdown.items():
    if dest is not None and dest not in nextdown:
        missingout["downstream_reach_outside_source_region"]+=1

# Site and WWTP original HYRIV coverage, limited strictly to EU network.
native_sites={x["reach"] for x in source_rows if x["reach"] is not None}
native_match={x for x in native_sites if x in nextdown}
site_match=sum(r["reach"] in nextdown for r in source_rows)
for r in source_rows:
    r["source_original_HYRIV_in_europe_network"]=bool(r["reach"] in nextdown)
if len(native_match)<500:raise RuntimeError("Original EEA station reach match unexpectedly sparse")

statuses=collections.Counter()
quality=collections.Counter()
local_count=collections.Counter()
local_flow=collections.defaultdict(float)
operational_count=collections.Counter()
operational_flow=collections.defaultdict(float)
qual1_count=collections.Counter()
qual1_flow=collections.defaultdict(float)
out_source=collections.Counter()
matched_original_facility_reaches=set()
for w in waste:
    status=str(w.get("STATUS","")).strip()
    statuses[status]+=1
    qual=str(w.get("QUAL_WASTE","")).strip()
    quality[qual]+=1
    reach=reachkey(w.get("HYRIV_ID"))
    f=num(w.get("WASTE_DIS"))
    if reach is None or reach not in nextdown:
        out_source["global_WWTP_not_on_EU_network_or_missing_original_HYRIV_ID"]+=1
        continue
    if f is None or f<0:
        out_source["matched_EU_reach_invalid_discharge_flow"]+=1
        continue
    matched_original_facility_reaches.add(reach)
    local_count[reach]+=1
    local_flow[reach]+=f
    if status.lower()=="operational":
        operational_count[reach]+=1
        operational_flow[reach]+=f
    if qual=="1":
        qual1_count[reach]+=1
        qual1_flow[reach]+=f
    out_source["WWTP_original_HYRIV_ID_matches_europe_network"]+=1
if out_source["WWTP_original_HYRIV_ID_matches_europe_network"]<4000:
    raise RuntimeError("HydroWASTE and HydroRIVERS original HYRIV_ID may use incompatible versions")

# Exact downstream DAG (each reach at most one NEXT_DOWN); topologically
# accumulate all registered model sources to downstream WITHOUT same reach
# source ambiguity. Source network cycles are a hard stop, not silently cut.
reverse=collections.defaultdict(list)
indegree=collections.Counter({x:0 for x in nextdown})
for id,dest in nextdown.items():
    if dest in nextdown:
        if dest==id:raise RuntimeError("Self-loop in original HydroRIVERS source")
        reverse[dest].append(id);indegree[dest]+=1
queue=collections.deque(k for k,d in indegree.items() if d==0)
c_count=collections.Counter(local_count);c_flow=collections.defaultdict(float,local_flow)
c_oper=collections.Counter(operational_count);c_operflow=collections.defaultdict(float,operational_flow)
c_q1=collections.Counter(qual1_count);c_q1flow=collections.defaultdict(float,qual1_flow)
processed=0
while queue:
    x=queue.popleft()
    processed+=1
    downstream=nextdown[x]
    if downstream in nextdown:
        c_count[downstream]+=c_count[x];c_flow[downstream]+=c_flow[x]
        c_oper[downstream]+=c_oper[x];c_operflow[downstream]+=c_operflow[x]
        c_q1[downstream]+=c_q1[x];c_q1flow[downstream]+=c_q1flow[x]
        indegree[downstream]-=1
        if indegree[downstream]==0:queue.append(downstream)
if processed!=n_reach:
    raise RuntimeError("HydroRIVERS Europe original NEXT_DOWN graph cycle; do not claim upstream exposure")
# For 50km source: a FULL reach length is added for every distinct upstream
# reach containing a WWTP; receiver reach is not counted. This means the
# 50km window is a **conservative source opportunity lower bound**.
def reachable50(receiver):
    remaining=collections.deque([(receiver,0.0)])
    seen={receiver}
    up_n=0;up_flow=0.;up_oper=0;up_operflow=0.;up_q1=0;up_q1flow=0.
    visited=0
    while remaining:
        child,km=remaining.popleft()
        for up in reverse.get(child,()):
            if up in seen:raise RuntimeError("Directed source graph has duplicated upstream reach path")
            seen.add(up)
            newdist=km+length[up]
            if newdist>50:continue
            remaining.append((up,newdist))
            visited+=1
            up_n+=local_count[up];up_flow+=local_flow[up]
            up_oper+=operational_count[up];up_operflow+=operational_flow[up]
            up_q1+=qual1_count[up];up_q1flow+=qual1_flow[up]
    # Directional-only downstream placebo: distinct reaches up to sum of
    # full reach lengths <= 50km. Not a matched randomized placebo.
    down_n=0;down_flow=0.;down_nodes=0;cur=receiver;d=0.;seen_down={cur}
    while True:
        down=nextdown[cur]
        if down is None or down not in nextdown:break
        if down in seen_down:raise RuntimeError("Original source graph downstream cycle")
        seen_down.add(down)
        d+=length[down]
        if d>50:break
        down_n+=local_count[down];down_flow+=local_flow[down]
        down_nodes+=1;cur=down
    return {
      "upstream_50km_distinct_reach_registered_WWTP_count":up_n,
      "upstream_50km_distinct_reach_registered_WWTP_WASTE_DIS_m3_day":up_flow,
      "upstream_50km_original_STATUS_operational_count":up_oper,
      "upstream_50km_original_STATUS_operational_model_flow_m3_day":up_operflow,
      "upstream_50km_original_QUAL_WASTE_1_count":up_q1,
      "upstream_50km_original_QUAL_WASTE_1_model_flow_m3_day":up_q1flow,
      "upstream_50km_complete_reaches_traversed":visited,
      "downstream_50km_distinct_reach_registered_WWTP_count_direction_control":down_n,
      "downstream_50km_distinct_reach_registered_WWTP_flow_m3_day_direction_control":down_flow,
      "downstream_50km_complete_reaches_traversed":down_nodes
    }
unique_site_reaches=sorted(native_match)
per_reach={}
for i,k in enumerate(unique_site_reaches):
    flux50=reachable50(k)
    per_reach[k]={
      "station_reach_original_HYRIV_ID":k,
      "source_original_same_reach_facilities_unresolved_COUNT":local_count[k],
      "source_original_same_reach_facilities_unresolved_WASTE_DIS_m3_day":local_flow[k],
      "strict_upstream_all_distinct_reach_source_registered_WWTP_count":c_count[k]-local_count[k],
      "strict_upstream_all_distinct_reach_WASTE_DIS_m3_day":max(0,c_flow[k]-local_flow[k]),
      "strict_upstream_all_distinct_reach_source_QUAL1_flow_m3_day":max(0,c_q1flow[k]-qual1_flow[k]),
      "strict_upstream_all_distinct_reach_SOURCE_STATUS_operational_flow_m3_day":max(0,c_operflow[k]-operational_flow[k]),
      "original_network_reach_DIS_AV_CMS":reachflow[k],
      "original_network_reach_LENGTH_KM":length[k],
      **flux50}
    if i and i%500==0:print("Actual EEA receptor reach source verified",i,"of",len(unique_site_reaches),flush=True)
positive=collections.Counter()
source_rows2=[]
by_H7=collections.defaultdict(list)
for r in source_rows:
    out=dict(r)
    id=r["reach"]
    if id not in per_reach:
        out["routed_source_valid"]=False
        source_rows2.append(out);by_H7[(r["H7"],r["PERIOD"])].append(out)
        continue
    x=per_reach[id]
    out.update(x)
    f=r["observed_site_DIS_AV_CMS_original"]
    if f is None or f<=0:
        out["routed_source_valid"]=False
    else:
        out["routed_source_valid"]=True
        out["strict_upstream_WWTP_fraction_of_original_model_river_discharge"]=x["strict_upstream_all_distinct_reach_WASTE_DIS_m3_day"]/(f*86400)
        out["strict_upstream_50km_WWTP_fraction_of_original_model_river_discharge"]=x["upstream_50km_distinct_reach_registered_WWTP_WASTE_DIS_m3_day"]/(f*86400)
        out["downstream_50km_WWTP_fraction_direction_control"]=x["downstream_50km_distinct_reach_registered_WWTP_flow_m3_day_direction_control"]/(f*86400)
        if x["strict_upstream_all_distinct_reach_source_registered_WWTP_count"]>0:
            positive["source_station_reach_strict_upstream_positive"]+=1
        if x["upstream_50km_distinct_reach_registered_WWTP_count"]>0:
            positive["source_station_reach_strict_upstream_50km_positive"]+=1
        if x["downstream_50km_distinct_reach_registered_WWTP_count_direction_control"]>0:
            positive["source_station_reach_DOWNSTREAM_50km_positive"]+=1
        if x["source_original_same_reach_facilities_unresolved_COUNT"]>0:
            positive["source_station_reach_same_reach_UNRESOLVED_positive"]+=1
    source_rows2.append(out);by_H7[(r["H7"],r["PERIOD"])].append(out)
if len(source_rows2)!=3525 or len(by_H7)!=454:raise RuntimeError("Source-only original station numbers changed after graph routing")
summary_rows=[]
for (h,period),group in sorted(by_H7.items()):
    valid=[x for x in group if x["routed_source_valid"]]
    d={"HYBAS_ID_L07":h,"PERIOD":period,
       "original_EEA_stations":len(group),
       "original_usable_network_receiver_stations":len(valid),
       "network_source_valid":bool(valid)}
    if valid:
        keys=[
            "strict_upstream_all_distinct_reach_source_registered_WWTP_count",
            "strict_upstream_all_distinct_reach_WASTE_DIS_m3_day",
            "strict_upstream_all_distinct_reach_SOURCE_STATUS_operational_flow_m3_day",
            "strict_upstream_all_distinct_reach_source_QUAL1_flow_m3_day",
            "strict_upstream_WWTP_fraction_of_original_model_river_discharge",
            "source_original_same_reach_facilities_unresolved_COUNT",
            "upstream_50km_distinct_reach_registered_WWTP_count",
            "upstream_50km_distinct_reach_registered_WWTP_WASTE_DIS_m3_day",
            "upstream_50km_original_STATUS_operational_model_flow_m3_day",
            "upstream_50km_original_QUAL_WASTE_1_model_flow_m3_day",
            "strict_upstream_50km_WWTP_fraction_of_original_model_river_discharge",
            "downstream_50km_distinct_reach_registered_WWTP_count_direction_control",
            "downstream_50km_distinct_reach_registered_WWTP_flow_m3_day_direction_control",
            "downstream_50km_WWTP_fraction_direction_control"
        ]
        for key in keys:
            d["station_equal_mean_"+key]=statistics.mean(x[key] for x in valid)
        d["n_station_with_strict_upstream_source"]=sum(x["strict_upstream_all_distinct_reach_source_registered_WWTP_count"]>0 for x in valid)
        d["n_station_with_50km_upstream_source"]=sum(x["upstream_50km_distinct_reach_registered_WWTP_count"]>0 for x in valid)
        d["n_station_with_downstream_50km_source"]=sum(x["downstream_50km_distinct_reach_registered_WWTP_count_direction_control"]>0 for x in valid)
    summary_rows.append(d)
num_valid=sum(x["routed_source_valid"] for x in source_rows2)
valid_basinperiod=[r for r in summary_rows if r["network_source_valid"]]
positive_basinperiod=[r for r in summary_rows if r["network_source_valid"] and r["n_station_with_strict_upstream_source"]>0]
positive50_basinperiod=[r for r in summary_rows if r["network_source_valid"] and r["n_station_with_50km_upstream_source"]>0]
positivebasin={x["HYBAS_ID_L07"] for x in positive_basinperiod}
gate=bool(
    num_valid/3525>=0.95 and len(valid_basinperiod)/454>=.90
    and len(positive_basinperiod)>=150 and len(positivebasin)>=75)
report={
 "schema":"GESUI_EU454_ORIGINAL_HYDRO_RIVERS_STATION_RECEPTOR_WWTP_SOURCE_GATE_v1",
 "network_public_source_original_url":RIVERS_URL,
 "network_source_resolved_url":resolved,
 "source_HydroRIVERS_v10_EU_ZIP_SHA256":rsha,
 "source_HydroRIVERS_v10_EU_attribute_DBf_SHA256":rdbf_sha,
 "source_HydroRIVERS_v10_EU_archive_file_bytes":len(rbytes),
 "source_HydroRIVERS_v10_EU_reaches":n_reach,
 "source_HydroRIVERS_v10_EU_original_DBf_fields":heads,
 "source_HydroRIVERS_original_reach_dbf_size":len(rraw),
 "original_graph_all_reaches_acyclic":processed==n_reach,
 "original_graph_source_missing_downstream_links":dict(missingout),
 "original_network_original_bad_attribute_count":dict(badsource),
 "original_HydroWASTE_zipped_global_58k_SHA256":WASTE_ZIP_SHA,
 "original_HydroWASTE_original_CSV_SHA256":WASTE_CSV_SHA,
 "original_HydroWASTE_status_all_facilities":dict(statuses),
 "original_HydroWASTE_discharge_quality_all_facilities":dict(quality),
 "original_HydroWASTE_EU_network_reach_matched_sources":dict(out_source),
 "original_HydroWASTE_EU_unique_reaches_with_facility":len(matched_original_facility_reaches),
 "source_original_EEA_benthic_2018_2022_github_blobs":{k:v[1] for k,v in EEA_BLOBS.items()},
 "source_original_EEA_to_HYBAS7_native_geo_manifest_SHA256":H7_SITE_SHA,
 "same_original_fixed_H7_study_rows":454,
 "same_original_fixed_distinct_H7":232,
 "source_original_3525_EEA_site_years":len(source_rows),
 "source_original_distinct_EEA_HYRIV_ID":len(native_sites),
 "source_original_distinct_reaches_matched_EU_HydroRIVERS":len(native_match),
 "source_original_EEA_site_years_network_reach_valid":site_match,
 "source_original_EEA_site_years_network_and_valid_river_flow":num_valid,
 "source_H7_periods_at_least1_receiver_station_routable":len(valid_basinperiod),
 "source_H7_periods_at_least1_strict_upstream_registered_WWTP":len(positive_basinperiod),
 "source_H7_periods_at_least1_50km_full_reach_UPSTREAM_WWTP":len(positive50_basinperiod),
 "source_distinct_H7_basins_with_strict_upstream_registered_WWTP":len(positivebasin),
 "source_H7_periods_with_50km_DOWNSTREAM_WWTP_direction_control":sum(x["network_source_valid"] and x["n_station_with_downstream_50km_source"]>0 for x in summary_rows),
 "source_EEA_site_3d_distance_direction_controls":dict(positive),
 "source_original_network_reach_same_reach_sources_excluded":True,
 "source_50km_full_reach_distance_upper_bound_conservative":True,
 "source_flow_model_m3_day_divided_by_original_EEA_station_cms_times86400":True,
 "source_reported_QUAL_WASTE_1_and_STATUS_operational_separate":True,
 "original_EEA_454_receptor_SOURCE_GATE_pass":gate,
 "predeclared_gate_limits":{"station_year_match_fraction_min":0.95,"H7_year_any_source_valid_fraction_min":0.90,
    "H7_periods_at_least1_strict_upstream_facility_min":150,"distinct_H7_with_upstream_facility_min":75},
 "source_no_benthic_EQR_plant_S50_or_lycaenidae_values_examined":True,
 "original_HydroRIVERS_raw_ZIP_not_distributed_due_to_publisher_license":True,
 "source_original_wastewater_registration_does_not_imply_active":True,
 "no_causal_pollution_or_biodiversity_claims":True
}
# Only derived original WWTP/monitoring station counts and source reader IDs,
# never distribute raw HydroRIVERS licensed standalone source.
p=ROOT/"EU454_ORIGINAL_EEA_MONITOR_RECEIVER_HYRIV_WWTP_NETWORK_SOURCE_ONLY.json"
d={"schema":report["schema"],"H7_year_original_derived_receiver_site_mean_exposures":summary_rows,
   "source_original_HydroRIVERS_EU_ZIP_SHA256":rsha,
   "source_original_HydroWASTE_ZIP_SHA256":WASTE_ZIP_SHA,
   "source_only_no_EQR_biology_read":True}
raw=(json.dumps(d,sort_keys=True,separators=(",",":"))+"\n").encode()
p.write_bytes(raw)
report["source_derived_receiver_WWTP_454_basinperiod_manifest_SHA256"]=hashlib.sha256(raw).hexdigest()
(ROOT/"EU454_HYDRO_RIVERS_TRUE_STATION_UPSTREAM_WWTP_REPORT.json").write_text(json.dumps(report,indent=2,sort_keys=True))
print("EU_HYDRO_RIVERS_RECEPTOR_SOURCE_GATE_START",flush=True)
print(json.dumps(report,sort_keys=True),flush=True)
print("EU_HYDRO_RIVERS_RECEPTOR_SOURCE_GATE_END",flush=True)
if not gate:raise RuntimeError("Official HydroRIVERS to original EEA monitoring reach source support gate failed, no ecological model allowed")
