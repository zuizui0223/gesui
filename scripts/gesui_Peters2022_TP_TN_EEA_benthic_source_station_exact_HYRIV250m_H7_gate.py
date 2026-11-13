#!/usr/bin/env python3
"""SOURCE ONLY original measured TP2022 (Peters WQ registry and original real
2022 measured river TP) x actual EEA2022 observed benthic EQR station position
same native HYRIV reach, within 250m primary, no outcome values accessed.

Population locks original 454 H7×period macro cohort and original 6330
station/HYBAS7 source, BEFORE measured TP vs EQR or other 3biota coefficient.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,statistics

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_2022_EEA_Peters_same_site_source"
ROOT.mkdir(parents=True,exist_ok=True)
SOURCES={
 "TP":("results/eea2022_measured_nutrient_load_river_20261107/TP_2022_reported_wwtp_load_and_measured_river_station.csv","a75a4cf140de64fd6038079a3fac6756a461f3f8"),
 "TN":("results/eea2022_measured_nutrient_load_river_20261107/TN_2022_reported_wwtp_load_and_measured_river_station.csv","05273bc667040dc484a074c823e60ab51c378fdf"),
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
H7SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
def blobsha(raw):return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
def read_source(k):
    path,sha=SOURCES[k]
    raw=pathlib.Path(path).read_bytes()
    if blobsha(raw)!=sha:raise RuntimeError("Original source Git blob hash changed "+k)
    rr=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
    if len(rr)!={"TP":1291,"TN":159,"EEA22":2475}[k]:
        raise RuntimeError("Original source station count not source-verified "+k)
    return rr

cpath=pathlib.Path(os.environ["ORIGINAL_454_H7_SAME_THREE_LAYER_COHORT"])
with cpath.open(encoding="utf-8-sig",newline="") as f:
    old=csv.DictReader(f)
    originalH7={r["HYBAS_ID_L07"] for r in old if r["PERIOD"]=="2020-2022"}
if len(originalH7)!=227:
    raise RuntimeError("Original native macro fixed 2022 study not exactly 227 H7 units")
original_geo_path=pathlib.Path(os.environ["ORIGINAL_6330_EEA_EPA_H7_STATION_IDS"])
if hashlib.sha256(original_geo_path.read_bytes()).hexdigest()!=H7SHA:
    raise RuntimeError("Original geographic source 6330 station SHA256 changed")
station_geo=json.loads(original_geo_path.read_text())["candidate_source_site_period_rows"]
if len(station_geo)!=6330:raise RuntimeError("Original 6330 source year/ID and H7 region missing")
h7_2022={}
for row in station_geo:
    if row["panel"]=="EEA22" and row["HYBAS_ID_L07"] is not None:
        key=row["native_id"]
        if key in h7_2022:raise RuntimeError("Duplicate original EEA22 native station ID H7")
        h7_2022[key]=row["HYBAS_ID_L07"]
if len(h7_2022)!=2224:
    raise RuntimeError("Expected 2224 original EEA2022 H7 source-positive geometry matches")

eqar=read_source("EEA22")
# This stage reads ONLY source metadata, NEVER numeric EQR value, its sign or
# reported WWTP load. Source+year EQR station ID and original HYRIV match.
qa=[]
by_reach=collections.defaultdict(list)
for r in eqar:
    d={
      "site_key":r["site_key"],"HYRIV_ID":r["HYRIV_ID"],
      "country":r["countryCode"],"lat":float(r["lat"]),"lon":float(r["lon"]),
      "BQE_method":r["procedureClassificationSystem"],
      "H7_id":h7_2022.get(r["site_key"])
    }
    qa.append(d)
    by_reach[d["HYRIV_ID"]].append(d)
assert len(qa)==2475
def km(a,b,c,d):
    R=6371.0088
    p1,p2=map(math.radians,(a,c))
    deltaP=math.radians(c-a)
    deltaL=math.radians(d-b)
    xx=math.sin(deltaP/2)**2+math.cos(p1)*math.cos(p2)*math.sin(deltaL/2)**2
    return 2*R*math.asin(min(1,math.sqrt(xx)))
if any(not(-90<=r["lat"]<=90 and -180<=r["lon"]<=180) for r in qa):
    raise RuntimeError("Invalid source original 2022 EEA coordinate")
results={}
for group in ("TP","TN"):
    water=read_source(group)
    matches_all=[]
    eligible=[]
    all_near=[]
    for r in water:
        WQ=r["WQ.ID"]
        lat=float(r["WQ.lat"]);lon=float(r["WQ.lon"])
        reach=r["HYRIV_ID"]
        original_visit_count=int(r["actual_2022_observation_count"])
        nmonths=int(r["distinct_2022_months"])
        if original_visit_count<4 or nmonths<3:
            raise RuntimeError("Original measured nutrient 2022 metadata 4+ visits/3+months gate source changed")
        nearest=sorted([(1000*km(lat,lon,x["lat"],x["lon"]),x) for x in by_reach.get(reach,[])],
                       key=lambda q:(q[0],q[1]["site_key"]))
        if nearest:
            d,x=nearest[0]
            item={"nutrient":group,"WQ_ID":WQ,"EEA2022_site_key":x["site_key"],
                  "EEA2022_BQE_original_method":x["BQE_method"],
                  "EEA2022_country":x["country"],"HYRIV_ID":reach,
                  "EEA_H7_L07":x["H7_id"],
                  "spherical_source_coordinate_separation_m":d,
                  "real_2022_WQ_valid_sample_count":original_visit_count,
                  "real_2022_WQ_distinct_months":nmonths,
                  "same_2022_calendar_year":True,
                  "2018_2022_EEA_EQR_not_read":True,
                  "2022_river_TP_TN_median_not_read":True}
            all_near.append(item)
            if d<=250:
                matches_all.append(item)
                if x["H7_id"] in originalH7:eligible.append(item)
    # Enforce 1:1 site-matches based on source geometry not biological values.
    def stable_1to1(v):
        seenWQ=set();seenEEA=set();out=[]
        for x in sorted(v,key=lambda z:(z["spherical_source_coordinate_separation_m"],z["WQ_ID"],z["EEA2022_site_key"])):
            if x["WQ_ID"] in seenWQ or x["EEA2022_site_key"] in seenEEA:continue
            out.append(x);seenWQ.add(x["WQ_ID"]);seenEEA.add(x["EEA2022_site_key"])
        return sorted(out,key=lambda z:(z["EEA2022_site_key"],z["WQ_ID"]))
    matched=stable_1to1(matches_all)
    macros=stable_1to1(eligible)
    def stats(src):
        return {
         "n_WQ_EEA_unique_pairs":len(src),
         "n_unique_H7_with_geometry":len({z["EEA_H7_L07"] for z in src if z["EEA_H7_L07"] is not None}),
         "n_unique_EEA_physical_stations":len({z["EEA2022_site_key"] for z in src}),
         "source_countries":dict(collections.Counter(z["EEA2022_country"] for z in src)),
         "source_country_count":len({z["EEA2022_country"] for z in src}),
         "source_EEA_BQE_methods":dict(collections.Counter(z["EEA2022_BQE_original_method"] for z in src)),
         "min_geodesic_source_station_offset_m":min((z["spherical_source_coordinate_separation_m"] for z in src),default=None),
         "median_geodesic_source_station_offset_m":statistics.median(z["spherical_source_coordinate_separation_m"] for z in src) if src else None,
         "n_exact_zero_coordinate_source_pairs":sum(z["spherical_source_coordinate_separation_m"]<0.01 for z in src)
        }
    results[group]={
     "original_all_EEA2022_stations":2475,
     "original_all_nutrient_2022_WQ_stations":len(water),
     "WQ_with_exact_same_native_HYRIV_ID_in_2022_EEA":len(all_near),
     "WQ_with_same_native_HYRIV_and_nearest_EEA_250m_before_1to1":len(matches_all),
     "source_all_EEA2022_native_reach_250m_1to1":stats(matched),
     "source_FIXED_original_227_H7_2022_cohort_native_reach_250m_1to1":stats(macros),
     "source_reference_same_native_reach_500m":sum(z["spherical_source_coordinate_separation_m"]<=500 for z in all_near),
     "source_reference_same_native_reach_2000m":sum(z["spherical_source_coordinate_separation_m"]<=2000 for z in all_near),
     "fixed_H7_2022_cohort_WQ_to_EEA_station_year_source_pairs":macros,
     "matched_WQ_EEA_all_EEA2022_not_limited_to_H7_cohort":matched,
     "overlap_is_2022_original_measured_nutrient_station_to_2022_original_river_benthic_station_not_GBIF_plant_sampling_site":True
    }
mac=results["TP"]["source_FIXED_original_227_H7_2022_cohort_native_reach_250m_1to1"]
countries=mac["source_countries"]
fr_n=countries.get("FR",0)
fr_H7=len({z["EEA_H7_L07"] for z in results["TP"]["fixed_H7_2022_cohort_WQ_to_EEA_station_year_source_pairs"] if z["EEA2022_country"]=="FR"})
primary_pass=mac["n_WQ_EEA_unique_pairs"]>=100 and mac["n_unique_H7_with_geometry"]>=40 and mac["source_country_count"]>=3
fr_only_pass=fr_n>=80 and fr_H7>=30
out={
 "schema":"GESUI_REAL_MEASURED_2022_PETERS_TP_TN_TO_SAME_EEA_EQR_HYRIV_250M_AND_FIXED227H7_SOURCE_ONLY_v1",
 "source_original_H7_EEA_station_manifest_SHA256":H7SHA,
 "source_EEA2022_WQ_GitBlob_SHA1":{k:v[1] for k,v in SOURCES.items()},
 "original_original_454_H7_period_cohort_n":454,
 "original_fixed_2022_H7_population":227,
 "original_original_2022_EEA_BQE_stations":2475,
 "source_fixed_same_native_reach_and_250m_threshold_before_new_ecological_association":True,
 "source_fixed_nonoverlapping_nearest_pair_matching_by_distance_and_string_ties":True,
 "source_original_WQ_observation_count_and_date_requirements_at_least4_samples_3months":True,
 "all_source_2022_site_overlap":{
    k:{x:v for x,v in z.items() if x not in ("fixed_H7_2022_cohort_WQ_to_EEA_station_year_source_pairs","matched_WQ_EEA_all_EEA2022_not_limited_to_H7_cohort")}
    for k,z in results.items()},
 "fixed_2022_macro_H7_TP_FR_station_count":fr_n,
 "fixed_2022_macro_H7_TP_FR_unique_H7_basins":fr_H7,
 "original_source_pan_EU_primary_gate_pass":primary_pass,
 "original_source_FR_only_true_2022_same_station_gate_pass":fr_only_pass,
 "original_source_ANY_predeclared_new_analysis_gate_pass":bool(primary_pass or fr_only_pass),
 "no_original_TP_TN_concentrations_or_numeric_EQR_read_into_selection":True,
 "do_not_test_TN_macro_if_30sites_and_15H7_not_met":True,
 "no_causal_wwtp_effect_or_realized_pollination":True
}
manifest={
 "schema":out["schema"],
 "original_454_H7_cohort_2022_genuine_TP_250m_source_matched_pairs":results["TP"]["fixed_H7_2022_cohort_WQ_to_EEA_station_year_source_pairs"],
 "original_454_H7_cohort_2022_genuine_TN_250m_source_matched_pairs":results["TN"]["fixed_H7_2022_cohort_WQ_to_EEA_station_year_source_pairs"],
 "original_native_EEA_EQR_response_and_Peters_measured_river_nutrients_not_used_to_match":True
}
raw=(json.dumps(manifest,sort_keys=True,separators=(",",":"))+"\n").encode()
out["source_2022_genuine_EEA_Peters_same_reach_station_250m_cohort_SHA256"]=hashlib.sha256(raw).hexdigest()
(ROOT/"EEA2022_PETERS_TP_TN_HYRIV_250M_TO_FIXED_227_H7_SOURCE_MATCH.json").write_bytes(raw)
(ROOT/"EEA2022_PETERS_TP_TN_EEA_STATION_250M_SOURCE_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("EEA2022_PETERS_TP_TN_STATION_OVERLAP_SOURCE_START",flush=True)
print(json.dumps(out,sort_keys=True),flush=True)
print("EEA2022_PETERS_TP_TN_STATION_OVERLAP_SOURCE_END",flush=True)
if not (primary_pass or fr_only_pass):
    raise RuntimeError("Fixed same-site water-chemistry/benthic source overlap gate FAILED. No ecological response models")
