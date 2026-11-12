#!/usr/bin/env python3
"""Outcome-blind EEA 2018→2022 exact same station/assessment method source pair
quality gate on original, fixed 234 HydroBASINS L07 EU repeat source cohort.
Do NOT read EQR numerical outcome or GBIF new facet/5-family values.
"""
import collections,csv,hashlib,io,json,os,pathlib

SITES_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
EU234_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189"
BLOBS={"2018":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
       "2022":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")}
path=pathlib.Path(os.environ["HYBAS7_ORIGINAL_NATIVE_RIVER_SOURCE"])
f=pathlib.Path(os.environ["H7_234_PRE_NEW_RESPONSE_SOURCE"])
if hashlib.sha256(path.read_bytes()).hexdigest()!=SITES_SHA or hashlib.sha256(f.read_bytes()).hexdigest()!=EU234_SHA:
    raise RuntimeError("Original same-basins EQR source pre-outcome identity changed")
orig=json.loads(path.read_text())["candidate_source_site_period_rows"]
eu234=json.loads(f.read_text())["original_two_year_EEA_H7_subbasins"]
basins={x["HYBAS_ID_L07"] for x in eu234}
assert len(basins)==234
joined={}
for z in orig:
    if z["panel"] not in ("EEA18","EEA22"):continue
    hid=z["HYBAS_ID_L07"]
    if hid is not None and hid in basins:
        key=(z["panel"],z["native_id"])
        if key in joined:raise RuntimeError("Source original river monitoring ID was duplicated")
        joined[key]=hid
def blobsha(raw):return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
data={}
for period,(name,expected) in BLOBS.items():
    raw=pathlib.Path(name).read_bytes()
    if blobsha(raw)!=expected:raise RuntimeError("Original EEA year source Git blob mismatch")
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    d={}
    for x in reader:
        key=x["site_key"]
        if key in d:raise RuntimeError("Duplicate original source station key in year "+period)
        # All actual EQR numeric outcomes are deliberately NOT read/retained at source stage
        d[key]={"site_key":key,"method":x.get("procedureClassificationSystem",""),
                "source_bqe_code":x.get("observedPropertyDeterminandBiologyEQRCode",""),
                "country":x["countryCode"],
                "lat":x["lat"],"lon":x["lon"]}
    data[period]=d
same=set(data["2018"])&set(data["2022"])
if len(same)!=1167:raise RuntimeError("Original 2018→22 exact same station identifiers changed")
source_pairs=[]
method_mismatch=collections.Counter()
n_h7=collections.Counter()
for key in sorted(same):
    h18=joined.get(("EEA18",key))
    h22=joined.get(("EEA22",key))
    if h18 is None or h22 is None or h18!=h22:
        method_mismatch["source_not_both_same_H7"]+=1
        continue
    p=data["2018"][key];q=data["2022"][key]
    if p["country"]!=q["country"]:
        method_mismatch["country_changed"]+=1;continue
    if p["method"]!=q["method"]:
        method_mismatch["procedure_method_changed"]+=1;continue
    if p["source_bqe_code"]!=q["source_bqe_code"]:
        method_mismatch["observed_bqe_code_changed"]+=1;continue
    n_h7[h18]+=1
    source_pairs.append({"HYBAS_ID_L07":h18,"original_site_key":key,
        "country":p["country"],"method":p["method"],"observed_bqe_code":p["source_bqe_code"],
        "source_same_station_method_both_2018_2022":True})
n_unique=len(n_h7)
gate=(n_unique>=120 and len(source_pairs)>=500)
out={
 "schema":"H7_234_EU_2018_2022_EXACT_STATION_BQE_METHOD_MATCH_PRE_NEW_CO_RESPONSE_v1",
 "original_L7_source_manifest_SHA256":SITES_SHA,
 "original_234_repeat_basin_WKT_taxonomy_manifest_SHA256":EU234_SHA,
 "original_EEA_Git_Blob_SHA":{k:v[1] for k,v in BLOBS.items()},
 "original_EEA_same_site_keys_all_countries":len(same),
 "original_234_EU_H7_basins_source_selection":len(basins),
 "matched_same_stations_and_same_BQE_method_and_country_n":len(source_pairs),
 "matched_same_stations_and_same_BQE_method_unique_H7_basins":n_unique,
 "n_by_original_country":dict(collections.Counter(x["country"] for x in source_pairs)),
 "n_by_original_2018_22_method":dict(collections.Counter(x["method"] for x in source_pairs)),
 "excluded_reason_source_counts":dict(method_mismatch),
 "source_gate_pass":gate,
 "never_read_numerical_EQR_biological_outcomes":True,
 "new_butterfly_family_or_plant_species_results_NOT_used_for_site_selection":True,
 "not_a_mechanistic_causal_wastewater_effect":True
}
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_234_EU_same_station_source"
root.mkdir(parents=True,exist_ok=True)
manifest=(json.dumps({"schema":out["schema"],"original_source_sha256":SITES_SHA,
                      "source_station_pairs":source_pairs,
                      "no_original_numeric_EQR_read":True},sort_keys=True,separators=(",",":"))+"\n").encode()
out["original_same_station_source_manifest_SHA256"]=hashlib.sha256(manifest).hexdigest()
(root/"H7_EU_234_ORIGINAL_EXACT_BQE_METHOD_SAME_STATION_PAIRS_2018_2022.json").write_bytes(manifest)
(root/"H7_EU_234_ORIGINAL_EXACT_BQE_METHOD_PAIR_SOURCE_QA.json").write_text(json.dumps(out,sort_keys=True,indent=2))
print("H7_234_EU_EXACT_SAME_STATION_SOURCE_QA_START")
print(json.dumps(out,sort_keys=True))
print("H7_234_EU_EXACT_SAME_STATION_SOURCE_QA_END")
if not gate:raise RuntimeError("Native 234 EQR repeat basin same method station gate STOP")
