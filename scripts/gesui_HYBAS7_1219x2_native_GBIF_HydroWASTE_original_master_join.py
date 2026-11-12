#!/usr/bin/env python3
"""Source-only assemble authentic Level7 × two-period master with all 1219
basin WWTP outlet pressure and real GBIF species-GROUP record counts.

Check exactly 4 immutable source shards, original source URLs/SHA, original
1219 study basins ×2 time bins. No ecological regression or species richness.
"""
import collections,csv,hashlib,json,os,pathlib
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_full_multitaxon_master"
ROOT.mkdir(parents=True,exist_ok=True)
MASTER_SHA="a9a8eaa5c328026e46118810f9fbde29cc24451a594ad9437ba2222201bf5e0a"
WWTP_SHA="a311127b04e0a081a32a764b183240c8642aa8e4e71f81016bde9a56e537fcac"
WKT_SHA="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
masterfile=pathlib.Path(os.environ["ORIGINAL_H7_MASTER"])
wastefile=pathlib.Path(os.environ["ORIGINAL_H7_WWTP_PRESSURE"])
if hashlib.sha256(masterfile.read_bytes()).hexdigest()!=MASTER_SHA:
    raise RuntimeError("Original H7 1219 × 2 source-only master SHA changed")
if hashlib.sha256(wastefile.read_bytes()).hexdigest()!=WWTP_SHA:
    raise RuntimeError("Original HydroWASTE 1219 static watershed source CSV changed")
def reader(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
master=reader(masterfile);wwtp=reader(wastefile)
if len(master)!=2438 or len(wwtp)!=1219:raise RuntimeError("Fixed source count unexpected")
by_waste={x["HYBAS_ID_L07"]:x for x in wwtp}
if len(by_waste)!=1219:raise RuntimeError("Repeated local WWTP basin source")
runs=[]
by={}
shard_sizes={}
for i in range(4):
    folder=pathlib.Path(os.environ["GBIF_H7_SHARD_"+str(i)])
    rfile=folder/"HYBAS7_GBIF_ORIGINAL_SOURCE_COUNTS_URL_SHA.json"
    qfile=folder/"HYBAS7_GBIF_ORIGINAL_SOURCE_SUPPORT_QA.json"
    data=json.loads(rfile.read_text())
    qa=json.loads(qfile.read_text())
    if qa["source_wkt_manifest_SHA256"]!=WKT_SHA or qa["mode"]!="full" or qa["shard"]!=i:
        raise RuntimeError("Full original GBIF shard source identity incorrect")
    if len(data)!=qa["fixed_source_requests"] or len(data)!=6*qa["fixed_source_shard_basin_count"]:
        raise RuntimeError("Original shard request count inconsistent")
    shard_sizes[str(i)]=qa["fixed_source_shard_basin_count"]
    runs.append(qa)
    for x in data:
        hid=x["HYBAS_ID_L07"];time=x["period"];tax=x["taxon"]
        k=(hid,time,tax)
        if k in by:raise RuntimeError("Same original H7×period×taxon input appears in several source shards")
        by[k]=x
if len(by)!=7314 or sum(shard_sizes.values())!=1219:
    raise RuntimeError("GBIF request missing from original 1219 ×2 ×3 source cohort")
taxa=["Plantae","Lycaenidae","Animalia"]
periods=["2017-2019","2020-2022"]
complete=0
well=0
countries=collections.Counter();states=set();both=collections.defaultdict(set)
source_errors=collections.Counter()
counts={q:collections.Counter() for q in taxa}
out=[]
for original in master:
    hid=original["HYBAS_ID_L07"];period=original["PERIOD"]
    if period not in periods or hid not in by_waste:raise RuntimeError("Master row not in original source")
    row=dict(original)
    for key in ["HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT",
        "HYDROWASTE_STATIC_LOCAL_MODELLED_WASTEWATER_DISCHARGE_M3_PER_DAY",
        "HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT",
        "HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY"]:
        row[key]=by_waste[hid][key]
    vals={}
    for tax in taxa:
        x=by[(hid,period,tax)]
        if x["status"]=="ok":
            c=x["GBIF_index_record_count"]
            if not isinstance(c,int) or c<0:raise RuntimeError("Genuine GBIF index count invalid")
            row["GBIF_INDEX_"+tax.upper()+"_2017_2022_FIXED_PERIOD_RECORD_COUNT"]=c
            vals[tax]=c
            counts[tax]["observed_source_queries"]+=1
            counts[tax]["true_numeric_zero_record_index"]+=c==0
        else:
            row["GBIF_INDEX_"+tax.upper()+"_2017_2022_FIXED_PERIOD_RECORD_COUNT"]=""
            vals[tax]=None
            source_errors[tax]+=1
        row["GBIF_"+tax.upper()+"_SOURCE_PRESENT_FLAG"]=int(x["status"]=="ok")
    all_present=all(v is not None for v in vals.values())
    row["GBIF_ALL3_SOURCE_PRESENT_FLAG"]=int(all_present)
    good=all_present and vals["Plantae"]>=50 and vals["Lycaenidae"]>=10 and vals["Animalia"]>=50
    row["SOURCE_INDEXED_PLANT50_LYC10_ANIMAL50_EFFORT_TIER"]=int(good)
    row["NEW_ORIGINAL_SPECIES_RICHNESS_NOT_YET_COMPUTED"]=True
    row["REAL_POLLINATION_OR_OUTLET_BIOLOGICAL_SOURCE_NOT_MEASURED"]=True
    if all_present:
        complete+=1
        both[hid].add(period)
    if good:
        well+=1
        codes=row["OBSERVED_COUNTRY_OR_STATE_CODES"].split(";")
        for code in codes:
            if code:
                countries[code[:2] if code.startswith("US-") else code]+=1
                if code.startswith("US-"):states.add(code[3:])
    out.append(row)
assert len(out)==2438
two=len([x for x,v in both.items() if len(v)==2])
passgate=(complete>=.90*2438 and well>=200 and len(countries)>=8
          and len({z["HYBAS_ID_L07"] for z in out if z["SOURCE_INDEXED_PLANT50_LYC10_ANIMAL50_EFFORT_TIER"]})>=50
          and two>=100)
source_summary={
 "schema":"GESUI_HYBAS7_1219_GLOBAL_MULTITAXA_GBIF_COUNT_AND_STATIC_WWTP_MASTER_v1",
 "original_basin_period_master_SHA256":MASTER_SHA,
 "original_HydroWASTE1219_crosswalk_SHA256":WWTP_SHA,
 "original_HYBAS7_true_polygon_WKT_SHA256":WKT_SHA,
 "original_full_GBIF_shard_QA":runs,
 "original_1219_four_shard_sizes":shard_sizes,
 "original_source_1219_basin_count":1219,
 "original_two_periods":periods,
 "original_expected_basin_periods":2438,
 "original_expected_source_taxon_queries":7314,
 "source_complete_all_three_taxa_basin_periods":complete,
 "source_complete_all_three_taxa_two_periods_same_basin":two,
 "source_plant50_lyc10_animal50_qualifying_basin_periods":well,
 "source_eligible_unique_level7_basins":len({z["HYBAS_ID_L07"] for z in out if z["SOURCE_INDEXED_PLANT50_LYC10_ANIMAL50_EFFORT_TIER"]}),
 "source_eligible_distinct_countries":len(countries),
 "source_eligible_US_states":len(states),
 "source_qualifying_country_codes":dict(countries),
 "GBIF_source_missing_by_taxon":dict(source_errors),
 "GBIF_original_source_count_QA":{k:dict(c) for k,c in counts.items()},
 "original_source_GATE_pass":passgate,
 "predeclared_gate":{"all3_complete_at_least_fraction":.90,
     "plants50_lycae10_anim50_basinperiods":200,"minimum_countries":8,
     "minimum_unique_basins":50,"minimum_two_period_complete_basins":100},
 "pure_recording_opportunity_counts_NO_species_richness_measured":True,
 "future_full_plant_speciesKey_facet_and_Papilionoidea_background_required":True,
 "no_ecological_coefficients_yet":True,
 "WWTP_pressure_is_static_estimated_basin_outlet_not_EEA_EPA_station_dose":True
}
outfile=ROOT/"HYBAS7_1219x2_GBIF_PLANTAE_LYCAENIDAE_ANIMALIA_HYDROWASTE_SOURCE_ONLY_MASTER.csv"
with outfile.open("w",newline="",encoding="utf8") as f:
    writer=csv.DictWriter(f,fieldnames=list(out[0]))
    writer.writeheader();writer.writerows(out)
source_summary["new_master_csv_SHA256"]=hashlib.sha256(outfile.read_bytes()).hexdigest()
(ROOT/"HYBAS7_1219x2_ORIGINAL_GBIF_WWTP_SOURCE_SUPPORT_QA.json").write_text(json.dumps(source_summary,indent=2,sort_keys=True))
print("HYBAS7_ALL_TAXA_WWTP_MASTER_RESULT_START",flush=True)
brief={k:v for k,v in source_summary.items() if k not in ("original_full_GBIF_shard_QA","source_qualifying_country_codes")}
print(json.dumps(brief,sort_keys=True),flush=True)
print("HYBAS7_ALL_TAXA_WWTP_MASTER_RESULT_END",flush=True)
if not passgate:
    print("SOURCE STOP: original GBIF coverage gate not met; no observational ecosystem layer regression permitted",flush=True)
