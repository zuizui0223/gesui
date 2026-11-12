#!/usr/bin/env python3
"""Outcome blind source reconciliation of 4 authentic original source shards
(234 H7 repeated EEA catchments ×2 source windows × genuine plant S50 and five
butterfly family counts) with the independently frozen original 1047 same-
station/method EQR source key pair. No ecology associations yet.
"""
import collections,csv,hashlib,json,math,os,pathlib,statistics
SOURCE_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189"
STATION_SHA="373a3a96ea790a7592cc4d10b0e4f9244fbdddab2c1af4901118792674874930"
MACRO_SHA="3f95d4825455cfab67a41de3ec4dab8778acf97d2804f2148108b65c95225c52"
originalfile=pathlib.Path(os.environ["SOURCE_EU234_COLORECT_NATIVE"])
stationfile=pathlib.Path(os.environ["SOURCE_EU_SAME_METHOD_STATIONS"])
macrofile=pathlib.Path(os.environ["SOURCE_H7_1219X2_GBIF_WWTP"])
def verify(path,sha):
    if hashlib.sha256(path.read_bytes()).hexdigest()!=sha:
        raise RuntimeError("Original immutable source file SHA256 changed: "+str(path))
verify(originalfile,SOURCE_SHA)
verify(stationfile,STATION_SHA)
verify(macrofile,MACRO_SHA)
orig=json.loads(originalfile.read_text())
stations=json.loads(stationfile.read_text())["source_station_pairs"]
if len(stations)!=1047 or len({r["original_site_key"] for r in stations})!=1047:
    raise RuntimeError("Native EEA identical-site/assessment-method 2018/2022 source freeze drift")
source_basin={x["HYBAS_ID_L07"] for x in orig["original_two_year_EEA_H7_subbasins"]}
if len(source_basin)!=234:raise RuntimeError("Original 234 common river basins changed")
rows=[]
shard_counts={}
for i in range(4):
    directory=pathlib.Path(os.environ["SOURCE_EU234_FACETS_SHARD_"+str(i)])
    r=directory/"HYBAS7_234EU_NATIVE_PLANT_S50_BUTTERFLY5_ORIGINAL_SOURCE_RESULTS.json"
    q=directory/"HYBAS7_234EU_NATIVE_SOURCE_JOINT_ECOLOGY_QA.json"
    z=json.loads(r.read_text()); qa=json.loads(q.read_text())
    if qa["fixed_original_source_and_taxon_manifest_sha256"]!=SOURCE_SHA or qa["shard"]!=i or qa["mode"]!="full":
        raise RuntimeError("Wrong original 234 source shard")
    if len(z)!=qa["original_expected_field_queries"] or len(z)!=qa["original_expected_fixed_basin_periods"]*6:
        raise RuntimeError("Original source shard size changed")
    shard_counts[str(i)]={"basins":qa["original_fixed_EU_basins_this_shard"],
      "events":qa["original_expected_fixed_basin_periods"],"complete_plant":qa["source_facet_plant_complete"],
      "S50_eligible":qa["source_plant_rarefied_50_eligible"],
      "all_five_family":qa["all_five_butterfly_families_request_complete"],
      "eligible_joint":qa["joint_plant_S50_and_butterfly20_genuine_source_support"]}
    rows+=z
if len(rows)!=468*6 or sum(v["basins"] for v in shard_counts.values())!=234:
    raise RuntimeError("Original 468x6 target source not fully obtained")
ind=collections.defaultdict(dict)
for x in rows:
    key=(x["HYBAS_ID_L07"],x["period"])
    fld=x["source_field"]
    if fld in ind[key]:raise RuntimeError("Duplicate original basin-period-facet key")
    ind[key][fld]=x
if len(ind)!=468 or {k[0] for k in ind}!=source_basin:
    raise RuntimeError("Unexpected original river basin identifier")
def csvrows(path):
    with path.open(encoding="utf8",newline="") as f:return list(csv.DictReader(f))
macro=csvrows(macrofile);macromap={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in macro}
if len(macromap)!=2438:raise RuntimeError("Original full biological source macro table changed")
by_station=collections.defaultdict(list)
for q in stations:by_station[q["HYBAS_ID_L07"]].append(q)
panel=[]
missing=collections.Counter();counts=collections.Counter();hist=[]
for (h,period),val in sorted(ind.items()):
    pl=val["Plantae_speciesKey_facet"]
    families=["Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae"]
    booleans=[val[f]["status"]=="ok" for f in families]
    butterflies=sum(val[f]["source_record_index_count"] for f in families) if all(booleans) else None
    lyc=val["Lycaenidae"]["source_record_index_count"] if val["Lycaenidae"]["status"]=="ok" else None
    rich=pl.get("species_50_draw_rarefied_richness")
    eligible=bool(rich is not None and butterflies is not None and butterflies>=20)
    if not eligible:
        if rich is None:missing["plant_S50_missing_or_insufficient"]+=1
        if butterflies is None:missing["some_5_family_API_missing"]+=1
        elif butterflies<20:missing["fewer_than_20_butterfly_family_records"]+=1
    m=macromap[(h,period)]
    if int(m["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"])<=0:
        raise RuntimeError("A source-frozen repeated EQR basin has no original recorded EEA study station")
    row={
      "HYBAS_ID_L07":h,"PERIOD":period,
      "ORIGINAL_STATION_METHOD_REPEAT_ELIGIBLE":bool(h in by_station),
      "N_ORIGINAL_SAME_STATION_METHOD_EEA2018_2022":len(by_station[h]),
      "EEA_SOURCE_MONITOR_STATION_N_IN_THIS_PERIOD":int(m["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"]),
      "source_EEA_country_or_state_codes":m["OBSERVED_COUNTRY_OR_STATE_CODES"],
      "original_source_plant_S50":rich,
      "original_source_plant_S50_eligible":rich is not None,
      "original_source_plant_speciesKey_count":pl.get("indexed_original_speciesKey_distinct"),
      "original_source_plant_facetted_record_coverage_fraction":pl.get("speciesKey_index_record_coverage_fraction"),
      "original_source_5family_butterfly_records":butterflies,
      "original_source_lycaenidae_records":lyc,
      "original_source_lycaenidae_share_of_five_families":lyc/butterflies if butterflies and lyc is not None else None,
      "original_5family_plant_S50_and_butterfly20_source_eligible":eligible,
      "original_GBIF_Animalia_record_opportunity":int(m["GBIF_INDEX_ANIMALIA_2017_2022_FIXED_PERIOD_RECORD_COUNT"]) if m["GBIF_INDEX_ANIMALIA_2017_2022_FIXED_PERIOD_RECORD_COUNT"] else None,
      "original_GBIF_Plantae_record_opportunity":int(m["GBIF_INDEX_PLANTAE_2017_2022_FIXED_PERIOD_RECORD_COUNT"]) if m["GBIF_INDEX_PLANTAE_2017_2022_FIXED_PERIOD_RECORD_COUNT"] else None,
      "registered_static_local_WWTP_count":int(m["HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT"]),
      "registered_static_upstream_outlet_WWTP_count":int(m["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT"]),
      "registered_static_WWTP_outlet_modelled_flow_m3perday":float(m["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY"]),
      "source_local_area_km2":float(m["LOCAL_SUBBASIN_AREA_KM2"]),
      "source_UP_AREA_KM2":float(m["UPSTREAM_TOTAL_AREA_KM2"]),
      "true_EEA_EOL_EQR_numeric_NOT_OPENED_YET":True,
      "plant_species_facets_not_direct_field_RICHNESS":True,
      "butterfly_5_family_index_count_not_abundance_or_pollination":True
    }
    panel.append(row)
    counts["all_complete_plant_facets"]+=pl["status"]=="ok"
    counts["all_complete_5_butterfly_families"]+=all(booleans)
    counts["plant_S50_computable"]+=rich is not None
    counts["all_response_source_eligible"]+=eligible
    if rich is not None:hist.append(rich)
bypair=collections.defaultdict(dict)
for x in panel:bypair[x["HYBAS_ID_L07"]][x["PERIOD"]]=x
pairs=[q for q in bypair.values() if all(q[p]["original_5family_plant_S50_and_butterfly20_source_eligible"] for p in ("2017-2019","2020-2022"))]
method_basins=sum(h in by_station for h in source_basin)
method_eligible=sum(h in by_station for h,q in bypair.items() if all(q[p]["original_5family_plant_S50_and_butterfly20_source_eligible"] for p in ("2017-2019","2020-2022")))
# Countries must be confirmed from original native monitoring station metadata,
# source IDs only; no EQR outcome or within-country effect.
country_code=collections.Counter()
for h,q in bypair.items():
    if not all(q[p]["original_5family_plant_S50_and_butterfly20_source_eligible"] for p in ("2017-2019","2020-2022")):continue
    for c in set(";".join(q[p]["source_EEA_country_or_state_codes"] for p in q).split(";")):
        if c:country_code[c]+=1
sourcegate=(counts["all_response_source_eligible"]>=300 and len(pairs)>=150 and len(country_code)>=12
   and len(stations)>=80 and method_basins>=100)
out={
 "schema":"GESUI_H7_234EU_original_S50_5family_EEA_SOURCE_ONLY_FINAL_QA_v1",
 "original_234_EU_WKT_and_taxonomy_manifest_SHA256":SOURCE_SHA,
 "original_EEA_same_method_physical_monitor_source_SHA256":STATION_SHA,
 "original_all_H7_source_GBIF_HydroWASTE_master_SHA256":MACRO_SHA,
 "source_shards":shard_counts,
 "source_basins":len(source_basin),
 "source_periods":len(ind),
 "new_raw_original_plant_facets_and_butterfly_queries":len(rows),
 "source_coverage":dict(counts),
 "source_missing_QA":dict(missing),
 "source_basins_both_periods_S50_and_5butterfly_eligible":len(pairs),
 "source_basins_both_periods_S50_and_5butterfly_and_original_exact_EQR_method_eligible":method_eligible,
 "source_original_same_physical_station_same_EQR_BQE_method_basins":method_basins,
 "source_original_EEA_country_codes_in_eligible_two_year_basins":dict(country_code),
 "source_unique_countries_in_eligible_two_periods":len(country_code),
 "source_predeclared_gate":{"source_eligible_468_periods_at_least":300,
                           "both_period_source_eligible_basins_at_least":150,
                           "different_EU_countries_at_least":12,
                           "same_method_stations_at_least":80},
 "source_predeclared_GATE_PASS":sourcegate,
 "source_predeclared_ecological_regression_ALLOWED":sourcegate,
 "reason_source_only_stop_if_failed":"Original source gate must not be relaxed after looking at plant/5family records, even if observed count tiers are numerically favorable",
 "no_EQR_or_environmental_butterfly_ecological_effect_fitted":True
}
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_234EU_SOURCE_RECONCILIATION"
root.mkdir(parents=True,exist_ok=True)
p=root/"EU234_H7_TRUE_ORIGINAL_PLANT_S50_BUTTERFLY5_EEA_WWTP_SOURCE_ONLY_CANDIDATES.json"
payload=(json.dumps({"source_stage":out["schema"],"original_source_eligible_basin_period_candidates":panel,
    "source_only_no_EQR_numeric_used":True},sort_keys=True,separators=(",",":"))+"\n").encode()
p.write_bytes(payload)
out["source_fixed_234_H7_source_candidate_SHA256"]=hashlib.sha256(payload).hexdigest()
(root/"EU234_H7_PLANT_S50_BUTTERFLY5_EEA_SOURCE_JOIN_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("H7_EU234_SOURCE_RECONCILIATION_START")
print(json.dumps(out,sort_keys=True))
print("H7_EU234_SOURCE_RECONCILIATION_END")
