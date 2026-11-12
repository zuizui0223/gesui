#!/usr/bin/env python3
"""Merge original source-locked 1219 H7 x2 time bins with original HydroWASTE
static outlet pressure NOW, without waiting for butterfly/plant GBIF records.

Each 2015-era static HydroWASTE pressure is labelled STATIC and repeated in
both source periods; NOT 2018-vs-2022 change or delivered dose to station.
"""
import collections,csv,hashlib,json,os,pathlib
MASTER_SHA="a9a8eaa5c328026e46118810f9fbde29cc24451a594ad9437ba2222201bf5e0a"
WWTP_SHA="a311127b04e0a081a32a764b183240c8642aa8e4e71f81016bde9a56e537fcac"
master=pathlib.Path(os.environ["HYBAS7_ORIGINAL_MASTER"])
wwtp=pathlib.Path(os.environ["HYDROWASTE_ORIGINAL_H7_SOURCE"])
if hashlib.sha256(master.read_bytes()).hexdigest()!=MASTER_SHA or hashlib.sha256(wwtp.read_bytes()).hexdigest()!=WWTP_SHA:
    raise RuntimeError("Original macro source or modeled HydroWASTE pressure changed")
def table(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
a=table(master);b=table(wwtp)
assert len(a)==2438 and len(b)==1219
d={x["HYBAS_ID_L07"]:x for x in b}
if len(d)!=1219:raise RuntimeError("Original HydroWASTE one-to-one Level7 basin key broken")
x=[]
for row in a:
    k=row["HYBAS_ID_L07"]
    if k not in d:raise RuntimeError("Original river-monitored basin missing real HydroWASTE model")
    w=d[k]
    z=dict(row)
    for col in [
      "HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT",
      "HYDROWASTE_STATIC_LOCAL_MODELLED_WASTEWATER_DISCHARGE_M3_PER_DAY",
      "HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT",
      "HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY"
    ]:
        if col in z:raise RuntimeError("Do not overwrite original source response")
        z[col]=w[col]
    z["HYDROWASTE_PRESSURE_SOURCE_YEAR"]="STATIC_ORIGINAL_V1_NOT_YEAR_RESOLVED"
    z["HYDROWASTE_OUTLET_PRESSURE_IS_MONITOR_STATION_UPSTREAM_DOSE"]=False
    z["SOURCE_ACTUAL_2018_22_NP_MASS_MISSING"]=True
    x.append(z)
assert len(x)==2438 and len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in x})==2438
p1={k:r for r in x for k in [(r["HYBAS_ID_L07"],r["PERIOD"])]}
for hid in d:
    q=p1[(hid,"2017-2019")]
    r=p1[(hid,"2020-2022")]
    if q["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY"]!=r["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_MODELED_WASTE_DIS_M3_PER_DAY"]:
        raise RuntimeError("Static modelled wastewater source differs by year incorrectly")
outdir=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HYBAS7_2PERIOD_STATIC_HYDROWASTE"
outdir.mkdir(parents=True,exist_ok=True)
csvpath=outdir/"HYBAS7_1219BASINS_2PERIOD_REAL_RIVER_MONITORING_STATIC_WWTP_PRESSURE_MASTER.csv"
with csvpath.open("w",encoding="utf8",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(x[0]))
    writer.writeheader();writer.writerows(x)
source={
 "schema":"gesui_1219_original_H7_basin_by_2period_true_HydroWASTE_static_join_v1",
 "original_aquatic_1219x2_master_SHA256":MASTER_SHA,
 "original_1219_HYBAS7_source_WWTP_SHA256":WWTP_SHA,
 "whole_row_count":len(x),"unique_HYBAS_IDs":len(d),
 "unique_H7_times_period_pairs":len(x),
 "valid_WWTP_outlet_local_basin_count":sum(int(v["HYDROWASTE_STATIC_LOCAL_REGISTERED_WWTP_COUNT"])>0 for v in b),
 "valid_WWTP_outlet_upstream_basins_count":sum(int(v["HYDROWASTE_STATIC_UPSTREAM_TO_BASIN_OUTLET_REGISTERED_WWTP_COUNT"])>0 for v in b),
 "EU_real_benthic_monitor_basin_periods":sum(int(q["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"])>0 for q in x),
 "US_real_benthic_monitor_basin_periods":sum(int(q["EPA_NRSA_BENTHIC_PHYSICAL_STATIONS_SOURCE_N"])>0 for q in x),
 "still_unavailable": ["GBIF basin×period plants","GBIF basin×period Lycaenidae","GBIF baseline survey opportunity taxa","accepted plant species facets","butterfly standardized monitoring", "actual source-specific WWTP nutrient concentrations per station", "realized pollination or adult emergence per basin"],
 "hydrowaste_upstream_to_outlet_not_to_any_EEA_or_USA_interior_station":True,
 "source_year_static_repeated_on_two_rows_no_temporal_change":True,
 "no_ecological_response_model_fitted":True,
 "source_only_master_SHA256":hashlib.sha256(csvpath.read_bytes()).hexdigest()
}
(outdir/"HYBAS7_ORIGINAL_2438_AQUATIC_STATIC_WWTP_MASTER_QA.json").write_text(json.dumps(source,indent=2,sort_keys=True))
print("HYBAS7_2438_STATIC_WWTP_MASTER_RESULT_START")
print(json.dumps(source,sort_keys=True))
print("HYBAS7_2438_STATIC_WWTP_MASTER_RESULT_END")
