#!/usr/bin/env python3
"""Create the first REAL common native (HydroBASINS L7 × period) macro table.

Uses source-only original EEA18/22 and US EPA 18/19 basin crosswalk from
GH Actions 43165326244. MISSING response layers stay null, not zero.
No EQR/EPT/plant/butterfly/ecosystem-service effect signs examined.
"""
import collections,csv,hashlib,json,os,pathlib
IN=pathlib.Path(os.environ["HYBAS7_FROZEN_NATIVE_SOURCE_MANIFEST"])
SOURCE_SHA256="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
b=IN.read_bytes()
if hashlib.sha256(b).hexdigest()!=SOURCE_SHA256:
    raise RuntimeError("Original HydroBASINS native river station mapping SHA mismatch")
meta=json.loads(b)
source=meta["candidate_source_site_period_rows"]
if len(source)!=6330:raise RuntimeError("Original 6330 source river-visit entries missing")
areas=collections.defaultdict(set);upareas=collections.defaultdict(set)
basins=set();bininfo=collections.defaultdict(lambda:{
    "eu_station_keys":set(),"us_station_keys":set(),"us_Uids":set(),
    "countries":set(),"source_types":set()})
for row in source:
    hid=row["HYBAS_ID_L07"]
    if hid is None:continue
    binid=(hid,row["period"])
    basins.add(hid)
    areas[hid].add(float(row["sub_area_km2"]))
    upareas[hid].add(float(row["up_area_km2"]))
    u=bininfo[binid]
    u["countries"].add(row["country"])
    u["source_types"].add(row["panel"])
    if row["panel"] in ("EEA18","EEA22"):u["eu_station_keys"].add(row["physical_id"])
    else:
        u["us_station_keys"].add(row["physical_id"])
        u["us_Uids"].add(row["native_id"])
if len(basins)<1000:raise RuntimeError("Unexpectedly small global Europe+USA basin cohort")
if any(len(q)!=1 for q in areas.values()) or any(len(q)!=1 for q in upareas.values()):
    raise RuntimeError("Same HydroBASINS ID has inconsistent source polygon area")
periods=["2017-2019","2020-2022"]
out=[]
for hid in sorted(basins):
    for period in periods:
        o=bininfo.get((hid,period))
        ee=len(o["eu_station_keys"]) if o is not None else 0
        us=len(o["us_station_keys"]) if o is not None else 0
        visits=len(o["us_Uids"]) if o is not None else 0
        out.append({
          "HYBAS_ID_L07":hid,"PERIOD":period,
          "LOCAL_SUBBASIN_AREA_KM2":next(iter(areas[hid])),
          "UPSTREAM_TOTAL_AREA_KM2":next(iter(upareas[hid])),
          "OBSERVED_COUNTRY_OR_STATE_CODES":";".join(sorted(o["countries"])) if o else "",
          "EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N":ee,
          "EEA_BENTHIC_EQR_SOURCE_SUPPORT":"observed_2018_or_2022" if ee else "not_measured_here",
          "EPA_NRSA_BENTHIC_PHYSICAL_STATIONS_SOURCE_N":us,
          "EPA_NRSA_NATIVE_UID_VISITS_SOURCE_N":visits,
          "EPA_NRSA_BENTHIC_SOURCE_SUPPORT":"observed_2018_2019" if us else "not_measured_here",
          "WWTP_EU_MEASURED_UPSTREAM_ANNUAL_N_TPY":None,
          "WWTP_EU_MEASURED_UPSTREAM_ANNUAL_P_TPY":None,
          "WWTP_GLOBAL_HYDROWASTE_STATIC_EFFLUENT_FLOW":None,
          "RIVER_EQR_WITHIN_METHOD_STANDARDIZED":None,
          "RIVER_EPT_FRACTION_WITHIN_SOURCE":None,
          "PLANTAE_GBIF_2017_2019_OR_2020_2022_RECORDING_QC":None,
          "PLANTAE_RAREFIED_SPECIES_RICHNESS":None,
          "BUTTERFLY_RECORDED_SOURCE_EFFORT":None,
          "LYCAENIDAE_OBSERVED_COMPOSITION":None,
          "MEASURED_EMERGENT_ADULT_BIOMASS":None,
          "REALIZED_POLLINATION_OR_PREDATOR_DIET":None,
          "CROSS_TAXON_COMPLETE_CASE_NOW":False
        })
summary={
 "schema":"GESUI_HYBAS7_X_2PERIOD_GLOBAL_MULTILAYER_EMPTY_SOURCE_FIRST_PANEL_v1",
 "original_pre_source_HYBAS7_geography_sha256":SOURCE_SHA256,
 "H7_basin_count":len(basins),"periods":periods,"rows_total":len(out),
 "basin_period_rows_with_EEA_monitoring":sum(r["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"]>0 for r in out),
 "basin_period_rows_with_EPA_monitoring":sum(r["EPA_NRSA_BENTHIC_PHYSICAL_STATIONS_SOURCE_N"]>0 for r in out),
 "basin_period_rows_with_either_source":sum(r["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"]>0 or r["EPA_NRSA_BENTHIC_PHYSICAL_STATIONS_SOURCE_N"]>0 for r in out),
 "basin_period_rows_with_both_US_EPT_and_EU_EQR":sum(r["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"]>0 and r["EPA_NRSA_BENTHIC_PHYSICAL_STATIONS_SOURCE_N"]>0 for r in out),
 "basins_with_EEA_both_periods":sum(
   any(r["HYBAS_ID_L07"]==b and r["PERIOD"]==period and r["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"]>0 for r in out for period in []) for b in []),
 "basin_period_with_comparable_plant_butterfly_response_from_original_coastal_ports":0,
 "all_outcome_cells_NULL_not_true_biological_zero":True,
 "source_nonnumeric_YET":True,
 "all_coefficients_unestimated":True,
 "hydrological_upstream_loads_NOT_computed":True,
 "source_time_or_sampling_bias_NOT_yet_validated":True,
 "missing_P2_USA_is_structural_no_observation_NOT_extinction":True,
 "scientific_next_gate":"Get native basin-period original GBIF plant/host and butterfly occurrence with effort normalization; route HYDROWASTE/EEA reported source along rivers. Only then compute common-case response-layer models."
}
# Fix explicit repeated basin count separately without hidden filtering.
eubys=collections.defaultdict(set)
for row in out:
    if row["EEA_BENTHIC_NUMERIC_EQR_ORIGINAL_STATIONS_SOURCE_N"]>0:
        eubys[row["HYBAS_ID_L07"]].add(row["PERIOD"])
summary["basins_with_EEA_both_periods"]=sum(len(z)==2 for z in eubys.values())
assert summary["basins_with_EEA_both_periods"]==234
assert summary["rows_total"]==summary["H7_basin_count"]*2
assert summary["basin_period_rows_with_both_US_EPT_and_EU_EQR"]==0
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HYBAS7_master"
ROOT.mkdir(parents=True,exist_ok=True)
csv_path=ROOT/"HYBAS7_2017_2019_2020_2022_ALL_RESPONSE_LAYERS_STATUS_MASTER_SOURCE_ONLY.csv"
with csv_path.open("w",encoding="utf8",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(out[0]))
    writer.writeheader();writer.writerows(out)
summary["source_only_master_csv_sha256"]=hashlib.sha256(csv_path.read_bytes()).hexdigest()
(ROOT/"HYBAS7_2PERIOD_GLOBAL_AQUATIC_PLANTS_BUTTERFLY_SOURCE_TABLE_QA.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
print("MACRO_HYBAS7_TIME_LAYER_INDEX_RESULT_START",flush=True)
print(json.dumps(summary,sort_keys=True),flush=True)
print("MACRO_HYBAS7_TIME_LAYER_INDEX_RESULT_END",flush=True)
