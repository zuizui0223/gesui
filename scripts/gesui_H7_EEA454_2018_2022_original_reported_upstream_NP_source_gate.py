#!/usr/bin/env python3
"""Original EEA2018/2022 actual ANNUAL WWTP reported point-source N/P
upstream of each monitored river station, aggregated to the fixed 454
real H7 basin×periods. SOURCE ONLY: does not analyze EQR or GBIF outcomes.

All 454 fixed original macro geographic rows are from the independently
successful 232 H7×period true cross-taxon source study #43215754868.
"""
import collections,csv,hashlib,io,json,math,os,pathlib,statistics

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_EEA_2018_2022_NP_SOURCE"
ROOT.mkdir(parents=True,exist_ok=True)
H7_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
BLOBS={
 "EEA18":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
SOURCE_STUDY_SCHEMA="GESUI_HYBAS7_234EU_source_matched_three_response_static_WWTP_gradient_posthoc_v1"
source_report=json.loads(pathlib.Path(os.environ["SOURCE_454_H7_MACRO_REPORT"]).read_text())
if source_report["schema"]!=SOURCE_STUDY_SCHEMA:raise RuntimeError("Original saved 454 H7-period source study changed")
with open(os.environ["SOURCE_454_H7_MACRO_COHORT"],encoding="utf-8-sig",newline="") as f:
    macro=list(csv.DictReader(f))
if len(macro)!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in macro})!=454:
    raise RuntimeError("Original 454 source-fixed Level7 period cohort changed")
basin_periods={(r["HYBAS_ID_L07"],r["PERIOD"]) for r in macro}
h7=pathlib.Path(os.environ["SOURCE_NATIVE_STATION_H7_MANIFEST"])
if hashlib.sha256(h7.read_bytes()).hexdigest()!=H7_SHA:
    raise RuntimeError("Original original 6330 station native L7 source SHA mismatch")
site_data=json.loads(h7.read_text())["candidate_source_site_period_rows"]
if len(site_data)!=6330:raise RuntimeError("Source native stations not original n6330")
site_by_period=collections.defaultdict(list)
for r in site_data:
    if r["panel"] not in BLOBS or r["HYBAS_ID_L07"] is None:continue
    period="2017-2019" if r["panel"]=="EEA18" else "2020-2022"
    k=(r["HYBAS_ID_L07"],period)
    if k in basin_periods:site_by_period[k].append(r["native_id"])
if set(site_by_period)!=basin_periods:
    raise RuntimeError("Some original fixed EEA H7 period missing native station IDs")
def gitsha(blob):return hashlib.sha1(("blob "+str(len(blob))+"\0").encode()+blob).hexdigest()
original={}
for panel,(path,sha) in BLOBS.items():
    b=pathlib.Path(path).read_bytes()
    if gitsha(b)!=sha:raise RuntimeError("Original 2018/2022 source EEA Git blob mismatch")
    rr=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"),newline="")))
    if len(rr)!={"EEA18":2634,"EEA22":2475}[panel]:raise RuntimeError("Original source EEA site count changed")
    original[panel]={r["site_key"]:r for r in rr}
    if len(original[panel])!=len(rr):raise RuntimeError("Duplicate original EEA site key within year")

def num(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None
# 365.25 days/year (same explicit conversion in all original station-year sources).
TONS_YEAR_PER_CMS_TO_MGL=1e6/(365.25*24*60*60)
P_KEYS={"EEA18":"2018_measured_upstream_P_t_per_year","EEA22":"upstream_2022_EEA_P_measured_t_per_yr"}
N_KEYS={"EEA18":"2018_measured_upstream_N_t_per_year","EEA22":"upstream_2022_EEA_N_measured_t_per_yr"}
bins=collections.Counter();positive_country=collections.Counter();positive_basins=set()
rows=[]
for hid,period in sorted(basin_periods):
    panel="EEA18" if period=="2017-2019" else "EEA22"
    record=original[panel]
    sites=site_by_period[(hid,period)]
    station=[];missingflow=0;negative=0
    for site in sites:
        if site not in record:raise RuntimeError("Monitoring station in GIS but not EEA original source")
        z=record[site]
        p=num(z[P_KEYS[panel]]);n=num(z[N_KEYS[panel]]);dis=num(z["DIS_AV_CMS"])
        if p is None or n is None or p<0 or n<0:
            negative+=1;continue
        # Native source environmental flow used as denominator only where >0.
        if dis is None or dis<=0:
            missingflow+=1;continue
        sourceyes=(p>0 or n>0)
        pressureP=p/dis*TONS_YEAR_PER_CMS_TO_MGL
        pressureN=n/dis*TONS_YEAR_PER_CMS_TO_MGL
        if sourceyes:
            positive_country[z["countryCode"]]+=1
            positive_basins.add(hid)
        station.append({
         "original_site_key":site,"P_tpy_reported_upstream":p,"N_tpy_reported_upstream":n,
         "DIS_AV_CMS":dis,"P_mgL_equiv_full_delivery":pressureP,
         "N_mgL_equiv_full_delivery":pressureN,
         "source_reported_positive":sourceyes,
         "source_2018_upstream_registered_WWTP_count":num(z.get("2018_upstream_qualifying_WWTP_count")) if panel=="EEA18" else None,
         "country":z["countryCode"],"method":z["procedureClassificationSystem"]
        })
    if not station:
        bins["no_station_with_valid_reported_NP_and_river_flow"]+=1
        rows.append({"HYBAS_ID_L07":hid,"PERIOD":period,"source_valid":False,
                     "n_original_EEA_stations_raw":len(sites),
                     "n_original_EEA_stations_source_qualified":0})
        continue
    npos=sum(z["source_reported_positive"] for z in station)
    bins["basin_period_with_positive_source"]+=npos>0
    bins["basin_period_with_no_positive_reported_source"]+=npos==0
    bins["source_positive_stations"]+=npos
    bins["source_qualified_stations"]+=len(station)
    if npos and period=="2017-2019":bins["source_positive_basins_period_early"]+=1
    if npos and period=="2020-2022":bins["source_positive_basins_period_late"]+=1
    country_name=collections.Counter(z["country"] for z in station)
    local={
      "HYBAS_ID_L07":hid,"PERIOD":period,
      "source_valid":True,"n_original_EEA_stations_raw":len(sites),
      "n_original_EEA_stations_source_qualified":len(station),
      "source_station_annual_upstream_NP_positive_n":npos,
      "source_station_annual_upstream_NP_positive_fraction":npos/len(station),
      "source_year_reported_upstream_P_tpy_site_equal_mean":statistics.mean(z["P_tpy_reported_upstream"] for z in station),
      "source_year_reported_upstream_N_tpy_site_equal_mean":statistics.mean(z["N_tpy_reported_upstream"] for z in station),
      "source_year_reported_P_mgL_equivalent_site_equal_mean":statistics.mean(z["P_mgL_equiv_full_delivery"] for z in station),
      "source_year_reported_N_mgL_equivalent_site_equal_mean":statistics.mean(z["N_mgL_equiv_full_delivery"] for z in station),
      "source_year_station_DIS_AV_CMS_mean":statistics.mean(z["DIS_AV_CMS"] for z in station),
      "source_source_positive_reported_not_proof_of_only_WWTP":True,
      "source_modelled_river_DIS_AV_CMS_not_reach_delivered_water_chemistry":True,
      "source_2022_actual_matched_upstream_plant_count_not_available":True,
      "source_country_modes":";".join(sorted(country_name,key=lambda k:(-country_name[k],k))),
      "source_country_n":len(country_name),"source_method_n":len({z["method"] for z in station}),
    }
    rows.append(local)
if len(rows)!=454 or len({(x["HYBAS_ID_L07"],x["PERIOD"]) for x in rows})!=454:raise RuntimeError("Original source grain violated")
positive=[r for r in rows if r.get("source_valid") and r["source_station_annual_upstream_NP_positive_n"]>0]
pbasins={z["HYBAS_ID_L07"] for z in positive}
sourcecountries={c for r in positive for c in r["source_country_modes"].split(";") if c}
gate=bool(
 len(positive)>=30 and len(pbasins)>=20
 and sum(z["PERIOD"]=="2017-2019" for z in positive)>=10
 and sum(z["PERIOD"]=="2020-2022" for z in positive)>=10
 and len(sourcecountries)>=3
 and sum(r["source_valid"] for r in rows)>=430
)
out={
 "schema":"GESUI_454_EU_H7_ONE_WATERSHED_TWO_PERIOD_YEAR_SPECIFIC_REPORTED_NP_SOURCE_ONLY_v1",
 "same_fixed_original_454_basin_periods":len(rows),
 "fixed_original_basin_ids":len({r["HYBAS_ID_L07"] for r in rows}),
 "original_station_H7_source_SHA256":H7_SHA,
 "EEA_original_source_GitBlob_SHA1":{a:b[1] for a,b in BLOBS.items()},
 "original_H7_saved_3layer_result_schema":SOURCE_STUDY_SCHEMA,
 "number_with_source_valid_nutrient_flow":sum(r["source_valid"] for r in rows),
 "number_with_positively_matched_reported_upstream_WWTP_load":len(positive),
 "positive_source_unique_H7_basins":len(pbasins),
 "positive_source_H7_period_EARLY":sum(z["PERIOD"]=="2017-2019" for z in positive),
 "positive_source_H7_period_LATE":sum(z["PERIOD"]=="2020-2022" for z in positive),
 "positive_source_distinct_countries":len(sourcecountries),
 "positive_source_original_country_codes":sorted(sourcecountries),
 "same_original_site_source_counts":dict(bins),
 "original_source_positive_station_country_frequency":dict(positive_country),
 "t_per_year_div_modelled_cms_to_mgL_equivalent_conversion":TONS_YEAR_PER_CMS_TO_MGL,
 "source_reported_t_per_year_only_not_actual_receiver_station_TNTP_or_WWTP_causation":True,
 "original_reported_positive_not_all_original_2018_2022_facilities_had_complete_reporting":True,
 "no_2018_2022_EQR_plant_S50_or_butterfly_results_inspected":True,
 "locked_source_only_gate_minimums":{
   "source_valid_basinperiods":430,"positive_basinperiods":30,
   "unique_source_positive_H7_basins":20,
   "positive_EARLY":10,"positive_LATE":10,
   "positive_reported_countries":3},
 "source_only_gate_pass":gate,
 "ecological_NP_coresponse_fit_ALLOWED":gate,
 "if_source_insufficient":"STOP_BEFORE_ECOLOGICAL_NUMERICAL_MODEL, preserve this source summary only"
}
source_manifest={"schema":out["schema"],"original_source_reconciliation_pre_ecology":True,
    "source_H7_period_reported_N_P_upstream_to_true_original_station":rows,
    "true_original_EEA_numeric_EQR_never_used_to_select":True}
raw=(json.dumps(source_manifest,sort_keys=True,separators=(",",":"))+"\n").encode()
outfile=ROOT/"EEA2018_EEA2022_454_H7_ORIGINAL_REPORTED_NP_SITE_FLOW_EQUIVALENT_SOURCE_ONLY.json"
outfile.write_bytes(raw)
out["source_positive_mass_to_H7_manifest_SHA256"]=hashlib.sha256(raw).hexdigest()
(ROOT/"EEA454_H7_2018_2022_ACTUAL_REPORTED_NP_POSITIVE_COVERAGE_PRE_EFFECT_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("EEA_H7_ACTUAL_ANNUAL_NP_SOURCE_ONLY_START",flush=True)
print(json.dumps(out,sort_keys=True),flush=True)
print("EEA_H7_ACTUAL_ANNUAL_NP_SOURCE_ONLY_END",flush=True)
