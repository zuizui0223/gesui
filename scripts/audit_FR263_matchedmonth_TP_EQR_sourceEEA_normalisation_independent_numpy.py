#!/usr/bin/env python3
"""Independent manual NumPy Ridge numeric verification of original French
263 identical 2018/2022 TP+I2M2 sites, 73 H7, 24 H6, month-matched
TP measurement sensitivity; plus EEA resultNormalisedEQRValue source lineage.

No original biological source re-selection; all analyses POSTHOC.
"""
import collections,csv,hashlib,io,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR263_TP_EQR_MONTH_INDEPENDENT"
ROOT.mkdir(parents=True,exist_ok=True)
report=json.loads(pathlib.Path(os.environ["ORIGINAL_FR263_MONTHMATCHED_RESULT"]).read_text())
if report["schema"]!="GESUI_FR263_REAL_TP_2018_22_MONTHMATCHED_I2M2_EQR_73H7_24H6_v1":
    raise RuntimeError("Original month-matched main study source result type changed")
with open(os.environ["ORIGINAL_FR73_H7_MONTHMATCHED_DATA"],encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
with open(os.environ["ORIGINAL_FR263_STATION_MONTHMATCHED_DATA"],encoding="utf-8-sig",newline="") as f:sites=list(csv.DictReader(f))
if len(rows)!=73 or len(sites)!=263 or \
    len({r["HYBAS_ID_L07"] for r in rows})!=73 or len({r["HYBAS_ID_L06"] for r in rows})!=24:
    raise RuntimeError("Original French n263 same-site original 73 H7 source grain changed")
oldsha={
 "2018":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "2022":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
def gitsha(b):return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()
orig={}
for year,(path,sha) in oldsha.items():
    b=pathlib.Path(path).read_bytes()
    if gitsha(b)!=sha:raise RuntimeError("Original EEA "+year+" physical station source SHA changed")
    R=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"),newline="")))
    if len(R)!=({"2018":2634,"2022":2475}[year]):
        raise RuntimeError("Original EEA source national station-year counts changed")
    orig[year]={r["site_key"]:r for r in R}
    if len(orig[year])!=len(R):raise RuntimeError("Duplicate EEA station in "+year)
# Original two committed EEA pipeline scripts clearly use the same original
# EEA 'resultNormalisedEQRValue' as both *2018 numeric* and *2022 normalized*.
script18=pathlib.Path("scripts/gesui_2018_independent_bqe_EQR_outcomes.py").read_text()
script22=pathlib.Path("scripts/gesui_model_2022_eea_aquatic_EQR.py").read_text()
if "data[\"numeric_EQR_2018\"]=observed.to_numpy" not in script18 or \
   "observed=pd.to_numeric(data.resultNormalisedEQRValue" not in script18 or \
   "df[\"normalised_EQR\"]=pd.to_numeric(df.resultNormalisedEQRValue" not in script22:
    raise RuntimeError("Original 2018 or 2022 score source does not use normalized EEA EQR!")
if 'data["raw_EQR_2018"]=raw.to_numpy' not in script18 or \
    "raw=pd.to_numeric(data.resultEQRValue" not in script18:
    raise RuntimeError("2018 raw EEA EQR provenance changed; cannot claim original raw not used")
method_count=collections.Counter()
original_raw2018_diff=0;original2018_raw_missing=0
for r in sites:
    id=r["EEA_identical_station_key"]
    a=orig["2018"].get(id);b=orig["2022"].get(id)
    if not a or not b:raise RuntimeError("Previously source-frozen original EEA physical station ID missing from EEA site-year file")
    if a["countryCode"]!="FR" or b["countryCode"]!="FR" or \
       a["procedureClassificationSystem"]!="I2M2" or b["procedureClassificationSystem"]!="I2M2":
        raise RuntimeError("Original source EEA normalized score procedure/country changed")
    if a["HYRIV_ID"]!=b["HYRIV_ID"]:raise RuntimeError("Original EEA native reach not identical across years")
    method_count["2018_2022_I2M2_FR_same_station"]+=1
    orig18=float(a["numeric_EQR_2018"]);orig22=float(b["normalised_EQR"])
    # Main source month-matched CSV preserves ΔI2M2, compare to original site index.
    olddelta=float(r["actual_2018_to_2022_same_site_delta_I2M2_EQR"])
    if not math.isclose(orig22-orig18,olddelta,rel_tol=1e-9,abs_tol=1e-11):
        raise RuntimeError("Original 2018 2022 exact EQR field site numerator and source month study disagree")
    try:
        raw=float(a["raw_EQR_2018"])
        if not math.isclose(raw,orig18,rel_tol=1e-9,abs_tol=1e-10):
            original_raw2018_diff+=1
    except (TypeError,ValueError):
        original2018_raw_missing+=1

H6=np.array([r["HYBAS_ID_L06"] for r in rows])
fold=np.array([int(r["original_H6_source_fold"]) for r in rows])
if len(set(H6))!=24 or set(fold)!={0,1,2,3,4}:
    raise RuntimeError("Original 24 H6 parent geography/folds changed")
keys=[
 "log_original_H7_local_area_km2",
 "source_original_upstream_cropland_pct",
 "source_original_upstream_urban_pct",
 "log_original_H7_outlet_natural_DIS_m3s",
 "mean_original_site_modelled_RIVER_Q_log1p",
 "log1p_original_same_site_2018_22_pair_n",
 "original_COMMON_month_delta_log1p_TP_sample_effort"]
Y=np.array([float(r["observed_delta_EQR_same_site_H7"]) for r in rows])
Tm=np.array([float(r["observed_original_common_month_logTP_delta_H7"]) for r in rows])
Ta=np.array([float(r["observed_original_annual_logTP_delta_H7_same_sites"]) for r in rows])
X=np.array([[float(r[k]) for k in keys] for r in rows])
Xl=np.column_stack([X,[float(r["observed_original_common_month_equal_logTP_2018"]) for r in rows]])
def MSE(A,y,group):
    p=np.full(len(y),np.nan)
    for g in sorted(set(group)):
        test=group==g;train=~test
        if not np.any(test) or not np.any(train):
            raise RuntimeError("Original H6 group heldout subset empty")
        m=A[train].mean(axis=0);s=A[train].std(axis=0);s=np.where(s>1e-12,s,1)
        Z=(A[train]-m)/s
        Q=(A[test]-m)/s
        design=np.column_stack([np.ones(sum(train)),Z])
        regression=np.eye(design.shape[1])*10.
        regression[0,0]=0.
        coefficients=np.linalg.solve(design.T@design+regression,design.T@y[train])
        p[test]=np.column_stack([np.ones(sum(test)),Q])@coefficients
    if not np.all(np.isfinite(p)):raise RuntimeError("independent month-matched model not fully predicted")
    return float(np.mean((y-p)**2))
check={
 "primary_B0_H6_5fold":(MSE(X,Y,fold),report["new_monthmatch_primary_h6_geo5fold"]["original_baseline_MSE"]),
 "primary_monthmatchedTP_H6_5fold":(MSE(np.column_stack([X,Tm]),Y,fold),report["new_monthmatch_primary_h6_geo5fold"]["original_additional_MSE"]),
 "rawannualTP_same_sample_H6_5fold":(MSE(np.column_stack([X,Ta]),Y,fold),report["same_source_population_rawannual_TPdelta_and_same_month_effort_baseline_h6_geo5fold"]["original_additional_MSE"]),
 "monthmatchedTP_H6_leaveoneparentout":(MSE(np.column_stack([X,Tm]),Y,H6),report["posthoc_new_monthmatch_h6_leaveoneparentout"]["original_additional_MSE"]),
 "rawannualTP_H6_leaveoneparentout":(MSE(np.column_stack([X,Ta]),Y,H6),report["posthoc_same_source_rawannual_TPdelta_h6_leaveoneparentout"]["original_additional_MSE"]),
 "source_2018month_TP_baseline_H6_5fold":(MSE(Xl,Y,fold),report["sensitivity_B0_plus_original_2018_monthmatched_baselineTP_then_month_delta"]["original_baseline_MSE"]),
 "source_2018month_TP_baseline_plusΔTP_H6_5fold":(MSE(np.column_stack([Xl,Tm]),Y,fold),report["sensitivity_B0_plus_original_2018_monthmatched_baselineTP_then_month_delta"]["original_additional_MSE"])
}
for label,(actual,previous) in check.items():
    if not math.isclose(actual,float(previous),rel_tol=2e-10,abs_tol=1e-12):
        raise RuntimeError("Original main model monthly sampling numeric MSE not reproduced "+label+":"+str(actual)+"/"+str(previous))
sp=float(spearmanr(Tm,Y).statistic)
if not math.isclose(sp,report["original_H7_monthmatched_delta_logTP_vs_delta_I2M2_EQR_Spearman"],rel_tol=1e-12):
    raise RuntimeError("Original H7 month TP/EQR Spearman not reproduced")
spearsite=float(spearmanr(
 [float(r["actual_2018_to_2022_MONTHMATCHED_delta_logTP"]) for r in sites],
 [float(r["actual_2018_to_2022_same_site_delta_I2M2_EQR"]) for r in sites]).statistic)
if not math.isclose(spearsite,report["source_sitelevel_monthmatched_delta_TP_vs_I2M2_EQR_Spearman_descriptive"],rel_tol=1e-12):
    raise RuntimeError("Original 263 site monthly ranks changed")
# Cross-layer exact source regression replay: all 2 source terrestrial outcomes,
# 2 source MSE each, as separate genuine whole-H7 GBIF recording indices.
other={}
for label,source_out,extra in [
 ("GBIF_Plantae_source_speciesKey_S50_change","delta_GBIF_original_plant_S50",
   ["delta_GBIF_plant_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p"]),
 ("GBIF_Lycaenidae_fivefamily_index_logodds_change","delta_GBIF_original_Lycaenidae_logodds",
   ["delta_GBIF_fivefamily_record_opportunity_log1p","delta_GBIF_Animalia_record_opportunity_log1p"])
]:
    y=np.array([float(z[source_out]) for z in rows])
    Xother=np.array([[float(z[k]) for k in keys+extra] for z in rows])
    a=MSE(Xother,y,fold);b=MSE(np.column_stack([Xother,Tm]),y,fold)
    item=report["same_H7_2017_19_2020_22_GBIF_plant_butterfly_change"][label]["source_H6_group5fold_after_GBIF_opportunity_controls"]
    if not math.isclose(a,float(item["original_baseline_MSE"]),rel_tol=2e-10,abs_tol=1e-12) or \
       not math.isclose(b,float(item["original_additional_MSE"]),rel_tol=2e-10,abs_tol=1e-12):
        raise RuntimeError("GBIF original recording-proxy whole-H7 monthTP source model MSE numerical drift "+label)
    other[label]={"index_source_only":True,"baseline":a,"plus_monthTP":b}
output={
 "schema":"GESUI_FR263_REAL_TP_I2M2_SAME_CALENDAR_MONTH_INDEPENDENT_NUMPY_EEA_LINEAGE_QA_v1",
 "original_main_monthmodel_schema":report["schema"],
 "original_source_263_same_WQ_EEA_station_ids":len(sites),
 "original_source_hydrologic_73_H7":len(rows),
 "original_H6_parent_group_count":len(set(H6)),
 "original_EEA_biological_source_pipelines_both_use_resultNormalisedEQRValue":True,
 "2018_EEA_numeric_EQR_and_2022_EEA_normalised_EQR_both_from_publisher_normalized_field":True,
 "2018_EEA_raw_EQR_value_is_separate_unmixed_original_resultEQRValue":True,
 "2018_raw_and_normalized_EQR_original_source_sites_differ_n":original_raw2018_diff,
 "2018_original_raw_EQR_source_missing_or_non_numeric_n":original2018_raw_missing,
 "original_methods_and_country_exact_site_repeated":dict(method_count),
 "all_seven_original_primary_source_MSE_2018_monthlevel_H6_and_LOHO_reproduced":True,
 "independent_numeric_checks":{k:{"independent_NumPy_MSE":float(a),"original_sklearn_MSE":float(b)} for k,(a,b) in check.items()},
 "independent_original_same_site_73H7_monthTP_rho":sp,
 "independent_original_263_EEA_site_monthTP_rho":spearsite,
 "independent_source_GBIF_two_observation_index_models_reproduced":other,
 "not_independent_field_observation_or_crosscountry_validation":True,
 "limitations":[
    "EEA resultNormalisedEQRValue is the original published normalized EQR in both years, but identical I2M2 method/column does not prove identical field sampling timing or annual national calibration.",
    "Numeric agreement between separate Ridge solvers replays the SAME original field observations; it is not ecological replication.",
    "2018/2022 TP measurements on matching calendar months still differ in actual sample days, discharge regime, weather and biotic response-lag context.",
    "Whole-H7 GBIF plant/butterfly source indices do not observe real riparian host availability or pollination.",
    "Wastewater reporting source P n4/273 in this matched French TP/EQR sample cannot identify WWTP emission effects."
 ]}
path=ROOT/"FR263_MONTHMATCHED_TP_I2M2_EEA_SOURCE_NORMALIZATION_INDEPENDENT_RIDGE_QA.json"
path.write_text(json.dumps(output,sort_keys=True,indent=2,allow_nan=False))
print("FR263_MONTHTP_INDEPENDENT_EEA_LINEAGE_QA_START")
print(json.dumps({k:v for k,v in output.items() if k!="limitations"},sort_keys=True,allow_nan=False))
print("FR263_MONTHTP_INDEPENDENT_EEA_LINEAGE_QA_END")
