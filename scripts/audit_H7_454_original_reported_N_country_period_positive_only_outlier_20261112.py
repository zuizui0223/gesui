#!/usr/bin/env python3
"""Independent (from saved source original rows) 454-basin-period upstream
reported WWTP-N source result re-calculation. Posthoc country/year/winsor QA.

Source nonreported index != clean water. No fresh source retrieval and
NO independent field experiment. All source original run IDs are immutable.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_454_reported_N_independent"
ROOT.mkdir(parents=True,exist_ok=True)
NP_SHA="012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"
src=pathlib.Path(os.environ["ORIGINAL_REPORTED_NP_SOURCE_MANIFEST"])
if hashlib.sha256(src.read_bytes()).hexdigest()!=NP_SHA:raise RuntimeError("Original source-only reported P/N immutable sample SHA mismatch")
v=json.loads(src.read_text())
byNP={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in v["source_H7_period_reported_N_P_upstream_to_true_original_station"]}
if len(byNP)!=454:raise RuntimeError("Original 454 source units changed")
full=json.loads(pathlib.Path(os.environ["ORIGINAL_454_TRI_LAYER_REPORTED_NP_RESULT"]).read_text())
if full["schema"]!="GESUI_ORIGINAL_H7_EU454_REPORTED_2018_2022_UPSTREAM_NP_VS_3_GBIF_EEA_RESPONSE_v1":
    raise RuntimeError("First original 454 N/P biological result schema changed")
with open(os.environ["ORIGINAL_454_TRI_LAYER_SOURCE_ROWS"],encoding="utf-8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
if len(rows)!=454:raise RuntimeError("First original 454 source rows changed")
for r in rows:
    if (r["HYBAS_ID_L07"],r["PERIOD"]) not in byNP:
        raise RuntimeError("Nonmatching native H7 original station/nutrient period")
    source=byNP[(r["HYBAS_ID_L07"],r["PERIOD"])]
    if not source["source_valid"]:raise RuntimeError("Original source-zero and missing data confused")
    r["original_report_N"]=math.log1p(float(source["source_year_reported_N_mgL_equivalent_site_equal_mean"]))
    r["original_report_P"]=math.log1p(1000*float(source["source_year_reported_P_mgL_equivalent_site_equal_mean"]))
    r["source_report_positive"]=int(source["source_station_annual_upstream_NP_positive_n"]>0)
assert sum(r["source_report_positive"] for r in rows)==56
def stage(sub,censor=None):
    n=len(sub)
    basin={r["HYBAS_ID_L07"] for r in sub}
    if n<100 or len(basin)<50:return {"status":"insufficient","n_periods":n,"n_basins":len(basin)}
    cats=sorted({r["country"] for r in sub})
    methods=sorted({r["original_majority_BQE_method"] for r in sub})
    keys=["log_local_H7_area_km2","log1p_original_EQR_site_n",
     "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
     "log1p_GBIF_5family_butterfly_source_records","log1p_original_EEA_site_DIS_AV_CMS"]
    X=np.column_stack([
      np.ones(n),
      np.array([[float(r[k]) for k in keys] for r in sub]),
      np.array([[float(r["country"]==c) for c in cats[1:]] for r in sub]),
      np.array([[float(r["original_majority_BQE_method"]==m) for m in methods[1:]] for r in sub]),
      np.array([float(r["PERIOD"]=="2020-2022") for r in sub]),
      np.array([float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in sub])
    ])
    raw=np.array([r["original_report_N"] for r in sub],dtype=float)
    if censor is not None:raw=np.minimum(raw,censor)
    XN=np.column_stack([X,raw])
    fold=np.array([int(hashlib.sha256(r["HYBAS_ID_L07"].encode()).hexdigest()[:8],16)%5 for r in sub])
    if min(np.bincount(fold,minlength=5))<4:return {"status":"fold_low","n_periods":n}
    out={"n_periods":n,"n_H7_basins":len(basin),"positive_source_periods":sum(r["source_report_positive"] for r in sub),
         "country_count":len(cats),"FR_source_periods":sum(r["country"]=="FR" for r in sub),"model":{}}
    for name,ykey,fullkey in [
       ("EQR","EQR_original_site_equal_basin_period_mean","original_measured_EEA_benthic_EQR"),
       ("plantS50","plant_S50_GBIF_speciesKey","GBIF_original_plant_S50"),
       ("Lycaenidae5","lycaenidae_fivefamily_index_share_logodds","GBIF_five_family_Lycaenidae_log_odds")]:
        Y=np.array([float(r[ykey]) for r in sub],dtype=float)
        def fit(A):
            pp=np.empty(n)
            for f in range(5):
                test=fold==f;train=~test
                Z=A.copy()
                mu=Z[train,1:].mean(axis=0)
                sd=Z[train,1:].std(axis=0)
                sd=np.where(sd>1e-12,sd,1.)
                Z[:,1:]=(Z[:,1:]-mu)/sd
                a=Z[train];b=Z[test]
                reg=np.eye(a.shape[1])*10.;reg[0,0]=0.
                beta=np.linalg.solve(a.T@a+reg,a.T@Y[train])
                pp[test]=b@beta
            return float(np.mean((Y-pp)**2))
        b=fit(X);w=fit(XN)
        out["model"][name]={
           "baseline_ridge10_5fold_MSE":b,
           "additional_reported_N_MSE":w,
           "additional_reported_N_gain_pct":100*(b-w)/b,
           "raw_rho_reported_N_with_response":float(spearmanr(raw,Y).statistic)
        }
    return out
original=stage(rows)
for key,modelkey in [
    ("EQR","original_measured_EEA_benthic_EQR"),
    ("plantS50","GBIF_original_plant_S50"),
    ("Lycaenidae5","GBIF_five_family_Lycaenidae_log_odds")
]:
    old=full["models"][modelkey]
    source_MSE=old["original_same_basins_source_P_N_CV_MSE"]
    if not math.isclose(original["model"][key]["baseline_ridge10_5fold_MSE"],source_MSE["B_original_static_HydroWASTE"],rel_tol=1e-12,abs_tol=1e-12):
        raise RuntimeError("Manual independent model original H7-ridge baseline disagree "+key)
    if not math.isclose(original["model"][key]["additional_reported_N_MSE"],source_MSE["B_plus_annual_reported_N_flow_converted"],rel_tol=1e-12,abs_tol=1e-12):
        raise RuntimeError("Manual independent model original H7-ridge reported N mismatch "+key)
bygroup={
 "full_original_source_454":original,
 "year_2017_19":[r for r in rows if r["PERIOD"]=="2017-2019"],
 "year_2020_22":[r for r in rows if r["PERIOD"]=="2020-2022"],
 "FR_only":[r for r in rows if r["country"]=="FR"],
 "without_FR":[r for r in rows if r["country"]!="FR"],
 "same_single_country_and_BQE_method":[r for r in rows if
   int(r["original_country_count_in_basin_year"])==1 and int(r["original_method_count_in_basin_year"])==1]
}
# All listed are outcome-blind source-key subsets, but looked at after main result.
for key in list(bygroup):
    if key!="full_original_source_454":
        bygroup[key]=stage(bygroup[key])
reportedPositive=np.array([r["original_report_N"] for r in rows if r["original_report_N"]>0],dtype=float)
if len(reportedPositive)<25:raise RuntimeError("Original positive reported N source too sparse for predeclared capping")
cap90=float(np.quantile(reportedPositive,.90))
cap95=float(np.quantile(reportedPositive,.95))
bygroup["source_positive_N_intensity_cap_95pct_posthoc"]=stage(rows,censor=cap95)
bygroup["source_positive_N_intensity_cap_90pct_posthoc"]=stage(rows,censor=cap90)
positive=[r for r in rows if r["source_report_positive"]>0]
conditional={}
for response,key in [
 ("EQR","EQR_original_site_equal_basin_period_mean"),
 ("plant_S50","plant_S50_GBIF_speciesKey"),
 ("Lycae5","lycaenidae_fivefamily_index_share_logodds")
]:
    for exposure in ["original_report_N","original_report_P"]:
        v=np.array([r[exposure] for r in positive])
        yy=np.array([float(r[key]) for r in positive])
        conditional[response+"_"+exposure]={
         "n_original_positive_source_basin_periods":len(positive),
         "spearman_r":float(spearmanr(v,yy).statistic)}
nonzero_country=collections.Counter(r["country"] for r in positive)
out={
 "schema":"GESUI_ORIGINAL_H7_EU454_REPORTED_N_FLOW_SOURCE_POSTHOC_ROBUSTNESS_v1",
 "source_only_original_NP_manifest_SHA256":NP_SHA,
 "source_original_sparsity":{"positive_basin_periods":len(positive),
   "positive_unique_H7_basins":len({r["HYBAS_ID_L07"] for r in positive}),
   "positive_country_counts":dict(nonzero_country),
   "N_positive_unique_records":len(reportedPositive)},
 "independent_full_3response_5fold_ridge10_reproduction":True,
 "posthoc_cap_original_nonzero_quantiles":{"P90_reported_N_log1p_mgL":cap90,"P95_reported_N_log1p_mgL":cap95},
 "original_same_H7_basin_period_subset_sensitivities":bygroup,
 "same_source_reported_positive_ONLY_raw_rank_diagnostics_NOT_causal":conditional,
 "limitations":[
  "All subgroup/cap checks selected after the initial source N sign was observed, so only descriptive sensitivity evidence.",
  "A very small minority of basin-periods have matched positive source records; source-zero is not a verified no-WWTP receptor/control group.",
  "Source P/N event distribution may vary by country and year, with unknown 2022 WWTP point-site reporting completeness.",
  "Site flow denominator is modeled long-run DIS_AV_CMS; P/N sources are originally annual reported tons/year, not observed N/P at the river monitoring stations.",
  "This is the same 454 original EEA/GBIF source sample, not a new field-site replication."
 ]}
path=ROOT/"H7_454_REPORTED_N_TRUE_SOURCE_COUNTRY_PERIOD_POSITIVE_ONLY_INDEPENDENT_QA.json"
path.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False))
print("H7_454_N_SOURCE_INDEPENDENT_ROBUSTNESS_START",flush=True)
print(json.dumps({k:v for k,v in out.items() if k!="limitations"},sort_keys=True,allow_nan=False),flush=True)
print("H7_454_N_SOURCE_INDEPENDENT_ROBUSTNESS_END",flush=True)
