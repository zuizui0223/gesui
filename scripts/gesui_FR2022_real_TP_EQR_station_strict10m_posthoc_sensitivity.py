#!/usr/bin/env python3
"""Posthoc strict <=10m exact station coordinate sensitivity of original
273 real WQ TP–EEA I2M2 EQR field sites; main 250m n273 unchanged.

Reaggregate only from original previously saved 273 physical station matched
source and original 74 native H7 sources, no additional TP/EQR outcome fetch.
Predeclared >=200 strict source pairs, >=60 H7 and >=15 H6 for QA sensitivity.
"""
import csv,json,math,os,pathlib,collections
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2022_TP_EQR_strict10m_QA"
ROOT.mkdir(parents=True,exist_ok=True)
main=json.loads(pathlib.Path(os.environ["ORIGINAL_FR2022_TP_EQR_RESULTS"]).read_text())
with open(os.environ["ORIGINAL_FR2022_TP_EQR_273_SITES"],encoding="utf-8-sig",newline="") as f:sites=list(csv.DictReader(f))
with open(os.environ["ORIGINAL_FR2022_TP_EQR_74_H7"],encoding="utf-8-sig",newline="") as f:h7=list(csv.DictReader(f))
if len(sites)!=273 or len(h7)!=74 or main["schema"]!="GESUI_FR_2022_273_SAME_RIVER_SITE_TP_I2M2_EQR_74_H7_24_H6_3BIOTA_RESULTS_v1":
    raise RuntimeError("Original true 2022 colocated water TP/EQR source sample changed")
h7src={r["HYBAS_ID_L07"]:r for r in h7}
strict=[s for s in sites if float(s["original_source_station_coordinate_distance_m"])<=10.]
by=collections.defaultdict(list)
for item in strict:by[item["HYBAS_ID_L07"]].append(item)
parents={h7src[h]["HYBAS_ID_L06"] for h in by}
if len(strict)<200 or len(by)<60 or len(parents)<15:
    raise RuntimeError("Original strictly-coincident station 10m source robustness gate FAILED, do not fit")
r10=[]
for h,ss in sorted(by.items()):
    old=h7src[h]
    mean=lambda col:float(np.mean([float(r[col]) for r in ss]))
    row=dict(old)
    row["original_matched_same_reach_TP_EQR_station_n"]=len(ss)
    row["mean_log_real_2022_river_TP_mgL_at_exact_EQR_sites"]=mean("log_real_measured_river_TP_2022_mgL")
    row["mean_same_site_original_2022_EEA_I2M2_EQR"]=mean("original_EEA_benthic_I2M2_2022_normalized_EQR")
    row["log_mean_original_site_DIS_AV_CMS"]=math.log1p(mean("original_EEA_site_modelled_river_DIS_AV_CMS"))
    row["log1p_n_actual_matched_TP_EQR_river_sites"]=math.log1p(len(ss))
    r10.append(row)
n=len(r10)
tp=np.array([float(x["mean_log_real_2022_river_TP_mgL_at_exact_EQR_sites"]) for x in r10])
eqr=np.array([float(x["mean_same_site_original_2022_EEA_I2M2_EQR"]) for x in r10])
rho=spearmanr(tp,eqr)
models=main["same_native_H7_2020_2022_response_results"]["true_same_2022_I2M2_benthic_EQR"]
cols=models["source_baseline_covariates"]
X0=np.array([[float(r[col]) for col in cols] for r in r10])
X1=np.column_stack([X0,tp])
fold=np.array([int(r["original_source_fold"]) for r in r10])
parents_arr=np.array([r["HYBAS_ID_L06"] for r in r10])
def fit(X,y):
    err=np.empty(len(y))
    for k in range(5):
        train=fold!=k;test=~train
        if sum(train)<30 or sum(test)<4:raise RuntimeError("Original geo-fold heldout too small at <=10m")
        A=X[train];B=X[test]
        mu=A.mean(axis=0);sd=A.std(axis=0)
        sd=np.where(sd>1e-12,sd,1)
        re=Ridge(alpha=10,fit_intercept=True,solver="svd")
        re.fit((A-mu)/sd,y[train])
        err[test]=(y[test]-re.predict((B-mu)/sd))**2
    return err
e0=fit(X0,eqr);e1=fit(X1,eqr)
G,inv=np.unique(parents_arr,return_inverse=True)
rng=np.random.default_rng(20261113)
idx=rng.integers(0,len(G),(2000,len(G)))
s0=np.bincount(inv,weights=e0,minlength=len(G))[idx].sum(axis=1)
s1=np.bincount(inv,weights=e1,minlength=len(G))[idx].sum(axis=1)
boot=100*(s0-s1)/s0
within=[]
for h,ss in by.items():
    if len(ss)<2:continue
    x=np.array([float(r["log_real_measured_river_TP_2022_mgL"]) for r in ss])
    y=np.array([float(r["original_EEA_benthic_I2M2_2022_normalized_EQR"]) for r in ss])
    x-=x.mean();y-=y.mean()
    within.append((float(x@y),float(x@x)))
slope=sum(x[0] for x in within)/sum(x[1] for x in within)
out={
 "schema":"GESUI_FR2022_TRUE_TP_EQR_10M_COINCIDENT_2022_EEA_WQ_POSTHOC_SPATIAL_QA_v1",
 "original_n273_h7_n74_pre_source_match_sha256":main["matched_original_250m_TP_EQR_station_source_SHA256"],
 "original_source_only_coincident_distance_threshold_m":10.,
 "strict_original_TP_EQR_same_reach_site_pairs":len(strict),
 "strict_original_H7_catchments":n,
 "strict_original_H6_parent_groups":len(G),
 "source_original_country_I2M2_France_only":True,
 "observed_same_H7_rho_TP_vs_EQR":float(rho.statistic),
 "observed_same_H7_descriptive_p_not_independent_h6":float(rho.pvalue),
 "within_H7_site_demeaned_TP_vs_EQR_slope":slope,
 "within_H7_with_multiple_sites":len(within),
 "H6_fold_holdout_MSE_baseline":float(e0.mean()),
 "H6_fold_holdout_MSE_plus_actual_river_TP":float(e1.mean()),
 "H6_fold_holdout_MSE_incremental_gain_pct":100*float(e0.mean()-e1.mean())/float(e0.mean()),
 "H6_parent_2000_bootstrap_descriptive95_gain_pct":np.percentile(boot,[2.5,97.5]).tolist(),
 "original_250m_baseline_source_for_comparison":{
   "matched_stations":len(sites),"basins":len(h7),
   "H7_spearman":main["same_native_H7_2020_2022_response_results"]["true_same_2022_I2M2_benthic_EQR"]["raw_basin_level_spearman_rho_log_measured_TP_vs_bio"],
   "H6_fold_MSE_gain_pct":main["same_native_H7_2020_2022_response_results"]["true_same_2022_I2M2_benthic_EQR"]["after_landcover_and_river_flow_MEASURED_TP_incremental_H6_CV"]["H6_parent_original_GroupKFold5_MSE_gain_pct"]
 },
 "posthoc_10m_source_subset_not_independent_replication":True,
 "no_effluent_NP_wwtp_pollination_or_causal_attribution":True}
(ROOT/"FR2022_10M_PHYSICAL_TP_I2M2_EQR_SOURCE_GEO_SENSITIVITY.json").write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False))
print("FR2022_STRICT_10M_COLOCATION_SENSITIVITY_START",flush=True)
print(json.dumps(out,sort_keys=True,allow_nan=False),flush=True)
print("FR2022_STRICT_10M_COLOCATION_SENSITIVITY_END",flush=True)
