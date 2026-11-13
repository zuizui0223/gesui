#!/usr/bin/env python3
"""Independent NUMPY Ridge vs original sklearn genuine 2018–22 same-station
263 physical WQ TP / FR EEA I2M2 site pairs, 73 H7 /24 H6 group holdout.

NEW diagnostic (posthoc to original within-site Δ near-null):
both original years' SAME-site cross-sectional TP ↔ benthic EQR links,
and within-H7/within-time contrasts. No new site or year selection.
"""
import collections,csv,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2018_2022_TP_EQR_SPATIAL_TEMPORAL_INDEPENDENT"
ROOT.mkdir(parents=True,exist_ok=True)
summary=json.loads(pathlib.Path(os.environ["ORIGINAL_REAL_2018_22_TP_EQR_MAIN_REPORT"]).read_text())
if summary["schema"]!="GESUI_FR_TRUE_263_SAME_SITE_I2M2_MEASURED_TP_2018_22_DELTA_73H7_24H6_3BIOTA_v1":
    raise RuntimeError("Original numeric longitudinal study type changed")
with open(os.environ["ORIGINAL_REAL_2018_22_TP_EQR_263_SITES"],encoding="utf-8-sig",newline="") as f:
    sites=list(csv.DictReader(f))
with open(os.environ["ORIGINAL_REAL_2018_22_TP_EQR_73_H7"],encoding="utf-8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
if len(sites)!=263 or len(rows)!=73 or len({x["EEA_identical_native_site_key_2018_2022"] for x in sites})!=263:
    raise RuntimeError("Original 263 genuine repeatedly monitored TP/EQR sites or 73 H7 changed")
if len({r["HYBAS_ID_L07"] for r in rows})!=73 or len({r["HYBAS_ID_L06"] for r in rows})!=24:
    raise RuntimeError("Original 73 H7 /24 H6 source region inference grain changed")
H6=np.array([x["HYBAS_ID_L06"] for x in rows])
fold=np.array([int(x["original_H6_source_fold"]) for x in rows],dtype=int)
H6uniq,inv=np.unique(H6,return_inverse=True)
if len(H6uniq)!=24 or min(np.bincount(fold,minlength=5))<10:
    raise RuntimeError("Original 24 real H6 groups and source-only five folds changed")
def numeric(X,Y,groups):
    X=np.asarray(X,dtype=float);Y=np.asarray(Y,dtype=float)
    pred=np.full(len(Y),np.nan,dtype=float)
    for group in sorted(set(groups)):
        hold=groups==group;train=~hold
        if not np.any(hold) or not np.any(train):
            raise RuntimeError("Empty source group in independent CV")
        A=X[train];B=X[hold]
        mean=A.mean(axis=0);sd=A.std(axis=0)
        sd=np.where(sd>1e-12,sd,1)
        a=(A-mean)/sd;b=(B-mean)/sd
        a=np.column_stack([np.ones(len(a)),a])
        b=np.column_stack([np.ones(len(b)),b])
        pen=np.eye(a.shape[1])*10.;pen[0,0]=0.
        fit=np.linalg.solve(a.T@a+pen,a.T@Y[train])
        pred[hold]=b@fit
    if not np.all(np.isfinite(pred)):raise RuntimeError("Independent NumPy Ridge withheld prediction missing")
    return (Y-pred)**2
rng=np.random.default_rng(20261113)
sample_id=rng.integers(0,len(H6uniq),(2000,len(H6uniq)))
def gain(err0,err1):
    mse0=float(np.mean(err0));mse1=float(np.mean(err1))
    a=np.bincount(inv,weights=err0,minlength=len(H6uniq))
    b=np.bincount(inv,weights=err1,minlength=len(H6uniq))
    boot0=np.sum(a[sample_id],axis=1)
    boot1=np.sum(b[sample_id],axis=1)
    g=100*(boot0-boot1)/boot0
    return {"model_base_MSE":mse0,"model_plus_TP_MSE":mse1,
       "incremental_TP_predictive_gain_pct":100*(mse0-mse1)/mse0,
       "H6_parent_2000bootstrap_95_descriptive_gain_pct":np.percentile(g,[2.5,97.5]).tolist(),
       "original_bootstrap_gain_positive_fraction":float(np.mean(g>0))}
keys=summary["independent_geographic_real_H6_parent_out_of_region"]["source_original_static_landuse_and_riverflow_covariates"]
if len(keys)!=7:raise RuntimeError("Original source-locked environmental covariates expected seven")
X=np.array([[float(r[k]) for k in keys] for r in rows])
Y=np.array([float(r["mean_site_delta_original_EQR_2022_minus_2018"]) for r in rows])
T=np.array([float(r["mean_site_delta_log_original_measured_TP_2022_minus_2018"]) for r in rows])
old_e0=numeric(X,Y,fold)
old_e1=numeric(np.column_stack([X,T]),Y,fold)
old_LO0=numeric(X,Y,H6)
old_LO1=numeric(np.column_stack([X,T]),Y,H6)
target=summary["independent_geographic_real_H6_parent_out_of_region"]
original_checks={
 "H6_fivefold_base":(old_e0.mean(),target["original_fixed_H6_5fold_delta_TP_added_EQR_prediction"]["original_grouped_predictive_MSE_baseline"]),
 "H6_fivefold_plus_deltaTP":(old_e1.mean(),target["original_fixed_H6_5fold_delta_TP_added_EQR_prediction"]["original_grouped_predictive_MSE_with_delta_TP"]),
 "H6_24_parent_leaveout_base":(old_LO0.mean(),target["posthoc_leave_1_of_24_H6_parent_out_delta_TP_added_EQR_prediction"]["original_grouped_predictive_MSE_baseline"]),
 "H6_24_parent_leaveout_plus_deltaTP":(old_LO1.mean(),target["posthoc_leave_1_of_24_H6_parent_out_delta_TP_added_EQR_prediction"]["original_grouped_predictive_MSE_with_delta_TP"])
}
for name,(new,old) in original_checks.items():
    if not math.isclose(float(new),float(old),rel_tol=1e-10,abs_tol=1e-12):
        raise RuntimeError("Independent NumPy Ridge cannot reproduce source-original longitudinal model "+name)
rho=float(spearmanr(T,Y).statistic)
if not math.isclose(rho,summary["primary_same_73_H7_basin_delta_TP_vs_delta_EQR_Spearman"]["rho"],rel_tol=1e-12):
    raise RuntimeError("Original n73 actual ΔTP/EQR rank did not reproduce")
rho_site=float(spearmanr(
 [float(z["real_original_delta_log_TP_2022_minus_2018"]) for z in sites],
 [float(z["real_original_delta_I2M2_EQR_2022_minus_2018"]) for z in sites]).statistic)
if not math.isclose(rho_site,summary["secondary_site_raw_delta_Spearman_not263_independent_basins"]["rho"],rel_tol=1e-12):
    raise RuntimeError("Original n263 observed ΔTP/EQR site rank did not reproduce")
byH7=collections.defaultdict(list)
for x in sites:byH7[x["HYBAS_ID_L07"]].append(x)
slopes=[];top=0.;bottom=0.
for x in byH7.values():
    if len(x)<2:continue
    a=np.array([float(z["real_original_delta_log_TP_2022_minus_2018"]) for z in x])
    b=np.array([float(z["real_original_delta_I2M2_EQR_2022_minus_2018"]) for z in x])
    a-=a.mean();b-=b.mean()
    top+=float(a@b);bottom+=float(a@a)
sl=float(top/bottom)
if not math.isclose(sl,summary["within_H7_site_demeaned_delta_logTP_delta_EQR"]["slope_delta_EQR_per_delta_logTP_within_H7"],rel_tol=1e-12):
    raise RuntimeError("Original within-H7 2018→22 physical station demeaned TP–I2M2 slope did not independently reproduce")
# New strictly SAME original sample per year, not 2018 sites with no 2022
# followup or 2022 cohorts selected by a particular TP/EQR magnitude.
spatial={}
base_static=[k for k in keys if k!="mean_original_delta_log1p_TP_positive_sample_opportunities"]
if len(base_static)!=6:raise RuntimeError("Expected source original static env six predictors")
for year in (2018,2022):
    term=str(year)
    xTP=np.array([float(z[f"mean_original_{term}_logTP"]) for z in rows])
    y=np.array([float(z[f"mean_original_{term}_EQR"]) for z in rows])
    sample_eff=np.array([np.mean([math.log1p(float(q[f"original_WQ_valid_{term}_TP_observation_count"])) for q in byH7[z["HYBAS_ID_L07"]]]) for z in rows])
    XX=np.column_stack([
        np.array([[float(z[k]) for k in base_static] for z in rows]),
        sample_eff
    ])
    err0=numeric(XX,y,fold)
    err1=numeric(np.column_stack([XX,xTP]),y,fold)
    errL0=numeric(XX,y,H6)
    errL1=numeric(np.column_stack([XX,xTP]),y,H6)
    stationTP=[math.log(float(z[f"real_original_{term}_measured_TP_median_mg_L"])) for z in sites]
    stationE=[float(z[f"real_original_{term}_EEA_I2M2_{'numeric' if year==2018 else 'normalized'}_EQR"]) for z in sites]
    spatial[str(year)]={
        "year":year,
        "same_original_FR_I2M2_physical_station_pairs":263,
        "same_original_unique_H7":73,
        "same_original_H6_parents":24,
        "original_H7_source_TP_vs_EQR_Spearman_rho":float(spearmanr(xTP,y).statistic),
        "original_same_physical_station_TP_vs_EQR_Spearman_rho_descriptive_only":float(spearmanr(stationTP,stationE).statistic),
        "source_true_H6_fivefold_measured_TP_predictive_gain":gain(err0,err1),
        "posthoc_24_H6_parent_leaveout_measured_TP_predictive_gain":gain(errL0,errL1),
        "original_source_same_calendar_year_sampling_effort_adjustment":True
    }
TP18=np.array([float(r["mean_original_2018_logTP"]) for r in rows])
TP22=np.array([float(r["mean_original_2022_logTP"]) for r in rows])
EQ18=np.array([float(r["mean_original_2018_EQR"]) for r in rows])
EQ22=np.array([float(r["mean_original_2022_EQR"]) for r in rows])
between_rho=float(spearmanr((TP18+TP22)/2,(EQ18+EQ22)/2).statistic)
summaryout={
 "schema":"GESUI_FR_263_IDENTICAL_2018_22_TP_I2M2_SPATIAL_VS_TEMPORAL_INDEPENDENT_SOURCE_QA_v1",
 "original_source_main_study_schema":summary["schema"],
 "original_official_Peters_TP_archive_SHA256":summary["source_publisher_measured_TP_zip_SHA256"],
 "original_source_exact_263_station_id_73_H7_24_H6":True,
 "independent_NUMPY_reproduces_original_4_primary_H6_5fold_and_24leaveout_MSE":True,
 "independently_reproduced_4_original_source_metrics":{
    key:{"original":float(old),"independent_numPy":float(new)} for key,(new,old) in original_checks.items()},
 "independent_spearman_delta73_same_H7":rho,
 "independent_spearman_delta263_same_sites":rho_site,
 "independently_recomputed_original_within_H7_demeaned_deltaTP_delta_EQR_slope":sl,
 "posthoc_original_same_site_year_independent_spatial_TP_EQR_relationships":spatial,
 "posthoc_original_same_H7_between_geographic_meanTP_vs_meanEQ_Spearman":between_rho,
 "posthoc_same_H7_original_cross_year_stability_TP2018_vs_TP2022_Spearman":float(spearmanr(TP18,TP22).statistic),
 "posthoc_same_H7_original_cross_year_stability_I2M2_EQR2018_vs_EQR2022_Spearman":float(spearmanr(EQ18,EQ22).statistic),
 "delta2018_22_source_73_H7_rho_kept_as_temporal_result":rho,
 "not_a_new_independent_field_site_or_EEA_country":True,
 "underreported_WWTP_P_mass_still_positive_only_4_of_273_2022":True,
 "warnings":[
    "These diagnostics reuse all original 263 French I2M2 river site pairs and 73 H7/24 H6 physical geography from main actual TP time series.",
    "Year-separated spatial covariance does NOT imply causal TP harm or 2018/2022 temporal changes.",
    "2018 and 2022 EQR numeric source columns share I2M2 method but may differ in sampling dates or normalization reference.",
    "TP yearly medians are not simultaneous with EEA EQR benthic sample dates; lag or seasonal ecological responses may not align with simple 4-year change.",
    "Neither registered wastewater source P inventory nor whole H7 GBIF plant/butterfly indexed response measures real WWTP nutrient dose or pollination."
 ]}
root=ROOT/"FR2018_22_263_SAME_SITE_TP_EQR_SPATIAL_NEGATIVE_VS_TEMPORAL_NUL_NUMPY_INDEPENDENT.json"
root.write_text(json.dumps(summaryout,sort_keys=True,indent=2,allow_nan=False))
print("FR2018_22_TP_EQR_SPATIAL_VS_TEMPORAL_INDEPENDENT_START")
print(json.dumps({k:v for k,v in summaryout.items() if k!="warnings"},sort_keys=True,allow_nan=False))
print("FR2018_22_TP_EQR_SPATIAL_VS_TEMPORAL_INDEPENDENT_END")
