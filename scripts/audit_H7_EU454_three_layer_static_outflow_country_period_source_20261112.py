#!/usr/bin/env python3
"""Independent source-saved native H7 454-row posthoc diagnostics:
rebuild original Ridge CV without importing primary analysis, compare signs
in P1 vs P2, FR vs nonFR, and exclude cross-border/mixed-method basins.

DOES NOT treat French catchments as 299 independent country replications,
and does NOT use real 2018/2022 wastewater emissions, emergence or pollination.
"""
import csv,hashlib,json,math,os,pathlib,collections
import numpy as np
from scipy.stats import spearmanr

src=pathlib.Path(os.environ["H7_ORIGINAL_THREE_LAYER_SOURCE_REPORT"])
tab=pathlib.Path(os.environ["H7_ORIGINAL_THREE_LAYER_ROWS"])
data=json.loads(src.read_text())
expected="GESUI_HYBAS7_234EU_source_matched_three_response_static_WWTP_gradient_posthoc_v1"
if data.get("schema")!=expected:raise RuntimeError("Original saved source analysis type changed")
with tab.open(encoding="utf8",newline="") as f:rows=list(csv.DictReader(f))
if len(rows)!=454 or len({(x["HYBAS_ID_L07"],x["PERIOD"]) for x in rows})!=454:
    raise RuntimeError("Genuine source 454 H7-period rows changed")
if len({x["HYBAS_ID_L07"] for x in rows})!=232:
    raise RuntimeError("Original source 232 unique H7 basins changed")
if not all(x["static_flow_is_not_actual_annual_discharge_or_receptor_specific_dose"]=="True" for x in rows):
    raise RuntimeError("Original static source interpretation changed")
by=collections.defaultdict(list)
for z in rows:by[z["HYBAS_ID_L07"]].append(z)
for h,sub in by.items():
    if len(sub)==2 and abs(float(sub[0]["W_static_upstream_MODEL_flow_per_UpArea"])-float(sub[1]["W_static_upstream_MODEL_flow_per_UpArea"]))>1e-11:
        raise RuntimeError("HydroWASTE static annual flow artificially varies inside one H7 basin")
if sorted(data["source_grain"]["countries"])!=sorted({r["country"] for r in rows}):
    raise RuntimeError("Original source country strata changed")
base_cols=["log_local_H7_area_km2","log1p_original_EQR_site_n",
  "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
  "log1p_GBIF_5family_butterfly_source_records"]
if data["source_grain"]["river_flow_covariate_complete_all_rows"]:
    base_cols+=["log1p_original_EEA_site_DIS_AV_CMS"]
def build(sub,pressure=False):
    countries=sorted({r["country"] for r in sub})
    methods=sorted({r["original_majority_BQE_method"] for r in sub})
    v=[]
    for r in sub:
        row=[1.]+[float(r[k]) for k in base_cols]+\
            [float(r["country"]==c) for c in countries[1:]]+\
            [float(r["original_majority_BQE_method"]==m) for m in methods[1:]]+\
            [float(r["PERIOD"]=="2020-2022")]
        if pressure:row.append(float(r["W_static_upstream_MODEL_flow_per_UpArea"]))
        v.append(row)
    return np.array(v,dtype=float)
def fold(h):
    return int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5
def oos(X,Y,groups):
    folds=np.array([fold(r) for r in groups])
    if min(np.bincount(folds,minlength=5))<2:
        return None
    pp=np.empty(len(Y))
    for f in range(5):
        train=folds!=f;test=~train
        Z=X.copy()
        mu=Z[train,1:].mean(axis=0)
        sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        A=Z[train]
        reg=np.eye(A.shape[1])*10;reg[0,0]=0
        coef=np.linalg.solve(A.T@A+reg,A.T@Y[train])
        pp[test]=Z[test]@coef
    return float(np.mean((Y-pp)**2))
outcome_keys={
 "EQR":"EQR_original_site_equal_basin_period_mean",
 "PlantS50":"plant_S50_GBIF_speciesKey",
 "Lycaenidae5":"lycaenidae_fivefamily_index_share_logodds"
}
def cohort_metric(sub):
    if len(sub)<25 or len({r["HYBAS_ID_L07"] for r in sub})<15:
        return {"n_basin_period":len(sub),"status":"insufficient_source"}
    b=build(sub);w=build(sub,True)
    ids=[r["HYBAS_ID_L07"] for r in sub]
    z={"n_basin_period":len(sub),"n_basins":len(set(ids)),
       "n_countries":len({r["country"] for r in sub}),
       "n_FR":sum(r["country"]=="FR" for r in sub)}
    for k,key in outcome_keys.items():
        Y=np.array([float(r[key]) for r in sub])
        a=oos(b,Y,ids);c=oos(w,Y,ids)
        corr=spearmanr(Y,[float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in sub])
        z[k]={"OOS_baseline_MSE":a,"OOS_with_W_MSE":c,
           "OOS_gain_pct":100*(a-c)/a if a is not None and c is not None and a>0 else None,
           "raw_rho":float(corr.statistic)}
    return z
subsets={
 "original_full_454":rows,
 "2017_19_only":[r for r in rows if r["PERIOD"]=="2017-2019"],
 "2020_22_only":[r for r in rows if r["PERIOD"]=="2020-2022"],
 "France_only":[r for r in rows if r["country"]=="FR"],
 "exclude_France":[r for r in rows if r["country"]!="FR"],
 "exclude_mixed_BQE_and_border_country":[r for r in rows if int(r["original_method_count_in_basin_year"])==1 and int(r["original_country_count_in_basin_year"])==1],
}
metrics={label:cohort_metric(sub) for label,sub in subsets.items()}
# Verify the original saved exact full-sample out-of-basin MSE values; any
# reimplementation discrepancy must STOP rather than displaying results.
for label,original in [
    ("EQR","actual_EEA_EQR_original_source"),
    ("PlantS50","GBIF_Plantae_speciesKey_S50_50rarefaction"),
    ("Lycaenidae5","GBIF_Lycaenidae_log_odds_5_butterfly_families")
]:
    m=metrics["original_full_454"][label]
    given=data["results_three_response_separate_models"][original]
    for field,other in [
        ("OOS_baseline_MSE","baseline_5fold_hash_basin_disjoint_MSE"),
        ("OOS_with_W_MSE","plus_static_modelled_WWTP_gradient_5fold_MSE")]:
        if not math.isclose(m[field],given[other],rel_tol=1e-10,abs_tol=1e-12):
            raise RuntimeError("Independent original source CV does not reproduce "+label+" "+field)
out={
 "schema":"gesui_same_H7_454_EU_three_layer_outlet_gradient_independent_posthoc_subset_QA_v1",
 "original_study_source_report_schema":expected,
 "original_study_basin_count":232,
 "original_source_basinperiod_n":454,
 "original_four_source_SHA256":data["original_source_candidates_SHA256"],
 "main_5fold_OOS_exact_independently_reproduced":True,
 "source_repeat_basins_both_periods":sum(len(v)==2 for v in by.values()),
 "all_outcomes_are_separate_source_proxies_not_realized_services":True,
 "subgroup_sensitivity_POSTHOC":metrics,
 "constraints":[
   "This is a country and period sensitivity after the full model sign was observed, not an independent holdout study.",
   "France and nonFrance share distinct original EEA biological methods; source method/country dummies remain imperfect proxies.",
   "H7 hash CV is basin-disjoint, not a country-out or contiguous-ecoregion spatial CV.",
   "No source N/P delivered to original river stations, and HydroWASTE static modeled pressure does not vary by the two periods.",
   "Measured plant flowering, butterfly activity, predator diet and pollination are absent; even a negative EQR result cannot establish a mediated ecosystem service."
 ]}
folder=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_454_OUTLET_GRADIENT_AUDIT"
folder.mkdir(parents=True,exist_ok=True)
(folder/"H7_454_EU_3LAYER_SOURCE_SAME_BASIN_WWTP_INDEPENDENT_PERIOD_COUNTRY_QA.json").write_text(json.dumps(out,indent=2,sort_keys=True))
print("H7_454_LAYER_INDEPENDENT_COUNTRY_PERIOD_QA_START")
print(json.dumps({k:v for k,v in out.items() if k!="constraints"},sort_keys=True))
print("H7_454_LAYER_INDEPENDENT_COUNTRY_PERIOD_QA_END")
