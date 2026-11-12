#!/usr/bin/env python3
"""Source-frozen cross-sectional EEA original 2018 and 2022 *reported* upstream
WWTP P/N source indicators vs EQR, original plant S50, butterfly5 index.

Original 454 H7 basin-periods, source annual N/P was audited and frozen
separately WITHOUT ecological outcomes. This is posthoc and NOT a sewage
causal effect; most station-basin source indices are unreported/zero.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_EU_454_REPORTED_NP_EFFECTS"
ROOT.mkdir(parents=True,exist_ok=True)
NP_SHA="012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"
NP_REPORT_SCHEMA="GESUI_454_EU_H7_ONE_WATERSHED_TWO_PERIOD_YEAR_SPECIFIC_REPORTED_NP_SOURCE_ONLY_v1"
BASE_REPORT_SCHEMA="GESUI_HYBAS7_234EU_source_matched_three_response_static_WWTP_gradient_posthoc_v1"
source_path=pathlib.Path(os.environ["SOURCE_EEA_2018_2022_NP_MANIFEST"])
if hashlib.sha256(source_path.read_bytes()).hexdigest()!=NP_SHA:
    raise RuntimeError("Original outcome-blind 454 H7 reported NP source SHA differs")
np_source=json.loads(source_path.read_text())["source_H7_period_reported_N_P_upstream_to_true_original_station"]
if len(np_source)!=454:raise RuntimeError("Original source-only NP sample N is not 454")
S={(q["HYBAS_ID_L07"],q["PERIOD"]):q for q in np_source}
if len(S)!=454:raise RuntimeError("Source-only annual NP origin H7 period duplicated")
if not all(r["source_valid"] and r["n_original_EEA_stations_source_qualified"]>=1 for r in np_source):
    raise RuntimeError("Source-only reported N/P/river flow station denominator missing")

original_report=json.loads(pathlib.Path(os.environ["SOURCE_PREVIOUS_THREE_OUTCOME_RESULT"]).read_text())
if original_report["schema"]!=BASE_REPORT_SCHEMA:
    raise RuntimeError("Original true 454 three layer source results do not match")
with open(os.environ["SOURCE_PREVIOUS_THREE_OUTCOME_COHORT"],encoding="utf-8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
if len(rows)!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in rows})!=454:
    raise RuntimeError("Previously verified same H7 cross-layer source period frame changed")
if set((r["HYBAS_ID_L07"],r["PERIOD"]) for r in rows)!=set(S):
    raise RuntimeError("Previous same H7 source period and actual reported loads not same native row units")
N=len(rows)
D=[]
for r in rows:
    s=S[(r["HYBAS_ID_L07"],r["PERIOD"])]
    present=1 if s["source_station_annual_upstream_NP_positive_n"]>0 else 0
    p=float(s["source_year_reported_P_mgL_equivalent_site_equal_mean"])
    n=float(s["source_year_reported_N_mgL_equivalent_site_equal_mean"])
    if not all(math.isfinite(x) and x>=0 for x in (p,n)):
        raise RuntimeError("Invalid source annual P/N converted through original river DIS")
    if (p>0 or n>0)!=bool(present):
        raise RuntimeError("Source positive detection indicator not equivalent to original report P/N mass")
    D.append({
      "HYBAS_ID_L07":r["HYBAS_ID_L07"],"PERIOD":r["PERIOD"],"country":r["country"],
      "source_positive_ind":present,
      "reported_P_mgL_equivalent_from_original_annual_tpy":p,
      "reported_N_mgL_equivalent_from_original_annual_tpy":n,
      "log1p_reported_P_ugL_equivalent":math.log1p(1000*p),
      "log1p_reported_N_mgL_equivalent":math.log1p(n),
      "n_source_station_positive":int(s["source_station_annual_upstream_NP_positive_n"])
    })
report_count=sum(x["source_positive_ind"] for x in D)
if report_count!=56 or len({x["HYBAS_ID_L07"] for x in D if x["source_positive_ind"]})!=38:
    raise RuntimeError("Original source-only positive-reporting 56 periods 38 basins changed")
if sum(x["PERIOD"]=="2017-2019" and x["source_positive_ind"] for x in D)!=28 or \
   sum(x["PERIOD"]=="2020-2022" and x["source_positive_ind"] for x in D)!=28:
    raise RuntimeError("Original source support in each year changed")
base_num=["log_local_H7_area_km2","log1p_original_EQR_site_n",
 "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
 "log1p_GBIF_5family_butterfly_source_records"]
if not original_report["source_grain"]["river_flow_covariate_complete_all_rows"]:
    raise RuntimeError("Exact original common baseline had valid river DIS_AV_CMS; do not change model")
base_num+=["log1p_original_EEA_site_DIS_AV_CMS"]
countries=sorted({r["country"] for r in rows})
methods=sorted({r["original_majority_BQE_method"] for r in rows})
x0=np.column_stack([
 np.ones(N),
 np.array([[float(r[k]) for k in base_num] for r in rows],dtype=float),
 np.array([[float(r["country"]==c) for c in countries[1:]] for r in rows],dtype=float),
 np.array([[float(r["original_majority_BQE_method"]==m) for m in methods[1:]] for r in rows],dtype=float),
 np.array([float(r["PERIOD"]=="2020-2022") for r in rows],dtype=float),
 np.array([float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in rows],dtype=float)
])
X={
 "B_original_static_HydroWASTE":x0,
 "B_plus_any_positive_original_reported_N_or_P":np.column_stack([x0,[r["source_positive_ind"] for r in D]]),
 "B_plus_annual_reported_P_flow_converted":np.column_stack([x0,[r["log1p_reported_P_ugL_equivalent"] for r in D]]),
 "B_plus_annual_reported_N_flow_converted":np.column_stack([x0,[r["log1p_reported_N_mgL_equivalent"] for r in D]]),
 "B_plus_both_reported_P_and_N":np.column_stack([x0,[r["log1p_reported_P_ugL_equivalent"] for r in D],[r["log1p_reported_N_mgL_equivalent"] for r in D]])
}
basins=[r["HYBAS_ID_L07"] for r in rows]
folds=np.array([int(hashlib.sha256(b.encode()).hexdigest()[:8],16)%5 for b in basins])
if min(np.bincount(folds,minlength=5))<35:
    raise RuntimeError("Original H7 basin hash fold support changed")
groups,inv=np.unique(basins,return_inverse=True)
G=len(groups)
if G!=232:raise RuntimeError("Expected exactly 232 source basin groups")
rng=np.random.default_rng(20261112)
draws=rng.integers(0,G,size=(2000,G))
def cvpred(A,Y):
    pred=np.empty(len(Y))
    for k in range(5):
        test=folds==k;train=~test
        Z=A.copy()
        mu=Z[train,1:].mean(axis=0)
        sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        a=Z[train];b=Z[test]
        reg=np.eye(a.shape[1])*10.;reg[0,0]=0.
        coef=np.linalg.solve(a.T@a+reg,a.T@Y[train])
        pred[test]=b@coef
    return pred
def h7clustercoef(A,Y,last_columns=1,weights=None):
    """Posthoc linear coefficient and H7 clustered descriptive interval."""
    Z=A.copy()
    mu=Z[:,1:].mean(axis=0)
    sd=Z[:,1:].std(axis=0)
    sd=np.where(sd>1e-12,sd,1.)
    Z[:,1:]=(Z[:,1:]-mu)/sd
    if weights is None:weights=np.ones(N)
    W=np.sqrt(weights)
    AA=Z*W[:,None]
    zz=Y*W
    beta=np.linalg.lstsq(AA,zz,rcond=None)[0]
    u=zz-AA@beta
    bread=np.linalg.pinv(AA.T@AA)
    meat=np.zeros((AA.shape[1],AA.shape[1]))
    for k in range(G):
        filt=inv==k
        v=AA[filt].T@u[filt]
        meat+=np.outer(v,v)
    rank=np.linalg.matrix_rank(AA)
    if N-rank<=0:raise RuntimeError("Too many nuisance predictors for original source sample")
    V=bread@meat@bread*G/(G-1)*(N-1)/(N-rank)
    se=np.sqrt(np.maximum(V.diagonal(),0))
    out=[]
    for idx in range(A.shape[1]-last_columns,A.shape[1]):
        out.append({"beta_per_1SD_reported_source_predictor":float(beta[idx]),
          "H7_cluster_descriptive95":[float(beta[idx]-1.96*se[idx]),float(beta[idx]+1.96*se[idx])],
          "original_source_predictor_SD":float(sd[idx-1])})
    return out
countrycounts=collections.Counter(r["country"] for r in rows)
countryweights=np.array([1/countrycounts[r["country"]] for r in rows],dtype=float)
countryweights/=countryweights.mean()
outcomes={
 "original_measured_EEA_benthic_EQR":"EQR_original_site_equal_basin_period_mean",
 "GBIF_original_plant_S50":"plant_S50_GBIF_speciesKey",
 "GBIF_five_family_Lycaenidae_log_odds":"lycaenidae_fivefamily_index_share_logodds"}
models={}
for label,column in outcomes.items():
    Y=np.array([float(r[column]) for r in rows])
    errors={model:(Y-cvpred(A,Y))**2 for model,A in X.items()}
    mse={k:float(e.mean()) for k,e in errors.items()}
    base=mse["B_original_static_HydroWASTE"]
    original_key={
     "original_measured_EEA_benthic_EQR":"actual_EEA_EQR_original_source",
     "GBIF_original_plant_S50":"GBIF_Plantae_speciesKey_S50_50rarefaction",
     "GBIF_five_family_Lycaenidae_log_odds":"GBIF_Lycaenidae_log_odds_5_butterfly_families"}[label]
    prev=original_report["results_three_response_separate_models"][original_key]["plus_static_modelled_WWTP_gradient_5fold_MSE"]
    if not math.isclose(base,prev,rel_tol=1e-10,abs_tol=1e-12):
        raise RuntimeError("Previous original 454 H7 baseline cannot be reproduced "+label)
    extra={}
    for k in list(X)[1:]:
        e0=np.bincount(inv,weights=errors["B_original_static_HydroWASTE"],minlength=G)
        e1=np.bincount(inv,weights=errors[k],minlength=G)
        sum0=np.sum(e0[draws],axis=1)
        sum1=np.sum(e1[draws],axis=1)
        gain_dist=100*(sum0-sum1)/sum0
        extra[k]={
          "H7_basin_disjoint_5fold_MSE_gain_pct_from_original_static_WWTP_baseline":100*(base-mse[k])/base,
          "H7_basin_bootstrap2000_descriptive95_gain_pct":np.percentile(gain_dist,[.025*100,.975*100]).tolist(),
          "H7_basin_bootstrap_positive_fraction":float(np.mean(gain_dist>0)),
          "H7_basin_cluster_OLS_descriptive_coefficients":h7clustercoef(X[k],Y,2 if k.endswith("both_reported_P_and_N") else 1),
          "country_balanced_H7_basin_cluster_OLS_coefficients":h7clustercoef(X[k],Y,2 if k.endswith("both_reported_P_and_N") else 1,countryweights)
        }
    raw=float(np.corrcoef(np.array([r["source_positive_ind"] for r in D]),Y)[0,1])
    models[label]={
        "original_same_basins_by_two_period_n":N,"H7_unique_basins":G,
        "original_baseline_static_WWTP_fivefold_MSE_exactly_reproduced":True,
        "mean_original_response":float(Y.mean()),
        "reported_positive_binary_raw_pearson_r_descriptive":raw,
        "original_same_basins_source_P_N_CV_MSE":mse,
        "original_added_reported_source_models":extra
    }
out={
 "schema":"GESUI_ORIGINAL_H7_EU454_REPORTED_2018_2022_UPSTREAM_NP_VS_3_GBIF_EEA_RESPONSE_v1",
 "original_study_inference":"original same 232 unique river basins, repeated 2018/2022 454 native period rows",
 "original_reported_NP_source_only_SHA256":NP_SHA,
 "previous_original_three_source_EEA_H7_result_schema":BASE_REPORT_SCHEMA,
 "source_support":{"basin_periods":N,"independent_H7_basins":G,
      "original_source_positive_periods":report_count,
      "original_source_positive_unique_H7_basins":38,
      "positive_source_EARLY":28,"positive_source_LATE":28,
      "original_source_positive_country_codes":sorted(set(z["country"] for z in D if z["source_positive_ind"])),
      "river_flow_converted_source_equivalent_not_receiver_actual_TNTP":True,
      "source_zero_is_no_positive_report_not_absence_of_sewage":True,
      "country_observations":dict(countrycounts),
      "source_P_and_N_nonzero_paired_log_correlation":float(np.corrcoef(
         np.log1p(np.array([d["reported_P_mgL_equivalent_from_original_annual_tpy"] for d in D])),
         np.log1p(np.array([d["reported_N_mgL_equivalent_from_original_annual_tpy"] for d in D])))[0,1])},
 "model_design":{"same_original_454_H7_rows_and_outcomes":True,"original_country_method_time_flow_observer_static_outlet_baseline":True,
       "5fold_hash_basin_disjoint_not_contiguous_geographic_independence":True,
       "original_source_P_and_N_tested_alone_and_joint_noncausal":True},
 "models":models,
 "limitations":[
   "Original EEA point-source mass reporting positive on only 56/454 H7 period rows, 38 unique H7 basins; other source matched zero is not true no sewage.",
   "HydroRIVERS nearest-upstream reporting linkage is a conservative inventory and 50km window; missing original plant reporting/proximity confounds.",
   "DIS_AV_CMS is climatological mean flow, annual loads are reported tons/year; synthetic mg/L equivalent is full-delivery upper-input concentration, NOT sampled river chemistry.",
   "P and N source indicators are strongly collinear, and their separate model predictive gains cannot be converted to independent nutrient effects.",
   "Original country and method and observed GBIF recording effort controls do not remove agriculture, urbanization, environmental gradients or survey collider bias.",
   "Plant S50 original GBIF speciesKey rarefied records and butterfly five-family record share are not botanical or pollination services.",
   "All analyses are posthoc after static HydroWASTE signal known, and observed periods are repeated basin IDs. No causal ecosystem service chain established."
 ]}
outpath=ROOT/"EEA454_SAME_H7_PERIOD_YEAR_REPORTED_NP_SITE_FLOW_EQUIVALENT_3LAYER_POSTHOC.json"
outpath.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False))
print("H7_454_ACTUAL_REPORTED_UPSTREAM_NP_MACRO_RESULT_START",flush=True)
brief={k:v for k,v in out.items() if k not in ("limitations",)}
print(json.dumps(brief,sort_keys=True,allow_nan=False),flush=True)
print("H7_454_ACTUAL_REPORTED_UPSTREAM_NP_MACRO_RESULT_END",flush=True)
