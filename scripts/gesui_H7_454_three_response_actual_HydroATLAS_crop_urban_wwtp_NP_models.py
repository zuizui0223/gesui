#!/usr/bin/env python3
"""Posthoc original-data macro SOURCE: same 232 native H7 EU basins ×454
period rows, real HydroATLAS v1 native ID-match static cropland and urban,
BasinATLAS natural discharge. Keep EEA EQR, GBIF botanical S50 and Lycaenidae
response observation models distinct, compare static modeled WWTP outlet flow
and real site-reported annual upstream N/P *indicators* after controls.

HydroATLAS true publisher original L07 DBF CRC and original 1219-ID SHA were
independently source verified BEFORE opening the new landuse-adjusted effects.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr

OUT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_454_LANDCOVER_REPORT"
OUT.mkdir(parents=True,exist_ok=True)
ATLAS_SHA="ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686"
NP_SHA="012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"
OLD_SPATIAL="GESUI_HYBAS7_234EU_source_matched_three_response_static_WWTP_gradient_posthoc_v1"
OLD_ANNUAL="GESUI_ORIGINAL_H7_EU454_REPORTED_2018_2022_UPSTREAM_NP_VS_3_GBIF_EEA_RESPONSE_v1"

def read_json_source(name,expectedSHA=None):
    p=pathlib.Path(os.environ[name]);raw=p.read_bytes()
    if expectedSHA is not None and hashlib.sha256(raw).hexdigest()!=expectedSHA:
        raise RuntimeError("Previous independently verified native source SHA changed "+name)
    return json.loads(raw)
atlas=read_json_source("SOURCE_ORIGINAL_H7_HYDROATLAS_1219",ATLAS_SHA)
np_source=read_json_source("SOURCE_ORIGINAL_H7_REPORT_NP_454",NP_SHA)
prior=read_json_source("SOURCE_ORIGINAL_H7_THREE_LAYER_SPATIAL_REPORT")
prior_annual=read_json_source("SOURCE_ORIGINAL_H7_THREE_LAYER_NP_REPORT")
if prior["schema"]!=OLD_SPATIAL or prior_annual["schema"]!=OLD_ANNUAL:
    raise RuntimeError("Previously verified original source record expected schema mismatch")
with open(os.environ["SOURCE_ORIGINAL_H7_THREE_LAYER_ROWS"],encoding="utf8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
if len(rows)!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in rows})!=454:
    raise RuntimeError("Original 454 native basin-period rows changed")
h7={k:v for k,v in atlas["original_1219_H7_covariate_by_real_ID"].items()}
if len(h7)!=1219:raise RuntimeError("Original authentic native BasinATLAS study level7 coverage not exactly 1219")
np_by={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in np_source["source_H7_period_reported_N_P_upstream_to_true_original_station"]}
if len(np_by)!=454:raise RuntimeError("Original 2018/22 upstream source station H7 period join changed")
landuse=[]
np_parts=[]
rawflows=[]
for x in rows:
    key=x["HYBAS_ID_L07"]
    if key not in h7:raise RuntimeError("Native source H7 ID missing from original genuine BasinATLAS DBF")
    l=h7[key]
    req=["crp_pc_use","urb_pc_use","crp_pc_sse","urb_pc_sse","dis_m3_pyr"]
    values={}
    for k in req:
        if l[k] is None:raise RuntimeError("Original source native attribute field missing "+k)
        v=float(l[k])
        if not math.isfinite(v):raise RuntimeError("Original HydroATLAS native field nonfinite "+k)
        values[k]=v
    for k in ("crp_pc_use","urb_pc_use","crp_pc_sse","urb_pc_sse"):
        if not 0<=values[k]<=100:
            raise RuntimeError("Original hydroatlas percentage coding not comparable 0..100 "+k)
    if values["dis_m3_pyr"]<0:raise RuntimeError("Original natural basin outlet discharge negative")
    landuse.append(values)
    rawflows.append(values["dis_m3_pyr"])
    source=np_by[(key,x["PERIOD"])]
    if not source["source_valid"]:raise RuntimeError("Previously frozen actual EEA reported upstream source missing")
    p=float(source["source_year_reported_P_mgL_equivalent_site_equal_mean"])
    n=float(source["source_year_reported_N_mgL_equivalent_site_equal_mean"])
    if not all(math.isfinite(i) and i>=0 for i in (p,n)):
        raise RuntimeError("Authentic original annual N/P source nonfinite")
    np_parts.append({
       "N_log1p_real_reported_mgL_full_source_delivery":math.log1p(n),
       "P_log1p_real_reported_ugL_full_source_delivery":math.log1p(p*1000),
       "source_reported_any_positive":int(source["source_station_annual_upstream_NP_positive_n"]>0)
    })
if len({x["HYBAS_ID_L07"] for x in rows})!=232:
    raise RuntimeError("Original 232 native H7 study basins changed")
if sum(x["source_reported_any_positive"] for x in np_parts)!=56:
    raise RuntimeError("Original actual annual reported source positive H7 subset changed")
# Per-basin static landcover exactly identical across periods
per_key={}
for i,r in enumerate(rows):
    key=r["HYBAS_ID_L07"]
    if key in per_key and per_key[key]!=landuse[i]:
        raise RuntimeError("Static original BasinATLAS landcover changing by observational year")
    per_key[key]=landuse[i]
N=len(rows)
countries=sorted({r["country"] for r in rows})
methods=sorted({r["original_majority_BQE_method"] for r in rows})
base_num=["log_local_H7_area_km2","log1p_original_EQR_site_n",
 "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
 "log1p_GBIF_5family_butterfly_source_records","log1p_original_EEA_site_DIS_AV_CMS"]
B=np.column_stack([
 np.ones(N),
 np.array([[float(r[k]) for k in base_num] for r in rows]),
 np.array([[float(r["country"]==k) for k in countries[1:]] for r in rows]),
 np.array([[float(r["original_majority_BQE_method"]==k) for k in methods[1:]] for r in rows]),
 np.array([float(r["PERIOD"]=="2020-2022") for r in rows])
])
W=np.array([float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in rows])
crop_up=np.array([r["crp_pc_use"] for r in landuse])
urban_up=np.array([r["urb_pc_use"] for r in landuse])
crop_local=np.array([r["crp_pc_sse"] for r in landuse])
urban_local=np.array([r["urb_pc_sse"] for r in landuse])
natural_flow=np.log1p(np.array(rawflows))
P=np.array([r["P_log1p_real_reported_ugL_full_source_delivery"] for r in np_parts])
R=np.array([r["N_log1p_real_reported_mgL_full_source_delivery"] for r in np_parts])
env_up=np.column_stack([crop_up,urban_up,natural_flow])
env_local=np.column_stack([crop_local,urban_local,natural_flow])
X={
 "B0_source_recording_country_method_flow":B,
 "B1_plus_STATIC_HydroWASTE":np.column_stack([B,W]),
 "B2_plus_UPSTREAM_crop_urban_natural_flow":np.column_stack([B,env_up]),
 "B3_UPSTREAM_env_plus_STATIC_HydroWASTE":np.column_stack([B,env_up,W]),
 "B4_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N":np.column_stack([B,env_up,W,R]),
 "B5_UPSTREAM_env_STATIC_plus_ANNUAL_reported_P":np.column_stack([B,env_up,W,P]),
 "B6_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N_and_P":np.column_stack([B,env_up,W,R,P]),
 "B2local_plus_LOCAL_crop_urban_natural_flow":np.column_stack([B,env_local]),
 "B3local_LOCAL_env_plus_STATIC_HydroWASTE":np.column_stack([B,env_local,W]),
 "B4local_LOCAL_env_STATIC_plus_ANNUAL_N":np.column_stack([B,env_local,W,R]),
 # Prior PR66 original source N baseline (B1+R), needs exact numerical check
 "PR66_B1_plus_ANNUAL_reported_N":np.column_stack([B,W,R])
}
basins=[r["HYBAS_ID_L07"] for r in rows]
folds=np.array([int(hashlib.sha256(str(b).encode()).hexdigest()[:8],16)%5 for b in basins])
if min(np.bincount(folds,minlength=5))<35:raise RuntimeError("H7 original fivefold study selection changed")
group_names,inv=np.unique(basins,return_inverse=True)
G=len(group_names)
if G!=232:raise RuntimeError("Original 232 independent H7 units changed")
rng=np.random.default_rng(20261112)
boot_ids=rng.integers(0,G,size=(2000,G))
def predict(A,Y):
    pred=np.empty(len(Y))
    for fold in range(5):
        test=folds==fold;train=~test
        Z=A.copy()
        mu=Z[train,1:].mean(axis=0)
        sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        tr=Z[train];te=Z[test]
        reg=np.eye(tr.shape[1])*10.;reg[0,0]=0.
        beta=np.linalg.solve(tr.T@tr+reg,tr.T@Y[train])
        pred[test]=te@beta
    return pred
def ols(A,Y,key_weights,which):
    Z=A.copy();mu=Z[:,1:].mean(axis=0);sd=Z[:,1:].std(axis=0)
    sd=np.where(sd>1e-12,sd,1.)
    Z[:,1:]=(Z[:,1:]-mu)/sd
    AA=Z*np.sqrt(key_weights)[:,None]
    YY=Y*np.sqrt(key_weights)
    coefficients=np.linalg.lstsq(AA,YY,rcond=None)[0]
    u=YY-AA@coefficients
    bread=np.linalg.pinv(AA.T@AA)
    meat=np.zeros((AA.shape[1],AA.shape[1]))
    for g in range(G):
        sel=inv==g
        v=AA[sel].T@u[sel]
        meat+=np.outer(v,v)
    rank=np.linalg.matrix_rank(AA)
    V=bread@meat@bread*G/(G-1)*(N-1)/(N-rank)
    se=np.sqrt(np.maximum(V.diagonal(),0))
    idx=A.shape[1]-1
    return {
      "source_variable":which,
      "beta_response_units_per_1SD_predictor":float(coefficients[idx]),
      "H7_cluster_descriptive95":[float(coefficients[idx]-1.96*se[idx]),float(coefficients[idx]+1.96*se[idx])],
      "source_original_predictor_SD":float(sd[idx-1]),
      "rank":int(rank),"number_X_columns":int(AA.shape[1])
    }
countrycounts=collections.Counter(r["country"] for r in rows)
cw=np.array([1/countrycounts[r["country"]] for r in rows],dtype=float)
cw/=cw.mean()
pairs=[
 ("B0_source_recording_country_method_flow","B1_plus_STATIC_HydroWASTE"),
 ("B0_source_recording_country_method_flow","B2_plus_UPSTREAM_crop_urban_natural_flow"),
 ("B1_plus_STATIC_HydroWASTE","B3_UPSTREAM_env_plus_STATIC_HydroWASTE"),
 ("B2_plus_UPSTREAM_crop_urban_natural_flow","B3_UPSTREAM_env_plus_STATIC_HydroWASTE"),
 ("B3_UPSTREAM_env_plus_STATIC_HydroWASTE","B4_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N"),
 ("B3_UPSTREAM_env_plus_STATIC_HydroWASTE","B5_UPSTREAM_env_STATIC_plus_ANNUAL_reported_P"),
 ("B3_UPSTREAM_env_plus_STATIC_HydroWASTE","B6_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N_and_P"),
 ("B2local_plus_LOCAL_crop_urban_natural_flow","B3local_LOCAL_env_plus_STATIC_HydroWASTE"),
 ("B3local_LOCAL_env_plus_STATIC_HydroWASTE","B4local_LOCAL_env_STATIC_plus_ANNUAL_N")
]
responses={
 "EEA_original_benthic_EQR":("EQR_original_site_equal_basin_period_mean","actual_EEA_EQR_original_source","original_measured_EEA_benthic_EQR"),
 "GBIF_original_plant_speciesKey_S50":("plant_S50_GBIF_speciesKey","GBIF_Plantae_speciesKey_S50_50rarefaction","GBIF_original_plant_S50"),
 "GBIF_original_Lycaenidae_fivefamily_logodds":("lycaenidae_fivefamily_index_share_logodds","GBIF_Lycaenidae_log_odds_5_butterfly_families","GBIF_five_family_Lycaenidae_log_odds")
}
final={}
for k,(target,prior_key,n_key) in responses.items():
    Y=np.array([float(row[target]) for row in rows])
    err={model:(Y-predict(A,Y))**2 for model,A in X.items()}
    MSE={x:float(np.mean(v)) for x,v in err.items()}
    # Mandatory exact numeric previous original model reproduction
    originalB=prior["results_three_response_separate_models"][prior_key]
    if not math.isclose(MSE["B0_source_recording_country_method_flow"],
      originalB["baseline_5fold_hash_basin_disjoint_MSE"],rel_tol=1e-11,abs_tol=1e-12):
        raise RuntimeError("Original PR65 baseline not numerically reproduced "+k)
    if not math.isclose(MSE["B1_plus_STATIC_HydroWASTE"],
      originalB["plus_static_modelled_WWTP_gradient_5fold_MSE"],rel_tol=1e-11,abs_tol=1e-12):
        raise RuntimeError("Original PR65 static model not numerically reproduced "+k)
    originalNP=prior_annual["models"][n_key]["original_same_basins_source_P_N_CV_MSE"]
    if not math.isclose(MSE["PR66_B1_plus_ANNUAL_reported_N"],
      originalNP["B_plus_annual_reported_N_flow_converted"],rel_tol=1e-11,abs_tol=1e-12):
        raise RuntimeError("Original PR66 annual N model not reproduced "+k)
    comparisons={}
    for a,b in pairs:
        q0=np.bincount(inv,weights=err[a],minlength=G)
        q1=np.bincount(inv,weights=err[b],minlength=G)
        boots0=np.sum(q0[boot_ids],axis=1)
        boots1=np.sum(q1[boot_ids],axis=1)
        pct=100*(boots0-boots1)/boots0
        comparisons[a+"→"+b]={
           "original_same_basin_5fold_MSE_gain_pct":100*(MSE[a]-MSE[b])/MSE[a],
           "H7_cluster2000_bootstrap_descriptive95_gain_pct":np.percentile(pct,[2.5,97.5]).tolist(),
           "bootstrap_positive_gain_fraction":float(np.mean(pct>0))}
    slopes={
      "B3_static_after_UPSTREAM_env":ols(X["B3_UPSTREAM_env_plus_STATIC_HydroWASTE"],Y,np.ones(N),"W_static"),
      "B4_actual_reported_N_after_UPSTREAM_env_and_static_W":ols(X["B4_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N"],Y,np.ones(N),"N_reported"),
      "B4_N_country_equal_total_weight":ols(X["B4_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N"],Y,cw,"N_reported"),
      "B5_reported_P_after_UPSTREAM_env_and_static_W":ols(X["B5_UPSTREAM_env_STATIC_plus_ANNUAL_reported_P"],Y,np.ones(N),"P_reported"),
      "B3local_static_after_LOCAL_env":ols(X["B3local_LOCAL_env_plus_STATIC_HydroWASTE"],Y,np.ones(N),"W_static"),
      "B4local_reported_N_after_LOCAL_env":ols(X["B4local_LOCAL_env_STATIC_plus_ANNUAL_N"],Y,np.ones(N),"N_reported")
    }
    final[k]={"N":N,"basins":G,
              "baseline_original_PR65_PR66_MSE_values_exactly_reproduced":True,
              "response_mean":float(Y.mean()),
              "all_11_fixed_model_fivefold_MSE":MSE,
              "predetermined_nested_model_comparisons":comparisons,
              "land_adjusted_source_gradients_H7_cluster_descriptive_OLS":slopes}
report={
 "schema":"GESUI_EU454_H7_ORIGINAL_THREE_RESPONSE_REPORTED_NP_LANDCOVER_ADJUSTED_v1",
 "genuine_source_original_HydroATLAS_1219_SHA256":ATLAS_SHA,
 "genuine_source_original_EEA2018_22_NP_454_SHA256":NP_SHA,
 "source_original_454_spatial_model_type":OLD_SPATIAL,
 "source_original_454_annual_NP_model_type":OLD_ANNUAL,
 "source_cohort":{"original_EU_H7_basin_periods":N,"original_independent_H7":G,
   "original_2017_2019_periods":sum(r["PERIOD"]=="2017-2019" for r in rows),
   "original_2020_2022_periods":sum(r["PERIOD"]=="2020-2022" for r in rows),
   "countries":dict(countrycounts),
   "positive_reported_annual_source_periods":sum(r["source_reported_any_positive"] for r in np_parts),
   "original_all_HydroATLAS_H7_IDs_matched":True,
   "original_HydroATLAS_static_upstream_cropland_percent_range":[float(crop_up.min()),float(crop_up.max())],
   "original_HydroATLAS_static_upstream_urban_percent_range":[float(urban_up.min()),float(urban_up.max())],
   "original_HydroATLAS_static_natural_outlet_flow_range_m3s":[float(min(rawflows)),float(max(rawflows))],
   "recorded_positive_source_absence_means_no_positive_MATCHED_REPORT_only":True,
   "HydroATLAS_landuse_NA_count":0},
 "true_environment_covariates":{
    "upstream_cropland_crp_pc_use_percent_static":True,
    "upstream_urban_urb_pc_use_percent_static":True,
    "BasinATLAS_natural_dis_m3_pyr_m3s_original_static":True,
    "local_cropland_urban_source_secondary_only":True,
    "environment_is_historically_static_not_real_2018_2022_landuse_change":True},
 "separate_layer_result_models":final,
 "limitations":[
  "HydroATLAS landcover and natural discharge are historical static sources and have different source periods than 2018/2022 original EEA monitoring and contemporary GBIF.",
  "HydroWASTE reported registered WWTP flow to original H7 BASIN OUTLET is not actual upstream exposure at an interior EEA field station.",
  "EEA report N/P annual matched-source mass is only positive for 56/454 H7-periods; 398 no-positive recorded inventory is not verified negative controls.",
  "All 3 outcomes are different measurements: actual EEA quality index, GBIF source-indexed rarefied plant speciesKey richness, GBIF five-family Lycaenidae composition. None is real pollination or ecosystem service.",
  "Nearby H7 basins can appear in different H7-hash CV folds and country/method FE are spatially nonindependent. Bootstrap intervals descriptive and source geography may not transfer.",
  "This is a POSTHOC landuse robustness study after original WWTP–EQR negative sign was seen, no causal nutrient conclusion regardless of coefficient sign.",
  "Source original BasinATLAS full 4.276GB ZIP publisher whole-file MD5 is metadata only, not independently checked; original level07 DBF compressed member CRC32 and 80MB uncompressed SHA fully verified."
 ]}
outpath=OUT/"H7_454_EU_THREE_LAYER_STATIC_WWTP_REPORTED_NP_HYDROATLAS_CROP_URBAN_NATURAL_FLOW_RESULTS.json"
outpath.write_text(json.dumps(report,sort_keys=True,indent=2,allow_nan=False))
print("H7_HYDROATLAS_454_LANDUSE_ADJUSTED_RESULTS_START",flush=True)
brief={k:v for k,v in report.items() if k!="limitations"}
print(json.dumps(brief,sort_keys=True,allow_nan=False),flush=True)
print("H7_HYDROATLAS_454_LANDUSE_ADJUSTED_RESULTS_END",flush=True)
