#!/usr/bin/env python3
"""Independent NUMPY normal-equation ridge repro and strict H6-parent leaveout
of saved original real 2022 France 273 co-located measured TP / benthic I2M2
source sites and same 74 H7 GBIF plant/butterfly source observations.

Posthoc sensitivity after positive original TP↔EQR co-response already seen.
No novel field data or causal WWTP/realized service inference.
"""
import collections,csv,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2022_TP_EQR_INDEPENDENT"
ROOT.mkdir(parents=True,exist_ok=True)
original=json.loads(pathlib.Path(os.environ["ORIGINAL_FR2022_TP_EQR_MAIN_RESULTS"]).read_text())
if original["schema"]!="GESUI_FR_2022_273_SAME_RIVER_SITE_TP_I2M2_EQR_74_H7_24_H6_3BIOTA_RESULTS_v1":
    raise RuntimeError("Original TP↔EQR 2022 France source test schema changed")
with open(os.environ["ORIGINAL_FR2022_TP_EQR_H7_ROWS"],encoding="utf-8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
with open(os.environ["ORIGINAL_FR2022_TP_EQR_STATION_ROWS"],encoding="utf-8-sig",newline="") as f:
    sites=list(csv.DictReader(f))
if len(rows)!=74 or len(sites)!=273:
    raise RuntimeError("Original 74-H7/273 same-station source counts changed")
if len({r["HYBAS_ID_L07"] for r in rows})!=74 or len({r["HYBAS_ID_L06"] for r in rows})!=24:
    raise RuntimeError("Original same-basin source region group changed")
if set(z["original_source_fold"] for z in rows)!=set("01234"):
    raise RuntimeError("Original five source-only folds changed")
fold=np.array([int(r["original_source_fold"]) for r in rows])
group=np.array([r["HYBAS_ID_L06"] for r in rows])
originalTP=np.array([float(r["mean_log_real_2022_river_TP_mgL_at_exact_EQR_sites"]) for r in rows])
originalW=np.array([float(r["modelled_HydroWASTE_source_to_H7_BASIN_OUTLET_index"]) for r in rows])
groups,groupinv=np.unique(group,return_inverse=True)
if len(groups)!=24:raise RuntimeError("Independent original H6 geometry source changed")
base0=original["same_native_H7_2020_2022_response_results"]
response={
 "true_same_2022_I2M2_benthic_EQR":"mean_same_site_original_2022_EEA_I2M2_EQR",
 "H7_GBIF_plant_speciesKey_S50":"GBIF_original_2020_2022_Plant_speciesKey_rarefied_S50",
 "H7_GBIF_Lycaenidae_fivefamily_logodds":"GBIF_original_2020_2022_Lycaenidae_fivefamily_index_logodds"}
def predict(X,Y,folds):
    n=len(Y);pred=np.empty(n)
    for k in sorted(set(folds)):
        test=(folds==k);train=~test
        Z=X.copy()
        mu=Z[train].mean(axis=0);std=Z[train].std(axis=0)
        std=np.where(std>1e-12,std,1.0)
        a=(Z[train]-mu)/std;b=(Z[test]-mu)/std
        aa=np.column_stack([np.ones(a.shape[0]),a])
        bb=np.column_stack([np.ones(b.shape[0]),b])
        pen=np.eye(aa.shape[1])*10.;pen[0,0]=0
        beta=np.linalg.solve(aa.T@aa+pen,aa.T@Y[train])
        pred[test]=bb@beta
    if not np.all(np.isfinite(pred)):raise RuntimeError("Independent original source CV prediction nan")
    return (Y-pred)**2
bootrng=np.random.default_rng(20261113)
ids=bootrng.integers(0,len(groups),(2000,len(groups)))
def gains(e0,e1):
    bas=float(e0.mean());plus=float(e1.mean())
    e0h=np.bincount(groupinv,weights=e0,minlength=len(groups))
    e1h=np.bincount(groupinv,weights=e1,minlength=len(groups))
    total0=e0h[ids].sum(axis=1);total1=e1h[ids].sum(axis=1)
    dist=100*(total0-total1)/total0
    return {"relative_OOS_MSE_gain_percent":100*(bas-plus)/bas,
            "first_original_MSE":bas,"additional_model_MSE":plus,
            "H6_parent_2000_bootstrap_descriptive95_percent":np.percentile(dist,[2.5,97.5]).tolist(),
            "fraction_group_bootstrap_positive":float(np.mean(dist>0))}
output={}
for label,key in response.items():
    prior=base0[label]
    covars=prior["source_baseline_covariates"]
    y=np.array([float(z[key]) for z in rows])
    X0=np.array([[float(z[c]) for c in covars] for z in rows])
    Xtp=np.column_stack([X0,originalTP])
    Xw=np.column_stack([X0,originalW])
    Xboth=np.column_stack([X0,originalTP,originalW])
    r0=predict(X0,y,fold)
    rtp=predict(Xtp,y,fold)
    rw=predict(Xw,y,fold)
    for a,b in [
       (r0.mean(),prior["after_landcover_and_river_flow_MEASURED_TP_incremental_H6_CV"]["original_baseline_MSE"]),
       (rtp.mean(),prior["after_landcover_and_river_flow_MEASURED_TP_incremental_H6_CV"]["original_additional_MSE"]),
       (rw.mean(),prior["after_landcover_and_river_flow_MODELED_WWTP_H7_OUTLET_incremental_H6_CV"]["original_additional_MSE"])]:
        if not math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-12):
            raise RuntimeError("Independent NUMPY ridge failed original scikit-learn 2022 TP/I2M2 ecology source replay "+label)
    x=group
    loho0=predict(X0,y,x);lohoTP=predict(Xtp,y,x);lohoW=predict(Xw,y,x);lohoBoth=predict(Xboth,y,x)
    output[label]={
       "independent_exact_Numpy_Ridge_REPRO_3_original_H6_5fold_MSE":True,
       "source_original_hydrological_5fold_gain_measured_TP":gains(r0,rtp),
       "source_original_hydrological_5fold_gain_modelled_registered_WWTP":gains(r0,rw),
       "posthoc_24fold_leave_each_REAL_H6_parent_OUT_measured_TP_increment":gains(loho0,lohoTP),
       "posthoc_24fold_leave_each_REAL_H6_parent_OUT_modelled_H7_outlet_increment":gains(loho0,lohoW),
       "posthoc_24fold_leave_each_REAL_H6_parent_OUT_add_WWTP_after_measured_TP":gains(lohoTP,lohoBoth),
       "posthoc_24fold_leave_each_REAL_H6_parent_OUT_add_measured_TP_after_WWTP":gains(lohoW,lohoBoth),
       "posthoc_24fold_group_out_MSE":{
           "B0":float(loho0.mean()),"B_TP":float(lohoTP.mean()),
           "B_WWTP":float(lohoW.mean()),"B_TP_WWTP":float(lohoBoth.mean())}
    }
# Measured in-stream TP source as RESPONSE to modelled WWTP registries, not
# chemical effluent causation.
sourceCovars=base0["true_same_2022_I2M2_benthic_EQR"]["source_baseline_covariates"]
XB0=np.array([[float(z[c]) for c in sourceCovars] for z in rows])
Xp=np.column_stack([XB0,originalW])
OOS0=predict(XB0,originalTP,fold);OOS1=predict(Xp,originalTP,fold)
LOHO0=predict(XB0,originalTP,group);LOHO1=predict(Xp,originalTP,group)
s=float(spearmanr(originalTP,originalW).statistic)
if not math.isclose(s,original["original_same_H7_measured_TP_vs_modelled_old_HydroWASTE_outlet_Spearman"]["rho"],rel_tol=1e-12):
    raise RuntimeError("Original TP-source index Spearman not reproducible")
siteTP=np.array([float(z["log_real_measured_river_TP_2022_mgL"]) for z in sites])
siteEQ=np.array([float(z["original_EEA_benthic_I2M2_2022_normalized_EQR"]) for z in sites])
rawR=float(spearmanr(siteTP,siteEQ).statistic)
if not math.isclose(rawR,original["same_physical_river_station_TP_EEA_EQR_association"]["native_real_2022_site_pair_273_unadjusted_Spearman_rho"],rel_tol=1e-12):
    raise RuntimeError("Original 273 genuine physical same-site TP/EQR correlation changed")
fromcollections=collections.defaultdict(list)
for item in sites:fromcollections[item["HYBAS_ID_L07"]].append(item)
sum_xy=0.;sum_xx=0.
for hid,g in fromcollections.items():
    if len(g)<2:continue
    v=np.array([float(z["log_real_measured_river_TP_2022_mgL"]) for z in g])
    w=np.array([float(z["original_EEA_benthic_I2M2_2022_normalized_EQR"]) for z in g])
    v=v-v.mean();w=w-w.mean()
    sum_xy+=float(v@w);sum_xx+=float(v@v)
slope=sum_xy/sum_xx
if not math.isclose(slope,original["same_physical_river_station_TP_EEA_EQR_association"]["within_H7_site_demeaned_unadjusted_logTP_to_EQR_slope"],rel_tol=1e-12):
    raise RuntimeError("Original within-H7 true 2022 TP/EQR station-paired coefficient changed")
sourceTP_W={
 "measured_TP_log_and_modelled_static_H7_outlet_Spearman_rho":s,
 "posthoc_original_2022_FR_TP_measured_response_to_WWTP_modelled_outlet_24H6_group_5fold":gains(OOS0,OOS1),
 "posthoc_original_2022_FR_TP_measured_response_to_WWTP_modelled_outlet_24H6_leaveout":gains(LOHO0,LOHO1),
 "no_actual_reported_WWTP_P_flow_at_these_all_source_stations":True}
out={
 "schema":"GESUI_FR2022_REAL_273_MEASURED_TP_I2M2_EQR_74_H7_24_H6_3BIOTA_NUMPY_INDEPENDENT_LOHO_v1",
 "original_source_numeric_study_type":original["schema"],
 "source_frozen_TP_same_station_2022_h7_sha256":original["matched_original_250m_TP_EQR_station_source_SHA256"],
 "source_frozen_h6_74H7_fold_sha256":original["native_74H7_H6_preoutcome_folds_SHA256"],
 "source_original_273_same_station_pairs":len(sites),
 "source_original_74_H7_catchments":len(rows),
 "source_original_24_H6_hydrologic_independent_blocks":len(groups),
 "original_source_273_site_Spearman_and_within_H7_demeaned_slope_independently_reproduced":True,
 "independent_3biota_3prior_H6_fivefold_MSE_exactly_reproduced":True,
 "posthoc_checks_per_three_observation_models":output,
 "posthoc_observed_TP_link_to_modelled_WWTP_source_geography":sourceTP_W,
 "limitations":[
   "Independent NUMPY numeric code replay, no independent country or field study.",
   "Posthoc 24-H6-leaveout evaluations and WTP+TP two-source contrasts were defined after the original TP↔EQR negative observation.",
   "All 74 source basins in France, one field biological method I2M2. No panEuropean or global transfer.",
   "Peters observed 2022 river TP is a 2022 annual median across sampled months, not same-day benthic chemistry.",
   "Modeled registered HydroWASTE wastewater outlet source is not actual operative effluent P load at the EEA station and cannot isolate agricultural P or landuse effects.",
   "GBIF indexed plant S50 and Lycaenidae 5-family relative counts do not measure floral visitors, native host, pollen transfer or ecosystem services."
 ]
}
path=ROOT/"FR2022_TP_I2M2_EQR_74H7_INDEPENDENT_NUMPY_24H6_LOHO_AND_WWTP_SOURCE_QA.json"
path.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False,default=lambda v:float(v) if isinstance(v,np.floating) else (int(v) if isinstance(v,np.integer) else v.tolist())))
small={k:v for k,v in out.items() if k!="limitations"}
print("FR2022_TP_I2M2_INDEPENDENT_LOHO_START",flush=True)
print(json.dumps(small,sort_keys=True,allow_nan=False,default=lambda v:float(v) if isinstance(v,np.floating) else (int(v) if isinstance(v,np.integer) else v.tolist())),flush=True)
print("FR2022_TP_I2M2_INDEPENDENT_LOHO_END",flush=True)
