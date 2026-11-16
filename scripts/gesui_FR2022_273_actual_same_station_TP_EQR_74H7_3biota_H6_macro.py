#!/usr/bin/env python3
"""Real same-2022 monitored receiving river chemistry TP vs EEA benthic EQR
on genuine TP↔EEA river site pairs, plus same HydroBASINS7 GBIF plant
indexed rarefied S50 and Lycaenidae five-family observation composition.

France-only source-locked 273 field station pairs /74 H7 /24 H6 parents.
Original source selection and H6 fold assignment frozen BEFORE numeric
TP/EQR models. This is observational posthoc, NOT WWTP-causal effects nor
measured ecosystem services.
"""
import collections,csv,hashlib,io,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2022_TRUE_TP_EQR_3BIOTA_MACRO"
ROOT.mkdir(parents=True,exist_ok=True)
S={
 "source_match":"448add2d6850cbc876c3e6b90bca278c6b711e64f9d363071e93bc0c18b84611",
 "h6_source_fold":"294eb69b76cd4cf40f1022a5e83278529b1d2bbceac9a38931004610a23eb052",
 "atlas":"ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686",
 "TP_git":"a75a4cf140de64fd6038079a3fac6756a461f3f8",
 "EQ_git":"325f199e985d9862ffee5b0edb205c61a7be3833"}
def file_json(env,sha):
    path=pathlib.Path(os.environ[env])
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=sha:raise RuntimeError("Original source SHA mismatch "+env)
    return json.loads(raw)
match=file_json("ORIGINAL_TP_EQR_SAME_REACH_250M_SOURCE",S["source_match"])
foldjson=file_json("ORIGINAL_TP_EQR_H6_FROZEN_FOLDS",S["h6_source_fold"])
atlas=file_json("ORIGINAL_H7_HYDROATLAS_1219_ATTRIBUTES",S["atlas"])["original_1219_H7_covariate_by_real_ID"]
pairs=match["original_454_H7_cohort_2022_genuine_TP_250m_source_matched_pairs"]
if len(pairs)!=273 or len({z["EEA_H7_L07"] for z in pairs})!=74:
    raise RuntimeError("Genuine 273 same-station pairs source cohort drift")
if any(z["EEA2022_country"]!="FR" or z["EEA2022_BQE_original_method"]!="I2M2" for z in pairs):
    raise RuntimeError("Original country and assessment method fixed single-source frame changed")
folds=foldjson["matched_original_74_H7_fold_source"]
if len(folds)!=74:raise RuntimeError("Genuine original H6 source fold n74 differs")
fold_by={z["H7"]:z["fold"] for z in folds}
parent_by={z["H7"]:z["H6"] for z in folds}
paircount_by={z["H7"]:z["n_2022_same_station_pairs"] for z in folds}
if len(set(parent_by.values()))!=24:
    raise RuntimeError("Genuine original 24 independent H6 watershed parent groups changed")
def gitsha(blob):
    return hashlib.sha1(("blob "+str(len(blob))+"\0").encode()+blob).hexdigest()
def datafile(name,original_sha):
    raw=pathlib.Path(name).read_bytes()
    if gitsha(raw)!=original_sha:raise RuntimeError("Original truly observed site Git blob changed "+name)
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline="")))
wq=datafile("results/eea2022_measured_nutrient_load_river_20261107/TP_2022_reported_wwtp_load_and_measured_river_station.csv",S["TP_git"])
eea=datafile("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",S["EQ_git"])
if len(wq)!=1291 or len(eea)!=2475:raise RuntimeError("Actual measured TP or original EQR registry sizes changed")
W={r["WQ.ID"]:r for r in wq}; E={r["site_key"]:r for r in eea}
if len(W)!=1291 or len(E)!=2475:raise RuntimeError("Duplicate original WQ/EEA monitoring ID")
with open(os.environ["ORIGINAL_FIXED_454_H7_BIODIVERSITY_ROWS"],encoding="utf-8-sig",newline="") as f:
    biology=list(csv.DictReader(f))
if len(biology)!=454:raise RuntimeError("Original 454 study basin/period frame changed")
macro={r["HYBAS_ID_L07"]:r for r in biology if r["PERIOD"]=="2020-2022"}
if len(macro)!=227:raise RuntimeError("Original source macro 2020–22 H7 population changed")
site_rows=[]; invalid=collections.Counter()
for source in pairs:
    ww=W[source["WQ_ID"]]
    eq=E[source["EEA2022_site_key"]]
    if ww["HYRIV_ID"]!=source["HYRIV_ID"] or eq["HYRIV_ID"]!=source["HYRIV_ID"]:
        raise RuntimeError("True matched original TP/EQR river reach source changed")
    if eq["procedureClassificationSystem"]!="I2M2" or eq["countryCode"]!="FR":
        raise RuntimeError("Single French I2M2 method unexpectedly changed")
    tp=float(ww["real_2022_station_median_mgL"])
    eqr=float(eq["normalised_EQR"])
    flow=float(eq["DIS_AV_CMS"])
    if not all(math.isfinite(x) for x in (tp,eqr,flow)) or tp<=0 or flow<=0:
        invalid["invalid_real_2022_nutrient_biological_or_flow"]+=1
        continue
    if int(ww["actual_2022_observation_count"])<4 or int(ww["distinct_2022_months"])<3:
        invalid["original_WQ_year_measured_sample_support_mismatched"]+=1
        continue
    site_rows.append({
      "WQ_ID":source["WQ_ID"],"EEA_original_site_key":source["EEA2022_site_key"],
      "HYBAS_ID_L07":source["EEA_H7_L07"],"source_H6_parent":parent_by[source["EEA_H7_L07"]],
      "source_fold":fold_by[source["EEA_H7_L07"]],
      "original_WQ_station_real_TP_2022_mgL":tp,
      "log_real_measured_river_TP_2022_mgL":math.log(tp),
      "original_EEA_benthic_I2M2_2022_normalized_EQR":eqr,
      "original_EEA_site_modelled_river_DIS_AV_CMS":flow,
      "original_observed_TP_sample_count_2022":int(ww["actual_2022_observation_count"]),
      "original_observed_TP_distinct_months_2022":int(ww["distinct_2022_months"]),
      "original_source_station_coordinate_distance_m":source["spherical_source_coordinate_separation_m"]
    })
if len(site_rows)<250 or len({r["HYBAS_ID_L07"] for r in site_rows})<65:
    raise RuntimeError("Predeclared original numerical station-pair gate n>=250 and nH7>=65 failed. Do not fit")
if len({r["WQ_ID"] for r in site_rows})!=len(site_rows) or len({r["EEA_original_site_key"] for r in site_rows})!=len(site_rows):
    raise RuntimeError("Genuine same-station 2022 TP-EQR source not one-to-one")

by=collections.defaultdict(list)
for r in site_rows:by[r["HYBAS_ID_L07"]].append(r)
if set(by)!=set(fold_by):raise RuntimeError("Original source 74 H7 physically co-located study frame drift")
rows=[]
for hid in sorted(by):
    ss=by[hid]
    if len(ss)!=paircount_by[hid]:raise RuntimeError("Source original TP/EQR site membership unexpectedly dropped")
    orig=macro[hid];land=atlas[hid]
    if not all(land.get(k) is not None for k in ("crp_pc_use","urb_pc_use","dis_m3_pyr")):
        raise RuntimeError("Genuine atlas landcover/source flow original missing")
    modelW=float(orig["W_static_upstream_MODEL_flow_per_UpArea"])
    r={
      "HYBAS_ID_L07":hid,"HYBAS_ID_L06":parent_by[hid],"original_source_fold":fold_by[hid],
      "study_country":"FR","biological_method":"I2M2","year":2022,
      "original_matched_same_reach_TP_EQR_station_n":len(ss),
      "mean_log_real_2022_river_TP_mgL_at_exact_EQR_sites":float(np.mean([x["log_real_measured_river_TP_2022_mgL"] for x in ss])),
      "median_original_real_TP_2022_mgL_at_EQR_sites":float(np.median([x["original_WQ_station_real_TP_2022_mgL"] for x in ss])),
      "mean_actual_2022_TP_sample_count":float(np.mean([x["original_observed_TP_sample_count_2022"] for x in ss])),
      "mean_same_site_original_2022_EEA_I2M2_EQR":float(np.mean([x["original_EEA_benthic_I2M2_2022_normalized_EQR"] for x in ss])),
      "log_mean_original_site_DIS_AV_CMS":math.log1p(float(np.mean([x["original_EEA_site_modelled_river_DIS_AV_CMS"] for x in ss]))),
      "log_H7_original_local_area_km2":float(orig["log_local_H7_area_km2"]),
      "log1p_n_actual_matched_TP_EQR_river_sites":math.log1p(len(ss)),
      "source_original_upstream_cropland_percent":float(land["crp_pc_use"]),
      "source_original_upstream_urban_percent":float(land["urb_pc_use"]),
      "log1p_H7_natural_modelled_outlet_DIS_m3s":math.log1p(float(land["dis_m3_pyr"])),
      "GBIF_original_2020_2022_Plant_speciesKey_rarefied_S50":float(orig["plant_S50_GBIF_speciesKey"]),
      "GBIF_original_2020_2022_Lycaenidae_fivefamily_index_logodds":float(orig["lycaenidae_fivefamily_index_share_logodds"]),
      "GBIF_original_log_Plantae_record_opportunity":float(orig["log1p_GBIF_Plantae_source_records"]),
      "GBIF_original_log_Animalia_record_opportunity":float(orig["log1p_GBIF_Animalia_source_records"]),
      "GBIF_original_log_fivefamily_butterfly_records":float(orig["log1p_GBIF_5family_butterfly_source_records"]),
      "modelled_HydroWASTE_source_to_H7_BASIN_OUTLET_index":modelW
    }
    rows.append(r)
assert len(rows)==74 and len({r["HYBAS_ID_L06"] for r in rows})==24

Ykeys={
 "true_same_2022_I2M2_benthic_EQR":"mean_same_site_original_2022_EEA_I2M2_EQR",
 "H7_GBIF_plant_speciesKey_S50":"GBIF_original_2020_2022_Plant_speciesKey_rarefied_S50",
 "H7_GBIF_Lycaenidae_fivefamily_logodds":"GBIF_original_2020_2022_Lycaenidae_fivefamily_index_logodds"
}
base=[
 "log_H7_original_local_area_km2",
 "source_original_upstream_cropland_percent",
 "source_original_upstream_urban_percent",
 "log1p_H7_natural_modelled_outlet_DIS_m3s",
 "log_mean_original_site_DIS_AV_CMS",
 "log1p_n_actual_matched_TP_EQR_river_sites"]
xTP=np.array([z["mean_log_real_2022_river_TP_mgL_at_exact_EQR_sites"] for z in rows])
xW=np.array([z["modelled_HydroWASTE_source_to_H7_BASIN_OUTLET_index"] for z in rows])
fold=np.array([z["original_source_fold"] for z in rows])
H6,inv=np.unique([z["HYBAS_ID_L06"] for z in rows],return_inverse=True)
if len(H6)!=24 or min(np.bincount(fold,minlength=5))<4:
    raise RuntimeError("Original pre-outcome 24 H6 geographic/source fold support changed")
rng=np.random.default_rng(20261113)
bootids=rng.integers(0,len(H6),(2000,len(H6)))
def evaluate(X,y):
    pred=np.empty(len(y))
    for f in range(5):
        train=fold!=f;test=fold==f
        tr=X[train];ts=X[test]
        mu=tr.mean(axis=0);sd=tr.std(axis=0);sd=np.where(sd>1e-12,sd,1.)
        model=Ridge(alpha=10.,fit_intercept=True,solver="svd")
        model.fit((tr-mu)/sd,y[train])
        pred[test]=model.predict((ts-mu)/sd)
    if not np.all(np.isfinite(pred)):raise RuntimeError("Original source H6 geographic CV has missing predictions")
    return (y-pred)**2
def report_gain(e0,e1):
    g=100*(float(np.mean(e0))-float(np.mean(e1)))/float(np.mean(e0))
    block0=np.bincount(inv,weights=e0,minlength=len(H6))
    block1=np.bincount(inv,weights=e1,minlength=len(H6))
    v0=np.sum(block0[bootids],axis=1);v1=np.sum(block1[bootids],axis=1)
    d=100*(v0-v1)/v0
    return {"H6_parent_original_GroupKFold5_MSE_gain_pct":g,
        "H6_parent_2000_bootstrap_descriptive95_MSE_gain_pct":np.percentile(d,[2.5,97.5]).tolist(),
        "H6_parent_2000_bootstrap_positive_fraction":float(np.mean(d>0)),
        "original_baseline_MSE":float(e0.mean()),"original_additional_MSE":float(e1.mean())}
def H6_OLS(A,y,which):
    n=len(y);Z=np.column_stack([np.ones(n),A])
    sd=Z[:,1:].std(axis=0);mu=Z[:,1:].mean(axis=0)
    sd=np.where(sd>1e-12,sd,1)
    Z[:,1:]=(Z[:,1:]-mu)/sd
    coef=np.linalg.lstsq(Z,y,rcond=None)[0];resid=y-Z@coef
    bread=np.linalg.pinv(Z.T@Z);meat=np.zeros_like(bread)
    for j in range(len(H6)):
        mask=inv==j
        v=Z[mask].T@resid[mask];meat+=np.outer(v,v)
    rank=np.linalg.matrix_rank(Z)
    cov=bread@meat@bread*(len(H6)/(len(H6)-1))*(n-1)/(n-rank)
    se=float(np.sqrt(max(0,cov[-1,-1])))
    return {"predictor":which,
       "adjusted_beta_per_1SD_log_real_2022_TP":float(coef[-1]),
       "H6_parent_cluster_descriptive95":[float(coef[-1]-1.96*se),float(coef[-1]+1.96*se)],
       "source_predictor_1SD_logTP":float(sd[-1]),
       "rank":rank}
res={}
for label,key in Ykeys.items():
    y=np.array([r[key] for r in rows],dtype=float)
    cols=list(base)
    if label=="H7_GBIF_plant_speciesKey_S50":
        cols+=["GBIF_original_log_Plantae_record_opportunity","GBIF_original_log_Animalia_record_opportunity"]
    elif label=="H7_GBIF_Lycaenidae_fivefamily_logodds":
        cols+=["GBIF_original_log_fivefamily_butterfly_records","GBIF_original_log_Animalia_record_opportunity"]
    X0=np.array([[float(z[c]) for c in cols] for z in rows],dtype=float)
    e0=evaluate(X0,y)
    e1=evaluate(np.column_stack([X0,xTP]),y)
    eW=evaluate(np.column_stack([X0,xW]),y)
    spr=spearmanr(xTP,y)
    res[label]={
       "source_original_number_74_H7":74,
       "n_original_same_2022_EEA_and_Peters_sites":len(site_rows),
       "raw_basin_level_spearman_rho_log_measured_TP_vs_bio":float(spr.statistic),
       "raw_basin_level_descriptive_p_not_independent_H6":float(spr.pvalue),
       "response_SD":float(y.std()),
       "source_baseline_covariates":cols,
       "after_landcover_and_river_flow_MEASURED_TP_incremental_H6_CV":report_gain(e0,e1),
       "after_landcover_and_river_flow_MODELED_WWTP_H7_OUTLET_incremental_H6_CV":report_gain(e0,eW),
       "measured_TP_OLS_after_environment_H6_cluster":H6_OLS(np.column_stack([X0,xTP]),y,"original_colocated_2022_measured_river_TP_log")
    }
# physically matched station-level TP ↔ EQR, distinct from H7 macro outcomes
s_tp=np.array([z["log_real_measured_river_TP_2022_mgL"] for z in site_rows])
s_eqr=np.array([z["original_EEA_benthic_I2M2_2022_normalized_EQR"] for z in site_rows])
r_site=spearmanr(s_tp,s_eqr)
pairs_by=collections.defaultdict(list)
for r in site_rows:pairs_by[r["HYBAS_ID_L07"]].append(r)
multi={k:items for k,items in pairs_by.items() if len(items)>=2}
demeanX=[];demeanY=[];within=[]
for h,items in multi.items():
    xs=np.array([x["log_real_measured_river_TP_2022_mgL"] for x in items])
    ys=np.array([x["original_EEA_benthic_I2M2_2022_normalized_EQR"] for x in items])
    vx=xs-xs.mean();vy=ys-ys.mean()
    demeanX.extend(vx);demeanY.extend(vy)
    within.append({"h7":h,"cross":float(vx@vy),"xx":float(vx@vx),"n":len(items)})
within_slope=float(sum(z["cross"] for z in within)/sum(z["xx"] for z in within))
H7_ids=sorted(pairs_by)
byh={z["h7"]:z for z in within}
within_rng=np.random.default_rng(20261113)
rboot=[]
for i in range(2000):
    sampled=within_rng.choice(H7_ids,size=len(H7_ids),replace=True)
    num=sum(byh[h]["cross"] for h in sampled if h in byh)
    denom=sum(byh[h]["xx"] for h in sampled if h in byh)
    if denom>0:rboot.append(num/denom)
within_src={
 "native_real_2022_site_pair_273_unadjusted_Spearman_rho":float(r_site.statistic),
 "native_real_2022_site_pair_unadjusted_descriptive_p_not_independent_H7":float(r_site.pvalue),
 "original_H7_with_2plus_exact_site_pairs":len(multi),
 "original_same_site_pairs_in_H7_with_2plus":sum(len(v) for v in multi.values()),
 "within_H7_site_demeaned_unadjusted_logTP_to_EQR_slope":within_slope,
 "within_H7_cluster2000_bootstrap_slope_descriptive95":np.percentile(rboot,[2.5,97.5]).tolist() if rboot else None,
 "no_independent_station_pvalue_claim":True
}
tp_modelled_relation=spearmanr(xTP,xW)
summary={
 "schema":"GESUI_FR_2022_273_SAME_RIVER_SITE_TP_I2M2_EQR_74_H7_24_H6_3BIOTA_RESULTS_v1",
 "matched_original_250m_TP_EQR_station_source_SHA256":S["source_match"],
 "native_74H7_H6_preoutcome_folds_SHA256":S["h6_source_fold"],
 "original_HydroATLAS_1219_real_source_SHA256":S["atlas"],
 "original_2022_Peters_TP_source_GitBlob_SHA1":S["TP_git"],
 "original_2022_EEA_I2M2_EQR_source_GitBlob_SHA1":S["EQ_git"],
 "original_FR_I2M2_2022_same_site_unique_pairs":len(site_rows),
 "original_FR_H7_basins":len(rows),
 "original_H6_hydrologic_parent_regions":len(H6),
 "source_predeclared_fold_74_H7":np.bincount(fold,minlength=5).tolist(),
 "new_TP_EQR_source_invalid_rows":dict(invalid),
 "actual_measured_2022_TP_mgL_range_at_same_sites":[float(np.exp(s_tp).min()),float(np.exp(s_tp).max())],
 "I2M2_EQR_range_at_same_sites":[float(s_eqr.min()),float(s_eqr.max())],
 "same_physical_river_station_TP_EEA_EQR_association":within_src,
 "original_same_H7_measured_TP_vs_modelled_old_HydroWASTE_outlet_Spearman":{
     "rho":float(tp_modelled_relation.statistic),"descriptive_p":float(tp_modelled_relation.pvalue)},
 "same_native_H7_2020_2022_response_results":res,
 "limitations":[
  "Only France I2M2 EQR supports true close-station TP overlap; NO pan-European or country-heldout transfer can be claimed.",
  "Original real measured TP year 2022 median has 4+ sampled visits 3+ months but not same sample day as EEA benthic assessment.",
  "GBIF flora and butterfly source units cover whole H7 polygons, not same physical aquatic station, true flowering or realized pollination.",
  "Registered HydroWASTE outlet wastewater volume is static modeled geography, not observed annual plant P emission, actual dilution or verified activity.",
  "Both measured TP and I2M2 EQR may reflect common catchment urbanization, agriculture, morphology, disturbance, substrate and monitoring selection; no randomized or WWTP-causal inference.",
  "The original 74 H7 H6 geographic parent groups were source-frozen with 16 of 232 original Europe L7 parent mapping dependent on explicit GEOS geometry repairs.",
  "This is a newly discovered same-site observational analysis after earlier related wastewater–EQR results were known; inference is exploratory POSTHOC."
 ],
 "recorded_noncausal_water_quality_and_three_biota_proxy_relationships":True
}
result=ROOT/"FR2022_EEA_I2M2_BENTHIC_REAL_PETERS_TP_273SITES_74H7_3BIOTA_H6_GROUPCV_RESULTS.json"
def num_default(v):
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,np.floating):return float(v)
    if isinstance(v,np.ndarray):return v.tolist()
    raise TypeError("Unknown numeric result "+str(type(v)))
result.write_text(json.dumps(summary,sort_keys=True,indent=2,default=num_default,allow_nan=False))
with (ROOT/"FR2022_REAL_TP_EQR_74_H7_BASIN_MEDIAN_3BIOTA_MACRO_DATA.csv").open("w",newline="",encoding="utf8") as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader();writer.writerows(rows)
with (ROOT/"FR2022_REAL_TP_I2M2_EQR_273_SAME_RIVER_STATION_MATCHED_DATA.csv").open("w",newline="",encoding="utf8") as f:
    writer=csv.DictWriter(f,fieldnames=list(site_rows[0]))
    writer.writeheader();writer.writerows(site_rows)
print("FR2022_TP_I2M2_EQR_3BIOTA_MACRO_RESULT_START",flush=True)
print(json.dumps({k:v for k,v in summary.items() if k!="limitations"},sort_keys=True,default=num_default,allow_nan=False),flush=True)
print("FR2022_TP_I2M2_EQR_3BIOTA_MACRO_RESULT_END",flush=True)
