#!/usr/bin/env python3
"""Comparative posthoc three-outcome macro analysis on authentic same
EU HydroBASINS L7 basin×period: original EEA 2018/2022 EQR, GBIF plant
species-key S50 and 5-family Lycaenidae share. ONE static HydroWASTE source.

Source-only 234 basin x2 WKT, 455 eligible periods precommitted. No chain
mediation or sewage-caused service claim. H7 hash CV is basin-disjoint
(NOT geographically blocked across neighboring basin IDs).
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np
from scipy.stats import spearmanr

OUT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_EU234_TRILAYER_STATIC_WWTP"
OUT.mkdir(parents=True,exist_ok=True)
SOURCE_SHA="07fa9744cd4e9574587a7b25c4ed4fa7f3183e377cac61c28723a4f24a86c875"
H7_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
BLOBS={
 "EEA18":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "EEA22":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
def fromsha(p,sha):
    b=pathlib.Path(p).read_bytes()
    if hashlib.sha256(b).hexdigest()!=sha:raise RuntimeError("Original source SHA changed: "+str(p))
    return json.loads(b)
C=fromsha(os.environ["EU234_FROZEN_CANDIDATES"],SOURCE_SHA)
H=fromsha(os.environ["HYBAS7_ORIGINAL_NATIVE_RIVER_SOURCE"],H7_SHA)
candidates=C["original_source_eligible_basin_period_candidates"]
h7map=H["candidate_source_site_period_rows"]
if len(candidates)!=468 or len(h7map)!=6330:raise RuntimeError("Original source-frozen cohort mismatched")
cmap={(a["HYBAS_ID_L07"],a["PERIOD"]):a for a in candidates}
if len(cmap)!=468:raise RuntimeError("Duplicate original fixed basin-period")
def gsha(b):
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()
E={}
for year,(path,sha) in BLOBS.items():
    b=pathlib.Path(path).read_bytes()
    if gsha(b)!=sha:raise RuntimeError("Original EEA source Git blob changed")
    rows=list(csv.DictReader(b.decode("utf-8-sig").splitlines()))
    if len(rows)!={"EEA18":2634,"EEA22":2475}[year]:
        raise RuntimeError("Original EEA station-year row count changed")
    E[year]={a["site_key"]:a for a in rows}
    if len(E[year])!=len(rows):raise RuntimeError("Duplicate station ID within a biology year")

sourcesites=collections.defaultdict(list)
for a in h7map:
    panel=a["panel"];hid=a["HYBAS_ID_L07"]
    if hid is None or panel not in ("EEA18","EEA22"):continue
    k=(str(hid),"2017-2019" if panel=="EEA18" else "2020-2022")
    if k in cmap:
        sourcesites[k].append(a["native_id"])
def mostcommon(vals):
    cc=collections.Counter(vals)
    if not cc:raise RuntimeError("Source method or country absent")
    return sorted(cc, key=lambda k:(-cc[k],k))[0]
def finite(v):
    try:
        x=float(v)
        return x if math.isfinite(x) else None
    except (ValueError,TypeError):return None
rows=[];excluded=collections.Counter()
for (hid,period),s in sorted(cmap.items()):
    if not s["original_5family_plant_S50_and_butterfly20_source_eligible"]:
        excluded["original_plantS50_or_fivefamily20_source_support_missing"]+=1;continue
    required=["original_source_plant_S50","original_source_lycaenidae_records",
      "original_source_5family_butterfly_records","original_GBIF_Animalia_record_opportunity",
      "original_GBIF_Plantae_record_opportunity","source_local_area_km2","source_UP_AREA_KM2",
      "registered_static_WWTP_outlet_modelled_flow_m3perday"]
    if any(s.get(k) is None for k in required):
        excluded["source_numerical_original_background_missing"]+=1;continue
    key="EEA18" if period=="2017-2019" else "EEA22"
    original_E=E[key]
    idx=sourcesites.get((hid,period),[])
    if not idx:
        excluded["no_original_EEA_EQR_site_source_in_H7"]+=1;continue
    stations=[]
    for site in idx:
        if site not in original_E:raise RuntimeError("Original EEA station source mapped but native record absent")
        z=original_E[site]
        eqr=finite(z["numeric_EQR_2018" if key=="EEA18" else "normalised_EQR"])
        if eqr is None:
            excluded["non_numeric_original_EQR_one_station"]+=1;continue
        if abs(eqr)>5:raise RuntimeError("Original EQR source not same numeric domain")
        flow=finite(z["DIS_AV_CMS"])
        stations.append((eqr,flow,z["countryCode"],z["procedureClassificationSystem"],site))
    if not stations:
        excluded["no_source_original_EQR_after_numeric_check"]+=1;continue
    nbut=float(s["original_source_5family_butterfly_records"])
    nly=float(s["original_source_lycaenidae_records"])
    if not (nbut>=20 and 0<=nly<=nbut):raise RuntimeError("Original butterfly family denominator invalid")
    wa=float(s["registered_static_WWTP_outlet_modelled_flow_m3perday"])
    area=float(s["source_local_area_km2"]);up=float(s["source_UP_AREA_KM2"])
    if area<=0 or up<=0 or not all(math.isfinite(a) for a in (area,up,wa)):
        raise RuntimeError("Original basin unit or static WWTP source impossible")
    row={
     "HYBAS_ID_L07":hid,"PERIOD":period,
     "country":mostcommon(z[2] for z in stations),
     "original_majority_BQE_method":mostcommon(z[3] for z in stations),
     "original_method_count_in_basin_year":len({z[3] for z in stations}),
     "original_country_count_in_basin_year":len({z[2] for z in stations}),
     "original_EQR_physical_site_n":len(stations),
     "EQR_original_site_equal_basin_period_mean":float(np.mean([z[0] for z in stations])),
     "plant_S50_GBIF_speciesKey":float(s["original_source_plant_S50"]),
     "lycaenidae_fivefamily_index_share_logodds":math.log((nly+0.5)/(nbut-nly+0.5)),
     "lycaenidae_fivefamily_index_share":nly/nbut,
     "W_static_upstream_MODEL_flow_per_UpArea":math.log1p(max(0,wa)/(up+1)),
     "log_local_H7_area_km2":math.log(area),
     "log1p_original_EQR_site_n":math.log1p(len(stations)),
     "log1p_GBIF_Animalia_source_records":math.log1p(float(s["original_GBIF_Animalia_record_opportunity"])),
     "log1p_GBIF_Plantae_source_records":math.log1p(float(s["original_GBIF_Plantae_record_opportunity"])),
     "log1p_GBIF_5family_butterfly_source_records":math.log1p(nbut),
     "source_valid_EQR_DIS_AV_CMS_all_stations":all(z[1] is not None and z[1]>=0 for z in stations),
     "source_mean_DIS_AV_CMS":float(np.mean([z[1] for z in stations])) if all(z[1] is not None and z[1]>=0 for z in stations) else None,
     "static_flow_is_not_actual_annual_discharge_or_receptor_specific_dose":True
    }
    rows.append(row)
basins={r["HYBAS_ID_L07"] for r in rows}
countries=collections.Counter(z["country"] for z in rows)
if len(rows)<320 or len(basins)<170 or len(countries)<12:
    raise RuntimeError("Predeclared source cross-sectional original 320 rows/170 basins/12 countries gate failed")
rowkeys={(z["HYBAS_ID_L07"],z["PERIOD"]) for z in rows}
if len(rows)!=len(rowkeys):raise RuntimeError("Duplicate original basin-period inference grain")
hasflow=all(z["source_valid_EQR_DIS_AV_CMS_all_stations"] for z in rows)
base_num=[
 "log_local_H7_area_km2",
 "log1p_original_EQR_site_n",
 "log1p_GBIF_Animalia_source_records",
 "log1p_GBIF_Plantae_source_records",
 "log1p_GBIF_5family_butterfly_source_records"]
if hasflow:
    for z in rows:z["log1p_original_EEA_site_DIS_AV_CMS"]=math.log1p(z["source_mean_DIS_AV_CMS"])
    base_num+=["log1p_original_EEA_site_DIS_AV_CMS"]
country_list=sorted(countries)
method_list=sorted(set(z["original_majority_BQE_method"] for z in rows))
# Cross-sectional source is largely same method per country, so method and
# country FE may be rank-deficient. Ridge stable; OLS pseudoinverse with diagnostics.
N=len(rows)
fields=base_num
def Xmake(addWWTP=False):
    xx=np.column_stack([
     np.ones(N),
     np.array([[float(z[k]) for k in fields] for z in rows],dtype=float),
     np.array([[float(z["country"]==c) for c in country_list[1:]] for z in rows],dtype=float),
     np.array([[float(z["original_majority_BQE_method"]==m) for m in method_list[1:]] for z in rows],dtype=float),
     np.array([float(z["PERIOD"]=="2020-2022") for z in rows],dtype=float)
    ])
    if addWWTP:xx=np.column_stack([xx,np.array([z["W_static_upstream_MODEL_flow_per_UpArea"] for z in rows])])
    return xx
X0=Xmake(False);X1=Xmake(True)
group=[z["HYBAS_ID_L07"] for z in rows]
folds=np.array([int(hashlib.sha256(str(v).encode()).hexdigest()[:8],16)%5 for v in group])
if min(np.bincount(folds,minlength=5))<35:raise RuntimeError("Original H7 fold support too weak")
group_id,group_inv=np.unique(group,return_inverse=True)
G=len(group_id)
rng=np.random.default_rng(20261112)
idx=rng.integers(0,G,size=(2000,G))
def predict(X,Y):
    preds=np.empty(N)
    for f in range(5):
        train=folds!=f;test=~train
        Z=X.copy()
        mu=Z[train,1:].mean(axis=0);sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        A=Z[train];B=Z[test]
        reg=np.eye(A.shape[1])*10.;reg[0,0]=0.
        beta=np.linalg.solve(A.T@A+reg,A.T@Y[train])
        preds[test]=B@beta
    return preds
def clustered_beta(Y,weights=None):
    Z=X1.copy()
    mu=Z[:,1:].mean(axis=0);sd=Z[:,1:].std(axis=0)
    sd=np.where(sd>1e-12,sd,1.)
    Z[:,1:]=(Z[:,1:]-mu)/sd
    if weights is None:weights=np.ones(N)
    A=Z*np.sqrt(weights)[:,None];b=Y*np.sqrt(weights)
    coef=np.linalg.lstsq(A,b,rcond=None)[0]
    resid=b-A@coef
    bread=np.linalg.pinv(A.T@A)
    meat=np.zeros((A.shape[1],A.shape[1]))
    for ix in range(G):
        groupx=group_inv==ix
        v=A[groupx].T@resid[groupx]
        meat+=np.outer(v,v)
    rank=np.linalg.matrix_rank(A)
    variance=bread@meat@bread*(G/(G-1)*(N-1)/(N-rank))
    se=np.sqrt(np.maximum(variance.diagonal(),0.))
    j=A.shape[1]-1
    beta=float(coef[j]); sdY=float(np.std(Y))
    return {
      "beta_EQR_or_S50_or_logit_response_per_1SD_static_WWTP":beta,
      "H7_basin_cluster_descriptive95":[float(beta-1.96*se[j]),float(beta+1.96*se[j])],
      "beta_outcome_SD_per_1SD_WWTP":float(beta/sdY) if sdY else None,
      "rank":rank,"source_columns":A.shape[1],"original_WWTP_log_gradient_SD":float(sd[-1])
    }
weights=np.array([1/countries[z["country"]] for z in rows],dtype=float)
weights/=weights.mean()
results={}
for name,key in [
 ("actual_EEA_EQR_original_source","EQR_original_site_equal_basin_period_mean"),
 ("GBIF_Plantae_speciesKey_S50_50rarefaction","plant_S50_GBIF_speciesKey"),
 ("GBIF_Lycaenidae_log_odds_5_butterfly_families","lycaenidae_fivefamily_index_share_logodds")
]:
    Y=np.array([float(z[key]) for z in rows])
    P0=predict(X0,Y);P1=predict(X1,Y)
    e0=(Y-P0)**2;e1=(Y-P1)**2
    mse0=float(e0.mean());mse1=float(e1.mean())
    grouped0=np.bincount(group_inv,weights=e0,minlength=G)
    grouped1=np.bincount(group_inv,weights=e1,minlength=G)
    ba=np.sum(grouped0[idx],axis=1)
    bb=np.sum(grouped1[idx],axis=1)
    boots=100*(ba-bb)/ba
    raw=spearmanr([z["W_static_upstream_MODEL_flow_per_UpArea"] for z in rows],Y)
    results[name]={
      "n_original_basinperiods":N,"n_independent_H7_basins":G,
      "mean_response":float(Y.mean()),
      "response_SD":float(Y.std()),
      "raw_spearman_static_WWTP_gradient":float(raw.statistic),
      "raw_spearman_p_descriptive_only":float(raw.pvalue),
      "baseline_5fold_hash_basin_disjoint_MSE":mse0,
      "plus_static_modelled_WWTP_gradient_5fold_MSE":mse1,
      "out_of_basin_5fold_MSE_gain_percent":100*(mse0-mse1)/mse0,
      "H7_cluster_2000_bootstrap_gain_descriptive95":np.percentile(boots,[2.5,97.5]).tolist(),
      "fraction_bootstrap_improves_prediction":float(np.mean(boots>0)),
      "all_source_country_and_method_FE_adjusted_H7_cluster_beta":clustered_beta(Y),
      "country_equal_total_weight_H7_cluster_beta":clustered_beta(Y,weights),
    }
# Source-saturation checks; wastewater status unknown and modeled discharge
# are NOT environmental intervention.
out={
 "schema":"GESUI_HYBAS7_234EU_source_matched_three_response_static_WWTP_gradient_posthoc_v1",
 "not_causal":True,"not_realized_ecosystem_services":True,
 "original_source_candidates_SHA256":SOURCE_SHA,
 "original_EEA_station_to_HYBAS7_source_SHA256":H7_SHA,
 "original_EEA_GitBlob_source_ids":{k:v[1] for k,v in BLOBS.items()},
 "source_grain":{"basinperiod_rows":N,"basins":G,"countries":dict(countries),
   "majority_BQE_methods":dict(collections.Counter(z["original_majority_BQE_method"] for z in rows)),
   "country_ambiguous_periods":sum(z["original_country_count_in_basin_year"]>1 for z in rows),
   "methods_mixed_periods":sum(z["original_method_count_in_basin_year"]>1 for z in rows),
   "original_record_QA_excluded_reasons":dict(excluded),
   "periods":dict(collections.Counter(z["PERIOD"] for z in rows)),
   "five_fold_basins":{str(k):len({group[i] for i in range(N) if folds[i]==k}) for k in range(5)},
   "year_to_year_basin_repeated_count":sum(len({z["PERIOD"] for z in rows if z["HYBAS_ID_L07"]==hid})==2 for hid in basins),
   "river_flow_covariate_complete_all_rows":hasflow,
   "modelled_wwtp_outlet_pressure_present_locally_fraction":sum(z["W_static_upstream_MODEL_flow_per_UpArea"]>0 for z in rows)/N},
 "fixed_common_covariates":fields+["country_FE","modal_original_BQE_method_FE","year_bin_FE"],
 "single_static_predictor":"log1p(static original HydroWASTE WASTE_DIS / (upstream catchment area +1)); no receptor-specific effluent dose",
 "results_three_response_separate_models":results,
 "limitations":[
   "Not an observed downriver WWTP dose to original EEA monitoring stations, original 58k source status mostly not reported, and no flow-standardized true concentrations.",
   "Full H7 basin static source and each of three 2018/2022 index/quality response variables share geography, but some environmental landuse and meteorology are unadjusted.",
   "EEA quality country/method adjusted but original source EQR column normalizations/year differences may be nonidentical.",
   "GBIF S50 is sampled speciesKey index richness and five-family Lycaenidae share is observation composition, neither direct true biodiversity, abundance or services.",
   "The baseline observational counts overlap numerator/denominator of some outcomes; collider and survey bias remain.",
   "Two-period repeated H7 is kept in the same hashed fold; nearby river basins could be in both train/test. Not spatially/region-independent external prediction.",
   "Comparing beta across response families means 1 SD source gradient versus 1 SD of each response, not equal original measurement units.",
   "This is POSTHOC exploration after seeing same-basin temporal co-change fail; must not claim a new preregistered confirmation."
 ]}
path=OUT/"EU234_H7_SAME_BASINPLANT_BUTTERFLY_EQR_STATIC_HYDROWASTE_GRADIENT_RESULTS.json"
def normalize_numpy_number(v):
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,np.floating):return float(v)
    if isinstance(v,np.ndarray):return v.tolist()
    raise TypeError("Non-serializable original source result type "+str(type(v)))
path.write_text(json.dumps(out,indent=2,sort_keys=True,default=normalize_numpy_number,allow_nan=False))
with (OUT/"EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv").open("w",newline="") as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]))
    writer.writeheader();writer.writerows(rows)
print("HYBAS7_234EU_THREE_LAYER_SPATIAL_WWTP_RESULT_START",flush=True)
print(json.dumps({"source_grain":out["source_grain"],"results":results},sort_keys=True,default=normalize_numpy_number,allow_nan=False),flush=True)
print("HYBAS7_234EU_THREE_LAYER_SPATIAL_WWTP_RESULT_END",flush=True)
