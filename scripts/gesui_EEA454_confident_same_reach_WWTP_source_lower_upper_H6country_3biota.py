#!/usr/bin/env python3
"""Same original EU232 H7 /454 periods, 90 true H6 parent spatial blocked
and stable-country-out CV. Physical reach-specific upstream HydroWASTE
WASTE_DIS compared to prior H7-basin outlet source and 50km downstream control.

Original native 938544 reach HydroRIVERS river network already processed in
an independent source-only SUCCESS before this ecological comparison. Only
derived receptor×period source data, no standalone licensed river network.
POSTHOC exploratory, no nutrient/delivery or realized service causal inference.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HYRIV_EU454_SAME_REACH_CONFIDENT_MACRO"
ROOT.mkdir(parents=True,exist_ok=True)
SHA={
 "receiver":"d50eba0f934278cb28eb403af7b84d4da5a19af209a6d05be312a2ee44ad6d6d",
 "H6":"372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc",
 "hydroatlas":"ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686",
 "annualNP":"012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"
}
def read(env,sha=None):
    p=pathlib.Path(os.environ[env]);raw=p.read_bytes()
    if sha and hashlib.sha256(raw).hexdigest()!=sha:
        raise RuntimeError("Previously frozen original source identity changed "+env)
    return json.loads(raw)
receiver=read("TRUE_RECEIVER_WWTP_SOURCE",SHA["receiver"])
if receiver["schema"]!="GESUI_EU454_ORIGINAL_HYDRO_RIVERS_STATION_RECEPTOR_WWTP_SOURCE_GATE_v1":
    raise RuntimeError("Original 454×period receptor source schema differs")
Xrec=receiver["H7_year_original_derived_receiver_site_mean_exposures"]
if len(Xrec)!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in Xrec})!=454:
    raise RuntimeError("Original HydroRIVERS derived receptor source has duplicate/missing H7×period keys")
rmap={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in Xrec}
h6=read("TRUE_SOURCE_H6_PARENT_MAPPING",SHA["H6"])
if h6["schema"]!="GESUI_232_EU_H7_TO_AUTHENTIC_HYBAS6_EXPLICIT_INVALID_GEOMETRY_REPAIR_SOURCE_SENSITIVITY_v2":
    raise RuntimeError("Previously explicit source-repaired parent-only geographic unit has changed")
parent={r["HYBAS_ID_L07"]:r["official_parent_HYBAS_ID_L06"] for r in h6["original_454_H7_study_basin_geographic_parent_mapping"]}
if len(parent)!=232 or len(set(parent.values()))!=90 or None in parent.values():
    raise RuntimeError("Original 232→90 H6 parent source mapping changed")
atlas=read("TRUE_SOURCE_H7_HYDROATLAS",SHA["hydroatlas"])["original_1219_H7_covariate_by_real_ID"]
annual=read("TRUE_SOURCE_EEA_REPORTED_NP",SHA["annualNP"])
np_source={(r["HYBAS_ID_L07"],r["PERIOD"]):r for r in annual["source_H7_period_reported_N_P_upstream_to_true_original_station"]}
prior=read("TRUE_SOURCE_PREVIOUS_H6_MODEL_RESULTS")
if prior["schema"]!="GESUI_454_EU_H7_90_H6_PARENT_15_COUNTRY_SPATIAL_TRANSFER_POSTHOC_v1":
    raise RuntimeError("Original 90 H6 parent or country model results changed")

with open(os.environ["TRUE_SOURCE_454_EEA_GBIF_ORIGINAL_ROWS"],encoding="utf-8-sig",newline="") as f:
    rawrows=list(csv.DictReader(f))
N=len(rawrows)
if N!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in rawrows})!=454:
    raise RuntimeError("Original 232 H7×454 period study cohort changed")
if len({r["HYBAS_ID_L07"] for r in rawrows})!=232:
    raise RuntimeError("Original 232 H7 study group identity changed")
if set((r["HYBAS_ID_L07"],r["PERIOD"]) for r in rawrows)!=set(rmap):
    raise RuntimeError("Original EEA receptor source not one-to-one match with real 454 study rows")

# Fixed original, source-only same-reach position-quality study carried out
# before viewing any NEW ecological result. Origin 471 original station-year
# records with 700 registered plant×site pairs; 352 confident UPSTREAM,
# 222 confident DOWNSTREAM, 126 source-orientation/position UNRESOLVED.
SOURCE_SAME_REACH_SHA="9bf89e7c372cfd294113f913955f48dc1b6262d9435b8918b208758e9005d302"
same_source=read("ORIGINAL_SAME_HYRIV_WWTP_POSITION_SOURCE",SOURCE_SAME_REACH_SHA)
if same_source["schema"]!="GESUI_EEA454_SOURCE_ONLY_TRUE_HYRIV_SAME_REACH_WWTP_POINT_ORDER_QA_20261113_v1":
    raise RuntimeError("Original SAME reach source-only quality gate is not the authentic study")
same= same_source["same_reach_EEA_station_year_source_decisions"]
if len(same)!=471 or sum(r["source_confident_upstream_plant_n"] for r in same)!=352 or \
   sum(r["source_confident_downstream_plant_n"] for r in same)!=222 or \
   sum(r["source_uncertain_same_reach_plant_n"] for r in same)!=126:
    raise RuntimeError("Original pre-biological same-reach WWTP source 471/352/222/126 changed")
def sourceblob(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()
original_eu_year={
 "2017-2019":("results/eea2018_independent_temporal_bqe_20261108/benthic_invertebrates_2018_WWTP_source_real_EQR.csv","0f30bfbc3dac13f6a182c8a910e9899c04778322"),
 "2020-2022":("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv","325f199e985d9862ffee5b0edb205c61a7be3833")
}
riverflows={}
for period,(f,digest) in original_eu_year.items():
    b=pathlib.Path(f).read_bytes()
    if sourceblob(b)!=digest:raise RuntimeError("Original published river source file hash changed")
    rr=list(csv.DictReader(b.decode("utf-8-sig").splitlines()))
    if len(rr)!={"2017-2019":2634,"2020-2022":2475}[period]:
        raise RuntimeError("Original EEA source station year count changed")
    for r in rr:
        x=float(r["DIS_AV_CMS"])
        if not math.isfinite(x) or x<=0:raise RuntimeError("Original EEA receiver flow unsupported")
        key=(period,r["site_key"])
        if key in riverflows:raise RuntimeError("Duplicate same-year EEA monitoring station flow")
        riverflows[key]=x
same_keys=set()
additional=collections.defaultdict(lambda:{"u":0.,"d":0.,"a":0.,"up_n":0,"down_n":0,"unknown_n":0,"source_sites":0})
for item in same:
    key=(item["HYBAS_ID_L07"],item["PERIOD"],item["original_EEA_station"])
    if key in same_keys:raise RuntimeError("Same EEA source station×H7×year duplicated in 471 geometry rows")
    same_keys.add(key)
    kp=(item["HYBAS_ID_L07"],item["PERIOD"])
    if kp not in rmap:raise RuntimeError("Original geometric source site not in fixed 454 basin-year cohort")
    f=riverflows[(item["PERIOD"],item["original_EEA_station"])]
    source=additional[kp]
    source["source_sites"]+=1
    source["u"]+=float(item["source_confident_upstream_WASTE_DIS_m3_day"])/(f*86400)
    source["d"]+=float(item["source_confident_downstream_WASTE_DIS_m3_day"])/(f*86400)
    source["a"]+=float(item["source_unresolved_same_reach_WASTE_DIS_m3_day"])/(f*86400)
    source["up_n"]+=int(item["source_confident_upstream_plant_n"])
    source["down_n"]+=int(item["source_confident_downstream_plant_n"])
    source["unknown_n"]+=int(item["source_uncertain_same_reach_plant_n"])
if len(same_keys)!=471 or sum(x["up_n"] for x in additional.values())!=352:
    raise RuntimeError("Source integrated same-reach original registered flow conservation failed")
for k in rmap:
    n=int(rmap[k]["original_usable_network_receiver_stations"])
    if n<1:raise RuntimeError("Previously source-validated original 454 receiver sites missing")
    vals=additional[k]
    for t in ("u","d","a"):vals[t]/=n
if sum(x["source_sites"] for x in additional.values())!=471:
    raise RuntimeError("Original same-reach station-year source accounting failed")

def nonnegative(v,name):
    f=float(v)
    if not math.isfinite(f) or f<0:
        raise RuntimeError("Original source negative/nonfinite "+name)
    return f
rows=[]
for base in rawrows:
    r=dict(base)
    key=r["HYBAS_ID_L07"];period=r["PERIOD"]
    src=rmap[(key,period)]
    if not src["network_source_valid"] or src["original_usable_network_receiver_stations"]<=0:
        raise RuntimeError("Original predeclared 454 study receptor source incomplete; do not impute")
    original=atlas[key]
    r["H6_parent_native"]=parent[key]
    r["source_crop_up"]=nonnegative(original["crp_pc_use"],"crop")
    r["source_urban_up"]=nonnegative(original["urb_pc_use"],"urban")
    r["source_natural_outlet_flow_log1p"]=math.log1p(nonnegative(original["dis_m3_pyr"],"outlet_flow"))
    fixed=lambda tag:nonnegative(src["station_equal_mean_"+tag],tag)
    r["source_W_receiver_strict_all_log1p"]=math.log1p(fixed("strict_upstream_WWTP_fraction_of_original_model_river_discharge"))
    r["source_W_receiver_strict50_log1p"]=math.log1p(fixed("strict_upstream_50km_WWTP_fraction_of_original_model_river_discharge"))
    r["source_W_receiver_DOWNSTREAM_50_negativecontrol_log1p"]=math.log1p(fixed("downstream_50km_WWTP_fraction_direction_control"))
    r["source_RECEPTOR_same_reach_UNRESOLVED_count"]=fixed("source_original_same_reach_facilities_unresolved_COUNT")

    # Original source strict-upstream numerator is an original site-equal MEAN
    # of per-site modelled WWTP m3/day / (original receiver cms ×86400).
    # The new confirmed SAME-reach source is the SUM of exact per-site
    # fractions divided by the SAME original number of valid site-years.
    # Unresolved same-reach source is a possible upper bound, NOT a confirmed
    # upstream or true absent effluent.
    extra=additional[(key,period)]
    strict_frac=fixed("strict_upstream_WWTP_fraction_of_original_model_river_discharge")
    downstream_frac=fixed("downstream_50km_WWTP_fraction_direction_control")
    r["source_W_receiver_CONFIDENT_same_lower_log1p"]=math.log1p(strict_frac+extra["u"])
    r["source_W_receiver_POSSIBLE_same_upper_log1p"]=math.log1p(strict_frac+extra["u"]+extra["a"])
    r["source_W_receiver_confident_same_downstream_log1p"]=math.log1p(downstream_frac+extra["d"])
    r["source_CONFIDENT_SAME_REACH_upstream_modelled_flow_perstation_ratio"]=extra["u"]
    r["source_CONFIDENT_SAME_REACH_downstream_modelled_flow_perstation_ratio"]=extra["d"]
    r["source_UNRESOLVED_SAME_REACH_worst_case_upstream_perstation_ratio"]=extra["a"]
    r["source_original_same_reach_up_source_plant_x_station_pairs"]=extra["up_n"]
    r["source_original_same_reach_remaining_ambiguous_source_pairs"]=extra["unknown_n"]

    # Status/quality original m3/day at site-mean, normalize by the *site-mean*
    # EEA modeled river DIS_AV_CMS, distinct from exact mean of per-site fractions.
    orig_flow=math.expm1(float(r["log1p_original_EEA_site_DIS_AV_CMS"]))
    if orig_flow<=0:raise RuntimeError("Original 454 row source modeled site flow nonpositive")
    r["source_W_operational50_log1p_mean_flow_ratio"]=math.log1p(
        fixed("upstream_50km_original_STATUS_operational_model_flow_m3_day")/(orig_flow*86400))
    r["source_W_QUAL1_50_log1p_mean_flow_ratio"]=math.log1p(
        fixed("upstream_50km_original_QUAL_WASTE_1_model_flow_m3_day")/(orig_flow*86400))
    r["source_50km_true_registered_positive"]=bool(fixed("upstream_50km_distinct_reach_registered_WWTP_count")>0)
    r["source_all_true_registered_positive"]=bool(fixed("strict_upstream_all_distinct_reach_source_registered_WWTP_count")>0)
    annualrec=np_source[(key,period)]
    if not annualrec["source_valid"]:raise RuntimeError("Original pre-source annual WWTP N sample invalid")
    r["source_reported_N_log1p_full_delivery_equiv"]=math.log1p(
       nonnegative(annualrec["source_year_reported_N_mgL_equivalent_site_equal_mean"],"reported_N"))
    r["source_annual_EEA_NP_reported_positive"]=int(annualrec["source_station_annual_upstream_NP_positive_n"]>0)
    rows.append(r)
assert len(rows)==454
source_positive=sum(z["source_all_true_registered_positive"] for z in rows)
if source_positive!=301 or sum(z["source_50km_true_registered_positive"] for z in rows)!=300:
    raise RuntimeError("Original 301/300 receptor-source positive basin-period count changed")
if sum(z["source_annual_EEA_NP_reported_positive"] for z in rows)!=56:
    raise RuntimeError("Original reported annual EEA source-positive n56 changed")

# Independent H6 parent holdout and stable country holdout, original schema.
CATS=sorted({r["country"] for r in rows})
METHODS=sorted({r["original_majority_BQE_method"] for r in rows})
base=["log_local_H7_area_km2","log1p_original_EQR_site_n",
      "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
      "log1p_GBIF_5family_butterfly_source_records","log1p_original_EEA_site_DIS_AV_CMS"]
def xmat(with_country=True):
    return np.column_stack([
     np.ones(N),
     np.array([[float(r[k]) for k in base] for r in rows]),
     np.array([[float(r["country"]==c) for c in CATS[1:]] for r in rows]) if with_country else np.zeros((N,0)),
     np.array([[float(r["original_majority_BQE_method"]==m) for m in METHODS[1:]] for r in rows]),
     np.array([float(r["PERIOD"]=="2020-2022") for r in rows]),
     np.array([[r["source_crop_up"],r["source_urban_up"],r["source_natural_outlet_flow_log1p"]] for r in rows])
    ])
def models(country=True):
    B=xmat(country)
    def col(k):return np.array([float(r[k]) for r in rows])
    W=col("W_static_upstream_MODEL_flow_per_UpArea")
    U=col("source_W_receiver_strict_all_log1p")
    U50=col("source_W_receiver_strict50_log1p")
    DN=col("source_W_receiver_DOWNSTREAM_50_negativecontrol_log1p")
    NN=col("source_reported_N_log1p_full_delivery_equiv")
    OP=col("source_W_operational50_log1p_mean_flow_ratio")
    Q1=col("source_W_QUAL1_50_log1p_mean_flow_ratio")
    UNC=col("source_W_receiver_CONFIDENT_same_lower_log1p")
    UPPER=col("source_W_receiver_POSSIBLE_same_upper_log1p")
    DADD=col("source_W_receiver_confident_same_downstream_log1p")
    return {
     "B2_environment_only":B,
     "B3_original_H7_basin_OUTLET_W":np.column_stack([B,W]),
     "B4_original_H7_basin_OUTLET_W_plus_EEA_annual_N":np.column_stack([B,W,NN]),
     "B2_plus_receiver_STRICT_all_upstream":np.column_stack([B,U]),
     "B2_plus_receiver_STRICT_50km":np.column_stack([B,U50]),
     "B2_plus_DOWNSTREAM_50km_direction_control":np.column_stack([B,DN]),
     "B2_plus_receiver_all_AND_old_basin_outlet":np.column_stack([B,U,W]),
     "B2_plus_receiver_all_AND_downstream_control":np.column_stack([B,U,DN]),
     "B2_plus_receiver_all_AND_annual_reported_N":np.column_stack([B,U,NN]),
     "B2_plus_receiver50_AND_operational_only50":np.column_stack([B,U50,OP]),
     "B2_plus_receiver50_AND_QUAL_WASTE1_reported_only50":np.column_stack([B,U50,Q1]),
     "B2_plus_old_basin_W_AND_DOWNSTREAM_control":np.column_stack([B,W,DN]),
     "B2_plus_CONFIDENT_SAME_reach_upstream_lower":np.column_stack([B,UNC]),
     "B2_plus_POSSIBLE_SAME_reach_upstream_upper":np.column_stack([B,UPPER]),
     "B2_plus_CONFIDENT_SAME_reach_lower_AND_outlet":np.column_stack([B,UNC,W]),
     "B2_plus_CONFIDENT_SAME_reach_lower_AND_downstream":np.column_stack([B,UNC,DN]),
     "B2_plus_DOWNSTREAM_confident_same_source":np.column_stack([B,DADD])
    }
A=models(True)
Ac=models(False)
same_h6=np.array([r["H6_parent_native"] for r in rows])
fgeo=np.array([int(hashlib.sha256(id.encode()).hexdigest()[:8],16)%5 for id in same_h6])
H6,ih6=np.unique(same_h6,return_inverse=True)
G=len(H6)
if G!=90 or min(np.bincount(fgeo,minlength=5))<15:
    raise RuntimeError("Original 90 native H6 hydrologic source geographic group fold broken")
def grouped_country():
    cnt=collections.defaultdict(collections.Counter)
    for r in rows:cnt[r["HYBAS_ID_L07"]][r["country"]]+=1
    group={h:sorted(v,key=lambda x:(-v[x],x))[0] for h,v in cnt.items()}
    return np.array([group[r["HYBAS_ID_L07"]] for r in rows])
country_group=grouped_country()
if len(set(country_group))!=11:
    raise RuntimeError("Original stable country 11 source blocks changed")
def predict(X,Y,groups):
    out=np.empty(N)
    for group in sorted(set(groups)):
        test=groups==group;train=~test
        if not np.any(test) or not np.any(train):raise RuntimeError("Empty source group fold")
        Z=X.copy()
        mu=Z[train,1:].mean(axis=0);sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1.0)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        a=Z[train];b=Z[test]
        reg=np.eye(a.shape[1])*10.;reg[0,0]=0.
        beta=np.linalg.solve(a.T@a+reg,a.T@Y[train])
        out[test]=b@beta
    if not np.all(np.isfinite(out)):raise RuntimeError("Nonfinite original source routed predictor response")
    return out
CMP=[
 ("B2_environment_only","B3_original_H7_basin_OUTLET_W"),
 ("B2_environment_only","B2_plus_receiver_STRICT_all_upstream"),
 ("B2_environment_only","B2_plus_receiver_STRICT_50km"),
 ("B2_environment_only","B2_plus_DOWNSTREAM_50km_direction_control"),
 ("B3_original_H7_basin_OUTLET_W","B2_plus_receiver_all_AND_old_basin_outlet"),
 ("B2_plus_receiver_STRICT_all_upstream","B2_plus_receiver_all_AND_downstream_control"),
 ("B2_plus_receiver_STRICT_all_upstream","B2_plus_receiver_all_AND_annual_reported_N"),
 ("B2_plus_receiver_STRICT_50km","B2_plus_receiver50_AND_operational_only50"),
 ("B2_plus_receiver_STRICT_50km","B2_plus_receiver50_AND_QUAL_WASTE1_reported_only50"),
 ("B3_original_H7_basin_OUTLET_W","B2_plus_old_basin_W_AND_DOWNSTREAM_control"),
 ("B3_original_H7_basin_OUTLET_W","B4_original_H7_basin_OUTLET_W_plus_EEA_annual_N"),
 ("B2_environment_only","B2_plus_CONFIDENT_SAME_reach_upstream_lower"),
 ("B2_environment_only","B2_plus_POSSIBLE_SAME_reach_upstream_upper"),
 ("B2_plus_receiver_STRICT_all_upstream","B2_plus_CONFIDENT_SAME_reach_upstream_lower"),
 ("B2_plus_CONFIDENT_SAME_reach_upstream_lower","B2_plus_POSSIBLE_SAME_reach_upstream_upper"),
 ("B3_original_H7_basin_OUTLET_W","B2_plus_CONFIDENT_SAME_reach_lower_AND_outlet"),
 ("B2_plus_CONFIDENT_SAME_reach_upstream_lower","B2_plus_CONFIDENT_SAME_reach_lower_AND_downstream"),
 ("B2_environment_only","B2_plus_DOWNSTREAM_confident_same_source")
]
rng=np.random.default_rng(20261113)
bootstrapH6=rng.integers(0,G,(2000,G))
def compare(e,weights,geo=False):
    out={}
    for x,y in CMP:
        mse0=float(np.mean(e[x]));mse1=float(np.mean(e[y]))
        item={"signed_heldout_MSE_gain_pct":100*(mse0-mse1)/mse0,
              "baseline_MSE":mse0,"plus_MSE":mse1}
        if geo:
            e0=np.bincount(ih6,weights=e[x],minlength=G)
            e1=np.bincount(ih6,weights=e[y],minlength=G)
            d0=np.sum(e0[bootstrapH6],axis=1);d1=np.sum(e1[bootstrapH6],axis=1)
            gains=100*(d0-d1)/d0
            item["H6_parent_bootstrap_2000_95_descriptive"]=np.percentile(gains,[2.5,97.5]).tolist()
            item["H6_parent_bootstrap_positive_fraction"]=float(np.mean(gains>0))
        out[x+"_TO_"+y]=item
    return out
responses={
 "original_EEA_benthic_EQR":("EQR_original_site_equal_basin_period_mean","EEA_benthic_EQR"),
 "GBIF_index_plant_S50":("plant_S50_GBIF_speciesKey","GBIF_plant_speciesKey_S50"),
 "GBIF_index_Lycaenidae_five_family_logodds":("lycaenidae_fivefamily_index_share_logodds","GBIF_Lycaenidae_five_family_logodds")
}
answers={}
for name,(field,oldkey) in responses.items():
    Y=np.array([float(r[field]) for r in rows])
    geo_predictions={model:predict(X,Y,fgeo) for model,X in A.items()}
    geo_err={model:(Y-p)**2 for model,p in geo_predictions.items()}
    # Exact baseline reproduction prior H6 original models before interpreting
    # new source-directed predictors.
    prev=prior["models"][oldkey]["NEW_H6_PARENT_BLOCKED_5fold_original_outcome_MSE"]
    for new,old in [
        ("B2_environment_only","B2_upstream_crop_urban_naturalflow"),
        ("B3_original_H7_basin_OUTLET_W","B3_upstream_plus_modelled_WWTP"),
        ("B4_original_H7_basin_OUTLET_W_plus_EEA_annual_N","B4_upstream_plus_WWTP_and_reported_N")
    ]:
        if not math.isclose(float(np.mean(geo_err[new])),float(prev[old]),rel_tol=1e-10,abs_tol=1e-12):
            raise RuntimeError("Old original source H6-block CV baseline not reproduced "+name+":"+new)
    country_predictions={model:predict(X,Y,country_group) for model,X in Ac.items()}
    country_err={model:(Y-p)**2 for model,p in country_predictions.items()}
    prev_country=prior["models"][oldkey]["NEW_country_LEFT_OUT_non_country_FE_MSE"]
    for new,old in [
        ("B2_environment_only","B2_upstream_crop_urban_naturalflow"),
        ("B3_original_H7_basin_OUTLET_W","B3_upstream_plus_modelled_WWTP"),
        ("B4_original_H7_basin_OUTLET_W_plus_EEA_annual_N","B4_upstream_plus_WWTP_and_reported_N")
    ]:
        if not math.isclose(float(np.mean(country_err[new])),float(prev_country[old]),rel_tol=1e-10,abs_tol=1e-12):
            raise RuntimeError("Old original source country-out CV baseline not reproduced "+name+":"+new)
    answers[name]={
      "old_source_H6_and_country_9_CV_baseline_values_exactly_reproduced":True,
      "H6_true_parent90_geo_block_source_receiver_new_predictive_comparisons":compare(geo_err,bootstrapH6,True),
      "H6_parent90_heldout_all_models_MSE":{k:float(v.mean()) for k,v in geo_err.items()},
      "country_leftout_new_receiver_predictive_comparisons":compare(country_err,None,False),
      "country_leftout_all_models_MSE":{k:float(v.mean()) for k,v in country_err.items()}
    }
report={
 "schema":"GESUI_454_H7_TRUE_HYRIV_CONFIDENT_SAME_REACH_WWTP_ADDITION_LOWER_UPPER_3BIOTA_POSTHOC_20261113_v1",
 "source_original_454_receptor_WWTP_geonetwork_SHA256":SHA["receiver"],
 "source_original_232_H7_to90_H6_explicit_repaired_source_SHA256":SHA["H6"],
 "source_original_1219_HydroATLAS_environment_SHA256":SHA["hydroatlas"],
 "source_original_reported_annual_2018_2022_upstream_NP_SHA256":SHA["annualNP"],
 "original_source_454_H7period_rows":N,
 "original_source_232_H7_distinct_basins":len({r["HYBAS_ID_L07"] for r in rows}),
 "original_source_90_H6_true_geo_hydrological_parent_groups":G,
 "source_11_country_stable_H7_blocks":len(set(country_group)),
 "source_2018_22_annual_reported_N_positive_H7_periods":sum(r["source_annual_EEA_NP_reported_positive"] for r in rows),
 "source_receptor_true_strict_upstream_original_454_period_positive":source_positive,
 "source_receptor_true_strict50_original_454_period_positive":sum(r["source_50km_true_registered_positive"] for r in rows),
 "source_original_same_reach_position_order_471_stationyear_SHA256":SOURCE_SAME_REACH_SHA,
 "source_original_SAME_reach_352_CONFIDENT_upstream_source_pairs":sum(r["source_original_same_reach_up_source_plant_x_station_pairs"] for r in rows),
 "source_original_SAME_reach_126_UNRESOLVED_pairs":sum(r["source_original_same_reach_remaining_ambiguous_source_pairs"] for r in rows),
 "source_original_SAME_reach_454_periods_any_CONFIDENT_upstream_new_flow":sum(r["source_CONFIDENT_SAME_REACH_upstream_modelled_flow_perstation_ratio"]>0 for r in rows),
 "source_original_STATION_receptor_upper_bound_not_true_source_presence":True,
 "source_original_EEA_same_reach_ambiguity_explicitly_EXCLUDED":False,
 "source_downstream_control_is_NOT_randomized_equal_reporting_placebo":True,
 "source_50km_conservative_all_original_source_reach_lengths_included":True,
 "original_H6_GEO_SOURCE_NO_REPAIR_GATE_failed_before_validity_corrected":True,
 "source_original_HydroWASTE_outfall_operation_and_flow_estimated_not_direct_effluent_dose":True,
 "source_new_receiver_effects_POSTHOC_not_prospectively_validated":True,
 "original_source_crosslayer_results":answers,
 "limitations":[
   "HydroRIVERS native graph identifies only source upstream of a monitored RIVER REACH, excluding same-reach ambiguous source/receiver order, not a measured source plume at the site's point.",
   "HydroWASTE original registered WWTP site statuses mostly not reported, and WASTE_DIS is static mixed-quality m3/day, not measured annual pollutant mass.",
   "The downstream 50 km control is much more frequently source-positive because downstream river/settlement geography differs and it is NOT a matched causal placebo.",
   "All original 454 rows reuse published 2018/2022 EEA EQR and GBIF indexed flora/butterfly records. Response-specific observation problems remain.",
   "Country-block validation lacks country fixed effects and cannot extrapolate unknown country-specific BQE methods or wastewater source reporting.",
   "Original HYBAS6 source invalid geometry explicitly repaired for 16 target H7 sites after strict initial source gate failed.",
   "Basin/H6 crossvalidation is source-driven posthoc; this analysis has no unseen external validation continent or decades.",
   "No measured river actual TN/TP at original EEA station, adult emergence, predator diet, flowering plant host, butterfly visit or realized pollination from this source."
 ]}
outpath=ROOT/"EU454_HYRIV_CONFIDENT_SAME_REACH_UPSTREAM_LOWER_UPPER_3LAYERS_SOURCE_RESULTS.json"
outpath.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False))
brief={k:v for k,v in report.items() if k not in ("limitations",)}
for val in brief["original_source_crosslayer_results"].values():
    val.pop("H6_parent90_heldout_all_models_MSE",None)
    val.pop("country_leftout_all_models_MSE",None)
print("HYRIV_SAME_REACH_CONFIDENT_MACRO_RESULT_START")
print(json.dumps(brief,sort_keys=True,allow_nan=False),flush=True)
print("HYRIV_SAME_REACH_CONFIDENT_MACRO_RESULT_END")
