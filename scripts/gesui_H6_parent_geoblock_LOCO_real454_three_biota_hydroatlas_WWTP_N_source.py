#!/usr/bin/env python3
"""Source-frozen 232 HydroBASINS-L7 basins × 454 genuine native period rows:
spatial GROUP CV by authentic repaired-only-when-invalid level6 parent, and
country-heldout transfer. Same original EEA benthic EQR, GBIF indexed plant
S50, Lycaenidae five family composition, HydroATLAS and HydroWASTE models.
No causal sewage, no ecosystem service observations. All results POSTHOC.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H6_GEO_COUNTRY_TRANSFER_RESULT"
ROOT.mkdir(parents=True,exist_ok=True)
SHA={
 "H6":"372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc",
 "HydroATLAS":"ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686",
 "Reported_NP":"012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"}
def rawjson(env,sha=None):
    file=pathlib.Path(os.environ[env]);raw=file.read_bytes()
    if sha and hashlib.sha256(raw).hexdigest()!=sha:
        raise RuntimeError("Previously source-frozen file SHA256 differs: "+env)
    return json.loads(raw)
geo=rawjson("ORIGINAL_TRUE_H6_GEOGRAPHY",SHA["H6"])
if geo["schema"]!="GESUI_232_EU_H7_TO_AUTHENTIC_HYBAS6_EXPLICIT_INVALID_GEOMETRY_REPAIR_SOURCE_SENSITIVITY_v2":
    raise RuntimeError("Original source geo v2 repaired-only-invalid not validated")
maps=geo["original_454_H7_study_basin_geographic_parent_mapping"]
group={x["HYBAS_ID_L07"]:x["official_parent_HYBAS_ID_L06"] for x in maps}
if len(group)!=232 or len(set(group.values()))!=90 or None in group.values():
    raise RuntimeError("Original 232 H7 → 90 H6 geography changed")
source_repaired={x["HYBAS_ID_L07"] for x in maps if x["parent_source_was_invalid_and_required_geometry_repair"]}
if len(source_repaired)!=16:raise RuntimeError("Original 16 repaired-parent study H7 source audit changed")
atlas=rawjson("ORIGINAL_H7_HYDROATLAS_SOURCE",SHA["HydroATLAS"])["original_1219_H7_covariate_by_real_ID"]
mass=rawjson("ORIGINAL_H7_ANNUAL_REPORTED_SOURCE",SHA["Reported_NP"])
masses={(x["HYBAS_ID_L07"],x["PERIOD"]):x for x in mass["source_H7_period_reported_N_P_upstream_to_true_original_station"]}
prior=rawjson("PREVIOUS_EXACT_H7_LAND_ADJUSTED_RESULT")
if prior["schema"]!="GESUI_EU454_H7_ORIGINAL_THREE_RESPONSE_REPORTED_NP_LANDCOVER_ADJUSTED_v1":
    raise RuntimeError("Previously saved original 454 source model type changed")
with open(os.environ["ORIGINAL_454_EEA_GBIF_SOURCE_ROWS"],encoding="utf-8-sig",newline="") as f:
    rows=list(csv.DictReader(f))
if len(rows)!=454 or len({(r["HYBAS_ID_L07"],r["PERIOD"]) for r in rows})!=454:
    raise RuntimeError("Prior original 454 source native H7 period cohort drift")
if len({r["HYBAS_ID_L07"] for r in rows})!=232 or \
    set(z["HYBAS_ID_L07"] for z in rows)!=set(group):
    raise RuntimeError("Source H7 new parent mapping not same original study 232 river basins")
periods={z["PERIOD"] for z in rows}
if periods!={"2017-2019","2020-2022"}:raise RuntimeError("Original source period drift")
for z in rows:
    h=z["HYBAS_ID_L07"]
    if h not in atlas:raise RuntimeError("Original HydroATLAS H7 missing")
    if (h,z["PERIOD"]) not in masses:raise RuntimeError("Original year-specific reported N missing")
    z["H6_parent_native"]=group[h]
    z["SOURCE_H6_PARENT_WAS_EXPLICITLY_REPAIRED"]=h in source_repaired
    u=atlas[h]
    z["crop_up"]=float(u["crp_pc_use"])
    z["urban_up"]=float(u["urb_pc_use"])
    z["naturalout"]=math.log1p(float(u["dis_m3_pyr"]))
    z["crop_local"]=float(u["crp_pc_sse"])
    z["urban_local"]=float(u["urb_pc_sse"])
    original=masses[(h,z["PERIOD"])]
    if not original["source_valid"]:raise RuntimeError("Source reported annual WWTP station flow unavailable")
    z["reportN"]=math.log1p(float(original["source_year_reported_N_mgL_equivalent_site_equal_mean"]))
    z["sourceN_positive"]=int(original["source_station_annual_upstream_NP_positive_n"]>0)
    if not (0<=z["crop_up"]<=100 and 0<=z["urban_up"]<=100 and math.isfinite(z["naturalout"])):
        raise RuntimeError("Source agricultural/urban/flow covariate bad")
if sum(z["sourceN_positive"] for z in rows)!=56:
    raise RuntimeError("Original 56 reported upstream source-positive H7 period index changed")
N=len(rows);group6=[r["H6_parent_native"] for r in rows]
folds=np.array([int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5 for x in group6])
uniqH6,invH6=np.unique(group6,return_inverse=True)
G=len(uniqH6)
assert G==90
if min(np.bincount(folds,minlength=5))<15:raise RuntimeError("Blocked H6 rows per fold insufficient")
groups_per_fold={k:len({group6[i] for i in range(N) if folds[i]==k}) for k in range(5)}
if min(groups_per_fold.values())<5:raise RuntimeError("H6 source parents per fold insufficient")
# Country holdout must keep BOTH periods of any one H7 together even if a
# border monitoring station's majority-country source code changes between years.
byH7=collections.defaultdict(collections.Counter)
for r in rows:byH7[r["HYBAS_ID_L07"]][r["country"]]+=1
H7stableCountry={h:sorted(counter,key=lambda c:(-counter[c],c))[0] for h,counter in byH7.items()}
for r in rows:r["heldout_H7_stable_country"]=H7stableCountry[r["HYBAS_ID_L07"]]
COUNTRIES=sorted({r["country"] for r in rows})
METHODS=sorted({r["original_majority_BQE_method"] for r in rows})
base_features=[
 "log_local_H7_area_km2","log1p_original_EQR_site_n",
 "log1p_GBIF_Animalia_source_records","log1p_GBIF_Plantae_source_records",
 "log1p_GBIF_5family_butterfly_source_records","log1p_original_EEA_site_DIS_AV_CMS"
]
def makeX(include_country=True):
    baseline=np.column_stack([
      np.ones(N),
      np.array([[float(r[k]) for k in base_features] for r in rows]),
      np.array([[float(r["country"]==c) for c in COUNTRIES[1:]] for r in rows]) if include_country else np.zeros((N,0)),
      np.array([[float(r["original_majority_BQE_method"]==m) for m in METHODS[1:]] for r in rows]),
      np.array([float(r["PERIOD"]=="2020-2022") for r in rows])
    ])
    env=np.array([[r["crop_up"],r["urban_up"],r["naturalout"]] for r in rows])
    envlocal=np.array([[r["crop_local"],r["urban_local"],r["naturalout"]] for r in rows])
    W=np.array([float(r["W_static_upstream_MODEL_flow_per_UpArea"]) for r in rows])
    R=np.array([r["reportN"] for r in rows])
    return {
      "B0":baseline,
      "B1":np.column_stack([baseline,W]),
      "B2_upstream_crop_urban_naturalflow":np.column_stack([baseline,env]),
      "B3_upstream_plus_modelled_WWTP":np.column_stack([baseline,env,W]),
      "B4_upstream_plus_WWTP_and_reported_N":np.column_stack([baseline,env,W,R]),
      "B2local":np.column_stack([baseline,envlocal]),
      "B3local_plus_modelled_WWTP":np.column_stack([baseline,envlocal,W]),
      "B4local_plus_WWTP_and_reported_N":np.column_stack([baseline,envlocal,W,R])
    }
X6=makeX(True);X_country=makeX(False)
source_columns={
 "B0":"B0_source_recording_country_method_flow",
 "B1":"B1_plus_STATIC_HydroWASTE",
 "B2_upstream_crop_urban_naturalflow":"B2_plus_UPSTREAM_crop_urban_natural_flow",
 "B3_upstream_plus_modelled_WWTP":"B3_UPSTREAM_env_plus_STATIC_HydroWASTE",
 "B4_upstream_plus_WWTP_and_reported_N":"B4_UPSTREAM_env_STATIC_plus_ANNUAL_reported_N",
 "B2local":"B2local_plus_LOCAL_crop_urban_natural_flow",
 "B3local_plus_modelled_WWTP":"B3local_LOCAL_env_plus_STATIC_HydroWASTE",
 "B4local_plus_WWTP_and_reported_N":"B4local_LOCAL_env_STATIC_plus_ANNUAL_N"}
OUTCOMES={
 "EEA_benthic_EQR":("EQR_original_site_equal_basin_period_mean","EEA_original_benthic_EQR"),
 "GBIF_plant_speciesKey_S50":("plant_S50_GBIF_speciesKey","GBIF_original_plant_speciesKey_S50"),
 "GBIF_Lycaenidae_five_family_logodds":("lycaenidae_fivefamily_index_share_logodds","GBIF_original_Lycaenidae_fivefamily_logodds")
}
def pred(A,Y,holdout):
    out=np.full(N,np.nan)
    for fold in sorted(set(holdout)):
        test=holdout==fold;train=~test
        if not np.any(train) or not np.any(test):raise RuntimeError("Empty evaluation/source training group")
        Z=A.copy()
        mu=Z[train,1:].mean(axis=0);sd=Z[train,1:].std(axis=0)
        sd=np.where(sd>1e-12,sd,1)
        Z[:,1:]=(Z[:,1:]-mu)/sd
        aa=Z[train];bb=Z[test]
        reg=np.eye(aa.shape[1])*10.
        reg[0,0]=0
        beta=np.linalg.solve(aa.T@aa+reg,aa.T@Y[train])
        out[test]=bb@beta
    if not np.all(np.isfinite(out)):raise RuntimeError("Some source H6 or country groups were not held out correctly")
    return out
rng=np.random.default_rng(20261113)
bootH6=rng.integers(0,G,size=(2000,G))
def summarize(err,modelnames,byH6=True):
    summary={}
    for a,b in [
      ("B2_upstream_crop_urban_naturalflow","B3_upstream_plus_modelled_WWTP"),
      ("B3_upstream_plus_modelled_WWTP","B4_upstream_plus_WWTP_and_reported_N"),
      ("B2local","B3local_plus_modelled_WWTP"),
      ("B3local_plus_modelled_WWTP","B4local_plus_WWTP_and_reported_N")
    ]:
        baseline=float(np.mean(err[a]));improved=float(np.mean(err[b]))
        v={
          "signed_source_out_of_sample_MSE_gain_pct":100*(baseline-improved)/baseline,
          "baseline_OOS_MSE":baseline,"additional_source_OOS_MSE":improved
        }
        if byH6:
            e0=np.bincount(invH6,weights=err[a],minlength=G)
            e1=np.bincount(invH6,weights=err[b],minlength=G)
            bs0=np.sum(e0[bootH6],axis=1);bs1=np.sum(e1[bootH6],axis=1)
            gain=100*(bs0-bs1)/bs0
            v["H6_parent_bootstrap2000_descriptive95_gain_pct"]=np.percentile(gain,[2.5,97.5]).tolist()
            v["H6_parent_bootstrap_positive_fraction"]=float(np.mean(gain>0))
        summary[a+"_TO_"+b]=v
    return summary
results={}
country_hold=np.array([z["heldout_H7_stable_country"] for z in rows])
uniqueCountry=sorted(set(country_hold))
countryIndex={k:i for i,k in enumerate(uniqueCountry)}
countrygroups=np.array([countryIndex[x] for x in country_hold])
country_group_membership_counts=collections.Counter(country_hold)
reassigned_country_rows=sum(x["country"]!=x["heldout_H7_stable_country"] for x in rows)
for outcome,(column,old) in OUTCOMES.items():
    Y=np.array([float(z[column]) for z in rows])
    originalMSE=prior["separate_layer_result_models"][old]["all_11_fixed_model_fivefold_MSE"]
    originalH7=np.array([int(hashlib.sha256(z["HYBAS_ID_L07"].encode()).hexdigest()[:8],16)%5 for z in rows])
    # Independently REPRODUCE source original all model H7-hash MSE before reading new geography effects.
    originalPred={label:pred(A,Y,originalH7) for label,A in X6.items()}
    originErrors={label:(Y-v)**2 for label,v in originalPred.items()}
    for label,key in source_columns.items():
        got=float(np.mean(originErrors[label]))
        want=float(originalMSE[key])
        if not math.isclose(got,want,rel_tol=1e-11,abs_tol=1e-12):
            raise RuntimeError("H7 original source numerical check failed "+outcome+":"+label)
    geoPred={label:pred(A,Y,folds) for label,A in X6.items()}
    geoError={label:(Y-v)**2 for label,v in geoPred.items()}
    countryPred={label:pred(A,Y,countrygroups) for label,A in X_country.items()}
    countryError={label:(Y-v)**2 for label,v in countryPred.items()}
    perCountry={}
    for country in uniqueCountry:
        mask=country_hold==country
        perCountry[country]={
         "original_period_rows_heldout":int(sum(mask)),
         "original_H7_basins_heldout":len({rows[i]["HYBAS_ID_L07"] for i in range(N) if mask[i]}),
         "original_reported_NP_positive_H7_period_rows":sum(rows[i]["sourceN_positive"] for i in range(N) if mask[i]),
         "outcomes_original_biodiversity":outcome,
         "B2_to_B3_static_WWTP_gain_pct":float(100*(countryError["B2_upstream_crop_urban_naturalflow"][mask].mean()-countryError["B3_upstream_plus_modelled_WWTP"][mask].mean())/countryError["B2_upstream_crop_urban_naturalflow"][mask].mean()),
         "B3_to_B4_reported_N_gain_pct":float(100*(countryError["B3_upstream_plus_modelled_WWTP"][mask].mean()-countryError["B4_upstream_plus_WWTP_and_reported_N"][mask].mean())/countryError["B3_upstream_plus_modelled_WWTP"][mask].mean())
        }
    eligible=[q for q in uniqueCountry if perCountry[q]["original_period_rows_heldout"]>=10 and perCountry[q]["original_H7_basins_heldout"]>=5]
    def eqCountryGain(a,b,countries):
        bas=np.mean([float(np.mean(countryError[a][country_hold==country])) for country in countries])
        new=np.mean([float(np.mean(countryError[b][country_hold==country])) for country in countries])
        return 100*(bas-new)/bas
    results[outcome]={
      "original_H7_hashed_5fold_every_source_model_MSE_reproduced":True,
      "original_H7_hashed_source_model_MSE":{k:float(v.mean()) for k,v in originErrors.items()},
      "NEW_H6_PARENT_BLOCKED_5fold_original_outcome_MSE":{k:float(v.mean()) for k,v in geoError.items()},
      "NEW_H6_spatial_parent_comparisons":summarize(geoError,list(X6)),
      "NEW_country_LEFT_OUT_non_country_FE_MSE":{k:float(v.mean()) for k,v in countryError.items()},
      "NEW_country_LEFT_OUT_row_weighted_comparisons":summarize(countryError,list(X_country),False),
      "NEW_country_LEFT_OUT_major_country_equal_MSE_gains":{
        "countries_included_only_with_>=10_period_and_>=5_H7":eligible,
        "B2_TO_B3_WWTP":eqCountryGain("B2_upstream_crop_urban_naturalflow","B3_upstream_plus_modelled_WWTP",eligible),
        "B3_TO_B4_annual_N":eqCountryGain("B3_upstream_plus_modelled_WWTP","B4_upstream_plus_WWTP_and_reported_N",eligible)
      },
      "NEW_country_LEFT_OUT_per_country_original_source_cv":perCountry
    }
geoQA={
 "schema":"GESUI_454_EU_H7_90_H6_PARENT_15_COUNTRY_SPATIAL_TRANSFER_POSTHOC_v1",
 "original_H6_parent_mapping_source_SHA256":SHA["H6"],
 "original_HydroATLAS_source_SHA256":SHA["HydroATLAS"],
 "original_2018_22_reported_NP_source_SHA256":SHA["Reported_NP"],
 "original_source_level6_archive_SHA256":geo["source_original_level6_archive_SHA256"],
 "source_original_strict_no_repair_H6_GATE_failed_then_explicit_geometry_sensitivity":True,
 "same_original_454_basin_periods":N,"same_original_232_H7_basins":len(group),
 "original_H6_parent_count":G,
 "source_H6_parent_fold_sizes_H7_periods":np.bincount(folds,minlength=5).tolist(),
 "source_H6_parent_fold_sizes_groups":groups_per_fold,
 "source_original_16_H7_assigned_via_repair_group":len(source_repaired),
 "source_stable_country_group_n":dict(country_group_membership_counts),
 "source_country_reassignment_across_two_original_period_rows":reassigned_country_rows,
 "original_reported_N_source_positive_years":sum(z["sourceN_positive"] for z in rows),
 "country_leave_out_definition":"Group by fixed H7 stable country across both periods, omit country dummies because country is unseen, retain method controls",
 "new_geographic_results_are_POSTHOC_geographic_sensitivity_not_preregistered_independent_validation":True,
 "evaluation_folds_H6_parent_nested_basin_disjoint_more_strict_than_H7_ID_hash":True,
 "source_no_true_receiver_annual_NP_concentration":True,
 "source_plant_butterfly_are_GBIF_recording_indices_not_services":True,
 "models":results,
 "limitations":[
  "The original no-repair HydroBASINS6 gate FAILED 216/232; 207 originally invalid Level6 polygons were explicitly repaired in separate source-only sensitivity to recover all 232, not hidden as official valid native polygons.",
  "Level6 groups are hydrologically connected spatial blocks and are a stronger generalization test than H7 hash, but adjacent Level6 groups/countries are still not independent ecological regions.",
  "Leave-one-country-out excludes country fixed effects and may encounter unseen assessment method. This is a genuinely different estimand and cannot be compared as an identical covariate-adjusted model to the earlier within-country H7 CV.",
  "National river biological EQR assessment methods and original reporting are heterogeneous and France is heavily represented; many countries have <10 H7-period records.",
  "Only 56/454 original source reported annual N-positive records, 4 in France, 52 outside France; no matched report is not proof of absent WWTP.",
  "HydroWASTE static modeled outflow is accumulated to H7 basin OUTLET, not the interior EEA sampling station. Annual source N report converted by modeled flow is a hypothetical equivalent, not real TN.",
  "Source cropland/urban and natural outlet flow are original historically static HydroATLAS attributes, not measured 2018/2022 land-use trends.",
  "GBIF indexed plant S50 and Lycaenidae share of five families do not measure plant flowering, insect abundance, pollination or realized ecosystem functions.",
  "The exploratory geography/CV checks are posterior to previous same-source N/EQR associations; no external site or untouched independent continental replication."
 ]
}
outfile=ROOT/"H7_454_EU_H6_REPAIRED_SOURCE_SPATIAL_AND_COUNTRY_BENTHIC_PLANT_BUTTERFLY_RESULTS.json"
outfile.write_text(json.dumps(geoQA,sort_keys=True,indent=2,allow_nan=False))
print("H6_GEOBLOCK_COUNTRY_MACRO_RESULT_START")
compact={k:v for k,v in geoQA.items() if k not in ("limitations","models")}
compact["results"]={
 k:{
  "H6":v["NEW_H6_spatial_parent_comparisons"],
  "LOCO_weighted":v["NEW_country_LEFT_OUT_row_weighted_comparisons"],
  "LOCO_country_equal":v["NEW_country_LEFT_OUT_major_country_equal_MSE_gains"],
  "LOCO_country_rows":{country:{"period_rows":d["original_period_rows_heldout"],"sourcePositive":d["original_reported_NP_positive_H7_period_rows"],"WWTPgain":d["B2_to_B3_static_WWTP_gain_pct"],"N_gain":d["B3_to_B4_reported_N_gain_pct"]} for country,d in v["NEW_country_LEFT_OUT_per_country_original_source_cv"].items()}
  } for k,v in results.items()}
print(json.dumps(compact,sort_keys=True,allow_nan=False),flush=True)
print("H6_GEOBLOCK_COUNTRY_MACRO_RESULT_END")
