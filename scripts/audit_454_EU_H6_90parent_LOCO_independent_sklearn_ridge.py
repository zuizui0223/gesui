#!/usr/bin/env python3
"""Independent sklearn Ridge CV verification of original H6-parent grouped and
stable-country leave-out macro biodiversity analyses from REAL source artifacts.

Numerically separate package (sklearn) from original manual numpy
linear-equation ridge. This is no new site evidence; all results POSTHOC.
"""
import collections,csv,hashlib,json,math,os,pathlib
import numpy as np
from sklearn.linear_model import Ridge

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H6_sklearn_independent"
ROOT.mkdir(parents=True,exist_ok=True)
SHA={
 "geo":"372b800b85abe2f395241987cea2000ff901cf5aff31bf1d1fb825304f2495fc",
 "atlas":"ce6780e01d2b98ed8e852574adb827cb3125be792a164e8be9e2ebb917b3a686",
 "source":"012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26"}
def read(k,sha=None):
    b=pathlib.Path(os.environ[k]).read_bytes()
    if sha and hashlib.sha256(b).hexdigest()!=sha:raise RuntimeError("Original published source identity mismatch "+k)
    return json.loads(b)
geo=read("ORIGINAL_H6_REPAIRED_SOURCE_MAP",SHA["geo"])["original_454_H7_study_basin_geographic_parent_mapping"]
h6={z["HYBAS_ID_L07"]:z["official_parent_HYBAS_ID_L06"] for z in geo}
atlas=read("ORIGINAL_H7_ATLAS_1219",SHA["atlas"])["original_1219_H7_covariate_by_real_ID"]
np_src=read("ORIGINAL_ANNUAL_NP_454",SHA["source"])["source_H7_period_reported_N_P_upstream_to_true_original_station"]
src={(z["HYBAS_ID_L07"],z["PERIOD"]):z for z in np_src}
prior=read("ORIGINAL_H6_MACRO_MODEL_RESULT")
if prior["schema"]!="GESUI_454_EU_H7_90_H6_PARENT_15_COUNTRY_SPATIAL_TRANSFER_POSTHOC_v1":raise RuntimeError("Original model changed")
with open(os.environ["ORIGINAL_454_SOURCE_BIOLOGICAL_ROWS"],newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
if len(rows)!=454 or len(h6)!=232:raise RuntimeError("Original 454 and 232 source grain changed")
COUNTRIES=sorted({r["country"] for r in rows})
METHODS=sorted({r["original_majority_BQE_method"] for r in rows})
h7countries=collections.defaultdict(collections.Counter)
for r in rows:h7countries[r["HYBAS_ID_L07"]][r["country"]]+=1
countryGroup={h:sorted(d,key=lambda k:(-d[k],k))[0] for h,d in h7countries.items()}
original_H6=[h6[r["HYBAS_ID_L07"]] for r in rows]
origCountry=[countryGroup[r["HYBAS_ID_L07"]] for r in rows]
folds=np.array([int(hashlib.sha256(x.encode()).hexdigest()[:8],16)%5 for x in original_H6])
basenum=len(rows)
def base(withcountry):
    rows0=[]
    for r in rows:
        features=[
        float(r["log_local_H7_area_km2"]),
        float(r["log1p_original_EQR_site_n"]),
        float(r["log1p_GBIF_Animalia_source_records"]),
        float(r["log1p_GBIF_Plantae_source_records"]),
        float(r["log1p_GBIF_5family_butterfly_source_records"]),
        float(r["log1p_original_EEA_site_DIS_AV_CMS"])]
        if withcountry:
            features += [float(r["country"]==c) for c in COUNTRIES[1:]]
        features += [float(r["original_majority_BQE_method"]==method) for method in METHODS[1:]]
        features += [float(r["PERIOD"]=="2020-2022")]
        rows0.append(features)
    return np.asarray(rows0,dtype=float)
def mat(which,withcountry):
    a=base(withcountry)
    additional=[]
    for r in rows:
        h=r["HYBAS_ID_L07"];v=atlas[h]
        if which.startswith("B2local") or which.startswith("B3local") or which.startswith("B4local"):
            env=[float(v["crp_pc_sse"]),float(v["urb_pc_sse"])]
        else:
            env=[float(v["crp_pc_use"]),float(v["urb_pc_use"])]
        env += [math.log1p(float(v["dis_m3_pyr"]))]
        if which.startswith("B3") or which.startswith("B4"):
            env += [float(r["W_static_upstream_MODEL_flow_per_UpArea"])]
        if which.startswith("B4"):
            q=src[(h,r["PERIOD"])]
            env += [math.log1p(float(q["source_year_reported_N_mgL_equivalent_site_equal_mean"]))]
        additional.append(env)
    return np.column_stack([a,np.asarray(additional,dtype=float)])
MODELS=["B2_upstream_crop_urban_naturalflow","B3_upstream_plus_modelled_WWTP",
        "B4_upstream_plus_WWTP_and_reported_N","B2local",
        "B3local_plus_modelled_WWTP","B4local_plus_WWTP_and_reported_N"]
def run(X,Y,groups):
    pred=np.empty(basenum)
    for group in sorted(set(groups)):
        mask=np.asarray(groups)==group
        train=~mask
        trX=X[train];teX=X[mask]
        mu=trX.mean(axis=0);sd=trX.std(axis=0)
        sd=np.where(sd>1e-12,sd,1.)
        A=(trX-mu)/sd;B=(teX-mu)/sd
        # Sklearn's fit_intercept=True is equivalent to the original manual
        # Ridge design with an unpenalized intercept and standardized X.
        reg=Ridge(alpha=10.0,fit_intercept=True,solver="svd")
        reg.fit(A,Y[train])
        pred[mask]=reg.predict(B)
    if not np.all(np.isfinite(pred)):raise RuntimeError("Independent model prediction nonfinite")
    return float(np.mean((Y-pred)**2))
audit={}
for label,col in [
 ("EEA_benthic_EQR","EQR_original_site_equal_basin_period_mean"),
 ("GBIF_plant_speciesKey_S50","plant_S50_GBIF_speciesKey"),
 ("GBIF_Lycaenidae_five_family_logodds","lycaenidae_fivefamily_index_share_logodds")
]:
    Y=np.array([float(r[col]) for r in rows])
    check={}
    for groupname,(withcountry,groups) in {
       "H6_parent_geo_block":(True,folds),
       "Country_heldout_with_country_FEs_removed":(False,origCountry)
    }.items():
        X={m:mat(m,withcountry) for m in MODELS}
        metrics={m:run(X[m],Y,groups) for m in MODELS}
        old=prior["models"][label][
            "NEW_H6_PARENT_BLOCKED_5fold_original_outcome_MSE" if groupname=="H6_parent_geo_block"
            else "NEW_country_LEFT_OUT_non_country_FE_MSE"]
        for m in MODELS:
            if not math.isclose(metrics[m],old[m],rel_tol=2e-8,abs_tol=1e-10):
                raise RuntimeError("Independent sklearn model differs from original manual Ridge: "+label+":"+groupname+":"+m+":"+str(metrics[m])+"/"+str(old[m]))
        def gain(a,b):
            return float(100*(metrics[a]-metrics[b])/metrics[a])
        check[groupname]={
           "all_original_6_model_MSE_exactly_matched":True,
           "source_GEOGRAPHIC_blocks_or_country_blocks":len(set(groups)),
           "EQR_or_plant_or_lyc_b2_to_b3_static_WWTP_gain_pct":gain(MODELS[0],MODELS[1]),
           "EQR_or_plant_or_lyc_b3_to_b4_annual_reported_N_gain_pct":gain(MODELS[1],MODELS[2]),
           "local_cropland_urban_b2_to_b3_static_WWTP_gain_pct":gain(MODELS[3],MODELS[4]),
           "local_cropland_urban_b3_to_b4_annual_N_gain_pct":gain(MODELS[4],MODELS[5]),
           "original_model_baseline_and_new_source_MSE_values":metrics
        }
    audit[label]=check
output={
 "schema":"GESUI_H6_90PARENT_454EU_3LAYER_SKLEARN_INDEPENDENT_CV_QA_20261113_v1",
 "source_H6_geo_SHA256":SHA["geo"],"source_H7_HydroATLAS_SHA256":SHA["atlas"],
 "source_EEA_reported_annual_N_SHA256":SHA["source"],
 "original_exact_real_H7_basins":len(h6),
 "original_same_basin_period_rows":len(rows),
 "original_H6_parent_groups":len(set(original_H6)),
 "original_country_heldout_stable_groups":len(set(origCountry)),
 "all_original_manual_Ridge10_H6_parent_and_country_CV_model_values_reproduced_with_independent_sklearn":True,
 "models_three_separate_biodiversity_responses":audit,
 "warnings":[
  "Independent numerical implementation, not independent biological study or prospectively locked ecological predictions.",
  "Country heldout drops country FE; original country methods can be unseen and each source group differs in sampling/reporting.",
  "16 original H7-to-H6 parents required explicitly corrected invalid source geometry after strict original H6 source gate failed.",
  "Neither native H6 group split nor country omit identifies causal sewage impact or cross-boundary ecosystem service."
 ]}
path=ROOT/"H6_90PARENT_EU454_EEA_PLANT_LYCAENIDAE_SKLEARN_INDEPENDENT_REPRODUCTION_QA.json"
path.write_text(json.dumps(output,indent=2,sort_keys=True))
brief={k:v for k,v in output.items() if k not in ("warnings",)}
for z in brief["models_three_separate_biodiversity_responses"].values():
    for x in z.values():x.pop("original_model_baseline_and_new_source_MSE_values",None)
print("H6_INDEPENDENT_SKLEARN_CV_QA_START")
print(json.dumps(brief,sort_keys=True))
print("H6_INDEPENDENT_SKLEARN_CV_QA_END")
