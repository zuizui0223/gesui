#!/usr/bin/env python3
"""Independent, zero-network 2022 EEA aquatic EQR <-> reported WWTP N/P audit.

Replays all declared country/procedure-adjusted, spatial-fold-heldout Ridge
models with fixed source-only fold IDs. Uses *no* code imported from the original
model implementation, no outcome-driven model selection, no new EEA downloads.
Reports geographic block and country heterogeneity, source N/P collinearity,
and source-site attrition. A positive predictive gain is NOT causal sewage effect.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sys
from pathlib import Path
import numpy as np,pandas as pd,scipy,sklearn
from scipy.stats import spearmanr
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge

SOURCE_FILES={
 "BQE_SOURCE_EXPOSURE_COVERAGE.json":"a6d7a95b634ea796cefd091a2c0510a769a683358db86ed8b4d9e1bd01ea3f17",
 "EQR_2022_river_source_only_station_registry.csv":"e12df511d64c4c580ef2422dccdbbd11a1309c84d29d467fd4527206e41fa636"}
OUTCOME_FILES={
 "SOURCE_LOCKED_2022_AQUATIC_EQR_RESULTS.json":"7e647daa7f797740bd970e14079b6338add4a4faf04bf67f8eaac1ddd9e3f374",
 "benthic_invertebrates_EEA2022_EQR_source_matched_monitoring.csv":"5e1f3c31b389d0bd769984ab055d3e2c54b63ccaa4f5eedfb2d82bcdc83c28e6",
 "phytobenthos_EEA2022_EQR_source_matched_monitoring.csv":"7fa5dfff9023337ac0a62a87fe406eb71783979858d9d924ab97019c6e4e5fc3"}
EXPOSURES={
 "P":("upstream_2022_EEA_P_measured_t_per_yr","downstream_2022_EEA_P_t_per_yr_placebo"),
 "N":("upstream_2022_EEA_N_measured_t_per_yr","downstream_2022_EEA_N_t_per_yr_placebo")}
CLASSES={"benthic_invertebrates":"EEA_13-01-4",
         "phytobenthos":"EEA_124-04-9"}
NAMES={"M0":"M0_country_method_hydrogeo",
       "M1":"M1_plus_upstream_measured_load",
       "M2":"M2_plus_downstream_negative_control",
       "M3":"M3_both_upstream_and_downstream"}
MAX_DIFFERENCE=4e-6
BOOTSTRAP_DRAWS=5000
BOOTSTRAP_SEED=111

def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load_with_hash(folder,manifest):
    for name,sha in manifest.items():
        p=folder/name
        if hash_file(p)!=sha:
            raise ValueError(f"Source artifact changed: {name}")
    return True

def get_matrix(frame,up,down,variant):
    x=pd.DataFrame(index=frame.index)
    x["lat"]=frame.lat.astype(float)
    x["lon"]=frame.lon.astype(float)
    x["logQ"]=np.log1p(frame.DIS_AV_CMS.astype(float))
    x["logDrainage"]=np.log1p(frame.UPLAND_SKM.astype(float))
    cat=frame.countryCode.astype("string").fillna("UNKNOWN").astype(str)
    meth=frame.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)
    x["country"]=cat
    x["country_method"]=cat+"|"+meth
    nums=["lat","lon","logQ","logDrainage"]
    if variant in ("M1","M3"):
        x["upstream"]=np.log1p(frame[up].astype(float))
        nums.append("upstream")
    if variant in ("M2","M3"):
        x["downstream"]=np.log1p(frame[down].astype(float))
        nums.append("downstream")
    if not np.isfinite(x[nums].to_numpy(float)).all():
        raise ValueError("Unreliable source numerical predictor")
    return x,nums

def fit_original_folds(frame,nutrient,variant):
    up,down=EXPOSURES[nutrient]
    x,nums=get_matrix(frame,up,down,variant)
    y=frame.normalised_EQR.astype(float).to_numpy()
    fold=frame.fixed_spatial_cv_fold.to_numpy(dtype=int)
    if set(fold)!={1,2,3,4,5}:raise ValueError("Original spatial folds missing")
    pred=np.full(len(frame),np.nan,dtype=float)
    for i in range(1,6):
        train=np.where(fold!=i)[0]
        test=np.where(fold==i)[0]
        pipe=Pipeline([
          ("transformer",ColumnTransformer([
             ("numerical",StandardScaler(),nums),
             ("categorical",OneHotEncoder(handle_unknown="ignore"),
              ["country","country_method"])
          ],remainder="drop")),
          ("ridge",Ridge(alpha=10.0))])
        pipe.fit(x.iloc[train],y[train])
        pred[test]=pipe.predict(x.iloc[test])
    if not np.isfinite(pred).all():raise ValueError("Heldout predictions incomplete")
    return (y-pred)**2,abs(y-pred)

def source_integrity(frame,source,code,original_support):
    registry=source[source.observedPropertyDeterminandBiologyEQRCode.eq(code)].copy()
    if registry.duplicated("site_key").any() or frame.duplicated("site_key").any():
        raise ValueError("Repeated biological monitor site in frozen table")
    if not set(frame.site_key)<=set(registry.site_key):
        raise ValueError("Some source-only EQR monitoring sites replaced after observing outcomes")
    if len(registry)!=original_support["source_sites_before_numeric_quality"]:
        raise ValueError("Source denominator drifted")
    if len(frame)!=original_support["numeric_original_EQR_sites"]:
        raise ValueError("Original EQR measured sample size drifted")
    if frame.region_3deg_block.astype(str).isna().any():
        raise ValueError("Source region missing")
    for k in ("upstream_2022_EEA_N_measured_t_per_yr",
              "upstream_2022_EEA_P_measured_t_per_yr",
              "downstream_2022_EEA_N_t_per_yr_placebo",
              "downstream_2022_EEA_P_t_per_yr_placebo"):
        old=registry.set_index("site_key")[k]
        joined=frame.set_index("site_key")[k]
        if not np.allclose(old.reindex(joined.index).to_numpy(float),
                           joined.to_numpy(float),atol=1e-6,rtol=0):
            raise ValueError("WWTP source loading changed after outcome inspection: "+k)
    method=registry.set_index("site_key").procedureClassificationSystem
    if not method.reindex(frame.site_key).fillna("UNKNOWN").tolist() == frame.procedureClassificationSystem.fillna("UNKNOWN").tolist():
        raise ValueError("Biological procedure selection changed after outcome inspection")
    source_fold=registry.set_index("site_key").fixed_spatial_cv_fold
    if not np.array_equal(source_fold.reindex(frame.site_key).astype(int),
                          frame.fixed_spatial_cv_fold.astype(int)):
        raise ValueError("Original spatial fold altered in outcome data")
    missing=registry.loc[~registry.site_key.isin(frame.site_key)]
    return {
      "source_all_2022_river_snap_sites":len(registry),
      "real_outcome_sites":len(frame),
      "source_sites_not_in_outcome_panel":len(missing),
      "not_in_panel_by_reporting_country":{
          str(k):int(v) for k,v in missing.countryCode.value_counts().items()},
      "not_in_panel_by_source_observation_status":{
          str(k):int(v) for k,v in missing.resultObservationStatus.fillna("MISSING").value_counts().items()},
      "important":"Missing panel rows are NOT labeled invalid EQR without checking raw EEA values and source key joins; this is a post-outcome coverage audit."
    }

def summarize_one(df,original,source,code,label):
    if df.normalised_EQR.lt(0).any() or df.normalised_EQR.isna().any():
        raise ValueError("Non-finite or negative EQR in source outcome")
    if df.countryCode.nunique()!=original["support"]["country_count"]:
        raise ValueError("Original geographic denominator drift")
    srcqc=source_integrity(df,source,code,original["support"])
    blocks=df.region_3deg_block.astype(str).to_numpy()
    country=df.countryCode.astype(str).to_numpy()
    y=df.normalised_EQR.to_numpy(float)
    result={"support":srcqc,"models":{}}
    for nutrient,(up,down) in EXPOSURES.items():
        observations={}
        errors={}
        saved=original["numerical_prediction"][nutrient]["fixed_geographic_block_predictive_comparisons"]
        for name,targetname in NAMES.items():
            square,absolute=fit_original_folds(df,nutrient,name)
            mse=float(np.mean(square))
            old=float(saved["CV_MSE_by_model"][targetname]["overall_heldout_MSE"])
            if abs(mse-old)>MAX_DIFFERENCE:
                raise ValueError(f"Unreproduced heldout model {label}/{nutrient}/{name}: {mse-old}")
            observations[name]={"new_MSE":mse,"original_MSE":old,"absolute_replay_difference":abs(mse-old),
                                "new_MAE":float(np.mean(absolute))}
            errors[name]=square
        source_gain=errors["M0"]-errors["M1"]
        negative=errors["M0"]-errors["M2"]
        byblock=[]
        for k in sorted(set(blocks)):
            idx=np.flatnonzero(blocks==k)
            byblock.append({"block":k,"n":len(idx),"M0_minus_M1_mean_squared_error":float(source_gain[idx].mean()),
                            "M0_minus_M2_mean_squared_error":float(negative[idx].mean())})
        values=np.asarray([r["M0_minus_M1_mean_squared_error"] for r in byblock])
        weights=np.asarray([r["n"] for r in byblock])
        rng=np.random.default_rng(BOOTSTRAP_SEED)
        draw=rng.integers(0,len(values),size=(BOOTSTRAP_DRAWS,len(values)))
        estimates=(values[draw]*weights[draw]).sum(axis=1)/weights[draw].sum(axis=1)
        # A block bootstrap of FIXED model predictions: not an independent
        # refitting/bootstrap estimator for ecological causal effect.
        bycountry=[]
        for c in sorted(set(country)):
            s=country==c
            bycountry.append({"country":c,"station_n":int(s.sum()),
                              "upstream_measured_load_positive_n":int((df.loc[s,up]>0).sum()),
                              "MSE_baseline_minus_upstream":float(source_gain[s].mean())})
        positive=df[up]>0
        rho_s=spearmanr(df.loc[positive,"upstream_2022_EEA_N_measured_t_per_yr"],
                        df.loc[positive,"upstream_2022_EEA_P_measured_t_per_yr"]).statistic
        result["models"][nutrient]={
          "confirmed_model_MSE":observations,
          "raw_EQR_source_rank_rho_original":original["numerical_prediction"][nutrient]["source_EQR_spearman_rho_descriptive"],
          "replayed_MSE_upstream_improvement_pct":
             float(100*source_gain.mean()/observations["M0"]["new_MSE"]),
          "replayed_MSE_downstream_placebo_improvement_pct":
             float(100*negative.mean()/observations["M0"]["new_MSE"]),
          "outcome_independent_spatial_blocks":len(byblock),
          "blocks_where_upstream_model_improves":int((values>0).sum()),
          "countries_where_upstream_model_improves":sum(a["MSE_baseline_minus_upstream"]>0 for a in bycountry),
          "reporting_countries":len(bycountry),
          "posthoc_fixed_pred_block_bootstrap_95pct_MSE_reduction":
             [float(v) for v in np.quantile(estimates,[0.025,0.975])],
          "block_bootstrap_is_exploratory_not_causal":True,
          "source_P_and_N_collinearity_among_positive_upstream_sites":{
             "n":int(positive.sum()),"spearman_N_P":float(rho_s)},
          "per_country_MSE_diagnostics":bycountry,
          "per_block_MSE_diagnostics":byblock
        }
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-dir",required=True,type=Path)
    ap.add_argument("--outcome-dir",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    load_with_hash(args.source_dir,SOURCE_FILES)
    load_with_hash(args.outcome_dir,OUTCOME_FILES)
    srcqc=json.loads((args.source_dir/"BQE_SOURCE_EXPOSURE_COVERAGE.json").read_text())
    origin=json.loads((args.outcome_dir/"SOURCE_LOCKED_2022_AQUATIC_EQR_RESULTS.json").read_text())
    if origin["source_only_original_run"]!=42432030617 or origin["study_reference_year"]!=2022:
        raise ValueError("Original source/run identifier changed")
    source=pd.read_csv(args.source_dir/"EQR_2022_river_source_only_station_registry.csv",
                       dtype={"site_key":str},low_memory=False)
    if len(source)!=8245 or source[["site_key","observedPropertyDeterminandBiologyEQRCode"]].duplicated().any():
        raise ValueError("Unanticipated source registry duplicate/denominator")
    results={}
    for label,code in CLASSES.items():
        df=pd.read_csv(args.outcome_dir/(label+"_EEA2022_EQR_source_matched_monitoring.csv"),
                       dtype={"site_key":str},low_memory=False)
        if set(df.fixed_spatial_cv_fold.astype(int))!={1,2,3,4,5}:
            raise ValueError("Source-fold incomplete")
        results[label]=summarize_one(df,origin["results"][label],source,code,label)
        print("GESUI_INDEPENDENT_EQR_BLOCK_AUDIT",json.dumps({
          "group":label,"n":len(df),
          "source_missing_original":results[label]["support"]["source_sites_not_in_outcome_panel"],
          "P_improvement":results[label]["models"]["P"]["replayed_MSE_upstream_improvement_pct"],
          "N_improvement":results[label]["models"]["N"]["replayed_MSE_upstream_improvement_pct"],
          "P_block_positive":results[label]["models"]["P"]["blocks_where_upstream_model_improves"],
          "N_block_positive":results[label]["models"]["N"]["blocks_where_upstream_model_improves"]},
          sort_keys=True),flush=True)
    assert origin["results"]["fish"]["status"].startswith("NO_POOLED_ECOLOGICAL_TEST")
    report={
      "schema":"gesui_eea2022_real_aquatic_EQR_independent_source_holdout_audit_v01",
      "status":"PASS_INDEPENDENT_OFFLINE_RECONSTRUCTION_2022_BQE",
      "source_only_run":42432030617,
      "original_measured_BQE_run":42432384196,
      "source_files_sha256":SOURCE_FILES,
      "outcome_files_sha256":OUTCOME_FILES,
      "runtime":{"python":sys.version.split()[0],"numpy":np.__version__,
                 "pandas":pd.__version__,"scipy":scipy.__version__,
                 "scikit_learn":sklearn.__version__},
      "N_P_2022_source_only_same_facility_subset":True,
      "source_Biology_original_numeric_EQR_not_plant_species_abundance":True,
      "downstream_2022_reported_WWTP_mass_is_negative_control":True,
      "tested_all_predeclared_GEQ_BQE_classes":["benthic_invertebrates","phytobenthos","fish"],
      "fish_failed_original_geographic_coverage_gate":True,
      "hypothesis_exploratory_not_confirmed_causal_effect":True,
      "note":"Bootstrap based on fixed heldout errors is a descriptive heterogeneity check, not a causal CI nor a nested model-refit bootstrap.",
      "audited_groups":results}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_INDEPENDENT_EQR_2022_FINAL_AUDIT",
        json.dumps({"status":report["status"],"groups":list(results),
                    "source_provenance_verified":True,
                    "fish_original_test_skipped":True,
                    "external_network_calls":0},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
