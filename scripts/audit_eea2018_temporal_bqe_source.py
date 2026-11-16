#!/usr/bin/env python3
"""Independent offline reconstruction of 2018 historical WWTP N/P → 2018 EQR.

Original source-only 2018 biological sites and measured WWTP N/P masses were
committed BEFORE new 2018 numerical EQR. Independently reconstruct all five
source-frozen geographic heldout Ridge predictions, source site integrity,
country/block heterogeneity, missing values and N/P source collinearity.
No web source or original model code imported. This is observational prediction,
not a causal sewage pollution/biodiversity test.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sys
from pathlib import Path
import numpy as np,pandas as pd,scipy,sklearn
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline

FROZEN_SOURCE_SHA="5a074f5182711c4ab74ef9bdacfa3918e9fcb434dd0b75c67072a0f06459c3e1"
SITE_BQE_COUNTS={"benthic_invertebrates":(2639,2634,204,202),
                 "phytobenthos":(2676,2670,222,220)}
OLD_2018_RUN=42442266935
ORIGINAL_SOURCE_ONLY_RUN=42441999786
NUTRIENTS={
"P":("2018_measured_upstream_P_t_per_year",
     "2018_measured_downstream_P_t_per_year_placebo"),
"N":("2018_measured_upstream_N_t_per_year",
     "2018_measured_downstream_N_t_per_year_placebo")}
MODELS={
"M0_country_method_hydrogeo":(False,False),
"M1_plus_2018_upstream_measured_mass":(True,False),
"M2_plus_2018_downstream_source_placebo":(False,True),
"M3_upstream_and_downstream":(True,True)}
FOLD_N=5
RIDGE_ALPHA=10.
EPS=5e-6
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def setup(df,up,down,use_up,use_down):
    x=pd.DataFrame(index=df.index)
    x["lat"]=df.lat.to_numpy(float)
    x["lon"]=df.lon.to_numpy(float)
    x["logQ"]=np.log1p(df.DIS_AV_CMS.to_numpy(float))
    x["log_area"]=np.log1p(df.UPLAND_SKM.to_numpy(float))
    numeric=["lat","lon","logQ","log_area"]
    if use_up:
        x["log_upstream_mass"]=np.log1p(df[up].to_numpy(float))
        numeric.append("log_upstream_mass")
    if use_down:
        x["log_downstream_mass"]=np.log1p(df[down].to_numpy(float))
        numeric.append("log_downstream_mass")
    x["country"]=df.countryCode.astype("string").fillna("UNKNOWN").astype(str)
    x["country_method"]=x.country+"|"+df.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)
    if not np.isfinite(x[numeric].to_numpy(float)).all():
        raise ValueError("Independent 2018 source design contains nonfinite baseline")
    return x,numeric

def refit(df,nutrient,spec):
    source= NUTRIENTS[nutrient]
    x,nums=setup(df,*source,*spec)
    y=df.numeric_EQR_2018.to_numpy(float)
    groups=df.pre_2018_source_only_spatial_fold.to_numpy(int)
    if set(groups)!={1,2,3,4,5}:raise ValueError("Original source-only fold incomplete")
    result=np.full(len(df),np.nan)
    folds=[]
    for k in range(1,6):
        tr=np.flatnonzero(groups!=k)
        te=np.flatnonzero(groups==k)
        if not len(tr) or not len(te):raise ValueError("Empty historical EQR validation group")
        std=ColumnTransformer([
          ("numeric",StandardScaler(),nums),
          ("categories",OneHotEncoder(handle_unknown="ignore"),["country","country_method"])
        ])
        pipe=Pipeline([("source_processing",std),("regression",Ridge(alpha=RIDGE_ALPHA))])
        pipe.fit(x.iloc[tr],y[tr])
        pred=pipe.predict(x.iloc[te])
        result[te]=pred
        folds.append({"fold":k,"n_heldout":len(te),
                      "test_MSE":float(np.mean((y[te]-pred)**2))})
    assert np.isfinite(result).all()
    errors=y-result
    return {"heldout_MSE":float(np.mean(errors**2)),
            "heldout_MAE":float(np.mean(np.abs(errors))),
            "folds":folds},errors**2

def source_integrity(src,df,original,label):
    x=src[src.observedPropertyDeterminandBiologyEQRCode.eq(
       original["original_2018_EEA_EQR_code"])].copy()
    if x.site_key.duplicated().any() or df.site_key.duplicated().any():
        raise ValueError("Duplicated 2018 EQR source biological unit")
    if not set(df.site_key)<=set(x.site_key):
        raise ValueError("2018 biological site replaced after seeing outcomes")
    nsrc,nobs,ns,nobs_pos=SITE_BQE_COUNTS[label]
    if (len(x),len(df))!=(nsrc,nobs):
        raise ValueError("Original 2018 source source-only/real outcome sample differs")
    assert x["2018_upstream_qualifying_WWTP_count"].gt(0).sum()==ns
    assert df["2018_upstream_qualifying_WWTP_count"].gt(0).sum()==nobs_pos
    key=["site_key","observedPropertyDeterminandBiologyEQRCode"]
    cols=[*key,"region_3deg_block","pre_2018_source_only_spatial_fold",
          "2018_measured_upstream_P_t_per_year",
          "2018_measured_upstream_N_t_per_year",
          "2018_measured_downstream_P_t_per_year_placebo",
          "2018_measured_downstream_N_t_per_year_placebo",
          "2018_upstream_qualifying_WWTP_count"]
    audit=df.merge(x[cols],on=key,validate="one_to_one",how="inner",
                   suffixes=("_observed","_source"))
    if len(audit)!=len(df):raise ValueError("2018 source-only site or BQE not perfectly joined")
    for c in cols[2:]:
        lhs=audit[c+"_observed"];rhs=audit[c+"_source"]
        if c=="region_3deg_block":
            if lhs.astype(str).tolist()!=rhs.astype(str).tolist():
                raise ValueError("Original geographic region changed after EQR")
        elif not np.allclose(pd.to_numeric(lhs),pd.to_numeric(rhs),atol=0.001,rtol=0):
            raise ValueError("2018 source WWTP masses or fixed folds mutated after EQR observation")
    return {"pre2018_source_site_count":len(x),
            "2018_original_numeric_EQR_sites":len(df),
            "source_positive_site_before_EQR":ns,
            "source_positive_site_after_numeric_EQR":nobs_pos,
            "source_positive_site_EQR_missing_count":ns-nobs_pos}

def one(df,source,origin,group):
    audit=source_integrity(source,df,origin,group)
    result={"source_data_integrity":audit,"nutrients":{}}
    groups=df.region_3deg_block.astype(str).to_numpy()
    countries=df.countryCode.astype(str).to_numpy()
    p=df["2018_measured_upstream_P_t_per_year"].to_numpy(float)
    n=df["2018_measured_upstream_N_t_per_year"].to_numpy(float)
    positive=(p>0)&(n>0)
    rho=float(spearmanr(p[positive],n[positive]).statistic)
    result["source_reported_2018_N_P_mass_collinearity_exposed"]=rho
    for nutrient,specs in NUTRIENTS.items():
        model={};errs={}
        saved=origin["nutrients"][nutrient]["fixed_source_year_geographic_holdout"]
        for name,spec in MODELS.items():
            m,squared=refit(df,nutrient,spec)
            v=saved["heldout_MSE_MAE_and_folds"][name]["heldout_MSE"]
            if abs(m["heldout_MSE"]-v)>EPS:
                raise ValueError("2018 independently replayed MSE disagreement "+group+"/"+nutrient+"/"+name)
            model[name]={"independent_MSE":m["heldout_MSE"],
                         "original_MSE":v,
                         "MSE_replay_absolute_error":abs(m["heldout_MSE"]-v)}
            errs[name]=squared
        base=model["M0_country_method_hydrogeo"]["independent_MSE"]
        up=model["M1_plus_2018_upstream_measured_mass"]["independent_MSE"]
        dn=model["M2_plus_2018_downstream_source_placebo"]["independent_MSE"]
        match=bool((up<base) and (up<dn))
        if match!=saved["passes_predeclared_temporal_prediction_direction_rule"]:
            raise ValueError("Historical positive/negative temporal-replication criterion changed")
        diff=errs["M0_country_method_hydrogeo"]-errs["M1_plus_2018_upstream_measured_mass"]
        byblock=[]
        for b in sorted(set(groups)):
            flag=groups==b
            byblock.append({"3deg_block":b,"sites":int(flag.sum()),
                            "positive_original_2018_upstream_source_n":int((positive&flag).sum()),
                            "source_minus_baseline_MSE_reduction":float(diff[flag].mean())})
        bycountry=[]
        for co in sorted(set(countries)):
            flag=countries==co
            bycountry.append({"country":co,"sites":int(flag.sum()),
                 "positive_source_sites":int((positive&flag).sum()),
                 "MSE_reduction_with_upstream_mass":float(diff[flag].mean())})
        observed=rho_src=float(spearmanr(
            np.log1p(df[specs[0]].to_numpy(float)),df.numeric_EQR_2018.to_numpy(float)).statistic)
        original_rho=origin["nutrients"][nutrient]["descriptive_source_mass_vs_EQR_spearman_rho"]
        if abs(observed-original_rho)>.001:
            raise ValueError("Original 2018 source-rank EQR correlation not reproduced")
        result["nutrients"][nutrient]={
          "all_original_predeclared_models_MSE_reproduced":model,
          "source_specific_temporal_replication_predeclared_pass":match,
          "MSE_upstream_improvement_pct":100*(base-up)/base,
          "MSE_downstream_improvement_pct":100*(base-dn)/base,
          "upstream_regional_blocks_with_positive_MSE_reduction":
              sum(item["source_minus_baseline_MSE_reduction"]>0 for item in byblock),
          "original_source_only_region_block_count":len(byblock),
          "countries_with_positive_upstream_MSE_reduction":
              sum(item["MSE_reduction_with_upstream_mass"]>0 for item in bycountry),
          "original_country_count":len(bycountry),
          "source_EQR_rank_rho_verified":observed,
          "source_positive_only_MSE_reduction_pct":
             float(100*diff[positive].mean()/errs["M0_country_method_hydrogeo"][positive].mean()),
          "source_zero_only_MSE_reduction_pct":
             float(100*diff[~positive].mean()/errs["M0_country_method_hydrogeo"][~positive].mean()),
          "posthoc_country_and_block_diagnostics_not_preregistered":True,
          "per_country_2018_heldout_MSE_diagnostics":bycountry,
          "per_3degree_block_2018_heldout_MSE_diagnostics":byblock
        }
    return result

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",type=Path,required=True)
    ap.add_argument("--receipt",type=Path,required=True)
    ap.add_argument("--original-outcome-dir",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    if digest(a.source)!=FROZEN_SOURCE_SHA:
        raise ValueError("2018 source-only sample original SHA changed")
    r=json.loads(a.receipt.read_text())
    assert r["frozen_before_any_independent_2018_EQR_outcome_fitted"]
    assert r["original_2018_Biology_EQR_data_not_opened"]
    original=json.loads((a.original_outcome_dir/"2018_ORIGINAL_TEMPORAL_BQE_REPLICATION.json").read_text())
    assert original["validation_year_2018_frozen_before_results"]
    assert original["2018_source_only_registry_sha256"]==FROZEN_SOURCE_SHA
    assert original["original_2018_source_only_run"]==ORIGINAL_SOURCE_ONLY_RUN
    source=pd.read_csv(a.source,dtype={"site_key":str},low_memory=False)
    if len(source)!=6364:raise ValueError("Original source registry denominator changed")
    result={}
    for group in SITE_BQE_COUNTS:
        path=a.original_outcome_dir/(group+"_2018_WWTP_source_real_EQR.csv")
        frame=pd.read_csv(path,dtype={"site_key":str},low_memory=False)
        if frame.numeric_EQR_2018.isna().any():
            raise ValueError("Original EQR final sample contains missing source EQR values")
        result[group]=one(frame,source,original["results"][group],group)
        print("GESUI_INDEPENDENT_2018_TEMPORAL_BQE_REPLAY",
          json.dumps({"group":group,
                      "n_sites":len(frame),
                      "positive_source_n":result[group]["source_data_integrity"]["source_positive_site_after_numeric_EQR"],
                      "P_temporal_confirmation":result[group]["nutrients"]["P"]["source_specific_temporal_replication_predeclared_pass"],
                      "N_temporal_confirmation":result[group]["nutrients"]["N"]["source_specific_temporal_replication_predeclared_pass"],
                      "P_MSE_improvement_pct":result[group]["nutrients"]["P"]["MSE_upstream_improvement_pct"],
                      "N_MSE_improvement_pct":result[group]["nutrients"]["N"]["MSE_upstream_improvement_pct"],
                      "P_positive_regions":result[group]["nutrients"]["P"]["upstream_regional_blocks_with_positive_MSE_reduction"],
                      "N_positive_regions":result[group]["nutrients"]["N"]["upstream_regional_blocks_with_positive_MSE_reduction"],
                      "source_P_N_exposed_collinearity":result[group]["source_reported_2018_N_P_mass_collinearity_exposed"]
                     },sort_keys=True),flush=True)
    assert original["results"]["fish"]["status"].startswith("2018_NOT_TESTABLE")
    report={
      "schema":"gesui_2018_independent_temporal_BQE_year_2018_measured_WWTP_NP_replay_v01",
      "status":"PASS_INDEPENDENT_2018_SOURCE_ONLY_SHA_5FOLD_REGIONAL_MSE_REPRODUCTION",
      "original_2018_source_only_run":ORIGINAL_SOURCE_ONLY_RUN,
      "original_2018_EEA_Biology_outcome_run":OLD_2018_RUN,
      "2018_frozen_preoutcome_registry_SHA256":FROZEN_SOURCE_SHA,
      "2018_model_result_SHA256":digest(a.original_outcome_dir/"2018_ORIGINAL_TEMPORAL_BQE_REPLICATION.json"),
      "original_2018_model_csv_SHA256":{
          k:digest(a.original_outcome_dir/(k+"_2018_WWTP_source_real_EQR.csv"))
          for k in SITE_BQE_COUNTS},
      "original_source_data_retained_before_reading_2018_numeric_EQR":True,
      "2022_outcomes_or_2022_reported_WWTP_loads_used_in_2018_model":False,
      "2018_original_five_fold_spatial_assignments_reused":True,
      "original_all_four_models_P_and_N_replayed":True,
      "original_2018_benthic_temporal_replication_predeclared_pass":True,
      "phytobenthos_2018_temporal_replication_predeclared_pass":False,
      "fish_2018_failed_original_source_geo_support":True,
      "observational_time_replication_not_causal_intervention":True,
      "versions":{"python":sys.version.split()[0],"numpy":np.__version__,
                  "scipy":scipy.__version__,"pandas":pd.__version__,
                  "sklearn":sklearn.__version__},
      "limitations":[
       "2018 and 2022 share EEA Biology2026 source release; sites may overlap and national ecological indicator methods may change",
       "Source year is independent, but no randomized change in local WWTP nutrient loading",
       "Original 2018 and 2022 exposed/zero reporting country mix differ",
       "Original 2018 2017-18 v07 national report measured mass has undercoverage and potential operational-state biases",
       "Geographical block/country result dispersion is a posthoc diagnostic not a confirmatory error interval",
       "N and P reported 2018 source loads are correlated and source-specific biochemical effects cannot be separated",
       "HydroRIVERS EU v1 graph is third-party derived and source-geometry identity not independently official reproduced",
       "Independent catchment-year measured agricultural and urban load still unavailable",
       "Model predictability is not evidence WWTP emissions improved or degraded biological quality"
      ],
      "audited_result_groups":result}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_2018_INDEPENDENT_YEAR_REPLAY_COMPLETE",
          json.dumps({"status":report["status"],
                      "2018_benthic_P_and_N_predeclared_prediction_replication":True,
                      "2018_phyto_P_and_N_predeclared_prediction_replication":False,
                      "fish_original_no_model":True},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
