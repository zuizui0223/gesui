#!/usr/bin/env python3
"""Posthoc source-load magnitude permutation falsification for EEA aquatic 2022 EQR.

Fixed original EEA 2022 biological sites, country/method baseline and source-only
five geographic folds; no outcome-driven feature/threshold selection here.
Keep locations of registered upstream WWTP sources, country, and optionally
Q strata FIXED; shuffle only annual reported N/P load magnitudes among source-
positive sites. No ecological/wastewater causal significance claims.
"""
from __future__ import annotations
import argparse,hashlib,json,math,sys
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.pipeline import Pipeline
import sklearn

ORIG={
 "benthic_invertebrates":{"path":"BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv",
                          "sha":"5e1f3c31b389d0bd769984ab055d3e2c54b63ccaa4f5eedfb2d82bcdc83c28e6",
                          "n":2475,"exposed":147,
                          "expected_MSE":{"P":0.05769459552059828,"N":0.05752327324679533}},
 "phytobenthos":{"path":"PHYTOBENTHOS_2022_SOURCE_EQR.csv",
                          "sha":"7fa5dfff9023337ac0a62a87fe406eb71783979858d9d924ab97019c6e4e5fc3",
                          "n":2784,"exposed":110,
                          "expected_MSE":{"P":0.035875931550570526,"N":0.03587108895518606}},
}
SEED=20261108
N_PERMUTATIONS=199
ALPHA=10
SOURCES={"P":"upstream_2022_EEA_P_measured_t_per_yr",
         "N":"upstream_2022_EEA_N_measured_t_per_yr"}

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def input_matrix(df):
    x=pd.DataFrame({
      "lat":df.lat.to_numpy(float),
      "lon":df.lon.to_numpy(float),
      "logQ":np.log1p(df.DIS_AV_CMS.to_numpy(float)),
      "log_area":np.log1p(df.UPLAND_SKM.to_numpy(float)),
      "country":df.countryCode.astype("string").fillna("UNKNOWN").astype(str),
      "country_method":(df.countryCode.astype("string").fillna("UNKNOWN").astype(str)
                        +"|"+df.procedureClassificationSystem.astype("string").fillna("UNKNOWN").astype(str)),
    })
    return x

def get_group_indices(df,source,flow_stratified):
    country=df.countryCode.fillna("UNKNOWN").astype(str).to_numpy()
    logQ=np.log1p(df.DIS_AV_CMS.to_numpy(float))
    result=[]
    for c in sorted(set(country)):
        idx=np.flatnonzero((country==c)&(source>0))
        if not len(idx):continue
        if not flow_stratified or len(idx)<8:
            if len(idx)>=2:result.append(idx)
            continue
        # only source-exposed sites; country + 2 or 3 pre-specified Q strata.
        nbins=3 if len(idx)>=18 else 2
        ranks=pd.Series(logQ[idx]).rank(method="first")
        bins=pd.qcut(ranks,q=nbins,labels=False).to_numpy()
        for j in range(nbins):
            chunk=idx[bins==j]
            if len(chunk)>=2:result.append(chunk)
    all_indices=[i for group in result for i in group]
    if len(all_indices)!=len(set(all_indices)):
        raise ValueError("Exposed sites reused across permutation cells")
    return result

def heldout_MSE(x_base,y,fold_indices,source):
    x=x_base.copy()
    x["source_load"]=source
    numerical=["lat","lon","logQ","log_area","source_load"]
    categorical=["country","country_method"]
    preds=np.full(len(y),np.nan)
    for k in range(1,6):
        test=np.flatnonzero(fold_indices==k)
        train=np.flatnonzero(fold_indices!=k)
        if not len(test) or not len(train):
            raise ValueError("One of the pre-registered regional folds is empty")
        model=Pipeline([
          ("features",ColumnTransformer([
            ("numeric",StandardScaler(),numerical),
            ("categorical",OneHotEncoder(handle_unknown="ignore"),categorical)
          ])),
          ("ridge",Ridge(alpha=ALPHA))
        ])
        model.fit(x.iloc[train],y[train])
        preds[test]=model.predict(x.iloc[test])
    if not np.isfinite(preds).all():
        raise ValueError("Missing heldout validation prediction")
    return float(np.mean((y-preds)**2))

def permutation_qc(df,group,nutrient):
    nsource=ORIG[group]["exposed"]
    original=np.log1p(df[SOURCES[nutrient]].to_numpy(float))
    assert int((original>0).sum())==nsource
    base=input_matrix(df)
    response=df.normalised_EQR.to_numpy(float)
    fold=df.fixed_spatial_cv_fold.to_numpy(int)
    assert set(fold)=={1,2,3,4,5}
    observed=heldout_MSE(base,response,fold,original)
    expected=ORIG[group]["expected_MSE"][nutrient]
    if abs(observed-expected)>4e-6:
        raise ValueError("Original 2022 predeclared model does not numerically replay")
    studies={}
    for design in ("within_country_source_positive",
                   "within_country_and_source_positive_river_flow_strata"):
        cells=get_group_indices(df,original,flow_stratified=(design.endswith("_strata")))
        rng=np.random.default_rng(SEED)
        simulated=[]
        for b in range(N_PERMUTATIONS):
            modified=original.copy()
            for idx in cells:
                modified[idx]=original[rng.permutation(idx)]
            # The binary observed WWTP source locations and country/river Q
            # remain untouched; ONLY the source mass amount is randomized.
            assert np.array_equal(modified>0,original>0)
            simulated.append(heldout_MSE(base,response,fold,modified))
        vals=np.asarray(simulated,dtype=float)
        better=int((vals<=observed+1e-12).sum())
        studies[design]={
          "actual_2022_observed_MSE":observed,
          "null_median_MSE":float(np.median(vals)),
          "null_05th_MSE":float(np.quantile(vals,.05)),
          "null_95th_MSE":float(np.quantile(vals,.95)),
          "null_MSE_min":float(np.min(vals)),
          "null_MSE_max":float(np.max(vals)),
          "n_permutations":N_PERMUTATIONS,
          "null_randomized_better_or_equal_count":better,
          "monte_carlo_one_sided_extremeness_probability":
                 float((better+1)/(N_PERMUTATIONS+1)),
          "number_of_exchangeable_cells":len(cells),
          "number_of_exposed_sites_reassigned":sum(len(t) for t in cells),
          "source_positive_presence_preserved":True,
          "country_membership_preserved":True,
          "river_flow_strata_preserved":design.endswith("_strata"),
          "null_is_not_a_valid_causal_or_confirmatory_exchangeability_test":True
        }
        print("GESUI_2022_BQE_POSTHOC_CONDITIONAL_SOURCE_MASS_PERMUTATION",
             json.dumps({"BQE":group,"nutrient":nutrient,"design":design,
                          **studies[design]},sort_keys=True),flush=True)
    return studies

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--results-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result={}
    for group,source in ORIG.items():
        path=args.results_dir/source["path"]
        if sha(path)!=source["sha"]:
            raise ValueError("Original EQR source file digest does not match original retained data")
        df=pd.read_csv(path,low_memory=False,dtype={"site_key":str})
        if len(df)!=source["n"] or df.site_key.duplicated().any():
            raise ValueError("Original biological site units changed")
        if df.normalised_EQR.isna().any():
            raise ValueError("Original sample contains unreported EQR values")
        result[group]={nutrient:permutation_qc(df,group,nutrient)
                       for nutrient in ("P","N")}
    output={
      "schema":"gesui_2022_posthoc_reported_WWTP_mass_within_exposed_country_and_Q_permutation_v01",
      "status":"PASS_SOURCE_LOCKED_CONDITIONAL_PERMUTATION_EXPLORATORY",
      "original_EEA_BQE_source_run":42432030617,
      "original_2022_EQR_outcome_run":42432384196,
      "original_fixed_Ridge_alpha":ALPHA,
      "original_5_geographic_folds_reused":True,
      "scientific_scope":"annual WWTP reported N/P mass conditional prediction, not direct pollutant dose",
      "posthoc_not_preregistered":True,
      "random_seed":SEED,
      "source_mass_randomized_only_among_registered_positive_source_sites":True,
      "source_zero_sites_not_reassigned":True,
      "original_EEA_reported_N_and_P_mass_are_collinear":True,
      "country_and_Q_stratification_not_full_agriculture_or_urban_adjustment":True,
      "numpy_version":np.__version__,
      "sklearn_version":sklearn.__version__,
      "date":datetime.now(timezone.utc).isoformat(),
      "limitations":[
        "This null relies on unproven conditional exchangeability of upstream WWTP loads within the selected country/river flow stratum",
        "River network and catchment-source errors, population density and agricultural pressures are unadjusted",
        "Country-specific EQR method scales and original source missingness could produce spurious predictive ordering",
        "Repeated 2022 EQR sites are observational, not randomized sewage load treatments",
        "This permutation is posthoc because original EQR predictive results were already known; p is descriptive only",
        "N and P loads strongly correlate in the original EEA both-positive-source sample, so nutrient-specific mechanism unidentified"
      ],
      "results":result}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_BQE_PERMUTATION_SOURCE_LOCKED_FINAL",json.dumps({
       "status":output["status"],
       "tests":{g:{n:{k:v["monte_carlo_one_sided_extremeness_probability"]
          for k,v in r.items()} for n,r in d.items()} for g,d in result.items()}
       },sort_keys=True),flush=True)
if __name__=="__main__":
    main()
