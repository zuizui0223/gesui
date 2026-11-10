#!/usr/bin/env python3
"""Macro robustness of original global 2015 WW production × indexed plants.

Exact original 1596 WPI/WorldCover/Jones source join; original positive n1055
country-FE model MUST match before any posthoc geographic/country transfer
evaluation. Four directions fixed in docs; results regardless of sign.
Country FE usable for within-country spatial holdout, NEVER transferred to
new unseen countries. No all-flora or 2015 observed sewage nutrient claim.
"""
from __future__ import annotations
import argparse,hashlib,importlib.util,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import statsmodels.formula.api as smf

ORIGINAL_SCRIPT="archive/chocho_pr122/scripts/exploratory_global_wwtp_jones_worldcover.py"
ORIGINAL_AUDIT="archive/chocho_pr122/provenance/reviewer_defenses/results/global_hydrowaste_port_hosts_20261102/JONES2015_WORLDCOVER2021_GLOBAL_WW_SOURCE_AUDIT_20261103.json"
MODEL_CONTROL="abs_lat+sinlon+coslon+ln_total+ln_wwtp+built_fraction_2021+crop_fraction_2021"
NEW_PREDICTOR="ln_produced"
RNG_SEED=20261110
N_COUNTRY_BS=2000

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load_archived_script(root):
    f=root/ORIGINAL_SCRIPT
    spec=importlib.util.spec_from_file_location("historical_global_wwtp_original",str(f))
    if spec is None or spec.loader is None:raise RuntimeError("Archived original Global WWp parser unavailable")
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def build_original_sample(root,jones_path,worldcover_path):
    archive=load_archived_script(root)
    raw=archive.read_source(jones_path,worldcover_path)
    if len(raw)!=1596 or raw.port_id.nunique()!=1596 or raw.country.nunique()!=142:
        raise ValueError("Original globally sampled 1596 WPI source ports changed")
    if int((raw.quality=="COMPLETE").sum())!=1241:
        raise ValueError("Original ESA2021 source 1241 valid port samples changed")
    df=raw.loc[raw.in_within_country_comparison].copy()
    df=df.loc[(df.quality=="COMPLETE")&(~df.incomplete_or_missing_grid)].copy()
    df["abs_lat"]=df.latitude.abs()
    df["sinlon"]=np.sin(np.deg2rad(df.longitude))
    df["coslon"]=np.cos(np.deg2rad(df.longitude))
    df["ln_total"]=np.log1p(df.all_taxonomically_indexed.astype(float))
    df["ln_wwtp"]=np.log1p(df.WWTP_count_25km.astype(float))
    df["ln_produced"]=np.log1p(df.produced_25km_millionm3yr.astype(float))
    df["ratio"]=np.log((df.alien_candidates.astype(float)+.5)/(df.native_candidates.astype(float)+.5))
    df["ln_alien"]=np.log1p(df.alien_candidates.astype(float))
    df["ln_native"]=np.log1p(df.native_candidates.astype(float))
    df=df.replace([np.inf,-np.inf],np.nan).dropna(subset=[
        "abs_lat","sinlon","coslon","ln_total","ln_wwtp","ln_produced",
        "built_fraction_2021","crop_fraction_2021","ratio"])
    if (df[["built_fraction_2021","crop_fraction_2021"]].min().min()<0 or
        df[["built_fraction_2021","crop_fraction_2021"]].max().max()>1):
        raise ValueError("Unexpected original ESA2021 covariate scale")
    group=df.groupby("country").agg(n=("port_id","size"),wvari=("ln_wwtp","nunique"))
    valid=group.index[(group.n>=5)&(group.wvari>=2)]
    df=df.loc[df.country.isin(valid)].sort_values(["country","port_id"]).copy()
    if len(df)!=1055 or df.country.nunique()!=42 or df.port_id.nunique()!=1055:
        raise ValueError("Original source-locked 1055/42 country FE WWp population not recovered")
    if any(df.produced_25km_millionm3yr<0):
        raise ValueError("Modeled WWp production negative at original matched source")
    audit=json.loads((root/ORIGINAL_AUDIT).read_text())
    expected=audit["city_raster_2021_matched_models"]["WWTP_plus_WWproduction_plus_urban_crop"]
    form="ratio~"+MODEL_CONTROL+"+"+NEW_PREDICTOR+"+C(country)"
    full=smf.ols(form,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
    beta=float(full.params[NEW_PREDICTOR])
    ci=list(map(float,full.conf_int().loc[NEW_PREDICTOR]))
    if abs(beta-expected["WWproduction_beta"])>1e-9 or any(abs(a-b)>1e-6 for a,b in zip(ci,expected["WWproduction_CI95"])):
        raise ValueError("Original positive β_WWp or 42-country clustered CI not independently reproduced")
    if abs(float(full.params["ln_wwtp"])-expected["WWTP_beta"])>1e-9:
        raise ValueError("Original 2022 HydroWASTE source regression did not reconcile")
    return df,{
        "source_original_WPI_ports":1596,
        "source_original_WPI_countries":142,
        "matched_source_worldcover_valid_ports":1241,
        "source_original_countryFE_sample_ports":1055,
        "source_original_countryFE_countries":42,
        "original_β_WWp_after_urban_crop_countryFE":beta,
        "original_country_cluster95":ci,
        "original_country_cluster_nominal_p":float(full.pvalues[NEW_PREDICTOR]),
        "original_β_WWTP_count":float(full.params["ln_wwtp"]),
        "original_R2":float(full.rsquared),
        "prior_result_archive_numeric_reproduced":True,
    }

def positive_source_country_diagnostics(df):
    vals=[]
    for country,subset in df.groupby("country",sort=True):
        if len(subset)>=8 and subset.ln_produced.nunique()>=3 and subset.ratio.nunique()>=3:
            r=spearmanr(subset.ln_produced,subset.ratio).statistic
            vals.append({"country":str(country),"ports":int(len(subset)),"rho":float(r)})
    return {"n_eligible":len(vals),"n_pos":sum(x["rho"]>0 for x in vals),
            "n_neg":sum(x["rho"]<0 for x in vals),
            "n_exact_zero":sum(x["rho"]==0 for x in vals),
            "source_WWp_v_host_photo_fraction_within_country_spearman":vals}

def country_deletion(df):
    rows=[]
    model="ratio~"+MODEL_CONTROL+"+"+NEW_PREDICTOR+"+C(country)"
    for leave in sorted(df.country.unique()):
        d=df.loc[df.country!=leave]
        r=smf.ols(model,data=d).fit(cov_type="cluster",cov_kwds={"groups":d.country})
        rows.append({"left_out_country":str(leave),"remaining_ports":int(len(d)),
                     "remaining_countries":int(d.country.nunique()),
                     "re_fitted_β_WWp":float(r.params[NEW_PREDICTOR]),
                     "country_cluster95":list(map(float,r.conf_int().loc[NEW_PREDICTOR]))})
    return {"n_country_holdouts":len(rows),
            "positive_WWp_regression_after_one_country_deleted_n":sum(x["re_fitted_β_WWp"]>0 for x in rows),
            "min_beta":min(x["re_fitted_β_WWp"] for x in rows),
            "max_beta":max(x["re_fitted_β_WWp"] for x in rows),
            "models":rows,
            "not_unseen_country_prediction_because_countryFE_refitted_in_known_countries":True}

def repeated_country_geo_2degree(df):
    df=df.copy()
    df["original_country2deg_block"]=(
        df.country.astype(str)+"|"+
        np.floor(df.latitude.astype(float)/2).astype(int).astype(str)+"|"+
        np.floor(df.longitude.astype(float)/2).astype(int).astype(str))
    folds={}
    for country,g in df.groupby("country",sort=True):
        blocks=sorted(g.original_country2deg_block.unique())
        for i,name in enumerate(blocks):folds[name]=int(i%5)
    df["geo_fold"]=df.original_country2deg_block.map(folds).astype(int)
    formula0="ratio~"+MODEL_CONTROL+"+C(country)"
    formula1="ratio~"+MODEL_CONTROL+"+"+NEW_PREDICTOR+"+C(country)"
    fixed=[{"port_id":int(x.port_id),"country":str(x.country),
            "original_geo2deg_block":x.original_country2deg_block,
            "fold":int(x.geo_fold)}
           for x in df[["port_id","country","original_country2deg_block","geo_fold"]].itertuples(index=False)]
    errors=[]
    skipped=Counter()
    for fold in range(5):
        train=df[df.geo_fold!=fold]
        test=df[df.geo_fold==fold].copy()
        countries=set(train.groupby("country").size().loc[lambda s:s>=5].index)
        skip=test.loc[~test.country.isin(countries)]
        skipped["country_5_training_port_requirement"]+=int(len(skip))
        test=test.loc[test.country.isin(countries)]
        if test.empty:continue
        model0=smf.ols(formula0,data=train).fit()
        model1=smf.ols(formula1,data=train).fit()
        p0=model0.predict(test).to_numpy()
        p1=model1.predict(test).to_numpy()
        y=test.ratio.to_numpy()
        for i,(row,e0,e1) in enumerate(zip(test.itertuples(),(y-p0)**2,(y-p1)**2)):
            errors.append({"port_id":int(row.port_id),"country":str(row.country),
                           "fold":fold,"2deg_block":row.original_country2deg_block,
                           "baseline_squared_error":float(e0),
                           "plus_WWp_squared_error":float(e1)})
    if len(errors)<500 or len({x["country"] for x in errors})<20:
        return {"status":"SPATIAL_SOURCE_COUNTRY_TRAINING_SUPPORT_FAIL",
                "tested_ports":len(errors),"tested_countries":len({x["country"] for x in errors}),
                "original_country2degree_block_fold_assignments":fixed,
                "excluded_source_test_port_counts":dict(skipped),
                "DO_NOT_CLAIM_POSITIVE_OR_NEGATIVE_VALIDATION":True}
    return {"status":"COUNTRY_FIXED_EFFECT_WITHIN_COUNTRY_SPATIAL_HOLDOUT_ESTIMATED",
            "original_2degree_country_spatial_blocks_n":len(folds),
            "tested_2025_country_port_records":len(errors),
            "tested_original_countries":len({x["country"] for x in errors}),
            "excluded_source_test_port_counts":dict(skipped),
            "original_country2degree_block_fold_assignments":fixed,
            "squared_error_by_original_port":errors,
            "estimates":summarize_errors(errors),
            "known_country_geographic_transfer_only_NOT_new_country_inference":True}

def summarize_errors(rows):
    if not rows:return None
    sums=defaultdict(lambda:[0.,0.,0])
    b=sum(x["baseline_squared_error"] for x in rows)
    a=sum(x["plus_WWp_squared_error"] for x in rows)
    for x in rows:
        v=sums[x["country"]]
        v[0]+=x["baseline_squared_error"]
        v[1]+=x["plus_WWp_squared_error"]
        v[2]+=1
    coun=[
        {"country":c,"n_ports":int(x[2]),
         "MSE_baseline":float(x[0]/x[2]),
         "MSE_plus_WWp":float(x[1]/x[2]),
         "delta_MSE_base_minus_plus_WWp":float((x[0]-x[1])/x[2])}
        for c,x in sorted(sums.items())]
    base_country=sum(r["MSE_baseline"] for r in coun)/len(coun)
    add_country=sum(r["MSE_plus_WWp"] for r in coun)/len(coun)
    return {"weighted_by_PORT":{
             "MSE_baseline":b/len(rows),"MSE_plus_WWp":a/len(rows),
             "relative_MSE_gain_percent":100*(b-a)/b if b else None},
            "weighted_equal_COUNTRY":{
             "MSE_baseline":base_country,"MSE_plus_WWp":add_country,
             "relative_MSE_gain_percent":100*(base_country-add_country)/base_country if base_country else None},
            "countries_with_positive_WWp_predictive_gain":sum(r["delta_MSE_base_minus_plus_WWp"]>0 for r in coun),
            "countries_with_negative_WWp_predictive_gain":sum(r["delta_MSE_base_minus_plus_WWp"]<0 for r in coun),
            "per_country":coun}

def unseen_country_transfer(df):
    form0="ratio~"+MODEL_CONTROL
    form1="ratio~"+MODEL_CONTROL+"+"+NEW_PREDICTOR
    rows=[]
    for country in sorted(df.country.unique()):
        train=df.loc[df.country!=country]
        test=df.loc[df.country==country]
        if len(train)<500 or len(test)<5:raise ValueError("Country support changed")
        m0=smf.ols(form0,data=train).fit()
        m1=smf.ols(form1,data=train).fit()
        pred0=m0.predict(test).to_numpy()
        pred1=m1.predict(test).to_numpy()
        yy=test.ratio.to_numpy()
        for idx,(pid,e0,e1) in enumerate(zip(test.port_id,(yy-pred0)**2,(yy-pred1)**2)):
            rows.append({"port_id":int(pid),"country":str(country),"baseline_squared_error":float(e0),
                         "plus_WWp_squared_error":float(e1),
                         "external_fit_coefficient_WWp":float(m1.params[NEW_PREDICTOR])})
    if len(rows)!=len(df):raise ValueError("Unseen-country fit failed to score every original 1055 port")
    summary=summarize_errors(rows)
    # Resample COUNTRY (not individual ports), holding their OOS predictions fixed.
    countries=summary["per_country"]
    rng=np.random.default_rng(RNG_SEED)
    samples=[]
    for i in range(N_COUNTRY_BS):
        indices=rng.integers(0,len(countries),size=len(countries))
        draw=[countries[k] for k in indices]
        b=np.mean([x["MSE_baseline"] for x in draw])
        a=np.mean([x["MSE_plus_WWp"] for x in draw])
        samples.append(100*(b-a)/b if b else np.nan)
    samples=[float(x) for x in samples if math.isfinite(x)]
    assert len(samples)==N_COUNTRY_BS
    qs=np.quantile(samples,[.025,.975])
    summary["country_cluster_descriptive_2000_bootstrap95_relative_MSE_gain_percent"]=list(map(float,qs))
    summary["probabilistic_fraction_country_resamples_with_positive_gain"]=sum(x>0 for x in samples)/len(samples)
    return {"status":"STRICT_UNSEEN_COUNTRY_TRANSFER_NO_COUNTRYFE",
            "source_countries_evaluated":42,"out_of_country_scored_ports":len(rows),
            "country_FE_explicitly_excluded_as_unavailable_for_new_country":True,
            "squared_error_by_original_port":rows,"estimates":summary,
            "posthoc_generalization_not_confirmatory_preregistration":True}

def negative_controls(df):
    out={}
    for y in ("ln_alien","ln_native"):
        form=y+"~abs_lat+sinlon+coslon+ln_wwtp+ln_produced+built_fraction_2021+crop_fraction_2021+C(country)"
        m=smf.ols(form,data=df).fit(cov_type="cluster",cov_kwds={"groups":df.country})
        out[y]={"source_beta_WWp":float(m.params[NEW_PREDICTOR]),
                "country_cluster95":list(map(float,m.conf_int().loc[NEW_PREDICTOR])),
                "n_ports":len(df),"n_countries":int(df.country.nunique()),
                "source_metric_is_component_GBIF_indexed_candidate_richness_not_butterfly_or_all_plants":True}
    return out

def main():
    args=argparse.ArgumentParser()
    args.add_argument("--jones",type=Path,required=True)
    args.add_argument("--worldcover",type=Path,required=True)
    args.add_argument("--out",type=Path,required=True)
    a=args.parse_args()
    root=Path(__file__).resolve().parents[1]
    sample,original=build_original_sample(root,a.jones,a.worldcover)
    countrydiag=positive_source_country_diagnostics(sample)
    influence=country_deletion(sample)
    within=repeated_country_geo_2degree(sample)
    unseen=unseen_country_transfer(sample)
    controls=negative_controls(sample)
    out={"schema":"gesui_2015WWp_original_1055ports42country_2026_11_macro_spatial_country_transfer_posthoc_v01",
         "UTC":datetime.now(timezone.utc).isoformat(),
         "study_priority":"GLOBAL wastewater production vs GBIF-indexed introduced/native candidate larval host PLANT RECORD RATIO, not Lycaenidae",
         "source_PANGAEA2015_WWproduction_original_portCSV_SHA256":sha(a.jones),
         "source_original_WorldCover2021_portCSV_SHA256":sha(a.worldcover),
         "source_original_chocho_nested_model_archive_SHA256":sha(root/ORIGINAL_AUDIT),
         "original_source_cohort_replay":original,
         "original_countries_WWp_within_country_rank_pattern":countrydiag,
         "old_positive_model_one_country_removed_sensitivity":influence,
         "new_geo_country2degree_5fold_within_country_transfer":within,
         "new_leave_one_country_out_without_country_fixed_effects":unseen,
         "component_alien_and_native_candidate_indexing_controls":controls,
         "null_and_noData_never_imputed_as_WWp_zero":True,
         "result_crosssectional_noncausal_not_all_flora_or_direct_NP_mass":True,
         "study_all_new_outofsample_tests_posthoc_AFTER_original_β_positive_already_known":True,
         "sampling_heterogeneity_and_2015_2019_2021_2001_2024_date_mismatch":True,
         "all_source_1596_port_level_input_tables_NOT_republished":True}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print(json.dumps({
        "original":original,
        "country_WWp_rank":{k:v for k,v in countrydiag.items() if k!="source_WWp_v_host_photo_fraction_within_country_spearman"},
        "one_country_leave":{k:v for k,v in influence.items() if k!="models"},
        "within_country_geographical_transfer":within["estimates"] if "estimates" in within else within,
        "unseen_country_transfer":unseen["estimates"],
        "outcome_component_controls":controls,
    },indent=2))

if __name__=="__main__":
    main()
