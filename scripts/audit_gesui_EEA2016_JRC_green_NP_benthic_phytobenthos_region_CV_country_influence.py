#!/usr/bin/env python3
"""Independent offline EEA2016 GREEN 2014–18 regional source posthoc QA.

Uses ONLY already-archived 2016 original region-heldout aggregate EQR
prediction errors plus verified pre-outcome 2016 NUTS2 source manifest.
Recompute each MSE gain by site, region, country, fold, leave-country-out;
N/P regional predictor collinearity and country influence. No new EQR or
sewage measurements accessed. POSTHOC descriptive influence, no claim of
prospectively registered independent 2016 WWTP treatment effect.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

D=Path("results/eea2016_thirdyear_source_20261111")
RESULT=D/"2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json"
SOURCE=D/"2016_EEA_BENTHIC_GISCO_NUTS2_GREEN_POINT_NP_SOURCE_GATE_BEFORE_EQR.json"
FOLDS=D/"2016_EEA_BENTHIC_JRC_NUTS2_FROZEN_5FOLDS_NO_EQR.json"
MODELS=("baseline","upstream_modelPointN","upstream_modelPointP")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gain(a,b):return 100*(a-b)/a
def wmean(rows,key,weight):
    num=sum(float(x[key])*float(x[weight]) for x in rows)
    den=sum(float(x[weight]) for x in rows)
    if den<=0:raise ValueError("Source NUTS2 region case weight must be positive")
    return num/den
def unweighted(rows,key):
    if not rows:raise ValueError("Source missing")
    return sum(float(z[key]) for z in rows)/len(rows)
def check(root):
    o=json.loads((root/RESULT).read_text())
    source=json.loads((root/SOURCE).read_text())
    folds=json.loads((root/FOLDS).read_text())
    if o["original_2016_sites_GISCO2016_JRC_modelled_SOURCE_before_EQR_SHA256"]!=sha(root/SOURCE) or o["original_original5fold_SOURCE_before_EQR_SHA256"]!=sha(root/FOLDS):
        raise ValueError("Independent 2016 outcome result does not match frozen source cohort and folds")
    outcome={}
    for name in ("benthic_invertebrates","phytobenthos"):
        bio=o["results"][name]
        if not bio["precommitted_numeric_EQR_support_gate_PASS"]:
            raise ValueError("No complete source-supported real BQE model to audit "+name)
        cv=bio["original_NUTS2_disjoint_holdout_Ridge10"]
        rr=cv["source_JRC_green_region_sampled_EEA_EQR_summary"]
        codes={x["nuts2"] for x in rr}
        expected=set(folds["source_biological_classes"][name]["fixed_region_to_fivefold_source_only"])
        if not codes.issubset(expected) or len(codes)!=cv["source_test_independent_NUTS2_regions"]:
            raise ValueError("Original 2016 2014–18 JRC NUTS2 BQE grouped source mismatch")
        if len({x["nuts2"] for x in rr})!=len(rr):raise ValueError("Duplicate original NUTS2 result")
        for row in rr:
            if row["source_frozen_test_fold"]!=folds["source_biological_classes"][name]["fixed_region_to_fivefold_source_only"][row["nuts2"]]:
                raise ValueError("Same NUTS2 2016 source region changed frozen fold after biological outcome")
        n=sum(x["n_original_2016_BQE_sites"] for x in rr)
        if n!=bio["precommitted_source_country_region_outcome_QC"]["original_EEA2016_numeric_EQR_valid_sites"]:
            raise ValueError("Original 2016 numeric EQR site sample changed after CV result")
        err={}
        for key in MODELS:
            accessor=lambda x: x["CV_region_MSE"][key]
            site=sum(x["n_original_2016_BQE_sites"]*accessor(x) for x in rr)/n
            reg=sum(accessor(x) for x in rr)/len(rr)
            if abs(site-cv["MSE_site_equal"][key])>1e-10 or abs(reg-cv["MSE_NUTS2_equal"][key])>1e-10:
                raise ValueError("Independent MSE site/region raw source arithmetic mismatch: "+name+"/"+key)
            err[key]={"site_equal":site,"NUTS2_equal":reg}
        bycountry=defaultdict(list)
        byfold=defaultdict(list)
        for x in rr:
            bycountry[x["country"]].append(x)
            byfold[x["source_frozen_test_fold"]].append(x)
        stats={}
        for country,part in sorted(bycountry.items()):
            val={key:sum(x["n_original_2016_BQE_sites"]*x["CV_region_MSE"][key] for x in part)/
                 sum(x["n_original_2016_BQE_sites"] for x in part) for key in MODELS}
            stats[country]={
                "original_2016_EEA_benthic_numeric_sites_in_country":sum(x["n_original_2016_BQE_sites"] for x in part),
                "original_green_NUTS2_independent_region_count":len(part),
                "NUP_baseline_MSE":val["baseline"],
                "NUP_N_2014to18_pointsource_OOS_gain_percent":gain(val["baseline"],val["upstream_modelPointN"]),
                "NUP_P_2014to18_pointsource_OOS_gain_percent":gain(val["baseline"],val["upstream_modelPointP"]),
                "not_country_randomized_WWTP_treatment":True
            }
        lof={}
        for country in sorted(bycountry):
            rest=[x for x in rr if x["country"]!=country]
            val={}
            for key in MODELS:
                val[key]=sum(x["n_original_2016_BQE_sites"]*x["CV_region_MSE"][key] for x in rest)/sum(x["n_original_2016_BQE_sites"] for x in rest)
            lof[country]={"remaining_site_n":sum(x["n_original_2016_BQE_sites"] for x in rest),
                "N_source_prediction_gain_percent":gain(val["baseline"],val["upstream_modelPointN"]),
                "P_source_prediction_gain_percent":gain(val["baseline"],val["upstream_modelPointP"])}
        foldgain={}
        for fold,part in sorted(byfold.items()):
            val={key:sum(x["n_original_2016_BQE_sites"]*x["CV_region_MSE"][key] for x in part)/sum(x["n_original_2016_BQE_sites"] for x in part) for key in MODELS}
            foldgain[str(fold)]={
                "test_sites":sum(x["n_original_2016_BQE_sites"] for x in part),
                "test_independent_NUTS2_regions":len(part),
                "N_gain_pct":gain(val["baseline"],val["upstream_modelPointN"]),
                "P_gain_pct":gain(val["baseline"],val["upstream_modelPointP"]),
            }
        source_N=np.array([x["mean_logPointN_modelled2014to18_density"] for x in rr])
        source_P=np.array([x["mean_logPointP_modelled2014to18_density"] for x in rr])
        meanEQ=np.array([x["mean_EQR_original_numeric"] for x in rr])
        spearman_NP=float(spearmanr(source_N,source_P).statistic)
        reg_loci={"modelled_regionequal_N_EQR_spearman":float(spearmanr(source_N,meanEQ).statistic),
                  "modelled_regionequal_P_EQR_spearman":float(spearmanr(source_P,meanEQ).statistic)}
        outcome[name]={
            "original_2016_EEA_numeric_BQE_sites":n,
            "original_2016_EEA_benthic_NUTS2_region_count":len(rr),
            "original_EEA_COUNTRY_count":len(bycountry),
            "original_NUTS2_frozen_folds_must_not_overlap":True,
            "independent_siteequal_and_NUTS2equal_MSE":err,
            "country_2016_modelled_load_forecast_influence":stats,
            "leave_one_country_out_on_ALREADY_FROZEN_CV_PREDICTIONS_not_retrained":lof,
            "N_predictor_LOO_country_positive_count":sum(x["N_source_prediction_gain_percent"]>0 for x in lof.values()),
            "P_predictor_LOO_country_positive_count":sum(x["P_source_prediction_gain_percent"]>0 for x in lof.values()),
            "N_predictor_country_with_positive_gain_count":sum(x["NUP_N_2014to18_pointsource_OOS_gain_percent"]>0 for x in stats.values()),
            "P_predictor_country_with_positive_gain_count":sum(x["NUP_P_2014to18_pointsource_OOS_gain_percent"]>0 for x in stats.values()),
            "five_frozen_NUTS2_fold_honest_source_prediction_gains":foldgain,
            "regional_original_2014to18_JRC_modelled_PsN_vs_PsP_spearman_collinearity":spearman_NP,
            "region_equal_raw_ecological_EQR_correlations_with_model_point_NP":reg_loci,
            "no_2016_measured_WWTP_load_or_causal_claim":True
        }
    return {
        "schema":"gesui_EEA2016_JRC2014to18_GREEN_pointNP_79NUTS2_2016realBQE_independent_offline_region_country_influence_QA_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_first_original_2016_EEA_BQE_GISCO_JRC_NP_SHA256":sha(root/SOURCE),
        "original_2016_EEA_NUTS2_disjoint_preoutcome_fold_SHA256":sha(root/FOLDS),
        "original_2016_actual_BQE_EQR_2014_18_JRC_source_result_SHA256":sha(root/RESULT),
        "source_results_AUDIT_independently_from_archived_region_OOS_error_data":outcome,
        "source_only_original_2016_BQE_correct_year_and_2014to18_JRC_modelled_source":True,
        "never_call_2016_modelled_POINT_NP_a_2016_actual_WWTP_effluent":True,
        "new_effects_are_posthoc_predictive_and_observational_not_causal":True,
    }

def main():
    a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True);p=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    r=check(root)
    p.out.parent.mkdir(parents=True,exist_ok=True)
    p.out.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print(json.dumps({bio:{
        "n":res["original_2016_EEA_numeric_BQE_sites"],
        "regions":res["original_2016_EEA_benthic_NUTS2_region_count"],
        "countries_positive_N":res["N_predictor_country_with_positive_gain_count"],
        "countries_positive_P":res["P_predictor_country_with_positive_gain_count"],
        "source_NP_collinearity_rho":res["regional_original_2014to18_JRC_modelled_PsN_vs_PsP_spearman_collinearity"],
        "folds":res["five_frozen_NUTS2_fold_honest_source_prediction_gains"]
    } for bio,res in r["source_results_AUDIT_independently_from_archived_region_OOS_error_data"].items()},indent=2))

if __name__=="__main__":main()
