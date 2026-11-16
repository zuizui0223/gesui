#!/usr/bin/env python3
"""Final adversarial EU454 tests: fixed-record plant rarefaction and urban-vs-W.

Original cohort, original trait scores, original 5-fold TRUE-H6 hash, original
source GBIF facets and original 4011 EIV keys ONLY. No network requery and no
local wastewater or observed larval feeding claims. This is a POSTHOC
sensitivity analysis; Monte Carlo quantiles are NOT independent confidence
intervals or evidence of plant-to-butterfly causal mediation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "20261115_eu454_plant_nutri"))
from adversarial_host_exclusion_country import prepare, source_records
from h6_model_if_gate import (
    B_GEOG, P_GEOG, SEED, N_BOOT, fold_index, predict_ridge,
    source_cluster_regression, cv_h6,
)

RSEED = 20261116
DRAWS = 120
DEPTHS = (20, 50, 100)
PRIMARY_DEPTH = 50
MIN_H7_PERIOD = 100
MIN_H6 = 40


def number(v):
    v = float(v)
    if not np.isfinite(v):
        raise RuntimeError("Nonfinite original source numerical result")
    return v


def combine_source(df, facets, numeric):
    """Retrieve scored original occurrence-count vectors; no zero imputation."""
    sites = []
    for row in df.itertuples(index=False):
        facet = facets.get((str(row.h), str(row.period)))
        if facet is None or facet.get("source_status") != "ok":
            raise RuntimeError("Missing original source facet: "+repr((row.h, row.period)))
        counts = facet.get("source_speciesKey_counts")
        if not isinstance(counts, dict):
            raise RuntimeError("Original species-key counts missing")
        scored = sorted((str(k), int(v)) for k, v in counts.items()
                        if str(k) in numeric and int(v) > 0)
        if len(scored) < 10:
            raise RuntimeError("Existing original scored-species hard gate broken")
        keys, amounts = zip(*scored)
        weights = np.asarray(amounts, dtype=np.int64)
        vals = np.asarray([numeric[k] for k in keys], dtype=np.float64)
        if not np.isclose(vals.mean(), float(row.Plant_N_incidence), rtol=0, atol=1e-9):
            raise RuntimeError("Rarefaction source does not exactly reproduce original species-equal score")
        if np.any(weights <= 0) or not np.isfinite(vals).all():
            raise RuntimeError("Original records malformed")
        sites.append((weights, vals, int(weights.sum())))
    if len(sites) != 451:
        raise RuntimeError("Original source cohort drift")
    return sites


def fixed_errors(df, target, models):
    """Same 5 original genuine-H6 folds for every model, no leakage."""
    names = sorted(models)
    h6 = df.h6.astype(str).to_numpy()
    y = df[target].to_numpy(float)
    folds = np.array([fold_index(g) for g in h6], int)
    predictions = {name: np.full(len(df), np.nan) for name in names}
    for fold in range(5):
        ix = np.flatnonzero(folds == fold)
        jj = np.flatnonzero(folds != fold)
        if len(ix) < 1 or len(jj) < 30 or set(h6[ix]) & set(h6[jj]):
            raise RuntimeError("Invalid original H6 fold split")
        train = df.iloc[jj].reset_index(drop=True)
        test = df.iloc[ix].reset_index(drop=True)
        for name in names:
            predictions[name][ix] = predict_ridge(train, test, target, models[name])
    by = {}
    for name in names:
        if not np.isfinite(predictions[name]).all():
            raise RuntimeError("Missing genuine H6 holdout predictions")
        error = (y-predictions[name])**2
        by[name] = np.asarray(
            [error[h6 == g].mean() for g in sorted(set(h6))], dtype=float
        )
    if not all(np.isfinite(x).all() and x.mean() > 1e-12 for x in by.values()):
        raise RuntimeError("Invalid source equal-H6 MSE")
    return by


def contrast(err0, err1, rng):
    a = np.asarray(err0, dtype=float)
    b = np.asarray(err1, dtype=float)
    if a.shape != b.shape:
        raise RuntimeError("Different geographic holdout units")
    pct = number(100*(a.mean()-b.mean())/a.mean())
    draws = rng.integers(0, len(a), size=(N_BOOT,len(a)))
    aa = a[draws].mean(axis=1)
    bb = b[draws].mean(axis=1)
    boot = 100*(aa-bb)/aa
    return {
        "H6_equal_out_of_fold_MSE_improvement_percent": pct,
        "same_fixed_H6_prediction_errors_bootstrap95": [number(x) for x in np.quantile(boot,[.025,.975])],
        "n_genuine_H6": len(a),
        "not_model_refit_bootstrap_or_causal_effect": True,
    }


def urban_vs_w(df, target, base):
    if "urban" not in base or "W" in base:
        raise RuntimeError("Incorrect original baseline model specification")
    base = [x for x in base if x != "urban"]
    models = {
        "neither": base,
        "W_only": base+["W"],
        "urban_only": base+["urban"],
        "both": base+["W","urban"],
    }
    errors = fixed_errors(df, target, models)
    rng = np.random.default_rng(RSEED)
    return {
        "target": target,
        "original_source_cohort": {"n_H7period": len(df), "n_H6": int(df.h6.nunique())},
        "original_constant_covariates": base,
        "design_models": models,
        "source_W_urban_pearson_r": number(df[["W","urban"]].corr().iloc[0,1]),
        "baseline_to_W": contrast(errors["neither"], errors["W_only"], rng),
        "baseline_to_urban": contrast(errors["neither"], errors["urban_only"], rng),
        "W_increment_after_urban": contrast(errors["urban_only"], errors["both"], rng),
        "urban_increment_after_W": contrast(errors["W_only"], errors["both"], rng),
        "W_joint_H6_adjusted": source_cluster_regression(df,target,base+["urban","W"],"W"),
        "urban_joint_H6_adjusted": source_cluster_regression(df,target,base+["W","urban"],"urban"),
        "observational_only": True,
    }


def ra_refaction(df, sites, depth, n_draws=DRAWS):
    eligible = np.array([total >= depth for _,_,total in sites], dtype=bool)
    sub = df.loc[eligible].reset_index(drop=True)
    ess = [x for x, e in zip(sites, eligible) if e]
    status = {
        "depth_scored_original_plant_occurrence_records": depth,
        "n_eligible_H7period": int(len(sub)),
        "n_eligible_H6": int(sub.h6.nunique()),
        "n_original_H7period": len(df),
        "source_support_exclusion_only": True,
    }
    if len(sub) < MIN_H7_PERIOD or sub.h6.nunique() < MIN_H6:
        return {**status, "status": "STOP_SOURCE_OPPORTUNITY",
                "no_outcome_fitted": True}
    # Counterfactual full-depth comparison on the EXACT SAME eligible subcohort.
    full = source_cluster_regression(
        sub, "primary_incidence_mean", B_GEOG+["Plant_N_incidence"], "Plant_N_incidence")
    fullcv = cv_h6(sub, "primary_incidence_mean", B_GEOG, "Plant_N_incidence")
    rng = np.random.default_rng(RSEED + depth)
    results = []
    for i in range(n_draws):
        draw = sub.copy()
        means = []
        for weights, score, total in ess:
            selected = rng.multivariate_hypergeometric(weights, depth)
            occupied = selected > 0
            if not occupied.any() or int(selected.sum()) != depth:
                raise RuntimeError("Fixed original-count rarefaction draw invalid")
            # Preserve the original species-equal estimand on the sampled
            # detected species, not the count-weighted mean.
            means.append(number(score[occupied].mean()))
        draw["rarefied_plant_N"] = np.asarray(means, dtype=float)
        if draw.rarefied_plant_N.std(ddof=0) < 1e-12:
            raise RuntimeError("Degenerate randomized original score")
        beta = source_cluster_regression(
            draw, "primary_incidence_mean", B_GEOG+["rarefied_plant_N"], "rarefied_plant_N")
        cv = cv_h6(draw, "primary_incidence_mean", B_GEOG, "rarefied_plant_N")
        results.append({
            "iteration": i+1,
            "H6_cluster_beta": beta["adjusted_beta_per_1SD_focal"],
            "H6_descriptive95_low": beta["H6_cluster_descriptive95"][0],
            "H6_descriptive95_high": beta["H6_cluster_descriptive95"][1],
            "H6_equal_out_MSE_gain_percent": cv["H6_equal_added_MSE_gain_pct"],
            "CV_fixed_H6_bootstrap_low": cv["fixed_out_of_fold_error_H6_resampling95_pct"][0],
            "correlation_full_N_rarefied_N": number(
                np.corrcoef(draw.Plant_N_incidence, draw.rarefied_plant_N)[0,1]),
        })
    res = pd.DataFrame(results)
    def summary(column):
        a = res[column].to_numpy(float)
        return {"median": number(np.median(a)),
                "p2p5_p97p5_MC_sensitivity": [number(v) for v in np.quantile(a,[.025,.975])],
                "all_draws_positive_fraction": number((a > 0).mean())}
    return {
        **status, "status": "SOURCE_GATE_PASS_POSTHOC_MC_SENSITIVITY",
        "n_source_fixed_record_monte_carlo_draws": n_draws,
        "baseline_original_plant_N_H6_beta": full["adjusted_beta_per_1SD_focal"],
        "baseline_original_plant_N_H6_beta95": full["H6_cluster_descriptive95"],
        "baseline_original_plant_N_H6_MSE_gain_pct": fullcv["H6_equal_added_MSE_gain_pct"],
        "fixed_record_rarefied_beta": summary("H6_cluster_beta"),
        "fixed_record_rarefied_H6_out_gain_pct": summary("H6_equal_out_MSE_gain_percent"),
        "full_vs_rarefied_N_pearson": summary("correlation_full_N_rarefied_N"),
        "fraction_draws_beta_cluster95_lower_above_zero": number((res.H6_descriptive95_low > 0).mean()),
        "fraction_draws_fixed_H6_bootstrap_lower_above_zero": number((res.CV_fixed_H6_bootstrap_low > 0).mean()),
        "source_detection_bias_not_eliminated": True,
        "MC_quantiles_are_not_independent_between-H6_confidence_intervals": True,
    }


def main(args):
    df, hosts, unionkeys, numeric, gate = prepare(args.facets, args.source)
    if len(df) != 451 or df.h6.nunique() != 90 or len(numeric) != 4011:
        raise RuntimeError("Unchanged original 451/90/4011 source freeze violated")
    facets = source_records(args.facets)
    sites = combine_source(df, facets, numeric)
    out = {
        "schema": "GESUI_POSTHOC_FINAL_EU454_FIXED_RECORD_AND_W_URBAN_FALSIFIERS_V1",
        "status": "SOURCE_ORIGINAL_POSTHOC_FALSIFIER",
        "frozen_main_panel_sha256": hashlib.sha256(
            Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv").read_bytes()).hexdigest(),
        "source_original_h6_cohort": int(df.h6.nunique()),
        "source_4011_exact_EIV_scored_accepted_vascular_plant_species": len(numeric),
        "source_occurrence_facets_for_original_451_H7period": True,
        "original_fixed_genuine_H6_hash_5fold_and_ridge_alpha10": True,
        "already_existing_host_exclusion_result_remains_separate": True,
        "urban_vs_W": {
            "plant_NUTRI_primary": urban_vs_w(df,"Plant_N_incidence",P_GEOG),
            "butterfly_HNS_secondary": urban_vs_w(
                df,"primary_incidence_mean", [k for k in B_GEOG if k != "W"]),
        },
        "rarefaction": {
            str(k): ra_refaction(df, sites, k) for k in DEPTHS
        },
        "primary_fixed_scored_record_depth": PRIMARY_DEPTH,
        "original_other_EIV_single_indicator_controls_and_nulls_unchanged": True,
        "noncausal_not_measured_wastewater_N_P": True,
        "posthoc_test_many_existing_analyses_not_confirmatory_independent_replication": True,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output/"GESUI_FINAL_RAREFACTION_W_URBAN_SOURCE_REPORT.json"
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    brief = {"status":out["status"], "original_H6":out["source_original_h6_cohort"],
             "rarefaction":{k:{"status":v["status"],"H6":v["n_eligible_H6"],
                             "beta":v.get("fixed_record_rarefied_beta"),
                             "gain":v.get("fixed_record_rarefied_H6_out_gain_pct")}
                            for k,v in out["rarefaction"].items()},
             "urbanW":{k:{"W_after_urban":v["W_increment_after_urban"],
                          "urban_after_W":v["urban_increment_after_W"]}
                       for k,v in out["urban_vs_W"].items()}}
    print("GESUI_FINAL_FALSIFIER="+json.dumps(brief,ensure_ascii=False,sort_keys=True),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--facets",type=Path,required=True)
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
