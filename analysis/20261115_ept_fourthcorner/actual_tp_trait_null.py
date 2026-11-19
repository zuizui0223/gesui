#!/usr/bin/env python3
"""Exploratory EPA three-cycle *ambient measured TP* vs observed EPT traits.

Extension of the exact source-locked 2026-11-15 fourthcorner cohort. Unlike the
previous XEMBED null, this routine treats log_tp as the focal environment in a
within-EPT-order taxon-trait permutation, and removes log_tp from the nuisance
matrix. It then checks whether TP adds predictive information *after* physical
habitat + TN on survey-cycle- and HUC8-disjoint test data.

Everything is observational; measured river TP is NOT measured WWTP-derived TP.
Taxon incidence is EPA observed TAXA_ID presence, not occupancy-corrected
presence/absence or demonstrated taxon turnover through time.
"""
from __future__ import annotations

import argparse
import collections
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.preprocessing import OneHotEncoder, StandardScaler

import fourthcorner as fc


def source_only_strict(full: pd.DataFrame, matrix: tuple, kind: str, era: str):
    """Posthoc stress subset from prior fourthcorner, retain original row order."""
    mask = ((full["known_taxa_n"] >= 3)
            & (full["trait_taxon_coverage"] >= 0.8)
            & (full["original_EPT_counted_individuals"] >= 10))
    site = full.loc[mask].copy()
    expected = {
        "tolerance": {"2008_09": (352, 296), "2013_14": (679, 483),
                      "2018_19": (567, 397)},
        "clinger": {"2008_09": (377, 313), "2013_14": (767, 526),
                    "2018_19": (624, 423)},
    }
    if (len(site), site.huc8.nunique()) != expected[kind][era]:
        raise ValueError("Prior strict source cohort changed: " + era + kind)

    id_to_idx = {name: i for i, name in enumerate(full.source_UID.astype(str))}
    rows = np.asarray([id_to_idx[x] for x in site.source_UID.astype(str)], int)
    abundance, incidence, values, orders, keys = matrix
    # Restrict shuffling to scored taxa actually present in the strict cohort.
    used = abundance[rows].sum(axis=0) > 0
    a = abundance[rows][:, used]
    b = incidence[rows][:, used]
    v = values[used]
    o = orders[used]
    names = np.asarray(keys, dtype=object)[used]
    if not np.allclose(a @ v, site.CWM_count, rtol=0, atol=1e-10):
        raise ValueError("Strict taxon count matrix CWM mismatch")
    if not np.allclose(b @ v, site.CWM_incidence, rtol=0, atol=1e-10):
        raise ValueError("Strict taxon incidence matrix CWM mismatch")
    return site, (a, b, v, o, names)


def label_null(site: pd.DataFrame, matrix: tuple, target: str, B: int, seed: int):
    """Conditional within-order permutation for measured log TP or TN.

    The same exact full-rank WLS nuisance design as fc.model is used.
    In particular, *focal target must be absent* from the control matrix.
    Per-order fixed taxon-score permutation keeps local abundance, scored
    taxon missingness, site environment and geography unchanged.
    """
    if target not in ("log_tp", "log_tn"):
        raise ValueError("This experiment is only for measured ambient NP")
    abundance, incidence, scores, orders, keys = matrix
    x = site[target].to_numpy(float)
    x = (x - x.mean()) / x.std(ddof=0)
    base = [v for v in fc.GEO + fc.NP + fc.PHAB if v != target]
    assert target not in base
    ctrl = pd.concat([site[base].astype(float),
                      pd.get_dummies(site[fc.CATEGORICAL], drop_first=True,
                                     dtype=float)], axis=1)
    C = sm.add_constant(ctrl.to_numpy(float), has_constant="add")
    if np.linalg.matrix_rank(C) != C.shape[1]:
        raise ValueError("Nuisance model rank deficient")
    group = site.huc8.astype(str).to_numpy()
    freq = collections.Counter(group)
    weights = np.asarray([1.0 / freq[g] for g in group], float)
    sw = np.sqrt(weights)
    predicted = C @ np.linalg.lstsq(sw[:, None] * C, sw * x, rcond=None)[0]
    residual = x - predicted
    den = float(np.dot(weights * residual, residual))
    if den <= 1e-10:
        raise ValueError("Measured ambient target essentially explained by nuisance")
    contributions = {
        "CWM_count": (weights * residual) @ abundance / den,
        "CWM_incidence": (weights * residual) @ incidence / den,
    }
    observed = {m: float(np.dot(weights_by_taxon, scores))
                for m, weights_by_taxon in contributions.items()}
    # Independent regression reproduction is required BEFORE any label p value.
    for metric in contributions:
        ref = fc.model(site, metric)[target]["beta"]
        if not np.isclose(observed[metric], ref, rtol=0, atol=1e-9):
            raise ValueError("Permutation WLS slope disagrees with original regression")
    rng = np.random.default_rng(seed)
    groups = [np.flatnonzero(orders == label) for label in np.unique(orders)]
    draws = {metric: np.zeros(B, dtype=float) for metric in contributions}
    for k in range(B):
        perm = scores.copy()
        for idx in groups:
            if len(idx) > 1:
                perm[idx] = rng.permutation(scores[idx])
        for metric, v in contributions.items():
            draws[metric][k] = float(v @ perm)
    out = {}
    for metric, dist in draws.items():
        center = float(dist.mean())
        obs = observed[metric]
        exceed = int(np.sum(np.abs(dist - center)
                            >= np.abs(obs - center) - 1e-12))
        out[metric] = {
            "observed_beta": obs,
            "within_order_twosided_p": (exceed + 1) / (B + 1),
            "null_center": center,
            "null_percentile95": [float(z) for z in np.quantile(dist, [.025, .975])],
            "scored_taxon_n": int(len(scores)),
            "EPT_order_n": int(len(groups)),
            "B": B,
        }
    return out


def heldout_nutrient_gain(data: dict, source_huc: dict, kind: str,
                          metric: str, target: str, era: str):
    """True original survey-cycle and whole HUC8 nonoverlap; no CV leakage."""
    test = data[(era, kind)].copy()
    train = pd.concat([data[(c, kind)] for c in fc.CYCLES if c != era],
                      ignore_index=True)
    train = train.loc[~train.huc8.isin(source_huc[era])].copy()
    if set(train.huc8) & set(test.huc8):
        raise ValueError("Original HUC8 leakage")
    if len(test) < 250 or test.huc8.nunique() < 150:
        return {"status": "SOURCE_GATE_STOP_TEST"}
    if len(train) < 300 or train.huc8.nunique() < 200:
        return {"status": "SOURCE_GATE_STOP_TRAIN",
                "train_n": len(train), "train_HUC8": train.huc8.nunique()}

    other = "log_tn" if target == "log_tp" else "log_tp"
    baseline = fc.GEO + fc.PHAB + [other]
    extended = baseline + [target]

    def predict(cols):
        processor = ColumnTransformer(
            [("numeric", StandardScaler(), cols),
             ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False),
              fc.CATEGORICAL)],
            remainder="drop", sparse_threshold=0,
        )
        x = processor.fit_transform(train)
        z = processor.transform(test)
        counts = train.huc8.value_counts()
        w = (len(train) / len(counts)) / train.huc8.map(counts)
        return Ridge(alpha=10, fit_intercept=True, solver="svd").fit(
            x, train[metric].to_numpy(float),
            sample_weight=w.to_numpy(float),
        ).predict(z)

    y = test[metric].to_numpy(float)
    base_pred = predict(baseline)
    target_pred = predict(extended)
    huc, inverse = np.unique(test.huc8.astype(str), return_inverse=True)
    huc_n = np.bincount(inverse)
    err0 = np.bincount(inverse, weights=(y - base_pred) ** 2) / huc_n
    err1 = np.bincount(inverse, weights=(y - target_pred) ** 2) / huc_n
    rng = np.random.default_rng(20261119 + fc.CYCLES.index(era)
                                + (0 if target == "log_tp" else 10))
    ix = rng.integers(0, len(huc), size=(2000, len(huc)))
    gains = 100 * (err0[ix].mean(axis=1) - err1[ix].mean(axis=1)) / err0[ix].mean(axis=1)
    return {
        "status": "ERA_AND_HUC8_DISJOINT",
        "target": target,
        "train_n": len(train), "train_HUC8": train.huc8.nunique(),
        "test_n": len(test), "test_HUC8": len(huc),
        "baseline_MSE": float(err0.mean()),
        "extended_MSE": float(err1.mean()),
        "added_MSE_gain_percent": float(100 * (err0.mean() - err1.mean()) / err0.mean()),
        "fixed_holdout_HUC8_2000_bootstrap95": [float(z) for z in np.quantile(gains, [.025, .975])],
    }


def run(B: int, out: Path):
    original, counts, sites, qa = fc.read()
    data, matrices = fc.compute(counts, sites, qa)
    result = {
        "schema": "GESUI_EPA3ERA_AMBIENT_TP_TOLERANCE_COUNT_INCIDENCE_V1",
        "status": "EXPLORATORY_POSTHOC",
        "same_original_EPA_3cycles_2209_source_sites": True,
        "same_original_EPA_source_UID_taxa": True,
        "source_native_original_2648_visits": int(len(original)),
        "measured_ambient_TP_NOT_source_specific_WWTP": True,
        "posthoc_strict_coverage_not_preregistered": True,
        "label_null_preserves_order_not_phylogeny": True,
        "TP_primary_TN_context_control": True,
        "models_adjust_for_TN_TP_and_6_PHab_variables": True,
        "strong_temporal_inference_from_same_sites": False,
        "permutation_B": B,
        "strict_by_era": {},
        "original_by_era": {},
        "external_era_HUC8": {},
    }
    source_huc = {
        era: set(original.loc[original.source_cycle == era, "huc8"])
        for era in fc.CYCLES
    }
    pvals = []
    for era in fc.CYCLES:
        result["original_by_era"][era] = {}
        result["strict_by_era"][era] = {}
        for kind in ("tolerance", "clinger"):
            full = data[(era, kind)]
            mat = matrices[(era, kind)]
            strict, strict_mat = source_only_strict(full, mat, kind, era)
            for cohort, site, L in (
                ("original_by_era", full, mat),
                ("strict_by_era", strict, strict_mat),
            ):
                slopes = {}
                for metric in ("CWM_count", "CWM_incidence", "count_minus_incidence"):
                    mod = fc.model(site, metric)
                    slopes[metric] = {"log_tp": mod["log_tp"],
                                      "log_tn": mod["log_tn"]}
                permutations = {
                    "log_tp": label_null(site, L, "log_tp", B, 20261119 + 31 * fc.CYCLES.index(era)
                                         + (1 if kind == "tolerance" else 2) + (100 if cohort == "strict_by_era" else 0)),
                    "log_tn": label_null(site, L, "log_tn", B, 20261129 + 31 * fc.CYCLES.index(era)
                                         + (1 if kind == "tolerance" else 2) + (100 if cohort == "strict_by_era" else 0)),
                }
                result[cohort][era][kind] = {
                    "n_sites": len(site), "HUC8": int(site.huc8.nunique()),
                    "slope": slopes, "trait_label_null": permutations,
                }
                if cohort == "strict_by_era":
                    for metric in ("CWM_count", "CWM_incidence"):
                        pvals.append((era, kind, metric, mod if False else slopes[metric]["log_tp"]["p_nominal"]))
        result["external_era_HUC8"][era] = {
            kind: {metric: {target: heldout_nutrient_gain(data, source_huc, kind, metric, target, era)
                            for target in ("log_tp", "log_tn")}
                   for metric in ("CWM_count", "CWM_incidence")}
            for kind in ("tolerance", "clinger")
        }

    # Family of twelve TP associations: 2 traits x 2 weighting schemes x 3 eras.
    corrected = fc.bh([p for _, _, _, p in pvals])
    result["TP_strict_BH12"] = [
        {"era": era, "kind": kind, "metric": metric, "p_nominal": p, "BH12_q": q}
        for (era, kind, metric, p), q in zip(pvals, corrected)
    ]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    print("GESUI_AMBIENT_TP_TRAIT_FALSIFIER", json.dumps({
        "strict": {era: {
            kind: {
                "count": result["strict_by_era"][era][kind]["slope"]["CWM_count"]["log_tp"],
                "incidence": result["strict_by_era"][era][kind]["slope"]["CWM_incidence"]["log_tp"],
                "TP_incidence_label_p": result["strict_by_era"][era][kind]["trait_label_null"]["log_tp"]["CWM_incidence"]["within_order_twosided_p"],
            } for kind in ("tolerance", "clinger")}
            for era in fc.CYCLES},
        "external": result["external_era_HUC8"],
        "BH12": result["TP_strict_BH12"],
    }, sort_keys=True, allow_nan=False), flush=True)
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--permutations", type=int, default=999)
    ap.add_argument("--out", type=Path, default=Path(os.environ.get("RUNNER_TEMP", "/tmp"))
                    / "gesui_tp_falsifier" / "AMBIENT_TP_TRAIT_FALSIFIER.json")
    args = ap.parse_args()
    run(args.permutations, args.out)
