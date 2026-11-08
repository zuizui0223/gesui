#!/usr/bin/env python3
"""Posthoc cross-year transfer audit restricted to country × EQR-method compatibility.

The original 2018 and 2022 WWTP source-year biological EQR datasets, source
checksums, model specification and site-disjoint testing remain unchanged.
This diagnostic was planned AFTER the positive year-specific/transport results.
It is not a causal inference or a new independent confirmation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

import audit_2018_2022_site_disjoint_temporal_model_transfer as frozen

MODEL_CHOICES = {
    "M0_baseline": (False, False),
    "M1_upstream_reported_mass": (True, False),
    "M2_downstream_negative_control": (False, True),
}
N_BOOT = 1999
BOOT_SEED = 20261108


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_locked_data(folder: Path, year: str, meta: dict) -> pd.DataFrame:
    path = folder / meta[year]
    actual_sha = digest(path)
    if actual_sha != meta[year + "sha"]:
        raise ValueError(f"Original EEA site/EQR CSV SHA changed: {path}")
    df = pd.read_csv(path, dtype={"site_key": str}, low_memory=False)
    if len(df) != meta[year + "n"] or df.site_key.duplicated().any():
        raise ValueError(f"Original EEA site population/identity changed: {path}")
    return df


def predict_error_vectors(train: pd.DataFrame, test: pd.DataFrame,
                          source_year: str, target_year: str,
                          nutrient: str) -> dict[str, np.ndarray]:
    response_train = train[frozen.YEAR_RESPONSE[source_year]].to_numpy(float)
    response_test = test[frozen.YEAR_RESPONSE[target_year]].to_numpy(float)
    if not np.isfinite(response_train).all() or not np.isfinite(response_test).all():
        raise ValueError("Invalid EQR biological outcome")
    errors = {}
    for name, (upstream, downstream) in MODEL_CHOICES.items():
        xtrain, nums = frozen.matrix(train, source_year, nutrient,
                                     upstream, downstream)
        xtest, nums_test = frozen.matrix(test, target_year, nutrient,
                                         upstream, downstream)
        if nums != nums_test:
            raise AssertionError("Original feature specification drift")
        model = Pipeline([
            ("source_normalizer", ColumnTransformer([
                ("numeric", StandardScaler(), nums),
                ("categorical", OneHotEncoder(handle_unknown="ignore"),
                 ["country", "country_method"]),
            ], remainder="drop")),
            ("ridge_fixed_10", Ridge(alpha=frozen.ALPHA)),
        ])
        model.fit(xtrain, response_train)
        errors[name] = (response_test - model.predict(xtest)) ** 2
    return errors


def summarise(errors: dict[str, np.ndarray], subset: np.ndarray,
              country: np.ndarray, source_positive: np.ndarray) -> dict:
    n = int(subset.sum())
    if not n:
        return {"n": 0, "source_positive_n": 0, "countries": 0,
                "mse": None, "upstream_gain_pct": None,
                "downstream_placebo_gain_pct": None,
                "upstream_better_than_baseline_and_placebo": None}
    mse = {name: float(values[subset].mean()) for name, values in errors.items()}
    m0, m1, m2 = (mse[x] for x in MODEL_CHOICES)
    gain = 100 * (m0 - m1) / m0
    placebo = 100 * (m0 - m2) / m0
    return {
        "n": n, "source_positive_n": int(np.sum(subset & source_positive)),
        "countries": int(np.unique(country[subset]).size),
        "mse": mse, "upstream_gain_pct": float(gain),
        "downstream_placebo_gain_pct": float(placebo),
        "upstream_better_than_baseline_and_placebo": bool(m1 < m0 and m1 < m2),
    }


def country_diagnostic(errors: dict[str, np.ndarray],
                       subset: np.ndarray, country: np.ndarray) -> dict:
    labels = sorted(set(country[subset]))
    if not labels:
        return {"countries": 0, "country_rows": [],
                "country_balanced_gain_pct": None,
                "descriptive_country_bootstrap_gain_95pct": None}
    rows = []
    for label in labels:
        mask = subset & (country == label)
        q = summarise(errors, mask, country, np.zeros(len(country), dtype=bool))
        rows.append({
            "country": str(label), "n": q["n"],
            "baseline_mse": q["mse"]["M0_baseline"],
            "upstream_mse": q["mse"]["M1_upstream_reported_mass"],
            "placebo_mse": q["mse"]["M2_downstream_negative_control"],
            "upstream_gain_pct": q["upstream_gain_pct"],
            "upstream_outperforms_placebo": bool(
                q["mse"]["M1_upstream_reported_mass"]
                < q["mse"]["M2_downstream_negative_control"]),
        })
    balanced = {
        name: float(np.mean([r[k] for r in rows]))
        for name, k in [
            ("M0", "baseline_mse"), ("M1", "upstream_mse"),
            ("M2", "placebo_mse")
        ]
    }
    balanced_gain = 100 * (balanced["M0"] - balanced["M1"]) / balanced["M0"]

    # Cluster resampling uses countries, NEVER individual monitoring stations.
    country_n = np.array([r["n"] for r in rows], dtype=float)
    totals = {
        "M0": country_n * np.array([r["baseline_mse"] for r in rows]),
        "M1": country_n * np.array([r["upstream_mse"] for r in rows]),
        "M2": country_n * np.array([r["placebo_mse"] for r in rows]),
    }
    rng = np.random.default_rng(BOOT_SEED)
    draws = rng.integers(0, len(rows), size=(N_BOOT, len(rows)))
    drawn_n = country_n[draws].sum(axis=1)
    m0 = totals["M0"][draws].sum(axis=1) / drawn_n
    m1 = totals["M1"][draws].sum(axis=1) / drawn_n
    m2 = totals["M2"][draws].sum(axis=1) / drawn_n
    boot_gain = 100 * (m0 - m1) / m0
    boot_direction = 100 * (m2 - m1) / m0
    return {
        "countries": len(rows), "country_rows": rows,
        "countries_with_positive_upstream_gain": sum(
            r["upstream_gain_pct"] > 0 for r in rows),
        "countries_upstream_outperforms_downstream": sum(
            r["upstream_outperforms_placebo"] for r in rows),
        "country_balanced_gain_pct": float(balanced_gain),
        "country_balanced_mse": balanced,
        "descriptive_country_bootstrap_replicates": N_BOOT,
        "descriptive_country_bootstrap_gain_95pct": [
            float(v) for v in np.quantile(boot_gain, [0.025, 0.975])],
        "descriptive_country_bootstrap_upstream_over_placebo_95pp": [
            float(v) for v in np.quantile(boot_direction, [0.025, 0.975])],
        "descriptive_bootstrap_probability_positive_gain": float(
            np.mean(boot_gain > 0)),
    }


def audit_one(train: pd.DataFrame, test: pd.DataFrame, source_year: str,
              target_year: str, nutrient: str, published: dict) -> dict:
    if set(train.site_key) & set(test.site_key):
        raise AssertionError("Training and held-out monitoring site overlap")
    ctrain = train.countryCode.fillna("UNKNOWN").astype(str)
    ctest = test.countryCode.fillna("UNKNOWN").astype(str)
    key_train = ctrain + "|" + train.procedureClassificationSystem.fillna(
        "UNKNOWN").astype(str)
    key_test = ctest + "|" + test.procedureClassificationSystem.fillna(
        "UNKNOWN").astype(str)
    known = key_test.isin(set(key_train)).to_numpy(bool)
    unknown = ~known
    country = ctest.to_numpy(str)
    source_positive = (
        test[frozen.EXPOSURES[target_year][nutrient][0]].to_numpy(float) > 0)
    err = predict_error_vectors(train, test, source_year, target_year, nutrient)

    # Independent, case-for-case equivalence to the already-committed
    # unstratified site-disjoint transport. No model-selection permitted.
    for name in MODEL_CHOICES:
        replay = float(err[name].mean())
        locked = published["MSE_and_MAE_by_model"][name]["test_MSE"]
        if not np.isclose(replay, locked, rtol=1e-7, atol=1e-9):
            raise AssertionError(f"Original transport MSE not reproduced: {name}")

    mask_all = np.ones(len(test), dtype=bool)
    result = {
        "training_year": source_year, "testing_year": target_year,
        "site_disjoint": True,
        "original_MSE_reproduced": True,
        "known_method_definition": "test countryCode|procedureClassificationSystem present in training",
        "all_sites": summarise(err, mask_all, country, source_positive),
        "known_method": summarise(err, known, country, source_positive),
        "unknown_method": summarise(err, unknown, country, source_positive),
        "known_method_source_positive": summarise(
            err, known & source_positive, country, source_positive),
        "known_method_source_zero": summarise(
            err, known & ~source_positive, country, source_positive),
        "known_method_country_analysis": country_diagnostic(err, known, country),
        "unknown_method_country_analysis": country_diagnostic(err, unknown, country),
        "known_test_country_method_share": float(np.mean(known)),
        "unchanged_source_and_EQR_years": True,
        "posthoc_not_confirmatory": True,
    }
    assert result["known_method"]["n"] + result["unknown_method"]["n"] == len(test)
    return result


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--year2018", required=True, type=Path)
    p.add_argument("--year2022", required=True, type=Path)
    p.add_argument("--original-transport", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    published = json.loads(args.original_transport.read_text())
    if published.get("scientific_status") != (
            "POSTHOC_TRANSFER_STRESS_TEST__NOT_NEW_PREREGISTERED_REPLICATION"):
        raise ValueError("Unrecognized original transport provenance")
    if published.get("source_model_Ridge_alpha_fixed") != 10.0:
        raise ValueError("Unexpected historical source model")

    groups = {}
    for group, meta in frozen.FILES.items():
        df18 = load_locked_data(args.year2018, "2018", meta)
        df22 = load_locked_data(args.year2022, "2022", meta)
        kept22 = df22[~df22.site_key.isin(set(df18.site_key))].copy()
        kept18 = df18[~df18.site_key.isin(set(df22.site_key))].copy()
        choices = [
            ("2018_training_to_2022_NEW_SITES", df18, kept22, "2018", "2022"),
            ("2022_training_to_2018_NEW_SITES", df22, kept18, "2022", "2018"),
        ]
        groups[group] = {}
        for key, training, testing, source_year, target_year in choices:
            groups[group][key] = {}
            for nutrient in ("P", "N"):
                q = audit_one(training, testing, source_year, target_year,
                              nutrient, published["results"][group][key][nutrient])
                groups[group][key][nutrient] = q
                print("GESUI_COMPATIBILITY_TRANSPORT", json.dumps({
                    "biology": group, "direction": key, "nutrient": nutrient,
                    "known_n": q["known_method"]["n"],
                    "known_source_n": q["known_method"]["source_positive_n"],
                    "known_gain_pct": q["known_method"]["upstream_gain_pct"],
                    "known_placebo_gain_pct": q["known_method"][
                        "downstream_placebo_gain_pct"],
                    "known_countries_positive": q[
                        "known_method_country_analysis"][
                            "countries_with_positive_upstream_gain"],
                }, sort_keys=True), flush=True)

    outcome = {
        "schema": "gesui_original_source_locked_country_method_compatible_transport_posthoc_v01",
        "utc_created": datetime.now(timezone.utc).isoformat(),
        "posthoc": True, "not_sewage_causality": True,
        "not_nutrient_specific": True,
        "bootstrap_country_cluster_descriptive_only": True,
        "original_two_year_site_disjoint_transport_sha256": digest(
            args.original_transport),
        "source_original_EQR_sha256": {
            f"{g}_{yr}": meta[yr + "sha"]
            for g, meta in frozen.FILES.items()
            for yr in ("2018", "2022")
        },
        "results": groups,
        "limitations": [
            "Designed after original two-year positive EQR prediction",
            "EQR assessment methods and calibration are national",
            "Unreported WWTP source is not proof of absence of sewage",
            "Country-cluster bootstrap is posthoc/descriptive and not causal",
            "WWTP P and N source masses strongly collinear",
            "HydroRIVERS topology, agricultural/urban runoff, river dilution and reporting completeness not causally validated"
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(outcome, ensure_ascii=False,
                                   sort_keys=True, indent=2) + "\n")
    print("GESUI_METHOD_COMPATIBLE_TEMPORAL_AUDIT_COMPLETE", flush=True)


if __name__ == "__main__":
    main()
