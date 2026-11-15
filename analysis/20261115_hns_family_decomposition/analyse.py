#!/usr/bin/env python3
"""Same-H7 butterfly host-N affinity (HNS): between- vs within-family composition.

Posthoc SOURCE-LOCKED observational decomposition. Never label these source
record scores butterfly abundance, actual feeding, wastewater dose or mediation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

SOURCE_SHA256 = "1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
COMMON = ["crop", "urban", "riverflow", "naturalout", "area",
          "nEEA", "nAnimal", "nPlant", "nButterfly",
          "primary_record_coverage_raw", "primary_host_completeness"]
FAMILY = ["primary_rf_Lycaenidae", "primary_rf_Pieridae",
          "primary_rf_Papilionidae", "primary_rf_Hesperiidae"]
OUTCOMES = {"record_weighted_HNS": "primary_record_mean",
            "between_family_HNS": "primary_between_family",
            "within_family_HNS": "primary_within_family",
            "incidence_equal_HNS": "primary_incidence_mean",
            "record_minus_incidence_HNS": "_record_minus_incidence"}


def verify_source(path: Path) -> pd.DataFrame:
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    if sha != SOURCE_SHA256:
        raise ValueError(f"Original 454-row panel SHA mismatch: {sha}")
    df = pd.read_csv(path, dtype={"h": str, "h6": str})
    if (len(df), df.h.nunique(), df.h6.nunique(), df.country.nunique()) != (454, 232, 90, 15):
        raise ValueError("Original source grain/count drift")
    if df.duplicated(["h", "period"]).any():
        raise ValueError("Repeated H7 x period")
    if df.groupby("h").h6.nunique().max() != 1:
        raise ValueError("H7 switches H6 parent")
    if set(df.period) != {"2017-2019", "2020-2022"}:
        raise ValueError("Original periods changed")
    err = (df.primary_record_mean -
           df.primary_between_family -
           df.primary_within_family).abs().max()
    if err > 1e-10:
        raise ValueError(f"Original additive identity failed: {err}")
    df["_record_minus_incidence"] = df.primary_record_mean - df.primary_incidence_mean
    df["country_method"] = (df.country.fillna("__UNKNOWN__").astype(str) + "::" +
                            df.method.fillna("__UNKNOWN__").astype(str))
    return df


def design(df: pd.DataFrame, include_family: bool):
    num = COMMON + (FAMILY if include_family else []) + ["W"]
    cols = [np.ones(len(df))]
    names = ["intercept"]
    for k in num:
        vals = df[k].to_numpy(dtype=float)
        if not np.isfinite(vals).all():
            raise ValueError(f"Missing numeric source: {k}")
        sd = vals.std(ddof=0)
        cols.append((vals - vals.mean()) / (sd if sd > 1e-12 else 1.0))
        names.append(k)
    for k in ["country_method", "period"]:
        categories = sorted(df[k].astype(str).unique().tolist())
        for v in categories[1:]:
            cols.append((df[k].astype(str) == v).to_numpy(dtype=float))
            names.append(f"{k}={v}")
    X = np.column_stack(cols)
    rank = int(np.linalg.matrix_rank(X))
    if rank != X.shape[1]:
        raise ValueError(f"Model design not full rank: {rank}/{X.shape[1]}")
    return X, names


def estimate(df: pd.DataFrame, include_family: bool):
    X, names = design(df, include_family)
    group = df.h6.astype(str)
    group_sizes = group.map(group.value_counts()).to_numpy(dtype=float)
    weights = len(df) / (group.nunique() * group_sizes)
    j = names.index("W")
    estimates = {}
    for label, column in OUTCOMES.items():
        y = df[column].to_numpy(dtype=float)
        if not np.isfinite(y).all():
            raise ValueError(f"Missing response: {column}")
        fit = sm.WLS(y, X, weights=weights).fit(
            cov_type="cluster",
            cov_kwds={"groups": group.to_numpy(), "use_correction": True},
            use_t=False)
        b = float(fit.params[j])
        se = float(fit.bse[j])
        estimates[label] = {
            "beta_W_per_1sd": b, "H6_cluster_se": se,
            "H6_cluster_descriptive95": [b - 1.96 * se, b + 1.96 * se],
            "n_rows": len(df), "n_H6": int(group.nunique())
        }
    a = estimates["record_weighted_HNS"]["beta_W_per_1sd"]
    b = estimates["between_family_HNS"]["beta_W_per_1sd"]
    c = estimates["within_family_HNS"]["beta_W_per_1sd"]
    if abs(a - b - c) > 1e-9:
        raise ValueError(f"Linear decomposition coefficient identity failed: {a},{b},{c}")
    if include_family:
        if abs(b) > 1e-9:
            raise ValueError("Between-family score must be exactly explained by family shares")
        if abs(a - 0.144823738741) > 0.000003:
            raise ValueError("Fails independent PR88 source/model baseline coefficient cross-check")
    return {
        "model": "W + country×method + period + common effort/habitat, " +
                 ("plus 4 family shares" if include_family else "without family shares"),
        "n_design_columns": len(names),
        "design_rank": X.shape[1],
        "family_adjusted_between_component_is_algebraically_zero": include_family,
        "estimates": estimates,
        "point_ratio_within_to_total": float(c / a) if abs(a) > 1e-12 else None
    }


def run(panel: Path, output: Path):
    df = verify_source(panel)
    no_family = estimate(df, False)
    with_family = estimate(df, True)
    result = {
        "schema": "GESUI_20261115_H7_HNS_BETWEEN_WITHIN_SOURCE_DECOMPOSITION_V1",
        "source_sha256": SOURCE_SHA256,
        "n_H7_period": 454, "n_H7": 232, "n_H6": 90, "n_countries": 15,
        "design_a_without_family_shares": no_family,
        "design_b_with_four_family_shares": with_family,
        "posthoc_exploratory": True,
        "causal_SEM_or_mediation_fitted": False,
        "W_is_modelled_source_geography_not_measured_wastewater_dose": True,
        "record_weighted_GBIF_source_is_not_butterfly_abundance": True,
        "within_family_is_not_observed_larval_host_feeding": True,
        "do_not_interpret_between_family_adjusted_zero_as_ecological_null": True,
        "no_added_independent_ecology_sites": True,
        "H6_cluster_CI_descriptive_unadjusted_multiple_tests": True
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "SOURCE_LOCKED_HNS_FAMILY_DECOMPOSITION.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    report = []
    for model_key in ["design_a_without_family_shares", "design_b_with_four_family_shares"]:
        model = result[model_key]
        for name, estimate_ in model["estimates"].items():
            report.append({"design": model_key, "response": name,
                           "beta_W": estimate_["beta_W_per_1sd"],
                           "lo95": estimate_["H6_cluster_descriptive95"][0],
                           "hi95": estimate_["H6_cluster_descriptive95"][1]})
    pd.DataFrame(report).to_csv(output / "SOURCE_LOCKED_HNS_FAMILY_DECOMPOSITION.csv", index=False)
    print(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False))
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", type=Path, default=Path(
        "results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    ap.add_argument("--outdir", type=Path, default=Path(
        "results/20261115_hns_family_decomposition"))
    args = ap.parse_args()
    run(args.panel, args.outdir)
