#!/usr/bin/env python3
"""Cross-layer residual co-response diagnostics on the frozen gesui H7 panel.

This is observational, exploratory, and NOT a fitted ecological mediation SEM.
Model W: static, modelled wastewater source at H7 outlet; EQR != emergence
biomass; GBIF plant S50 != floristic survey; butterfly HNS != abundance.

Usage:
 python analysis/20261115_crosslayer_residuals/recompute.py \
   --panel results/20261115_sync/hns_macro/H7_454_HNS_panel.csv \
   --output results/20261115_crosslayer_residuals/recomputed.json
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

PANEL_SHA256 = "1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
COMMON = ["crop", "urban", "riverflow", "naturalout", "area",
          "nEEA", "nAnimal", "nPlant", "nButterfly", "W"]
HNS_EXTRA = ["primary_record_coverage_raw", "primary_host_completeness",
             "primary_rf_Lycaenidae", "primary_rf_Pieridae",
             "primary_rf_Papilionidae", "primary_rf_Hesperiidae"]
OUTCOMES = {
    "EQR": ("EQR_original_site_equal_basin_period_mean", COMMON),
    "plant_S50": ("plant_S50", COMMON),
    "HNS_record": ("primary_record_mean", COMMON + HNS_EXTRA),
    "HNS_incidence": ("primary_incidence_mean", COMMON + HNS_EXTRA),
}
PAIRINGS = [("EQR", "HNS_record"), ("plant_S50", "HNS_record"),
            ("EQR", "plant_S50"), ("EQR", "HNS_incidence")]


def source_gate(path: Path) -> pd.DataFrame:
    data = path.read_bytes()
    sha = hashlib.sha256(data).hexdigest()
    if sha != PANEL_SHA256:
        raise ValueError(f"Source SHA mismatch: {sha}")
    df = pd.read_csv(path, dtype={"h": str, "h6": str})
    if len(df) != 454 or df.h.nunique() != 232 or df.h6.nunique() != 90:
        raise ValueError("Original 454/232/90 sample changed")
    if df.duplicated(["h", "period"]).any() or df.groupby("h").h6.nunique().max() != 1:
        raise ValueError("Wrong source grain / H6 lineage")
    if set(df.period) != {"2017-2019", "2020-2022"}:
        raise ValueError("Original source years changed")
    needed = set(COMMON + HNS_EXTRA + [v[0] for v in OUTCOMES.values()])
    if not needed.issubset(df.columns):
        raise ValueError(f"Missing variables: {needed - set(df.columns)}")
    if not np.isfinite(df[list(needed)].to_numpy(dtype=float)).all():
        raise ValueError("Missing or infinite original data")
    return df


def residualize(df: pd.DataFrame, outcome: str, numeric: list[str]):
    """Exactly H6-equal within-sample weighted least squares.

    Numeric features are z-standardized on ALL 454 original rows. Categorical
    fixed effects are one full-rank country-by-EQR-method stratum and period.
    This reproduces the design used in original gesui multilayer diagnostics.
    """
    x = [np.ones(len(df))]
    for name in numeric:
        value = df[name].to_numpy(dtype=float)
        sd = float(value.std(ddof=0))
        x.append((value - value.mean()) / (sd if sd > 1e-10 else 1.0))
    strata = (df.country.fillna("UNKNOWN").astype(str) + "::" +
              df.method.fillna("UNKNOWN").astype(str)).to_numpy()
    for level in sorted(set(strata))[1:]:
        x.append((strata == level).astype(float))
    periods = df.period.astype(str).to_numpy()
    for level in sorted(set(periods))[1:]:
        x.append((periods == level).astype(float))
    X = np.column_stack(x)
    rank = int(np.linalg.matrix_rank(X))
    if rank != X.shape[1]:
        raise ValueError(f"Collinear design: rank {rank}/{X.shape[1]}")
    count = df.h6.value_counts().to_dict()
    weight = df.h6.map(lambda h: 1.0 / count[h]).to_numpy(dtype=float)
    sw = np.sqrt(weight)
    y = df[outcome].to_numpy(dtype=float)
    beta = np.linalg.lstsq(sw[:, None] * X, sw * y, rcond=None)[0]
    return y - X @ beta, {"columns": int(X.shape[1]), "rank": rank,
                         "W_coefficient_per_full_panel_SD": float(beta[numeric.index("W") + 1])}


def group_means(df: pd.DataFrame, a: np.ndarray, b: np.ndarray, mask: np.ndarray):
    sub = pd.DataFrame({
        "h6": df.loc[mask, "h6"].astype(str).to_numpy(),
        "a": a[mask], "b": b[mask]
    })
    means = sub.groupby("h6", sort=True)[["a", "b"]].mean()
    if len(means) < 3:
        raise ValueError("Too few H6 groups")
    r = float(np.corrcoef(means.a.to_numpy(), means.b.to_numpy())[0, 1])
    return r, means.to_numpy(dtype=float), len(means)


def bootstrap_fixed_residuals(means: np.ndarray, repetitions: int = 1000):
    """Cluster-resample H6 GROUP MEAN fixed residuals; no model refitting.

    Reproduces the independent 2026-11-15 JS diagnostic LCG exactly.
    This interval omits coefficient-estimation and model-selection uncertainty.
    """
    state = 20261115
    n = len(means)
    values = []
    for _ in range(repetitions):
        selected = []
        for _j in range(n):
            state = (1664525 * state + 1013904223) & 0xFFFFFFFF
            idx = int((state / 4294967296.0) * n)
            selected.append(idx)
        a, b = means[np.asarray(selected)].T
        values.append(float(np.corrcoef(a, b)[0, 1]))
    values.sort()
    return [values[25], values[975]]


def row_H6_equal_r(df: pd.DataFrame, a: np.ndarray, b: np.ndarray, mask: np.ndarray):
    subset = df.loc[mask]
    counts = subset.h6.value_counts().to_dict()
    weights = subset.h6.map(lambda k: 1.0 / counts[k]).to_numpy(dtype=float)
    aa, bb = a[mask], b[mask]
    wa = np.average(aa, weights=weights)
    wb = np.average(bb, weights=weights)
    return float(np.sum(weights * (aa-wa) * (bb-wb)) /
                 np.sqrt(np.sum(weights * (aa-wa)**2) *
                         np.sum(weights * (bb-wb)**2)))


def calculate(path: Path) -> dict:
    df = source_gate(path)
    residuals = {}
    info = {}
    for label, (column, numbers) in OUTCOMES.items():
        residuals[label], info[label] = residualize(df, column, numbers)
    out = {
        "source_sha256": PANEL_SHA256,
        "sample": {"H7_period": 454, "H7": 232, "H6": 90, "countries": 15},
        "design": info,
        "exposure": "modelled static basin-outlet wastewater geography W, NOT measured effluent at the receptor",
        "status": "exploratory, posthoc, descriptive fixed-residual bootstrap; no multiplicity correction",
        "causal_mediation_or_SEM_fitted": False,
        "pairs": {}
    }
    all_mask = np.ones(len(df), dtype=bool)
    for label1, label2 in PAIRINGS:
        a, b = residuals[label1], residuals[label2]
        r, gmeans, n = group_means(df, a, b, all_mask)
        rec = {
            "H6_group_mean_Pearson_r": r,
            "H6_fixed_residual_resample95": bootstrap_fixed_residuals(gmeans),
            "H6_equal_row_residual_Pearson_r": row_H6_equal_r(df, a, b, all_mask),
            "H6": n, "rows": len(df), "sensitivity": {}
        }
        for name, mask in [
            ("France", (df.country == "FR").to_numpy(dtype=bool)),
            ("outside_France", (df.country != "FR").to_numpy(dtype=bool)),
        ]:
            rr, _, nn = group_means(df, a, b, mask)
            rec["sensitivity"][name] = {
                "rows": int(mask.sum()), "H6": nn, "H6_group_mean_r": rr,
                "H6_equal_row_r": row_H6_equal_r(df, a, b, mask)
            }
        out["pairs"][label1+"__"+label2] = rec
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = calculate(args.panel)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, ensure_ascii=False, indent=2, sort_keys=True))
    print(json.dumps({
        "source_sha256": output["source_sha256"],
        "pairs": {p: {"r": v["H6_group_mean_Pearson_r"],
                       "95": v["H6_fixed_residual_resample95"]}
                  for p, v in output["pairs"].items()}
    }, indent=2))


if __name__ == "__main__":
    main()
