#!/usr/bin/env python3
"""POSTHOC geographical scale diagnostic for Lycaenidae photo share vs WWTP P.

Uses already-frozen aggregate public sources only; 2022 EEA P is the SAME spatial
source covariate applied to all three photo years, not a year-varying exposure.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

from gesui_lycaenidae_2021_2023_posthoc_photo_composition_differencing import (
    source_json, sha256, solve, build_rows,
)

SEED = 20261109
BOOT_N = 4000
SHUFFLE_N = 3999
VARIABLES = ("logP", "logN", "logSourceCount")
RESPONSES = ("Lycaenidae2021", "Lycaenidae2022", "Lycaenidae2023", "LycaenidaeMean")
MODELS = {
    "P_only": (),
    "P_plus_land": ("logPlant", "built", "crop", "logArea"),
    "P_plus_land_and_cell_geo": ("logPlant", "built", "crop", "logArea", "lon", "lat"),
}


def lcg_uniform(state):
    state = (1664525 * state + 1013904223) & 0xFFFFFFFF
    return state, state / (2 ** 32)


def contrast(rows, source, response, controls=()):
    """Within 2-degree geographic blocks; singletons do not identify slope."""
    if not rows:
        return None
    x, y = [], []
    for a in rows:
        fields = (source, *controls, response)
        means = {field: sum(r[field] for r in a) / len(a) for field in fields}
        for r in a:
            x.append([r[f] - means[f] for f in (source, *controls)])
            y.append(r[response] - means[response])
    b = solve(x, y)
    return None if b is None else {"slope": b[0], "n": len(x), "predictor_count": len(x[0])}


def between(groups, source, response):
    """Unweighted between-block mean OLS with country, geography and land controls."""
    means = []
    for a in groups:
        fields = (source, response, "logPlant", "built", "crop", "logArea", "lon", "lat")
        m = {f: sum(r[f] for r in a) / len(a) for f in fields}
        m["country"] = a[0]["country"]
        means.append(m)
    x = [
        [1., float(r["country"] == "DE"), r["logPlant"], r["built"],
         r["crop"], r["logArea"], r["lon"], r["lat"], r[source]]
        for r in means
    ]
    b = solve(x, [r[response] for r in means])
    return None if b is None else {"slope": b[-1], "regions": len(means), "predictor_count": len(x[0])}


def pct(v, p):
    v = sorted(v)
    return v[int(math.floor((len(v) - 1) * p))]


def build_inputs(root):
    files = {
        "2021_family": root / "results/lycaenidae_2021_prior_butterfly_photo_source/2021_HISTORICAL_3FAMILY_PHOTO_COUNTS_ORIGINAL43_SOURCE_ONLY.json",
        "2022_family": root / "results/lycaenidae_de_at_2022_macro_source_photo/2022_DE_AT_LYCAENIDAE_PHOTO_SHARE_VS_WWTP_SOURCE_EXPLORATORY.json",
        "2023_family": root / "results/lycaenidae_2023_photo_temporal_source/2023_INDEPENDENT_YEAR_FAMILY_PHOTO_TOTALS_ORIGINAL43_SOURCE_QC.json",
        "2022_land": root / "results/lycaenidae_worldcover2021_source_only/2022_LYCAENIDAE_ESAWC2021_DE_AT_URBAN_CROP_LANDAREA_POSTHOC_BLOCKCV.json",
    }
    src = {k: source_json(path) for k, path in files.items()}
    rows = build_rows(
        src["2021_family"], src["2022_family"], src["2023_family"],
        src["2022_land"], sha256(files["2022_land"])
    )
    index = {
        (r["countryCode"], int(float(r["lat_floor"])), int(float(r["lon_floor"]))): r
        for r in src["2022_land"]["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    }
    for r in rows:
        s = index[(r["country"], int(r["lat"]), int(r["lon"]))]
        if abs(r["logP"] - math.log1p(float(s["registered_P_t_year"]))) > 1e-10:
            raise ValueError("Mismatch in frozen WWTP P source")
        r["logN"] = math.log1p(float(s["registered_N_t_year"]))
        r["logSourceCount"] = math.log1p(float(s["registered_2022_EEA_both_measured_FW_outfalls"]))
        for field in ("logP", "logPlant", "built", "crop", "logArea", "lat", "lon"):
            r[field] = float(r[field])
        for name, val in r["_y"].items():
            r[name] = float(val)
        r["LycaenidaeMean"] = sum(r["Lycaenidae" + str(y)] for y in (2021, 2022, 2023)) / 3.0
    grouped = {}
    for r in rows:
        grouped.setdefault(r["block"], []).append(r)
    groups = list(grouped.values())
    usable = [a for a in groups if len(a) > 1]
    assert len(rows) == 42 and len(groups) == 22 and len(usable) == 11
    assert sum(len(a) for a in usable) == 31
    assert Counter(r["country"] for r in rows) == {"DE": 27, "AT": 15}
    return files, rows, groups, usable


def analysis(root, bootstrap=BOOT_N, shuffles=SHUFFLE_N):
    files, rows, groups, usable = build_inputs(root)
    out = {
        "schema": "gesui_lycaenidae_posthoc_within_vs_between_2deg_P_source_photo_v01",
        "data_provenance_sha256": {k: sha256(p) for k, p in files.items()},
        "posthoc_noncausal_geographical_diagnostic": True,
        "source": "EEA 2022 registered annual reported WWTP P/N mass, not measured butterfly substrate or 2021/2023 exposure",
        "response": "iNaturalist family-photo logit Lycaenidae / (Pieridae + Papilionidae), not butterfly abundance",
        "sample": {
            "country_1deg_cells": len(rows), "country_2deg_groups": len(groups),
            "groups_with_at_least_two_cells": len(usable),
            "within_region_identifying_cells": sum(len(a) for a in usable),
            "singleton_groups_not_estimating_within_slope": len(groups) - len(usable),
            "country_cell_counts": dict(Counter(r["country"] for r in rows)),
            "three_family_photo_years": [2021, 2022, 2023],
        },
        "within": {}, "between": {}, "sensitivity_bootstrap": {},
        "withinP_shuffle_diagnostic": {},
    }
    for source in VARIABLES:
        out["within"][source] = {}
        out["between"][source] = {}
        for response in RESPONSES:
            out["within"][source][response] = {
                name: contrast(usable, source, response, controls)
                for name, controls in MODELS.items()
            }
            out["between"][source][response] = between(groups, source, response)

    state = SEED
    for response in ("Lycaenidae2022", "LycaenidaeMean"):
        out["sensitivity_bootstrap"][response] = {}
        out["withinP_shuffle_diagnostic"][response] = {}
        for name, controls in MODELS.items():
            results = []
            for _ in range(bootstrap):
                sample = []
                for _ in range(len(usable)):
                    state, p = lcg_uniform(state)
                    sample.append(usable[int(p * len(usable))])
                b = contrast(sample, "logP", response, controls)
                if b is not None:
                    results.append(b["slope"])
            orig = out["within"]["logP"][response][name]["slope"]
            out["sensitivity_bootstrap"][response][name] = {
                "draws": bootstrap, "valid": len(results), "seed": SEED,
                "descriptive_cluster95": [pct(results, .025), pct(results, .975)],
                "positive_fraction": sum(v > 0 for v in results) / len(results),
                "effect_matches_original": abs(orig) < 1e20,
            }

            permutation_slopes = []
            for _ in range(shuffles):
                shuffled = []
                for group in usable:
                    cloned = [dict(row) for row in group]
                    ps = [r["logP"] for r in cloned]
                    for j in range(len(ps) - 1, 0, -1):
                        state, u = lcg_uniform(state)
                        k = int(u * (j + 1))
                        ps[j], ps[k] = ps[k], ps[j]
                    for r, v in zip(cloned, ps):
                        r["logP"] = v
                    shuffled.append(cloned)
                b = contrast(shuffled, "logP", response, controls)
                if b is not None:
                    permutation_slopes.append(b["slope"])
            tail = sum(abs(b) >= abs(orig) for b in permutation_slopes)
            out["withinP_shuffle_diagnostic"][response][name] = {
                "permutations_valid": len(permutation_slopes),
                "two_sided_abs_slope_fraction_addone": (1 + tail) / (1 + len(permutation_slopes)),
                "not_confirmatory_p_due_to_unproved_urban_land_source_exchangeability": True,
            }
    # Main within-area result numerics are a diagnostic not a causal threshold.
    assert abs(out["within"]["logP"]["Lycaenidae2022"]["P_plus_land"]["slope"] - (-0.043)) < 0.002
    assert abs(out["within"]["logP"]["LycaenidaeMean"]["P_plus_land"]["slope"] - (-0.107)) < 0.002
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = analysis(root)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "sample": out["sample"],
        "between_P_2022": out["between"]["logP"]["Lycaenidae2022"],
        "within_P_2022": out["within"]["logP"]["Lycaenidae2022"],
        "within_P_3year_mean": out["within"]["logP"]["LycaenidaeMean"],
        "bootstrap_2022": out["sensitivity_bootstrap"]["Lycaenidae2022"],
        "shuffle_2022": out["withinP_shuffle_diagnostic"]["Lycaenidae2022"],
    }, indent=2))


if __name__ == "__main__":
    main()
