#!/usr/bin/env python3
"""POSTHOC: 2021-2023 Lycaenidae photographic composition vs fixed 2022 WWTP P geography.

Read only frozen public aggregate JSON already committed in this draft branch.
No butterfly census, wastewater exposure time series, or causal effect inference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

FAMILIES = ("Lycaenidae", "Pieridae", "Papilionidae")
YEARS = (2021, 2022, 2023)
SEED = 20261109
BOOTSTRAPS = 3000


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(row):
    return row["countryCode"], int(float(row["lat_floor"])), int(float(row["lon_floor"]))


def source_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def standardize(values):
    mean = sum(values) / len(values)
    sd = math.sqrt(sum((x - mean) ** 2 for x in values) / len(values))
    if sd < 1e-12:
        raise ValueError("Zero-variance source feature")
    return [(x - mean) / sd for x in values]


def solve(x, y):
    """Full-rank unweighted OLS normal equation, Gauss-Jordan with pivoting."""
    n, p = len(x), len(x[0])
    a = [[0.0] * (p + 1) for _ in range(p)]
    for i in range(n):
        for j in range(p):
            a[j][p] += x[i][j] * y[i]
            for k in range(p):
                a[j][k] += x[i][j] * x[i][k]
    for k in range(p):
        pivot = max(range(k, p), key=lambda j: abs(a[j][k]))
        if abs(a[pivot][k]) < 1e-8:
            return None
        a[k], a[pivot] = a[pivot], a[k]
        div = a[k][k]
        for j in range(k, p + 1):
            a[k][j] /= div
        for i in range(p):
            if i == k:
                continue
            factor = a[i][k]
            for j in range(k, p + 1):
                a[i][j] -= factor * a[k][j]
    return [r[p] for r in a]


def slope(rows, response):
    b = solve([r["_design"] for r in rows], [r["_y"][response] for r in rows])
    return None if b is None else b[-1]


def build_rows(photo21, photo22, photo23, land, land_hash):
    if photo21["original_2021_WorldCover2022_photo_source_full_SHA256"] != land_hash:
        raise ValueError("2021 family sources do not match original land/2022 source")
    if photo23["original_2022_WorldCover_land_photo_source_sha256"] != land_hash:
        raise ValueError("2023 family sources do not match original land/2022 source")
    if not photo23["2023_outcome_year_source_gate_pass"]:
        raise ValueError("2023 source coverage not adequate")
    l = land["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    p21 = {key(r): r for r in photo21["source_2021_same43_cell_historical_butterfly_photo_records"]}
    p22 = {key(r): r for r in photo22["all_fixed_cells_2022_family_photo_source_receipts"]}
    p23 = {key(r): r for r in photo23["source_2023_all_original43_butterfly_photo_count_and_receipts"]}
    if any(len(p) != 43 for p in (p21, p23)) or len(l) != 43:
        raise ValueError("Unexpected source-area cohort size")
    rows = []
    for a in l:
        ky = key(a)
        if ky not in p21 or ky not in p22 or ky not in p23:
            raise ValueError("Missing original source-area cell")
        h, j, k = p21[ky], p22[ky], p23[ky]
        n = {
            2021: int(h["2021_three_family_photo_total"]),
            2022: int(j["three_family_total_2022"]),
            2023: int(k["2023_three_family_photo_total"]),
        }
        if n[2022] != int(a["2022_total_3family_photos"]):
            raise ValueError("2022 photo source disagrees with land original")
        if n[2021] < 15 or n[2022] < 20 or n[2023] < 20:
            continue
        r = {
            "country": a["countryCode"],
            "lat": float(a["lat_floor"]),
            "lon": float(a["lon_floor"]),
            "block": a["countryCode"] + "|" + a["country_2deg_block"],
            "logP": math.log1p(float(a["registered_P_t_year"])),
            "logPlant": math.log1p(float(a["2021_Plantae_photo_n"])),
            "built": float(a["WC2021_512_built_land_fraction"]),
            "crop": float(a["WC2021_512_crop_land_fraction"]),
            "logArea": math.log(float(a["WC2021_512_country_terrestrial_area_km2"])),
            "_y": {},
        }
        family_yr = {}
        for f in FAMILIES:
            family_yr[(f, 2021)] = int(h["2021_" + f + "_photo_n"])
            family_yr[(f, 2022)] = int(j[f + "_2022_photo_n"])
            family_yr[(f, 2023)] = int(k["2023_" + f + "_photo_n"])
        for yr in YEARS:
            if sum(family_yr[(f, yr)] for f in FAMILIES) != n[yr]:
                raise ValueError("Photo family counts do not sum to same-year total")
        for f in FAMILIES:
            for yr in YEARS:
                v = family_yr[(f, yr)]
                r["_y"][f + str(yr)] = math.log((v + 0.5) / (n[yr] - v + 0.5))
        rows.append(r)
    if len(rows) != 42 or Counter(r["country"] for r in rows) != {"DE": 27, "AT": 15}:
        raise ValueError("Unexpected balanced original geographic family-photo sample")
    if len({r["block"] for r in rows}) != 22:
        raise ValueError("Unexpected 2-degree source-area block count")
    feature_names = ("lat", "lon", "logPlant", "built", "crop", "logArea")
    zfeatures = {f: standardize([r[f] for r in rows]) for f in feature_names}
    pz = [0.0] * len(rows)
    for c in ("DE", "AT"):
        ids = [i for i, r in enumerate(rows) if r["country"] == c]
        p = standardize([rows[i]["logP"] for i in ids])
        for i, z in zip(ids, p):
            pz[i] = z
    for i, r in enumerate(rows):
        r["_design"] = [1.0, float(r["country"] == "DE")] + [zfeatures[f][i] for f in feature_names] + [pz[i]]
    return rows


def contrasts(b):
    return {
        "2021": b["2021"],
        "2022": b["2022"],
        "2023": b["2023"],
        "delta_2022_minus_2021": b["2022"] - b["2021"],
        "delta_2023_minus_2022": b["2023"] - b["2022"],
        "year2022_minus_mean_2021_2023": b["2022"] - (b["2021"] + b["2023"]) / 2,
    }


def slopes_all(rows):
    result = {}
    for f in FAMILIES:
        b = {str(y): slope(rows, f + str(y)) for y in YEARS}
        if any(v is None for v in b.values()):
            return None
        result[f] = contrasts(b)
    return result


def percentile(arr, frac):
    arr = sorted(arr)
    return arr[int(math.floor((len(arr) - 1) * frac))]


def run(rows):
    original = slopes_all(rows)
    groups = list(dict.fromkeys(r["block"] for r in rows))
    bygroup = {g: [r for r in rows if r["block"] == g] for g in groups}
    draws = {f: {k: [] for k in original[f]} for f in FAMILIES}
    state, valid = SEED, 0
    for _ in range(BOOTSTRAPS):
        sample = []
        for __ in groups:
            state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
            group = groups[int((state / (2 ** 32)) * len(groups))]
            sample.extend(bygroup[group])
        fit = slopes_all(sample)
        if fit is None:
            continue
        valid += 1
        for f in FAMILIES:
            for name in fit[f]:
                draws[f][name].append(fit[f][name])
    out = {
        "schema": "gesui_2021_2023_lycaenid_sourceP_photo_composition_fixed_temporal_posthoc_v01",
        "analysis_is_posthoc_after_positive_2022_source_photo": True,
        "N_country_1degree_cells": len(rows),
        "N_country_2degree_geographic_blocks": len(groups),
        "model": "OLS cell photo-family logit; source is 2022 original reported WWTP P mass log1p scaled within DE/AT; intercept, country, latitude, longitude, 2021 Plantae photo count, WorldCover built, crop, terrestrial land area (all standardized). Two-degree cluster resampling. Outcomes are photographic reporting shares, not population counts.",
        "source_covariate_is_fixed_2022_reported_mass_in_every_photo_year": True,
        "not_a_difference_in_treatment_exposure": True,
        "seed": SEED, "bootstrap_target_draws": BOOTSTRAPS, "bootstrap_valid_draws": valid,
        "family_results": {},
    }
    for f in FAMILIES:
        out["family_results"][f] = {}
        for name in original[f]:
            d = draws[f][name]
            los = []
            for group in groups:
                fitted = slopes_all([r for r in rows if r["block"] != group])
                if fitted:
                    los.append(fitted[f][name])
            out["family_results"][f][name] = {
                "coefficient_per_within_country_1SD_logP": original[f][name],
                "descriptive_cluster_bootstrap_95pct": [percentile(d, .025), percentile(d, .975)],
                "cluster_bootstrap_fraction_positive": sum(v > 0 for v in d) / len(d),
                "leave_one_2degree_block_out_positive_count": sum(v > 0 for v in los),
                "leave_one_2degree_block_out_total": len(los),
                "leave_one_2degree_block_out_range": [min(los), max(los)],
            }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    paths = {
        "2021_family": root / "results/lycaenidae_2021_prior_butterfly_photo_source/2021_HISTORICAL_3FAMILY_PHOTO_COUNTS_ORIGINAL43_SOURCE_ONLY.json",
        "2022_family": root / "results/lycaenidae_de_at_2022_macro_source_photo/2022_DE_AT_LYCAENIDAE_PHOTO_SHARE_VS_WWTP_SOURCE_EXPLORATORY.json",
        "2023_family": root / "results/lycaenidae_2023_photo_temporal_source/2023_INDEPENDENT_YEAR_FAMILY_PHOTO_TOTALS_ORIGINAL43_SOURCE_QC.json",
        "2022_land": root / "results/lycaenidae_worldcover2021_source_only/2022_LYCAENIDAE_ESAWC2021_DE_AT_URBAN_CROP_LANDAREA_POSTHOC_BLOCKCV.json",
    }
    raw = {k: source_json(v) for k, v in paths.items()}
    rows = build_rows(raw["2021_family"], raw["2022_family"], raw["2023_family"], raw["2022_land"], sha256(paths["2022_land"]))
    out = run(rows)
    out["source_sha256"] = {k: sha256(v) for k, v in paths.items()}
    out["photo_observation_data_not_independent_random_census"] = True
    out["no_direct_WWTP_to_butterfly_substrate_exposure_confirmed"] = True
    test = out["family_results"]["Lycaenidae"]
    # Numeric replay guardrails, not a confirmatory significance threshold.
    assert abs(test["2022"]["coefficient_per_within_country_1SD_logP"] - 0.15407574448824474) < 1e-7
    assert abs(test["delta_2022_minus_2021"]["coefficient_per_within_country_1SD_logP"] - 0.09308411494709429) < 1e-7
    assert abs(test["delta_2023_minus_2022"]["coefficient_per_within_country_1SD_logP"] - (-0.03853175738006566)) < 1e-7
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"n": out["N_country_1degree_cells"], "blocks": out["N_country_2degree_geographic_blocks"], "valid": out["bootstrap_valid_draws"], "Lycaenidae": test}, indent=2))


if __name__ == "__main__":
    main()
