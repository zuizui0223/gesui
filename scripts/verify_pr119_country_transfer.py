#!/usr/bin/env python3
"""Verify original PR119 signed ZIP and recompute country-transfer error aggregates.

Scientific scope: posthoc reaggregation of SAVED country-heldout prediction
errors only, not a new biological model, causal test, confidence interval,
independent replication, or any source relaxation. Stdlib-only.
Source: gesui Actions #43901046421 / artifact #14240413667.
"""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ZIP_SHA256 = "a22ea230c8e6beb4533bc696504194879927e9fd6a553377fba5255ce1464bc2"
MEMBER = "GESUI_STRICT_SAME_HOST_NF_H7_PLANTF_EQR_CROSSDOMAIN_RESULTS.json"
SCHEMA = "GESUI_20261117_H7_BUTTERFLY_SAME_ACTUAL_FERRERPARIS_LARVALHOST_EIV_NF_MATCHED_PLANT_AND_WATER_EQR_V1"
COUNTRIES = ["FR", "HR", "IT", "LT", "NL", "NO"]


def gain(base, added):
    if base <= 0:
        raise ValueError("MSE baseline <=0")
    return 100 * (1 - added / base)


def aggregate(rows):
    bn = sum(x["baseline_MSE"] * x["n_H7period"] for x in rows)
    an = sum(x["added_MSE"] * x["n_H7period"] for x in rows)
    return {
        "n_H7period": sum(x["n_H7period"] for x in rows),
        "n_H6_from_source_country_rows": sum(x["n_true_H6"] for x in rows),
        "H7weighted_MSE_gain_pct": gain(bn, an),
        "country_equal_MSE_gain_pct": gain(
            sum(x["baseline_MSE"] for x in rows),
            sum(x["added_MSE"] for x in rows)
        ),
        "arithmetic_mean_six_country_percent_gains":
            sum(x["country_MSE_improvement_pct"] for x in rows) / len(rows),
    }


def audit(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != ZIP_SHA256:
        raise ValueError("Original signed ZIP digest mismatch")
    with zipfile.ZipFile(path) as z:
        if z.namelist() != [MEMBER]:
            raise ValueError("Original member list drift")
        doc = json.loads(z.read(MEMBER))
    if doc["schema"] != SCHEMA or not doc["status"].startswith("PASS_SAME_PUBLISHED_HOST"):
        raise ValueError("Original source schema / status mismatch")
    if (doc["source_siteperiods_with_geographic_species_support"],
        doc["source_true_H6_of_supported_sites"]) != (447, 89):
        raise ValueError("Original 447 source H7period / 89 H6 drift")
    results = doc["source_same_host_family_score_water_land_2022_observational_results"]
    pairs = {
        "plant_moisture_after_plantN_and_butterfly_hostN":
            "original_paired_plant_MOIST_added_after_N_and_bfly_samehostN",
        "benthic_EQR_after_plantN_F_and_butterfly_hostN":
            "actual_aquatic_EQR_after_BOTH_plant_traits_and_bfly_samehostN",
    }
    output = {
        "original_actions_run": 43901046421,
        "original_artifact_id": 14240413667,
        "original_zip_sha256": ZIP_SHA256,
        "source_H7period": 447,
        "source_H6": 89,
        "posthoc_aggregate_only_no_new_ecological_fit": True,
        "six_country_eligible_subset": {},
    }
    for label, key in pairs.items():
        original = results[key]["source_country_H6_purged_prediction"]
        rows = original["source_country_results"]
        if len(rows) != 6 or [x["country_holdout"] for x in rows] != COUNTRIES:
            raise ValueError("Six-country source/order mismatch")
        for x in rows:
            if abs(gain(x["baseline_MSE"], x["added_MSE"]) -
                   x["country_MSE_improvement_pct"]) > 1e-9:
                raise ValueError("Country heldout MSE mismatch: " + x["country_holdout"])
        full = aggregate(rows)
        if full["n_H7period"] != 424 or full["n_H6_from_source_country_rows"] != 81:
            raise ValueError("Unexpected six-country 424 / 81 supported source")
        if abs(full["country_equal_MSE_gain_pct"] -
               original["country_equal_added_MSE_gain_pct"]) > 1e-9:
            raise ValueError("Original country-equal result mismatch")
        output["six_country_eligible_subset"][label] = {
            "country_summary": full,
            "country_holdout_MSE_gain_pct": {
                x["country_holdout"]: gain(x["baseline_MSE"], x["added_MSE"])
                for x in rows
            },
            "leave_one_country_OUT_of_original_errors_H7weighted_gain_pct": {
                x["country_holdout"]: aggregate(
                    [y for y in rows if y["country_holdout"] != x["country_holdout"]]
                )["H7weighted_MSE_gain_pct"]
                for x in rows
            },
            "number_of_positive_country_holdouts":
                sum(x["country_MSE_improvement_pct"] > 0 for x in rows),
            "original_fixed_errors_country_CI":
                original["country_equal_fixed_errors_bootstrap95"],
        }
    original_time = results["two_period_source_correlations"]
    output["original_two_period_H7_temporal_falsifier"] = {
        "n_H7_repeat": original_time["source_two_periods_H7_same_match_n"],
        "delta_hostF_vs_plantF_Spearman":
            original_time["delta_fly_hostF_vs_plantF_Spearman"],
        "delta_hostF_vs_EQR_Spearman":
            original_time["delta_fly_hostF_vs_EQR_Spearman"],
    }
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-zip", type=Path, required=True)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    result = audit(args.artifact_zip)
    payload = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(payload, encoding="utf-8")
    print(payload)
