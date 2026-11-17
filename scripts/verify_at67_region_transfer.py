#!/usr/bin/env python3
"""Check original Austria67 N-v-P and WorldCover ZIP SHA; reweight original errors.

No new observations, model fits, p-values or causal identification.
Original Actions #43890958636 and #43888461430.
"""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

SHA = {
    "np": "f92c7eb204bf86b83a4da9a1ec2cc37dd199a1eb324cdb176f175a5a6d2d732f",
    "land": "09246a4b11cd638fcc5382c00cf6ced9b97668709b9374fd87f827102074837d",
}
def open_source(path, label):
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA[label]:
        raise ValueError(label + " original ZIP SHA256 mismatch")
    with zipfile.ZipFile(path) as z:
        if len(z.namelist()) != 1:
            raise ValueError("Source ZIP must have exactly one JSON")
        return json.loads(z.read(z.namelist()[0]))

def audit(np_source, land_source):
    counts = land_source["source_valid_sites_by_original_3deg_geo_region"]
    regions = sorted(counts)
    if sum(counts.values()) != 67 or len(counts) != 5:
        raise ValueError("Original five region/67 site support drift")
    original = np_source[
        "source_real_2022_benthic_EQR_five_geographic_block_models"
    ]["MSE_by_region"]
    contrasts = np_source[
        "source_real_2022_benthic_EQR_N_v_P_contrasts_adjusted_for_real_TP"
    ]
    tests = {
        "upstream_P_Q": ("B", "P_UP", "up_reported_P_over_Q"),
        "upstream_N_Q": ("B", "N_UP", "up_reported_N_over_Q"),
        "upstream_P_Q_after_presence": (
            "N_PRESENCE", "P_UP_AFTER_N_PRESENCE",
            "up_reported_P_Q_after_source_presence"),
        "source_presence": ("B", "N_PRESENCE", "source_reporting_presence_only"),
        "downstream_P_Q_placebo": ("B", "P_DOWN", "down_P_Q_placebo"),
        "downstream_N_Q_placebo": ("B", "N_DOWN", "down_N_Q_placebo"),
        "N_Q_after_P_Q": ("P_UP", "P_N_UP", "up_N_Q_added_after_P_Q"),
    }
    def delta(b, a, selected, site_weighted=False):
        weight = lambda s: counts[s] if site_weighted else 1
        bn = sum(original[b][s] * weight(s) for s in selected)
        an = sum(original[a][s] * weight(s) for s in selected)
        if bn <= 0: raise ValueError("Invalid heldout baseline MSE")
        return 100 * (1 - an / bn)
    output = {
        "status": "POSTHOC_ORIGINAL_FIVE_REGION_ERRORS_ONLY_NO_NEW_FIT",
        "source_original_signed_actions": [43890958636, 43888461430],
        "source_original_zip_SHA256": SHA,
        "source_original_sites_per_region": counts,
        "same_source_upstream_N_Q_P_Q_Spearman":
            np_source["source_P_N_observational_correlation"][
                "upstream_reported_P_and_N_potential_Q_diluted_signal_spearman"],
        "derived_checks": {},
        "no_2022_riparian_biological_response_or_causal_sewage_effect": True,
    }
    for name, (base, added, source_key) in tests.items():
        original_five_region = contrasts[source_key]["source_5block_equal_MSE_added_gain_pct"]
        computed = delta(base, added, regions)
        if abs(computed - original_five_region) > 1e-9:
            raise ValueError("Cannot reproduce original model score: "+name)
        output["derived_checks"][name] = {
            "original_region_equal_gain_pct": computed,
            "new_67_site_weighted_gain_pct": delta(base, added, regions, True),
            "leave_one_of_five_regions_out_of_saved_errors_siteweighted_pct": {
                region: delta(base, added, [s for s in regions if s != region], True)
                for region in regions
            },
            "source_original_region_specific_gains":
                contrasts[source_key]["source_5block_individually_added_gain_pct"],
            "region_equal_without_16_5_saved_error_gain_pct":
                delta(base, added, [s for s in regions if s != "16_5"]),
        }
    return output

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--np-zip", type=Path, required=True)
    p.add_argument("--land-zip", type=Path, required=True)
    p.add_argument("--out", type=Path)
    a = p.parse_args()
    result = audit(open_source(a.np_zip,"np"), open_source(a.land_zip,"land"))
    payload = json.dumps(result,ensure_ascii=False,indent=2) + "\n"
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(payload,encoding="utf-8")
    print(payload)
