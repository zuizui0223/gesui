"""Check whether stricter published DMR support changes existing model cohorts."""
import argparse
import hashlib
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--cache", type=Path, required=True)
args = p.parse_args()
files = {
    "crosswalk": "NWCA2021_568WETLAND_160POTW_EXACT_PERMIT_FEATURE_SOURCE_BEFORE_2021DMR.json",
    "cohorts": "NWCA2021_CAL2021_ACTUAL_DMR_TNTP_NEARFAR_STATEECO_SOURCE_ONLY_GATE.json",
    "strict": "EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"}
data = {k: json.loads((args.cache / v).read_text()) for k, v in files.items()}
strict = data["strict"]["source_original_2021_exact_permit_feature_DMR_reporting_summary"]
n = {r["NPDES_exact_permit_feature_SHA256"] for r in strict if r["reported_original_numeric_2021_TN00600_period_end_dates"]}
ph = {r["NPDES_exact_permit_feature_SHA256"] for r in strict if r["reported_original_numeric_2021_TP00665_period_end_dates"]}
pairs = data["crosswalk"]["registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology"]
sites = data["crosswalk"]["original_568_NWCA_source_site_municipal_permit_geometries"]
regions = {r["source_Wetland_SITE_ID_VISIT_SHA256"]: (r["original_state"], r["original_ecoregion"]) for r in sites}
out = {}
for name, features in [("DMR_actual2021_BOTH_TN_TP_same_EXO", n & ph), ("DMR_actual2021_TN_or_TP", n | ph)]:
    near = set()
    for r in pairs:
        raw = r["original_public_NPDES_permit_id"].strip().upper() + "|" + r["original_public_NPDES_feature_number"].strip().upper()
        if hashlib.sha256(raw.encode()).hexdigest() in features:
            near.add(r["source_Wetland_SITE_ID_VISIT_SHA256"])
    original = data["cohorts"]["original_source_DMR_2021_ANY_NP_and_BOTH_NP_strict_same_state_ecoregion_controls"][name]
    old_n = set(original["original_NWCA2021_NEAR_source_site_SHA256s"])
    old_f = set(original["original_NWCA2021_FAR_source_site_SHA256s"])
    assert (old_n | old_f | near) <= regions.keys()
    far_regions = {regions[s] for s in old_f}
    matched_n = {s for s in near if regions[s] in far_regions}
    near_regions = {regions[s] for s in matched_n}
    matched_f = {s for s in old_f if regions[s] in near_regions}
    out[name] = {"strict_near_before_region_matching": len(near),
                 "old_near": len(old_n), "strict_matched_near": len(matched_n),
                 "old_far": len(old_f), "strict_matched_far": len(matched_f),
                 "near_removed": len(old_n - matched_n), "near_added": len(matched_n - old_n),
                 "far_removed": len(old_f - matched_f),
                 "same_exact_near_and_far_membership": old_n == matched_n and old_f == matched_f,
                 "matched_states": len({regions[s][0] for s in matched_n}),
                 "matched_state_ecoregions": len(near_regions),
                 "both_nutrients_30_site_original_gate_pass": len(near) >= 30 if "BOTH" in name else None}
receipt = {"scope": "Same registered geometry, stricter published monitoring support; no plant outcomes read",
           "input_sha256": {v: hashlib.sha256((args.cache / v).read_bytes()).hexdigest() for v in files.values()},
           "cohorts": out,
           "limitations": ["Source receipts not original national DMR replay", "No hydrologic delivery established",
                           "No raw plant model replay", "Old source thresholds retained"]}
target = Path(__file__).resolve().parents[2] / "docs/manuscript_riparian_interactions/evidence/NWCA_STRICT_COHORT_CHECK_20261120.json"
target.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
print(json.dumps(out))
