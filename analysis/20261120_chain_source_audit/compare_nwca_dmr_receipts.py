"""Compare existing public source receipts without publishing native permit IDs."""
import argparse
import hashlib
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--cache", type=Path, required=True)
args = parser.parse_args()
names = ["NWCA2021_160POTW_ORIGINAL_CAL2021_DMR_TN_TP_EXO_FEATURE_SOURCE_GATE.json",
         "EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"]
a, b = [json.loads((args.cache / name).read_text()) for name in names]
rows = a["full_source_permit_feature_parameter_monitoring_receipts"]
def key(r):
    raw = r["original_public_NPDES_permit_id"].strip().upper() + "|" + r["original_public_NPDES_feature_number"].strip().upper()
    return hashlib.sha256(raw.encode()).hexdigest()

target = {key(r) for r in rows}
sets_a, sets_b = {}, {}
summary = {}
for code, field in [("00600", "TN00600"), ("00665", "TP00665")]:
    sets_a[code] = {key(r) for r in rows if r["original_pollutant_code"] == code
                    and r["original_2021_valid_numeric_monitor_period_end_dates_unique"]}
    sets_b[code] = {r["NPDES_exact_permit_feature_SHA256"] for r in b["source_original_2021_exact_permit_feature_DMR_reporting_summary"]
                    if r["NPDES_exact_permit_feature_SHA256"] in target
                    and r[f"reported_original_numeric_2021_{field}_period_end_dates"]}
    summary[code] = {"PR53": len(sets_a[code]), "PR54_on_PR53_reported_feature_keys": len(sets_b[code]),
                     "only_PR53": len(sets_a[code] - sets_b[code]), "only_PR54": len(sets_b[code] - sets_a[code])}
summary["both"] = {"PR53": len(sets_a["00600"] & sets_a["00665"]),
                   "PR54_on_PR53_reported_feature_keys": len(sets_b["00600"] & sets_b["00665"])}
codes = {}
for r in rows:
    for code, count in r["original_DMR_qualifier_codes"].items():
        codes[code] = codes.get(code, 0) + count
receipt = {"scope": "Existing aggregate receipts; no nationwide DMR or biological model replay",
           "input_sha256": {name: hashlib.sha256((args.cache / name).read_bytes()).hexdigest() for name in names},
           "feature_counts": summary,
           "PR53_all_parameter_qualifier_counts_not_deduplicated": codes,
           "interpretation": "Code predicates differ in monitoring location, NODI and estimated-value acceptance. Aggregate receipts do not assign each discrepant feature to a unique cause.",
           "native_source_ids_omitted": True}
out = Path(__file__).resolve().parents[2] / "docs/manuscript_riparian_interactions/evidence/NWCA_DMR_RECEIPT_COMPARISON_20261120.json"
out.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
print(json.dumps(summary))
