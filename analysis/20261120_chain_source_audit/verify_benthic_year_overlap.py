"""Audit fixed EEA station tables; no new ecological model or API requests."""
import collections
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
BLOBS = {2018: "0f30bfbc3dac13f6a182c8a910e9899c04778322",
         2022: "325f199e985d9862ffee5b0edb205c61a7be3833"}
P_FIELDS = {2018: "2018_measured_upstream_P_t_per_year",
            2022: "upstream_2022_EEA_P_measured_t_per_yr"}
N_FIELDS = {2018: "2018_measured_upstream_N_t_per_year",
            2022: "upstream_2022_EEA_N_measured_t_per_yr"}


def main():
    tables, receipts = {}, {}
    for year, blob in BLOBS.items():
        raw = subprocess.check_output(["git", "cat-file", "blob", blob], cwd=ROOT)
        git_hash = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        assert git_hash == blob, "Unexpected Git object"
        rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        tables[year] = {r["site_key"]: r for r in rows}
        assert len(tables[year]) == len(rows), "Duplicate station keys"
        assert all((float(r[P_FIELDS[year]]) > 0) == (float(r[N_FIELDS[year]]) > 0)
                   for r in rows), "N/P presence differs; do not silently pool"
        receipts[str(year)] = {"git_blob": blob, "sha256": hashlib.sha256(raw).hexdigest(),
                               "rows": len(rows), "unique_site_keys": len(tables[year])}
    old, new = tables[2018], tables[2022]
    shared = sorted(old.keys() & new.keys())
    transitions = collections.Counter(
        f"{int(float(old[k][P_FIELDS[2018]]) > 0)}_to_{int(float(new[k][P_FIELDS[2022]]) > 0)}"
        for k in shared)
    result = {
        "scope": "Rejoin fixed station-summary CSVs; not raw EEA download or model replay",
        "inputs": receipts,
        "shared_sites": len(shared), "2018_only": len(old.keys() - new.keys()),
        "2022_only": len(new.keys() - old.keys()),
        "shared_by_country": dict(sorted(collections.Counter(old[k]["countryCode"] for k in shared).items())),
        "shared_same_method": sum(old[k]["procedureClassificationSystem"] == new[k]["procedureClassificationSystem"] for k in shared),
        "shared_same_reach": sum(old[k]["HYRIV_ID"] == new[k]["HYRIV_ID"] for k in shared),
        "reported_source_transitions": dict(sorted(transitions.items())),
        "limitations": ["Reporting changes are not verified WWTP interventions.",
                        "Station-disjoint sets need not have independent upstream sources or catchments.",
                        "Both years use the same EEA Biology release; no independent-study count estimated."]}
    assert [len(old), len(new), len(shared)] == [2634, 2475, 1167]
    assert transitions == {"0_to_0": 1139, "0_to_1": 7, "1_to_0": 4, "1_to_1": 17}
    assert result["shared_same_method"] == 1148
    out = ROOT / "docs/manuscript_riparian_interactions/evidence/BENTHIC_YEAR_OVERLAP_RECHECK_20261120.json"
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
