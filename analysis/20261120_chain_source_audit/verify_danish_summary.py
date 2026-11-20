"""Check published aggregate arithmetic; does not replay licensed eBMS records."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "docs/manuscript_riparian_interactions"
snapshot = DOC / "PR_SOURCE_SNAPSHOT_20261120.json"
prs = json.loads(snapshot.read_text(encoding="utf-8"))["prs"]
selected = [p for p in prs if 37 <= p["number"] <= 43]
assert sorted(p["number"] for p in selected) == list(range(37, 44))
# PR43 fixed report, c960df1064a43671e74dda7eaae1445a4316a644.
# Order: Lycaenidae, Pieridae, Papilionidae, Nymphalidae, Hesperiidae.
near = [382, 680, 0, 3822, 372]
far = [644, 2602, 0, 3544, 192]


def odds_ratio(a, b, families):
    return (a[0] / sum(a[1:families])) / (b[0] / sum(b[1:families]))


results = {"rate_ratio": (382 / 539) / (644 / 757),
           "three_family_OR": odds_ratio(near, far, 3),
           "five_family_OR": odds_ratio(near, far, 5)}
assert results["three_family_OR"] > 1
assert results["five_family_OR"] < 1
assert results["rate_ratio"] < 1
deletions = {}
for name, family, n, f in [("Pieris brassicae", 1, 145, 1439),
                            ("Maniola jurtina", 3, 1828, 870),
                            ("Polyommatus icarus", 0, 142, 363),
                            ("Plebejus argus", 0, 63, 0)]:
    a, b = near.copy(), far.copy()
    a[family] -= n
    b[family] -= f
    deletions[name] = {"three_family_OR": odds_ratio(a, b, 3),
                       "five_family_OR": odds_ratio(a, b, 5)}
assert deletions["Maniola jurtina"]["five_family_OR"] > 1
assert deletions["Polyommatus icarus"]["five_family_OR"] > 1
receipt = {
    "scope": "Published aggregate arithmetic only; no new ecological model or uncertainty estimate",
    "snapshot_sha256": hashlib.sha256(snapshot.read_bytes()).hexdigest(),
    "series": "Denmark eBMS v1.4; PRs are dependent analyses, not seven studies",
    "prs": [{k: p[k] for k in ("number", "headRefOid", "url")} for p in selected],
    "source_report": "https://github.com/zuizui0223/gesui/blob/c960df1064a43671e74dda7eaae1445a4316a644/docs/RESULTS_DK2025_NEW399_SPECIES_3VS5FAMILY_COMPOSITION_2021_HABITAT_20261110.md",
    "published_aggregate_counts": {"near": near, "remote": far,
                                   "near_linked_days": 539, "remote_linked_days": 757},
    "recomputed": results,
    "algebraic_species_exclusions_not_interventions": deletions,
    "limitations": ["No raw event/occurrence rejoin in this check",
                    "No independent variance derivable from pooled individual totals",
                    "Observation dates without occurrence rows are not verified zero counts",
                    "Recorded adults and larval-host HNS are different outcomes"]}
output = DOC / "evidence/DANISH_AGGREGATE_RECHECK_20261120.json"
output.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps(results))
