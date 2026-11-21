"""Real pure audit boundary with synthetic full baseline; no source observations."""
import json
import tempfile
import unittest
from pathlib import Path

import ea2018_permit_outlet_location_audit as audit


def synthetic_baseline():
    candidates = {"AN-SYNTHETIC" + str(i) for i in range(930)}
    pairs = {("AN-SYNTHETIC" + str(i), ("AN", "FAKE" + str(i), "1", "A", "1")) for i in range(346)}
    confirmed = {p for p in pairs if int(p[0].replace("AN-SYNTHETIC", "")) < 295}
    active = {p for p in confirmed if int(p[0].replace("AN-SYNTHETIC", "")) < 16}
    rows = [dict(EA_REGION=k[0], PERMIT_NUMBER=k[1], PERMIT_VERSION=k[2],
                 OUTLET_NUMBER=k[3], EFFLUENT_NUMBER=k[4],
                 OUTLET_GRID_REF="" if k[1] == "FAKE0" else "SU0000000000",
                 DISCHARGE_NGR="SU1000010000", EFFLUENT_GRID_REF="SU2000020000") for _, k in confirmed]
    matching = {p: {"outlet_coord": (400000, 100000)} for p in candidates}
    registry = {p: (400600, 100000) for p in candidates}
    return candidates, pairs, confirmed, active, {k for _, k in confirmed}, rows, matching, registry


class AuditTests(unittest.TestCase):
    def test_real_baseline_to_location_summary_preserves_missing(self):
        result = audit.audit_fixed_candidates(*synthetic_baseline())
        self.assertIsInstance(result, dict)
        self.assertEqual(result["prior_identity_partition"]["with_ea_all_corroborated_pair"], 295)
        result = result["fixed_295_location_support"]
        self.assertEqual(result["candidate_points"], 295)
        self.assertEqual(result["outlet_status_counts"], {"missing": 1, "valid": 294})
        self.assertEqual(result["distances"]["original_EEA_outlet"]["250m"]["within"], 294)
        self.assertNotIn("FAKE", json.dumps(result))
        self.assertNotIn("SYNTHETIC", json.dumps(result))

    def test_baseline_mismatch_stops_before_new_location_rows(self):
        args = list(synthetic_baseline())
        args[3] = set()
        def private_rows():
            raise AssertionError("PRIVATE_SOURCE_NOT_TO_BE_READ")
            yield
        args[5] = private_rows()
        with self.assertRaises(ValueError):
            audit.audit_fixed_candidates(*args)

    def test_controlled_stop_discards_partial_counts_and_exception_text(self):
        def failure(out):
            out["stage"] = "outlet_locations"
            out["partial"] = "PRIVATE_SOURCE_VALUE"
            raise ValueError("PRIVATE_SOURCE_VALUE")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "out.json"
            self.assertEqual(audit.guarded_run(path, failure), 1)
            result = json.loads(path.read_text())
            self.assertEqual(result["status"], "PERMIT_OUTLET_LOCATION_AUDIT_STOP")
            self.assertEqual(result["error_stage"], "outlet_locations")
            self.assertNotIn("PRIVATE_SOURCE_VALUE", path.read_text())
            self.assertNotIn("partial", result)

    def test_success_serializes_only_completed_producer(self):
        def success(out):
            out["status"] = audit.SUCCESS
            out["aggregate"] = {"n": 1}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "out.json"
            self.assertEqual(audit.guarded_run(path, success), 0)
            self.assertEqual(json.loads(path.read_text())["aggregate"], {"n": 1})


if __name__ == "__main__":
    unittest.main()
