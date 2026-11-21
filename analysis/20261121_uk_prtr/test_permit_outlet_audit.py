"""Real pure audit boundary with synthetic full baseline; no source observations."""
import csv
import io
import json
import tempfile
import unittest
from pathlib import Path

import ea2018_permit_outlet_location_audit as audit

SYNTHETIC_PAIR_SHA = "f4096167ee35e2bd458512605e9f1d70abfce713d2a9a3bd6fe0345f94b8d44a"


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
        result = audit.audit_fixed_candidates(*synthetic_baseline(), expected_pair_sha256=SYNTHETIC_PAIR_SHA)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["fixed_295_native_permit_pair_sha256"], SYNTHETIC_PAIR_SHA)
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
            audit.audit_fixed_candidates(*args, expected_pair_sha256=SYNTHETIC_PAIR_SHA)

    def test_count_preserving_permit_pair_swap_stops_before_location_read(self):
        args = list(synthetic_baseline())
        old = next(pair for pair in args[2] if pair[0] == "AN-SYNTHETIC100")
        new = (old[0], ("AN", "FAKE-REPLACEMENT", "1", "A", "1"))
        for index in (1, 2):
            args[index] = args[index] - {old} | {new}
        args[4] = args[4] - {old[1]} | {new[1]}
        def unread_locations():
            raise AssertionError("Changed identity must stop before location consumption")
            yield
        args[5] = unread_locations()
        with self.assertRaises(ValueError):
            audit.audit_fixed_candidates(*args, expected_pair_sha256=SYNTHETIC_PAIR_SHA)

    def test_same_size_duplicate_replacing_required_key_stops(self):
        args = list(synthetic_baseline())
        rows = args[5]
        self.assertEqual(len(rows), 295)
        args[5] = rows[1:] + [rows[1]]
        self.assertEqual(len(args[5]), 295)
        with self.assertRaises(ValueError):
            audit.audit_fixed_candidates(*args, expected_pair_sha256=SYNTHETIC_PAIR_SHA)

    def test_location_csv_rejects_truncated_and_extra_cells(self):
        header = ",".join(audit.KEY_FIELDS["consents"] + audit.FIELDS) + "\n"
        for row in ("AN,FAKE,1,A,1,SU0000000000,SU1000010000\n",
                    "AN,FAKE,1,A,1,SU0000000000,SU1000010000,SU2000020000,EXTRA\n"):
            with self.subTest(width=len(row.split(","))):
                with self.assertRaises(ValueError):
                    list(audit.location_csv_rows(io.StringIO(header + row)))

    def test_location_csv_preserves_explicit_empty_coordinate(self):
        header = ",".join(audit.KEY_FIELDS["consents"] + audit.FIELDS) + "\n"
        rows = list(audit.location_csv_rows(io.StringIO(header + "AN,FAKE,1,A,1,,SU1000010000,SU2000020000\n")))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["OUTLET_GRID_REF"], "")

    def test_location_csv_rejects_unclosed_quoted_record(self):
        header = ",".join(audit.KEY_FIELDS["consents"] + audit.FIELDS) + "\n"
        with self.assertRaises((ValueError, csv.Error)):
            list(audit.location_csv_rows(io.StringIO(header + 'AN,FAKE,1,A,1,"SU0000000000')))

    def test_registry_pages_pin_coordinate_bearing_original_bytes(self):
        baseline = audit.identity.ROOT / "results/uk2018_permit_identity_audit_20261121/run_44907990658.json"
        registry = json.loads(baseline.read_text())["current_WQE_registry"]
        result = audit.verify_registry_pages(registry)
        self.assertEqual(result["verified_pages"], 23)
        registry["page_receipts"][0]["source_sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            audit.verify_registry_pages(registry)

    def test_registry_pages_reject_incomplete_or_duplicate_receipts(self):
        baseline = audit.identity.ROOT / "results/uk2018_permit_identity_audit_20261121/run_44907990658.json"
        registry = json.loads(baseline.read_text())["current_WQE_registry"]
        for receipts in (registry["page_receipts"][:-1],
                         registry["page_receipts"][1:] + [registry["page_receipts"][1]]):
            with self.subTest(pages=len(receipts)):
                with self.assertRaises(ValueError):
                    audit.verify_registry_pages({**registry, "page_receipts": receipts})

    def test_registry_reference_file_is_itself_pinned(self):
        baseline = audit.identity.ROOT / "results/uk2018_permit_identity_audit_20261121/run_44907990658.json"
        registry = json.loads(baseline.read_text())["current_WQE_registry"]
        with tempfile.TemporaryDirectory() as folder:
            altered = Path(folder) / "changed.json"
            altered.write_bytes(baseline.read_bytes() + b" ")
            with self.assertRaises(ValueError):
                audit.verify_registry_pages(registry, baseline_path=altered)

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
