"""Synthetic identity/date tests; no original permit or sampling-point rows."""
import datetime as dt
import json
import unittest
import tempfile
from pathlib import Path

import ofwat2018_permit_sample_identity_audit as runner

import permit_sample_identity_support as audit


def row(**changes):
    value = dict(year=2018, sample_samplingpoint_notation="AN-SYNTH_A",
                 EA_REGION="AN", permit_number="00SYNTH/P", version="01",
                 outlet_number="001", effluent_number="1", _sheet="first")
    value.update(changes)
    return value


class IdentityTests(unittest.TestCase):
    def test_key_normalizes_numeric_components_but_preserves_permit(self):
        self.assertEqual(audit.permit_key(row(), "ofwat"),
                         ("AN", "00SYNTH/P", "1", "1", "1"))

    def test_missing_version_is_not_filled_from_another_record(self):
        self.assertIsNone(audit.permit_key(row(version=None), "ofwat"))

    def test_alphanumeric_outlet_identifiers_are_not_discarded_as_non_numeric(self):
        self.assertEqual(audit.permit_key(row(outlet_number="A"), "ofwat"),
                         ("AN", "00SYNTH/P", "1", "A", "1"))
        pairs, qa = audit.ofwat_pairs([row(outlet_number="A")])
        self.assertEqual(len(pairs), 1)
        self.assertEqual(qa["complete_identity_rows"], 1)

    def test_noninteger_identifier_text_is_preserved_without_rounding(self):
        self.assertEqual(audit.permit_key(row(version="1.5"), "ofwat"),
                         ("AN", "00SYNTH/P", "1.5", "1", "1"))

    def test_pair_losses_remain_distinct_from_point_losses(self):
        pairs, _ = audit.ofwat_pairs([row(), row(version=2)])
        out = audit.summarize_candidates({"AN-SYNTH_A"}, pairs, pairs, set())
        self.assertEqual(out["with_corroboration_only_in_all_not_active"], 1)
        self.assertEqual(out["ea_active_corroborated_distinct_pairs"], 0)
        self.assertEqual(out["ea_all_corroborated_distinct_pairs_absent_from_active"], 2)

    def test_complete_five_key_aliases_join_without_losing_version(self):
        source = dict(EA_REGION="AN", PERMIT_NUMBER="00SYNTH/P",
                      PERMIT_VERSION="02", OUTLET_NUMBER="1", EFFLUENT_NUMBER="1")
        self.assertEqual(audit.permit_key(source, "consents"),
                         ("AN", "00SYNTH/P", "2", "1", "1"))
        self.assertNotEqual(audit.permit_key(source, "consents"),
                            audit.permit_key(row(), "ofwat"))

    def test_duplicate_rows_and_sheets_do_not_create_extra_sites(self):
        pairs, qa = audit.ofwat_pairs([row(), row(_sheet="second"), row(year=2017)])
        self.assertEqual(pairs, {("AN-SYNTH_A", ("AN", "00SYNTH/P", "1", "1", "1"))})
        self.assertEqual(qa["year_rows"], 2)
        self.assertEqual(qa["duplicate_identity_rows"], 1)
        self.assertEqual(qa["distinct_points"], 1)

    def test_region_local_collisions_do_not_match(self):
        pairs, qa = audit.ofwat_pairs([
            row(), row(EA_REGION="NW"),
            row(EA_REGION="NW", sample_samplingpoint_notation="NW-SYNTH_A")])
        self.assertEqual(len(pairs), 2)
        self.assertEqual(qa["excluded"]["point_region_mismatch"], 1)

    def test_midlands_spelling_is_explicit_and_not_guessed(self):
        pairs, qa = audit.ofwat_pairs([
            row(EA_REGION="MI", sample_samplingpoint_notation="MD-SYNTH_A"),
            row(EA_REGION="MI", sample_samplingpoint_notation="MI-SYNTH_A"),
            row(EA_REGION="AN", sample_samplingpoint_notation="SYNTH_A")])
        self.assertEqual(pairs, {("MD-SYNTH_A", ("MI", "00SYNTH/P", "1", "1", "1"))})
        self.assertEqual(qa["complete_identity_rows"], 1)

    def test_two_versions_or_outlets_remain_distinct(self):
        pairs, qa = audit.ofwat_pairs([row(), row(version=2), row(outlet_number=2)])
        self.assertEqual(len(pairs), 3)
        self.assertEqual(qa["points_with_multiple_keys"], 1)
        self.assertEqual(qa["points_with_multiple_keys_ignoring_version"], 1)

    def test_candidate_counts_and_public_output_preserve_ambiguity(self):
        pairs, _ = audit.ofwat_pairs([row(), row(version=2),
            row(sample_samplingpoint_notation="AN-SYNTH_B", permit_number="SYNTH/B")])
        corroborated = {p for p in pairs if p[0] == "AN-SYNTH_B"}
        out = audit.summarize_candidates({"AN-SYNTH_A", "AN-SYNTH_B", "AN-SYNTH_C"},
                                         pairs, corroborated, set())
        self.assertEqual(out["candidate_points"], 3)
        self.assertEqual(out["with_multiple_ofwat_keys"], 1)
        self.assertEqual(out["without_ofwat_pair"], 1)
        self.assertEqual(out["with_unique_ofwat_key_and_ea_corroboration"], 1)
        self.assertEqual(out["with_corroboration_only_in_all_not_active"], 1)
        self.assertNotIn("SYNTH", json.dumps(out))

    def test_corroboration_cannot_introduce_an_unreported_pair(self):
        unexpected = {("AN-SYNTH_X", ("AN", "SYNTH/X", "1", "1", "1"))}
        with self.assertRaises(ValueError):
            audit.summarize_candidates({"AN-SYNTH_X"}, set(), unexpected, set())


class DatesTests(unittest.TestCase):
    def test_overlapping_recorded_interval_is_only_a_date_screen(self):
        self.assertEqual(audit.date_overlap_status(dict(EFFECTIVE_DATE="2010-01-01",
                            REVOCATION_DATE="2019-01-01")), "overlaps_by_recorded_dates")
        self.assertEqual(audit.date_overlap_status(dict(EFFECTIVE_DATE=dt.date(2018, 12, 31),
                            REVOCATION_DATE="")), "overlaps_by_recorded_dates")

    def test_revocation_on_first_day_excludes_and_future_effective_excludes(self):
        self.assertEqual(audit.date_overlap_status(dict(EFFECTIVE_DATE="2010-01-01",
                            REVOCATION_DATE="2018-01-01")), "revoked_before_year")
        self.assertEqual(audit.date_overlap_status(dict(EFFECTIVE_DATE="2019-01-01",
                            REVOCATION_DATE="")), "effective_after_year")

    def test_missing_invalid_or_reversed_dates_cannot_confirm_activity(self):
        examples = [({}, "missing_effective_date"),
            ({"EFFECTIVE_DATE":"not a date"}, "invalid_effective_date"),
            ({"EFFECTIVE_DATE":"2010-01-01", "REVOCATION_DATE":"bad"}, "invalid_revocation_date"),
            ({"EFFECTIVE_DATE":"2018-01-02", "REVOCATION_DATE":"2018-01-01"}, "reversed_or_empty_interval")]
        for value, expected in examples:
            with self.subTest(expected=expected):
                self.assertEqual(audit.date_overlap_status(value), expected)


class SourceBoundaryTests(unittest.TestCase):
    def test_joint_corroboration_excludes_conflicting_rows_and_preserves_historical_only(self):
        sources = [row(sample_samplingpoint_notation="AN-SYNTH_" + x,
                       permit_number="SYNTH/" + x) for x in "ABCD"]
        pairs, _ = audit.ofwat_pairs(sources)
        def consent(letter, **changes):
            item = dict(EA_REGION="AN", PERMIT_NUMBER="SYNTH/" + letter,
                        PERMIT_VERSION="1", OUTLET_NUMBER="1", EFFLUENT_NUMBER="1",
                        EFFECTIVE_DATE="2010-01-01", REVOCATION_DATE="",
                        EFF_TYPE_DESCRIPTION=runner.FINAL_TYPE)
            item.update(changes)
            return item
        all_index = runner.consent_index([consent("A"), consent("B"), consent("C"),
                         consent("C", EFFECTIVE_DATE="2019-01-01"), consent("D")])
        active_index = runner.consent_index([consent("B")])
        effluents = [dict(EA_REGION="AN", PERMIT_REF="SYNTH/" + x, VERSION="1",
                          OUTLET_NUMBER="1", EFFLUENT_NUMBER="1", EFF_SAMPLE_POINT="SYNTH_" + x,
                          SPT_DESC=runner.FINAL_TYPE) for x in "ABCD"]
        effluents.append(dict(effluents[-1], EFF_SAMPLE_POINT="SYNTH_CHANGED"))
        confirmed, active, qa = runner.corroborate(pairs, set(), all_index, active_index, effluents)
        self.assertEqual(len(confirmed), 2)
        self.assertEqual(len(active), 1)
        self.assertEqual(qa["ofwat_pair_fates"]["conflicting_required_source_fields"], 1)
        self.assertEqual(qa["ofwat_pair_fates"]["conflicting_Effluents_point_or_type"], 1)
        out = audit.summarize_candidates({"AN-SYNTH_A", "AN-SYNTH_B"}, pairs, confirmed, active)
        self.assertEqual(out["with_corroboration_only_in_all_not_active"], 1)


    def test_repeated_reordered_header_keeps_keys_and_drops_concentration(self):
        header = list(runner.OFWAT_FIELDS) + ["MEAN"]
        first = row()
        values = [first.get(k) for k in header[:-1]] + [98765.4]
        reordered = list(reversed(header))
        second = row(sample_samplingpoint_notation="AN-SYNTH_B")
        other = [second.get(k) if k != "MEAN" else 12345.6 for k in reordered]
        rows, qa = runner.project_sheet([header, values, reordered, other])
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[1]["sample_samplingpoint_notation"], "AN-SYNTH_B")
        self.assertEqual(qa["header_count"], 2)
        self.assertNotIn("MEAN", rows[0])
        self.assertNotIn("98765", json.dumps(rows))

    def test_missing_key_column_is_source_stop(self):
        with self.assertRaises(ValueError):
            runner.project_sheet([["year", "sample_samplingpoint_notation"], [2018, "AN-SYNTH_A"]])

    def test_changed_source_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "source"
            p.write_bytes(b"changed source")
            with self.assertRaises(ValueError):
                runner.verify_file(p, "0" * 64)

    def test_errors_publish_only_controlled_class_and_stage(self):
        def fail(out):
            out["stage"] = "permit_database"
            raise ValueError("private permit SYNTH-PRIVATE and coordinates")
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "result.json"
            self.assertEqual(runner.guarded_run(p, fail), 1)
            value = p.read_text()
            self.assertNotIn("SYNTH-PRIVATE", value)
            self.assertEqual(json.loads(value)["error_stage"], "permit_database")
            self.assertEqual(json.loads(value)["error_type"], "ValueError")

    def test_serialization_failure_is_not_success(self):
        def invalid(out):
            out.update(status=runner.SUCCESS, stage="serialization", invalid=float("nan"))
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder) / "result.json"
            self.assertEqual(runner.guarded_run(p, invalid), 1)
            self.assertEqual(json.loads(p.read_text())["status"], "SOURCE_IDENTITY_AUDIT_STOP")


if __name__ == "__main__":
    unittest.main()
