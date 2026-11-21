"""Synthetic geometry/source regressions; all identifiers are invented."""
import json
import math
import unittest

import permit_outlet_location_support as loc
from permit_sample_identity_support import permit_key


def row(**updates):
    value = dict(EA_REGION="AN", PERMIT_NUMBER="SYNTHETIC-PERMIT",
                 PERMIT_VERSION="01", OUTLET_NUMBER="A", EFFLUENT_NUMBER="1",
                 OUTLET_GRID_REF="SU000000", DISCHARGE_NGR="SU100100",
                 EFFLUENT_GRID_REF="SU200200")
    value.update(updates)
    return value


class GridTests(unittest.TestCase):
    def test_tile_and_axis_are_not_swapped(self):
        self.assertEqual(loc.grid_cell("SU1234512345"), (412345, 112345, 1))
        self.assertEqual(loc.grid_cell("NZ123456"), (412300, 545600, 100))
        self.assertEqual(loc.grid_cell("TQ000000"), (500000, 100000, 100))
        self.assertEqual(loc.grid_cell("SV00"), (0, 0, 10000))

    def test_spacing_case_and_precision_are_retained(self):
        self.assertEqual(loc.grid_cell(" su 1234 1234 "), (412340, 112340, 10))
        self.assertEqual(loc.grid_cell("SU123456"), (412300, 145600, 100))

    def test_bad_or_out_of_domain_references_rejected(self):
        for raw in (None, "", "SU", "SU12345", "SI123456", "AA000000", "ZZ123456",
                    "SU123456789012", "SU12/34", "500000,100000"):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                loc.grid_cell(raw)

    def test_distance_encloses_entire_cell(self):
        result = loc.distance_bounds((900, 1000), (1000, 1000, 100))
        self.assertIsInstance(result, tuple)
        self.assertEqual(result[0], 100)
        self.assertAlmostEqual(result[1], math.sqrt(50000))
        self.assertEqual(loc.distance_bounds((1050, 1050), (1000, 1000, 100))[0], 0)

    def test_coarse_grid_boundary_not_falsely_precise(self):
        self.assertEqual(loc.distance_class((900, 1000), (1000, 1000, 100), 250), "within")
        self.assertEqual(loc.distance_class((900, 1000), (1000, 1000, 100), 100), "boundary")
        self.assertEqual(loc.distance_class((900, 1000), (1000, 1000, 100), 99), "outside")

    def test_nonfinite_geometry_and_invalid_radius_rejected(self):
        for point, cell in [((float("nan"), 1), (1, 1, 1)), ((1, 1), (1, 1, 0))]:
            with self.assertRaises(ValueError):
                loc.distance_bounds(point, cell)
        with self.assertRaises(ValueError):
            loc.distance_class((0, 0), (1, 1, 1), -1)


class LocationTests(unittest.TestCase):
    def index(self, rows, eligible):
        result = loc.location_index(rows, eligible)
        self.assertIsInstance(result, tuple)
        self.assertEqual(len(result), 2)
        return result

    def test_duplicate_and_whitespace_do_not_inflate_geometry(self):
        a = row(); key = permit_key(a, "consents")
        index, qa = self.index([a, row(OUTLET_GRID_REF=" su 000 000 ")], {key})
        self.assertEqual(index[key]["OUTLET_GRID_REF"], {"status": "valid", "cell": (400000, 100000, 100)})
        self.assertEqual(qa["eligible_distinct_keys"], 1)
        self.assertEqual(qa["versioned_outlet_status_counts"], {"valid": 1})

    def test_missing_outlet_is_not_filled_from_site_or_effluent(self):
        a = row(OUTLET_GRID_REF=""); key = permit_key(a, "consents")
        index, _ = self.index([a], {key})
        self.assertEqual(index[key]["OUTLET_GRID_REF"], {"status": "missing"})
        self.assertEqual(index[key]["DISCHARGE_NGR"]["status"], "valid")

    def test_invalid_outlet_is_separate_from_missing(self):
        a = row(OUTLET_GRID_REF="XX?"); key = permit_key(a, "consents")
        index, _ = self.index([a], {key})
        self.assertEqual(index[key]["OUTLET_GRID_REF"], {"status": "invalid"})

    def test_missing_row_not_silently_skipped(self):
        a = row(); key = permit_key(a, "consents")
        index, _ = self.index([], {key})
        self.assertEqual(index[key]["OUTLET_GRID_REF"], {"status": "absent_key"})

    def test_same_key_conflict_does_not_select_first(self):
        a = row(); key = permit_key(a, "consents")
        for second in (row(OUTLET_GRID_REF="SU001001"), row(OUTLET_GRID_REF="")):
            index, _ = self.index([a, second], {key})
            self.assertEqual(index[key]["OUTLET_GRID_REF"]["status"], "conflicting")

    def test_other_effluent_on_same_versioned_outlet_checked(self):
        a = row(); key = permit_key(a, "consents")
        other = row(EFFLUENT_NUMBER="2", OUTLET_GRID_REF="SU001001")
        index, qa = self.index([a, other], {key})
        self.assertEqual(index[key]["OUTLET_GRID_REF"]["status"], "conflicting_versioned_outlet")
        self.assertEqual(qa["versioned_outlet_status_counts"], {"conflicting": 1})

    def test_different_version_region_or_outlet_not_collapsed(self):
        a = row(); key = permit_key(a, "consents")
        others = [row(PERMIT_VERSION="2", OUTLET_GRID_REF="SU001001"),
                  row(EA_REGION="NE", OUTLET_GRID_REF="SU001001"),
                  row(OUTLET_NUMBER="B", OUTLET_GRID_REF="SU001001")]
        index, _ = self.index([a] + others, {key})
        self.assertEqual(index[key]["OUTLET_GRID_REF"]["status"], "valid")

    def test_summary_uses_outlet_and_preserves_missing_denominator(self):
        a = row(OUTLET_GRID_REF="SU0000000000")
        b = row(PERMIT_NUMBER="SYNTHETIC-MISSING", OUTLET_GRID_REF="")
        ka, kb = permit_key(a, "consents"), permit_key(b, "consents")
        index, _ = self.index([a, b], {ka, kb})
        pairs = {("AN-SYNTHETIC1", ka), ("AN-SYNTHETIC2", kb)}
        matching = {p: {"outlet_coord": (400000, 100000)} for p, _ in pairs}
        registry = {p: (400600, 100000) for p, _ in pairs}
        result = loc.summarize_locations(pairs, index, matching, registry)
        self.assertIsInstance(result, dict)
        self.assertEqual(result["candidate_points"], 2)
        self.assertEqual(result["outlet_status_counts"], {"missing": 1, "valid": 1})
        self.assertEqual(result["distances"]["original_EEA_outlet"]["250m"], {"within": 1, "outside": 0, "boundary": 0})
        self.assertEqual(result["distances"]["current_WQE_samplepoint"]["250m"], {"within": 0, "outside": 1, "boundary": 0})
        self.assertEqual(result["outlet_vs_site_grid"]["different_cells"], 1)
        self.assertNotIn("SYNTHETIC", json.dumps(result))
        self.assertNotIn("400000", json.dumps(result))

    def test_missing_matching_point_stops(self):
        a = row(); key = permit_key(a, "consents")
        index, _ = self.index([a], {key})
        with self.assertRaises(ValueError):
            loc.summarize_locations({("AN-SYNTHETIC", key)}, index, {}, {})

    def test_multiple_keys_for_one_point_stops(self):
        a, b = row(), row(OUTLET_NUMBER="B")
        ka, kb = permit_key(a, "consents"), permit_key(b, "consents")
        index, _ = self.index([a, b], {ka, kb})
        with self.assertRaises(ValueError):
            loc.summarize_locations({("AN-SYNTHETIC", ka), ("AN-SYNTHETIC", kb)}, index,
                                    {"AN-SYNTHETIC": {"outlet_coord": (400000, 100000)}},
                                    {"AN-SYNTHETIC": (400000, 100000)})


if __name__ == "__main__":
    unittest.main()
