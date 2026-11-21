"""Hand-checked tests for aggregate-only candidate/TP-availability auditing."""
from __future__ import annotations

import copy
import json
import unittest

try:
    from premoval_source_support import build_source_support
except ModuleNotFoundError as exc:
    if exc.name != "premoval_source_support":
        raise
    build_source_support = None


class MembershipOnlyTPs(dict):
    """Values exist, but the audit may inspect only keys and membership."""

    def __getitem__(self, key):
        raise AssertionError("TP payload must not be read")

    def get(self, *args, **kwargs):
        raise AssertionError("TP payload must not be read")

    def values(self):
        raise AssertionError("TP payload must not be read")

    def items(self):
        raise AssertionError("TP payload must not be read")


def fixture():
    # A/B: <2000; C/D: 2000..9999; E/F: 10000..99999;
    # G/H: >=100000; I/J: missing/nonpositive. Each pair has P=0/P=1.
    # Eligible TP at A, B, D, E, H: P=0 has 2/5; P=1 has 3/5.
    rows = [
        ("A", 0, 0, "UKC", 1000, 0, "eligible_numeric_TP"),
        ("B", 1, 0, "UKC", 1999, 10, "eligible_numeric_TP"),
        ("C", 0, 1, "UKC", 2000, 20, "no_eligible_numeric_TP"),
        ("D", 1, 1, "UKC", 9999, 30, "eligible_numeric_TP"),
        ("E", 0, 0, "UKD", 10000, 40, "eligible_numeric_TP"),
        ("F", 1, 0, "UKD", 99999, 50, "no_returned_TP_record"),
        ("G", 0, 0, "UKE", 100000, 60, "no_returned_TP_record"),
        ("H", 1, 1, "UKE", 200000, 70, "eligible_numeric_TP"),
        ("I", 0, 0, None, None, 80, "no_eligible_numeric_TP"),
        ("J", 1, 1, None, 0, 90, "no_returned_TP_record"),
    ]
    matches, plants, tps, fates = {}, {}, MembershipOnlyTPs(), {}
    for index, (letter, pflag, nflag, region, capacity, distance, fate) in enumerate(rows):
        sa, plant = "private_sa_" + letter, "private_plant_" + letter
        matches[sa] = {
            "outlet_id": ("UK5", plant, "private_outlet_" + letter),
            "outlet_coord": (400001.125 + index, 300001.875 + index),
            "original2018_FW_outlet_to_SA_distance_m": distance,
        }
        plants[plant] = {
            "PRemoval": pflag, "NRemoval": nflag,
            "capacity": capacity, "original_region": region,
            "name": "private_facility_name_" + letter,
        }
        fates[sa] = {"fate": fate, "returned_source_rows": 0}
        if fate == "eligible_numeric_TP":
            tps[sa] = {"annual_median_TP_mgL": 987654.321}
    return matches, plants, tps, fates


class SourceSupportTests(unittest.TestCase):
    def setUp(self):
        self.assertTrue(callable(build_source_support), "source-support helper not implemented")
        self.matches, self.plants, self.tps, self.fates = fixture()

    def audit(self, with_fates=False):
        return build_source_support(
            self.matches, self.plants, self.tps,
            self.fates if with_fates else None,
        )

    def test_availability_counts_use_all_candidates_as_denominator(self):
        result = self.audit()
        self.assertEqual(result["totals"], {
            "candidate_n": 10, "eligible_numeric_TP_n": 5,
            "not_eligible_numeric_TP_n": 5, "eligible_numeric_TP_rate": 0.5,
        })
        self.assertEqual(result["by_treatment_flag"]["PRemoval"]["by_flag"], {
            "0": {"candidate_n": 5, "eligible_numeric_TP_n": 2,
                  "not_eligible_numeric_TP_n": 3, "eligible_numeric_TP_rate": 0.4},
            "1": {"candidate_n": 5, "eligible_numeric_TP_n": 3,
                  "not_eligible_numeric_TP_n": 2, "eligible_numeric_TP_rate": 0.6},
        })
        self.assertEqual(result["by_treatment_flag"]["NRemoval"]["by_flag"], {
            "0": {"candidate_n": 6, "eligible_numeric_TP_n": 3,
                  "not_eligible_numeric_TP_n": 3, "eligible_numeric_TP_rate": 0.5},
            "1": {"candidate_n": 4, "eligible_numeric_TP_n": 2,
                  "not_eligible_numeric_TP_n": 2, "eligible_numeric_TP_rate": 0.5},
        })

    def test_region_and_joint_tables_preserve_missing_and_unobserved_cells(self):
        group = self.audit()["by_treatment_flag"]["PRemoval"]
        self.assertEqual(group["by_flag_and_NUTS1"]["UKC"], {
            "0": {"candidate_n": 2, "eligible_numeric_TP_n": 1,
                  "not_eligible_numeric_TP_n": 1, "eligible_numeric_TP_rate": 0.5},
            "1": {"candidate_n": 2, "eligible_numeric_TP_n": 2,
                  "not_eligible_numeric_TP_n": 0, "eligible_numeric_TP_rate": 1.0},
        })
        joint = group["by_flag_and_NUTS1_and_capacity_band"]
        self.assertEqual(joint["UKD"]["10000_to_99999"]["1"], {
            "candidate_n": 1, "eligible_numeric_TP_n": 0,
            "not_eligible_numeric_TP_n": 1, "eligible_numeric_TP_rate": 0.0,
        })
        self.assertEqual(joint["missing"]["missing_or_nonpositive"]["0"]["candidate_n"], 1)
        self.assertEqual(joint["missing"]["missing_or_nonpositive"]["1"]["candidate_n"], 1)
        self.assertIsNone(joint["UKC"]["100000plus"]["0"]["eligible_numeric_TP_rate"])

    def test_capacity_boundaries_have_fixed_bands(self):
        table = self.audit()["by_treatment_flag"]["PRemoval"]["by_flag_and_capacity_band"]
        observed = {
            band: (row["0"]["candidate_n"], row["1"]["candidate_n"],
                   row["0"]["eligible_numeric_TP_n"], row["1"]["eligible_numeric_TP_n"])
            for band, row in table.items()
        }
        self.assertEqual(observed, {
            "<2000": (1, 1, 1, 1), "2000_to_9999": (1, 1, 0, 1),
            "10000_to_99999": (1, 1, 1, 0), "100000plus": (1, 1, 0, 1),
            "missing_or_nonpositive": (1, 1, 0, 0),
        })

    def test_negative_and_nonfinite_capacity_remain_missing(self):
        for capacity in (-1, float("nan"), float("inf")):
            with self.subTest(capacity=capacity):
                self.plants["private_plant_J"]["capacity"] = capacity
                group = self.audit()["by_treatment_flag"]["PRemoval"]
                self.assertEqual(group["by_flag_and_capacity_band"]["missing_or_nonpositive"]["1"]["candidate_n"], 1)
                self.assertEqual(group["capacity_PE_summary"]["candidates"]["1"]["missing_n"], 1)

    def test_support_is_separate_for_candidates_and_observed_sites(self):
        support = self.audit()["by_treatment_flag"]["PRemoval"]["both_flag_support"]
        self.assertEqual(support["NUTS1"]["candidates"], {
            "cells_with_any_site": 4, "cells_with_both_flags": 4,
            "covered_n": 10, "lacking_support_n": 0,
            "covered_by_flag": {"0": 5, "1": 5},
            "lacking_support_by_flag": {"0": 0, "1": 0},
        })
        self.assertEqual(support["NUTS1"]["eligible_numeric_TP"], {
            "cells_with_any_site": 3, "cells_with_both_flags": 1,
            "covered_n": 3, "lacking_support_n": 2,
            "covered_by_flag": {"0": 1, "1": 2},
            "lacking_support_by_flag": {"0": 1, "1": 1},
        })
        self.assertEqual(support["NUTS1_x_capacity_band"]["candidates"]["cells_with_both_flags"], 5)
        self.assertEqual(support["capacity_band"]["eligible_numeric_TP"], {
            "cells_with_any_site": 4, "cells_with_both_flags": 1,
            "covered_n": 2, "lacking_support_n": 3,
            "covered_by_flag": {"0": 1, "1": 1},
            "lacking_support_by_flag": {"0": 1, "1": 2},
        })
        self.assertEqual(support["NUTS1_x_capacity_band"]["eligible_numeric_TP"],
                         support["capacity_band"]["eligible_numeric_TP"])
        nitrogen = self.audit()["by_treatment_flag"]["NRemoval"]["both_flag_support"]
        self.assertEqual(nitrogen["capacity_band"]["candidates"]["covered_n"], 4)
        self.assertEqual(nitrogen["capacity_band"]["eligible_numeric_TP"]["covered_n"], 0)
        self.assertEqual(nitrogen["capacity_band"]["eligible_numeric_TP"]["lacking_support_n"], 5)

    def test_summaries_use_valid_capacity_without_zero_imputation(self):
        group = self.audit()["by_treatment_flag"]["PRemoval"]
        self.assertEqual(group["capacity_PE_summary"]["candidates"], {
            "0": {"count": 4, "missing_n": 1, "min": 1000.0, "median": 6000.0, "max": 100000.0},
            "1": {"count": 4, "missing_n": 1, "min": 1999.0, "median": 54999.0, "max": 200000.0},
        })
        self.assertEqual(group["capacity_PE_summary"]["eligible_numeric_TP"]["0"], {
            "count": 2, "missing_n": 0, "min": 1000.0, "median": 5500.0, "max": 10000.0,
        })
        self.assertEqual(group["outlet_to_SA_distance_m_summary"]["candidates"]["1"], {
            "count": 5, "missing_n": 0, "min": 10.0, "median": 50.0, "max": 90.0,
        })
        self.assertEqual(group["outlet_to_SA_distance_m_summary"]["eligible_numeric_TP"]["1"], {
            "count": 3, "missing_n": 0, "min": 10.0, "median": 30.0, "max": 70.0,
        })

    def test_optional_fates_are_aggregated_by_flag(self):
        result = self.audit(with_fates=True)
        self.assertEqual(result["by_treatment_flag"]["PRemoval"]["site_fates_by_flag"], {
            "0": {"no_returned_TP_record": 1, "no_eligible_numeric_TP": 2, "eligible_numeric_TP": 2},
            "1": {"no_returned_TP_record": 2, "no_eligible_numeric_TP": 0, "eligible_numeric_TP": 3},
        })
        self.assertNotIn("site_fates_by_flag", self.audit()["by_treatment_flag"]["PRemoval"])

    def test_tp_payloads_and_private_source_fields_do_not_escape(self):
        result = self.audit(with_fates=True)  # MembershipOnlyTPs rejects every payload read.
        serialized = json.dumps(result, sort_keys=True, allow_nan=False)
        for forbidden in ("private_sa_", "private_plant_", "private_outlet_",
                          "private_facility_name_", "400001.125", "300001.875", "987654.321"):
            self.assertNotIn(forbidden, serialized)
        self.assertEqual(result["source_key_counts"], {
            "matches": 10, "plants": 10, "tps": 5, "site_fates": 10,
            "candidate_native_plants": 10, "candidate_native_outlets": 10,
        })

    def test_mapping_hash_is_order_independent_and_tracks_identity_and_coordinates(self):
        original = self.audit()["pair_mapping_sha256"]
        self.assertRegex(original, r"^[0-9a-f]{64}$")
        self.matches = dict(reversed(list(self.matches.items())))
        self.assertEqual(self.audit()["pair_mapping_sha256"], original)
        self.tps = MembershipOnlyTPs()
        self.assertEqual(self.audit()["pair_mapping_sha256"], original)
        self.matches["private_sa_A"]["outlet_coord"] = (400001.375, 300001.875)
        self.assertNotEqual(self.audit()["pair_mapping_sha256"], original)
        self.matches, self.plants, self.tps, self.fates = fixture()
        self.matches["private_sa_A"], self.matches["private_sa_B"] = (
            self.matches["private_sa_B"], self.matches["private_sa_A"])
        self.assertNotEqual(self.audit()["pair_mapping_sha256"], original)

    def test_rejects_observed_sites_outside_candidate_mapping(self):
        self.tps["outside_candidate"] = object()
        with self.assertRaisesRegex(ValueError, "outside candidate"):
            self.audit()

    def test_rejects_reused_native_outlet(self):
        self.matches["private_sa_B"]["outlet_id"] = self.matches["private_sa_A"]["outlet_id"]
        with self.assertRaisesRegex(ValueError, "duplicate native outlet"):
            self.audit()

    def test_rejects_two_outlets_of_one_native_plant(self):
        self.matches["private_sa_B"]["outlet_id"] = ("UK5", "private_plant_A", "second_outlet")
        with self.assertRaisesRegex(ValueError, "duplicate native plant"):
            self.audit()

    def test_rejects_missing_native_plant_and_unknown_flags(self):
        del self.plants["private_plant_A"]
        with self.assertRaisesRegex(ValueError, "missing native plant"):
            self.audit()
        self.matches, self.plants, self.tps, self.fates = fixture()
        for field, value in (("PRemoval", 2), ("NRemoval", None), ("PRemoval", "1")):
            with self.subTest(field=field, value=value):
                plants = copy.deepcopy(self.plants)
                plants["private_plant_A"][field] = value
                with self.assertRaisesRegex(ValueError, "unrecognized treatment flag"):
                    build_source_support(self.matches, plants, self.tps)

    def test_rejects_incomplete_or_inconsistent_fates(self):
        del self.fates["private_sa_J"]
        with self.assertRaisesRegex(ValueError, "fate keys"):
            self.audit(with_fates=True)
        self.matches, self.plants, self.tps, self.fates = fixture()
        self.fates["private_sa_A"]["fate"] = "no_returned_TP_record"
        with self.assertRaisesRegex(ValueError, "fate.*membership"):
            self.audit(with_fates=True)
        self.fates["private_sa_A"]["fate"] = "unknown"
        with self.assertRaisesRegex(ValueError, "unrecognized site fate"):
            self.audit(with_fates=True)

    def test_rejects_nonfinite_or_negative_distance_and_nonfinite_coordinates(self):
        for distance in (-1, float("nan"), float("inf")):
            with self.subTest(distance=distance):
                matches = copy.deepcopy(self.matches)
                matches["private_sa_A"]["original2018_FW_outlet_to_SA_distance_m"] = distance
                with self.assertRaisesRegex(ValueError, "match distance"):
                    build_source_support(matches, self.plants, self.tps)
        self.matches["private_sa_A"]["outlet_coord"] = (float("nan"), 300001.875)
        with self.assertRaisesRegex(ValueError, "outlet coordinates"):
            self.audit()

    def test_rejects_invalid_region_instead_of_emitting_a_source_label(self):
        self.plants["private_plant_A"]["original_region"] = "private_free_text"
        with self.assertRaisesRegex(ValueError, "NUTS1"):
            self.audit()

    def test_empty_input_reports_undefined_rates_and_empty_summaries(self):
        result = build_source_support({}, {}, MembershipOnlyTPs())
        self.assertEqual(result["totals"], {
            "candidate_n": 0, "eligible_numeric_TP_n": 0,
            "not_eligible_numeric_TP_n": 0, "eligible_numeric_TP_rate": None,
        })
        group = result["by_treatment_flag"]["PRemoval"]
        self.assertEqual(group["capacity_PE_summary"]["candidates"]["0"], {
            "count": 0, "missing_n": 0, "min": None, "median": None, "max": None,
        })
        self.assertEqual(group["both_flag_support"]["capacity_band"]["candidates"]["cells_with_any_site"], 0)
        json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
