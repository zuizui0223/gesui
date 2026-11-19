import copy
import csv
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"analysis/20261119_de_sludge_saarland/de_sludge_origin_2016to22_SOURCE_trend.py"
CSV=ROOT/"analysis/20261119_de_sludge_saarland/DE_GERMANY_SAARLAND_2016TO22_SLUDGE_ORIGIN_ALLOCATION_SOURCE.csv"
spec=importlib.util.spec_from_file_location("de_sludge_origin",SCRIPT)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def source_rows():
    with CSV.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

class GermanySludgeRoutesTest(unittest.TestCase):
    def test_seven_annual_source_rows_total_mass_close(self):
        r=m.origin_routing(source_rows())
        self.assertEqual(len(r["annual_rows"]),7)
        self.assertEqual([x["year"] for x in r["annual_rows"]],list(range(2016,2023)))
        self.assertEqual(r["annual_rows"][1]["DE_agricultural_dry_t"],306230)
        self.assertEqual(r["annual_rows"][6]["DE_agricultural_dry_t"],233065)
        self.assertEqual(r["annual_rows"][1]["SL_agricultural_dry_t"],2140)
        self.assertEqual(r["annual_rows"][6]["SL_agricultural_dry_t"],3502)

    def test_national_and_state_originals_have_opposing_descriptive_changes(self):
        r=m.origin_routing(source_rows())
        self.assertAlmostEqual(r["2017to2022"]["DE_agri_percent_change"],-23.892,places=3)
        self.assertAlmostEqual(r["2017to2022"]["SL_agri_percent_change"],63.645,places=3)
        self.assertAlmostEqual(r["2017to2022"]["SL_minus_RestDE_change_in_agri_share_percentage_points_DESCRIPTIVE_ONLY"],12.374,places=3)
        self.assertTrue(r["national_includes_Saarland_avoid_pseudoreplication"])
        self.assertTrue(r["non_overlapping_comparison_uses_Saarland_vs_rest_DE"])

    def test_no_field_exposure_or_biological_inference(self):
        r=m.origin_routing(source_rows())
        self.assertTrue(r["regional_allocation_by_WWTP_origin_not_actually_receiving_farmland"])
        self.assertTrue(r["thermal_share_not_evidence_of_more_aqueous_N_P_discharge"])
        self.assertFalse(r["source_GO_for_biological_model"])
        self.assertTrue(r["no_EQR_NUTRI_HNS_used"])
        self.assertTrue(r["no_hypothesis_test_or_p_values"])

    def test_one_stated_mass_bad_fails_closure(self):
        rows=source_rows()
        rows[2]["SL_agricultural_dry_t"]="3884"
        with self.assertRaises(ValueError):
            m.origin_routing(rows)

    def test_double_counting_state_national_not_allowed(self):
        rows=source_rows()
        rows[2]["SL_total_dry_t"]="1744507"
        with self.assertRaises(ValueError):
            m.origin_routing(rows)

    def test_years_not_shifted(self):
        rows=source_rows()
        rows[1]["year"]="2016"
        with self.assertRaises(ValueError):
            m.origin_routing(rows)

    def test_missing_measurements_not_zero(self):
        rows=source_rows()
        rows[2]["SL_agricultural_dry_t"]=""
        with self.assertRaises(ValueError):
            m.origin_routing(rows)

    def test_negative_category_rejected(self):
        rows=source_rows()
        rows[1]["DE_other_disposal_dry_t"]="-1"
        with self.assertRaises(ValueError):
            m.origin_routing(rows)

if __name__=="__main__":
    unittest.main()
