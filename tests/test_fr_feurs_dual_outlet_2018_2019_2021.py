import copy
import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "analysis/20261119_fr_sludge_routeC/fr_feurs_dual_outlet_source_mass_gate.py"
SOURCE = ROOT / "analysis/20261119_fr_sludge_routeC/FR_FEURS_PRIMARY_SOURCE_YEARS_2018_2019_2021.json"
spec = importlib.util.spec_from_file_location("feurs_source", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class FeursDualOutletSourceGateTest(unittest.TestCase):
    def setUp(self):
        self.source = json.loads(SOURCE.read_text(encoding="utf-8"))

    def test_feurs_2018_field_nutrient_proxies_and_P_oxide_conversion(self):
        result = m.audit(self.source)
        r = result["source_year_records"]["2018"]
        self.assertAlmostEqual(r["sludge_N_total_assay_mean_times_applied_dry_mass_t_SCREENING_ONLY"], 22.910356)
        self.assertAlmostEqual(r["sludge_P2O5_assay_mean_times_applied_dry_mass_t_SCREENING_ONLY"], 14.246776)
        self.assertAlmostEqual(r["sludge_elemental_P_equivalent_t_SCREENING_ONLY"], 14.246776 * m.P_OVER_P2O5, places=5)
        self.assertAlmostEqual(r["independent_sludge_N_area_times_mean_reported_units_t_SCREENING_ONLY"], 23.58338)
        self.assertTrue(r["P2O5_not_elemental_P"])
        self.assertTrue(r["do_not_interpret_as_measured_total_mass_or_plant_uptake"])

    def test_2019_2021_water_screening_indices_are_not_measured_loads(self):
        r = m.audit(self.source)["source_year_records"]
        self.assertAlmostEqual(r["2019"]["unweighted_mean_NTK_x_annual_flow_t_SCREENING_INDEX_NOT_N_LOAD"], 5.569174, places=5)
        self.assertAlmostEqual(r["2019"]["unweighted_mean_PT_x_annual_flow_t_SCREENING_INDEX_NOT_P_LOAD"], 0.486627, places=5)
        self.assertAlmostEqual(r["2021"]["unweighted_mean_NTK_x_annual_flow_t_SCREENING_INDEX_NOT_N_LOAD"], 8.966588, places=5)
        self.assertAlmostEqual(r["2021"]["unweighted_mean_PT_x_annual_flow_t_SCREENING_INDEX_NOT_P_LOAD"], 0.706458, places=5)
        for year in ("2019", "2021"):
            self.assertFalse(r[year]["sludge_actual_spreading_confirmed"])
            self.assertTrue(r[year]["NTK_is_not_total_N"])
            self.assertTrue(r[year]["unweighted_mean_not_flow_weighted"])

    def test_all_years_fail_true_paired_mass_and_biological_gates(self):
        r = m.audit(self.source)
        self.assertEqual(r["source_years_with_measured_both_aqueous_and_field_mass"], [])
        self.assertEqual(r["strict_paired_export_mass_GATE"], "STOP")
        self.assertEqual(r["strict_biological_path_GATE"], "STOP")
        self.assertTrue(r["nonmatching_2018_2019_2021_years_never_combined"])
        for record in r["source_year_records"].values():
            self.assertFalse(record["true_same_year_water_and_field_nutrient_mass_pair"])
            self.assertTrue(record["source_only_NO_biological_values_read"])

    def test_claimed_measured_2018_water_total_N_is_rejected(self):
        s = copy.deepcopy(self.source)
        s["years"][0]["true_measured_annual_effluent_total_N_t"] = 1
        with self.assertRaises(ValueError):
            m.audit(s)

    def test_claimed_2018_complete_pair_is_rejected(self):
        s = copy.deepcopy(self.source)
        s["years"][0]["same_year_outlet_mass_pair_observable"] = True
        with self.assertRaises(ValueError):
            m.audit(s)

    def test_fabricated_2019_land_application_is_rejected(self):
        s = copy.deepcopy(self.source)
        s["years"][1]["actual_sludge_spreading_confirmed"] = True
        with self.assertRaises(ValueError):
            m.audit(s)

    def test_2019_mean_must_not_be_declared_flow_weighted(self):
        s = copy.deepcopy(self.source)
        s["years"][1]["effluent_mean_is_flow_weighted"] = True
        with self.assertRaises(ValueError):
            m.audit(s)

    def test_station_change_is_rejected(self):
        s = copy.deepcopy(self.source)
        s["STEU_SANDRE"] = "fictional-id"
        with self.assertRaises(ValueError):
            m.audit(s)


if __name__ == "__main__":
    unittest.main()
