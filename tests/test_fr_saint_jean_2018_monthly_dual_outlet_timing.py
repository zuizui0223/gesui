import copy
import csv
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
P=ROOT/"analysis/20261119_fr_2wwtp_dualoutlet/fr_saint_jean_2018_monthly_dual_outlet_timing_gate.py"
CSV=ROOT/"analysis/20261119_fr_2wwtp_dualoutlet/FR_SAINT_JEAN_2018_MONTHLY_TREATED_A4_AND_REALIZED_SLUDGE_EVENT_COUNTS.csv"
spec=importlib.util.spec_from_file_location("sj2018",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def rows():
    with CSV.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

class SJMonthlyWaterLandTests(unittest.TestCase):
    def test_original_annual_A4_water_and_real_spread_months_reconcile(self):
        r=m.audit(rows())
        self.assertEqual(r["annual_A4_treated_effluent_m3_reported"],343221)
        self.assertEqual(r["positive_A4_treated_outlet_flow_months"],12)
        self.assertEqual(r["months_with_realized_sludge_application_logs"],["03","09"])
        self.assertEqual(r["realized_agricultural_sludge_ledger_events"],8)
        self.assertEqual(r["water_A4_volume_m3_in_months_with_sludge_events"],58199)
        self.assertAlmostEqual(r["fraction_annual_A4_treated_water_in_sludge_event_months_DESCRIPTIVE_ONLY"],58199/343221,places=5)
        self.assertEqual(r["independent_station_years"],1)

    def test_not_a_causal_aqueous_terrestrial_nutrient_trajectory(self):
        r=m.audit(rows())
        self.assertTrue(r["monthly_actual_sludge_N_P_mass_unidentified"])
        self.assertTrue(r["monthly_aqueous_NGL_PT_export_mass_unidentified"])
        self.assertTrue(r["A5_bypass_meter_malfunction_preserved"])
        self.assertTrue(r["12_months_are_not_independent_stations"])
        self.assertEqual(r["monthly_causal_coupled_nutrient_transfer_gate"],"STOP")
        self.assertFalse(r["no_EQR_NUTRI_HNS_values_read"] is False)

    def test_wrong_month_flow_fails(self):
        d=rows()
        d[2]["reported_treated_A4_volume_m3"]="32801"
        with self.assertRaises(ValueError):m.audit(d)

    def test_shifted_sludge_event_month_fails(self):
        d=rows()
        d[2]["actual_sludge_application_log_events_in_month"]="0"
        d[3]["actual_sludge_application_log_events_in_month"]="5"
        with self.assertRaises(ValueError):m.audit(d)

    def test_annual_mass_not_in_source(self):
        d=rows()
        d[8]["monthly_outlet_TN_mass_measured_kg"]="150"
        with self.assertRaises(ValueError):m.audit(d)

    def test_bypass_not_confirmed_zero(self):
        d=rows()
        d[0]["A5_bypass_volume_reliable"]="true"
        with self.assertRaises(ValueError):m.audit(d)

    def test_not_independent_months(self):
        d=rows()
        d.pop()
        with self.assertRaises(ValueError):m.audit(d)

    def test_dont_impute_monthly_sludge_mass_from_annual(self):
        d=rows()
        d[2]["actual_sludge_wet_mass_by_month_t"]="845"
        with self.assertRaises(ValueError):m.audit(d)

if __name__=="__main__":
    unittest.main()
