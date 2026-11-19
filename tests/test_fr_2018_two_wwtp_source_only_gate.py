import copy
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"analysis/20261119_fr_2wwtp_dualoutlet"
P=DIR/"fr_2018_two_wwtp_source_only_gate.py"
spec=importlib.util.spec_from_file_location("fr2wwtp",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class FranceTwoIndependentWWTP2018Tests(unittest.TestCase):
    def setUp(self):
        self.main=m.getcsv(DIR/"FR_2018_TWO_WWTP_SOURCE_ROWS.csv",m.MAIN_COLS)
        self.events=m.getcsv(
            DIR/"FR_SAINT_JEAN_ANGELY_2018_ANONYMIZED_EVENT_MONTH_COMMUNE_COUNTS.csv",
            m.EVENT_COLS)

    def test_two_independent_wwtp_one_year(self):
        r=m.audit(self.main,self.events)
        self.assertEqual(r["independent_identified_STEU_year_units"],2)
        self.assertEqual(r["distinct_receiving_rivers_from_original_sources"],2)
        self.assertEqual(r["2018_realized_plant_application_ledger_rows"],51)
        self.assertAlmostEqual(r["land_actual_area_ha_sum_TWO_SOURCES_DESCRIPTIVE_ONLY"],235.73)
        self.assertEqual(r["land_actual_wet_t_sum_TWO_SOURCES_DESCRIPTIVE_ONLY"],2763.0)
        self.assertEqual(r["independent_H6_assigned_and_audited"],0)

    def test_water_A4_same_year_one_station_not_complete_discharge(self):
        r=m.audit(self.main,self.events)
        self.assertEqual(r["one_plant_same_year_A4_treated_outlet_m3_reported"],1)
        self.assertEqual(r["Saint_Jean_2018_A4_m3"],343221)
        self.assertTrue(r["Saint_Jean_2018_by_pass_A5_instrumentation_unreliable"])
        self.assertEqual(r["2018_plant_years_with_both_measured_aqueous_TN_TP_and_realized_land_N_P"],0)
        sj=[x for x in r["source_rows"] if x["station_SANDRE"]=="0517347V004"][0]
        self.assertTrue(sj["water_A4_volume_is_not_confirmed_total_received_with_unreliable_A5"])

    def test_actual_Saint_Jean_ledger_counts_only_not_anonymous_fake_masses(self):
        r=m.audit(self.main,self.events)
        self.assertEqual(r["Saint_Jean_month_register_rows"],{"03":5,"09":3})
        self.assertTrue(r["sludge_total_dry_mass_across_STEUs_liming_mismatch_not_pooled"])

    def test_missing_measured_aqueous_and_field_NP_are_not_zero(self):
        r=m.audit(self.main,self.events)
        self.assertFalse(r["NUTRI_HNS_EQR_values_read"])
        self.assertEqual(r["actual_field_geometries_verified_and_linked"],0)
        self.assertEqual(r["coupled_N_P_mass_balance_gate"],"STOP")
        self.assertEqual(r["biotic_response_and_causal_path_gate"],"STOP")

    def test_wrongly_treating_influent_NGL_kg_as_effluent_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["annual_aqueous_total_N_measured_t"]="41.013"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_wrongly_treating_influent_PT_kg_as_effluent_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["annual_aqueous_total_P_measured_t"]="5.331"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_wrongly_making_flows_total_or_missing_by_pass_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["water_bypass_volume_reliability"]="confirmed_zero_A5"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_wrongly_substituting_other_year_A4_flow_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["water_A4_treated_outlet_m3"]="401689"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_liming_prevents_dry_mass_fake_equivalence(self):
        main=copy.deepcopy(self.main)
        main[0]["land_dry_matter_lime_status"]="unlimed"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_pseudoreplicated_station_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["station_SANDRE"]="0442094S0002"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_non2018_application_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["year"]="2019"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_private_farm_field_rejected(self):
        main=copy.deepcopy(self.main)
        main[1]["farmer_name"]="not allowed"
        with self.assertRaises(ValueError):m.audit(main,self.events)

    def test_event_group_year_and_count_tamper_rejected(self):
        rows=copy.deepcopy(self.events)
        rows[0]["register_event_rows"]="3"
        with self.assertRaises(ValueError):m.audit(self.main,rows)
        rows=copy.deepcopy(self.events)
        rows[0]["year"]="2021"
        with self.assertRaises(ValueError):m.audit(self.main,rows)


if __name__=="__main__":
    unittest.main()
