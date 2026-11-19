import copy
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
DIR=ROOT/"analysis/20261119_fr_siaap_2022"
P=DIR/"fr_siaap_2022_real_spreading_SOURCE_audit.py"
spec=importlib.util.spec_from_file_location("siaap",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class SIAAPSeineAval2022SourceTest(unittest.TestCase):
    def setUp(self):
        self.depts=m.csvread(DIR/"FR_SIAAP_SEINE_AVAL_2022_REALIZED_SLUDGE_BY_DEPARTMENT.csv",m.DEPT_HEADER)
        self.products=m.csvread(DIR/"FR_SIAAP_SEINE_AVAL_2022_REALIZED_SLUDGE_BY_PRODUCT.csv",m.PRODUCT_HEADER)

    def test_two_independent_original_2022_reconciliation_tables(self):
        r=m.audit(self.depts,self.products)
        self.assertTrue(r["independent_geographic_subtotals_check_pass"])
        self.assertTrue(r["independent_product_subtotals_check_pass"])
        a=r["actual_2022_sludge_spreading_recap"]
        self.assertAlmostEqual(a["applied_wet_t"],81195.15)
        self.assertAlmostEqual(a["applied_dry_t"],35102.66)
        self.assertAlmostEqual(a["actually_spread_area_ha"],10860.49)
        self.assertEqual(a["original_field_parcel_count"],1092)
        self.assertEqual(a["department_count"],13)
        self.assertEqual(a["sludge_originating_independent_stations"],1)

    def test_actual_application_not_delivery_or_2021_plant_exposure(self):
        r=m.audit(self.depts,self.products)
        self.assertEqual(r["official_agricultural_campaign_year"],2022)
        self.assertEqual(r["source_delivery_window_includes_prior_2021"],"2021-11-08 to 2022-09-30")
        self.assertTrue(r["sludge_delivered_in_2021_not_necessarily_applied_in_2021"])
        self.assertEqual(len(r["realized_field_spreading_windows_2022"]),2)

    def test_field_and_river_2022_NP_and_biology_stopped(self):
        r=m.audit(self.depts,self.products)
        self.assertTrue(r["SANDRE_year2022_direct_original_code_not_independently_fixed"])
        self.assertTrue(r["actual_parcel_geometries_not_downloaded"])
        self.assertTrue(r["2022_aqueous_annual_measured_TN_TP_export_not_identified"])
        self.assertEqual(r["independent_H6_and_H7_actual_field_exposure_counts_confirmed"],0)
        self.assertFalse(r["biology_EQR_NUTRI_HNS_accessed"])
        self.assertEqual(r["source_two_outlet_mass_pair_and_biological_causal_models"],"STOP")

    def test_adding_one_extra_department_wet_t_fails(self):
        d=copy.deepcopy(self.depts)
        d[0]["applied_wet_sludge_t"]="17315.29"
        with self.assertRaises(ValueError):m.audit(d,self.products)

    def test_adding_one_extra_department_parcel_fails(self):
        d=copy.deepcopy(self.depts)
        d[0]["actual_parcel_count"]="189"
        with self.assertRaises(ValueError):m.audit(d,self.products)

    def test_changing_product_dry_t_fails_independent_table(self):
        products=copy.deepcopy(self.products)
        products[0]["applied_dry_sludge_t"]="11905.55"
        with self.assertRaises(ValueError):m.audit(self.depts,products)

    def test_year_mismatch_and_zero_not_source(self):
        d=copy.deepcopy(self.depts)
        d[0]["source_campaign_year"]="2021"
        with self.assertRaises(ValueError):m.audit(d,self.products)
        d=copy.deepcopy(self.depts)
        d[0]["actual_applied_area_ha"]="0"
        with self.assertRaises(ValueError):m.audit(d,self.products)

    def test_privacy_field_in_disaggregated_data_rejected(self):
        d=copy.deepcopy(self.depts)
        d[0]["farm_name"]="do not publish"
        with self.assertRaises(ValueError):m.audit(d,self.products)

    def test_duplicate_department_rejected(self):
        d=copy.deepcopy(self.depts)
        d[1]["department_code"]=d[0]["department_code"]
        with self.assertRaises(ValueError):m.audit(d,self.products)

    def test_no_substitution_of_another_product(self):
        p=copy.deepcopy(self.products)
        p[0]["sludge_product"]="compost"
        with self.assertRaises(ValueError):m.audit(self.depts,p)

if __name__=="__main__":
    unittest.main()
