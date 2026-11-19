import importlib.util,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_eea2022_seine_aval_measured_NP_exact_station_gate.py"
spec=importlib.util.spec_from_file_location("seine_eea",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
class EEAExactPlantTest(unittest.TestCase):
    def plant(self,code="FR037800501000",rep="FR2022                  ",n=105.0,p=21.0):
        return {"uwwCode":code,"repCode":rep,"uwwName":"PARIS SEINE AVAL",
                "uwwNDischargeMeasured":n,"uwwPDischargeMeasured":p,
                "uwwNDischargeEstimated":1000,"uwwPDischargeCalculated":2000}
    def outlet(self,code="FR037800501000",rep="FR2022",wb="FW"):
        return {"uwwCode":code,"repCode":rep,"dcpWaterBodyType":wb}
    def test_verified_exact_join_both_measured_and_fw(self):
        a=m.match_exact([self.plant()],[self.outlet()])
        self.assertEqual(a["status"],"EXACT_EEA_2022_MEASURED_FRESHWATER_SOURCE_MATCH")
        self.assertEqual(a["reported_measured_annual_N_P_t"],{"N_t_per_2022":105.0,"P_t_per_2022":21.0})
        self.assertFalse(a["paired_measured_water_N_P_and_2022_actual_field_N_P_mass_GO"])
    def test_correct_identifiers_require_both_fr_prefix_and_exact_official_name(self):
        a=m.match_exact([self.plant(code="037800501000")],[self.outlet()])
        self.assertIsNone(a["reported_measured_annual_N_P_t"])
        b=self.plant()
        b["uwwName"]="Another Seine Aval WWTP"
        c=m.match_exact([b],[self.outlet()])
        self.assertIsNone(c["reported_measured_annual_N_P_t"])

    def test_padding_in_original_reporter_codes_normalized(self):
        a=m.match_exact([self.plant()],[self.outlet(rep="FR2022                ")])
        self.assertTrue(a["confirmed_one_direct_Fw_outfall"])

    def test_no_verified_exact_id_by_name_alone(self):
        a=m.match_exact([self.plant("FRENCH_AVAL_CODE")],[])
        self.assertEqual(a["status"],"EXACT_STEU_NOT_UNIQUELY_CONFIRMED")
        self.assertIsNone(a["reported_measured_annual_N_P_t"])
    def test_wrong_country_rejected(self):
        a=m.match_exact([self.plant(rep="DE_TEST_2022")],[self.outlet(rep="DE_TEST_2022")])
        self.assertEqual(a["status"],"EXACT_STEU_NOT_UNIQUELY_CONFIRMED")
    def test_estimated_load_does_not_become_measured(self):
        a=m.match_exact([self.plant(n=None)],[self.outlet()])
        self.assertEqual(a["status"],"EXACT_STEU_PRESENT_INCOMPLETE_SOURCE")
        self.assertIsNone(a["reported_measured_annual_N_P_t"])
    def test_multiple_outlets_undetermined(self):
        a=m.match_exact([self.plant()],[self.outlet(),self.outlet()])
        self.assertEqual(a["same_key_outfall_count"],2)
        self.assertIsNone(a["reported_measured_annual_N_P_t"])
    def test_estuary_not_direct_freshwater(self):
        a=m.match_exact([self.plant()],[self.outlet(wb="ES")])
        self.assertFalse(a["confirmed_one_direct_Fw_outfall"])
        self.assertIsNone(a["reported_measured_annual_N_P_t"])
    def test_missing_p_prevents_pair(self):
        a=m.match_exact([self.plant(p=None)],[self.outlet()])
        self.assertFalse(a["reported_P_measurement_available"])
        self.assertIsNone(a["reported_measured_annual_N_P_t"])
    def test_duplicate_exact_station_rejected(self):
        a=m.match_exact([self.plant(),self.plant()],[self.outlet()])
        self.assertEqual(a["status"],"EXACT_STEU_NOT_UNIQUELY_CONFIRMED")
if __name__=="__main__":unittest.main()
