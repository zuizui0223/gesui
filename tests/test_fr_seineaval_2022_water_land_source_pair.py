import copy,importlib.util,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_seineaval_2022_water_land_source_pair_audit.py"
spec=importlib.util.spec_from_file_location("crossmedia",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def source_pair():
    months=[]
    for i in range(1,13):
        months.append({"month":f"{i:02d}","A3_NGL_kg":1000.0*i,
            "A3_PT_kg":100*i,
            "A4_NGL_kg":400*i,"A4_PT_kg":20*i,
            "A5_NGL_kg":20*i,"A5_PT_kg":5*i,
            "A4_treated_water_m3":10000*i})
    annual={}
    for point in ("A3","A4","A5"):
        for chem in ("NGL","PT"):
            annual[f"{point}_{chem}_kg_2022"]=sum(v[f"{point}_{chem}_kg"] for v in months)
    annual["A4_treated_water_m3_2022"]=sum(v["A4_treated_water_m3"] for v in months)
    water={
      "status":"OFFICIAL_2022_VALIDATED_MONTHLY_A4_NGL_PT_MASS_RECONSTRUCTED",
      "expected_SHA256":m.PIN_AESN_SHA,"STEU":"037800501000","year":2022,
      "actual_source":{
        "original_source_discovery_run":44536921058,
        "original_reporting_year":2022,
        "monthly_2022_source_indicator_count":84,"original_STEU_rows":528,
        "A3_influent_not_effluent":True,
        "A5_bypass_separate_not_summed_into_A4":True,
        "no_EQR_NUTRI_HNS_used":True,
        "biological_effect_causal_path_GATE":"STOP",
        "monthly":months,
        "annual_source_flux_indicators":annual,
      }}
    soil={
      "schema":"gesui_SIAAP_SeineAval_2022_actual_applied_lot_weighted_NTK_P2O5_SOURCE_ESTIMATE_v1",
      "actually_spread_campaign_year":2022,
      "actual_applied_wet_t":81195.15,
      "independent_WWTP_sources":1,
      "NTK_is_not_total_N_for_NGL_comparison":True,
      "coupled_water_land_mass_balance_and_mediation_GO":False,
      "source_only_no_EQR_NUTRI_HNS":True,
      "estimated_applied_Kjeldahl_N_NTK_t":704.35879,
      "estimated_applied_P2O5_t":3402.411135,
      "estimated_applied_elemental_P_equivalent_t":1484.9034630074,
    }
    return water,soil

class PairedSourceAuditTest(unittest.TestCase):
    def test_same_2022_source_at_two_outlets_not_closed_mass_balance(self):
        water,soil=source_pair()
        result=m.combine(water,soil)
        self.assertEqual(result["matched_calendar_year"],2022)
        self.assertEqual(result["origin_ministry_SANDRE"],"037800501000")
        self.assertFalse(result["true_closed_2022_plant_N_P_source_mass_balance"])
        self.assertTrue(result["nitrogen_chemical_species_not_identical_NTK_vs_NGL"])
        self.assertEqual(result["causal_N_P_transfer_and_biotic_effect_gate"],"STOP")
        self.assertTrue(result["all_outcome_biology_EQR_NUTRI_HNS_unopened"])
    def test_separate_aqueous_A5_bypass_from_A4_treated_load(self):
        water,soil=source_pair()
        result=m.combine(water,soil)
        self.assertAlmostEqual(result["aquatic_official_validated_monthly_calculated_indicator_A4"]["NGL_t_2022"],31.2)
        self.assertAlmostEqual(result["aquatic_additional_bypass_A5_not_merged_into_A4"]["NGL_t_2022"],1.56)
        self.assertTrue(result["aquatic_additional_bypass_A5_not_merged_into_A4"]["NGL_t_2022"]>0)
    def test_original_monthly_water_totals_cannot_be_fabricated(self):
        water,soil=source_pair()
        water["actual_source"]["monthly"][0]["A4_NGL_kg"]=999999.0
        with self.assertRaises(ValueError):
            m.combine(water,soil)
    def test_wrong_calendar_year_rejected(self):
        water,soil=source_pair()
        soil["actually_spread_campaign_year"]=2021
        with self.assertRaises(ValueError):m.combine(water,soil)
    def test_missing_original_aquatic_SHA_rejected(self):
        water,soil=source_pair()
        water["expected_SHA256"]="not-pinned"
        with self.assertRaises(ValueError):m.combine(water,soil)
    def test_missing_one_month_rejected(self):
        water,soil=source_pair()
        water["actual_source"]["monthly"].pop()
        with self.assertRaises(ValueError):m.combine(water,soil)
    def test_wrong_elemental_P2O5_conversion_rejected(self):
        water,soil=source_pair()
        soil["estimated_applied_elemental_P_equivalent_t"]=3402.411135
        with self.assertRaises(ValueError):m.combine(water,soil)
    def test_aquatic_yearly_flux_comparison_does_not_create_plant_biotic_effect(self):
        water,soil=source_pair()
        out=m.combine(water,soil)
        self.assertTrue(out["same_year_source_outlet_elemental_P_magnitude_comparison"]["source_mass_comparison_not_2022_material_allocation_fraction"])
        self.assertFalse(out["field_specific_exposure_and_no_spreading_controls_verified"])
    def test_pseudoreplication_or_different_plant_is_refused(self):
        water,soil=source_pair()
        water["STEU"]="0442094S0002"
        with self.assertRaises(ValueError):m.combine(water,soil)

if __name__=="__main__":unittest.main()
