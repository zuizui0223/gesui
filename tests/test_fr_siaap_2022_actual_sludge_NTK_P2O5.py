import csv,importlib.util,copy,unittest
from pathlib import Path
root=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022"
spec=importlib.util.spec_from_file_location("siaapnutrient",root/"fr_siaap_2022_actual_sludge_NTK_elementalP_SOURCE_screen.py")
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with (root/"FR_SIAAP_SEINE_AVAL_2022_REALIZED_SLUDGE_BY_PRODUCT.csv").open(newline="") as f:PRODUCTS=list(csv.DictReader(f))
with (root/"FR_SIAAP_SEINE_AVAL_2022_APPLIED_SLUDGE_NTK_P2O5_AGRONOMIC.csv").open(newline="") as f:COMPOSITION=list(csv.DictReader(f))

class ActualAgronomicNutrientSourceTest(unittest.TestCase):
    def test_lot_weighted_applied_NTK_and_P2O5_calculation(self):
        s=m.estimate(PRODUCTS,COMPOSITION)
        self.assertAlmostEqual(s["estimated_applied_Kjeldahl_N_NTK_t"],704.35879,places=5)
        self.assertAlmostEqual(s["estimated_applied_P2O5_t"],3402.411135,places=5)
        self.assertAlmostEqual(s["estimated_applied_elemental_P_equivalent_t"],1484.9034630074,places=5)
        self.assertEqual(s["actually_spread_campaign_year"],2022)
        self.assertEqual(s["actual_applied_wet_t"],81195.15)
    def test_products_applied_vs_produced_must_not_be_confused(self):
        products=copy.deepcopy(PRODUCTS)
        products[0]["applied_wet_sludge_t"]="25141"
        with self.assertRaises(ValueError):m.estimate(products,COMPOSITION)
    def test_applied_vs_produced_composition_distinguished(self):
        comp=copy.deepcopy(COMPOSITION)
        comp[0]["applied_sludge_P2O5_kg_per_wet_t_lot_weighted_mean"]="48.3"
        with self.assertRaises(ValueError):m.estimate(PRODUCTS,comp)
    def test_P2O5_not_elemental_P_kg(self):
        s=m.estimate(PRODUCTS,COMPOSITION)
        self.assertLess(s["estimated_applied_elemental_P_equivalent_t"],s["estimated_applied_P2O5_t"])
        self.assertGreater(s["elemental_P_over_P2O5_factor"],0.436)
        self.assertLess(s["elemental_P_over_P2O5_factor"],0.437)
    def test_NTK_not_identical_to_aqueous_NGL(self):
        s=m.estimate(PRODUCTS,COMPOSITION)
        self.assertTrue(s["NTK_is_not_total_N_for_NGL_comparison"])
        self.assertFalse(s["coupled_water_land_mass_balance_and_mediation_GO"])
        self.assertTrue(s["agronomic_chemical_N_P_is_not_field_host_plant_uptake"])
    def test_same_year_even_if_some_delivered_2021(self):
        s=m.estimate(PRODUCTS,COMPOSITION)
        self.assertTrue(s["2021_production_or_delivery_lots_can_be_spread_in_2022"])
    def test_missing_lot_composition_fails(self):
        with self.assertRaises(ValueError):m.estimate(PRODUCTS,COMPOSITION[:1])
    def test_comp_bad_year_fails(self):
        comp=copy.deepcopy(COMPOSITION)
        comp[0]["source_campaign_year"]="2021"
        with self.assertRaises(ValueError):m.estimate(PRODUCTS,comp)
    def test_anonymous_only_no_private_parcel_fields(self):
        comp=copy.deepcopy(COMPOSITION)
        comp[0]["parcel_id"]="private"
        with self.assertRaises(ValueError):m.estimate(PRODUCTS,comp)

if __name__=="__main__":unittest.main()
