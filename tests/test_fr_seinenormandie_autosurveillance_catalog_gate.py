import importlib.util,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_seinenormandie_2022_autosurveillance_catalog_gate.py"
spec=importlib.util.spec_from_file_location("aesn",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
class CatalogTest(unittest.TestCase):
    def row(self,y,month):
        return {"title":f"autosurveillance-donnees-{y}-{month:02d}.csv","id":f"resource-{y}-{month:02d}",
                "format":"csv","filesize":10000000,"url":"https://example.invalid/data.csv"}
    def test_fixed_2022_2018_month_inventory(self):
        rows=[self.row(y,i) for y in (2018,2022) for i in range(1,13)]
        result=m.parse_catalog({"title":"official source","resources":rows})
        self.assertEqual(len(result["resource_names_2022"]),12)
        self.assertEqual(len(result["resource_names_2018"]),12)
        self.assertTrue(result["has_12_distinct_months_for_2022"])
        self.assertTrue(result["has_12_distinct_months_for_2018"])
        self.assertTrue(result["no_original_data_rows_opened"])
        self.assertTrue(result["reported_treated_effluent_NGL_PT_not_yet_identified"])
        self.assertFalse(result["paired_water_land_2022_mass_GO"])
    def test_2022_coverage_incomplete_never_fabricated(self):
        result=m.parse_catalog({"resources":[self.row(2022,7)]})
        self.assertFalse(result["has_12_distinct_months_for_2022"])
        self.assertEqual(result["source_years_with_monthly_resources"],[2022])
    def test_duplicate_month_is_not_independent_resource(self):
        with self.assertRaises(ValueError):
            m.parse_catalog({"resources":[self.row(2022,7),self.row(2022,7)]})
    def test_unknown_schema_rejected(self):
        with self.assertRaises(ValueError):
            m.parse_catalog({"title":"not a public dataset"})
    def test_2019_index_never_called_2022(self):
        r=m.parse_catalog({"resources":[self.row(2019,4)]})
        self.assertEqual(r["resource_names_2022"],[])
        self.assertFalse(r["has_12_distinct_months_for_2022"])
if __name__=="__main__":unittest.main()
