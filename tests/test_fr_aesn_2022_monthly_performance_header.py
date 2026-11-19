import importlib.util,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_aesn_2022_monthly_performance_indicator_HEADER_gate.py"
spec=importlib.util.spec_from_file_location("perf",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def catalog(url=None,size=53300000):
    if url is None:
        url="https://static.data.gouv.fr/resources/assainissement-stations-depuration-indicateurs-de-performance/20240531-145021/"+m.FILENAME
    return {"resources":[{"title":m.FILENAME,"filename":m.FILENAME,"url":url,"filesize":size,"id":"official-ID","format":"csv"}]}

class SourceMonthlyPerformanceTest(unittest.TestCase):
    def test_exact_official_resource_2022_confirmed(self):
        x=m.approved_resource(catalog())
        self.assertTrue(x["url"].endswith(m.FILENAME))
        self.assertEqual(x["expected_bytes"],53300000)
    def test_catalog_human_readable_title_with_original_official_filename_in_URL(self):
        d=catalog()
        d["resources"][0]["title"]="Résultats annuels des stations"
        d["resources"][0]["filename"]=None
        a=m.approved_resource(d)
        self.assertTrue(a["url"].endswith(m.FILENAME))

    def test_nonofficial_domain_rejected(self):
        with self.assertRaises(ValueError):
            m.approved_resource(catalog("https://example.com/"+m.FILENAME))
    def test_wrong_year_rejected(self):
        with self.assertRaises(ValueError):
            m.approved_resource({"resources":[{"title":"stepcol-indicateurs-mensuels-2021.csv"}]})
    def test_implausible_resource_size_rejected(self):
        with self.assertRaises(ValueError):
            m.approved_resource(catalog(size=10))
    def test_header_only_source_metadata_no_nutrient_estimation(self):
        h=m.header("ANNEE;MOIS;Station_Code SANDRE;Point;Parametre;Flux;Unite;Concentration;Rendement\n".encode())
        self.assertEqual(h["field_count"],9)
        self.assertIn("Flux",h["possible_semantic_columns"]["unit_mass_flow"])
        self.assertTrue(h["N_P_annual_mass_not_measured_from_header"])
        self.assertTrue(h["Seine_Aval_exact_station_match_NOT_yet_checked"])
    def test_corrupted_or_duplicate_column_schema_rejected(self):
        with self.assertRaises(ValueError):
            m.header(b"a;a;b;c;d;e\n")
if __name__=="__main__":unittest.main()
