import importlib.util,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_aesn_2022_autosurveillance_csv_header_source_probe.py"
spec=importlib.util.spec_from_file_location("aesn_header",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class AESNOfficialCSVHeader(unittest.TestCase):
    def source(self,url="https://static.data.gouv.fr/resources/demo/autosurveillance-donnees-2022-01.csv"):
        return {"resources":[{"title":m.MONTH_NAME,"format":"csv","url":url,"filesize":110379713,"id":"fixture"}]}
    def test_exact_original_month_url_only(self):
        url,meta=m.original_resource_url(self.source())
        self.assertTrue(url.startswith("https://static.data.gouv.fr/"))
        self.assertEqual(meta["official_title"],m.MONTH_NAME)
    def test_external_host_rejected(self):
        with self.assertRaises(ValueError):
            m.original_resource_url(self.source("https://example.org/autosurveillance-donnees-2022-01.csv"))
    def test_wrong_year_does_not_count(self):
        with self.assertRaises(ValueError):
            m.original_resource_url({"resources":[{"title":"autosurveillance-donnees-2021-01.csv","format":"csv","url":"https://static.data.gouv.fr/demo.csv"}]})
    def test_header_fields_are_candidates_not_source_flux(self):
        b="code_st e u;date_prelevement;code_point;code_parametre;valeur;unite\n".replace(" ","").encode("utf8")
        x=m.parse_first_header(b)
        self.assertEqual(x["column_count"],6)
        self.assertEqual(x["detected_delimiter"],";")
        self.assertIn("code_parametre",x["matching_candidate_column_names"]["potential_NGL_PT_parameter"])
        self.assertTrue(x["2022_SeineAval_specific_records_not_read"])
        self.assertTrue(x["annual_measured_NGL_PT_not_calculated"])
    def test_cp1252_header_accepted_but_no_value_read(self):
        x=m.parse_first_header("station;date_prélèvement;unité;valeur\n".encode("cp1252"))
        self.assertEqual(x["header_encoding"],"cp1252")
        self.assertEqual(x["column_count"],4)
    def test_duplicate_header_names_fail(self):
        with self.assertRaises(ValueError):
            m.parse_first_header(b"station;station;date\n")
if __name__=="__main__":unittest.main()
