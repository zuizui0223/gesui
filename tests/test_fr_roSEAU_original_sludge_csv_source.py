import importlib.util,unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_roSEAU_original_sludge_csv_source_audit.py"
spec=importlib.util.spec_from_file_location("fr_sludge_csv",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class OfficialROSEAUCSVTests(unittest.TestCase):
    def test_pipe_export_detects_application_commune(self):
        csvtext="annee|CodeSTEU|communes_epandage|production_boues_tMS|traitement|statut\n2018|X01|75056|18.2|boues|actual\n2018|X02|31069|0|boues|actual\n"
        parsed=m.sniff(csvtext.encode())
        self.assertEqual(parsed["delimiter"],"|")
        self.assertEqual(parsed["record_rows"],2)
        self.assertEqual(len(parsed["columns"]),6)
        self.assertIn("communes_epandage",[x["column"] for x in parsed["source_relevant_fields"]])
    def test_rejects_html_not_csv(self):
        with self.assertRaises(ValueError):
            m.sniff(b"<html><body>no data</body></html>")
    def test_rejects_malformed_rows(self):
        with self.assertRaises(ValueError):
            m.sniff(b"a|b|c|d|e\n1|2|3|4\n")
    def test_year_scope_frozen(self):
        self.assertEqual(m.PERMITTED_YEARS,{"2017","2018","2019","2020","2021","2022"})
        self.assertTrue(m.EXPORT.endswith("/controllers/GetDataExportBoueHandler.php"))
if __name__=="__main__":unittest.main()
