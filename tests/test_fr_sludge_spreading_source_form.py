import importlib.util,unittest
from pathlib import Path
from bs4 import BeautifulSoup

P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_sillage_spreading_form_source_probe.py"
spec=importlib.util.spec_from_file_location("fr_sludge",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class FormSourceTests(unittest.TestCase):
    def test_get_post_form_year_options(self):
        html='<form method="post" action="/export/boues"><select name="annee"><option value="2017">2017</option><option value="2018">2018</option><option value="2022">2022</option></select><select name="format"><option value="csv">CSV</option></select><input name="csrf_token" value="do-not-export-me"></form>'
        d=m.form_info(BeautifulSoup(html,"html.parser"),"https://host.gouv.fr/pages/boue")
        self.assertEqual(d["visible_2017_to_2022_options"],["2017","2018","2022"])
        self.assertEqual(d["forms"][0]["method"],"POST")
        self.assertEqual(d["forms"][0]["action"],"https://host.gouv.fr/export/boues")
        self.assertTrue(d["forms"][0]["controls"][-1]["default_value_redacted"])
        self.assertNotIn("do-not-export-me",str(d))
    def test_off_origin_is_not_same_origin(self):
        d=m.form_info(BeautifulSoup('<form action="https://unknown.test/post"></form>',"html.parser"),
                      "https://host.gouv.fr/source")
        self.assertFalse(d["forms"][0]["same_origin"])
    def test_extract_public_js_clues_no_secret(self):
        s='const password="supersecret"; function download_boues(){return "export csv"}'
        x=m.safe_js_excerpt(s)
        self.assertTrue(x)
        self.assertNotIn("supersecret"," ".join(x))
    def test_original_wanted_year_fixed(self):
        self.assertEqual(m.WANTED_YEARS,{"2017","2018","2019","2020","2021","2022"})
if __name__=="__main__":unittest.main()
