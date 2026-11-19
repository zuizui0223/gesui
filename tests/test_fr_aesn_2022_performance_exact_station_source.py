import csv,importlib.util,tempfile,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_aesn_2022_performance_exact_sandre_source_coverage.py"
spec=importlib.util.spec_from_file_location("perf2022",p)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def base():
    return {"Annee":"2022","Mois":"1","Station_Code SANDRE":"037800501000",
            "Station_Nom":"PARIS SEINE AVAL",
            "Point_Prelevement_Code":"A4",
            "Point_Prelevement_Libelle":"Sortie station",
            "Parametre_Code SANDRE":"1551",
            "Parametre_Libelle":"NGL","Resultat_Valeur":"9999999",
            "Resultat_Unite_Code":"168","Resultat_Unite_Symbole":"mg(N)/L",
            "Etat_calcul":"C"}

class MonthlySourceCoverageTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.p=Path(self.tmp.name)/"source.csv"
    def tearDown(self):
        self.tmp.cleanup()
    def make(self,rows,cols=m.EXPECTED_COLS):
        with self.p.open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=cols,delimiter=";")
            w.writeheader()
            for row in rows:
                w.writerow({x:row.get(x,"") for x in cols})
        return self.p
    def test_exact_sandre_station_and_NGL_PT_point_months(self):
        ngl=base()
        pt={**base(),"Parametre_Code SANDRE":"1350","Parametre_Libelle":"PT"}
        jan={**base(),"Mois":"2"}
        self.make([ngl,pt,jan])
        d=m.audit_rows(self.p)
        self.assertEqual(d["exact_station_record_rows"],3)
        self.assertEqual(d["NGL_PT_1551_1350_month_coverage_by_point"],
                         [{"point":"A4","SANDRE_parameter":"1350","distinct_month_count":1,"months":["01"]},
                          {"point":"A4","SANDRE_parameter":"1551","distinct_month_count":2,"months":["01","02"]}])
        self.assertTrue(d["original_Resultat_Valeur_NEVER_READ"])
        self.assertFalse(d["paired_water_land_nutrient_mediation_GO"])
    def test_entrance_A3_not_mixed_with_A4_outlet(self):
        d={**base(),"Point_Prelevement_Code":"A3"}
        self.make([d,base()])
        s=m.audit_rows(self.p)
        self.assertEqual(len(s["NGL_PT_1551_1350_month_coverage_by_point"]),2)
    def test_name_only_cannot_become_exact_station(self):
        d={**base(),"Station_Code SANDRE":"999000", "Station_Nom":"PARIS SEINE AVAL"}
        self.make([d])
        s=m.audit_rows(self.p)
        self.assertEqual(s["exact_station_record_rows"],0)
        self.assertEqual(len(s["nonmatching_name_candidates"]),1)
    def test_numeric_values_do_not_affect_source_coverage(self):
        self.make([base()])
        a=m.audit_rows(self.p)
        d={**base(),"Resultat_Valeur":"cannot-be-converted"}
        self.make([d])
        b=m.audit_rows(self.p)
        self.assertEqual(a,b)
    def test_wrong_year_is_not_2022_application(self):
        d={**base(),"Annee":"2018"}
        self.make([d])
        a=m.audit_rows(self.p)
        self.assertEqual(a["exact_station_record_year_mismatches"],1)
        self.assertEqual(a["exact_station_record_rows"],0)
    def test_private_or_unexpected_fields_rejected(self):
        self.make([base()],cols=m.EXPECTED_COLS+("farm_name",))
        with self.assertRaises(ValueError):m.audit_rows(self.p)
    def test_approved_distribution_only(self):
        url="https://static.data.gouv.fr/resources/a/"+m.NAME
        d={"resources":[{"title":"2022 production","url":url,"filesize":53500000}]}
        result,size=m.original_resource(d)
        self.assertEqual(size,53500000)
        self.assertEqual(result,url)
    def test_untrusted_resource_not_approved(self):
        d={"resources":[{"title":"2022 production","url":"https://example.org/"+m.NAME,"filesize":53500000}]}
        with self.assertRaises(ValueError):m.original_resource(d)
    def test_source_not_analyzed_from_header_only(self):
        self.assertIn("Resultat_Valeur",m.EXPECTED_COLS)
        self.assertEqual(m.PARAMS,{"1551":"NGL_azote_global","1350":"PT_phosphore_total"})

if __name__=="__main__":unittest.main()
