import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_aesn_2022_january_seineaval_site_point_source_audit.py"
spec=importlib.util.spec_from_file_location("aesn_2022_site",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
COLUMNS=[
 "Année","Mois","Finalité analyse","Type_Code","Type_Libellé",
 "Station_No Interne","Station_Code SANDRE","Station_Nom",
 "Departement_Numero","Commune_Code INSEE","Commune_Libellé",
 "XL93","YL93","Station_Capacite nominale",
 "Prelevement_Commentaire","Point_Prelevement","Prelevement_Date",
 "Parametre_Code SANDRE","Parametre_Libelle","Parametre_Type",
 "Support_Code","Support_Libelle","Fraction_Code","Fraction_Nom",
 "Resultat_Valeur","Resultat_Unite_Code","Resultat_Unite_Symbole",
 "Resultat_Unite_Libelle","Limite quantification","Qualification"]

def row(code="037800501000",name="PARIS SEINE AVAL",point="A4",param="1551",value="999"):
    r={k:"" for k in COLUMNS}
    r.update({"Année":"2022","Mois":"1","Station_Code SANDRE":code,
       "Station_Nom":name,"Point_Prelevement":point,
       "Prelevement_Date":"2022-01-11",
       "Parametre_Code SANDRE":param,
       "Parametre_Libelle":"Azote global" if param=="1551" else "Phosphore total",
       "Parametre_Type":"C","Resultat_Unite_Symbole":"mg/L",
       "Resultat_Valeur":value,"Qualification":"1"})
    return r

class AESNSeineAvalSourceTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.p=Path(self.tmp.name)/m.NAME
    def tearDown(self):
        self.tmp.cleanup()
    def write(self,rows,fields=COLUMNS):
        with self.p.open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter=";")
            w.writeheader()
            for r in rows:
                w.writerow({key:r.get(key,"") for key in fields})
    def test_original_official_URL_restricts_single_month_and_size(self):
        d={"resources":[{"title":m.NAME,"format":"csv","filesize":110379713,
              "url":"https://static.data.gouv.fr/resources/a/"+m.NAME}]}
        url,size=m.approved_url(d)
        self.assertEqual(size,110379713)
        self.assertIn(m.NAME,url)
    def test_dangerous_or_nonofficial_distribution_rejected(self):
        d={"resources":[{"title":m.NAME,"format":"csv","filesize":110379713,
              "url":"https://example.org/"+m.NAME}]}
        with self.assertRaises(ValueError):m.approved_url(d)
    def test_actual_A4_parameter_coverage_without_values(self):
        self.write([row(),row(param="1350",value="88888"),row(code="another",name="other")])
        z=m.count_site(self.p)
        self.assertEqual(z["exact_SeineAval_station_code_rows"],2)
        self.assertEqual(z["exact_code_station_point_row_counts"]["A4"],2)
        self.assertEqual(z["exact_code_point_parameter_NGL_PT_counts"]["A4"],{"1551":1,"1350":1})
        self.assertTrue(z["numeric_Resultat_Valeur_never_read"])
        self.assertFalse(z["paired_2022_water_land_mass_and_causal_response_GO"])
    def test_id_is_not_inferred_from_station_name(self):
        self.write([row(code="somethingelse"),row(code="FR037800501000")])
        z=m.count_site(self.p)
        self.assertEqual(z["exact_SeineAval_station_code_rows"],1)
        self.assertEqual(len(z["other_name_candidates_not_original_ID_confirmed"]),1)
    def test_non_2022_month_rows_are_not_2022_sources(self):
        r=row();r["Année"]="2021"
        self.write([r])
        z=m.count_site(self.p)
        self.assertEqual(z["original_source_year_month_mismatch_rows"],1)
        self.assertEqual(z["exact_SeineAval_station_code_rows"],0)
    def test_distinguish_A3_influent_from_A4_effluent(self):
        self.write([row(point="A3"),row(point="A4")])
        z=m.count_site(self.p)
        self.assertEqual(z["exact_code_point_parameter_NGL_PT_counts"]["A4"],{"1551":1})
        self.assertEqual(z["exact_code_point_parameter_NGL_PT_counts"]["A3"],{"1551":1})
    def test_unexpected_columns_fail_closed(self):
        self.write([row()],fields=COLUMNS+["farm_name"])
        with self.assertRaises(ValueError):m.count_site(self.p)
    def test_csv_schema_does_not_turn_means_into_annual_loads(self):
        self.write([row()])
        z=m.count_site(self.p)
        self.assertTrue(z["NGL_PT_quality_and_units_not_yet_sufficient_for_annual_load"])
        self.assertFalse(z["2022_annual_reported_water_N_P_load_confirmed"])
if __name__=="__main__":unittest.main()
