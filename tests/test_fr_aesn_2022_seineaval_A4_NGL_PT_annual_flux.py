import csv,copy,importlib.util,tempfile,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_siaap_2022/fr_aesn_SeineAval_2022_A4_NGL_PT_annual_flux_SOURCE.py"
spec=importlib.util.spec_from_file_location("seineflux",P)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def make_rows():
    out=[]
    for i in range(1,13):
        for point in ("A3","A4","A5"):
            for p in ("1551","1350"):
                out.append({"Annee":"2022","Mois":str(i),
                    "Station_Code SANDRE":"037800501000","Station_Nom":"SEINE AVAL",
                    "Point_Prelevement_Code":point,"Parametre_Code SANDRE":p,
                    "Resultat_Unite_Code":"67","Resultat_Unite_Symbole":"kg",
                    "Etat_calcul":"VALIDE","Resultat_Valeur":str((100 if p=="1551" else 10)*i)})
        out.append({"Annee":"2022","Mois":str(i),
                "Station_Code SANDRE":"037800501000","Station_Nom":"SEINE AVAL",
                "Point_Prelevement_Code":"A4","Parametre_Code SANDRE":"1552",
                "Resultat_Unite_Code":"115","Resultat_Unite_Symbole":"m3",
                "Etat_calcul":"VALIDE","Resultat_Valeur":str(5000*i)})
        for j in range(37):
            out.append({"Annee":"2022","Mois":str(i),
                "Station_Code SANDRE":"037800501000","Station_Nom":"SEINE AVAL",
                "Point_Prelevement_Code":"A2","Parametre_Code SANDRE":str(2000+j),
                "Resultat_Unite_Code":"67","Resultat_Unite_Symbole":"kg",
                "Etat_calcul":"VALIDE","Resultat_Valeur":"999"})
    assert len(out)==528
    return out

class AnnualFluxFromOfficialMonthlyIndicatorsTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.p=Path(self.tmp.name)/"2022.csv"
    def tearDown(self):
        self.tmp.cleanup()
    def save(self,rows):
        with self.p.open("w",newline="",encoding="utf8") as f:
            w=csv.DictWriter(f,fieldnames=m.EXPECTED_COLS,delimiter=";")
            w.writeheader()
            for r in rows:w.writerow({k:r.get(k,"") for k in m.EXPECTED_COLS})
        return self.p
    def test_12_months_kg_and_A4_water_export_close(self):
        d=m.scan(self.save(make_rows()))
        self.assertEqual(d["original_STEU_rows"],528)
        self.assertEqual(d["monthly_2022_source_indicator_count"],84)
        self.assertAlmostEqual(d["annual_source_flux_indicators"]["A4_NGL_kg_2022"],7800)
        self.assertAlmostEqual(d["annual_source_flux_indicators"]["A4_PT_kg_2022"],780)
        self.assertAlmostEqual(d["annual_source_flux_indicators"]["A4_treated_water_m3_2022"],390000)
        self.assertTrue(d["A5_bypass_separate_not_summed_into_A4"])
        self.assertEqual(d["biological_effect_causal_path_GATE"],"STOP")
    def test_missing_one_outlet_month_rejected(self):
        x=make_rows()
        x=[a for a in x if not (a["Mois"]=="2" and a["Point_Prelevement_Code"]=="A4" and a["Parametre_Code SANDRE"]=="1551")]
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_unvalidated_indicator_rejected(self):
        x=make_rows()
        x[0]["Etat_calcul"]="ESTIME"
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_wrong_unit_kg_per_day_not_kg_rejected(self):
        x=make_rows()
        x[0]["Resultat_Unite_Symbole"]="kg/j"
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_water_volume_not_NP_mass(self):
        x=make_rows()
        z=next(y for y in x if y["Point_Prelevement_Code"]=="A4" and y["Parametre_Code SANDRE"]=="1552")
        z["Resultat_Unite_Symbole"]="kg"
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_pseudoreplication_duplicate_source_record_rejected(self):
        x=make_rows()
        x[1]=copy.deepcopy(x[0])
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_source_station_name_mismatch_rejected(self):
        x=make_rows()
        x[0]["Station_Nom"]="ANOTHER WWTP"
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_year_mismatch_rejected(self):
        x=make_rows()
        x[0]["Annee"]="2021"
        with self.assertRaises(ValueError):m.scan(self.save(x))
    def test_numeric_invalid_rejected(self):
        for v in ("NaN","-1","Inf","imputed"):
            with self.assertRaises(ValueError):m.decimal_value(v)
    def test_decimal_comma_accepted(self):
        self.assertEqual(m.decimal_value("123,45"),m.Decimal("123.45"))
    def test_reporter_original_SHA_freeze_is_defined(self):
        self.assertEqual(len(m.PIN_SHA),64)
        self.assertEqual(m.PIN_BYTES,55938440)
    def test_biology_not_inspected_in_any_source(self):
        d=m.scan(self.save(make_rows()))
        self.assertTrue(d["no_EQR_NUTRI_HNS_used"])
        self.assertTrue(d["not_a_closed_WWTP_N_P_mass_balance"])
if __name__=="__main__":unittest.main()
