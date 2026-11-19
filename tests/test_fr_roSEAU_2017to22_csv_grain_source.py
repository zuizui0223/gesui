import importlib.util,unittest
from pathlib import Path

P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_roSEAU_2017to22_csv_grain_source_audit.py"
spec=importlib.util.spec_from_file_location("grain",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class ROSEAUGrainTests(unittest.TestCase):
    def makecsv(self,rows):
        cols=list(m.FIELD.values())+[f"extra_col_{i}" for i in range(32-len(m.FIELD))]
        self.assertGreaterEqual(len(cols),30)
        return ("|".join(cols)+"\n"+"\n".join("|".join(row.get(col,"") for col in cols) for row in rows)+"\n").encode()
    def test_duplicates_are_not_independent_spreading_events(self):
        rows=[]
        for parameter in ["N","P","Zn"]:
            rows.append({
                m.FIELD["station"]:"ST01",m.FIELD["commune"]:"12345",
                m.FIELD["spreading"]:"12,5",m.FIELD["produced"]:"20",
                m.FIELD["sample_date"]:"2018-06-01",m.FIELD["sample_year"]:"2018",
                m.FIELD["parameter"]:parameter,m.FIELD["material"]:"MT1",
                m.FIELD["dept"]:"12",m.FIELD["region"]:"Occitanie",
                m.FIELD["bassin"]:"Garonne"
            })
        info=m.classify(self.makecsv(rows))
        self.assertEqual(info["raw_records"],3)
        self.assertEqual(info["official_unique_stations_nonempty"],1)
        self.assertEqual(info["official_distinct_station_commune_literal_pairs"],1)
        self.assertEqual(info["distinct_explicit_INSEE_codes_only_not_validated_against_official_INSEE_geometry"],1)
        self.assertEqual(info["raw_quantity_positive_zero_invalid_or_blank"]["positive"],3)
        self.assertFalse(info["looks_like_one_row_per_independent_spreading_event"])
    def test_labeled_commune_not_exact_INSEE(self):
        self.assertEqual(m.type_commune("Lyon (69123)"),"label_plus_INSEE_code")
        self.assertEqual(m.type_commune("Paris"),"text_name_or_unknown")
        self.assertEqual(m.type_commune("75056"),"exact_INSEE_code")
        self.assertEqual(m.type_commune("2A004"),"exact_INSEE_code")
        self.assertEqual(m.type_commune(""),"empty")
        self.assertEqual(m.type_commune("City1, City2"),"multi_or_delimited_text_UNRESOLVED")
    def test_numeric_not_ambiguous(self):
        from decimal import Decimal
        self.assertEqual(m.numeric("1 234,5"),Decimal("1234.5"))
        self.assertEqual(m.numeric("0"),0)
        self.assertIsNone(m.numeric("-1"))
        self.assertIsNone(m.numeric("1,234.5"))
        self.assertIsNone(m.numeric(""))
    def test_years_are_fixed(self):
        self.assertEqual(m.YEARS,{"2017","2018","2019","2020","2021","2022"})
if __name__=="__main__":unittest.main()
