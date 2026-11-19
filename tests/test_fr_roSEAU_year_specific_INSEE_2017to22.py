import importlib.util,unittest,io,csv,zipfile,sys
from pathlib import Path
P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_roSEAU_annual_historic_INSEE_source_gate.py"
sys.path.insert(0,str(P.parent))
spec=importlib.util.spec_from_file_location("historicCog",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class HistoricalCOGTest(unittest.TestCase):
    def test_old_2018_txt_code_join(self):
        rows="CDC;CHEFLIEU;REG;DEP;COM;AR;CT;TNCC;ARTMAJ;NCC;ARTMIN;NCCENR\n0;0;53;29;049;3;02;0;;ELLIANT;;Elliant\n"
        out=io.BytesIO()
        with zipfile.ZipFile(out,"w") as z:z.writestr("comsimp2018.txt",rows.encode("iso-8859-1"))
        data,meta=m.official_cog("2018",out.getvalue(),min_count=1)
        self.assertEqual(data["29049"],{"elliant"})
        self.assertEqual(meta["official_unique_communes"],1)
    def test_current_csv_use_exact_COM_and_TYPECOM(self):
        raw="TYPECOM,COM,NCC,NCCENR,LIBELLE\nCOM,33033,BAURECH,Baurech,Baurech\nCOMD,33034,SOMETHING,Something,Something\n"
        data,meta=m.official_cog("2022",raw.encode(),min_count=1)
        self.assertIn("33033",data)
        self.assertNotIn("33034",data)
    def test_official_2019_csv_lowercase_fieldnames(self):
        txt='typecom,com,ncc,nccenr,libelle\\nCOM,01001,ABERGEMENT CLEMENCIAT,Abergement-Clemenciat,L Abergement-Clemenciat\\n'
        b=io.BytesIO()
        with zipfile.ZipFile(b,"w") as z:z.writestr("communes-01012019.csv",txt)
        result,meta=m.official_cog("2019",b.getvalue(),min_count=1)
        self.assertIn("01001",result)
        self.assertIn("abergementclemenciat",result["01001"])
    def test_frozen_source_hashes_for_all_years(self):
        self.assertEqual(set(m.SOURCE_SHA),{"2017","2018","2019","2020","2021","2022"})
        self.assertEqual(set(m.COG),set(m.SOURCE_SHA))
    def test_parse_commune_list_from_official_2018(self):
        c,e=m.parse_cell("29049|Elliant,29062|Gouézec")
        self.assertEqual(len(c),2)
        self.assertEqual(e,[])
if __name__=="__main__":unittest.main()
