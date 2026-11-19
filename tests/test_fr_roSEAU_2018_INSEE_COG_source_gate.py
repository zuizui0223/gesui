import importlib.util,unittest
from pathlib import Path
from decimal import Decimal

P=Path(__file__).resolve().parents[1]/"analysis/20261119_fr_sludge_routeC/fr_roSEAU_2018_original_INSEE_COG_source_gate.py"
spec=importlib.util.spec_from_file_location("cog_roSEAU",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class CommuneParserTests(unittest.TestCase):
    def test_original_2018_public_four_commune_example(self):
        sample="03208|Le Pin,03173|Molinet,71370|Rigny-sur-Arroux,71176|Digoin"
        vals,errors=m.parse_cell(sample)
        self.assertEqual(len(vals),4)
        self.assertEqual({x[0] for x in vals},{"03208","03173","71370","71176"})
        self.assertFalse(errors)
    def test_one_commune_field(self):
        vals,err=m.parse_cell("29049|Elliant")
        self.assertEqual(vals,[("29049","Elliant")])
        self.assertEqual(err,[])
    def test_unrecognized_name_is_not_made_into_code(self):
        self.assertEqual(m.parse_cell("Le Pin")[0],[])
        self.assertEqual(len(m.parse_cell("Le Pin")[1]),1)
    def test_distinguishes_empty_and_zero_sludge_station_total(self):
        self.assertIsNone(m.quantity(""))
        self.assertEqual(m.quantity("0"),Decimal(0))
        self.assertEqual(m.quantity("1 234,5"),Decimal("1234.5"))
        self.assertIsNone(m.quantity("-1"))
    def test_commune_names_normalize_without_guessing_code(self):
        self.assertEqual(m.commune_normalized("Scaër"),"scaer")
        self.assertEqual(m.commune_normalized("MARCILLAC"),"marcillac")
        self.assertNotEqual(m.commune_normalized("Scaër"),m.commune_normalized("Saint-Androny"))
if __name__=="__main__":unittest.main()
