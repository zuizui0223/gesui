import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "analysis/20261118_point_diffuse_source/jrc_soilN_pointN_nuts2_source_gate.py"
spec = importlib.util.spec_from_file_location("jrc_soilN_gate", SCRIPT)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)


class SoilNSourceOnlyTests(unittest.TestCase):
    def test_distinguish_valid_zero_and_missing(self):
        a,q=g.parse_rows(
            [{"cntrCode":"DE","Nuts2":"DE11","PsN":0},
             {"cntrCode":"DE","Nuts2":"DE12","PsN":None}],
            "cntrCode","Nuts2","PsN")
        self.assertEqual(a["DE|DE11"],0)
        self.assertIsNone(a["DE|DE12"])
        self.assertEqual(q["true_zero_regions"],1)

    def test_duplicate_region_fail(self):
        with self.assertRaises(ValueError):
            g.parse_rows([{"country":"DE","region":"DE11","x":1},
                          {"country":"DE","region":"DE11","x":2}],
                         "country","region","x")

    def test_different_models_do_not_fill_missing(self):
        a={key:{"FR|FR11":1.,"FR|FR12":2.} for key in g.SOURCES}
        del a["land_atmos_N_mg_m2"]["FR|FR12"]
        joined, gate=g.join(a)
        self.assertTrue(joined[0]["complete_all_six_modelled_sources"])
        self.assertFalse(joined[1]["complete_all_six_modelled_sources"])
        self.assertFalse(gate["crosswalk_source_only_pass"])

    def test_nonnegative_only(self):
        self.assertIsNone(g.numeric(-1))
        self.assertIsNone(g.numeric(float("nan")))
        self.assertEqual(g.numeric(0),0)


if __name__=="__main__":
    unittest.main()
