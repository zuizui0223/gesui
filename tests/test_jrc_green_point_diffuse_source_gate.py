import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "analysis/20261118_point_diffuse_source/jrc_green_point_diffuse_source_gate.py"
spec = importlib.util.spec_from_file_location("jrc_green_source", SCRIPT)
src = importlib.util.module_from_spec(spec)
spec.loader.exec_module(src)


class SourceOnlyTest(unittest.TestCase):
    def fixture(self, miss=False):
        pN = [dict(BasinID=1, PsN=4, AreaSqKm=10), dict(BasinID=2, PsN=0, AreaSqKm=20)]
        pP = [dict(BasinID=1, PsP=2, AreaSqKm=10), dict(BasinID=2, PsP=0, AreaSqKm=20)]
        dN = [dict(BasinID=1, HydroID=101, De=2), dict(BasinID=1, HydroID=102, De=3),
              dict(BasinID=2, HydroID=103, De=None if miss else 0)]
        dP = [dict(BasinID=1, HydroID=101, De=1), dict(BasinID=1, HydroID=102, De=2),
              dict(BasinID=2, HydroID=103, De=0)]
        return {k: src.normalize(rows, kind, nutrient) for k, rows, kind, nutrient in [
            ("point_N", pN, "point", "N"), ("point_P", pP, "point", "P"),
            ("diffuse_N", dN, "diffuse", "N"), ("diffuse_P", dP, "diffuse", "P")
        ]}

    def test_zero_is_not_missing(self):
        panel, qc = src.make_panel(self.fixture())
        self.assertEqual(panel[0]["diffuse_N_modelled_tonnes_2014_2018"], 5)
        self.assertEqual(panel[1]["diffuse_N_modelled_tonnes_2014_2018"], 0)
        self.assertTrue(panel[1]["full_four_source_support"])
        self.assertFalse(qc["source_only_native_BasinID_crosswalk_pass"])

    def test_missing_catchment_propagates(self):
        panel, _ = src.make_panel(self.fixture(miss=True))
        self.assertIsNone(panel[1]["diffuse_N_modelled_tonnes_2014_2018"])
        self.assertFalse(panel[1]["full_four_source_support"])

    def test_duplicate_point_fails(self):
        with self.assertRaises(ValueError):
            src.normalize([dict(BasinID=1, PsN=1, AreaSqKm=10),
                           dict(BasinID=1, PsN=2, AreaSqKm=10)], "point", "N")

    def test_duplicate_catchment_fails(self):
        with self.assertRaises(ValueError):
            src.normalize([dict(BasinID=1, HydroID=7, De=1),
                           dict(BasinID=1, HydroID=7, De=2)], "diffuse", "N")

    def test_area_disagreement_fails(self):
        layers = self.fixture()
        layers["point_P"][1][1] = 12
        with self.assertRaises(ValueError):
            src.make_panel(layers)


if __name__ == "__main__":
    unittest.main()
