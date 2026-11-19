import importlib.util
from pathlib import Path
import unittest

p=Path(__file__).resolve().parents[1]/"analysis/20261119_jrc_native_basins/jrc_native_full_source_area_census.py"
spec=importlib.util.spec_from_file_location("nativejrcfull",p)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class NativeFullSourceTests(unittest.TestCase):
    def test_integer_ID_alone_not_proof(self):
        self.assertFalse(m.compatible_areas(120,300))
        self.assertTrue(m.compatible_areas(100,100.00001))
        self.assertFalse(m.compatible_areas(None,100))
    def test_all_native_source_values_necessary(self):
        soilA={i:(1.,100.+i) for i in range(130)}
        soilM={i:(5.,100.+i) for i in range(130)}
        soilO={i:(4.,100.+i) for i in range(130)}
        soilF={i:(2.,100.+i) for i in range(130)}
        point={i:(10.,100.+i) for i in range(130)}
        cand,stat=m.check_land_and_candidates({
            "soil_atmosN_2015":soilA,"soil_mineralN_2010_19":soilM,
            "soil_organicN_2010_19":soilO,"soil_fixationN_2010_19":soilF,
            "water_pointN_2014_18":point
        })
        self.assertEqual(len(cand),130)
        self.assertEqual(stat["soil_3main_values_nonmissing_hydroID"],130)
        self.assertFalse(stat["point_and_soil_fully_geospatially_proved_same_HydroID"])
    def test_missing_not_zero_and_area_disagreement(self):
        a={i:(1.,100.) for i in range(4)}
        b={i:(2.,100.) for i in range(4)}
        c={i:(3.,100.) for i in range(4)}
        point={i:(10.,100.) for i in range(4)}
        a[2]=(None,100.)
        b[3]=(2.,200.)
        out,report=m.check_land_and_candidates({
          "soil_atmosN_2015":a,"soil_mineralN_2010_19":b,
          "soil_organicN_2010_19":c,"soil_fixationN_2010_19":{},
          "water_pointN_2014_18":point})
        self.assertEqual(report["soil_3main_values_nonmissing_hydroID"],3)
        self.assertEqual(report["soil_intersection_area_conflicts"],1)
        self.assertEqual(len(out),2)

if __name__=="__main__":unittest.main()
