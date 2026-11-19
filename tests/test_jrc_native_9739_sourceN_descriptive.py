import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parents[1]/"analysis/20261119_jrc_native_basins/jrc_native_9739_sourceN_descriptive.py"
spec=importlib.util.spec_from_file_location("jrcdesc",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class JRCNativeSourceDescriptiveTests(unittest.TestCase):
    def test_missing_is_not_zero(self):
        x=[1,2,3,4,np.nan]
        y=[1,2,3,4,999]
        self.assertAlmostEqual(m.pearson(x,y),1)
        self.assertAlmostEqual(m.rankcorr(x,y),1)
    def test_constant_predictor(self):
        self.assertIsNone(m.pearson([2,2,2,2],[1,2,3,4]))
    def test_all_source_fields_described(self):
        n=40
        source=pd.DataFrame({
          "candidate_ID_number_not_proof_of_same_basin":list(range(n)),
          "soil_area_sqkm":[1+i*10 for i in range(n)],
          "point_area_sqkm":[1+i*10 for i in range(n)],
          "water_point_N_modelled_t":[0 if i%4==0 else 2*i+1 for i in range(n)],
          "soil_atmos_N_mean_mg_m2":[100+i for i in range(n)],
          "soil_mineral_N_mean_kg_ha":[5+i%7 for i in range(n)],
          "soil_organic_N_mean_kg_ha":[20+i%9 for i in range(n)],
          "soil_fixation_N_mean_kg_ha":[np.nan if i%5==0 else 3+i%4 for i in range(n)],
        })
        result=m.summarize(source)
        self.assertEqual(result["native_candidate_count"],40)
        self.assertEqual(result["native_model_pointN_true_zero_n"],10)
        self.assertEqual(result["native_soil_fixation_missing_n"],8)
        self.assertEqual(len(result["descriptive_source_correlations"]),4)
        self.assertEqual(len(result["native_all_candidate_area_quartiles"]),4)
if __name__=="__main__":
    unittest.main()
