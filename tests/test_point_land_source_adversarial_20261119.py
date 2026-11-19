import importlib.util
from pathlib import Path
import unittest
import numpy as np
P=Path(__file__).resolve().parents[1]/"analysis/20261119_point_land_N_source/point_land_N_area_country_adversarial.py"
spec=importlib.util.spec_from_file_location("pointN",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class SourceOnlyTests(unittest.TestCase):
    def test_demean(self):
        c=["A","A","B","B","C","C"]
        x=[10,12,20,22,30,32]
        y=[50,48,100,98,150,148]
        self.assertAlmostEqual(m.within_country_r(x,y,c),-1.)
        self.assertEqual(len(m.moments_by_country(x,y,c)),3)
    def test_country_weights(self):
        w=m.country_equal_weights(["A","A","A","A","B","B"])
        self.assertAlmostEqual(sum(w[:4]),1.)
        self.assertAlmostEqual(sum(w[4:]),1.)
    def test_bootstrap_deterministic(self):
        c=["A","A","B","B","C","C"]
        x=[0,1,0,1,0,1]
        y=[1,0,0,1,1,1.1]
        mom=m.moments_by_country(x,y,c)
        self.assertEqual(m.bootstrap_country_r(mom,nboot=99),m.bootstrap_country_r(mom,nboot=99))
        self.assertEqual(len(m.country_loo(mom)),3)
    def test_model_source_partial_R2(self):
        x=np.array([1,2,4,5,7,8,9,10.])
        B=np.column_stack([np.ones(8),[1,1,1,1,0,0,0,0]])
        D=np.array([4,5,7,8,10,11,12,13.])
        out=m.fit_source_multivariate(x,[D],B)
        self.assertGreaterEqual(out["partial_R2_pointN_given_landN"],0)
        self.assertLessEqual(out["partial_R2_pointN_given_landN"],1)
    def test_constant_and_json(self):
        self.assertIsNone(m.pearson([2,2,2,2],[1,2,3,4]))
        self.assertIsNone(m.clean(float("nan")))
        self.assertEqual(m.clean(np.int64(2)),2)

if __name__=="__main__":unittest.main()
