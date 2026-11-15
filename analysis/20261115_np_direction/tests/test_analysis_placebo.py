import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze_placebo import group_error
class TestNegativeControl(unittest.TestCase):
    def test_identical_errors_zero_gain(self):
        a=np.array([1.,3.,2.,5.]);r=group_error(a,a,['g1','g1','g2','g2'])
        self.assertAlmostEqual(r['gain'],0.,places=12)
        self.assertEqual(r['n_H6'],2)
    def test_better_fits_positive_gain(self):
        a=np.array([2.,4.,6.,8.]);b=a*0.75
        r=group_error(a,b,['g1','g1','g2','g2'])
        self.assertAlmostEqual(r['gain'],25.,places=10)
    def test_cluster_bootstrap_consistent(self):
        a=np.array([2.,4.,6.,8.]);b=a*0.5
        r=group_error(a,b,['g1','g1','g2','g2'])
        self.assertTrue(all(abs(t-50)<1e-12 for t in r['95_interval']))
if __name__=='__main__':unittest.main()
