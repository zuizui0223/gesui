import sys
import unittest
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'analysis/20261122_butterfly_within_genus'))
try:
    from inference import fit_conditional, cluster_covariance, permute_trait_labels
except ImportError:
    fit_conditional=cluster_covariance=permute_trait_labels=None


class InferenceTests(unittest.TestCase):
    def test_mle_matches_closed_form_conditional_pair_odds(self):
        self.assertIsNotNone(fit_conditional,'conditional optimizer missing')
        x=np.tile(np.array([[[-1.0],[1.0]]]),(100,1,1))
        y=np.array([[0,1]]*75+[[1,0]]*25)
        fit=fit_conditional(x,y,np.repeat(2,100),np.ones(100))
        self.assertTrue(fit['converged'])
        self.assertAlmostEqual(fit['params'][0],0.5*np.log(3),places=6)

    def test_cluster_score_summed_before_outer_product(self):
        self.assertIsNotNone(cluster_covariance,'cluster covariance missing')
        # H=2; region scores are 1+2=3 and -3. Meat=18; CR1=2.
        cov=cluster_covariance(np.array([[2.0]]),np.array([[1.0],[2.0],[-3.0]]),['a','a','b'])
        self.assertAlmostEqual(float(cov[0,0]),9.0,places=12)

    def test_weighted_mle_agrees_with_integer_repetition(self):
        self.assertIsNotNone(fit_conditional,'conditional optimizer missing')
        x=np.tile(np.array([[[-1.0],[1.0]]]),(2,1,1))
        y=np.array([[0,1],[1,0]])
        fit=fit_conditional(x,y,np.array([2,2]),np.array([3.0,1.0]))
        self.assertAlmostEqual(fit['params'][0],0.5*np.log(3),places=6)

    def test_trait_labels_are_global_and_stay_within_butterfly_genus(self):
        self.assertIsNotNone(permute_trait_labels,'trait permutation missing')
        original={'a':1.0,'b':4.0,'c':8.0,'d':3.0,'e':7.0}
        genus={'a':'Alpha','b':'Alpha','c':'Alpha','d':'Beta','e':'Beta'}
        rng=np.random.default_rng(127)
        for _ in range(20):
            s=permute_trait_labels(original,genus,rng)
            self.assertEqual(sorted(s[k] for k in 'abc'),[1.0,4.0,8.0])
            self.assertEqual(sorted(s[k] for k in 'de'),[3.0,7.0])
            self.assertAlmostEqual((s['a']-s['b'])+(s['b']-s['c']),s['a']-s['c'])


if __name__=='__main__':
    unittest.main()
