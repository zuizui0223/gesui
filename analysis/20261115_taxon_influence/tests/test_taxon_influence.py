
import unittest
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from taxon_influence import cwm_from_matrix, ranked_taxa

class TaxonInfluenceTests(unittest.TestCase):
    def test_no_drop_recovers_weighted_mean(self):
        B=np.array([[2,1,1],[0,2,2]])
        score=np.array([2.,8.,np.nan])
        actual=cwm_from_matrix(B,score,[])
        np.testing.assert_allclose(actual,[4.,8.])
    def test_drop_one_renormalizes(self):
        B=np.array([[2,1,1],[1,3,0]])
        score=np.array([2.,8.,np.nan])
        actual=cwm_from_matrix(B,score,[0])
        np.testing.assert_allclose(actual,[8.,8.])
    def test_empty_known_denominator_not_zero(self):
        B=np.array([[1,0],[0,1]])
        score=np.array([5.,np.nan])
        with self.assertRaises(ValueError):
            cwm_from_matrix(B,score,[0])
    def test_ranked_taxa_absolute_descending_and_stable_ties(self):
        names=["b","a","c"]
        contribution=[-.2,.2,.05]
        self.assertEqual(ranked_taxa(names,contribution), ["a","b","c"])

if __name__=="__main__":
    unittest.main()
