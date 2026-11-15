import importlib.util, pathlib, sys, unittest
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))

class RobustnessTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('robustness_core'), 'Missing new robustness implementation')
        import robustness_core
        return robustness_core
    def test_h6_equal_weights(self):
        m=self.module();w=m.group_weights(np.array(['a','a','b']))
        np.testing.assert_allclose(w,[.25,.25,.5])
    def test_columns_never_restore_removed_host(self):
        m=self.module();rng=np.random.default_rng(1)
        keep=np.array([False,True,True,True]);s=np.array(['x','x','x','y'])
        q=m.retained_permutation(s,keep,rng)
        self.assertEqual(set(q),{1,2,3});self.assertNotIn(0,q)
        np.testing.assert_array_equal(s[q],s[keep])
    def test_family_partition(self):
        m=self.module();pb=np.array([[.5,.3],[.1,.8]]);pp=np.array([[.2,.7],[.4,.1]])
        A=np.array([[1,0],[1,1]])
        total=m.edge_contributions(pb,pp,A)
        separate=sum(m.edge_contributions(pb[:,[i]],pp,A[[i]]) for i in range(2))
        np.testing.assert_allclose(total,separate)
        np.testing.assert_allclose(total.sum(1),[.37,.44])
    def test_delete_host_keeps_rarefaction_probabilities(self):
        m=self.module();pb=np.array([[.5,.3]]);pp=np.array([[.2,.7]]);A=np.array([[1,0],[1,1]])
        v=m.edge_contributions(pb,pp[:,[1]],A[:,[1]])
        self.assertAlmostEqual(v.sum(),.21)
    def test_null_degrees_preserved(self):
        m=self.module();A=np.array([[1,0,1,0],[0,1,0,1]])
        k=np.array([True,False,True,True]);s=np.array(['x','x','x','x'])
        q=m.retained_permutation(s,k,np.random.default_rng(44))
        np.testing.assert_array_equal(A[:,q].sum(1),A[:,k].sum(1))
    def test_empty_subset_rejected(self):
        m=self.module()
        with self.assertRaises(ValueError):m.group_weights(np.array([]))

class CompositionTests(unittest.TestCase):
    def test_source_hellinger(self):
        import robustness_core as m
        self.assertTrue(hasattr(m,'hellinger'), 'Missing composition control')
        np.testing.assert_allclose(m.hellinger(np.array([[3,1],[0,2]])),[[np.sqrt(.75),.5],[0,1]])
    def test_absent_family_not_invented(self):
        import robustness_core as m
        self.assertTrue(hasattr(m,'hellinger'), 'Missing composition control')
        self.assertTrue(np.isnan(m.hellinger(np.zeros((1,2)))).all())

if __name__=='__main__':unittest.main()
