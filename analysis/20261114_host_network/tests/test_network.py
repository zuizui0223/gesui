import unittest,importlib.util,pathlib,math,itertools,numpy as np
P=pathlib.Path(__file__).parents[1]/'scripts'/'network.py'
class NetworkTests(unittest.TestCase):
 def load(self):
  self.assertTrue(P.exists(),'Network calculations not implemented')
  sp=importlib.util.spec_from_file_location('net',P);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
 def test_inclusion(self):
  m=self.load()
  self.assertAlmostEqual(float(m.inclusion(4,1,2)),.5)
  self.assertAlmostEqual(float(m.inclusion(4,2,2)),5/6)
  self.assertEqual(float(m.inclusion(4,0,2)),0.)
  self.assertEqual(float(m.inclusion(4,4,2)),1.)
  with self.assertRaises(ValueError):m.inclusion(4,2,5)
 def test_expected_links_and_support(self):
  m=self.load(); b=np.array([2,2]);p=np.array([1,1,2]);A=np.array([[1,1,0],[0,0,1]],float)
  out=m.metrics(b,p,A,2,2)
  self.assertAlmostEqual(out['expected_links'],(5/6)*(1/2+1/2+5/6))
  self.assertAlmostEqual(out['host_support'],5/6)
  self.assertAlmostEqual(out['catalogue_coverage'],1.)
  self.assertEqual(out['candidate_links'],3)
 def test_unknown_catalogue_not_zero_support(self):
  m=self.load();b=np.array([2,2]);p=np.array([4,0]);A=np.array([[1,0],[0,0]],float)
  out=m.metrics(b,p,A,2,2)
  self.assertAlmostEqual(out['host_support'],1.)
  self.assertAlmostEqual(out['catalogue_coverage'],.5)
 def test_degree_preservation(self):
  m=self.load();A=np.array([[1,0,1,0],[0,1,1,1]],float);strata=np.array(['a','a','b','b'])
  perm=m.shuffle_columns(strata,np.random.default_rng(5))
  self.assertTrue(np.array_equal(strata[perm],strata))
  self.assertTrue(np.array_equal(A[:,perm].sum(1),A.sum(1)))
  self.assertTrue(np.array_equal(np.sort(A[:,perm].sum(0)),np.sort(A.sum(0))))
 def test_rarefaction_bruteforce(self):
  m=self.load(); records=[0,0,1,2,2]
  for k in range(1,6):
   combos=list(itertools.combinations(range(5),k))
   brute=sum(len({records[t] for t in c}) for c in combos)/len(combos)
   self.assertAlmostEqual(m.inclusion(5,np.array([2,1,2]),k).sum(),brute,places=12)
if __name__=='__main__':unittest.main()
