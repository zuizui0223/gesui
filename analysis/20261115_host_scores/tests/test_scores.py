import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import score_core as m

class TestScores(unittest.TestCase):
 def test_missing_not_zero(self):
  x=m.summarize([('a',8.0),('b',None)])
  self.assertEqual(x['mean'],8.0); self.assertEqual(x['coverage'],.5)
  self.assertEqual(x['lower_bound'],4.5);self.assertEqual(x['upper_bound'],8.5)
 def test_empty(self):
  x=m.summarize([]);self.assertIsNone(x['mean']);self.assertIsNone(x['coverage'])
 def test_singleton_valid(self):
  x=m.summarize([('a',3.0)]);self.assertEqual(x['sd'],0);self.assertEqual(x['band'],'low')
 def test_dedup(self):
  x=m.summarize([('a',8.0),('a',8.0),('b',4.0)])
  self.assertEqual(x['n_hosts'],2);self.assertEqual(x['mean'],6.0)
 def test_conflicting_duplicate(self):
  with self.assertRaises(ValueError):m.summarize([('a',2.0),('a',8.0)])
 def test_invalid_zero(self):
  with self.assertRaises(ValueError):m.summarize([('a',0)])
 def test_missing_and_indifferent(self):
  self.assertIsNone(m.numeric_eiv('x'));self.assertIsNone(m.numeric_eiv('NA'))
  self.assertIsNone(m.numeric_eiv(''));self.assertEqual(m.numeric_eiv('7.4'),7.4)
 def test_alias_conflict(self):
  d={'taxon a':{'N':2.0},'taxon b':{'N':8.0}}
  self.assertEqual(m.match_names(['Taxon a','Taxon b'],d)['status'],'ambiguous')
 def test_exact_and_unknown(self):
  d={'taxon a':{'N':2.0}}
  self.assertEqual(m.match_names(['Taxon a'],d)['names'],['taxon a'])
  self.assertEqual(m.match_names(['Taxon b'],d)['status'],'unmatched')
 def test_middle_not_no_effect(self):
  x=m.summarize([('a',5.0)]);self.assertEqual(x['band'],'intermediate')
class TestUsedAliases(unittest.TestCase):
 def test_unique_used_names_are_not_ambiguous(self):
  self.assertEqual(m.ambiguous_used_keys([{'butterflyKey':'1','butterfly_names':['A a']},{'butterflyKey':'2','butterfly_names':['B b']}]),set())
 def test_repeated_edges_same_key_not_ambiguous(self):
  self.assertEqual(m.ambiguous_used_keys([{'butterflyKey':'1','butterfly_names':['A a']},{'butterflyKey':'1','butterfly_names':['a a']}]),set())
 def test_conflicting_used_name_marks_both(self):
  self.assertEqual(m.ambiguous_used_keys([{'butterflyKey':'1','butterfly_names':['A a']},{'butterflyKey':'2','butterfly_names':['a a']}]),{'1','2'})
if __name__=='__main__': unittest.main()
