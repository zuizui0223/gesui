import unittest
import numpy as np
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
try:
 from score_core import aggregate_traits, blocked_splits, ridge_fit_predict
except ImportError:
 aggregate_traits=blocked_splits=ridge_fit_predict=None

class ScoreTests(unittest.TestCase):
 def test_implementation_exists(self):
  self.assertIsNotNone(aggregate_traits, 'Trait aggregation has not been implemented')
 def calc(self,b=None):
  self.assertIsNotNone(aggregate_traits)
  return aggregate_traits(np.array([[2,1,7],[0,3,0]]) if b is None else b,np.array([2.,8.,np.nan]),np.array(['A','B','A']),np.array([1.,.5,np.nan]),np.array([12,3]))
 def test_weighted_mean_does_not_impute_missing(self):
  r=self.calc();np.testing.assert_allclose(r['record_mean'],[4.,8.]);np.testing.assert_allclose(r['incidence_mean'],[5.,8.])
 def test_coverage_keeps_unresolved_and_unscored_denominators(self):
  r=self.calc();np.testing.assert_allclose(r['record_coverage_raw'],[.25,1.]);np.testing.assert_allclose(r['species_coverage'],[2/3,1.])
 def test_no_scored_record_stays_nan(self):
  r=self.calc(np.array([[0,0,7],[0,0,0]]));self.assertTrue(np.isnan(r['record_mean']).all());np.testing.assert_equal(r['record_coverage_raw'],[0.,0.])
 def test_centered_decomposition_exact(self):
  r=self.calc();np.testing.assert_allclose(r['record_mean'],r['between_family']+r['within_family'])
 def test_negative_counts_rejected(self):
  with self.assertRaises(ValueError):self.calc(np.array([[-1,1,7],[0,3,0]]))
 def test_raw_denominator_too_small_rejected(self):
  with self.assertRaises(ValueError):aggregate_traits(np.array([[20,1,7],[0,3,0]]),np.array([2.,8.,np.nan]),np.array(['A','B','A']),np.array([1.,.5,np.nan]),np.array([12,3]))
 def test_family_shares_sum_one(self):
  r=self.calc();np.testing.assert_allclose(r['record_family_shares'].sum(1),[1.,1.])
 def test_permuted_species_order_invariant(self):
  r=self.calc();ix=[2,0,1];q=aggregate_traits(np.array([[2,1,7],[0,3,0]])[:,ix],np.array([2.,8.,np.nan])[ix],np.array(['A','B','A'])[ix],np.array([1.,.5,np.nan])[ix],np.array([12,3]));np.testing.assert_allclose(r['record_mean'],q['record_mean'])
class ModelTests(unittest.TestCase):
 def test_h6_is_always_disjoint(self):
  self.assertIsNotNone(blocked_splits);h=np.array(['a','a','b','c','d','e','f']);c=np.array(['FR','FR','DE','DE','FR','IT','IT'])
  for _,tr,te in blocked_splits(h,c,'H6'):
   self.assertFalse(set(h[tr])&set(h[te]))
 def test_country_purges_entire_shared_h6(self):
  self.assertIsNotNone(blocked_splits);h=np.array(['a','a','b','c','d','e']);c=np.array(['FR','DE','DE','FR','IT','IT'])
  for g,tr,te in blocked_splits(h,c,'country'):
   self.assertFalse(set(h[tr])&set(h[te]));self.assertEqual(set(c[te]),{g})
 def test_ridge_intercept_unpenalized(self):
  self.assertIsNotNone(ridge_fit_predict);X=np.column_stack([np.ones(4),np.arange(4)]);y=np.full(4,5.);p=ridge_fit_predict(X,y,X,np.ones(4));np.testing.assert_allclose(p,5.,atol=1e-12)
 def test_weighted_ridge_matches_direct_objective(self):
  self.assertIsNotNone(ridge_fit_predict);X=np.column_stack([np.ones(4),np.arange(4)]);y=np.array([2.,1.,5.,3.]);w=np.array([1.,2.,3.,4.]);expected=X@np.linalg.solve(X.T@(w[:,None]*X)+np.diag([0.,10.]),X.T@(w*y));np.testing.assert_allclose(ridge_fit_predict(X,y,X,w),expected)
if __name__=='__main__':unittest.main()

class OutputSchemaTests(unittest.TestCase):
 def test_group_count_and_model_metrics_have_distinct_keys(self):
  import pandas as pd
  from model_core import model,numeric_columns
  rng=np.random.default_rng(43);n=100
  d=pd.DataFrame(rng.normal(size=(n,len(numeric_columns()))),columns=numeric_columns())
  d['h']=[f'h{i}' for i in range(n)];d['h6']=[f'parent{i//5}' for i in range(n)]
  d['country']='FR';d['country_holdout']='FR';d['period']='p';d['method']='m';d['W']=rng.random(n);d['y']=rng.normal(size=n)
  result=model(d,'y',country=False)
  self.assertEqual(result.get('n_H6'),20)
  self.assertIsInstance(result['H6'],dict)
