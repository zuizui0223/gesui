import unittest
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from conditional_source import source_cells, paired_h6_gain, validate_equal_sample

class SourceConditionalTest(unittest.TestCase):
 def test_four_source_cells_are_exhaustive(self):
  u=np.array([0,0,1,1]);d=np.array([0,1,0,1])
  self.assertEqual(source_cells(u,d),{'neither':1,'up_only':1,'down_only':1,'both':1})
 def test_shared_denominator_model_gain(self):
  h=np.array(['a','a','b','b'])
  down=np.array([2.,6.,4.,8.]);both=np.array([1.,3.,2.,4.])
  ret=paired_h6_gain(down,both,h,seed=1,rep=80)
  self.assertAlmostEqual(ret['h6_mse_gain_pct'],50.)
  self.assertTrue(np.allclose(ret['bootstrap95'],[50,50]))
 def test_equal_model_errors_zero_gain(self):
  h=np.array(['a','a','b','b'])
  e=np.array([2.,6.,4.,8.]);r=paired_h6_gain(e,e,h)
  self.assertAlmostEqual(r['h6_mse_gain_pct'],0)
 def test_group_equal_weight_not_row_equal(self):
  # First group has 3 rows and the second only 1.
  h=np.array(['a','a','a','b'])
  e=np.array([1.,1.,1.,10.]);f=np.array([0.,0.,0.,10.])
  r=paired_h6_gain(e,f,h)
  self.assertAlmostEqual(r['h6_mse_gain_pct'],100*(5.5-5)/5.5)
 def test_any_nan_and_group_duplication_rejected(self):
  with self.assertRaises(ValueError):
   paired_h6_gain([2.,float('nan')],[1.,1.],['a','b'])
  with self.assertRaises(ValueError):
   validate_equal_sample([('1','A'),('1','A')],expected_n=2)
if __name__=='__main__':unittest.main()
