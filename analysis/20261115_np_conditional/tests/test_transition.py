import unittest
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from source_transition_audit import audit_report_transitions

class TestOriginalSourceTransition(unittest.TestCase):
 def test_duplicate_native_unit_rejected(self):
  x=pd.DataFrame([{'h':'1','h6':'a','period':'2017-2019','country':'X', 'upstream_report_positive':0,'downstream_report_positive':0}]*454)
  with self.assertRaises(ValueError): audit_report_transitions(x)
 def test_invalid_binary_rejected(self):
  x=pd.DataFrame([{'h':str(i),'h6':'a','period':'2017-2019','country':'X', 'upstream_report_positive':2,'downstream_report_positive':0} for i in range(454)])
  with self.assertRaises(ValueError): audit_report_transitions(x)
if __name__=='__main__':unittest.main()
