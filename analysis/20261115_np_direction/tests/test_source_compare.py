import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import source_compare as sc
class TestSource(unittest.TestCase):
    def test_conversion_unit(self):
        self.assertAlmostEqual(sc.FACTOR,0.03168808781402895,places=12)
    def test_licensing_source_sha_cannot_be_skipped(self):
        with self.assertRaises(ValueError):sc.source_data(b'{"candidate_source_site_period_rows":[]}','mapping')
    def test_missing_and_negative_dose_not_zero(self):
        self.assertIsNone(sc.bounded(''))
        self.assertIsNone(sc.bounded('-1'))
        self.assertIsNone(sc.bounded('NaN'))
        self.assertEqual(sc.bounded('0'),0)
    def test_git_blob_hash_known(self):
        self.assertEqual(sc.gitblob(b'hello'), 'b6fc4c620b67d95f953a5c1c1230aaab5db5a1b0')
if __name__=='__main__':unittest.main()
