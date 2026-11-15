import unittest
import tempfile
import zipfile
import os
import sys
from pathlib import Path
import pandas as pd
from types import SimpleNamespace
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from source_join import load_locked_source, join_source, derive_exposure, examine_source_eligibility

class SourceJoinTest(unittest.TestCase):
    def setUp(self):
        self.root=Path(__file__).resolve().parents[1]
        self.source=Path(os.environ.get('GESUI_NP_SOURCE_INPUT',str(self.root/'inputs/gesui_h7_actual_reported_np_source.zip')))
        self.panel=Path(os.environ.get('GESUI_HNS_PANEL',str(self.root/'inputs/H7_454_HNS_panel.csv')))

    def test_source_read_validates_original_source_sha(self):
        src,qa=load_locked_source(self.source)
        self.assertEqual(len(src),454)
        self.assertEqual(qa['source_sha256'],'012ad79feabc8c11e0d74f803730b5b7913b96a64a6af26fb0c97263d9aa7c26')

    def test_locked_join_is_one_to_one_and_keeps_unknown_reporting(self):
        src,qa=load_locked_source(self.source)
        panel=pd.read_csv(self.panel,dtype={'h':str,'h6':str})
        x,report=join_source(panel,src)
        self.assertEqual(len(x),454)
        self.assertEqual((x.reported_positive==1).sum(),56)
        self.assertEqual(x.loc[x.reported_positive==0,'reporting_class'].unique().tolist(),['no_matched_positive_report'])
        self.assertEqual(x.loc[x.reported_positive==1,'reporting_class'].unique().tolist(),['matched_positive_report'])
        self.assertEqual(report['country_mode_mismatch_n'],11)
        self.assertEqual(report['flow_value_max_abs_diff'],1.1368683772161603e-13)

    def test_source_raw_json_matches_downloaded_archive(self):
        if self.source.suffix.lower()=='.json':
            self.skipTest('Using raw original publisher artifact in CI')
        with zipfile.ZipFile(self.source) as z:
            payload=z.read('EEA2018_EEA2022_454_H7_ORIGINAL_REPORTED_NP_SITE_FLOW_EQUIVALENT_SOURCE_ONLY.json')
        with tempfile.TemporaryDirectory() as temp:
            raw=Path(temp)/'EEA2018_EEA2022_454_H7_ORIGINAL_REPORTED_NP_SITE_FLOW_EQUIVALENT_SOURCE_ONLY.json'
            raw.write_bytes(payload)
            fromraw,qa=load_locked_source(raw)
            fromzip,qz=load_locked_source(self.source)
            self.assertEqual(qa['source_sha256'],qz['source_sha256'])
            self.assertEqual(len(fromraw),len(fromzip))

    def test_missing_report_is_not_biological_zero(self):
        frame=pd.DataFrame([{'n_source':0,'n_mg_l':0.0,'p_mg_l':0.0},
                            {'n_source':1,'n_mg_l':0.25,'p_mg_l':0.001}])
        out=derive_exposure(frame)
        self.assertEqual(out.loc[0,'reporting_class'],'no_matched_positive_report')
        self.assertEqual(out.loc[1,'reporting_class'],'matched_positive_report')
        self.assertNotIn('true_wastewater_free',set(out.reporting_class))

    def test_sample_gate(self):
        src,_=load_locked_source(self.source)
        panel=pd.read_csv(self.panel,dtype={'h':str,'h6':str})
        x,_=join_source(panel,src)
        qc=examine_source_eligibility(x)
        self.assertEqual(qc['n_positive_periods'],{'2017-2019':28,'2020-2022':28})
        self.assertEqual(qc['positive_H7'],38)
        self.assertEqual(qc['positive_H6'],31)
        self.assertEqual(qc['positive_rows'],56)
        self.assertEqual(len(qc['positive_country_counts']),12)
        self.assertTrue(qc['exploratory_support_gate_pass'])

if __name__=='__main__':unittest.main()
