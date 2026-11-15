import importlib.util
import sys
import unittest
import os
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from source_join import join_source,load_locked_source
from analyze_reported_np import make_source_features, refit_old_static_W, grouped_mse_comparison, paired_source_change, country_h6_disjoint_folds

ROOT=Path(__file__).resolve().parents[1]

def old_core():
    path=Path(os.environ.get('GESUI_MULTILAYER_CORE_PATH','/mnt/data/gesui_sameh7_multilayer_20261115/scripts/multilayer.py'))
    spec=importlib.util.spec_from_file_location('reference_multilayer',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class ActualAnnualSourceModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        src,_=load_locked_source(Path(os.environ.get('GESUI_NP_SOURCE_INPUT',str(ROOT/'inputs/gesui_h7_actual_reported_np_source.zip'))))
        df=pd.read_csv(Path(os.environ.get('GESUI_HNS_PANEL',str(ROOT/'inputs/H7_454_HNS_panel.csv'))),dtype={'h':str,'h6':str})
        cls.df,_=join_source(df,src)
        cls.core=old_core()

    def test_exact_original_HNS_and_EQR_W_baselines(self):
        z=refit_old_static_W(self.df,self.core)
        self.assertAlmostEqual(z['HNS'],10.014323780152122,places=9)
        self.assertAlmostEqual(z['EQR'],6.724969789376491,places=9)

    def test_reporting_feature_does_not_propose_true_control(self):
        d=make_source_features(self.df)
        self.assertEqual(d.reported_positive.sum(),56)
        self.assertTrue(np.all(d.loc[d.reported_positive==0,'log_reported_N_mgL_equivalent']==0))
        self.assertTrue((d.reported_positive==0).sum()==398)
        self.assertTrue(set(d.reporting_class)=={'no_matched_positive_report','matched_positive_report'})

    def test_h6_equal_error_calibration(self):
        d=self.df.iloc[:4].copy();d['h6']=['X','X','Y','Y'];e0=np.array([1.,1.,4.,4.]);e1=np.array([.5,.5,3.,3.])
        z=grouped_mse_comparison(e0,e1,d.h6.to_numpy(),n_resamples=300,seed=21)
        self.assertAlmostEqual(z['H6_equal_MSE_gain_percent'],100*(2.5-1.75)/2.5)
        self.assertEqual(z['n_H6'],2)

    def test_country_validation_never_leaks_parent_H6(self):
        x=pd.DataFrame({'country':['A','A','B','B'],'h6':['1','2','2','3']})
        d=country_h6_disjoint_folds(x)
        heldB=next(t for t in d if t['target_country']=='B')
        self.assertEqual(heldB['test_n'],2)
        self.assertEqual(heldB['train_n'],1)
        self.assertEqual(heldB['overlap_H6'],0)

    def test_time_source_report_status_is_not_operation_change(self):
        z=paired_source_change(self.df)
        self.assertEqual(z['n_pairs'],222)
        self.assertEqual(z['report_status_transitions'],{'0_to_0':185,'0_to_1':10,'1_to_0':9,'1_to_1':18})
        self.assertTrue(z['not_equivalent_to_WWTP_opening_or_shutdown'])

if __name__=='__main__':unittest.main()
