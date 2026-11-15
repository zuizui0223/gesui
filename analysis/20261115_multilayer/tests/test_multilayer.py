import sys
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from multilayer import check_grain, make_fold_ids, ridge_oof, paired_deltas, add_country_method

class TestGrainAndSplit(unittest.TestCase):
    def test_duplicate_H7_period_is_rejected(self):
        d = pd.DataFrame({'h':['a','a'], 'h6':['x','x'], 'period':['2017-2019','2017-2019']})
        with self.assertRaises(ValueError):
            check_grain(d)
    def test_joint_country_method_factor_is_unique_and_full_rank(self):
        d=pd.DataFrame({'country':['FR','FR','IT','IT'], 'method':['I2M2','I2M2','Other','Other'],
                        'period':['2017-2019','2020-2022','2017-2019','2020-2022'],
                        'area':[1.0,2.0,3.0,4.5]})
        q=add_country_method(d)
        self.assertEqual(set(q['country_method']),{'FR::I2M2','IT::Other'})
        from multilayer import _matrix, CATEGORICAL
        X,_=_matrix(q,q,['area'],CATEGORICAL)
        self.assertEqual(int(np.linalg.matrix_rank(X)),X.shape[1])
    def test_same_H6_always_same_fold(self):
        d = pd.DataFrame({'h':['a','b','c','a'], 'h6':['x','x','y','x'],
                          'period':['2017-2019','2017-2019','2017-2019','2020-2022']})
        folds=make_fold_ids(d['h6'])
        self.assertEqual(folds[0],folds[1])
        self.assertEqual(folds[0],folds[3])
    def test_train_only_encoding_yields_complete_oof(self):
        rows=[]
        for i in range(25):
            for t in ('2017-2019','2020-2022'):
                rows.append(dict(h=f'h{i}',h6=f'g{i}',period=t,country='FR' if i%2 else 'IT',
                                 method='method1',crop=i/5,W=i/10,plant_S50=20+i/6,
                                 y=1+i/20+(t=='2020-2022')*.08))
        d=pd.DataFrame(rows)
        p=ridge_oof(d,'y',['crop'],['country','method','period'],['W'],alpha=10)
        self.assertEqual(p.shape,(50,))
        self.assertTrue(np.isfinite(p).all())
    def test_pairing_drops_nonrepeat_basins_and_preserves_method(self):
        d=pd.DataFrame({'h':['a','a','b','c','c'],'h6':['x','x','y','z','z'],
            'period':['2017-2019','2020-2022','2017-2019','2017-2019','2020-2022'],
            'method':['A','A','A','A','B'],'country':['FR','FR','FR','IT','IT'],'primary_record_mean':[1,2,1,2,4],
            'plant_S50':[2,3,1,4,3], 'EQR_original_site_equal_basin_period_mean':[.3,.4,.5,.6,.7]})
        x=paired_deltas(d)
        self.assertEqual(len(x),2)
        self.assertEqual(int(x.same_method.sum()),1)
        self.assertEqual(float(x.loc[x.h=='a','d_HNS'].iloc[0]),1.0)

if __name__=='__main__':unittest.main()
