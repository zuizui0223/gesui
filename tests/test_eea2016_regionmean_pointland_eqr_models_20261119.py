import importlib.util
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

P=Path(__file__).resolve().parents[1]/"analysis/20261119_point_land_N_source/eea2016_regionmean_eqr_point_vs_landN.py"
spec=importlib.util.spec_from_file_location("e2016",P)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def synthetic():
    n=80
    return pd.DataFrame({
        "country":["C"+str(i%8) for i in range(n)],
        "nuts2":["N"+str(i) for i in range(n)],
        "fold":[i%5+1 for i in range(n)],
        "lat":[40+.1*i for i in range(n)],
        "lon":[10+.03*(i%17) for i in range(n)],
        "logarea":[3+np.log1p(i%12) for i in range(n)],
        "log_source_n":[1+np.log1p(i%9) for i in range(n)],
        "pointN":[np.log1p((i%15)+1) for i in range(n)],
        "atmos":[np.log1p(5+i%11) for i in range(n)],
        "mineral":[np.log1p(1+(i*3)%19) for i in range(n)],
        "organic":[np.log1p(1+(i*5)%29) for i in range(n)],
        "EQR":[.3+.08*np.sin(i*.27) for i in range(n)],
        "original_valid_site_n":[2+(i%9) for i in range(n)],
    })

class AquaticOnlyFourModelTests(unittest.TestCase):
    def test_ridge_matches_independent_sklearn(self):
        d=synthetic()
        a=d.iloc[:55]
        b=d.iloc[55:]
        X,Z=m.build_design(a,b,("pointN",*m.D_VARS))
        scikit=Ridge(alpha=10.,fit_intercept=True).fit(X,a.EQR.to_numpy()).predict(Z)
        manual=m.fit_predict(a,b,("pointN",*m.D_VARS))
        self.assertLess(float(np.max(np.abs(manual-scikit))),1e-10)
    def test_full_out_of_fold_predictions_recovered(self):
        d=synthetic()
        result=m.cross_validate(d)
        self.assertEqual(len(result),80)
        self.assertFalse(result[["B0_OOF_predicted_mean_EQR","BP_OOF_predicted_mean_EQR",
                                 "BD_OOF_predicted_mean_EQR","BPD_OOF_predicted_mean_EQR"]].isna().any().any())
        metrics=m.metrics(result)
        self.assertEqual(len(metrics["country_mean_OOF_squared_error_by_country"]),8)
    def test_positive_prediction_gain_definition(self):
        d=synthetic()
        res=m.metrics(m.cross_validate(d))
        gain=res["source_unique_predictions_country_equal"]
        self.assertTrue(np.isfinite(list(gain.values())).all())

if __name__=="__main__":unittest.main()
