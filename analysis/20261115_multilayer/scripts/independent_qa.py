#!/usr/bin/env python3
"""Independent sklearn estimator check of original H6-out predictive diagnostics.

Crucially does not reuse multilayer._matrix(), ridge_oof() or CV errors;
only reads frozen input, published model specs and source results.
"""
import argparse
import hashlib
import json
import collections
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge

PINNED_SHA='1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca'
CATEGORICAL=['country_method','period']
TESTED=['W_to_EQR','W_to_PlantS50','W_to_HNS','Plant_and_EQR_after_W_to_HNS']


def independ_oof(df,target,numeric):
    y=df[target].to_numpy(dtype=float)
    h6=df.h6.astype(str).to_numpy()
    folds=np.array([int(hashlib.sha256(h.encode()).hexdigest()[:8],16)%5 for h in h6])
    pred=np.full(len(df),np.nan)
    for k in range(5):
        train=folds!=k;test=~train
        if not test.any():continue
        A=df.iloc[np.flatnonzero(train)];B=df.iloc[np.flatnonzero(test)]
        scaler=StandardScaler().fit(A[numeric].to_numpy(dtype=float))
        enc=OneHotEncoder(drop='first',handle_unknown='ignore',sparse_output=False).fit(A[CATEGORICAL])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            X=np.column_stack([scaler.transform(A[numeric]),enc.transform(A[CATEGORICAL])])
            Y=np.column_stack([scaler.transform(B[numeric]),enc.transform(B[CATEGORICAL])])
        freq=collections.Counter(h6[train]);w=np.array([1/freq[h] for h in h6[train]])
        w/=w.mean()
        est=Ridge(alpha=10,fit_intercept=True,solver='svd')
        est.fit(X,y[train],sample_weight=w)
        pred[test]=est.predict(Y)
    if not np.isfinite(pred).all():raise ValueError('Missing out-of-sample predictor')
    return pred


def main(panel,outjson,dest):
    if hashlib.sha256(panel.read_bytes()).hexdigest()!=PINNED_SHA:raise ValueError('Source SHA changed')
    df=pd.read_csv(panel,dtype={'h':str,'h6':str})
    df['country_method']=df['country'].fillna('__UNKNOWN__').astype(str)+'::'+df['method'].fillna('__UNKNOWN__').astype(str)
    results=json.loads(outjson.read_text())
    report={'source_SHA256':PINNED_SHA,'method':'independent sklearn StandardScaler/OneHotEncoder/Ridge/ H6 equal OOF','models':{},'passed':True}
    for name in TESTED:
        m=results['models'][name];base=m['baseline_numeric'];extra=m['added_numeric'];yname=m['outcome']
        a=independ_oof(df,yname,base);b=independ_oof(df,yname,base+extra)
        yy=df[yname].to_numpy(dtype=float);h6=df.h6.astype(str).to_numpy()
        groups=sorted(set(h6))
        ea=np.array([np.mean((yy[h6==g]-a[h6==g])**2) for g in groups])
        eb=np.array([np.mean((yy[h6==g]-b[h6==g])**2) for g in groups])
        gain=100*(ea.mean()-eb.mean())/ea.mean()
        earlier=m['predictive_added_value']['H6_equal_MSE_gain_percent']
        diff=abs(gain-earlier)
        report['models'][name]={'n_test_rows':len(df),'new_gain':float(gain),'archived_gain':float(earlier),'difference':float(diff)}
        if diff>1e-7:report['passed']=False
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(report,indent=2,sort_keys=True))
    print(json.dumps(report,sort_keys=True))
    if not report['passed']:raise SystemExit('Independent source diagnostic disagreement')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--panel',type=Path,required=True);ap.add_argument('--results',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args();main(a.panel,a.results,a.out)
