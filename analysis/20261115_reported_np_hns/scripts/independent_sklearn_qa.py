#!/usr/bin/env python3
"""Independent sklearn numerical verification of EEA-reporting × three H7 outcomes.

Does not import the original manual Ridge fitter or its feature matrix builder.
"""
import argparse,collections,hashlib,json,warnings
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
from source_join import load_locked_source,join_source,sha256,FROZEN_PANEL

BASE=['crop','urban','riverflow','naturalout','area','nEEA','nAnimal','nPlant','nButterfly','W']
HNS_EXTRAS=['primary_record_coverage_raw','primary_host_completeness','primary_rf_Lycaenidae','primary_rf_Pieridae','primary_rf_Papilionidae','primary_rf_Hesperiidae']
CATS=['country_method','period']
RECIPES={
 'B0_static_W_only':[],
 'B1_plus_matched_source_presence':['reported_positive'],
 'B2_plus_N_report_index':['log_reported_N_mgL_equivalent'],
 'B3_plus_P_report_index':['log_reported_P_ugL_equivalent'],
 'B4_plus_presence_and_N':['reported_positive','log_reported_N_mgL_equivalent'],
 'B5_plus_presence_and_P':['reported_positive','log_reported_P_ugL_equivalent'],
 'B6_plus_presence_and_N_and_P':['reported_positive','log_reported_N_mgL_equivalent','log_reported_P_ugL_equivalent'],
}

def sklearn_oof(df,target,features):
    h6=df.h6.astype(str).to_numpy()
    folds=np.array([int(hashlib.sha256(q.encode()).hexdigest()[:8],16)%5 for q in h6])
    y=df[target].to_numpy(dtype=float);pred=np.full(len(df),np.nan)
    for fold in range(5):
        tr=folds!=fold;te=~tr
        fit=df.iloc[np.flatnonzero(tr)];val=df.iloc[np.flatnonzero(te)]
        scaler=StandardScaler().fit(fit[features].to_numpy(dtype=float))
        enc=OneHotEncoder(drop='first',handle_unknown='ignore',sparse_output=False).fit(fit[CATS])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            X=np.column_stack((scaler.transform(fit[features]),enc.transform(fit[CATS])))
            Xtest=np.column_stack((scaler.transform(val[features]),enc.transform(val[CATS])))
        ct=collections.Counter(h6[tr]);w=np.array([1/ct[g] for g in h6[tr]],dtype=float);w/=w.mean()
        reg=Ridge(alpha=10,fit_intercept=True,solver='svd')
        reg.fit(X,y[tr],sample_weight=w)
        pred[te]=reg.predict(Xtest)
    if not np.isfinite(pred).all():raise RuntimeError('sklearn predictions incomplete')
    return pred

def main(panel,source_zip,source_results,outpath):
    if sha256(panel)!=FROZEN_PANEL:raise RuntimeError('Original HNS source SHA drift')
    source,metadata=load_locked_source(source_zip)
    orig=pd.read_csv(panel,dtype={'h':str,'h6':str})
    data,qa=join_source(orig,source)
    data['country_method']=data.country.fillna('__UNKNOWN__').astype(str)+'::'+data.method.fillna('__UNKNOWN__').astype(str)
    earlier=json.loads(source_results.read_text())
    report={'schema':'GESUI_20261115_EEA_NP_H7_INDEPENDENT_SKLEARN_QA_v1',
            'source_SHA':metadata['source_sha256'],'original_panel_SHA':FROZEN_PANEL,
            'n_source_rows':len(data),'n_H6':data.h6.nunique(),
            'method':'separate sklearn StandardScaler,OneHotEncoder,Ridge with H6 train-only fits, group-equal MSE',
            'results':{}, 'max_absolute_MSE_gain_deviation':None,'success':False}
    diffs=[]
    for label,(outcome,vars_) in {
         'HNS':('primary_record_mean',BASE+HNS_EXTRAS),
         'EQR':('EQR_original_site_equal_basin_period_mean',BASE),
         'PlantS50':('plant_S50',BASE),
    }.items():
        y=data[outcome].to_numpy(dtype=float);h6=data.h6.astype(str).to_numpy();names=sorted(set(h6))
        errors={}
        for model,extras in RECIPES.items():
            pred=sklearn_oof(data,outcome,vars_+extras)
            errors[model]=(y-pred)**2
        report['results'][label]={}
        for model in list(RECIPES)[1:]:
            baseline='B1_plus_matched_source_presence' if model in ('B4_plus_presence_and_N','B5_plus_presence_and_P','B6_plus_presence_and_N_and_P') else 'B0_static_W_only'
            aa=np.array([errors[baseline][h6==g].mean() for g in names]);bb=np.array([errors[model][h6==g].mean() for g in names])
            gain=100*(aa.mean()-bb.mean())/aa.mean()
            observed=earlier['models'][label]['comparisons'][model]['H6_equal_MSE_gain_percent']
            diff=abs(gain-observed)
            diffs.append(diff)
            report['results'][label][model]={'gain_pct':float(gain),'manual_original_pct':float(observed),'absolute_difference':float(diff)}
    report['max_absolute_MSE_gain_deviation']=max(diffs)
    report['success']=bool(max(diffs)<1e-7 and len(diffs)==18)
    outpath.parent.mkdir(parents=True,exist_ok=True)
    outpath.write_text(json.dumps(report,indent=2,sort_keys=True))
    print(json.dumps({'success':report['success'],'comparisons':len(diffs),'max_diff':max(diffs)},sort_keys=True))
    if not report['success']:raise SystemExit('Independent numerical models differ from prior frozen source computation')

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--panel',type=Path,required=True)
    ap.add_argument('--source-zip',type=Path,required=True)
    ap.add_argument('--results',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    main(args.panel,args.source_zip,args.results,args.out)
