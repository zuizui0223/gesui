#!/usr/bin/env python3
"""Independent sklearn H6-heldout verification of directional report placebo.

Verifies *predicted source-index outcome* numerical calculations only; no
independent ecological observations or biological mechanism identification.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,warnings
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
from analyze_placebo import join_direction
BASE=['crop','urban','riverflow','naturalout','area','nEEA','nAnimal','nPlant','nButterfly','W']
HNS=['primary_record_coverage_raw','primary_host_completeness','primary_rf_Lycaenidae','primary_rf_Pieridae','primary_rf_Papilionidae','primary_rf_Hesperiidae']
CATS=['country_method','period']
RECIPES={
 'static_W':[],
 'source_up_present':['uPresence'],
 'source_down_present':['dPresence'],
 'source_both_present':['uPresence','dPresence'],
 'source_up_N':['uN'],'source_down_N':['dN'],'source_both_N':['uN','dN'],
 'source_up_P':['uP'],'source_down_P':['dP'],'source_both_P':['uP','dP'],
}
def sk_oof(d,yvar,features):
    h6=d.h6.astype(str).to_numpy()
    fold=np.asarray([int(hashlib.sha256(v.encode()).hexdigest()[:8],16)%5 for v in h6])
    y=d[yvar].to_numpy(float);pred=np.full(len(d),np.nan)
    for k in range(5):
        tr=fold!=k;te=~tr
        x=d.iloc[np.flatnonzero(tr)];z=d.iloc[np.flatnonzero(te)]
        num=StandardScaler().fit(x[features].to_numpy(float))
        oh=OneHotEncoder(drop='first',handle_unknown='ignore',sparse_output=False).fit(x[CATS])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',UserWarning)
            X=np.column_stack((num.transform(x[features]),oh.transform(x[CATS])))
            Z=np.column_stack((num.transform(z[features]),oh.transform(z[CATS])))
        weight=collections.Counter(h6[tr]);w=np.array([1/weight[g] for g in h6[tr]],dtype=float);w/=w.mean()
        model=Ridge(alpha=10.,fit_intercept=True,solver='svd')
        model.fit(X,y[tr],sample_weight=w)
        pred[te]=model.predict(Z)
    if not np.isfinite(pred).all():raise RuntimeError('OOF source row missed')
    return (y-pred)**2

def run(panel,source,placebo,models,out):
    d,_,audit=join_direction(panel,source,placebo)
    d['country_method']=d.country.fillna('__UNKNOWN__').astype(str)+'::'+d.method.fillna('__UNKNOWN__').astype(str)
    saved=json.loads(models.read_text())
    results={};fail=[]
    for label,(target,cols) in {'HNS':('primary_record_mean',BASE+HNS),
      'EQR':('EQR_original_site_equal_basin_period_mean',BASE),
      'plant_S50':('plant_S50',BASE)}.items():
        errors={k:sk_oof(d,target,cols+features) for k,features in RECIPES.items()}
        h6=d.h6.astype(str).to_numpy();groups=sorted(set(h6))
        base=np.array([errors['static_W'][h6==q].mean() for q in groups]);
        compare={}
        for model in list(RECIPES)[1:]:
            e=np.array([errors[model][h6==q].mean() for q in groups]);gain=float(100*(base.mean()-e.mean())/base.mean())
            original=saved['outcomes'][label]['scores'][model]['gain'];diff=abs(gain-original)
            if diff>1e-7:fail.append([label,model,diff])
            compare[model]={'independent_sklearn_gain':gain,'original_manual_gain':original,'absolute_gain_diff':diff}
        results[label]=compare
    report={'schema':'GESUI_454_H7_DIRECTIONAL_PLACEBO_INDEPENDENT_NUMERICAL_QA_V1',
     'source_rows':len(d),'original_H6':d.h6.nunique(),'n_independent_specification_comparisons':sum(map(len,results.values())),
     'max_gain_abs_diff':max(x['absolute_gain_diff'] for s in results.values() for x in s.values()),
     'pass':not fail,'failed_specifications':fail,'original_upstream_source_H7_periods':audit['H7_period_upstream_positive'],
     'not_independent_ecology':True,'individual_butterfly_abundance_not_observed':True,'results':results}
    out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2,sort_keys=True))
    print(json.dumps({k:report[k] for k in ('pass','n_independent_specification_comparisons','max_gain_abs_diff')},sort_keys=True))
    if fail:raise SystemExit('Independent numerical source-model mismatch')
    return report
if __name__=='__main__':
    a=argparse.ArgumentParser()
    for k in ('panel','source','placebo','models','out'):a.add_argument('--'+k,type=pathlib.Path,required=True)
    d=a.parse_args();run(d.panel,d.source,d.placebo,d.models,d.out)
