#!/usr/bin/env python3
"""Whole-ecoregion OOD test of original river width beyond modeled flow and area.

Source/posthoc protocol:
docs/PROTOCOL_20261122_NRSA2324_LEAVE_ECO9_OUT_WIDTH_TRANSPORT.md
No sewage source dose, habitat restoration causal inference or ecological flux.
"""
import argparse,hashlib,json,math,pathlib
import numpy as np,pandas as pd,statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge
import gesui_nrsa2324_tp_canopy_benthos_interaction_20261122 as original
import gesui_nrsa2324_flow_area_width_stress_20261122 as stress

TEST_MODELS={
 "without_width":stress.CORE_MODELS["B_NC_FA"],
 "with_width":stress.CORE_MODELS["B_NC_FA_W"]
}
EXPECTED_HUC8_BASE={
 ("BERW","MMI_BENT"):(264.71084481132004,254.32878760787608),
 ("BERW","EPT_fraction"):(0.041372101876049616,0.040613377243087535),
 ("BETB","MMI_BENT"):(310.4172641143832,289.4003970467362),
 ("BETB","EPT_fraction"):(0.03824591904263095,0.03676116387349702)
}
def fit_oof(d,y):
 eco_labels=sorted(d.eco9.dropna().unique())
 if len(eco_labels)!=9:raise ValueError("Unexpected EPA nine ecological reporting regions "+str(eco_labels))
 p={k:np.full(len(d),np.nan) for k in TEST_MODELS}
 fold_report={}
 for eco in eco_labels:
  test=d.eco9.eq(eco).to_numpy()
  overlap_huc8=set(d.loc[test,"huc8"])
  train=(~test)&(~d.huc8.isin(overlap_huc8).to_numpy())
  ntest=int(test.sum());ntrain=int(train.sum())
  hs_test=int(d.loc[test,"huc8"].nunique());hs_train=int(d.loc[train,"huc8"].nunique())
  discarded=int((~test).sum()-ntrain)
  fold_report[eco]={"test_UID":ntest,"test_HUC8":hs_test,
      "train_UID":ntrain,"train_HUC8":hs_train,
      "train_UID_purged_shared_test_HUC8":discarded,
      "status":"SOURCE_SUPPORT"}
  if min(ntest,hs_test,ntrain,hs_train)<0:raise AssertionError("Invalid fold")
  if ntest<10 or hs_test<5 or ntrain<250 or hs_train<100:
   fold_report[eco]["status"]="SOURCE_STOP"
   continue
  for k,cols in TEST_MODELS.items():
   transform=ColumnTransformer([
    ("num",StandardScaler(),cols),
    ("cat",OneHotEncoder(handle_unknown="ignore",sparse_output=False),stress.CAT)
   ],sparse_threshold=0)
   xt=transform.fit_transform(d.loc[train,cols+stress.CAT])
   xv=transform.transform(d.loc[test,cols+stress.CAT])
   model=Ridge(alpha=10,solver="svd").fit(xt,y[train])
   p[k][test]=model.predict(xv)
  fold_report[eco]["status"]="FIT_OK"
 if not all(x["status"]=="FIT_OK" for x in fold_report.values()):
  return {"status":"REGIONAL_SUPPORT_STOP",
    "folds":fold_report,"reason":"At least one ecological region fails fixed source eligibility"}
 if not all(np.isfinite(v).all() for v in p.values()):
  return {"status":"OOF_PREDICTION_STOP","folds":fold_report}
 e0=(y-p["without_width"])**2
 e1=(y-p["with_width"])**2
 region0=pd.Series(e0,index=d.eco9.to_numpy()).groupby(level=0).mean().sort_index()
 region1=pd.Series(e1,index=d.eco9.to_numpy()).groupby(level=0).mean().reindex(region0.index)
 R0=region0.to_numpy();R1=region1.to_numpy()
 gain=100*(R0-R1)/R0
 rng=np.random.default_rng(20261122)
 draw=rng.integers(0,len(R0),size=(2000,len(R0)))
 bg=100*(R0[draw].mean(axis=1)-R1[draw].mean(axis=1))/R0[draw].mean(axis=1)
 return {"status":"COMPLETE_NINE_REGIONS",
   "heldout_ecoregions":9,
   "n_regions_positive":int(np.count_nonzero(gain>0)),
   "n_regions_negative":int(np.count_nonzero(gain<0)),
   "region_equal_MSE_reference":float(R0.mean()),
   "region_equal_MSE_plus_width":float(R1.mean()),
   "region_equal_MSE_gain_pct":float(100*(R0.mean()-R1.mean())/R0.mean()),
   "region_bootstrap_fixedOOF_95_pct":[float(q) for q in np.quantile(bg,[.025,.975])],
   "visit_weighted_MSE_gain_pct":float(100*(e0.mean()-e1.mean())/e0.mean()),
   "heldout_ecoregion_diagnostics":{
     k:{**fold_report[k],
        "baseline_MSE":float(region0[k]),"with_width_MSE":float(region1[k]),
        "MSE_gain_pct":float(100*(region0[k]-region1[k])/region0[k])}
     for k in region0.index},
   "zero_shared_HUC8_between_each_fold_train_and_test":True}
def descriptive_width(d,y):
 cols=TEST_MODELS["with_width"]
 z=d[cols].astype(float)
 sd=z.std(ddof=0).replace(0,1)
 z=(z-z.mean())/sd
 category=pd.get_dummies(d[stress.CAT],drop_first=True,dtype=float)
 X=pd.concat((z.reset_index(drop=True),category.reset_index(drop=True)),axis=1)
 X=sm.add_constant(X,has_constant="add")
 yy=(y-y.mean())/y.std(ddof=0)
 try:
  fit=sm.OLS(yy,X.to_numpy(float)).fit(cov_type="cluster",
   cov_kwds={"groups":d.huc8.to_numpy(),"use_correction":True})
  i=list(X.columns).index("log_width");ci=fit.conf_int(alpha=.05)[i]
  return {"standardized_width_beta":float(fit.params[i]),
   "HUC8_cluster_descriptive95":[float(ci[0]),float(ci[1])],
   "p_uncorrected":float(fit.pvalues[i]),
   "no_causal_effect_claim":True}
 except Exception as ex:
  return {"status":"DESCRIPTIVE_COEFFICIENT_FAILED","reason":str(ex)[:300]}
def main():
 ap=argparse.ArgumentParser()
 for k in original.SHA:ap.add_argument("--"+k,type=pathlib.Path,required=True)
 ap.add_argument("--out",type=pathlib.Path,required=True)
 a=ap.parse_args();source={}
 for k,src_hash in original.SHA.items():
  p=getattr(a,k)
  if original.sha(p)!=src_hash:raise ValueError("Original EPA source digest mismatch "+k)
  source[k]=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
 panel,checks,n=original.panel(source)
 if n!=2139 or sum(checks.values()):raise ValueError("Original source UID QA failed")
 panel=stress.numeric_embedded(stress.normalize_hydro(panel,source["siteinfo"]),source["physicalhabitat"])
 results={}
 for method in ("BERW","BETB"):
  results[method]={}
  for outcome in ("MMI_BENT","EPT_fraction"):
   original_d=original.cohort(panel,method,outcome)
   d,support=stress.covariate_support(original_d,method)
   if len(d)<250 or d.huc8.nunique()<100:
    results[method][outcome]={"status":"SOURCE_STOP","support":support};continue
   y=d[outcome].to_numpy(float)
   base,_,_=stress.oof(d,y,TEST_MODELS["without_width"],stress.CAT)
   full,_,_=stress.oof(d,y,TEST_MODELS["with_width"],stress.CAT)
   ref0,ref1=EXPECTED_HUC8_BASE[(method,outcome)]
   if not np.isclose(base.mean(),ref0,rtol=1e-8,atol=1e-9) or not np.isclose(full.mean(),ref1,rtol=1e-8,atol=1e-9):
    raise RuntimeError("Original 5fold HUC8 source-locked reference not numerically reproduced")
   report={"status":"RUN_FINISHED","cohort":{"UID":len(d),"physical_sites":int(d.site.nunique()),
        "HUC8":int(d.huc8.nunique()),"eco9":int(d.eco9.nunique())},
    "huc8_previous_ref_reproduced":True,
    "width_descriptive":descriptive_width(d,y),
    "region_extrapolation":fit_oof(d,y)}
   results[method][outcome]=report
   print("NRSA_ECO9_OOD_WIDTH_RESULT",method,outcome,json.dumps(report,sort_keys=True),flush=True)
 out={"schema":"gesui_NRSA2324_full_ecoregion_holdout_flowarea_width_extrapolation_v1",
  "EPA_original_SHA":original.SHA,"origin_5table_UID":n,
  "preoutcome_method_and_source_preserved":True,
  "posthoc":True,"no_sewage_dose":True,"results":results}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
 print("NRSA_ECO9_OOD_WIDTH_SUMMARY",json.dumps({m:{response:{
    "n":v["cohort"]["UID"],
    "coef":v["width_descriptive"].get("standardized_width_beta"),
    "extra_region_gain":v["region_extrapolation"].get("region_equal_MSE_gain_pct"),
    "regions_positive":v["region_extrapolation"].get("n_regions_positive")}
  for response,v in g.items()} for m,g in results.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
