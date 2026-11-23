#!/usr/bin/env python3
"""Independent manual source-locked EEA France376 LULC-year-delta/EQR-delta QA.

Original model sklearn Pipeline NOT used. Recomputes genuine 5 geo-fold Ridge10
via NumPy normal equations with training-only column means and population SD,
and original 14-block error resampling. Replicates ALL null/positive results.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np
import gesui_FR381_IO_LULC_2018_22_tree_class_delta_vs_EQR_delta_sourceGO_20261124 as source

COLS={
 "B_static_2018_status_geo_RZ":source.B,
 "B_plus_delta_tree_nearminusfar":source.B+["nearminusfar_tree_delta"],
 "B_plus_delta_tree_far_only":source.B+["far_tree_delta"]}
TARGET=source.TARGET
REGION=source.GEO
FOLD=source.FOLD
def close(x,y,label):
 if not np.allclose(np.asarray(x,float),np.asarray(y,float),rtol=2e-9,atol=2e-10):
  raise RuntimeError("Independent source-locked France376 manual Ridge QA FAILED: "+label+" "+str((x,y)))
def actual_source_EEA(a):
 signed=json.loads(a.original_result.read_text())
 if signed.get("schema")!="GESUI_FRANCE_RIPARIAN_IO_2018_22_TREE_CLASS_DELTA_TO_EEA_BENTHIC_EQR_DELTA_ACTUAL_v1":
  raise RuntimeError("Original source-locked ecological result absent")
 d,geo,stop=source.source_frame(a)
 if stop is not None or d is None or len(d)!=376:
  raise RuntimeError("Original pre-outcome 376 source quality eligibility failed")
 if hashlib.sha256(a.land_source.read_bytes()).hexdigest()!=signed["original_preoutcome_categorical_pixel_source_SHA256"]:
  raise RuntimeError("Signed original 376-site source-only annual land pixel census digest differs")
 close(d.EQR_delta.mean(),signed["actual_2018_2022_EQR_delta_response_mean"],"actual EQR delta")
 close(d.near_tree_delta.mean(),signed["2018_2022_IO_V2_original_tree_relative_land_fraction_change_mean"],"near delta trees")
 close(d.far_tree_delta.mean(),signed["2018_2022_IO_V2_original_outer_tree_fraction_change_mean"],"far delta trees")
 return d,signed
def manual_fold_source(df):
 y=df[TARGET].to_numpy(float)
 folds=df[FOLD].to_numpy(int)
 errors={};oof_fold={}
 for name,fields in COLS.items():
  x=df[fields].to_numpy(float)
  pred=np.full(len(df),np.nan)
  fold_stats={}
  for i in range(1,6):
   tr=folds!=i;te=folds==i
   A=x[tr];T=x[te]
   mu=A.mean(axis=0)
   sigma=A.std(axis=0,ddof=0)
   sigma=np.where(sigma==0,1.,sigma)
   train=(A-mu)/sigma
   test=(T-mu)/sigma
   yt=y[tr]
   inter=float(yt.mean())
   # Unpenalized intercept and native Euclidean 10 ridge penalty:
   # manually solve normal equations independent of original sklearn.
   beta=np.linalg.solve(train.T@train+10.*np.eye(train.shape[1]),
                        train.T@(yt-inter))
   pred[te]=inter+test@beta
   fold_stats[str(i)]={"n":int(te.sum()),"MSE":float(np.mean((y[te]-pred[te])**2))}
  if not np.isfinite(pred).all():raise RuntimeError("Source 376 manual 5fold not covering every site")
  errors[name]=(y-pred)**2
  oof_fold[name]=fold_stats
 return errors,oof_fold
def manual_mse_contrasts(d,e,signed):
 regions=d[REGION].to_numpy(str)
 codes=sorted(set(regions))
 if len(codes)!=14:raise RuntimeError("Original 14 EEA independent geography blocks changed")
 grouped={m:np.array([np.mean(r[regions==k]) for k in codes],float) for m,r in e.items()}
 rboot=np.random.default_rng(20261124).integers(0,len(codes),size=(2000,len(codes)))
 base="B_static_2018_status_geo_RZ"
 out={}
 for model in ("B_plus_delta_tree_nearminusfar","B_plus_delta_tree_far_only"):
  z=e[base];f=e[model]
  b=grouped[base];a=grouped[model]
  boot=100*(b[rboot].mean(axis=1)-a[rboot].mean(axis=1))/b[rboot].mean(axis=1)
  stats={"source_original_2018_2022_eligible_sites":len(d),
         "site_MSE_baseline":float(z.mean()),
         "site_MSE_added":float(f.mean()),
         "site_gain_percent":float(100*(z.mean()-f.mean())/z.mean()),
         "region_equal_gain_percent":float(100*(b.mean()-a.mean())/b.mean()),
         "14region_bootstrap95_percent":[float(t) for t in np.quantile(boot,[.025,.975])],
         "n_regions_improved":int((a<b).sum())}
  prior=signed["fixed_predictive_comparisons"][model+"_after_original_base"]
  for p,q in (("site_MSE_baseline","heldout_original_baseline_site_equal_MSE"),
              ("site_MSE_added","heldout_added_site_equal_MSE"),
              ("site_gain_percent","site_equal_MSE_improvement_pct"),
              ("region_equal_gain_percent","source_3degree_region_equal_MSE_improvement_pct"),
              ("14region_bootstrap95_percent","fixed_outoffold_2000_georegion_error_boot_descriptive95_pct"),
              ("n_regions_improved","n_regions_with_improved_errors")):
   close(stats[p],prior[q],model+" "+p)
  out[model]={"status":"INDEPENDENT_MANUAL_NUMPY_SOURCE_QA_PASS",**stats}
 return out
def main():
 parser=argparse.ArgumentParser()
 for k in ("rz","eq22","cohort_source","old_original_lagged","land_source","original_result","out"):
  parser.add_argument("--"+k,type=pathlib.Path,required=True)
 args=parser.parse_args()
 df,signed=actual_source_EEA(args)
 losses,folds=manual_fold_source(df)
 for model,groups in signed["original_fold_MSE_by_source"].items():
  for fold,rec in groups.items():
   close(folds[model][fold]["n"],rec["n"],model+" fold "+fold+" sample")
   close(folds[model][fold]["MSE"],rec["MSE"],model+" fold "+fold+" mse")
 res=manual_mse_contrasts(df,losses,signed)
 report={
  "schema":"GESUI_REAL_FRANCE376_IO_V2_2018_22_DELTA_CLASS_TO_EEA_DELTA_EQR_MANUAL_NUMPY_QA_v1",
  "original_preoutcome_376site_land_quality_source_SHA256":hashlib.sha256(args.land_source.read_bytes()).hexdigest(),
  "original_primary_source_ecology_result_SHA256":hashlib.sha256(args.original_result.read_bytes()).hexdigest(),
  "independent_original_source_matched_sites":376,"original_pre_frozen_spatial_source_folds":5,
  "original_3degree_spatial_blocks":14,"original_fixed_training_only_Ridge_alpha":10,
  "manual_numerical_solver_does_not_use_original_sklearn_ridge":True,
  "all_original_fivefold_MSEs_and_14block_bootstrap_nulls_reproduced":True,
  "same_observational_field_source_not_new_independent_field_sampling":True,
  "source_ETP_WWTP_inference_not_established":True,
  "results":res}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FRANCE376_ORIGINAL_LULC_TREE_DELTA_EQR_MANUAL_NUMPY_QA",json.dumps({
  "n":len(df),"all_source_models_and_region_CI_reproduced":True,
  "comparisons":res},sort_keys=True),flush=True)
if __name__=="__main__":main()
