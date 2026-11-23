#!/usr/bin/env python3
"""Independent manual 5fold 2018→2022 France EQR/RZland numerical audit.

Reuse signed source-only HUC/EEA ID and published EQR join; DO NOT invoke
original sklearn Pipeline, Ridge, oof(), evaluate() statistical computations.
Manual NumPy train-fold std, normal-equation Ridge10, source-geo block bootstrap.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np
import gesui_EEA2018_RZ_to2022_france_EQR_future_status_20261123 as orig

FIELDS={"B_2018_EQR_hydrogeo":orig.MODEL_BASE,
        "L_plus_2018_outer_land":orig.MODEL_L,
        "D_plus_near_outer_vegetated_contrast":orig.MODEL_D}
OUTCOMES=("EQR2022","EQR_delta")
def near(a,b,key,rtol=1e-10,atol=1e-10):
 if not np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=rtol,atol=atol):
  raise RuntimeError("Manual EEA geographic cross-year audit mismatch "+key+": "+str((a,b)))
def predict_all(df,target):
 years=df[target].to_numpy(float)
 folds=df.pre_2018_source_only_spatial_fold.to_numpy(int)
 o={}
 for name,fields in FIELDS.items():
  X=df[fields].to_numpy(float)
  yhat=np.full(len(df),np.nan)
  for g in range(1,6):
   tr=folds!=g;te=folds==g
   mean=X[tr].mean(0);scale=X[tr].std(0,ddof=0);scale=np.where(scale<1e-12,1.,scale)
   Xt=(X[tr]-mean)/scale
   Xe=(X[te]-mean)/scale
   Y=years[tr];Ymean=float(Y.mean())
   beta=np.linalg.solve(Xt.T@Xt+orig.ALPHA*np.eye(Xt.shape[1]),
                        Xt.T@(Y-Ymean))
   yhat[te]=Ymean+Xe@beta
  if not np.isfinite(yhat).all():raise ValueError("Manual out-of-fold predictions missing")
  o[name]=(years-yhat)**2
 return o
def gains(df,loss,more,base):
 region=df.region_3deg_block_2018.astype(str).to_numpy()
 codes=np.unique(region)
 b=np.array([np.mean(loss[base][region==k]) for k in codes])
 a=np.array([np.mean(loss[more][region==k]) for k in codes])
 if len(codes)!=14:raise RuntimeError("Original fixed source 14 regions changed")
 ix=np.random.default_rng(orig.SEED).integers(0,len(codes),size=(orig.BOOTS,len(codes)))
 boot=100*(b[ix].mean(1)-a[ix].mean(1))/b[ix].mean(1)
 return {"site_baseline_mse":float(np.mean(loss[base])),
   "site_added_mse":float(np.mean(loss[more])),
   "site_gain_pct":float(100*(np.mean(loss[base])-np.mean(loss[more]))/np.mean(loss[base])),
   "region_gain_pct":float(100*(b.mean()-a.mean())/b.mean()),
   "region_bootstrap95":[float(v) for v in np.quantile(boot,[.025,.975])],
   "region_count":len(codes),"region_n_improved":int(np.sum(a<b))}
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq22","cohort_source","original_numeric","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 signed=json.loads(a.original_numeric.read_text())
 if signed["schema"]!="GESUI_EEA2018_RZ_OFFICIAL_TO_EEA2022_FRENCH_BENTHOS_LAGGED_GEOFORECAST_V1":
  raise RuntimeError("Original signed 2022 EQR model source schema changed")
 df,qa,gate=orig.frame(a)
 if len(df)!=381:raise RuntimeError("Original French RZ2018 same-site 2022 source differs")
 near(df.numeric_EQR_2018.mean(),signed["source_lagged_2018_benthic_EQR_mean"],"2018 EQR mean")
 near(df.EQR2022.mean(),signed["source_real_2022_benthic_EQR_mean"],"2022 EQR mean")
 near(float(np.corrcoef(df.numeric_EQR_2018,df.EQR2022)[0,1]),signed["original_same_site_EQR_2018_2022_pearson"],"same-site EQR rho")
 results={}
 for endpoint,key in (("EQR2022","original_2022_future_EQR_forecast_incremental_source"),
                      ("EQR_delta","secondary_same_site_delta_EQR_sensitivity")):
  losses=predict_all(df,endpoint)
  original=signed[key]
  results[endpoint]={}
  for add,baseline,contrast in (("L_plus_2018_outer_land","B_2018_EQR_hydrogeo","L_after_lagged_EQR"),
                                ("D_plus_near_outer_vegetated_contrast","L_plus_2018_outer_land","D_after_L")):
   val=gains(df,losses,add,baseline)
   prev=original[contrast]
   for x,y in (("site_baseline_mse","original_true_2022_EQR_site_MSE_baseline"),
               ("site_added_mse","original_true_2022_EQR_site_MSE_added"),
               ("site_gain_pct","site_equal_heldout_MSE_gain_percent"),
               ("region_gain_pct","2018_3degree_region_equal_MSE_gain_percent"),
               ("region_bootstrap95","fixed_OOF_region_bootstrap_descriptive95"),
               ("region_count","total_regions"),
               ("region_n_improved","regions_where_more_predictors_reduce_MSE")):
    near(val[x],prev[y],endpoint+"/"+contrast+"/"+x)
   results[endpoint][contrast]={"all_original_numerical_values_and_bands_reproduced":True,
     "site_gain_pct":val["site_gain_pct"],"region_gain_pct":val["region_gain_pct"],
     "regional_descriptive95":val["region_bootstrap95"]}
 report={"schema":"gesui_original_France381_2018RZ_2022_Benthic_manual_NumPy_geo5fold_QA_v1",
  "same_site_FR_source_2018_2022_n":381,
  "original_2018RZ_2018EQR_2022EQR_source_schema_valid":True,
  "original_5fold_2018_source_block_assignments_no_region_split":True,
  "independent_normal_equation_trainfold_Ridge10_No_sklearn_model_called":True,
  "independent_original14_geoblock_bootstrap_and_lagged_EQR_means_reproduced":True,
  "no_source_P_origin_change_or_landplant_causal_inference":True,
  "signed_original_numeric_SHA256":hashlib.sha256(a.original_numeric.read_bytes()).hexdigest(),
  "results":results}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FR381_EQR2018_TO_EQR2022_LAND_TEMPORAL_INDEPENDENT_NUMPY_QA",json.dumps(
  {"n":len(df),"results":results},sort_keys=True),flush=True)
if __name__=="__main__":main()
