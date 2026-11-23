#!/usr/bin/env python3
"""Independent manual weighted Ridge10 5fold audit of France381 SAME SITES:
2018 same-year and 2022 no-lag EQR baseline vs 2018 RZ land and bank contrast.

No original sklearn Pipeline/StandardScaler/Ridge code called for model fit.
Recalculate all 4 fixed source-frozen contrasts and 14-block error bootstraps.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np
import gesui_EEA2018_RZ_to2022_france_EQR_future_status_20261123 as prior
import gesui_france381_RZ2018_same_cohort_2018_vs2022_spatial_status_20261123 as original

MODELS=original.MODELS
def near(a,b,label):
 if not np.allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),atol=2e-10,rtol=2e-10):
  raise RuntimeError("Independent EEA France source response/model mismatch "+label+" "+str((a,b)))
def predictions(frame,target):
 original_y=frame[target].to_numpy(float)
 fold=frame.pre_2018_source_only_spatial_fold.to_numpy(int)
 o={}
 for model,fields in MODELS.items():
  matrix=frame[fields].to_numpy(float)
  pred=np.full(len(frame),np.nan)
  for g in range(1,6):
   train=fold!=g;test=fold==g
   vals=matrix[train];predX=matrix[test]
   mean=vals.mean(axis=0)
   std=vals.std(axis=0,ddof=0)
   std=np.where(std<1e-12,1.,std)
   Z=(vals-mean)/std
   V=(predX-mean)/std
   ytrain=original_y[train]
   ym=ytrain.mean()
   beta=np.linalg.solve(Z.T@Z+10.*np.eye(Z.shape[1]),Z.T@(ytrain-ym))
   pred[test]=ym+V@beta
  if not np.isfinite(pred).all():raise RuntimeError("Manual original 2018 fold predictions failed")
  o[model]=(original_y-pred)**2
 return o
def contrasts(frame,loss):
 regional=frame.region_3deg_block_2018.astype(str).to_numpy()
 codes=np.unique(regional)
 if len(codes)!=14:raise RuntimeError("Original France source geographic blocks altered")
 grouped={name:np.array([np.mean(e[regional==r]) for r in codes]) for name,e in loss.items()}
 draws=np.random.default_rng(20261123).integers(0,len(codes),size=(2000,len(codes)))
 out={}
 for name,base,more in (("L0_after_original_hydrogeo","B0_2018_geo_hydro","L0_2018_geo_hydro_plus_land"),
                         ("D0_after_L0","L0_2018_geo_hydro_plus_land","D0_2018_nearouter_contrast")):
  a=loss[base];b=loss[more];ra=grouped[base];rb=grouped[more]
  boots=100*(ra[draws].mean(axis=1)-rb[draws].mean(axis=1))/ra[draws].mean(axis=1)
  out[name]={"baseline_MSE":float(a.mean()),"full_MSE":float(b.mean()),
    "site_gain":float(100*(a.mean()-b.mean())/a.mean()),
    "region_gain":float(100*(ra.mean()-rb.mean())/ra.mean()),
    "region_boot95":[float(v) for v in np.quantile(boots,[.025,.975])],
    "regional_14_improvement_n":int(np.sum(rb<ra))}
 return out
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq22","cohort_source","original_lag_numeric","original_spatial_numeric","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 orig=json.loads(a.original_spatial_numeric.read_text())
 base=json.loads(a.original_lag_numeric.read_text())
 if orig["schema"]!="GESUI_EEA_FR381_RZ2018_same_site_EQR2018_v_EQR2022_spatial_v_future_source_diagnostic_v1":
  raise RuntimeError("Original 381 France source-locked spatial diagnostic missing")
 if orig["source_lagged_original_numerical_SHA256"]!=hashlib.sha256(a.original_lag_numeric.read_bytes()).hexdigest():
  raise RuntimeError("Original primary France 2022 source biological result SHA changed")
 df,qa,gate=prior.frame(a)
 if len(df)!=381 or df.site_key.duplicated().any():raise RuntimeError("France381 source ID changed")
 assert abs(df.EQR2022.mean()-base["source_real_2022_benthic_EQR_mean"])<1e-12
 allres={}
 for target in ("numeric_EQR_2018","EQR2022"):
  independent=contrasts(df,predictions(df,target))
  saved=orig["outcomes"][target]["fixed_posthoc_source_covariate_gains"]
  for name,d in independent.items():
   expected=saved[name]
   for key,k in (("baseline_MSE","heldout_site_equal_baseline_MSE"),
                 ("full_MSE","heldout_site_equal_added_MSE"),
                 ("site_gain","original_5fold_site_equal_MSE_gain_pct"),
                 ("region_gain","2018_3deg_block_equal_MSE_gain_pct"),
                 ("region_boot95","fixed_OOF_block_descriptive_bootstrap95_pct"),
                 ("regional_14_improvement_n","source_regions_where_added_model_improves")):
    near(d[key],expected[k],target+"/"+name+"/"+key)
  allres[target]=independent
 result={"schema":"GESUI_FR381_2018_vs2022_RZ_land_spatial_NOlag_independent_NumPyQA_v1",
  "original_2018_to2022_FR_same_ID_method_reach_sites":381,
  "original_geoblock_count":14,
  "source_2018RZ_2018EQ_R_and_2022EQ_R_recap_verified":True,
  "original_signs_including_nulls_reproduced":True,
  "manual_closed_form_Ridge10_no_original_sklearn_prediction":True,
  "all_four_original_source_frozen_landcover_results_and_14block_bootstraps_reproduced":True,
  "original_signed_Spatial_2018_2022_SHA256":hashlib.sha256(a.original_spatial_numeric.read_bytes()).hexdigest(),
  "no_WWTP_causal_source_or_butterfly_plant_biomass":True,
  "results":allres}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n")
 print("FR381_SAMECOHORT_INDEPENDENT_2018_2022_SPATIAL_EQR_QA",json.dumps({
   outcome:{name:{"site_gain_pct":x["site_gain"],"region_gain_pct":x["region_gain"],
                  "region_boot95":x["region_boot95"]}
    for name,x in res.items()} for outcome,res in allres.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
