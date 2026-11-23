#!/usr/bin/env python3
"""Source-matched FR381 original RZ2018 cross-sectional vs 2022 EQR.

FOUR pre-fixed posthoc spatial-vs-lagged falsifiers from same exact 381
EEA native site IDs and original 2018 geography. No sewage exposure causal
claim; 2018 land is not 2022 change in land cover or actual leaf nitrogen.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import gesui_EEA2018_RZ_to2022_france_EQR_future_status_20261123 as original

NUM_BASE=["lat","lon","LogQ","LogArea"]
NUM_L=NUM_BASE+["NaturalFar","UrbanFar","CropFar"]
NUM_D=NUM_L+["NaturalNearMinusFar"]
MODELS={"B0_2018_geo_hydro":NUM_BASE,
        "L0_2018_geo_hydro_plus_land":NUM_L,
        "D0_2018_nearouter_contrast":NUM_D}
ORIG_FR="FR"
SEED=20261123;NBOOT=2000
def source_oof(df,target):
 if len(df)!=381:raise RuntimeError("Frozen same-site France cohort changed")
 y=df[target].to_numpy(float)
 fold=df.pre_2018_source_only_spatial_fold.to_numpy(int)
 blocks=df.region_3deg_block_2018.to_numpy(str)
 if len(set(blocks))!=14 or set(fold)!={1,2,3,4,5}:raise RuntimeError("Original source region/folds changed")
 for reg in set(blocks):
  if len(set(fold[blocks==reg]))!=1:raise RuntimeError("Original source blocks split over folds")
 losses={}
 fold_MSE={}
 for name,cols in MODELS.items():
  X=df[cols].to_numpy(float)
  pred=np.full(len(y),np.nan)
  fold_MSE[name]={}
  for f in range(1,6):
   tr=fold!=f;te=fold==f
   pipe=Pipeline([("scale",StandardScaler()),("ridge",Ridge(alpha=10,solver="svd"))])
   pipe.fit(X[tr],y[tr])
   pred[te]=pipe.predict(X[te])
   fold_MSE[name][str(f)]={"test_sites":int(np.sum(te)),
     "original_site_MSE":float(np.mean((y[te]-pred[te])**2))}
  if not np.isfinite(pred).all():raise RuntimeError("2018 spatial source original fold missing")
  losses[name]=(y-pred)**2
 return losses,fold_MSE
def contrasts(df,losses):
 blocks=df.region_3deg_block_2018.astype(str).to_numpy()
 regions=np.array(sorted(set(blocks)))
 if len(regions)!=14:raise RuntimeError("Original 2018 block count mismatch")
 same={k:np.array([np.mean(v[blocks==r]) for r in regions]) for k,v in losses.items()}
 rng=np.random.default_rng(SEED)
 ix=rng.integers(0,len(regions),size=(NBOOT,len(regions)))
 a={}
 for key,base,full in (
  ("L0_after_original_hydrogeo","B0_2018_geo_hydro","L0_2018_geo_hydro_plus_land"),
  ("D0_after_L0","L0_2018_geo_hydro_plus_land","D0_2018_nearouter_contrast")):
  b=losses[base];f=losses[full];bb=same[base];ff=same[full]
  boot=100*(bb[ix].mean(axis=1)-ff[ix].mean(axis=1))/bb[ix].mean(axis=1)
  a[key]={"heldout_site_equal_baseline_MSE":float(b.mean()),
    "heldout_site_equal_added_MSE":float(f.mean()),
    "original_5fold_site_equal_MSE_gain_pct":float(100*(b.mean()-f.mean())/b.mean()),
    "2018_3deg_block_equal_MSE_gain_pct":float(100*(bb.mean()-ff.mean())/bb.mean()),
    "fixed_OOF_block_descriptive_bootstrap95_pct":[float(v) for v in np.quantile(boot,[.025,.975])],
    "14_region_blocks":len(regions),
    "source_regions_where_added_model_improves":int(np.sum(ff<bb)),
    "posthoc_falsifier_not_causality":True}
 return a
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq22","cohort_source","original_lag_numeric","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 orig=json.loads(a.original_lag_numeric.read_text())
 if orig["schema"]!="GESUI_EEA2018_RZ_OFFICIAL_TO_EEA2022_FRENCH_BENTHOS_LAGGED_GEOFORECAST_V1":
  raise RuntimeError("Source signed primary lag result missing")
 df,qa,gate=original.frame(a)
 if len(df)!=381:raise RuntimeError("Original exact source-purged France381 changed")
 if gate["source_frozen_acceptance"]!="FRANCE_ONLY":raise RuntimeError("Source go status drift")
 if abs(float(df.EQR2022.mean())-orig["source_real_2022_benthic_EQR_mean"])>1e-12:
  raise RuntimeError("Original EQR2022 response revised")
 if abs(float(df.numeric_EQR_2018.mean())-orig["source_lagged_2018_benthic_EQR_mean"])>1e-12:
  raise RuntimeError("Original EQR2018 response revised")
 res={}
 for target in ("numeric_EQR_2018","EQR2022"):
  loss,folds=source_oof(df,target)
  res[target]={"source_true_year":2018 if target=="numeric_EQR_2018" else 2022,
    "true_original_source_sites":len(df),
    "source_fold_model_MSE":folds,
    "fixed_posthoc_source_covariate_gains":contrasts(df,loss)}
 ans={"schema":"GESUI_EEA_FR381_RZ2018_same_site_EQR2018_v_EQR2022_spatial_v_future_source_diagnostic_v1",
  "source_years":[2018,2022],
  "source_true_france_sites":381,
  "original_only_source_5geofolds":True,"source_2018_3deg_region_count":14,
  "original_primary_2022_after_2018_EQR_site_MSE_gain":orig["original_2022_future_EQR_forecast_incremental_source"]["L_after_lagged_EQR"],
  "source_fixed_alpha":10.,
  "source_same_2018_land_context_no_2022_land_change_measured":True,
  "no_potw_time_change_or_aquatic_terrestrial_material_transfer":True,
  "posthoc_not_confirmatory_spatial_results":True,
  "source_lagged_original_numerical_SHA256":hashlib.sha256(a.original_lag_numeric.read_bytes()).hexdigest(),
  "outcomes":res}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(ans,sort_keys=True,indent=2)+"\n")
 print("FRANCE381_SAME_COHORT_SPATIAL_V_FUTURE_LAND_FALSIFIER",json.dumps({
  k:{name:v["fixed_posthoc_source_covariate_gains"][name] for name in (
    "L0_after_original_hydrogeo","D0_after_L0")}
  for k,v in res.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
