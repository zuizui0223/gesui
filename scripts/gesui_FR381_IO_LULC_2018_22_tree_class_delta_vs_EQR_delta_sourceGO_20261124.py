#!/usr/bin/env python3
"""Only after strict source-quality GO: actual France EEA2018→2022 benthic
EQR change vs genuine annual 2018–2022 IO V2 land-cover CLASS change.

Scientific protocol fixed before 2018–2022 land pixel response inspection;
five original geographical folds, Ridge alpha10, 14 region-equal errors.
NOT canopy density, plant physiology, WWTP dose, or causal material transfer.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib,sys
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
sys.path.insert(0,"scripts")
import gesui_EEA2018_RZ_to2022_france_EQR_future_status_20261123 as original

SIG="GESUI_REAL_FR381_2018_2022_IO_10M_CLASS_SOURCE_YEARPAIRED_PIXELS_V1"
B=["numeric_EQR_2018","lat","lon","LogQ","LogArea","NaturalFar","UrbanFar","CropFar"]
TARGET="EQR_delta"
GEO="region_3deg_block_2018"
FOLD="pre_2018_source_only_spatial_fold"
NBOOT=2000
SEED=20261124
def source_frame(a):
 qa=json.loads(a.land_source.read_text())
 if qa.get("schema")!=SIG:
  raise RuntimeError("Original 381 source-only annual class/quality source missing")
 if qa["source_strict_high_quality_go"] is not True or qa["source_status"]!="SOURCE_GO_ORIGINAL_LULC2018_2022_SOURCE_QUALITY_PIXEL_PASS_ECOLOGY_NEXT":
  out={"schema":"GESUI_FRANCE381_2018_2022_LULC_TO_DELTA_EQR_SOURCE_NO_GO",
       "status":"SOURCE_STOP_PIXEL_SUPPORT_NO_ECOLOGY_FITTED",
       "original_lulc_eligible_sites":qa.get("source_france_original_site_quality_eligible"),
       "no_original_aquatic_EQR_values_read":True}
  return None,qa,out
 if (qa["source_original_EEA_FR381_sites_before_outcome"]!=381 or
     qa["source_original_2018_fixed_fold_valid_site_counts"] is None):
  raise RuntimeError("Source 381 original cohort not locked")
 signed=json.loads(a.old_original_lagged.read_text())
 if signed.get("schema")!="GESUI_EEA2018_RZ_OFFICIAL_TO_EEA2022_FRENCH_BENTHOS_LAGGED_GEOFORECAST_V1":
  raise RuntimeError("Earlier 381 real EEA source/folds verification missing")
 # Only NOW original EQR numeric gets opened, after independent source-only GO.
 frame,src,gt=original.frame(a)
 if len(frame)!=381 or len(qa["source_full_census_original_site_quality_records"])!=381:
  raise RuntimeError("Original 381 same-site source count unexpectedly changed")
 per={z["site_hash_sha256"]:z for z in qa["source_full_census_original_site_quality_records"]}
 if len(per)!=381:raise RuntimeError("Original 381 source hash duplicated")
 parts=[]
 for _,x in frame.iterrows():
  site_hash=hashlib.sha256(x.site_key.encode()).hexdigest()
  if site_hash not in per:raise RuntimeError("Source no original France EEA site hash")
  source=per[site_hash]
  if not source["source_quality_pass"]:continue
  near=source["real_2018_2022_category_fraction_source"]["near"]["class"]
  far=source["real_2018_2022_category_fraction_source"]["far"]["class"]
  parts.append({"site_key":x.site_key,"near_tree_delta":near["trees_delta_frac"],
     "far_tree_delta":far["trees_delta_frac"],"nearminusfar_tree_delta":near["trees_delta_frac"]-far["trees_delta_frac"]})
 if not parts:return None,qa,{"status":"SOURCE_STOP_NO_ACTUAL_QUALIFIED_SITES","no_EEA_model_fitted":True}
 extra=pd.DataFrame(parts)
 d=frame.merge(extra,on="site_key",how="inner",validate="one_to_one")
 if (len(d)!=qa["source_france_original_site_quality_eligible"] or
     set(d[FOLD].astype(int))!={1,2,3,4,5} or
     d[GEO].nunique()<10):
  raise RuntimeError("Source image-based eligibility join changed after biological GO")
 regionfold=d.groupby(GEO)[FOLD].nunique()
 if (regionfold>1).any():raise RuntimeError("Spatial block leaked across fixed original geofolds")
 if not np.isfinite(d[B+["near_tree_delta","far_tree_delta","nearminusfar_tree_delta",TARGET]].to_numpy(float)).all():
  raise RuntimeError("True 2018/22 source data missing")
 return d,qa,None
def cv(d):
 sourcefold=d[FOLD].to_numpy(int)
 y=d[TARGET].to_numpy(float)
 models={
  "B_static_2018_status_geo_RZ":B,
  "B_plus_delta_tree_nearminusfar":B+["nearminusfar_tree_delta"],
  "B_plus_delta_tree_far_only":B+["far_tree_delta"]}
 errors={};byfold={}
 for key,cols in models.items():
  preds=np.full(len(d),np.nan)
  byfold[key]={}
  for fold in range(1,6):
   tr=sourcefold!=fold;te=sourcefold==fold
   xtrain=d.loc[tr,cols].to_numpy(float)
   xvalid=d.loc[te,cols].to_numpy(float)
   z=StandardScaler()
   xt=z.fit_transform(xtrain)
   xv=z.transform(xvalid)
   m=Ridge(alpha=10,solver="svd")
   m.fit(xt,y[tr])
   preds[te]=m.predict(xv)
   byfold[key][str(fold)]={"n":int(te.sum()),"MSE":float(np.mean((y[te]-preds[te])**2))}
  if not np.isfinite(preds).all():raise RuntimeError("Original 5fold source outcome predictions incomplete")
  errors[key]=(y-preds)**2
 return errors,byfold
def improvements(d,e):
 regions=d[GEO].to_numpy(str)
 keys=sorted(set(regions))
 if len(keys)<10:raise RuntimeError("Original source region source gate disallowed")
 block={name:np.array([np.mean(a[regions==k]) for k in keys],float) for name,a in e.items()}
 base="B_static_2018_status_geo_RZ"
 rng=np.random.default_rng(SEED)
 idx=rng.integers(0,len(keys),size=(NBOOT,len(keys)))
 out={}
 for name in ("B_plus_delta_tree_nearminusfar","B_plus_delta_tree_far_only"):
  eb=e[base];em=e[name]
  bb=block[base];bm=block[name]
  gains=100*(bb[idx].mean(1)-bm[idx].mean(1))/bb[idx].mean(1)
  out[name+"_after_original_base"]={
   "source_real_sites":len(d),"source_regional_blocks":len(keys),
   "source_original_2018_geo_folds":5,
   "heldout_original_baseline_site_equal_MSE":float(eb.mean()),
   "heldout_added_site_equal_MSE":float(em.mean()),
   "site_equal_MSE_improvement_pct":float(100*(eb.mean()-em.mean())/eb.mean()),
   "source_3degree_region_equal_MSE_improvement_pct":float(100*(bb.mean()-bm.mean())/bb.mean()),
   "fixed_outoffold_2000_georegion_error_boot_descriptive95_pct":[float(v) for v in np.quantile(gains,[.025,.975])],
   "n_regions_with_improved_errors":int((bm<bb).sum()),
   "same_source_2018_2022_status_not_causal_or_WWTP":True}
 return out
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq22","cohort_source","old_original_lagged","land_source","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 d,source,stop=source_frame(a)
 if stop is not None:report=stop
 else:
  error,folds=cv(d)
  out=improvements(d,error)
  x=d["nearminusfar_tree_delta"].to_numpy(float)
  report={"schema":"GESUI_FRANCE_RIPARIAN_IO_2018_22_TREE_CLASS_DELTA_TO_EEA_BENTHIC_EQR_DELTA_ACTUAL_v1",
   "status":"POSTHOC_SOURCE_LOCKED_TRUE_2018_22_LULC_TREE_CLASS_CHANGE_TO_DELTA_EQR",
   "real_france_source_sites":len(d),
   "original_preoutcome_categorical_pixel_source_SHA256":hashlib.sha256(a.land_source.read_bytes()).hexdigest(),
   "source_prespecified_only_new_tree_delta_predictor":"nearminusfar_tree_delta",
   "2018_2022_IO_V2_original_tree_relative_land_fraction_change_mean":float(d.near_tree_delta.mean()),
   "2018_2022_IO_V2_original_outer_tree_fraction_change_mean":float(d.far_tree_delta.mean()),
   "nearminusfar_tree_delta_stats":{"mean":float(x.mean()),"min":float(x.min()),
      "q10_q90":[float(v) for v in np.quantile(x,[.1,.9])],"max":float(x.max())},
   "actual_2018_2022_3deg_spatial_holdout_model_alpha":10,
   "actual_2018_2022_EQR_delta_response_mean":float(d.EQR_delta.mean()),
   "same_site_original_EQR_2018_vs2022_not_independent":True,
   "original_geo_fold_site_counts":d[FOLD].value_counts().sort_index().to_dict(),
   "fixed_predictive_comparisons":out,
   "original_fold_MSE_by_source":folds,
   "source_time_interval_does_not_establish_temporal_lead":True,
   "no_actual_canopy_density_biomass_or_2018_land_LCLU_change":True,
   "not_independent_WWTP_material_source_or_insect_emergence":True,
   "posthoc_ecology_not_confirmatory":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FR381_IO_2018_2022_TREE_CLASS_DELTA_EQR_DELTA_SOURCELOCKED",json.dumps(
  {"status":report["status"],"source_n":report.get("real_france_source_sites"),
   "results":report.get("fixed_predictive_comparisons"),
   "source_quality_stop":stop},sort_keys=True),flush=True)
if __name__=="__main__":main()
