#!/usr/bin/env python3
"""Source-frozen actual France same monitoring station 2018 RZ habitat→2022 EQR.

This deliberately tests TEMPORAL predictive information beyond the 2018 actual
same-site EQR, rather than reusing the earlier 2018 spatial habitat-EQR signal.
All 381 French source sites and original geographical folds fixed before biology.
No sewage-specific wastewater flux, actual riparian biomass, or plant uptake.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,sys
import numpy as np,pandas as pd
from sklearn.linear_model import Ridge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
sys.path.insert(0,"analysis/20261117_RZ2018")
import eea2018_RZ2018_benthic_EQR_geoblock_model_if_fullsource as rz

MODEL_BASE=["numeric_EQR_2018","lat","lon","LogQ","LogArea"]
MODEL_L=MODEL_BASE+["NaturalFar","UrbanFar","CropFar"]
MODEL_D=MODEL_L+["NaturalNearMinusFar"]
ALPHA=10.
BOOTS=2000
PRELOCK_SCHEMA="GESUI_EEA2018_RZ2018_TO_EEA2022_BENTHIC_COHORT_SOURCE_ONLY_V1"
EXPECTED_FR=381
def source_quality(frame):
 # If an original 3-degree block is split across original 5folds, holdout
 # loses spatial independence, so STOP rather than repack folds after Y.
 blocks=frame.groupby("region_3deg_block_2018").pre_2018_source_only_spatial_fold.nunique()
 if (blocks>1).any():
  raise RuntimeError("2018 3-degree block leakage across original fixed folds")
 counts=frame.pre_2018_source_only_spatial_fold.value_counts()
 if len(counts)!=5 or min(counts)<20:raise RuntimeError("2018 frozen spatial folds not source-supported in France")
 return {"source_2018_original_five_fold_test_n":{
   str(int(k)):int(v) for k,v in counts.sort_index().items()},
   "2018_3degree_source_blocks":int(len(blocks)),
   "no_3deg_block_split_across_original_folds":True}
def read_2022(args):
 cols=["site_key","countryCode","procedureClassificationSystem","HYRIV_ID",
       "lat","lon","normalised_EQR","region_3deg_block"]
 t=pd.read_csv(args.eq22,usecols=cols,dtype={"site_key":str})
 if len(t)!=2475 or t.site_key.duplicated().any():raise RuntimeError("Original 2022 benthic numeric 2475 source changed")
 if not np.isfinite(t.normalised_EQR.to_numpy(float)).all() or (t.normalised_EQR<0).any():
  raise RuntimeError("2022 normalized original EQR source contains invalid values")
 return t
def frame(args):
 gate=json.loads(args.cohort_source.read_text())
 if gate["schema"]!=PRELOCK_SCHEMA or gate["source_frozen_acceptance"]!="FRANCE_ONLY":
  raise RuntimeError("Source-frozen France temporal cohort GO not valid")
 if gate["same2018RZ_2022_method_reach_coords_stable_n"]!=409 or gate["source_france_stable_site_n"]!=EXPECTED_FR:
  raise RuntimeError("Source-frozen France 381 site cohort has changed")
 if gate["source_original_2022_EQR_git_blob_SHA"]!="325f199e985d9862ffee5b0edb205c61a7be3833":
  raise RuntimeError("Signed native original 2022 Git blob changed")
 land,landmeta=rz.source_frame(args.rz)
 E,_=rz.source_eqr()
 t=read_2022(args)
 if len(E)!=2634 or len(land)!=978:raise RuntimeError("Original official land/benthic source changed")
 old=E.merge(land,on="site_key",how="inner",validate="one_to_one")
 if len(old)!=976:raise RuntimeError("Original 2018 EQR mapped RZ sample 976 changed")
 merged=old.merge(t,on="site_key",how="inner",suffixes=("_2018","_2022"),validate="one_to_one")
 proper=(merged.procedureClassificationSystem_2018.fillna("").ne("")&
         merged.procedureClassificationSystem_2018.eq(merged.procedureClassificationSystem_2022)&
         merged.HYRIV_ID_2018.eq(merged.HYRIV_ID_2022)&
         merged.countryCode_2018.eq(merged.countryCode_2022)&
         np.max(np.abs(merged[["lat_2018","lon_2018"]].to_numpy(float)-
                           merged[["lat_2022","lon_2022"]].to_numpy(float)),axis=1)<=.01)
 good=merged.loc[proper].copy()
 print("EEA_FR2018_22_BIOLOGICAL_RECAPTURE_PRECHECK",json.dumps({
   "2018_RZ_benthic_join":len(old),
   "2018_RZ_2022_all_revisit":len(merged),
   "method_equal_nonempty":int((merged.procedureClassificationSystem_2018.fillna("").ne("")&
               merged.procedureClassificationSystem_2018.eq(merged.procedureClassificationSystem_2022)).sum()),
   "same_reach":int(merged.HYRIV_ID_2018.eq(merged.HYRIV_ID_2022).sum()),
   "same_country":int(merged.countryCode_2018.eq(merged.countryCode_2022).sum()),
   "within_01deg_centroid":int((np.max(np.abs(merged[["lat_2018","lon_2018"]].to_numpy(float)-
               merged[["lat_2022","lon_2022"]].to_numpy(float)),axis=1)<=.01).sum()),
   "source_joint_prebiological_valid":len(good),
   "source_frozen_eligible":409,
   "source_frozen_FR_381":381},sort_keys=True),flush=True)
 if len(good)!=409:raise RuntimeError("Original source-frozen exact stable sites 409 changed")
 fr=good.loc[good.countryCode_2018.eq("FR")].copy()
 if len(fr)!=EXPECTED_FR or fr.site_key.duplicated().any():
  raise RuntimeError("Exact source-frozen France 381 stable-method sites changed")
 if not fr.procedureClassificationSystem_2018.eq(fr.procedureClassificationSystem_2022).all():
  raise RuntimeError("France BQE method changed after SOURCE cohort lock")
 fr["LogQ"]=np.log1p(fr.DIS_AV_CMS_2018.astype(float))
 fr["LogArea"]=np.log1p(fr.UPLAND_SKM_2018.astype(float))
 fr["EQR2022"]=fr.normalised_EQR.to_numpy(float)
 fr["EQR_delta"]=fr.EQR2022-fr.numeric_EQR_2018
 fr["lat"]=fr.lat_2018
 fr["lon"]=fr.lon_2018
 if not np.isfinite(fr[MODEL_D+["EQR2022","EQR_delta"]].to_numpy(float)).all():
  raise RuntimeError("Original 2018/2022 same-site EQR source data missing")
 return fr,source_quality(fr),gate
def oof(df,outcome,alpha=ALPHA):
 folds=df.pre_2018_source_only_spatial_fold.astype(int).to_numpy()
 y=df[outcome].to_numpy(float)
 predictions={}
 fold_info={}
 for name,col in (("B_2018_EQR_hydrogeo",MODEL_BASE),
                  ("L_plus_2018_outer_land",MODEL_L),
                  ("D_plus_near_outer_vegetated_contrast",MODEL_D)):
  p=np.full(len(df),np.nan)
  for i in range(1,6):
   tr=folds!=i;te=folds==i
   model=Pipeline([("std",StandardScaler()),("ridge",Ridge(alpha=alpha,solver="svd"))])
   model.fit(df.loc[tr,col].astype(float),y[tr])
   p[te]=model.predict(df.loc[te,col].astype(float))
  if not np.isfinite(p).all():raise RuntimeError("Geoblock-heldout model missing predictions")
  predictions[name]=p
  fold_info[name]={str(i):{"source_test_sites":int((folds==i).sum()),
    "MSE":float(np.mean((y[folds==i]-p[folds==i])**2))}
        for i in range(1,6)}
 return predictions,fold_info
def evaluate(df,target,pred):
 y=df[target].to_numpy(float)
 block=df.region_3deg_block_2018.astype(str).to_numpy()
 blocks=np.unique(block)
 losses={key:(y-v)**2 for key,v in pred.items()}
 if len(blocks)<5:raise RuntimeError("Original 3degree block count below source")
 rng=np.random.default_rng(20261123)
 draws=rng.integers(0,len(blocks),size=(BOOTS,len(blocks)))
 grouped={key:np.array([np.mean(v[block==z]) for z in blocks]) for key,v in losses.items()}
 comp={}
 for more,less,key in (("L_plus_2018_outer_land","B_2018_EQR_hydrogeo","L_after_lagged_EQR"),
                       ("D_plus_near_outer_vegetated_contrast","L_plus_2018_outer_land","D_after_L")):
  b=losses[less];f=losses[more]
  bb=grouped[less];ff=grouped[more]
  boot=100*(bb[draws].mean(axis=1)-ff[draws].mean(axis=1))/bb[draws].mean(axis=1)
  comp[key]={"original_true_2022_EQR_site_MSE_baseline":float(b.mean()),
    "original_true_2022_EQR_site_MSE_added":float(f.mean()),
    "site_equal_heldout_MSE_gain_percent":float(100*(b.mean()-f.mean())/b.mean()),
    "2018_3degree_region_equal_MSE_gain_percent":float(100*(bb.mean()-ff.mean())/bb.mean()),
    "fixed_OOF_region_bootstrap_descriptive95":[float(x) for x in np.quantile(boot,[.025,.975])],
    "regions_where_more_predictors_reduce_MSE":int((ff<bb).sum()),
    "total_regions":len(blocks),
    "no_causal_or_external_country_inference":True}
 return comp,{"source_3degree_region_count":len(blocks)}
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq22","cohort_source","out"):p.add_argument("--"+k,type=pathlib.Path,required=True)
 args=p.parse_args()
 d,qa,gate=frame(args)
 p,folds=oof(d,"EQR2022")
 result,_=evaluate(d,"EQR2022",p)
 # Internal state forecast using old EQR: direct no-fit persistence.
 persist=float(np.mean((d.EQR2022-d.numeric_EQR_2018)**2))
 # Sensitivity: fixed identical source/inputs but 2018→2022 change endpoint.
 # Do not compare this to raw EQR outcome as a separate independent study.
 pp,fold_change=oof(d,"EQR_delta")
 change,_=evaluate(d,"EQR_delta",pp)
 result_json={
   "schema":"GESUI_EEA2018_RZ_OFFICIAL_TO_EEA2022_FRENCH_BENTHOS_LAGGED_GEOFORECAST_V1",
   "status":"POSTHOC_FRANCE_ONLY_SOURCE_LOCKED_TEMPORAL_BENTHIC_FORECAST",
   "source_original_go":"FRANCE_ONLY",
   "2018_20122_true_matched_site_n":len(d),
   "original_2018_country":"FR",
   "original_2018_assessment_method_counts":d.procedureClassificationSystem_2018.value_counts().sort_index().to_dict(),
   "original_same_EEA_site_years_are_not_independent_observations":True,
   "source_frozen_RZ2018_true_land_site_2018_EQR2018_2022_match_qa":qa,
   "2018_land_mapped_NOT_2022_remeasured":True,
   "source_2018_upstream_reported_P_both_2022_status_not_eligible_for_source_change_effect":True,
   "source_2018_EEA_station_EQR_source_SHA256":"2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02",
   "source_2022_git_blob_SHA":"325f199e985d9862ffee5b0edb205c61a7be3833",
   "original_source_cohort_signed_FR_sites":gate["source_france_stable_site_n"],
   "source_lagged_2018_benthic_EQR_mean":float(d.numeric_EQR_2018.mean()),
   "source_real_2022_benthic_EQR_mean":float(d.EQR2022.mean()),
   "source_2018to2022_EQR_delta_mean":float(d.EQR_delta.mean()),
   "original_same_site_EQR_2018_2022_pearson":float(np.corrcoef(d.numeric_EQR_2018,d.EQR2022)[0,1]),
   "naive_EQR2018_persistence_predicts_2022_MSE":persist,
   "model_alpha":ALPHA,"source_train_only_scaling":True,
   "original_5fixed_2018_spatial_folds":qa["source_2018_original_five_fold_test_n"],
   "original_2022_future_EQR_forecast_regional_errors":folds,
   "original_2022_future_EQR_forecast_incremental_source":result,
   "secondary_same_site_delta_EQR_sensitivity":change,
   "posthoc_observational_not_vegetation_buffer_causation":True,
   "no_WWTP_effluent_exposure_change_or_insect_emergence_measured":True}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(result_json,indent=2,sort_keys=True)+"\n")
 print("EEA2018_RZ_TO_FRANCE2022_BENTHIC_TEMPORAL_FORECAST",json.dumps({
   "source_n":len(d),"folds":qa["source_2018_original_five_fold_test_n"],
   "EQ2018_2022_mean":[result_json["source_lagged_2018_benthic_EQR_mean"],result_json["source_real_2022_benthic_EQR_mean"]],
   "correlation":result_json["original_same_site_EQR_2018_2022_pearson"],
   "future":result,"change":change},sort_keys=True),flush=True)
if __name__=="__main__":main()
