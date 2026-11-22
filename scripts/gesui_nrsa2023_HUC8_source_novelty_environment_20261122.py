#!/usr/bin/env python3
"""Original EPA 2023 fourth-cycle SOURCE ENVIRONMENT, no EPT outcome test.

Compare HUC8 once represented in 2008/13/18 pre-selected source cohorts to
HUC8 never represented there; no biological numbers used in these calculations.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib
import numpy as np,pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as previous

FIELDS=["lat","lon","stream_order","log_cond","log_tn","log_tp",
        "XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
DISTANCE_BLOCKS={
 "geography":["lat","lon","stream_order"],
 "eco_environment":["lat","lon","stream_order","log_tn","log_tp","XEMBED","log_width"]
}
MATCH_FIELDS=DISTANCE_BLOCKS["eco_environment"]
CALIPER=3.5
def modal(series):
 counts=series.astype(str).value_counts()
 return sorted(counts[counts==counts.max()].index)[0]
def huc_aggregate(df):
 cols=["huc8","eco9","site_id"]+FIELDS
 if any(col not in df.columns for col in cols):raise RuntimeError("Missing exact EPA covariate")
 if not np.isfinite(df[FIELDS].to_numpy(float)).all():
  raise RuntimeError("Nonfinite source environmental covariate before group mean")
 z=df.groupby("huc8",sort=True)[FIELDS].mean()
 z["eco9"]=df.groupby("huc8").eco9.apply(modal)
 z["n_UID"]=df.groupby("huc8").size()
 z["n_SITE_ID"]=df.groupby("huc8").site_id.nunique()
 z["n_eco9_per_HUC8"]=df.groupby("huc8").eco9.nunique()
 return z
def standardized(train,a,b,fields):
 center=train[fields].mean(axis=0).to_numpy(float)
 sd=train[fields].std(axis=0,ddof=0).to_numpy(float)
 sd=np.where(sd>1e-9,sd,1.)
 return (a[fields].to_numpy(float)-center)/sd,(b[fields].to_numpy(float)-center)/sd
def source_geometry(train,target,fields):
 xa,xb=standardized(train,train,target,fields)
 ec_a=train.eco9.to_numpy(str);ec_b=target.eco9.to_numpy(str)
 n1=np.full(len(train),np.nan);n2=np.full(len(target),np.nan)
 for eco in sorted(set(ec_a)|set(ec_b)):
  tr=np.where(ec_a==eco)[0];te=np.where(ec_b==eco)[0]
  if len(tr)>1:
   mat=cdist(xa[tr],xa[tr],"euclidean")
   np.fill_diagonal(mat,np.inf)
   n1[tr]=mat.min(axis=1)
  if len(tr) and len(te):
   n2[te]=cdist(xb[te],xa[tr],"euclidean").min(axis=1)
 cutoff=float(np.quantile(n1[np.isfinite(n1)],.95))
 return n2,cutoff,int(np.sum(~np.isfinite(n2)))
def summary(v):
 arr=np.asarray(v,dtype=float)
 arr=arr[np.isfinite(arr)]
 if len(arr)==0:return {"n":0}
 return {"n":len(arr),"mean":float(np.mean(arr)),
   "sd":float(np.std(arr,ddof=0)),"median":float(np.median(arr)),
   "q10":float(np.quantile(arr,.1)),"q90":float(np.quantile(arr,.9))}
def source_profile(target,ever,train):
 records={}
 a=target.loc[ever];b=target.loc[~ever]
 for k in FIELDS:
  aa=a[k].to_numpy(float);bb=b[k].to_numpy(float)
  pooled=np.sqrt((aa.var(ddof=0)+bb.var(ddof=0))/2)
  records[k]={"old_source_ever":summary(aa),
    "never_old_source":summary(bb),
    "novel_minus_ever_standardized_mean_difference":
      float((bb.mean()-aa.mean())/pooled) if pooled>1e-10 else None,
    "historical_training":summary(train[k].to_numpy(float))}
 return records
def choose_pairs(train,target,ever):
 old=target.loc[ever].copy();novel=target.loc[~ever].copy()
 center=train[MATCH_FIELDS].mean(axis=0).to_numpy(float)
 sd=train[MATCH_FIELDS].std(axis=0,ddof=0).to_numpy(float)
 sd=np.where(sd>1e-9,sd,1.)
 matched=[];per_eco={}
 for eco in sorted(target.eco9.unique()):
  oo=old.loc[old.eco9==eco];nn=novel.loc[novel.eco9==eco]
  if len(oo)==0 or len(nn)==0:
   per_eco[eco]={"old":len(oo),"novel":len(nn),"pairs":0}
   continue
  a=(nn[MATCH_FIELDS].to_numpy(float)-center)/sd
  b=(oo[MATCH_FIELDS].to_numpy(float)-center)/sd
  dist=cdist(a,b,"euclidean")
  idx_a,idx_b=linear_sum_assignment(dist)
  for x,y in zip(idx_a,idx_b):
   dd=float(dist[x,y])
   if dd<=CALIPER:matched.append((str(nn.index[x]),str(oo.index[y]),eco,dd))
  per_eco[eco]={"old":len(oo),"novel":len(nn),
    "pairs":sum(1 for x in matched if x[2]==eco)}
 matched.sort()
 raw="\n".join(f"{n}|{o}|{eco}|{d:.12f}" for n,o,eco,d in matched)
 digest=hashlib.sha256(raw.encode()).hexdigest()
 return matched,per_eco,digest
def build(args):
 historic=previous.history(args.old_panel)
 # covariates() joins SHA-pinned SiteInfo, WaterChem, Phab and periphyton source
 # availability, not EPT values. The signed source-only 2023 cohort is n1165.
 env=previous.covariates(args,{}).copy()
 if len(env)!=1165 or env.huc8.nunique()!=693 or env.site_id.nunique()!=1046:
  raise RuntimeError("Original 2023 outcome-blind source n1165/693/1046 has changed")
 if len(historic)!=2648 or historic.huc8.nunique()!=1030:
  raise RuntimeError("Older three-cycle source census changed")
 old_union=set(historic.huc8)
 target=huc_aggregate(env);ever=target.index.isin(old_union)
 if (int(ever.sum()),int((~ever).sum()))!=(479,214):
  raise RuntimeError("Original HUC8 source novelty membership counts changed")
 if (int(target.loc[ever,"n_UID"].sum()),int(target.loc[~ever,"n_UID"].sum()))!=(879,286):
  raise RuntimeError("Original source novelty UID count changed")
 train=historic.loc[~historic.huc8.isin(set(target.index))].copy()
 if len(train)!=1185 or train.huc8.nunique()!=551:
  raise RuntimeError("Previously held-out 3-era HUC8 train source cohort changed")
 old_train=huc_aggregate(train)
 if set(target.index)&set(old_train.index):raise RuntimeError("Source training target HUC8 leakage")
 profiles=source_profile(target,ever,old_train)
 domain={}
 for name,cols in DISTANCE_BLOCKS.items():
  distances,cutoff,missing=source_geometry(old_train,target,cols)
  domain[name]={
   "original_train_HUC8_nearest_HUC8_LOO95_cutoff":cutoff,
   "no_same_ecoregion_source_train":missing,
   "old_source_ever":{"distance":summary(distances[ever]),
      "fraction_above_training_LOO95":float(np.mean(distances[ever]>cutoff))},
   "never_old_source":{"distance":summary(distances[~ever]),
      "fraction_above_training_LOO95":float(np.mean(distances[~ever]>cutoff))}}
 matched,by_eco,digest=choose_pairs(old_train,target,ever)
 pair_count=len(matched)
 source_go=pair_count>=70 and len([x for x in by_eco.values() if x["pairs"]>0])>=5 and pair_count>=int(np.ceil(.6*214))
 out={
  "schema":"gesui_NRSA2023_original_source_only_HUC8_novelty_env_domain_v1",
  "source_2023_original_UID":len(env),"source_2023_HUC8":len(target),
  "historical_source_original_UID":len(historic),
  "historical_source_HUC8_union":len(old_union),
  "historical_training_n_after_target_HUC8_removal":len(train),
  "historical_training_HUC8":len(old_train),
  "training_2023_shared_HUC8":0,
  "source_2023_groups":{
   "ever_prior_source_HUC8":{"UID":879,"HUC8":479},
   "never_prior_source_HUC8":{"UID":286,"HUC8":214}},
  "region_counts":{
   group:target.loc[mask].eco9.value_counts().sort_index().to_dict()
   for group,mask in (("ever_old_source",ever),("never_old_source",~ever))},
  "multiple_eco9_within_single_HUC8":{
   "ever":int((target.loc[ever,"n_eco9_per_HUC8"]>1).sum()),
   "never":int((target.loc[~ever,"n_eco9_per_HUC8"]>1).sum())},
  "environment_HUC8_equal_profiles":profiles,
  "distance_to_historical_train_HUC8":domain,
  "within_ecoregion_source_only_match":{
   "algorithm":"one_to_one_Hungarian_Euclidean_no_replacement_then_fixed_caliper",
   "numerical_features":MATCH_FIELDS,
   "training_HUC8_only_feature_scaling":True,
   "fixed_max_standardized_distance":CALIPER,
   "matched_pair_count":pair_count,
   "fraction_of_original_novel_HUC8_matched":pair_count/214,
   "matched_ecoregion_count":sum(x["pairs"]>0 for x in by_eco.values()),
   "ecoregion_match_source_counts":by_eco,
   "source_pair_membership_and_distance_SHA256":digest,
   "source_go":source_go,
   "source_stop_reason":None if source_go else "PRE_FIXED_MIN_70_PAIRS_5_ECO9_AND_60PCT_NOVEL_FAILED"},
  "source_only_no_EPT_outcomes_modeled":True,
  "HUC8_novelty_is_past_source_sample_membership_not_true_ecological_novelty":True}
 return out
def main():
 parser=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","old_panel","out"):
  parser.add_argument("--"+k,type=pathlib.Path,required=True)
 x=parser.parse_args()
 # Source SHA originals were already checked in original upstream workflow.
 # Check direct source digests independently here.
 from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha
 for key in ("sites","water","habitat","chla"):
  if sha(getattr(x,key))!=SHA[key]:
   raise RuntimeError("EPA original source SHA changed "+key)
 data=build(x)
 x.out.parent.mkdir(parents=True,exist_ok=True)
 x.out.write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
 print("NRSA_NOVEL_HUC8_ENV_SOURCE_ONLY",json.dumps({
  "groups":data["source_2023_groups"],
  "region_counts":data["region_counts"],
  "distance":data["distance_to_historical_train_HUC8"],
  "standardized_mean_differences":{k:v["novel_minus_ever_standardized_mean_difference"]
    for k,v in data["environment_HUC8_equal_profiles"].items()},
  "matched":data["within_ecoregion_source_only_match"]
  },sort_keys=True),flush=True)
if __name__=="__main__":main()
