#!/usr/bin/env python3
"""EPA source-only novelty HUC8 rematching on substrate or hydro-landscape.

Reconstructs old EPA source panels only to freeze source-HUC membership.
Does NOT READ any 2023 benthic counts or EPT responses, no effects fitted.
"""
import argparse,collections,hashlib,json,math,pathlib
import numpy as np,pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as prev
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as dom
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

BASE=dom.MATCH_FIELDS
BLOCKS={"grain":["PCT_FN","PCT_SAFN","LSUB_DMM"],
        "hydrology_landuse":["log_mean_flow","log_wsarea","urban_fraction"]}
BLOCK_CALIPER=3.5
def block_sources(args,env):
 h=pd.read_csv(args.habitat,dtype=str,keep_default_na=False,low_memory=False).set_index("UID")
 site=pd.read_csv(args.sites,dtype=str,keep_default_na=False,low_memory=False).set_index("UID")
 if any(k not in h.columns for k in BLOCKS["grain"]):
  raise RuntimeError("Original EPA source lacks frozen PCT_FN/PCT_SAFN/LSUB_DMM")
 if any(k not in site.columns for k in ("MEAN_FLOW","WSAREASQKM","URBN_NRS23")):
  raise RuntimeError("Original EPA source lacks frozen hydrology/urban fields")
 raw=pd.DataFrame(index=env.index)
 for k in BLOCKS["grain"]:raw[k]=pd.to_numeric(h.loc[env.index,k],errors="coerce")
 q=pd.to_numeric(site.loc[env.index,"MEAN_FLOW"],errors="coerce")
 a=pd.to_numeric(site.loc[env.index,"WSAREASQKM"],errors="coerce")
 raw["log_mean_flow"]=np.log10(q.where(q>0))
 raw["log_wsarea"]=np.log10(a.where(a>0))
 raw["urban_fraction"]=site.loc[env.index,"URBN_NRS23"].map({"Urban":1.0,"Non-Urban":0.0})
 for k in ("PCT_FN","PCT_SAFN"):raw.loc[~raw[k].between(0,100),k]=np.nan
 return raw
def normalizer(oldtrain,target,cols):
 centers=np.empty(len(cols));scale=np.empty(len(cols))
 for i,k in enumerate(cols):
  # the old three-cycle original files lack grain and modeled flow metadata;
  # these supplementary features are standardized from ORIGINAL combined
  # 2023 HUC8 environment ONLY, without any benthic outcome information.
  src=oldtrain[k] if k in oldtrain.columns else target[k]
  centers[i]=float(src.mean())
  scale[i]=float(src.std(ddof=0))
  if not (np.isfinite(scale[i]) and scale[i]>1e-8):scale[i]=1.
 return centers,scale
def choose(target,oldtrain,ever,features):
 a=target.loc[~ever].dropna(subset=features).copy()
 b=target.loc[ever].dropna(subset=features).copy()
 center,scale=normalizer(oldtrain,target,features)
 pairs=[];ecoregion_counts={}
 for eco in sorted(target.eco9.unique()):
  u=a[a.eco9.eq(eco)];v=b[b.eco9.eq(eco)]
  if len(u)==0 or len(v)==0:
   ecoregion_counts[eco]={"never":len(u),"ever":len(v),"matched":0}
   continue
  mu=(u[features].to_numpy(float)-center)/scale
  mv=(v[features].to_numpy(float)-center)/scale
  dist=cdist(mu,mv,"euclidean")
  rows,cols=linear_sum_assignment(dist)
  keep=[(str(u.index[i]),str(v.index[j]),eco,float(dist[i,j])) for i,j in zip(rows,cols)
        if dist[i,j]<=BLOCK_CALIPER]
  pairs.extend(keep)
  ecoregion_counts[eco]={"never":len(u),"ever":len(v),"matched":len(keep)}
 pairs.sort()
 raw="\n".join(f"{n}|{o}|{eco}|{d:.12f}" for n,o,eco,d in pairs)
 digest=hashlib.sha256(raw.encode()).hexdigest()
 return pairs,ecoregion_counts,digest
def smd(x,y):
 a=np.asarray(x,float);b=np.asarray(y,float)
 v=np.sqrt((a.var(ddof=0)+b.var(ddof=0))/2.)
 return float((a.mean()-b.mean())/v) if v>1e-10 else None
def pair_balance(target,pairs,features):
 if not pairs:return {}
 n=target.loc[[p[0] for p in pairs]]
 e=target.loc[[p[1] for p in pairs]]
 return {k:{"novel_minus_ever_SMD":smd(n[k],e[k]),
   "novel_mean":float(n[k].mean()),"ever_mean":float(e[k].mean())}
  for k in features}
def gradations(target):
 X=target[["PCT_FN","PCT_SAFN","LSUB_DMM"]].to_numpy(float)
 y=target.XEMBED.to_numpy(float)
 std=np.std(X,axis=0);std=np.where(std>1e-8,std,1.)
 z=(X-X.mean(0))/std
 Z=np.c_[np.ones(len(z)),z]
 coef=np.linalg.lstsq(Z,y,rcond=None)[0]
 pred=Z@coef
 err=((pred-y)**2).sum()
 den=((y-y.mean())**2).sum()
 return {"HUC8_equal_multivariable_R2_XEMBED_from_grain":float(1-err/den),
   "XEMBED_vs_PCT_FN_Pearson":float(np.corrcoef(y,X[:,0])[0,1]),
   "XEMBED_vs_PCT_SAFN_Pearson":float(np.corrcoef(y,X[:,1])[0,1]),
   "XEMBED_vs_LSUB_DMM_Pearson":float(np.corrcoef(y,X[:,2])[0,1]),
   "regression_is_source_only_descriptive_not_causal":True}
def run(args):
 for k in ("sites","water","habitat","chla"):
  if sha(getattr(args,k))!=SHA[k]:raise RuntimeError("Original EPA source SHA mismatch "+k)
 old=prev.history(args.old_panel)
 env=prev.covariates(args,{})
 if len(env)!=1165 or env.huc8.nunique()!=693 or len(old)!=2648 or old.huc8.nunique()!=1030:
  raise RuntimeError("Old or target original source population changed")
 target=dom.huc_aggregate(env)
 oldtrain=old.loc[~old.huc8.isin(set(target.index))].copy()
 oldh=dom.huc_aggregate(oldtrain)
 if len(oldtrain)!=1185 or len(oldh)!=551:raise ValueError("Older training source HUC group changed")
 if set(target.index)&set(oldh.index):raise RuntimeError("HUC8 training leakage")
 x=block_sources(args,env)
 # Same HUC8-equal environmental grain as the original preoutcome source-only match.
 ex=x.join(env[["huc8"]]).groupby("huc8").mean()
 extra_cols=[k for b in BLOCKS.values() for k in b]
 target=target.join(ex[extra_cols],how="left")
 ever=target.index.isin(set(old.huc8))
 if (ever.sum(),(~ever).sum())!=(479,214):raise ValueError("HUC8 source novelty grouping changed")
 basepairs,baseeco,basedigest=dom.choose_pairs(oldh,target,ever)
 if basedigest!="356105578fc8fe7c7deda89fda373b371fec0989bee3ad0be028fd3eee195366":
  raise RuntimeError("Original 214 HUC source matching digest changed")
 outcome={}
 for name,block in BLOCKS.items():
  cols=BASE+block
  ever_full=target.loc[ever,block].notna().all(axis=1)
  never_full=target.loc[~ever,block].notna().all(axis=1)
  source_go=(int(ever_full.sum())>=int(np.ceil(.8*479)) and
             int(never_full.sum())>=int(np.ceil(.8*214)))
  matches,ecocount,digest=choose(target,oldh,ever,cols)
  eco=int(sum(x["matched"]>0 for x in ecocount.values()))
  accepted=bool(source_go and len(matches)>=70 and len(matches)>=int(np.ceil(.6*214)) and eco>=5)
  outcome[name]={"original_fixed_environment_variables":cols,
    "source_ever_eligible_HUC8":int(ever_full.sum()),
    "source_novel_eligible_HUC8":int(never_full.sum()),
    "source_ever_complete_fraction":float(ever_full.mean()),
    "source_novel_complete_fraction":float(never_full.mean()),
    "source_only_match_pairs":len(matches),
    "source_match_eco9":eco,
    "source_match_ecoregion_counts":ecocount,
    "source_match_pair_membership_SHA256":digest,
    "source_match_distance_q25_median_q75":
      [float(np.quantile([p[3] for p in matches],q)) for q in (.25,.5,.75)] if matches else None,
    "matched_environment_balance":pair_balance(target,matches,cols),
    "source_go":accepted,
    "source_stop_reason":None if accepted else "frozen_80pct_both_HUC_70_pair_60pct_novel_and_5_eco9_failed"}
 valid=target.dropna(subset=BLOCKS["grain"])
 outcome_grain=gradations(valid)
 output={"schema":"NRSA2023_original_four_cycle_HUC_novelty_extra_source_domains_v1",
  "original_EPA_SHA":{k:SHA[k] for k in ("sites","water","habitat","chla")},
  "original_2023_UID":1165,"original_2023_HUC8":693,
  "original_old_source_three_cycles_UID":2648,"old_source_HUC8_union":1030,
  "original_training_UID":1185,"old_training_HUC8":551,
  "original_prior_7feature_match_pair_SHA256":basedigest,
  "source_only_blocks":outcome,
  "original_XEMBED_compositional_source_relationships":outcome_grain,
  "no_EPT_biological_outcomes_or_previous_errors_used_to_choose_match_pairs":True,
  "posthoc_matching_features_known_from_ecology_but_not_outcome_signs":True}
 return output
def main():
 ap=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","old_panel","out"):
  ap.add_argument("--"+k,type=pathlib.Path,required=True)
 a=ap.parse_args()
 result=run(a)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("NRSA2023_SUBSTRATE_HYDRO_HUC_NOVELTY_SOURCE_GATE",json.dumps({
  "blocks":result["source_only_blocks"],
  "embedded_grain":result["original_XEMBED_compositional_source_relationships"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
