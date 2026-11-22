#!/usr/bin/env python3
"""Evaluate SOURCE-ONLY signed additional original EPA substrate/hydro HUC matches.

Predictions remain fixed: original past-era EPA old-3 models, target 2023 unseen
from training. New HUC matches selected without target EPT outcomes.
"""
import argparse,hashlib,json,pathlib
import numpy as np,pandas as pd
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as prev
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as dom
import gesui_nrsa2023_huc_novelty_extra_source_matching_20261122 as extra
import gesui_nrsa2023_novel_HUC8_source_matched_EPT_prediction_20261122 as outcome
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

def close(a,b,where):
 if isinstance(a,bool) or isinstance(b,bool):
  if a!=b:raise RuntimeError("Source contract "+where)
 elif isinstance(a,(int,float)) and isinstance(b,(int,float)):
  if not np.isclose(a,b,rtol=1e-12,atol=1e-12):
   raise RuntimeError("Source contract "+where+str((a,b)))
 elif isinstance(a,list):
  if len(a)!=len(b):raise RuntimeError("Source list count "+where)
  for i,(x,y) in enumerate(zip(a,b)):close(x,y,where+str(i))
 elif isinstance(a,dict):
  if set(a)!=set(b):raise RuntimeError("Source keys "+where)
  for k in a:close(a[k],b[k],where+"."+str(k))
 elif a!=b:raise RuntimeError("Source exact "+where)
def main():
 p=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","old_panel","source_qa",
           "original_result","extension_source","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 signed=json.loads(a.extension_source.read_text())
 assert signed["schema"]=="NRSA2023_original_four_cycle_HUC_novelty_extra_source_domains_v1"
 for k in ("sites","water","habitat","chla","benthos"):
  if sha(getattr(a,k))!=SHA[k]:raise RuntimeError("Original published EPA source SHA mismatch "+k)
 current=extra.run(a)
 for k in ("source_only_blocks","original_XEMBED_compositional_source_relationships",
           "original_prior_7feature_match_pair_SHA256",
           "original_2023_UID","original_2023_HUC8",
           "original_old_source_three_cycles_UID","old_source_HUC8_union",
           "original_training_UID","old_training_HUC8"):
  close(current[k],signed[k],k)
 old=prev.history(a.old_panel)
 orig_source=json.loads(a.source_qa.read_text())
 oldenv=prev.covariates(a,orig_source)
 if len(oldenv)!=1165 or oldenv.huc8.nunique()!=693:raise RuntimeError("Old source panel target mismatch")
 ht=dom.huc_aggregate(oldenv)
 train=old.loc[~old.huc8.isin(set(ht.index))].copy()
 trainh=dom.huc_aggregate(train)
 if len(train)!=1185 or train.huc8.nunique()!=551:raise RuntimeError("Older 3-cycle training cohort differs")
 grain=extra.block_sources(a,oldenv)
 grp=grain.join(oldenv[["huc8"]]).groupby("huc8").mean()
 ht=ht.join(grp[[k for z in extra.BLOCKS.values() for k in z]],how="left")
 ever=ht.index.isin(set(old.huc8))
 if (ever.sum(),(~ever).sum())!=(479,214):raise RuntimeError("Native HUC groups changed")
 # Open original counted EPT only after source-only all-pairs identity verified
 # with previously signed source set and original 2023 target.
 targ=prev.canonical_target(prev.make_target(a,orig_source),orig_source,"full")
 ref=json.loads(a.original_result.read_text())
 assert ref["schema"]=="gesui_original_EPA_4thcycle_2023_24_EPT_NP_PHab_outofera_HUC8_holdout_v1"
 if len(targ)!=1165 or targ.huc8.nunique()!=693:raise RuntimeError("Original 2023 target changed")
 if set(train.huc8)&set(targ.huc8):raise RuntimeError("HUC8 train target leakage")
 pred=prev.all_predictions(train,targ)
 outcome.fixed_prediction_reference(pred,targ,ref["results"]["full"])
 losses={key:outcome.huc_loss(targ,z) for key,z in pred.items()}
 blocks={}
 for name,cols in extra.BLOCKS.items():
  source=signed["source_only_blocks"][name]
  block={"source_go":source["source_go"],
    "source_matched_HUC8_pairs":source["source_only_match_pairs"],
    "source_match_digest":source["source_match_pair_membership_SHA256"]}
  if not source["source_go"]:
   block["status"]="SOURCE_STOP_DO_NOT_FIT_OR_REPORT_OUTCOME";blocks[name]=block;continue
  pairs,_,digest=extra.choose(ht,trainh,ever,extra.BASE+cols)
  if digest!=source["source_match_pair_membership_SHA256"] or len(pairs)!=source["source_only_match_pairs"]:
   raise RuntimeError("Frozen pre-EPT source-only "+name+" matched HUC8 pair identity changed")
  r=outcome.paired_errors(losses,pairs)
  block["status"]="POSTHOC_FIXED_OOF_MATCHED_COMPLETED"
  block["HUC8_equal_matched_loss_contrasts"]=r
  block["matched_balance"]=extra.pair_balance(ht,pairs,extra.BASE+cols)
  block["covariate_standardization_source"]="original_old_training_HUC8_for_original7_and_all_2023_HUC8_for_new_features"
  blocks[name]=block
  print("NRSA_NEW_SUBSTRATE_OR_HYDRO_MATCHED_EPT",name,json.dumps(
   {"n_source_matched_HUC8_pairs":len(pairs),
    "beta_metric":"heldout percent error reduction, NOT causal effect",
    "contrasts":{k:{"novel_gain_pct":v["novel_gain_pct"],
      "seen_gain_pct":v["previously_seen_source_gain_pct"],
      "novel_minus_seen_gain_pctpoints":v["novel_minus_seen_gain_pct_points"],
      "difference_bootstrap95":v["matched_pair_bootstrap95_gain_difference"]}
     for k,v in r["comparisons"].items()}},sort_keys=True),flush=True)
 output={"schema":"NRSA2023_original_HUC8_novelty_extra_substrate_hydro_matched_EPT_posthoc_v1",
  "new_source_environment_SHA256":hashlib.sha256(a.extension_source.read_bytes()).hexdigest(),
  "source_only_2023_matching_blocks_frozen":True,
  "source_2023_UID":1165,"source_2023_HUC8":693,
  "historic_3cycle_train_UID":1185,"historic_train_HUC8":551,
  "no_target_2023_EPT_in_old_training":True,
  "original_full_target_MSE_and_gains_reproduced":True,
  "source_original_new_extra_match_blocks":blocks,
  "exploratory_posthoc_not_causal_or_sewage_or_pollution_conduit":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("NRSA_EXTRA_BLOCKS_HUC_NOVELTY_POSTHOC_SUMMARY",json.dumps({
  k:{"source_go":v["source_go"],"pairs":v["source_matched_HUC8_pairs"],
     "status":v["status"],
     "bed_alone":{key:q for key,q in v.get("HUC8_equal_matched_loss_contrasts",{}).get("comparisons",{}).get("bed_alone_after_NP",{}).items()
       if key in ["novel_gain_pct","previously_seen_source_gain_pct","novel_minus_seen_gain_pct_points","matched_pair_bootstrap95_gain_difference"]}}
  for k,v in blocks.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
