#!/usr/bin/env python3
"""Posthoc source-only HUC8 ecological-domain matching + ORIGINAL fixed EPT predictions.

Match source original 2023 HUC8 without EPT outcomes, then verify frozen
source SHA and prior full-target original three-era Ridge predictions.
No refits by 2023 outcomes, no claims of causal novelty or wastewater.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib
import numpy as np,pandas as pd
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as original
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as domain

PRIMARY=[("PHab_after_NP","B+NP","B+NP+PHab"),
         ("bed_after_NP_other5","B+NP+PHab_except_bed","B+NP+PHab"),
         ("bed_alone_after_NP","B+NP","B+NP+bed"),
         ("NP_after_PHab","B+PHab","B+NP+PHab")]
def fixed_prediction_reference(preds,target,reference):
 evaled=original.evaluate(target,preds)
 for name,r in reference["target"]["model_errors"].items():
  v=evaled["model_errors"][name]["MSE_equal_HUC8"]
  if abs(v-r["MSE_equal_HUC8"])>max(1e-10,1e-8*abs(r["MSE_equal_HUC8"])):
   raise RuntimeError("Original EPT full target model MSE failed: "+name)
 for name,r in reference["target"]["comparisons"].items():
  v=evaled["comparisons"][name]["HUC8_equal_gain_pct"]
  if abs(v-r["HUC8_equal_gain_pct"])>1e-9:
   raise RuntimeError("Frozen previous full 2023 response comparison failed: "+name)
 return evaled
def huc_loss(target,pred):
 y=target.ept_fraction.to_numpy(float)
 return pd.Series((y-pred)**2,index=target.huc8.to_numpy()).groupby(level=0).mean()
def pct(a,b):
 return float(100*(a.mean()-b.mean())/a.mean())
def summarize(a):
 z=np.asarray(a,dtype=float)
 return {"n":int(len(z)),"mean":float(z.mean()),"median":float(np.median(z)),
         "q05":float(np.quantile(z,.05)),"q95":float(np.quantile(z,.95))}
def paired_errors(losses,pairs):
 outputs={}
 m=len(pairs)
 group_eco=collections.Counter([p[2] for p in pairs])
 results={"pairs":m,"matched_by_eco9":dict(sorted(group_eco.items())),
         "distance":summarize([p[3] for p in pairs]),
         "source_pair_alignment_not_outcome_selected":True}
 for name,base,more in PRIMARY:
  nov=np.array([[losses[base][p[0]],losses[more][p[0]]] for p in pairs],dtype=float)
  ever=np.array([[losses[base][p[1]],losses[more][p[1]]] for p in pairs],dtype=float)
  if not (np.isfinite(nov).all() and np.isfinite(ever).all()):raise ValueError("Missing model loss in matched pair")
  g_n=pct(nov[:,0],nov[:,1]);g_o=pct(ever[:,0],ever[:,1])
  rng=np.random.default_rng(int(hashlib.sha256(name.encode()).hexdigest()[:8],16))
  idx=rng.integers(0,m,size=(2000,m))
  n0=nov[:,0][idx].mean(1);n1=nov[:,1][idx].mean(1)
  o0=ever[:,0][idx].mean(1);o1=ever[:,1][idx].mean(1)
  gbn=100*(n0-n1)/n0;gbo=100*(o0-o1)/o0
  paired_diff=gbn-gbo
  abs_reduction=(nov[:,0]-nov[:,1])-(ever[:,0]-ever[:,1])
  boot_abs=abs_reduction[idx].mean(1)
  outputs[name]={
   "novel_HUC8":m,"matched_ever_HUC8":m,
   "novel_reference_MSE":float(nov[:,0].mean()),
   "novel_added_MSE":float(nov[:,1].mean()),
   "novel_gain_pct":g_n,
   "novel_bootstrap_pair95":[float(z) for z in np.quantile(gbn,[.025,.975])],
   "previously_seen_source_reference_MSE":float(ever[:,0].mean()),
   "previously_seen_source_added_MSE":float(ever[:,1].mean()),
   "previously_seen_source_gain_pct":g_o,
   "previously_seen_source_bootstrap_pair95":[float(z) for z in np.quantile(gbo,[.025,.975])],
   "novel_minus_seen_gain_pct_points":float(g_n-g_o),
   "matched_pair_bootstrap95_gain_difference":[float(z) for z in np.quantile(paired_diff,[.025,.975])],
   "novel_minus_seen_absolute_MSE_reduction":float(abs_reduction.mean()),
   "matched_pair_bootstrap95_absolute_MSE_difference":[float(z) for z in np.quantile(boot_abs,[.025,.975])],
   "not_a_randomized_treatment_contrast":True}
 results["comparisons"]=outputs
 return results
def numeric_match_profiles(old_huc,target_huc,pairs):
 records={}
 matched_ever=target_huc.loc[[p[1] for p in pairs]]
 matched_novel=target_huc.loc[[p[0] for p in pairs]]
 for k in domain.MATCH_FIELDS:
  oo=matched_ever[k].to_numpy(float);nn=matched_novel[k].to_numpy(float)
  un_old=target_huc.loc[target_huc.index.isin(old_huc)][k].to_numpy(float)
  un_nov=target_huc.loc[~target_huc.index.isin(old_huc)][k].to_numpy(float)
  def smd(a,b):
   v=np.sqrt((a.var(ddof=0)+b.var(ddof=0))/2)
   return float((b.mean()-a.mean())/v) if v>1e-10 else None
  records[k]={"unmatched_novel_minus_seen_SMD":smd(un_old,un_nov),
              "matched_novel_minus_seen_SMD":smd(oo,nn)}
 return records
def main():
 parser=argparse.ArgumentParser()
 for k in ("sites","water","habitat","benthos","chla","old_panel",
           "source_qa","original_result","novelty_source","out"):
  parser.add_argument("--"+k,type=pathlib.Path,required=True)
 args=parser.parse_args()
 signed=json.loads(args.novelty_source.read_text())
 if signed["schema"]!="gesui_NRSA2023_original_source_only_HUC8_novelty_env_domain_v1":
  raise ValueError("Unexpected source-only HUC8 ecological support contract")
 match_state=signed["within_ecoregion_source_only_match"]
 if not match_state["source_go"]:
  print("NRSA_NOVELTY_MATCHED_OUTCOME_SOURCE_STOP",json.dumps(match_state,sort_keys=True))
  result={"schema":"NRSA2023_HUC8_novelty_matched_outcomes_STOP_v1",
          "source_matched_GO":False,"source_reason":match_state["source_stop_reason"]}
  args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
  return
 origsource=json.loads(args.source_qa.read_text())
 ref=json.loads(args.original_result.read_text())
 old=original.history(args.old_panel)
 if len(old)!=2648 or old.huc8.nunique()!=1030:raise RuntimeError("Original 3-cycle old sample changed")
 # Read the exactly same source-only environmental candidate and reproduce
 # the old SOURCE-ONLY signed matched pair hash before any EPT outcome use.
 env=original.covariates(args,origsource)
 huc_target=domain.huc_aggregate(env)
 train=old.loc[~old.huc8.isin(set(huc_target.index))].copy()
 train_h=domain.huc_aggregate(train)
 ever=huc_target.index.isin(set(old.huc8))
 pairs,byeco,pair_sha=domain.choose_pairs(train_h,huc_target,ever)
 if pair_sha!=match_state["source_pair_membership_and_distance_SHA256"] or len(pairs)!=match_state["matched_pair_count"]:
  raise RuntimeError("Original source-only pair membership or distance differs before EPT model")
 if len(train)!=1185 or train.huc8.nunique()!=551:raise ValueError("Original three old-cycle train source changed")
 # Original EPA source-only gate was signed in an entirely separate Actions run.
 target=original.canonical_target(original.make_target(args,origsource),origsource,"full")
 if len(target)!=1165 or target.huc8.nunique()!=693:raise ValueError("Original 2023 heldout cohort changed")
 if set(train.huc8)&set(target.huc8):raise RuntimeError("Target HUC8 leakage")
 preds=original.all_predictions(train,target)
 fixed_prediction_reference(preds,target,ref["results"]["full"])
 losses={k:huc_loss(target,v) for k,v in preds.items()}
 result=paired_errors(losses,pairs)
 result["original_2023_response_mean_per_HUC"]={
  "novel_matched":float(target.groupby("huc8").ept_fraction.mean().reindex([p[0] for p in pairs]).mean()),
  "seen_matched":float(target.groupby("huc8").ept_fraction.mean().reindex([p[1] for p in pairs]).mean())}
 result["matched_env_SMD"]=numeric_match_profiles(set(old.huc8),huc_target,pairs)
 report={
  "schema":"gesui_NRSA2023_HUC8_matched_ecoregion_environment_fixedEPTprediction_v1",
  "original_SOURCE_pair_digest":pair_sha,"source_matched_GO":True,
  "source_matched_n_HUC8":len(pairs),"no_training_overlap_any_2023_HUC8":True,
  "old_EPA_original_source_n_2648":True,"old_train_n_HUC8":551,
  "full_2023_HUC8":693,"source_EPT_prior_full_model_MSE_reproduced":True,
  "original_EPA_SHA":ref["2023_source_original_SHA"],
  "matched_results":result,
  "posthoc_ecological_support_analysis":True,
  "no_WWTP_source_or_ecological_mechanism_identified":True}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("NRSA_NOVELTY_ENVIRONMENT_MATCHED_EPT",json.dumps({
  "pairs":len(pairs),"source_pair_sha":pair_sha,
  "matched_smd":result["matched_env_SMD"],
  "outcome_mean":result["original_2023_response_mean_per_HUC"],
  "contrasts":result["comparisons"]
  },sort_keys=True),flush=True)
if __name__=="__main__":main()
