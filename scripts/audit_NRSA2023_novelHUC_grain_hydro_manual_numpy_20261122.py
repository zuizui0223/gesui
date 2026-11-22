#!/usr/bin/env python3
"""Independent EPA NARS fourth-era matched-source numerical QA.

Compared with sklearn Ridge original, independently refit previous 3 ERA
HUC8-equal original EPT models via weighted NumPy normal equations,
recompute exact within-HUC errors and pair bootstrap fixed predictions.
"""
import argparse,hashlib,json,pathlib
import numpy as np
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as prev
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as dom
import gesui_nrsa2023_huc_novelty_extra_source_matching_20261122 as ext
from audit_NRSA_fourthcycle_2023_manual_numpy_novelty import manual_fit_predict
from audit_NRSA2023_HUC8_novelty_matched_ManualNumPy_20261122 import (
 grouped_squared_error,bootstrap_gains,check_equal,CHECK)
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

def main():
 ap=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","old_panel",
           "source_qa","extension_source","original_matched","out"):
  ap.add_argument("--"+k,type=pathlib.Path,required=True)
 a=ap.parse_args()
 signed=json.loads(a.extension_source.read_text())
 initial=json.loads(a.original_matched.read_text())
 if initial["schema"]!="NRSA2023_original_HUC8_novelty_extra_substrate_hydro_matched_EPT_posthoc_v1":
  raise ValueError("Signed new matching result schema mismatch")
 if signed["schema"]!="NRSA2023_original_four_cycle_HUC_novelty_extra_source_domains_v1":
  raise ValueError("Signed source-only schema mismatch")
 for k in ("sites","water","habitat","chla","benthos"):
  if sha(getattr(a,k))!=SHA[k]:raise RuntimeError("EPA source original SHA changed "+k)
 old=prev.history(a.old_panel)
 env=prev.covariates(a,{})
 h=dom.huc_aggregate(env)
 grp=ext.block_sources(a,env).join(env[["huc8"]]).groupby("huc8").mean()
 h=h.join(grp[[k for q in ext.BLOCKS.values() for k in q]],how="left")
 ever=h.index.isin(set(old.huc8))
 tr=old.loc[~old.huc8.isin(set(h.index))].copy()
 trainh=dom.huc_aggregate(tr)
 if len(tr)!=1185 or trainh.shape[0]!=551 or h.shape[0]!=693:
  raise RuntimeError("Frozen original 2023 / prior source membership changed")
 source_qa=json.loads(a.source_qa.read_text())
 target=prev.canonical_target(prev.make_target(a,source_qa),source_qa,"full")
 if len(target)!=1165 or target.huc8.nunique()!=693:
  raise RuntimeError("Original EPT 2023 source 1165/693 changed")
 oldmodels={k for _,z,y in CHECK for k in (z,y)}
 errs={key:grouped_squared_error(target,manual_fit_predict(tr,target,prev.MODEL[key]))
       for key in sorted(oldmodels)}
 decisions={}
 for block,cols in ext.BLOCKS.items():
  source=signed["source_only_blocks"][block]
  actual=initial["source_original_new_extra_match_blocks"][block]
  if not source["source_go"]:
   if actual["status"]!="SOURCE_STOP_DO_NOT_FIT_OR_REPORT_OUTCOME":
    raise RuntimeError("Source STOP changed to numerical claim")
   decisions[block]={"status":"ORIGINAL_SOURCE_STOP_REPRODUCED"};continue
  pairs,_,digest=ext.choose(h,trainh,ever,ext.BASE+cols)
  if digest!=source["source_match_pair_membership_SHA256"] or len(pairs)!=source["source_only_match_pairs"]:
   raise RuntimeError("Original EPA source-only match pair identity changed "+block)
  if digest!=actual["source_match_digest"]:
   raise RuntimeError("Signed numeric match differs from source-only")
  comparisons={}
  for label,base,full in CHECK:
   z=actual["HUC8_equal_matched_loss_contrasts"]["comparisons"][label]
   n0=np.array([errs[base][p[0]] for p in pairs]);n1=np.array([errs[full][p[0]] for p in pairs])
   o0=np.array([errs[base][p[1]] for p in pairs]);o1=np.array([errs[full][p[1]] for p in pairs])
   now=bootstrap_gains(n0,n1,o0,o1,label)
   compare=[
    ("novel_gain","novel_gain_pct"),("seen_gain","previously_seen_source_gain_pct"),
    ("novel_boot95","novel_bootstrap_pair95"),("seen_boot95","previously_seen_source_bootstrap_pair95"),
    ("diff_gain","novel_minus_seen_gain_pct_points"),
    ("diff_boot95","matched_pair_bootstrap95_gain_difference"),
    ("diff_abs","novel_minus_seen_absolute_MSE_reduction"),
    ("diff_abs_boot95","matched_pair_bootstrap95_absolute_MSE_difference")]
   for akey,rkey in compare:check_equal(now[akey],z[rkey],block+"."+label+"."+akey)
   comparisons[label]={"independent_replay":True,
     "seen_gain":now["seen_gain"],"novel_gain":now["novel_gain"],
     "difference":now["diff_gain"],"difference95":now["diff_boot95"]}
  decisions[block]={"status":"ALL_FOUR_ORIGINAL_MATCHED_NUMERICAL_RESULTS_INDEPENDENTLY_REPLAYED",
    "signed_source_matched_pairs":len(pairs),"source_pair_digest":digest,
    "comparisons":comparisons}
 report={"schema":"EPA2023_original_grain_hydrology_novelHUC8_matches_independent_NumPyQA_v1",
  "original_new_source_SHA256":hashlib.sha256(a.extension_source.read_bytes()).hexdigest(),
  "original_new_numeric_SHA256":hashlib.sha256(a.original_matched.read_bytes()).hexdigest(),
  "original_previous_3_training_UID":1185,"training_HUC8":551,"new_2023_target_UID":1165,
  "zero_same_HUC8_old_training_new_target":True,
  "independent_NumPy_Ridge_not_original_sklearn":True,
  "independent_HUC8_pair_bootstrap_recomputed":True,
  "posthoc_not_causal_or_sewage_mechanism":True,
  "block_verdicts":decisions}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("NRSA2023_NEW_GRAIN_HYDRO_MATCHED_NUMPY_INDEPENDENT_PASS",json.dumps(
  {k:{"status":v["status"],"pairs":v.get("signed_source_matched_pairs"),
    "bed_alone":v.get("comparisons",{}).get("bed_alone_after_NP")}
    for k,v in decisions.items()},sort_keys=True))
if __name__=="__main__":main()
