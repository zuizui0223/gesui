#!/usr/bin/env python3
"""Independent NumPy original EPA fourth-cycle source-matched EPT audit.

Recalculates model predictions with weighted normal equations, not sklearn.
Recalculates per-HUC8 residuals, pair comparison and same-pair bootstraps
without importing the first matched-EPT result computation functions.
"""
import argparse,collections,hashlib,json,pathlib
import numpy as np
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as previous
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as domain
from audit_NRSA_fourthcycle_2023_manual_numpy_novelty import manual_fit_predict
CHECK=[
 ("PHab_after_NP","B+NP","B+NP+PHab"),
 ("bed_after_NP_other5","B+NP+PHab_except_bed","B+NP+PHab"),
 ("bed_alone_after_NP","B+NP","B+NP+bed"),
 ("NP_after_PHab","B+PHab","B+NP+PHab")]
def grouped_squared_error(df,yhat):
 d=collections.defaultdict(lambda:[0.,0])
 y=df.ept_fraction.to_numpy(float)
 for h,real,p in zip(df.huc8.to_numpy(),y,yhat):
  d[str(h)][0]+=(real-p)**2
  d[str(h)][1]+=1
 return {h:v[0]/v[1] for h,v in d.items()}
def bootstrap_gains(reference_n,full_n,reference_o,full_o,name):
 n=len(reference_n)
 seed=int(hashlib.sha256(name.encode()).hexdigest()[:8],16)
 rng=np.random.default_rng(seed)
 ix=rng.integers(0,n,size=(2000,n))
 n0=reference_n[ix].mean(1);n1=full_n[ix].mean(1)
 o0=reference_o[ix].mean(1);o1=full_o[ix].mean(1)
 gn=100*(n0-n1)/n0;go=100*(o0-o1)/o0
 diff=gn-go
 abs_diff=(reference_n-full_n)-(reference_o-full_o)
 bootstrap_abs=abs_diff[ix].mean(1)
 return {
  "novel_gain":float(100*(reference_n.mean()-full_n.mean())/reference_n.mean()),
  "seen_gain":float(100*(reference_o.mean()-full_o.mean())/reference_o.mean()),
  "novel_boot95":[float(x) for x in np.quantile(gn,[.025,.975])],
  "seen_boot95":[float(x) for x in np.quantile(go,[.025,.975])],
  "diff_gain":float(gn.mean()*0+100*(reference_n.mean()-full_n.mean())/reference_n.mean()-
                    100*(reference_o.mean()-full_o.mean())/reference_o.mean()),
  "diff_boot95":[float(x) for x in np.quantile(diff,[.025,.975])],
  "diff_abs":float(abs_diff.mean()),
  "diff_abs_boot95":[float(x) for x in np.quantile(bootstrap_abs,[.025,.975])]}
def check_equal(actual,expected,label,tol=1e-10):
 if not np.allclose(actual,expected,rtol=tol,atol=tol):
  raise RuntimeError("Independent matched EPA residual verification FAILED "+label+" "+str((actual,expected)))
def main():
 ap=argparse.ArgumentParser()
 for k in ("sites","water","habitat","benthos","chla","old_panel",
           "source_qa","novelty_source","matched_original","out"):
  ap.add_argument("--"+k,type=pathlib.Path,required=True)
 args=ap.parse_args()
 ref=json.loads(args.matched_original.read_text())
 src=json.loads(args.novelty_source.read_text())
 if not ref["source_matched_GO"] or not src["within_ecoregion_source_only_match"]["source_go"]:
  raise ValueError("Original matched source GO missing")
 old=previous.history(args.old_panel)
 env=previous.covariates(args,{})
 th=domain.huc_aggregate(env)
 tr=old.loc[~old.huc8.isin(set(th.index))].copy()
 trh=domain.huc_aggregate(tr)
 ever=th.index.isin(set(old.huc8))
 pairs,counts,digest=domain.choose_pairs(trh,th,ever)
 if digest!=ref["original_SOURCE_pair_digest"] or digest!=src["within_ecoregion_source_only_match"]["source_pair_membership_and_distance_SHA256"]:
  raise RuntimeError("Independent matching HUC8 source digest mismatch")
 original_src=json.loads(args.source_qa.read_text())
 t=previous.canonical_target(previous.make_target(args,original_src),original_src,"full")
 if len(t)!=1165 or t.huc8.nunique()!=693 or len(tr)!=1185 or tr.huc8.nunique()!=551:
  raise RuntimeError("Source counts no longer frozen")
 if set(tr.huc8)&set(t.huc8):raise RuntimeError("HUC8 data leakage")
 cols={q for _,a,b in CHECK for q in (a,b)}
 losses={key:grouped_squared_error(t,manual_fit_predict(tr,t,previous.MODEL[key]))
         for key in sorted(cols)}
 result={}
 for name,base,full in CHECK:
  v=ref["matched_results"]["comparisons"][name]
  n0=np.array([losses[base][p[0]] for p in pairs])
  n1=np.array([losses[full][p[0]] for p in pairs])
  o0=np.array([losses[base][p[1]] for p in pairs])
  o1=np.array([losses[full][p[1]] for p in pairs])
  fresh=bootstrap_gains(n0,n1,o0,o1,name)
  checks=[
   ("novel_gain","novel_gain_pct"),("seen_gain","previously_seen_source_gain_pct"),
   ("novel_boot95","novel_bootstrap_pair95"),
   ("seen_boot95","previously_seen_source_bootstrap_pair95"),
   ("diff_gain","novel_minus_seen_gain_pct_points"),
   ("diff_boot95","matched_pair_bootstrap95_gain_difference"),
   ("diff_abs","novel_minus_seen_absolute_MSE_reduction"),
   ("diff_abs_boot95","matched_pair_bootstrap95_absolute_MSE_difference")]
  for k,k2 in checks:check_equal(fresh[k],v[k2],name+"/"+k)
  result[name]={"posthoc_numeric_independent_replay":True,
   "novel_gain_pct":fresh["novel_gain"],"seen_gain_pct":fresh["seen_gain"],
   "difference_pctpoints":fresh["diff_gain"],
   "difference_boot95":fresh["diff_boot95"]}
 out={"schema":"NRSA2023_source_matched_HUC8_EPT_manual_NumPy_QA_v1",
  "original_novelty_matched_pairs":len(pairs),
  "exact_source_preoutcome_SHA":digest,
  "original_three_cycle_train_HUC8":551,
  "original_target_HUC8":693,
  "independent_weighted_10alpha_ridge_solver":True,
  "independent_grouped_HUC8_error_and_pair_bootstrap":True,
  "all_four_matched_biological_comparisons_verified":True,
  "posthoc_matched_source_observational_not_randomized":True,
  "results":result}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print("NRSA_NOVEL_HUC_MATCHED_INDEPENDENT_NUMPY_QA",json.dumps(
  {"match_source_hash":digest,"n_pairs":len(pairs),
   "results":result},sort_keys=True),flush=True)
if __name__=="__main__":main()
