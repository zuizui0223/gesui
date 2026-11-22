#!/usr/bin/env python3
"""Independent weighted Ridge normal equations audit of EPA fourth-era EPT.

The original source joining is intentionally reused from separately source-only
audited original EPA files, but ALL model fits, normalization, category encodings,
prediction and bootstrap calculations are independently coded with NumPy.
"""
import argparse,hashlib,json,pathlib
import numpy as np,pandas as pd
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as original
import gesui_nrsa2324_fourthcycle_source_only_20261122 as QA

def manual_fit_predict(train,target,columns):
 x_train=[];x_test=[]
 for k in columns:
  tr=train[k].to_numpy(float);te=target[k].to_numpy(float)
  mean=float(tr.mean());sd=float(tr.std(ddof=0))
  if not sd>0:sd=1.
  x_train.append(((tr-mean)/sd)[:,None]);x_test.append(((te-mean)/sd)[:,None])
 for cat in original.CATS:
  t=train[cat].astype(str).to_numpy()
  v=target[cat].astype(str).to_numpy()
  levels=sorted(set(t.tolist()))
  x_train.append((t[:,None]==np.asarray(levels)[None,:]).astype(float))
  x_test.append((v[:,None]==np.asarray(levels)[None,:]).astype(float))
 X=np.hstack(x_train);Z=np.hstack(x_test)
 if not (np.isfinite(X).all() and np.isfinite(Z).all()):
  raise RuntimeError("Invalid independent standardized source features")
 y=train.ept_fraction.to_numpy(float)
 counts=train.huc8.value_counts()
 n=len(train);g=len(counts)
 w=((n/g)/train.huc8.map(counts)).to_numpy(float)
 if abs(w.mean()-1)>1e-9:raise RuntimeError("Catchment-equal weights not normalized")
 meanx=np.average(X,axis=0,weights=w);meany=np.average(y,weights=w)
 xx=(X-meanx)*np.sqrt(w)[:,None]
 yy=(y-meany)*np.sqrt(w)
 gram=xx.T@xx+10*np.eye(X.shape[1])
 rhs=xx.T@yy
 coef=np.linalg.solve(gram,rhs)
 result=(Z-meanx)@coef+meany
 if not np.isfinite(result).all():raise RuntimeError("Nonfinite independent ridge prediction")
 return result

def group_loss(target,pred):
 err=(target.ept_fraction.to_numpy(float)-pred)**2
 return pd.Series(err,index=target.huc8.to_numpy()).groupby(level=0).mean().sort_index()

def ci_gain(groups0,groups1,label):
 a=groups0.to_numpy();b=groups1.reindex(groups0.index).to_numpy()
 seed=int(hashlib.sha256(label.encode()).hexdigest()[:8],16)
 ids=np.random.default_rng(seed).integers(0,len(a),size=(2000,len(a)))
 boot=100*(a[ids].mean(axis=1)-b[ids].mean(axis=1))/a[ids].mean(axis=1)
 return {"HUC8_gain_pct":float(100*(a.mean()-b.mean())/a.mean()),
  "fixed_OOF_boot95_pct":[float(x) for x in np.quantile(boot,[.025,.975])]}
def same(a,b,desc,tol=1e-8):
 if not np.allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),
                    rtol=tol,atol=tol,equal_nan=False):
  raise RuntimeError("Independent solver mismatch "+desc+": "+str((a,b)))
def check_one(test,old,case,signed):
 tr=old[~old.huc8.isin(set(test.huc8))].reset_index(drop=True).copy()
 assert not(set(tr.huc8)&set(test.huc8))
 assert len(tr)==signed["source"]["train_after_removing_every_shared_target_HUC8_UID"]
 assert tr.huc8.nunique()==signed["source"]["train_HUC8"]
 pred={key:manual_fit_predict(tr,test,features)
       for key,features in original.MODEL.items()}
 losses={key:group_loss(test,values) for key,values in pred.items()}
 worst=0.
 for name,loss in losses.items():
  err=float(loss.mean())
  reference=signed["target"]["model_errors"][name]["MSE_equal_HUC8"]
  worst=max(worst,abs(err-reference))
  same(err,reference,case+"/"+name+" original heldout MSE")
 main_comparison={}
 for label,base,full in original.COMPARE:
  v=ci_gain(losses[base],losses[full],label)
  expected=signed["target"]["comparisons"][label]
  same(v["HUC8_gain_pct"],expected["HUC8_equal_gain_pct"],case+"/"+label+" gain")
  same(v["fixed_OOF_boot95_pct"],expected["fixed_target_HUC8_bootstrap95"],case+"/"+label+" bootstrap")
  main_comparison[label]=v
 # HUC8 defined before opening the split's responses: whether ever appeared in
 # any original old-era SOURCE, not whether model trained on that catchment.
 union=set(old.huc8)
 partitions={"old_source_ever_sampled_but_removed_from_train":test.huc8.isin(union),
             "never_sampled_old_source_HUC8":~test.huc8.isin(union)}
 split_results={}
 for name,mask in partitions.items():
  target=test.loc[mask].reset_index(drop=True)
  if len(target)<150 or target.huc8.nunique()<100:
   split_results[name]={"status":"SOURCE_STOP","UID":len(target),"HUC8":target.huc8.nunique()}
   continue
  m=mask.to_numpy()
  err_groups={key:group_loss(target,pp[m]) for key,pp in pred.items()}
  summary={}
  for label,b,f in original.COMPARE:
   summary[label]=ci_gain(err_groups[b],err_groups[f],label)
  split_results[name]={"status":"POSTHOC_CATCHMENT_NOVELTY_PARTITION",
      "original_UID":len(target),"HUC8":target.huc8.nunique(),
      "physical_sites":target.site_id.nunique(),
      "eco9_distribution":target.eco9.value_counts().to_dict(),
      "EPS_original_response_subsample":"EPT counted fraction",
      "comparisons":summary}
 # Target-only previously signed result already checked via numeric aggregate.
 check_novel=signed["novel_target_only_HUC8"]
 if check_novel["status"]=="NOVEL_TARGET_HUC8_COMPLETE":
  n=split_results["never_sampled_old_source_HUC8"]
  assert n["original_UID"]==check_novel["results"]["source_test_UID"]
  assert n["HUC8"]==check_novel["results"]["source_test_HUC8"]
  for label,v in n["comparisons"].items():
   wanted=check_novel["results"]["comparisons"][label]
   same(v["HUC8_gain_pct"],wanted["HUC8_equal_gain_pct"],case+"/"+label+" novel")
   same(v["fixed_OOF_boot95_pct"],wanted["fixed_target_HUC8_bootstrap95"],case+"/"+label+" novel bootstrap")
 return {"source_target_UID":len(test),"source_target_HUC8":test.huc8.nunique(),
    "source_training_UID":len(tr),"source_training_HUC8":tr.huc8.nunique(),
    "maximum_absolute_MSE_error_against_original_Ridge":worst,
    "all_ten_original_model_MSEs_and_eight_bootstrap_comparisons_replayed":True,
    "main_full_target":main_comparison,
    "catchment_novelty":split_results}
def main():
 p=argparse.ArgumentParser()
 for k in QA.SHA:p.add_argument("--"+k,required=True,type=pathlib.Path)
 p.add_argument("--source-qa",type=pathlib.Path,required=True)
 p.add_argument("--old-panel",type=pathlib.Path,required=True)
 p.add_argument("--original-result",type=pathlib.Path,required=True)
 p.add_argument("--out",type=pathlib.Path,required=True)
 a=p.parse_args()
 source=json.loads(a.source_qa.read_text())
 ref=json.loads(a.original_result.read_text())
 if ref["schema"]!="gesui_original_EPA_4thcycle_2023_24_EPT_NP_PHab_outofera_HUC8_holdout_v1":
  raise ValueError("Original 2023 actual signed result source schema mismatch")
 if ref["2023_source_original_SHA"]!=QA.SHA or source["original_EPA_sha256"]!=QA.SHA:
  raise ValueError("EPA original 2023 original SHA set conflict")
 for k,pin in QA.SHA.items():
  if QA.sha(getattr(a,k))!=pin:raise RuntimeError("Official EPA source digest mismatch "+k)
 old=original.history(a.old_panel)
 target=original.make_target(a,source)
 result={}
 for case in ("full","periphyton"):
  t=original.canonical_target(target,source,case)
  result[case]=check_one(t,old,case,ref["results"][case])
 report={"schema":"gesui_NR​​SA_fourthcycle_independent_weighted_numpy_and_novelHUC_QA_v1",
   "original_EPA_source_SHA":QA.SHA,
   "original_3cycle_source_expected":original.OLD_EXPECT,
   "all_model_original_numeric_and_CI_independent_replay_pass":True,
   "no_target_year_used_in_training":True,
   "no_shared_train_and_test_HUC8":True,
   "posthoc_novelty_partition_not_further_independent_validation":True,
   "not_WWTP_or_cross_ecosystem_transfer":True,"results":result}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
 print("FOURCYCLE_INDEPENDENT_NUMPY_QA",json.dumps({
  "both_target_source_cohorts_verified":True,
  "full_target_MSE_error":result["full"]["maximum_absolute_MSE_error_against_original_Ridge"],
  "source_novelty_splits":{case:{k:{"UID":v.get("original_UID"),"HUC8":v.get("HUC8"),
    "PHab_gain":v.get("comparisons",{}).get("PHab_after_NP",{}).get("HUC8_gain_pct"),
    "bed_gain":v.get("comparisons",{}).get("bed_after_NP_other5",{}).get("HUC8_gain_pct")}
   for k,v in res["catchment_novelty"].items()} for case,res in result.items()}},sort_keys=True),flush=True)
if __name__=="__main__":main()
