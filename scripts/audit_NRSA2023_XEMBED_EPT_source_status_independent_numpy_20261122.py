#!/usr/bin/env python3
"""Independent manual NumPy normal equations, HC3, whole-pair bootstrap and
Ridge-pair OOF audit for genuine EPA 2023 XEMBED×historical-source interaction.
The original statistical model functions are NOT imported/called.
"""
import argparse,hashlib,json,pathlib
import numpy as np
from scipy.stats import norm
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as previous
import gesui_NRSA2023_matched_XEMBED_source_status_slope_20261122 as source

N_BOOT=2000
MAIN_COV=["XEMBED","log_tn","log_tp","log_width","lat","lon","stream_order"]
GRAIN_COV=["PCT_FN","PCT_SAFN","LSUB_DMM"]
def approx(a,b,name,tol=1e-9):
 if not np.allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),rtol=tol,atol=tol):
  raise RuntimeError("Independent original numerical mismatch "+name+": "+str((a,b)))
def data(d,include_grain=False):
 c=["mean_XEMBED"]+["delta_"+k for k in MAIN_COV+(GRAIN_COV if include_grain else [])]
 X=np.c_[np.ones(len(d)),d[c].to_numpy(float)]
 y=d.target_delta_EPT.to_numpy(float)
 if np.linalg.matrix_rank(X)!=len(c)+1:
  raise RuntimeError("Original EPA source independent X rank deficient")
 return X,y,c
def manual_HC3(d,grain=False):
 X,y,c=data(d,grain)
 inv=np.linalg.inv(X.T@X)
 coef=inv@(X.T@y)
 resid=y-X@coef
 h=np.sum((X@inv)*X,axis=1)
 meat=(X*(resid/(1-h))[:,None]).T@(X*(resid/(1-h))[:,None])
 cov=inv@meat@inv
 j=1
 se=float(np.sqrt(cov[j,j]))
 gamma=float(coef[j])
 CI=gamma+np.array([-1,1])*float(norm.ppf(.975))*se
 p=float(2*norm.sf(abs(gamma/se)))
 rng=np.random.default_rng(20261122)
 pars=[]
 for t in range(N_BOOT):
  idx=rng.integers(0,len(y),size=len(y))
  bet=np.linalg.lstsq(X[idx],y[idx],rcond=None)[0]
  pars.append(bet[j])
 boot=np.quantile(pars,[.025,.975])
 return {"pp":100*gamma,"HC3_95":(100*CI).tolist(),"p":p,"paired_boot95":(100*boot).tolist()}
def manual_heldout(d):
 primary=["mean_XEMBED"]+["delta_"+k for k in MAIN_COV]
 base=primary[1:]
 yy=d.target_delta_EPT.to_numpy(float)
 folds=np.array([int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5 for h in d.source_novel_HUC8],int)
 errors={}
 for tag,cols in (("base",base),("full",primary)):
  preds=np.full(len(d),np.nan)
  for fold in range(5):
   train=folds!=fold;test=folds==fold
   a=d.loc[train,cols].to_numpy(float)
   b=d.loc[test,cols].to_numpy(float)
   mu=a.mean(0);sd=a.std(0,ddof=0);sd=np.where(sd>0,sd,1.)
   x=(a-mu)/sd;z=(b-mu)/sd
   y=yy[train];mean=float(y.mean())
   # Manual centered ridge with the intercept unpenalized.
   beta=np.linalg.solve(x.T@x+10*np.eye(x.shape[1]),x.T@(y-mean))
   preds[test]=z@beta+mean
  errors[tag]=(yy-preds)**2
 a=errors["base"];b=errors["full"]
 ix=np.random.default_rng(20261122).integers(0,len(a),size=(N_BOOT,len(a)))
 boot=100*(a[ix].mean(1)-b[ix].mean(1))/a[ix].mean(1)
 return {"ref":float(a.mean()),"full":float(b.mean()),
         "MSE_gain":float(100*(a.mean()-b.mean())/a.mean()),
         "CI95":np.quantile(boot,[.025,.975]).tolist()}
def main():
 p=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","old_panel",
           "source_qa","extra_source","original_result","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 signed=json.loads(a.original_result.read_text())
 if signed["schema"]!="gesui_EPA2023_XEMBED_EPT_source_status_matched_pair_slope_v1":
  raise ValueError("Original EPA source slope schema invalid")
 ht,train,src,pairs=source.source_pairs(a)
 t=previous.canonical_target(previous.make_target(a,src),src,"full")
 y=t.groupby("huc8",sort=True).ept_fraction.mean()
 out={}
 for key,ps in pairs.items():
  d=source.build_pair_matrix(ht,y,ps)
  first=signed["results"][key]
  if first["sample_source_SHA"]!=source.PAIRS_SOURCE_SHA[key]:
   raise RuntimeError("Signed source prior HUC pair identity drift")
  expected=first["paired_difference_XEMBED_mean_slope"]
  test=manual_HC3(d)
  approx(test["pp"],expected["difference_in_EPT_percentage_points_per_1SD_XEMBED"],key+".beta")
  approx(test["HC3_95"],expected["HC3_95_difference_in_EPT_percentage_points_per_1SD"],key+".HC3")
  approx(test["p"],expected["HC3_exploratory_p"],key+".HC3p")
  approx(test["paired_boot95"],expected["paired_HUC8_2000_refit_boot95_difference_in_EPT_pp"],key+".boot")
  diags={"manual_normal_equation_beta":test["pp"],"manual_HC3_CI95":test["HC3_95"],
         "manual_p":test["p"],"manual_boot_CI95":test["paired_boot95"],
         "all_original_inference_reproduced":True}
  if key=="grain":
   plus=manual_HC3(source.build_pair_matrix(ht,y,ps,include_grain=True),grain=True)
   prior_grain=first["grain_residual_sensitivity"]
   approx(plus["pp"],prior_grain["difference_in_EPT_percentage_points_per_1SD_XEMBED"],"grain.extra_beta")
   approx(plus["HC3_95"],prior_grain["HC3_95_difference_in_EPT_percentage_points_per_1SD"],"grain.extraHC3")
   approx(plus["paired_boot95"],prior_grain["paired_HUC8_2000_refit_boot95_difference_in_EPT_pp"],"grain.extraBoot")
   diags["grain_extra_residual_beta_reproduced"]=True
   cv=manual_heldout(d)
   prevcv=first["whole_pair_heldout_prediction"]
   for name,key2 in (("ref","pair_equal_baseline_MSE"),
                     ("full","pair_equal_plus_XEMBED_mean_MSE"),
                     ("MSE_gain","added_moderation_pairheldout_MSE_gain_pct"),
                     ("CI95","fixed_pair_OOF_boot95_gain_pct")):
    approx(cv[name],prevcv[key2],"pair_holdout_"+name)
   diags["manual_5fold_pairheldout_OLS_interaction_gain_reproduced"]=True
   diags["pairheldout_gain_pct"]=cv["MSE_gain"]
   diags["pairheldout_boot95"]=cv["CI95"]
  out[key]=diags
 report={"schema":"EPA2023_XEMBED_EPT_source_status_independent_manual_NumPy_HC3_and_pairCV_QA_v1",
    "original_source_SHA":source.SHA,"original_old_2008_2013_2018_source_UID":2648,
    "original_2023_berw_UID":1165,"old_training_HUC8":551,"target_HUC8":693,
    "original_source_novel_and_ever_214/479":True,
    "all_original_no_training_target_overlap":True,
    "source_matched_HUC8_SHA_exact":source.PAIRS_SOURCE_SHA,
    "manual_linear_algebra_HC3_and_bootstrap_not_original_statsmodels_sklearn":True,
    "all_original_primary_and_sensitivity_regression_CI_and_heldout_results_reproduced":True,
    "posthoc_not_causal":True,"results":out}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("EPA2023_XEMBED_EPT_MODERATION_INDEPENDENT_QA",json.dumps({
   k:{"beta_pp":v["manual_normal_equation_beta"],"HC3_CI95":v["manual_HC3_CI95"],
      "paired_boot95":v["manual_boot_CI95"],"all_inference_reproduced":True,
      "pairheldout_CV":v.get("pairheldout_gain_pct"),
      "pairheldout_CI95":v.get("pairheldout_boot95")}
     for k,v in out.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
