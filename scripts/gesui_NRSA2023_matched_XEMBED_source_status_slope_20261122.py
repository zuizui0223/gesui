#!/usr/bin/env python3
"""2023 EPA source-only HUC8 pairs: one EPT × XEMBED status-slope falsifier.

Outcome is fixed original EPA counted EPT fraction, averaged by HUC8.
Source grouping means ever / never in prior three *eligible source samples*,
NOT a causal intervention or evidence river never sampled anywhere.
"""
import argparse,hashlib,json,pathlib
import numpy as np,pandas as pd,statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
import gesui_NRSA_fourthcycle_2023_EPT_heldout_transport as prior
import gesui_nrsa2023_HUC8_source_novelty_environment_20261122 as domain
import gesui_nrsa2023_huc_novelty_extra_source_matching_20261122 as extra
from gesui_nrsa2324_fourthcycle_source_only_20261122 import SHA,sha

PAIRS_SOURCE_SHA={
 "grain":"8c78c9060e0d512666062760f0467baaa54f6cfabeeaf76867ad609e3f1ea8e2",
 "hydrology_landuse":"4d21361aa86707039d4afd62cfa11c70f629c50ea01e1a0663e4c00930f254b2"}
DIFF_ORIGINAL=["XEMBED","log_tn","log_tp","log_width","lat","lon","stream_order"]
EXTRA_GRAIN=["PCT_FN","PCT_SAFN","LSUB_DMM"]
N_BOOT=2000
def source_pairs(args):
 src=json.loads(args.source_qa.read_text())
 if src["original_EPA_sha256"]!=SHA:raise RuntimeError("Source-locked original 2023 EPA SHA changed")
 frozen=json.loads(args.extra_source.read_text())
 if frozen["schema"]!="NRSA2023_original_four_cycle_HUC_novelty_extra_source_domains_v1":
  raise RuntimeError("Source-only frozen grain/hydro match contract changed")
 for k in SHA:
  if sha(getattr(args,k))!=SHA[k]:raise RuntimeError("EPA original file SHA mismatch "+k)
 history=prior.history(args.old_panel)
 env=prior.covariates(args,src)
 if len(history)!=2648 or history.huc8.nunique()!=1030 or len(env)!=1165 or env.huc8.nunique()!=693:
  raise RuntimeError("Original old and target source cohorts altered")
 ht=domain.huc_aggregate(env)
 oldtrain=history.loc[~history.huc8.isin(set(ht.index))].copy()
 oldh=domain.huc_aggregate(oldtrain)
 if len(oldtrain)!=1185 or oldh.shape[0]!=551 or set(oldh.index)&set(ht.index):
  raise RuntimeError("Original training catchment source not isolated from target")
 supplements=extra.block_sources(args,env).join(env[["huc8"]]).groupby("huc8").mean()
 ht=ht.join(supplements[[x for cols in extra.BLOCKS.values() for x in cols]],how="left")
 old_status=ht.index.isin(set(history.huc8))
 if (old_status.sum(),(~old_status).sum())!=(479,214):
  raise RuntimeError("Old-source ever/never target group membership changed")
 pairs={}
 for k,block in extra.BLOCKS.items():
  s=frozen["source_only_blocks"][k]
  pp,eco,hashval=extra.choose(ht,oldh,old_status,extra.BASE+block)
  if (hashval!=PAIRS_SOURCE_SHA[k] or hashval!=s["source_match_pair_membership_SHA256"] or
      len(pp)!=s["source_only_match_pairs"] or not s["source_go"]):
   raise RuntimeError("Source-only original signed paired HUC contract failed "+k)
  pairs[k]=pp
 return ht,oldtrain,src,pairs
def build_pair_matrix(ht,y,pairs,include_grain=False):
 n=np.array([p[0] for p in pairs]);o=np.array([p[1] for p in pairs])
 if len(np.unique(n))!=len(n) or len(np.unique(o))!=len(o) or set(n)&set(o):
  raise RuntimeError("HUC8 pair re-used, not independent")
 if (set(n)|set(o)) - set(ht.index):raise RuntimeError("Source pairing missing original HUC8")
 cols=DIFF_ORIGINAL+(EXTRA_GRAIN if include_grain else [])
 z={}
 for name in cols:
  center=float(ht[name].mean())
  sd=float(ht[name].std(ddof=0))
  if not np.isfinite(sd) or sd<1e-8:raise RuntimeError("Source-only predictor has no independent variation "+name)
  a=ht.loc[n,name].to_numpy(float);b=ht.loc[o,name].to_numpy(float)
  z["delta_"+name]=(a-b)/sd
  if name=="XEMBED":
   z["mean_XEMBED"]=(0.5*(a+b)-center)/sd
   z["mean_XEMBED_unscaled"]=0.5*(a+b)
   z["XEMBED_sd_2023_source"]=sd
 dx=pd.DataFrame(z)
 dx["target_delta_EPT"]=y.reindex(n).to_numpy(float)-y.reindex(o).to_numpy(float)
 dx["source_novel_HUC8"]=n
 dx["source_prior_HUC8"]=o
 dx["source_eco9"]=[p[2] for p in pairs]
 if not np.isfinite(dx.select_dtypes("number").to_numpy(dtype=float)).all():
  raise RuntimeError("Nonfinite original source matched EPT or covariates")
 return dx
def fit(d,more_grain=False):
 columns=["mean_XEMBED"]+["delta_"+k for k in DIFF_ORIGINAL+(EXTRA_GRAIN if more_grain else [])]
 X=sm.add_constant(d[columns].astype(float),has_constant="add")
 y=d.target_delta_EPT.to_numpy(float)
 if np.linalg.matrix_rank(X.to_numpy(float))!=X.shape[1]:
  return {"status":"RANK_STOP","matrix_rank":int(np.linalg.matrix_rank(X))}
 res=sm.OLS(y,X).fit(cov_type="HC3")
 p=float(res.pvalues["mean_XEMBED"])
 ci=res.conf_int().loc["mean_XEMBED"].to_numpy(float)
 gamma=float(res.params["mean_XEMBED"])
 # Repeat paired-HUC8 resampling without duplicating individual visits.
 seed=20261122
 rng=np.random.default_rng(seed)
 xv=X.to_numpy(float);pars=[]
 for i in range(N_BOOT):
  ix=rng.integers(0,len(y),size=len(y))
  coef=np.linalg.lstsq(xv[ix],y[ix],rcond=None)[0]
  pars.append(coef[columns.index("mean_XEMBED")+1])
 ci_boot=np.quantile(np.array(pars,dtype=float),[.025,.975])
 return {"status":"POSTHOC_PAIRED_GAMMA_COMPLETE",
   "n_matched_pairs":len(d),"n_pair_ecoregions":d.source_eco9.nunique(),
   "adjusted_width_main_controls":["delta_"+k for k in DIFF_ORIGINAL],
   "additional_grain_residual_differences":more_grain,
   "XEMBED_one_pooled2023_SD":float(d.XEMBED_sd_2023_source.iloc[0]),
   "mean_XEMBED_original_units_q10_q90":list(map(float,np.quantile(d.mean_XEMBED_unscaled,[.1,.9]))),
   "difference_in_EPT_fractions_per_1SD_XEMBED":gamma,
   "difference_in_EPT_percentage_points_per_1SD_XEMBED":100*gamma,
   "HC3_95_difference_in_EPT_percentage_points_per_1SD":[float(100*v) for v in ci],
   "HC3_exploratory_p":p,
   "paired_HUC8_2000_refit_boot95_difference_in_EPT_pp":[float(100*v) for v in ci_boot],
   "adjusted_R2_descriptive":float(res.rsquared_adj),
   "n_regressors_including_intercept":int(X.shape[1]),
   "not_causal_source_history_or_sediment":True}
def heldout(d):
 columns=["mean_XEMBED"]+["delta_"+k for k in DIFF_ORIGINAL]
 base=columns[1:]
 yy=d.target_delta_EPT.to_numpy(float)
 folds=np.array([int(hashlib.sha256(str(z).encode()).hexdigest()[:8],16)%5 for z in d.source_novel_HUC8])
 counts=[int((folds==f).sum()) for f in range(5)]
 if min(counts)<20:raise RuntimeError("Bad HUC8 matched-pair fold source")
 predictions={}
 for name,features in (("main_effects_only",base),("mean_XEMBED_moderation",columns)):
  pred=np.full(len(d),np.nan)
  for f in range(5):
   tr=folds!=f;te=folds==f
   scaler=StandardScaler()
   X=scaler.fit_transform(d.loc[tr,features].astype(float))
   V=scaler.transform(d.loc[te,features].astype(float))
   pred[te]=Ridge(alpha=10,solver="svd").fit(X,yy[tr]).predict(V)
  if not np.isfinite(pred).all():raise RuntimeError("Heldout full-pair fold missing predictions")
  predictions[name]=pred
 e0=(yy-predictions["main_effects_only"])**2
 e1=(yy-predictions["mean_XEMBED_moderation"])**2
 rng=np.random.default_rng(20261122)
 ix=rng.integers(0,len(d),size=(N_BOOT,len(d)))
 boots=100*(e0[ix].mean(1)-e1[ix].mean(1))/e0[ix].mean(1)
 return {"pair_huc8_disjoint_fivefold":True,"fixed_fold_pair_counts":counts,
   "pair_equal_baseline_MSE":float(e0.mean()),
   "pair_equal_plus_XEMBED_mean_MSE":float(e1.mean()),
   "added_moderation_pairheldout_MSE_gain_pct":float(100*(e0.mean()-e1.mean())/e0.mean()),
   "fixed_pair_OOF_boot95_gain_pct":[float(v) for v in np.quantile(boots,[.025,.975])],
   "prediction_target":"paired EPT fraction difference, not absolute EPT share",
   "not_independent_external_data":True}
def descriptive_slopes(ht,y,pairs):
 result={}
 for name,ids in (("source_novel",[p[0] for p in pairs]),("old_source_ever",[p[1] for p in pairs])):
  x=ht.loc[ids,"XEMBED"].to_numpy(float)
  yy=y.reindex(ids).to_numpy(float)
  sd=float(ht.XEMBED.std(ddof=0))
  xx=(x-float(ht.XEMBED.mean()))/sd
  X=sm.add_constant(xx)
  fit=sm.OLS(yy,X).fit(cov_type="HC3")
  ci=fit.conf_int()[1]
  result[name]={"unadjusted_slope_EPT_pp_per_1SD_XEMBED":float(100*fit.params[1]),
   "HC3_descriptive95_pp":[float(100*ci[0]),float(100*ci[1])],
   "no_adjustment_in_this_subgroup_diagnostic":True}
 return result
def main():
 p=argparse.ArgumentParser()
 for k in ("sites","water","habitat","chla","benthos","old_panel","source_qa","extra_source","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 ht,old,src,pairs=source_pairs(a)
 # Biological outcome opened only after original environmental source pair
 # identity and zero train/test HUC8 overlap have passed.
 target=prior.canonical_target(prior.make_target(a,src),src,"full")
 if len(target)!=1165 or target.huc8.nunique()!=693:
  raise RuntimeError("Original 2023 EPT response cohort changed")
 y=target.groupby("huc8",sort=True).ept_fraction.mean()
 if set(y.index)!=set(ht.index):raise RuntimeError("Source HUC and biological target HUC mismatch")
 result={}
 for key,pair in pairs.items():
  d=build_pair_matrix(ht,y,pair)
  ratio=(np.quantile(d.mean_XEMBED_unscaled,.9)-np.quantile(d.mean_XEMBED_unscaled,.1))/float(d.XEMBED_sd_2023_source.iloc[0])
  if len(d)<150 or d.source_eco9.nunique()<5 or ratio<1:
   result[key]={"status":"SOURCE_STOP","n_pairs":len(d),
     "source_eco9":d.source_eco9.nunique(),"q90_q10_bed_range_per_pooled_SD":float(ratio)}
   continue
  primary=fit(d)
  if primary["status"]!="POSTHOC_PAIRED_GAMMA_COMPLETE":
   result[key]={"status":"COLLINEARITY_STOP","primary":primary};continue
  if key=="grain": # single primary source-matched HUC cohort
   primary["source_frozen_primary_test"]=True
   results_cv=heldout(d)
   frag=fit(build_pair_matrix(ht,y,pair,include_grain=True),more_grain=True)
   rawslopes=descriptive_slopes(ht,y,pair)
  else:
   primary["source_frozen_primary_test"]=False
   results_cv=None;frag=None;rawslopes=None
  result[key]={"status":"POSTHOC_SOURCE_MODERATION_COMPLETE",
   "sample_source_SHA":PAIRS_SOURCE_SHA[key],
   "q90_q10_bed_range_per_pooled_SD":float(ratio),
   "paired_difference_XEMBED_mean_slope":primary,
   "whole_pair_heldout_prediction":results_cv,
   "grain_residual_sensitivity":frag,
   "unadjusted_old_vs_novel_slopes":rawslopes}
  print("NRSA2023_XEMBED_PAIR_SOURCE_MODERATION",key,json.dumps(result[key],sort_keys=True),flush=True)
 output={"schema":"gesui_EPA2023_XEMBED_EPT_source_status_matched_pair_slope_v1",
   "EPA_original_five_source_SHA":SHA,
   "original_2008_2013_2018_uid":2648,
   "original_2023_berw_uid":1165,
   "original_target_HUC8":693,
   "old_source_ever_HUC8":479,
   "old_source_never_HUC8":214,
   "oldtrain_only_1185_UID_551_HUC8_zero_target_overlap":True,
   "one_primary_difference_in_XEMBED_EPT_slope_test":"grain_214_pairs",
   "result_labels_posthoc_not_causal":True,
   "no_sewage_derived_NP_or_aquatic_emergence":True,
   "results":result}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("NRSA2023_PAIR_STATUS_XEMBED_SUMMARY",json.dumps({
   key:{"status":d["status"],"n":d.get("paired_difference_XEMBED_mean_slope",{}).get("n_matched_pairs"),
    "slope_pp":d.get("paired_difference_XEMBED_mean_slope",{}).get("difference_in_EPT_percentage_points_per_1SD_XEMBED"),
    "HC3_CI":d.get("paired_difference_XEMBED_mean_slope",{}).get("HC3_95_difference_in_EPT_percentage_points_per_1SD"),
    "p":d.get("paired_difference_XEMBED_mean_slope",{}).get("HC3_exploratory_p"),
    "pairCV_MSEgain":(d.get("whole_pair_heldout_prediction") or {}).get("added_moderation_pairheldout_MSE_gain_pct")}
    for key,d in result.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
