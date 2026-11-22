#!/usr/bin/env python3
"""Adversarial same-flow and same-width support test of NRSA wadeable/boatable sign.
Original source-only cutoffs frozen in independent artifact before biotic outcomes.
Observational, posthoc; no WWTP discharge or community-causal inference.
"""
import argparse,json,pathlib
import numpy as np,pandas as pd,statsmodels.api as sm
import gesui_nrsa2324_tp_canopy_benthos_interaction_20261122 as orig
import gesui_nrsa2324_flow_area_width_stress_20261122 as stress

SPECS={
 "RV_source_only":(100,70),
 "RV_common_flow":(80,50),
 "RV_common_flow_and_width":(60,40)}
REF_MODEL={"without_width":stress.CORE_MODELS["B_NC_FA"],
           "with_width":stress.CORE_MODELS["B_NC_FA_W"]}
CAT=stress.CAT
def descriptive(d,y):
 nums=REF_MODEL["with_width"]
 df=d[nums].astype(float)
 z=(df-df.mean())/df.std(ddof=0).replace(0,1)
 cats=pd.get_dummies(d[CAT],drop_first=True,dtype=float)
 X=pd.concat([z.reset_index(drop=True),cats.reset_index(drop=True)],axis=1)
 X=sm.add_constant(X,has_constant="add")
 yy=(y-np.mean(y))/np.std(y,ddof=0)
 try:
  fit=sm.OLS(yy,X.to_numpy(float)).fit(cov_type="cluster",
      cov_kwds={"groups":d.huc8.to_numpy(),"use_correction":True})
  i=list(X.columns).index("log_width");ci=fit.conf_int()[i]
  return {"status":"DESCRIPTIVE_COMPLETE",
   "beta_width_ySD_per_1SD_logwidth":float(fit.params[i]),
   "cluster_HUC8_CI95":[float(ci[0]),float(ci[1])],
   "p_uncorrected":float(fit.pvalues[i]),
   "OLS_rank":int(np.linalg.matrix_rank(X.to_numpy(float))),
   "OLS_covariate_columns":len(X.columns)}
 except Exception as e:
  return {"status":"DESCRIPTIVE_FAILED","error":str(e)[:300]}
def source_and_results(d,method,q):
 s=q["methods"][method]["survey_subsets"]
 if not all(k in s for k in SPECS):raise ValueError("Source-only frozen subsets missing")
 fmin,fmax=q["fixed_common_log10_flow_bounds"]
 wmin,wmax=q["fixed_common_log10_width_bounds"]
 results={}
 for k,(min_n,min_huc) in SPECS.items():
  sel=d[d.flow_size.eq("RV")].copy()
  if k!="RV_source_only":
   sel=sel[sel.log_flow.between(fmin,fmax,inclusive="both")].copy()
  if k=="RV_common_flow_and_width":
   sel=sel[sel.log_width.between(wmin,wmax,inclusive="both")].copy()
  subset={"source_method":method,"source_only":s[k],
      "observed_UID":len(sel),"physical_sites":int(sel.site.nunique()),
      "independent_HUC8":int(sel.huc8.nunique()),"states":int(sel.state.nunique())}
  if not s[k]["source_GO"]:
   results[k]={**subset,"status":"SOURCE_STOP_BY_FROZEN_ORIGINAL_SOURCE_GATE"};continue
  if len(sel)<min_n or sel.huc8.nunique()<min_huc:
   results[k]={**subset,"status":"SOURCE_STOP_AFTER_COMPLETE_CASE"};continue
  sel=sel.sort_index();y=sel["target_endpoint"].to_numpy(float)
  try:
   base,visit0,folds=stress.oof(sel,y,REF_MODEL["without_width"],CAT)
   full,visit1,_=stress.oof(sel,y,REF_MODEL["with_width"],CAT)
   metric=stress.prev.gain(base,full)
   metric["visit_weighted_MSE_gain_pct"]=float(100*(visit0-visit1)/visit0)
   subset["HUC8_fivefold_counts"]=folds
   subset["HUC8_prediction"]=metric
   subset["adjusted_width"]=descriptive(sel,y)
   subset["status"]="NUMERIC_COMPLETE"
  except Exception as e:
   subset["status"]="ANALYSIS_STOP"
   subset["error"]=type(e).__name__+": "+str(e)[:300]
  results[k]=subset
 return results
def main():
 ap=argparse.ArgumentParser()
 for k in orig.SHA:ap.add_argument("--"+k,type=pathlib.Path,required=True)
 ap.add_argument("--source-qa",type=pathlib.Path,required=True)
 ap.add_argument("--out",type=pathlib.Path,required=True)
 a=ap.parse_args();qa=json.loads(a.source_qa.read_text())
 if qa.get("schema")!="gesui_NRSA2324_original_RV_flow_width_overlap_source_only_v1":
  raise ValueError("Source-only preceding artifact schema changed")
 if qa["source_SHA"]!={k:orig.SHA[k] for k in ("siteinfo","physicalhabitat")}:
  raise ValueError("Numerical job did not use same EPA source-locked original tables")
 raw={}
 for k,pin in orig.SHA.items():
  p=getattr(a,k)
  if orig.sha(p)!=pin:raise ValueError("Original national EPA 5-table SHA mismatch")
  raw[k]=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
 origpanel,checks,n=orig.panel(raw)
 if n!=qa["common_UID_all_original"] or sum(checks.values()):
  raise RuntimeError("Source-only and numeric original UID source no longer equal")
 panel=stress.numeric_embedded(stress.normalize_hydro(origpanel,raw["siteinfo"]),raw["physicalhabitat"])
 results={}
 for method in ("BERW","BETB"):
  results[method]={}
  for endpoint in ("MMI_BENT","EPT_fraction"):
   d=orig.cohort(panel,method,endpoint)
   d,support=stress.covariate_support(d,method)
   d["target_endpoint"]=d[endpoint].to_numpy(float)
   analysis=source_and_results(d,method,qa)
   results[method][endpoint]={"original_complete_case_n":len(d),
      "cohort_HUC8":d.huc8.nunique(),"source_support":support,
      "fixed_source_overlap_subsets":analysis}
   print("NRSA_RV_OVERLAP_RESULT",method,endpoint,json.dumps(analysis,sort_keys=True),flush=True)
 report={"schema":"gesui_NRSA2324_RV_flow_width_overlap_hydrologic_sign_falsifier_v1",
  "original_EPA_SHA":orig.SHA,
  "source_only_preoutcome_RV_cutoffs":{k:v for k,v in qa.items() if k.startswith("fixed_common_")},
  "source_only_artifact_SHA256":__import__("hashlib").sha256(a.source_qa.read_bytes()).hexdigest(),
  "same_original_UID":n,"all_results_explicitly_posthoc":True,
  "no_sewage_exposure_reconstruction":True,
  "comparison_is_not_randomized":True,
  "results":results}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("NRSA_RV_OVERLAP_SUMMARY",json.dumps({m:{endpoint:
  {k:{"status":val["status"],"UID":val["observed_UID"],
     "HUC8":val["independent_HUC8"],
     "beta":val.get("adjusted_width",{}).get("beta_width_ySD_per_1SD_logwidth"),
     "MSE_gain":val.get("HUC8_prediction",{}).get("HUC8_equal_MSE_gain_pct")}
   for k,val in z["fixed_source_overlap_subsets"].items()}
  for endpoint,z in g.items()} for m,g in results.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
