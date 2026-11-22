#!/usr/bin/env python3
"""Posthoc EPA NRSA2324 nutrient vs habitat additional prediction, same 4 cohorts.

After the TP x canopy interaction analysis. Not sewage effect or causal mediation.
The source/outcome sample and fold match the earlier interaction source exactly.
"""
import argparse,hashlib,json,pathlib
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge
import gesui_nrsa2324_tp_canopy_benthos_interaction_20261122 as original

BASE=["lat","lon","order","year","season_sin","season_cos","log_cond"]
N=["log_tn","log_tp"]
H=["log_width","canopy"]
CAT=["state","eco9"]
MODELS={
 "B_geography":BASE,
 "B_plus_N_ambient":BASE+N,
 "B_plus_H_field":BASE+H,
 "B_plus_N_plus_H":BASE+N+H}
EXPECTED={
 ("BERW","MMI_BENT"):257.58943288675243,
 ("BERW","EPT_fraction"):0.040524098828369694,
 ("BETB","MMI_BENT"):290.3427640085802,
 ("BETB","EPT_fraction"):0.03666350258154949}
def oof(d,y,num):
 fold=np.asarray([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5 for x in d.huc8])
 pred=np.full(len(d),np.nan)
 for k in range(5):
  tr=fold!=k;te=~tr
  prep=ColumnTransformer([
    ("numerical",StandardScaler(),num),
    ("categories",OneHotEncoder(handle_unknown="ignore",sparse_output=False),CAT)],
    sparse_threshold=0)
  xt=prep.fit_transform(d.loc[tr,num+CAT])
  xs=prep.transform(d.loc[te,num+CAT])
  pred[te]=Ridge(alpha=10,solver="svd").fit(xt,y[tr]).predict(xs)
 if not np.isfinite(pred).all():raise ValueError("Nonfinite heldout prediction")
 errors=(y-pred)**2
 return (pd.Series(errors,index=d.huc8.to_numpy()).groupby(level=0).mean().sort_index(),
         float(np.mean(errors)))
def gain(a,b):
 A=a.to_numpy();B=b.reindex(a.index).to_numpy()
 ix=np.random.default_rng(20261122).integers(0,len(A),size=(2000,len(A)))
 arr=100*(A[ix].mean(axis=1)-B[ix].mean(axis=1))/A[ix].mean(axis=1)
 return {"HUC8_equal_MSE_gain_pct":float(100*(A.mean()-B.mean())/A.mean()),
   "HUC8_fixedOOF_bootstrap95_pct":[float(v) for v in np.quantile(arr,[.025,.975])]}
def main():
 ap=argparse.ArgumentParser()
 for k in original.SHA:ap.add_argument("--"+k,required=True,type=pathlib.Path)
 ap.add_argument("--out",required=True,type=pathlib.Path)
 a=ap.parse_args()
 t={}
 for k,pin in original.SHA.items():
  p=getattr(a,k)
  if original.sha(p)!=pin:raise RuntimeError("Original EPA SHA changed "+k)
  t[k]=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
 panel,qa,n=original.panel(t)
 if n!=2139 or sum(qa.values())!=0:raise RuntimeError("Five-source UID join changed")
 res={}
 for m in ("BERW","BETB"):
  res[m]={}
  for response in ("MMI_BENT","EPT_fraction"):
   d=original.cohort(panel,m,response)
   if len(d)<250 or d.huc8.nunique()<100:raise RuntimeError("Previously established cohort gate now fails")
   y=d[response].to_numpy(float)
   fits={name:oof(d,y,cols) for name,cols in MODELS.items()}
   ref=float(fits["B_plus_N_plus_H"][0].mean())
   exp=EXPECTED[(m,response)]
   if not np.isclose(ref,exp,rtol=1e-7,atol=1e-9):
    raise RuntimeError("Original HUC8 heldout both-main-effects reference MSE not reproduced "+str((m,response,ref,exp)))
   results={}
   for name,(group,visit) in fits.items():
    results[name]={"HUC8_equal_MSE":float(group.mean()),"visit_weighted_MSE":visit}
   comparison={}
   for label,a0,b0 in [
    ("ambient_N_after_geography","B_geography","B_plus_N_ambient"),
    ("habitat_H_after_geography","B_geography","B_plus_H_field"),
    ("ambient_N_after_habitat","B_plus_H_field","B_plus_N_plus_H"),
    ("habitat_H_after_ambient_N","B_plus_N_ambient","B_plus_N_plus_H")]:
    comparison[label]=gain(fits[a0][0],fits[b0][0])
    va=fits[a0][1];vb=fits[b0][1]
    comparison[label]["visit_weighted_MSE_gain_pct"]=float(100*(va-vb)/va)
   res[m][response]={"n_UID":len(d),"physical_sites":int(d.site.nunique()),
       "n_HUC8":int(d.huc8.nunique()),"model_MSE":results,"comparisons":comparison,
       "original_HUC8_ref_MSE_reproduced":True}
   print("NRSA_ADDITIVE_DRIVER_FOURWAY",m,response,json.dumps(res[m][response],sort_keys=True),flush=True)
 report={"schema":"gesui_NRSA2324_samecohort_additive_ambientNP_vs_physical_canopy_width_posthoc_v1",
  "source_5table_UID":n,"original_SHA":original.SHA,"reference_exact_reproduced":True,
  "models":MODELS,"outcomes":res,"posthoc_after_interaction_seen":True,
  "HUC8_group_disjoint":True,"not_WWTP_specific":True,
  "source_mediator_or_causality_unproven":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
 print("NRSA_ADDITIVE_DRIVER_SUMMARY",json.dumps({m:{r:
   {k:round(v["HUC8_equal_MSE_gain_pct"],3) for k,v in z["comparisons"].items()}
   for r,z in o.items()} for m,o in res.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
