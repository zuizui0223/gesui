#!/usr/bin/env python3
"""Posthoc EPA NRSA2023/24 flow/area stress test of physical width and nutrients.

Models original water TN/TP and benthic MMI/EPT, not sewage-derived P.
Follows frozen docs/PROTOCOL_20261122_NRSA2324_FLOW_AREA_WIDTH_IDENTIFICATION.md.
"""
import argparse,collections,hashlib,json,pathlib
import numpy as np,pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge
import gesui_nrsa2324_tp_canopy_benthos_interaction_20261122 as orig
import gesui_nrsa2324_ambient_np_habitat_additive_posthoc_20261122 as prev

NUM_GEO=prev.BASE
N=["log_tn","log_tp"]
C=["canopy"]
W=["log_width"]
FLOW=["log_flow"]
AREA=["log_area"]
CAT=["state","eco9","urban"]
CORE_MODELS={
 "B_NC_F":NUM_GEO+N+C+FLOW,
 "B_NC_F_W":NUM_GEO+N+C+FLOW+W,
 "B_NC_A":NUM_GEO+N+C+AREA,
 "B_NC_A_W":NUM_GEO+N+C+AREA+W,
 "B_NC_FA":NUM_GEO+N+C+FLOW+AREA,
 "B_NC_FA_W":NUM_GEO+N+C+FLOW+AREA+W,
 "B_C_FA_W":NUM_GEO+C+FLOW+AREA+W
}
COMPARE=[
 ("width_after_modeled_flow","B_NC_F","B_NC_F_W"),
 ("width_after_upstream_area","B_NC_A","B_NC_A_W"),
 ("width_after_flow_and_area","B_NC_FA","B_NC_FA_W"),
 ("nutrients_after_hydro_area_width","B_C_FA_W","B_NC_FA_W")]
EXTENDED={
 "full_size_sediment_without_width":NUM_GEO+N+C+FLOW+AREA+["log_flow_squared","embedded"],
 "full_size_sediment_with_width":NUM_GEO+N+C+FLOW+AREA+["log_flow_squared","embedded"]+W
}
ORIGINAL_EXPECTED=prev.EXPECTED

def normalize_hydro(panel,table):
 s=table.drop_duplicates("UID").set_index("UID")
 z=panel.copy()
 rawflow=pd.to_numeric(s.MEAN_FLOW,errors="coerce").reindex(z.index)
 rawarea=pd.to_numeric(s.WSAREASQKM,errors="coerce").reindex(z.index)
 z["log_flow"]=np.log10(rawflow.where(rawflow>0))
 z["log_area"]=np.log10(rawarea.where(rawarea>0))
 z["log_flow_squared"]=z.log_flow**2
 z["urban"]=s.URBN_NRS23.astype(str).str.strip().reindex(z.index)
 z["flow_size"]=s.FLOW_SIZE.astype(str).str.strip().reindex(z.index)
 return z
def numeric_embedded(z,physical):
 p=physical.drop_duplicates("UID").set_index("UID")
 z["embedded"]=pd.to_numeric(p.XEMBED,errors="coerce").reindex(z.index)
 z.loc[~z.embedded.between(0,100),"embedded"]=np.nan
 return z
def covariate_support(d,method):
 raw=d.copy()
 nonnull=raw[["log_flow","log_area","log_width","canopy"]].replace([np.inf,-np.inf],np.nan)
 core=raw.dropna(subset=["log_flow","log_area","urban"]).copy()
 core=core[core.urban.isin(["Urban","Non-Urban"])]
 core=core.sort_index()
 good=core[["log_flow","log_width","log_area"]].dropna()
 pearson={}
 for col in ("log_flow","log_area"):
  pearson["corr_log_width_"+col]=(float(good.log_width.corr(good[col])) if len(good)>3 else None)
 support={"original_complete":len(raw),"all_hydro_complete":len(core),
   "all_hydro_unique_physical_sites":int(core.site.nunique()),
   "all_hydro_HUC8":int(core.huc8.nunique()),
   "correlations_original_complete_hydro":pearson,
   "median_flow_cfs":float(10**core.log_flow.median()) if len(core) else None,
   "median_area_km2":float(10**core.log_area.median()) if len(core) else None,
   "median_width_m":float(10**core.log_width.median()) if len(core) else None}
 return core,support
def oof(d,y,numcols,catcols):
 fold=np.array([int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5 for h in d.huc8],int)
 ns={str(f):int(d.loc[fold==f,"huc8"].nunique()) for f in range(5)}
 if min(ns.values())<3:raise ValueError("Too few HUC8 in a heldout fold")
 pred=np.full(len(y),np.nan)
 for f in range(5):
  train=fold!=f;test=~train
  enc=ColumnTransformer([
    ("continuous",StandardScaler(),numcols),
    ("classes",OneHotEncoder(handle_unknown="ignore",sparse_output=False),catcols)
   ],sparse_threshold=0)
  a=enc.fit_transform(d.loc[train,numcols+catcols])
  b=enc.transform(d.loc[test,numcols+catcols])
  reg=Ridge(alpha=10,solver="svd",fit_intercept=True).fit(a,y[train])
  pred[test]=reg.predict(b)
 if not np.isfinite(pred).all():raise ValueError("Unpredicted original heldout rows")
 err=(y-pred)**2
 groups=pd.Series(err,index=d.huc8.to_numpy()).groupby(level=0).mean().sort_index()
 return groups,float(err.mean()),ns
def test_nested(d,y,numsets,cats,comparisons):
 scored={k:oof(d,y,nums,cats) for k,nums in numsets.items()}
 comp={}
 for name,m0,m1 in comparisons:
  g0,v0,_=scored[m0];g1,v1,_=scored[m1]
  comp[name]={**prev.gain(g0,g1),
   "visit_weighted_MSE_gain_pct":float(100*(v0-v1)/v0)}
 mses={k:{"HUC8_equal_MSE":float(v[0].mean()),"visit_MSE":v[1]} for k,v in scored.items()}
 return {"comparisons":comp,"HUC8_fold_counts":list(scored.values())[0][2],
   "model_MSEs":mses}
def main():
 p=argparse.ArgumentParser()
 for k in orig.SHA:p.add_argument("--"+k,required=True,type=pathlib.Path)
 p.add_argument("--out",required=True,type=pathlib.Path)
 x=p.parse_args()
 src={}
 for k,expected in orig.SHA.items():
  path=getattr(x,k)
  if orig.sha(path)!=expected:raise ValueError("Pinned original EPA national data changed "+k)
  src[k]=pd.read_csv(path,dtype=str,keep_default_na=False,low_memory=False)
 panel,qa,n=orig.panel(src)
 if n!=2139 or sum(qa.values()):raise ValueError("Original 5-source UID/date/site join changed")
 panel=normalize_hydro(panel,src["siteinfo"])
 panel=numeric_embedded(panel,src["physicalhabitat"])
 out={}
 for method in ("BERW","BETB"):
  out[method]={}
  for endpoint in ("MMI_BENT","EPT_fraction"):
   original=orig.cohort(panel,method,endpoint)
   original_y=original[endpoint].to_numpy(float)
   previous_ref=prev.oof(original,original_y,NUM_GEO+N+C+W)[0].mean()
   exp=ORIGINAL_EXPECTED[(method,endpoint)]
   if not np.isclose(previous_ref,exp,rtol=1e-8,atol=1e-10):
    raise ValueError("Pre-stress original source model MSE not reproduced for "+method+"/"+endpoint)
   d,support=covariate_support(original,method)
   status={"baseline_exact_reproduced":True,"source_support":support,
     "method":"field_midchannel_canopy" if method=="BERW" else "field_bank_canopy"}
   if len(d)<250 or d.huc8.nunique()<100:
    status.update({"status":"SOURCE_STOP","reason":"modeled flow+area complete case below source gate"})
    out[method][endpoint]=status;continue
   y=d[endpoint].to_numpy(float)
   status["core"]=test_nested(d,y,CORE_MODELS,CAT,COMPARE)
   status["status"]="CORE_NUMERIC_COMPLETE"
   ext=d.dropna(subset=["embedded","log_flow_squared"]).copy()
   ext=ext[~ext.flow_size.isin(["","-8","NA","None"])].sort_index()
   status["extended_source_support"]={
     "UID":len(ext),"physical_sites":int(ext.site.nunique()),"HUC8":int(ext.huc8.nunique())}
   if len(ext)>=250 and ext.huc8.nunique()>=100:
    status["extended"]=test_nested(ext,ext[endpoint].to_numpy(float),EXTENDED,
     CAT+["flow_size"],[("width_after_saturated_hydrology_sediment",
       "full_size_sediment_without_width","full_size_sediment_with_width")])
    status["extended_status"]="NUMERIC_COMPLETE"
   else:status["extended_status"]="SOURCE_STOP_MISSING_SEDIMENT_OR_FLOW_CATEGORY"
   out[method][endpoint]=status
   print("NRSA_HYDROLOGY_SIZE_STRESS",method,endpoint,json.dumps(status,sort_keys=True),flush=True)
 report={"schema":"gesui_NRSA2324_modeled_flow_area_width_benthos_posthoc_v1",
  "original_sha":orig.SHA,"original_joint_UID":n,
  "native_source_QA":qa,"model_specifications":CORE_MODELS,
  "original_method_cohorts_preserved":True,
  "posthoc":True,"no_WWTP_discharge_assigned":True,
  "modeled_mean_flow_not_field_measured_flow":True,
  "not_riparian_restoration_causal":True,
  "results":out}
 x.out.parent.mkdir(parents=True,exist_ok=True)
 x.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("NRSA_FLOW_AREA_WIDTH_STRESS_SUMMARY",json.dumps(
   {m:{r:{"status":z["status"],"n":z["source_support"]["all_hydro_complete"],
    "HUC8":z["source_support"]["all_hydro_HUC8"],
    "flow_plus_width_gain":z.get("core",{}).get("comparisons",{}).get("width_after_modeled_flow",{}).get("HUC8_equal_MSE_gain_pct"),
    "flowarea_plus_width_gain":z.get("core",{}).get("comparisons",{}).get("width_after_flow_and_area",{}).get("HUC8_equal_MSE_gain_pct"),
    "nutrient_gain":z.get("core",{}).get("comparisons",{}).get("nutrients_after_hydro_area_width",{}).get("HUC8_equal_MSE_gain_pct")}
    for r,z in rows.items()} for m,rows in out.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
