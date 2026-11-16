#!/usr/bin/env python3
"""Posthoc but source-frozen physical-habitat falsification of NRSA2023-24
same-visit ambient TN/TP, areal AFDM, Chla and Chla/AFDM responses.

Uses authentic EPA PHab original source SHA + earlier pre-outcome n2129 UID
map, original n1828 EPA biofilm strict source cohort. Wide and full tiers were
specified on source availability BEFORE reading PHab numeric effects.

No WWTP effluent, production, true direct sunlight, ecosystem services,
or benthic diversity response is estimated here.
"""
import hashlib, json, os, pathlib, urllib.request, time
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
import statsmodels.api as sm

DIR=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"
DIR.mkdir(parents=True,exist_ok=True)
SOURCE_URL="https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_physicalhabitat_metrics.csv"
SOURCE_SHA="8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4"
source_path=pathlib.Path(os.environ["FROZEN_PHAB_SOURCE_QA"])
source_q=json.loads(source_path.read_text())
if source_q["phab"]["source_sha256"]!=SOURCE_SHA or source_q["native_site_visits_matched"]!=2129:
 raise RuntimeError("Original pre-numeric PHab source SHA/UID and quality gate changed")
if not source_q["gate_phab_support_min400_100HUC8_20states"]:
 raise RuntimeError("Original PHab eligibility source gate did not pass")
frozen_uids=set(source_q["source_eligible_original_UIDs"])
if len(frozen_uids)!=2129:raise RuntimeError("Frozen source UID manifest not original cohort")
previous=json.loads((DIR/"OBSERVED_NRSA_NUTRIENT_PERIPHYTON_EPT_TWO_ERAS_NO_WWTP_CAUSAL.json").read_text())
original=previous["years"]["2023_24"]["strict"]
panel=pd.read_csv(DIR/"AUDIT_ONLY_2023_24_PANEL.csv",
                  dtype={"UID":str,"state":str,"eco9":str,"huc8":str,"site_id":str})
if len(panel)!=2129 or not set(panel.UID)==frozen_uids or panel.UID.duplicated().any():
 raise RuntimeError("Pre-PHab 2,129 source-only visits were changed")
p=DIR/"original_epa_phab.csv"
for i in range(3):
 try:
  req=urllib.request.Request(SOURCE_URL,headers={"User-Agent":"gesui-physical-habitat-sameUID-numeric/20261114"})
  with urllib.request.urlopen(req,timeout=100) as u,p.open("wb") as f:
   while x:=u.read(1048576):f.write(x)
  if hashlib.sha256(p.read_bytes()).hexdigest()!=SOURCE_SHA:raise RuntimeError("New EPA PHab source SHA mismatch")
  break
 except Exception:
  if i==2:raise
  time.sleep(2**i)
phab=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
p.unlink()
if phab.UID.duplicated().any():raise RuntimeError("PHab one-to-one UID corrupted")
if not frozen_uids.issubset(set(phab.UID)):raise RuntimeError("PHab frozen UID support is not full")
phab=phab.set_index("UID")
wide=["XPCAN","XGVEG","XFC_OHV","XWIDTH"]
extended=["XPCAN","XGVEG","XFC_OHV","XWIDTH","XCDENMID","XEMBED"]
for k in extended:
 if k not in phab.columns:raise RuntimeError("Precommitted PHab measure absent: "+k)
 phab[k]=pd.to_numeric(phab[k],errors="coerce")
phab.loc[~(phab.XWIDTH>0),"XWIDTH"]=np.nan
phab["log_width"]=np.log10(phab.XWIDTH)
# PHab source values are used untruncated; vegetation layers can overlap and
# their sum is not necessarily a mutually exclusive 0–1 fractional cover.
# XPCAN and XFC_OHV are reported areal proportions; no after-sign clipping.
qc_values={}
for k in extended:
 v=phab[k].dropna()
 qc_values[k]={"present_n":int(len(v)),"negative_n":int((v<0).sum()),
  "above_one_n":int((v>1).sum()) if k in ["XPCAN","XGVEG","XFC_OHV"] else None,
  "max":float(v.max()) if len(v) else None}
panel=panel.set_index("UID").join(phab[extended+["log_width"]],how="left")
panel.index.name="UID"
core=["log_tn","log_tp","log_cond","lat","lon","year","stream_order","season_sin",
      "season_cos","log_chla","log_afdm"]
flag=panel.chla_flag.fillna("").eq("") & panel.afdm_flag.fillna("").eq("") & panel.tn_mgL_rawflag.fillna("").eq("") & panel.tp_mgL_rawflag.fillna("").eq("")
strict=panel[flag].dropna(subset=core)
strict=strict[(strict.huc8!="")&(strict.state!="")&(strict.eco9!="")]
if len(strict)!=1828 or strict.huc8.nunique()!=946 or strict.state.nunique()!=48:
 raise RuntimeError("Original pre-PHab strict complete-case sample n1828/946/48 changed")
strict=strict.copy()
strict["log_pigment_per_AFDM"]=strict.log_chla-strict.log_afdm
for modelname,cols,minimum in [("wide",["XPCAN","XGVEG","XFC_OHV","log_width"],850),
                                ("dense_substrate",["XPCAN","XGVEG","XFC_OHV","log_width","XCDENMID","XEMBED"],450)]:
 cohort=strict.dropna(subset=cols).sort_index().copy()
 if len(cohort)<minimum or cohort.huc8.nunique()<100 or cohort.state.nunique()<20:
  raise RuntimeError(f"Source original physical metric support failure {modelname}: n{len(cohort)}")
 if modelname=="wide" and any(cohort[X].isna().any() for X in cols):
  raise RuntimeError("PHab wide predictor missing in supposed complete cases")
 if modelname=="wide":widecohort=cohort
 else:densecohort=cohort
def fold(src):
 return np.array([int(hashlib.sha256(str(v).encode()).hexdigest()[:8],16)%5 for v in src.huc8])
def features(d,cols):
 C=pd.concat([d[["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]].astype(float),
              pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 # Fixed baseline built once on same within-tier outcome sample.
 physical=d[cols].astype(float)
 NP=d[["log_tn","log_tp"]].astype(float)
 return {"baseline":C,
         "nutrients":pd.concat([C,NP],axis=1),
         "physical":pd.concat([C,physical],axis=1),
         "physical+nutrients":pd.concat([C,physical,NP],axis=1)}
def ridge_oos_mse(X,Y,folds):
 out=np.full(len(Y),np.nan)
 a=X.to_numpy(dtype=float)
 for k in range(5):
  tr=folds!=k;te=~tr
  if tr.sum()<100 or te.sum()<60:raise RuntimeError("Fold original HUC8 sample support insufficient")
  scalar=StandardScaler().fit(a[tr])
  z=scalar.transform(a)
  model=Ridge(alpha=10,fit_intercept=True,solver="svd")
  model.fit(z[tr],Y[tr])
  out[te]=model.predict(z[te])
 if not np.isfinite(out).all():raise RuntimeError("Heldout prediction missing")
 return float(np.mean((Y-out)**2))
def slope(d,X,Y):
 # Original response in native fraction or log10 units; TN/TP scaled by
 # full-sample 1 SD, with HUC8 clustered descriptive uncertainty.
 A=X.copy()
 for nm in ["log_tn","log_tp"]:
  A[nm]=(A[nm]-A[nm].mean())/A[nm].std(ddof=0)
 ymod=sm.OLS(Y,sm.add_constant(A.to_numpy(float),has_constant="add")).fit(
   cov_type="cluster",cov_kwds={"groups":d.huc8,"use_correction":True})
 return {nm:{"beta_per_1SD_log10_nutrient":float(ymod.params[A.columns.get_loc(nm)+1]),
    "HUC8_cluster_descriptive95":list(map(float,ymod.conf_int()[A.columns.get_loc(nm)+1]))}
   for nm in ("log_tn","log_tp")}
results={}
for name,d in [("wide",widecohort),("dense_substrate",densecohort)]:
 cols=["XPCAN","XGVEG","XFC_OHV","log_width"] + (["XCDENMID","XEMBED"] if name=="dense_substrate" else [])
 X=features(d,cols);folds=fold(d)
 ans={"strict_pre_phab_original_sample_n":1828,"complete_case_n":len(d),
      "n_source_HUC8":d.huc8.nunique(),"n_states":d.state.nunique(),"n_distinct_physical_SITE_ID":d.site_id.nunique(),
      "phys_covariates":cols,"responses":{}}
 for resp in ["log_pigment_per_AFDM","log_chla","log_afdm"]:
  y=d[resp].to_numpy(dtype=float)
  mse={m:ridge_oos_mse(mx,y,folds) for m,mx in X.items()}
  gain=lambda a,b:100*(mse[a]-mse[b])/mse[a]
  ans["responses"][resp]={
    "MSEs_HUC8_fivefold":mse,
    "nutrient_added_gain_to_geo_pct":gain("baseline","nutrients"),
    "nutrient_added_gain_to_phab_pct":gain("physical","physical+nutrients"),
    "PHab_added_gain_to_nutrients_pct":gain("nutrients","physical+nutrients"),
    "TP_TN_slope_geo_control":slope(d,X["nutrients"],y),
    "TP_TN_slope_PHab_control":slope(d,X["physical+nutrients"],y)}
 results[name]=ans
out={
 "schema":"gesui_NRSA2324_original_ambient_NP_PHAB_PIGMENT_AFDM_falsification_20261114",
 "original_2129_source_frozen_UIDs":len(frozen_uids),"original_strict_baseline_n":len(strict),
 "publisher_PHAB_original_sha256":SOURCE_SHA,
 "phab_data_semantic":"canopy/overhanging vegetation/ground vegetation/width; dense canopy and bed embedding only where native EPA measurements support",
 "phab_raw_metric_quality_non_numeric":qc_values,
 "new_selected_on_previous_ratio_result_posthoc":True,
 "no_EPA_WWTP_effluent_NP_attribution":True,
 "light_shading_physical_proxy_not_measured_RSDS":True,
 "AFDM_is_standing_organic_mass_not_net_primary_productivity":True,
 "ratio_is_algebraic_contrast_not_an_independent_field_endpoint":True,
 "no_benthic_2023_24_diversity_available":True,
 "models":results,
 "source_qa_artifact":"NRSA2324_PHAB_SOURCE_ONLY_QA_20261114.json",
 "source_bias_and_nulls_must_be_preserved":True}
f=DIR/"NRSA2324_PHAB_NUTRIENT_PIGMENT_AREAL_BIOMASS_POSTHOC.json"
f.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA2324_PHAB_NUMERIC_START")
print(json.dumps({k:{"n":v["complete_case_n"],"HUC8":v["n_source_HUC8"],
  "responses":{m:{"deltaN_geo":p["nutrient_added_gain_to_geo_pct"],
              "deltaN_after_PHab":p["nutrient_added_gain_to_phab_pct"],
              "delta_PHab_after_N":p["PHab_added_gain_to_nutrients_pct"],
              "slopes_geo":p["TP_TN_slope_geo_control"],
              "slopes_PHab":p["TP_TN_slope_PHab_control"]}
        for m,p in v["responses"].items()}}
  for k,v in results.items()},sort_keys=True))
print("NRSA2324_PHAB_NUMERIC_END")
