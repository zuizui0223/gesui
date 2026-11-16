#!/usr/bin/env python3
"""Exploratory NRSA2023/24 genuine same-UID nutrient × biofilm mass-pigment allocation.

Triggered only after original EPA 2023/24 benthic Count/MMI source-only failure.
Reuses EXACT original PR73 2129 site cohort and strict n1828 source+quality
without changing time window, selected sites or nutrient endpoint criteria.

Areal Chla in ug/cm2 and AFDM in mg/cm2. log10(chla/afdm) represents
pigment per AFDM dry organic mass [ug/mg], NOT photosynthetic rate, autotroph
fraction, wastewater effects or ecosystem service. The ratio is algebraically
related to Chla and AFDM; never claim independent validation of previous TP sign.
"""
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from scipy.stats import spearmanr

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"
SOURCE=ROOT/"AUDIT_ONLY_2023_24_PANEL.csv"
PREV=ROOT/"OBSERVED_NRSA_NUTRIENT_PERIPHYTON_EPT_TWO_ERAS_NO_WWTP_CAUSAL.json"
OUT=ROOT/"NRSA2324_PIGMENT_PER_ORGANIC_MASS_EXPLORATORY.json"

df=pd.read_csv(SOURCE,dtype={"UID":str,"site_id":str,"state":str,"eco9":str,"huc8":str})
prior=json.loads(PREV.read_text())
ref=prior["years"]["2023_24"]["strict"]
allowed= df.chla_flag.fillna("").eq("") & df.afdm_flag.fillna("").eq("") & df.tn_mgL_rawflag.fillna("").eq("") & df.tp_mgL_rawflag.fillna("").eq("")
basevars=["lat","lon","stream_order","year","season_sin","season_cos","log_cond","log_tn","log_tp","log_chla","log_afdm"]
d=df[allowed].dropna(subset=basevars).copy()
d=d[(d.huc8.notna())&(d.huc8!="")&(d.state.notna())&(d.state!="")&(d.eco9.notna())&(d.eco9!="")]
d=d.sort_values("UID").reset_index(drop=True)
if len(d)!=1828 or d.huc8.nunique()!=946 or d.state.nunique()!=48 or len(df)!=2129:
 raise RuntimeError("Original strict PR73 data source cohort changed: no estimation")

d["pigment_ug_per_mg"]=10.0**(d.log_chla-d.log_afdm)
d["log_pigment_per_afdm"]=d.log_chla-d.log_afdm
if not np.isfinite(d.log_pigment_per_afdm).all():raise RuntimeError("Nonpositive original Chla/AFDM was not correctly excluded")
for nm in ("log_chla","log_afdm"):
 prev=float(ref["responses"][nm]["huc8_disjoint_5fold"]["baseline"]["MSE"])
 if not np.isfinite(prev):raise RuntimeError("Old source ledger malformed")

# Whole HUC8 heldout; no new fold selection after observing ratio.
fold=np.asarray([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5 for x in d.huc8])
basecols=["lat","lon","stream_order","year","season_sin","season_cos","log_cond"]
X0=pd.concat([d[basecols].astype(float),pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
preds={
 "baseline":X0,
 "plus_NP":X0.assign(log_tn=d.log_tn.to_numpy(float),log_tp=d.log_tp.to_numpy(float)),
 "plus_NP_AFDM":X0.assign(log_tn=d.log_tn.to_numpy(float),log_tp=d.log_tp.to_numpy(float),log_afdm=d.log_afdm.to_numpy(float))
}
def heldout(y, X):
 pred=np.full(len(d),np.nan)
 a=X.to_numpy(dtype=float)
 for k in range(5):
  tr=fold!=k;te=fold==k
  if tr.sum()<100 or te.sum()<100:raise RuntimeError("Invalid HUC8 folds")
  scalar=StandardScaler().fit(a[tr])
  z=scalar.transform(a)
  reg=Ridge(alpha=10,solver="svd",fit_intercept=True).fit(z[tr],y[tr])
  pred[te]=reg.predict(z[te])
 if not np.isfinite(pred).all():raise RuntimeError("Original HUC8 fold prediction has missing units")
 return float(np.mean((y-pred)**2))
def cluster_beta(y,terms,extra=None):
 # Models are exploratory and independently implemented using statsmodels.
 X=preds["baseline"].copy()
 vars=list(terms)
 for nm in vars:X[nm]=(d[nm].to_numpy(float)-d[nm].mean())/d[nm].std(ddof=0)
 if extra is not None:
  X[extra]=(d[extra].to_numpy(float)-d[extra].mean())/d[extra].std(ddof=0)
 fit=sm.OLS(y,sm.add_constant(X.to_numpy(float))).fit(
     cov_type="cluster",cov_kwds={"groups":d.huc8.to_numpy(),"use_correction":True})
 return {name:{"beta_per_1SD_source_log10_predictor":float(fit.params[i+1]),
               "HUC8_cluster_descriptive95":list(map(float,fit.conf_int(alpha=.05)[i+1]))}
         for i,name in enumerate(X.columns) if name in vars or name==extra}
def evaluate(outcome,models,conditional=None):
 y=d[outcome].to_numpy(float)
 mses={m:heldout(y,preds[m]) for m in models}
 basems=mses["baseline"]
 slopes=cluster_beta(y,["log_tn","log_tp"],conditional)
 return {"n_visits":len(d),"HUC8":d.huc8.nunique(),"response_mean":float(np.mean(y)),
         "response_sd":float(np.std(y,ddof=0)),
         "response_q25_median_q75":list(map(float,np.quantile(y,[.25,.5,.75]))),
         "out_of_HUC8_fivefold_MSE":mses,
         "NP_heldout_MSE_improvement_vs_baseline_pct":float(100*(basems-mses["plus_NP"])/basems),
         "posthoc_within_huc8_clustered_slopes":slopes,
         "not_temporal_or_sewage_causal":True}
report={
 "schema":"gesui_NR​​SA2324_GENUINE_PIGMENT_PER_AREAL_AFDM_MACRO_EXPLORATORY_v1",
 "original_source_frozen_2023_24_2129":True,
 "original_strict_qualifying_UIDs":len(d),
 "original_states":d.state.nunique(),
 "original_HUC8":d.huc8.nunique(),
 "strict_original_source_preserved":True,
 "question_selected_after_old_chla_and_afdm_results_known":True,
 "original_chla_area_unit":"ug/cm2",
 "original_AFDM_area_unit":"mg/cm2",
 "ratio_unit":"ug chlorophyll-a per mg ash-free dry biofilm material",
 "ecological_meaning_limit":"Pigment content relative to AFDM, not primary productivity, autotroph mass fraction, or energy export",
 "no_2023_benthic_mmi_or_count_attached":True,
 "not_EPA_wastewater_or_DMR_NP_dose":True,
 "no_pollination_food_web_function_claim":True,
 "descriptive_same_visit": {
   "AFDM_vs_chla_source_log10_Spearman":float(spearmanr(d.log_afdm,d.log_chla).statistic),
   "pigment_per_mg_median":float(d.pigment_ug_per_mg.median()),
   "pigment_per_mg_q25_q75":list(map(float,d.pigment_ug_per_mg.quantile([.25,.75])))
 },
 "responses":{
   "log10_pigment_per_AFDM":evaluate("log_pigment_per_afdm",["baseline","plus_NP"]),
   "log10_Chla_original":evaluate("log_chla",["baseline","plus_NP","plus_NP_AFDM"],"log_afdm"),
   "log10_AFDM_original":evaluate("log_afdm",["baseline","plus_NP"]),
 },
 "notes":[
  "Ratio slopes algebraically linked to Chla and AFDM. The ratio result is a NEW estimand on known old data, not an independent replication.",
  "Adjusted Chla model additionally conditions on measured AFDM, not a causal mediation regression.",
  "Original NRSA native UID visits include repeat physical SITE_ID; HUC8 spatial grouped validation and geographic confounding remain.",
  "MMI could not be independently obtained from official original EPA links (source QA separately).",
  "All site-level and measured biological numeric rows remain ephemeral; only aggregated results exported."
 ]}
OUT.write_text(json.dumps(report,indent=2,sort_keys=True))
print("NRSA2324_PIGMENT_AFDM_SOURCE_RESULT")
print(json.dumps({"n":len(d),"HUC8":d.huc8.nunique(),"Spearman_AFDM_CHLA":report["descriptive_same_visit"]["AFDM_vs_chla_source_log10_Spearman"],
 "results":{k:{"mse_gain_NP_pct":v["NP_heldout_MSE_improvement_vs_baseline_pct"],
                "slopes":v["posthoc_within_huc8_clustered_slopes"]} for k,v in report["responses"].items()}},sort_keys=True))
