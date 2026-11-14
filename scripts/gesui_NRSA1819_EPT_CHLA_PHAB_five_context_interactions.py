#!/usr/bin/env python3
"""Posthoc, source-frozen national EPA NRSA 2018/19 EPT–periphyton Chla
interactions with real mid-channel canopy, embeddedness, vegetation and width.

Population strictly original n946 native UID/603 HUC8. All 5 moderators are
precommitted in docs before testing signs. Full exposure and physical covariates
forced in every model. Whole-HUC8 geographic heldout Ridge, HUC8 clustered OLS
Holm multiple interaction contrast. NO municipal effluent, intervention,
insect emergence flux, biomass or ecosystem-service causal claims.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_numeric"
ORIGINAL=ROOT/"NRSA1819_SAME_UID_BIOTIC_TRIAD_PHAB_NP_REAL_RESULTS.json"
PANEL=ROOT/"AUDIT_ONLY_2018_19_EPT_Chla_PHAB_SOURCE_PANEL.csv"
src=pd.read_csv(PANEL,dtype={"UID":str,"state":str,"eco9":str,"huc8":str,"site_id":str})
old=json.loads(ORIGINAL.read_text())
if old["schema"]!="gesui_NRSA1819_TRUE_SAME_UID_EPT_Chla_TAXONRICHNESS_PHAB_NP_OBSERVATIONAL_MACRO_v1":
 raise RuntimeError("Original source model/field schema changed")
if len(src)!=1221 or src.UID.duplicated().any():raise RuntimeError("Source 1221 native EPA UID changed")
COHORT=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
        "log_tn","log_tp","log_chla","ept_fraction","total_benthic_subsample_count",
        "XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width","XCDENMID","XEMBED"]
quality=(src.chla_flag.fillna("").eq("")&
         src.tn_mgL_rawflag.fillna("").eq("")&
         src.tp_mgL_rawflag.fillna("").eq(""))
d=src[quality].dropna(subset=COHORT).copy()
d=d[(d.state!="")&(d.eco9!="")&(d.huc8!="")].sort_values("UID").reset_index(drop=True)
if len(d)!=946 or d.huc8.nunique()!=603 or d.state.nunique()!=47 or d.site_id.nunique()!=880:
 raise RuntimeError("Original prechosen dense PHab strict source n946/HUC8=603/47 states/880 SITE_ID changed")
n=len(d)
MODS=["XCDENMID","XEMBED","XPCAN","XGVEG","log_width"]
CONTROL=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
         "log_tn","log_tp","XPCAN","XGVEG","XFC_OHV","W1_HALL",
         "log_width","XCDENMID","XEMBED","log_chla"]
if not all(k in d for k in CONTROL+MODS):raise RuntimeError("Original PHab source predictors incomplete")
# One-hot state and ecoregion always retained (same as frozen PR77).
onehot=pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)
Xmain=pd.concat([d[CONTROL].astype(float),onehot],axis=1)
Y=d.ept_fraction.to_numpy(float)
groups=d.huc8.astype(str).to_numpy()
fold=np.asarray([int(hashlib.sha256(str(g).encode()).hexdigest()[:8],16)%5 for g in groups])
if min(np.bincount(fold,minlength=5))<100:raise RuntimeError("HUC8 fixed 5fold support failure")
m_chla=float(d.log_chla.mean());sd_chla=float(d.log_chla.std(ddof=0))
if sd_chla<=1.e-10:raise RuntimeError("Source Chla effect variation collapsed")
qmod={name:{"mean":float(d[name].mean()),"sd":float(d[name].std(ddof=0)),
            "q25":float(d[name].quantile(.25)),"q75":float(d[name].quantile(.75))}
      for name in MODS}
if any(x["sd"]<=1.e-10 for x in qmod.values()):raise RuntimeError("Source habitat modifier variation collapsed")
zchla=(d.log_chla.to_numpy(float)-m_chla)/sd_chla

def add_extras(columns,train,moderators=()):
 # Original nuisance predictors + pre-specified curvilinear Chla + interactions.
 X=Xmain.to_numpy(float).copy()
 m=X[train].mean(axis=0);s=X[train].std(axis=0)
 z=(X-m)/np.where(s>1.e-10,s,1)
 c=z[:,CONTROL.index("log_chla")]
 extras=[c*c]
 for key in moderators:
  zi=z[:,CONTROL.index(key)]
  extras.append(c*zi)
 z=np.column_stack([z]+extras)
 # Ridge penalty sees each added product standardized on TRAIN only.
 mean=z[train].mean(axis=0)
 sd=z[train].std(axis=0)
 final=(z-mean)/np.where(sd>1e-10,sd,1)
 return final

# MAIN already has 2018/19 fixed stratified mixed geography. 'base' has all
# physical and nutrients + Chla, Q adds Chla^2, I adds exactly one moderator.
def predict(label,mods=(),curve=True,subset=None):
 use=np.ones(n,dtype=bool) if subset is None else subset
 if not use.all():raise RuntimeError("Primary CV must use identical complete cohort")
 pred=np.full(n,np.nan)
 for k in range(5):
  train=fold!=k;test=~train
  x=add_extras( CONTROL,train,mods if curve else ()) if curve else None
  if not curve:
   arr=Xmain.to_numpy(float)
   mu=arr[train].mean(0);sd=arr[train].std(0)
   x=(arr-mu)/np.where(sd>1e-10,sd,1)
  ztrain=np.column_stack([np.ones(train.sum()),x[train]])
  ztest=np.column_stack([np.ones(test.sum()),x[test]])
  ridge=np.eye(ztrain.shape[1])*10
  ridge[0,0]=0
  coefs=np.linalg.solve(ztrain.T@ztrain+ridge,ztrain.T@Y[train])
  pred[test]=ztest@coefs
 if not np.isfinite(pred).all():raise RuntimeError("Missing heldout predictions")
 return pred

P={"B":predict("B",curve=False),"Q":predict("Q",curve=True)}
for name in MODS:
 P["I_"+name]=predict("I_"+name,(name,),curve=True)
P["Iall"]=predict("Iall",tuple(MODS),curve=True)
mse={k:float(np.mean((Y-p)**2)) for k,p in P.items()}
def gain(base,alt):return float(100*(mse[base]-mse[alt])/mse[base])
def boot(b,a):
 u,inv=np.unique(groups,return_inverse=True)
 e0=np.bincount(inv,weights=(Y-P[b])**2,minlength=len(u))
 e1=np.bincount(inv,weights=(Y-P[a])**2,minlength=len(u))
 rng=np.random.default_rng(20261114)
 ids=rng.integers(0,len(u),size=(2000,len(u)))
 den=e0[ids].sum(1);num=e1[ids].sum(1)
 arr=100*(den-num)/den
 return [float(v) for v in np.quantile(arr,[.025,.5,.975])]

# Original full sample descriptive OLS, not an independent intervention.
z=(d.log_chla.to_numpy(float)-m_chla)/sd_chla
XOLS=Xmain.copy()
XOLS["logChla_z_squared"]=z*z
for name in MODS:
 zz=(d[name].to_numpy(float)-qmod[name]["mean"])/qmod[name]["sd"]
 XOLS["logChla_z_x_"+name]=z*zz
def fit_ols(fields):
 X=XOLS[["logChla_z_squared"]+[]] if False else XOLS.drop(columns=[
  ("logChla_z_x_"+j) for j in MODS if j not in fields
 ])
 model=sm.OLS(Y,sm.add_constant(X.to_numpy(float),has_constant="add")).fit(
  cov_type="cluster",cov_kwds={"groups":groups,"use_correction":True})
 return X,model
results={}
for name in MODS:
 X,mod=fit_ols((name,))
 col="logChla_z_x_"+name
 idx=X.columns.get_loc(col)+1
 v=float(mod.params[idx]);se=float(mod.bse[idx])
 p=float(2*norm.sf(abs(v/se))) if se>0 else 1.
 cis=list(map(float,mod.conf_int()[idx]))
 base_logChla_index=X.columns.get_loc("log_chla")+1
 base_sq_index=X.columns.get_loc("logChla_z_squared")+1
 mean_slope=float(mod.params[base_logChla_index]*sd_chla)
 # Squared curvature derivative zero at original Chla mean after z-centering.
 q1=(qmod[name]["q25"]-qmod[name]["mean"])/qmod[name]["sd"]
 q3=(qmod[name]["q75"]-qmod[name]["mean"])/qmod[name]["sd"]
 results[name]={
  "interaction_beta_per_1SD_Chla_and_1SD_environment":v,
  "HUC8_cluster_descriptive95":cis,
  "two_sided_cluster_normal_approx_p_posthoc":p,
  "source_mod_q25":qmod[name]["q25"],"source_mod_q75":qmod[name]["q75"],
  "predicted_EPT_per_1SD_logChla_at_low_mod":mean_slope+v*q1,
  "predicted_EPT_per_1SD_logChla_at_high_mod":mean_slope+v*q3,
  "heldout_MSE_extra_gain_over_curved_Q_pct":gain("Q","I_"+name),
  "HUC8_bootstrap_2000_fixedOOF95_extra_gain_pct":boot("Q","I_"+name),
 }
# Holm 5 interaction tests (each same n, full nuisance and Chla quadratic).
ranking=sorted(MODS,key=lambda n:results[n]["two_sided_cluster_normal_approx_p_posthoc"])
prior=0.
for i,name in enumerate(ranking):
 raw=results[name]["two_sided_cluster_normal_approx_p_posthoc"]
 prior=max(prior,min(1.,raw*(len(MODS)-i)))
 results[name]["Holm_5_family_adjusted_p_posthoc"]=float(prior)
modall=fit_ols(MODS)[1]
joint={
 "five_terms_heldout_gain_over_curve_pct":gain("Q","Iall"),
 "HUC8_bootstrap_fixedOOF95":boot("Q","Iall"),
 "global_interaction_coefficient_joint_test_not_preregistered":True,
}
# Explanatory only: source original one XCDENMID>100 might affect the
# whole-sample Chla×canopy-density coefficient. Audit its event count.
quality={
 "original_XCDENMID_gt100_rows":int((d.XCDENMID>100).sum()),
 "original_XGVEG_gt1_rows":int((d.XGVEG>1).sum()),
 "original_native_site_visits":len(d),
 "unique_physical_site_ID":d.site_id.nunique(),
 "states":d.state.nunique(),"HUC8":d.huc8.nunique(),
 "fivefold_groups":np.bincount(fold,minlength=5).tolist(),
 "initial_2018_year_only":int((d.year==2018).sum()),
 "initial_2019_year_only":int((d.year==2019).sum())
}
report={
 "schema":"gesui_NR​​SA1819_EPT_CHLA_fivePHAB_moderators_preattached_posthoc_v1",
 "source_original_benthic_chla_NP_PHAB_2018_sample_locked":True,
 "source_old_PR77_EPT_Chla_PHAB_QA_success":True,
 "original_parent_2018_EPT_Chla_PHAB_json":str(ORIGINAL.name),
 "question_and_holm_family_frozen_before_numeric_results":True,
 "exploratory_posthoc_not_confirmatory":True,
 "original_PHAB_measurement_not_external_solar_radiation":True,
 "no_wastewater_effluent_source_or_service_attribution":True,
 "not_species_richness_or_adult_insect_biomass":True,
 "source_quality":quality,
 "cohort":len(d),"HUC8":d.huc8.nunique(),
 "MSE_HUC8_disjoint_5fold":mse,
 "curve_Chla_squared_gain_over_linear_B_pct":gain("B","Q"),
 "single_moderators":results,
 "joint_five_interactions":joint,
 "interpretation_limits":[
  "Interaction slopes depend on full environmental adjustment; they are not treatments.",
  "Holm correction protects nominal familywise multiple testing but cannot undo posthoc hypothesis selection.",
  "Observed canopy is a spatial shading proxy, not RSDS or diurnal light exposure.",
  "Same-source observational cross-layer association is not biological mediation or ecosystem service.",
  "Grouped holdouts have fixed grouping; bootstrap is only a fixed out-of-fold error comparison."]
}
dest=ROOT/"NRSA1819_EPT_CHLA_PHAB_INTERACTION_SOURCE_SIGNED_EXPLORATORY.json"
dest.write_text(json.dumps(report,indent=2,sort_keys=True))
print("EPT_CHLA_PHAB_INTERACTION_FULL_RESULT_START")
print(json.dumps({"QA":quality,"mse":mse,"curvature_gain":report["curve_Chla_squared_gain_over_linear_B_pct"],
                   "interaction_results":results,"joint":joint},sort_keys=True))
print("EPT_CHLA_PHAB_INTERACTION_FULL_RESULT_END")
