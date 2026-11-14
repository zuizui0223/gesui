#!/usr/bin/env python3
"""INDEPENDENT FIELD CYCLE: EPA NRSA 2013/14 genuine BERW benthic EPT counted
fraction vs areal PCHL pigment and areal PBIO AFDM, original XEMBED (%).
Uses original exact-native UID, source-first 1203 probability BERW UIDs,
true measured river TN/TP, canopy and physical covariates. No data from
NRSA 2018/19 used to estimate 2013 coefficients. Pre-frozen test imported
from original source protocol. No sewage discharge, insect emergence or service.
"""
from __future__ import annotations
import hashlib,json,os,pathlib,urllib.request,time
import numpy as np,pandas as pd
from scipy.stats import spearmanr,norm
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1314_independent_numeric"
ROOT.mkdir(parents=True,exist_ok=True)
E=pathlib.Path(os.environ["NRSA1314_SEMANTIC_SOURCE_QA"])
source=json.loads(E.read_text())
if source["schema"]!="NRSA1314_ORIGINAL_EPA_2019_SIX_TABLE_PRENUMERIC_SEMANTIC_CONTRACT_v1":
 raise RuntimeError("Independently locked source 2013 semantic gate changed")
G=source["original_source_gate"]
if not G["eligible_for_numeric_if_n350_HUC200_states20"] or len(G["source_eligible_native_UIDs"])!=1203:
 raise RuntimeError("Original source BERW PROB 1203/754 HUC8/48 states gate changed")
F=source["original_six_source_receipts_sha256"]
HOST="https://www.epa.gov/sites/default/files/2019-04/"
FILES={
 "sites":"nrsa1314_siteinformation_wide_04292019.csv",
 "water":"nrsa1314_widechem_04232019.csv",
 "chla":"nrsa1314_widepchl_04232019.csv",
 "afdm":"nrsa1314_widepbio_04232019.csv",
 "phab":"nrsa1314_phabmed_04232019.csv",
 "benthic":"nrsa1314_bentcnts_04232019.csv"}
def load(name):
 f=ROOT/FILES[name]
 for i in range(4):
  try:
   request=urllib.request.Request(HOST+FILES[name],headers={"User-Agent":"gesui-NRSA1314-independent-raw-EPA-numeric/20261114"})
   with urllib.request.urlopen(request,timeout=120) as r,f.open("wb") as out:
    while buf:=r.read(2**20):out.write(buf)
   if hashlib.sha256(f.read_bytes()).hexdigest()!=F[name]:
    raise RuntimeError("Original 2013 EPA source bytes changed "+name)
   break
  except Exception:
   f.unlink(missing_ok=True)
   if i==3:raise
   time.sleep(2**i)
 d=pd.read_csv(f,dtype=str,encoding="latin1",keep_default_na=False,low_memory=False)
 f.unlink();d.columns=d.columns.str.upper().str.strip()
 return d
raw={k:load(k) for k in FILES}
allowed=set(G["source_eligible_native_UIDs"])
if len(allowed)!=1203:raise RuntimeError("Original source-blind 2013 sample changed")
for k in ("sites","water","chla","afdm","phab"):
 if raw[k].UID.duplicated().any():raise RuntimeError("Actual 2013 original native UID nonunique in "+k)
 raw[k]=raw[k][raw[k].UID.isin(allowed)].set_index("UID").sort_index()
 if len(raw[k])!=1203:raise RuntimeError("2013 source native UID missing in "+k)
ref=raw["sites"]
if not ref.SITETYPE.eq("PROB").all():
 raise RuntimeError("Original 2013 source-only probability site criteria changed")
source_id_discordance={}
for k in ("water","chla","afdm","phab"):
 n=int((ref.SITE_ID!=raw[k].SITE_ID).sum()) if "SITE_ID" in raw[k] else None
 source_id_discordance[k]=n
 if n:
  raise RuntimeError("True original EPA 2013 same site UID mapping has SITE_ID conflicts "+k+":"+str(n))
# Strict same DATE_COL if present (same visit; source original field dates)
date_discordance={}
for k in ("water","chla","afdm","phab"):
 if "DATE_COL" in raw[k]:
  t=pd.to_datetime(ref.DATE_COL,errors="coerce",format="mixed")
  v=pd.to_datetime(raw[k].DATE_COL,errors="coerce",format="mixed")
  bad=t.notna()&v.notna()&t.dt.normalize().ne(v.dt.normalize())
  date_discordance[k]=int(bad.sum())
  fmt_only=int((ref.DATE_COL.astype(str)!=raw[k].DATE_COL.astype(str)).sum()-bad.sum())
  print("NRSA1314_ORIGINAL_DATE_ALIGNMENT",k,
        "true_unmatched",date_discordance[k],"format_only_or_missing",fmt_only,
        "original_sample_date_pairs",list(zip(ref.DATE_COL.loc[bad].head(4).tolist(),raw[k].DATE_COL.loc[bad].head(4).tolist())),flush=True)
  if date_discordance[k]>0:raise RuntimeError("Native original SITE_ID same but original visit DATE_COL conflicts: "+k)
B=raw["benthic"]
B=B[B.UID.isin(allowed)&B.SAMPLE_TYPE.str.upper().eq("BERW")].copy()
B["N"]=pd.to_numeric(B.TOTAL300,errors="coerce")
if (B.N.dropna()<0).any():raise RuntimeError("Negative BERW individual counts not permitted")
B=B[B.N.gt(0)].copy()
B["EPT"]=B.ORDER.str.upper().str.strip().isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"))
if B.empty:raise RuntimeError("No positive original 2013 BERW sample source")
by=B.groupby("UID")
btot=by.N.sum()
eptcnt=(B.N*B.EPT.astype(float)).groupby(B.UID).sum()
ept_share=(eptcnt/btot).reindex(ref.index)
taxonID=B.groupby("UID").TAXA_ID.nunique().reindex(ref.index)
if not ept_share.dropna().between(0,1).all():
 raise RuntimeError("Invalid original 2013 counted EPT share")
def number(frame,column):
 return pd.to_numeric(frame[column],errors="coerce")
W=raw["water"];C=raw["chla"];A=raw["afdm"];H=raw["phab"]
if set(W.NTL_RESULT_UNITS.str.upper().str.strip())!={"MG/L"}:raise RuntimeError("Original 2013 total N units changed")
if set(W.PTL_RESULT_UNITS.str.upper().str.strip())!={"UG/L"}:raise RuntimeError("Original 2013 total P units changed")
if not set(C.CHLA_RESULT_UNITS.str.upper().str.strip()).issubset({"UG/CM2",""}):
 raise RuntimeError("Original 2013 periphyton Chla is not areal ug/cm2")
if not set(A.AFDM_RESULT_UNITS.str.upper().str.strip()).issubset({"MG/CM2","MG",""}):
 raise RuntimeError("Original 2013 periphyton biomass unit schema changed")
if not set(H.columns).issuperset({"XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","XWIDTH"}):
 raise RuntimeError("Original 2013 PHab common physical columns changed")
D=pd.DataFrame(index=ref.index)
D["UID"]=D.index
for old,new in [("SITE_ID","site_id"),("PSTL_CODE","state"),("AG_ECO9","eco9"),("HUC8","huc8")]:
 D[new]=ref[old].astype(str).str.strip()
D["lat"]=number(ref,"LAT_DD83");D["lon"]=number(ref,"LON_DD83")
D["stream_order"]=number(ref,"STRAH_ORD")
DATE=pd.to_datetime(ref.DATE_COL,errors="coerce")
D["year"]=DATE.dt.year
D["doy"]=DATE.dt.dayofyear
D["season_sin"]=np.sin(2*np.pi*D.doy/365.2425)
D["season_cos"]=np.cos(2*np.pi*D.doy/365.2425)
D["ept_fraction"]=ept_share
D["total_benthic_counted"]=btot.reindex(D.index)
D["unique_original_taxonID_n"]=taxonID
D["logcond"]=np.log10(number(W,"COND_RESULT").where(lambda v:v>0))
D["logTN"]=np.log10(number(W,"NTL_RESULT").where(lambda v:v>0))
D["logTP"]=np.log10((number(W,"PTL_RESULT")/1000).where(lambda v:v>0))
D["logchla"]=np.log10(number(C,"CHLA_RESULT").where(lambda v:v>0))
D["logafdm"]=np.log10(number(A,"AFDM_RESULT").where(lambda v:v>0))
for v in ("XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","XWIDTH"):
 D[v]=number(H,v)
D["logwidth"]=np.log10(D.XWIDTH.where(D.XWIDTH>0))
D["strict_original_water_flags"]=W.NTL_NARS_FLAG.str.strip().eq("") & W.PTL_NARS_FLAG.str.strip().eq("")
D["strict_original_chla_flag"]=C.CHLA_NARS_FLAG.str.strip().eq("")
D["strict_original_AFDM_flag"]=A.AFDM_NARS_FLAG.str.strip().eq("")
D["original_afdm_areal_unit"]=A.AFDM_RESULT_UNITS.str.upper().str.strip().eq("MG/CM2")
D["original_chla_areal_unit"]=C.CHLA_RESULT_UNITS.str.upper().str.strip().eq("UG/CM2")
cov=["lat","lon","year","season_sin","season_cos","stream_order","logcond",
     "logTN","logTP","XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","logwidth","logchla","ept_fraction"]
use=(D.strict_original_water_flags & D.strict_original_chla_flag & D.original_chla_areal_unit)
sample=D[use].dropna(subset=cov).copy()
sample=sample[(sample.huc8!="")&(sample.state!="")&(sample.eco9!="")&(sample.year.isin([2013,2014]))]
sample=sample.sort_index()
if len(sample)<350 or sample.huc8.nunique()<200 or sample.state.nunique()<20:
 raise RuntimeError(f"Independent actual 2013 genuine same-visit numeric gate failed n{len(sample)} HUC8{sample.huc8.nunique()} states{sample.state.nunique()}")
# Exact same available PHab common columns as the harmonized 2018 study,
# with no 2013-only XGV substitution for 2018 XGVEG.
physical=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","logwidth"]
geog=["lat","lon","year","season_sin","season_cos","stream_order","logcond","logTN","logTP"]+physical
def X_base(d,exposure):
 base=pd.concat([d[geog].astype(float),pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 base[exposure]=d[exposure].astype(float)
 return base
def ols_cluster(d,X,Y,exp,mod):
 z0=(d[exp].to_numpy(float)-d[exp].mean())/d[exp].std(ddof=0)
 z1=(d[mod].to_numpy(float)-d[mod].mean())/d[mod].std(ddof=0)
 full=X.assign(quadratic=z0*z0,interaction=z0*z1)
 fit=sm.OLS(Y,sm.add_constant(full.to_numpy(float),has_constant="add")).fit(
    cov_type="cluster",cov_kwds={"groups":d.huc8.to_numpy(),"use_correction":True})
 i=len(full.columns)
 return {"beta":float(fit.params[i]),"ci":list(map(float,fit.conf_int()[i])),
         "raw_normal_p":float(2*norm.sf(abs(fit.params[i]/fit.bse[i])))}
def CV(d,X,Y,exp,mod,include_interaction):
 a=X.to_numpy(float);res=np.full(len(Y),np.nan)
 cv=np.asarray([int(hashlib.sha256(str(g).encode()).hexdigest()[:8],16)%5 for g in d.huc8])
 expind=X.columns.get_loc(exp);modind=X.columns.get_loc(mod)
 for k in range(5):
  tr=cv!=k;te=~tr
  if tr.sum()<120 or te.sum()<20:raise RuntimeError("2013 within-cycle fold support insufficient")
  ss=StandardScaler().fit(a[tr]);z=ss.transform(a)
  extra=[z[:,expind]**2]
  if include_interaction:extra.append(z[:,expind]*z[:,modind])
  Z=np.column_stack([z]+extra)
  zz=StandardScaler().fit(Z[tr]).transform(Z)
  model=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(zz[tr],Y[tr])
  res[te]=model.predict(zz[te])
 if not np.isfinite(res).all():raise RuntimeError("Prediction leakage / missing holdouts")
 return res
def boot(y,q,i,groups):
 unq,inv=np.unique(groups,return_inverse=True)
 w0=np.bincount(inv,weights=(y-q)**2,minlength=len(unq))
 w1=np.bincount(inv,weights=(y-i)**2,minlength=len(unq))
 rng=np.random.default_rng(20261114)
 samples=rng.integers(0,len(unq),size=(2000,len(unq)))
 a=w0[samples].sum(1);b=w1[samples].sum(1)
 return list(map(float,np.quantile(100*(a-b)/a,[.025,.5,.975])))
def model(d,exposure="logchla",response="ept_fraction"):
 Y=d[response].to_numpy(float)
 X=X_base(d,exposure)
 b=ols_cluster(d,X,Y,exposure,"XEMBED")
 pred_q=CV(d,X,Y,exposure,"XEMBED",False)
 pred_i=CV(d,X,Y,exposure,"XEMBED",True)
 mse0=float(np.mean((Y-pred_q)**2));mse1=float(np.mean((Y-pred_i)**2))
 return {"n":len(d),"HUC8":d.huc8.nunique(),"states":d.state.nunique(),
         "interaction":b,"MSE_Q":mse0,"MSE_interaction":mse1,
         "HUC8_OOS_MSE_gain_pct":float(100*(mse0-mse1)/mse0),
         "fixed_OOF_HUC8_2000bootstrap95_gain":boot(Y,pred_q,pred_i,d.huc8),
         "observed_Spearman_EPT_vs_biofilm":float(spearmanr(d.ept_fraction,d[exposure]).statistic)}
m=model(sample)
mass=sample[sample.strict_original_AFDM_flag & sample.original_afdm_areal_unit].dropna(subset=["logafdm"]).copy()
m_mass=model(mass,"logafdm","ept_fraction") if len(mass)>=350 and mass.huc8.nunique()>=200 else {"gate":"NO_AREAL_MASS_SAMPLE_SUPPORT","n":len(mass)}
# Ephemeral raw field panel for independent review within the same GH runner;
# deliberately excluded from uploaded artifact and git repository.
sample.to_csv(ROOT/"AUDIT_ONLY_NRSA1314_BERW_EPT_CHLA_SOURCE_PANEL.csv",index_label="UID")
mass.to_csv(ROOT/"AUDIT_ONLY_NRSA1314_AREAL_AFDM_SOURCE_PANEL.csv",index_label="UID")
# Periphyt-on AFDM response vs ambient N/P and sediment is not response
# in a model deliberately optimized for positive sign; only status and a
# directly measured same-visit correlation reported as secondary.
qa={"original_source_wadeable_probability_UID":len(allowed),
 "strict_EPTChla_NP_PHAB_complete":len(sample),
 "strict_positive_areal_AFDM_same_source":len(mass),
 "original_source_year_counts":{str(k):int(v) for k,v in sample.year.value_counts().items()},
 "N_HUC8":sample.huc8.nunique(),"states":sample.state.nunique(),
 "physical_sites":sample.site_id.nunique(),
 "source_uid_join_site_conflicts":source_id_discordance,
 "source_uid_join_date_conflicts":date_discordance,
 "sample_benthos_positive_record_n":len(B),
 "sample_BERW_positive_native_UIDs":B.UID.nunique(),
 "original_sample_EPT_share_median":float(sample.ept_fraction.median()),
 "original_chla_log_median":float(sample.logchla.median()),
 "original_areal_AFDM_log_median":float(mass.logafdm.median()) if len(mass) else None,
 "source_AFDM_v_CHLA_log_spearman":float(spearmanr(mass.logafdm,mass.logchla).statistic) if len(mass)>=10 else None}
out={"schema":"gesui_EPA_NR​​SA1314_INDEPENDENT_2013_wade_BERW_EPT_CHLA_AFDM_EMBEDDEDNESS_v1",
 "original_preselected_1203_UID_source_manifest_SHA256":hashlib.sha256(E.read_bytes()).hexdigest(),
 "source_original_six_EPA_file_SHA":F,"genuine_independent_2013_14_survey_cycle":True,
 "same_native_UID_fields_only":True,"not_2018_field_replication":True,
 "different_physical_XGVEG_missing_so_harmonized_six_common_variables_only":True,
 "no_2013_outcome_based_source_selection":True,
 "posthoc_original_2018_19_hypothesis_known":True,
 "no_WWTP_discharge_exposure_or_ecosystem_services":True,
 "quality":qa,
 "2013_14_biofilmChla_EPT_embed_primary":m,
 "2013_14_arealAFDM_EPT_embed_secondary":m_mass}
p=ROOT/"NRSA1314_INDEPENDENT_BERW_AREAL_PCHL_AFDM_EMBED_RESULT.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA1314_INDEPENDENT_ECOLOGICAL_CYCLE_RESULT",json.dumps({k:v for k,v in out.items() if k in ("quality","2013_14_biofilmChla_EPT_embed_primary","2013_14_arealAFDM_EPT_embed_secondary")},sort_keys=True))
