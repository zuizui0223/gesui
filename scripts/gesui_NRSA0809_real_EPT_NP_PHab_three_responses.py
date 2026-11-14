#!/usr/bin/env python3
"""Third independent survey EPA NRSA 2008–09: SOURCE-LOCKED actual benthic
EPT percent, EPT_NTAX and benthic MMI vs real measured stream TN/TP, PHab.
No 2008 Chla data; no invented substitute/biomass or wastewater discharge.
"""
import hashlib,json,os,pathlib,urllib.request,time
import numpy as np,pandas as pd
import statsmodels.api as sm
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

R=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_real_numeric"
R.mkdir(parents=True,exist_ok=True)
SRC=json.loads(pathlib.Path(os.environ["NRSA0809_FROZEN_SOURCE_QA"]).read_text())
if SRC["schema"]!="gesui_NRSA0809_original_wadeable_CHEMW_BERWW_PHAB_MMI_pre_numeric_v1":
 raise RuntimeError("Original original-source-only protocol not verified")
G=SRC["original_wadeable_primary_source_gate"]
uids=set(G["original_UIDs"])
if len(uids)!=1299 or G["original_source_HUC8"]!=796 or G["original_source_states"]!=48:
 raise RuntimeError("Pre-outcome 2008 source 1299/796/48 frozen eligible UID population changed")
FILES=SRC["original_EPA_CSVs_SHA256"]
URL={
"siteinfo":"https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.csv",
"water":"https://www.epa.gov/sites/default/files/2015-09/chem.csv",
"phab":"https://www.epa.gov/sites/default/files/2015-09/phabmed.csv",
"bentcount":"https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv",
"bentmmi":"https://www.epa.gov/sites/default/files/2015-09/bentcond.csv",
}
def load(k):
 p=R/(k+".csv")
 for attempt in range(3):
  try:
   req=urllib.request.Request(URL[k],headers={"User-Agent":"gesui-nrsa0809-thirdcycle-actual-raw-nutrient-benthos"})
   with urllib.request.urlopen(req,timeout=120) as w,p.open("wb") as f:
    while b:=w.read(1<<20):f.write(b)
   if hashlib.sha256(p.read_bytes()).hexdigest()!=FILES[k]:
    raise RuntimeError("EPA 2008 source no longer matches SHA "+k)
   break
  except Exception:
   p.unlink(missing_ok=True)
   if attempt==2:raise
   time.sleep(2**attempt)
 x=pd.read_csv(p,dtype=str,encoding="latin1",keep_default_na=False,low_memory=False)
 p.unlink()
 x.columns=x.columns.str.upper().str.strip()
 return x
s=load("siteinfo")
w=load("water")
h=load("phab")
b=load("bentcount")
m=load("bentmmi")
if s.UID.duplicated().any() or h.UID.duplicated().any() or m.UID.duplicated().any():
 raise RuntimeError("Official 2008 source repeated original UID in one-to-one table")
ww=w[w.SAMPLE_TYPE.str.upper().str.strip().eq("CHEMW")&
     w.SAMPLE_CAT.str.upper().str.strip().eq("P")].copy()
bb=b[b.SAMPLE_TYPE.str.upper().str.strip().eq("BERWW")&
     b.SAMPLE_CAT.str.upper().str.strip().eq("P")].copy()
if ww.UID.duplicated().any():raise RuntimeError("Original CHEMW primary native UID repeated")
if not uids<=set(bb.UID):
 raise RuntimeError("Original native source probability wadeable benthic UID no longer present")
s=s[s.UID.isin(uids)].set_index("UID").sort_index()
ww=ww[ww.UID.isin(uids)].set_index("UID").sort_index()
h=h[h.UID.isin(uids)].set_index("UID").sort_index()
m=m[m.UID.isin(uids)].set_index("UID").sort_index()
if len(s)!=1299 or len(ww)!=1299 or len(h)!=1299 or len(m)!=1299:
 raise RuntimeError("Frozen 2008 source's complete five original native UID population changed")
for name,x in (("CHEMW",ww),("PHab",h),("BENT MMI",m)):
 if not (x.SITE_ID.eq(s.SITE_ID)).all():
  raise RuntimeError("EPA 2008 native original same UID but different SITE_ID in "+name)
# Original DATE_COL, may be different text formatting or PHab value empty.
date_q={}
for name,x in (("CHEMW",ww),("PHab",h),("BENT MMI",m)):
 if "DATE_COL" not in x:continue
 d0=pd.to_datetime(s.DATE_COL,errors="coerce",format="mixed")
 d1=pd.to_datetime(x.DATE_COL,errors="coerce",format="mixed")
 valid=d0.notna()&d1.notna()
 mismatch=valid&d0.dt.normalize().ne(d1.dt.normalize())
 date_q[name]={"both_dates_parseable":int(valid.sum()),"actual_date_discrepancy":int(mismatch.sum())}
 if int(mismatch.sum())>0:
  raise RuntimeError("Source native UID original visit date mismatch for "+name)
# Require published BENT indicator from actual BERWW method and observed sample.
mtype=m.BENT_SAMPLE_TYPE.str.upper().str.strip().eq("BERWW")
observed=m.SAMPLED_BENT.str.upper().str.strip().eq("YES")
source_mmi_ok=mtype&observed
D=pd.DataFrame(index=s.index)
D["state"]=s.STATE.astype(str).str.strip()
D["eco9"]=s.AGGR_ECO9_2015.astype(str).str.strip()
D["huc8"]=s.HUC8.astype(str).str.strip()
D["site_id"]=s.SITE_ID.astype(str)
D["lat"]=pd.to_numeric(s.LAT_DD83,errors="coerce")
D["lon"]=pd.to_numeric(s.LON_DD83,errors="coerce")
D["year"]=pd.to_numeric(s.YEAR,errors="coerce")
D["stream_order"]=pd.to_numeric(s.STRAHLERORDER,errors="coerce")
date=pd.to_datetime(s.DATE_COL,errors="coerce",format="mixed")
D["dayofyear"]=date.dt.dayofyear
D["season_sin"]=np.sin(2*np.pi*D.dayofyear/365.2425)
D["season_cos"]=np.cos(2*np.pi*D.dayofyear/365.2425)
def numeric(x,k):return pd.to_numeric(x[k],errors="coerce")
D["log_cond"]=np.log10(numeric(ww,"COND").where(lambda x:x>0))
D["log_tn"]=np.log10((numeric(ww,"NTL")/1000).where(lambda x:x>0))
D["log_tp"]=np.log10((numeric(ww,"PTL")/1000).where(lambda x:x>0))
D["source_nutrients_unalerted"]=ww.NTL_ALERT.str.strip().eq("") & ww.PTL_ALERT.str.strip().eq("")
D["source_cond_unalerted"]=ww.COND_ALERT.str.strip().eq("")
D["MMI_record_wade_observed"]=source_mmi_ok
D["ept_pct_original"]=numeric(m,"EPT_PIND")
D["ept_fraction"]=D.ept_pct_original/100.0
D["ept_native_taxon_count"]=numeric(m,"EPT_NTAX")
D["log1p_EPT_taxon_count"]=np.log1p(D.ept_native_taxon_count.where(D.ept_native_taxon_count>=0))
D["benthic_MMI"]=numeric(m,"MMI_BENT")
phab=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","XWIDTH"]
for col in phab:D[col]=numeric(h,col)
D["log_width"]=np.log10(D.XWIDTH.where(D.XWIDTH>0))
# Source-population diagnostics, avoiding false zero.
ranges={}
for k in ["ept_pct_original","ept_native_taxon_count","benthic_MMI","XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","XWIDTH"]:
 v=D[k]
 ranges[k]={"n_numeric":int(v.notna().sum()),"min":float(v.min()) if v.notna().any() else None,
  "max":float(v.max()) if v.notna().any() else None,"invalid_negative":int((v<0).sum()),
  "above_100":int((v>100).sum())}
if (D.ept_pct_original.dropna()>100).any() or (D.ept_pct_original.dropna()<0).any():
 raise RuntimeError("Source EPA EPT percent should be within documented 0..100")
if (D.ept_native_taxon_count.dropna()<0).any():
 raise RuntimeError("Source EPT taxon count is negative, unsafe to model")
if (D.benthic_MMI.dropna()>100).any() or (D.benthic_MMI.dropna()<0).any():
 raise RuntimeError("Source benthic MMI score outside documented 0..100")
flags=D.source_nutrients_unalerted & D.source_cond_unalerted & D.MMI_record_wade_observed
baseline=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
          "log_tn","log_tp","XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
# Keep original EPT=0 where a genuine published biological zero was measured.
# SOURCE measurement completeness diagnostics, before fitting; not a
# license to loosen original pre-committed sample/quality gates based on signs.
print("NRSA0809_PRENUMERIC_ORIGINAL_SOURCE_QA",json.dumps({
  "source_1299_UIDs":len(D),
  "BENT_SAMPLE_TYPE":m.BENT_SAMPLE_TYPE.str.strip().value_counts().head(12).to_dict(),
  "SAMPLED_BENT":m.SAMPLED_BENT.str.strip().value_counts().head(12).to_dict(),
  "MMI_wade_yes":int(source_mmi_ok.sum()),
  "original_CHEMW_nutrient_no_alert":int(D.source_nutrients_unalerted.sum()),
  "source_COND_no_alert":int(D.source_cond_unalerted.sum()),
  "joint_all_quality_flags":int(flags.sum()),
  "all_covariate_numeric_nonmissing":{col:int(D[col].notna().sum()) for col in baseline+["ept_fraction","log1p_EPT_taxon_count","benthic_MMI"]},
  "joint_good_flags_and_all_numeric":int(D.loc[flags].dropna(subset=baseline+["ept_fraction","log1p_EPT_taxon_count","benthic_MMI"]).shape[0]),
  "original_measurement_ranges":ranges
 },sort_keys=True),flush=True)
Q=D[flags].dropna(subset=baseline+["ept_fraction","log1p_EPT_taxon_count","benthic_MMI"]).copy()
Q=Q[(Q.state!="")&(Q.eco9!="")&(Q.huc8!="")&Q.year.isin([2008,2009])].sort_index()
if len(Q)<300 or Q.huc8.nunique()<150 or Q.state.nunique()<20:
 raise RuntimeError(f"Precommitted original numeric support failure: n{len(Q)} HUC8{Q.huc8.nunique()} state{Q.state.nunique()}")
if Q.site_id.nunique()<300:raise RuntimeError("Repeated survey visits pseudo-replicated")
G=Q.huc8.to_numpy()
fold=np.asarray([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5 for x in G])
if min(np.bincount(fold,minlength=5))<30:raise RuntimeError("Fivefold spatial HUC8 holdout lacks support")
geo=pd.concat([
 Q[["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]].astype(float),
 pd.get_dummies(Q[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
phys=Q[["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]].astype(float)
npchem=Q[["log_tn","log_tp"]].astype(float)
M={
 "geo":geo,
 "geo+N_P":pd.concat([geo,npchem],axis=1),
 "geo+PHab":pd.concat([geo,phys],axis=1),
 "geo+PHab+N_P":pd.concat([geo,phys,npchem],axis=1)}
def cvpred(X,y):
 arr=X.to_numpy(float);pred=np.full(len(y),np.nan)
 for k in range(5):
  tr=fold!=k;te=~tr
  scaler=StandardScaler().fit(arr[tr])
  z=scaler.transform(arr)
  model=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(z[tr],y[tr])
  pred[te]=model.predict(z[te])
 if not np.isfinite(pred).all():raise RuntimeError("2008 original heldout prediction missing")
 return pred
def slope(X,y):
 xx=X.copy()
 for col in ("log_tn","log_tp"):
  xx[col]=(xx[col]-xx[col].mean())/xx[col].std(ddof=0)
 fit=sm.OLS(y,sm.add_constant(xx.to_numpy(float),has_constant="add")).fit(
     cov_type="cluster",cov_kwds={"groups":G,"use_correction":True})
 return {col:{"beta_per_source_1SD_log10":float(fit.params[xx.columns.get_loc(col)+1]),
   "HUC8_cluster_descriptive95":list(map(float,fit.conf_int()[xx.columns.get_loc(col)+1]))}
   for col in ("log_tn","log_tp")}
def boot(y,a,b):
 u,iv=np.unique(G,return_inverse=True)
 q0=np.bincount(iv,weights=(y-a)**2,minlength=len(u))
 q1=np.bincount(iv,weights=(y-b)**2,minlength=len(u))
 rng=np.random.default_rng(20261114)
 reps=rng.integers(0,len(u),size=(2000,len(u)))
 s0=q0[reps].sum(1);s1=q1[reps].sum(1)
 return list(map(float,np.quantile(100*(s0-s1)/s0,[.025,.5,.975])))
results={}
for name in ("ept_fraction","log1p_EPT_taxon_count","benthic_MMI"):
 y=Q[name].to_numpy(float)
 preds={key:cvpred(x,y) for key,x in M.items()}
 e={key:float(np.mean((y-v)**2)) for key,v in preds.items()}
 gain=lambda a,b:float(100*(e[a]-e[b])/e[a])
 results[name]={
  "n_native_UID":len(Q),"n_HUC8":Q.huc8.nunique(),
  "observed_y_mean":float(y.mean()),
  "HUC8_5fold_MSE":e,
  "N_P_added_to_geo_gain_pct":gain("geo","geo+N_P"),
  "N_P_added_after_PHab_gain_pct":gain("geo+PHab","geo+PHab+N_P"),
  "PHab_added_after_N_P_gain_pct":gain("geo+N_P","geo+PHab+N_P"),
  "HUC8_fixedOOF_2000_bootstrap95_N_P_after_PHab":boot(y,preds["geo+PHab"],preds["geo+PHab+N_P"]),
  "source_original_instream_NP_cluster_associations":slope(M["geo+PHab+N_P"],y)}
qa={
 "source_exact_precomputed_frame_UIDs":len(uids),"SOURCE_probability_wadeable_eligible":len(D),
 "source_method_MMI_BERWW_and_YES":int(source_mmi_ok.sum()),
 "source_numeric_same_UID_strict_cases":len(Q),"HUC8":Q.huc8.nunique(),
 "states":Q.state.nunique(),"distinct_physical_SITE_ID":Q.site_id.nunique(),
 "original_year_visits":{str(k):int(v) for k,v in Q.year.value_counts().items()},
 "source_original_same_UID_site_identity_verified":True,"source_date_match_QA":date_q,
 "original_source_unit_converted_NTL_PTL_ug_L_to_mg_L":True,
 "source_chem_alerts_excluded_not_zero":True,
 "source_original_2008_MMI_percent_converted_to_fraction":True,
 "source_ranges":ranges}
OUT={
 "schema":"gesui_NRSA0809_true_original_3_benthic_endpoints_vs_stream_NP_PHAB_v1",
 "genuine_EPA_original_2008_SOURCE_QA_manifest_SHA256":hashlib.sha256(pathlib.Path(os.environ["NRSA0809_FROZEN_SOURCE_QA"]).read_bytes()).hexdigest(),
 "source_2008_numerical_frame_was_selected_before_signs":True,
 "older_2008_MMI_and_EPT_are_not_independent_biological_replicates":True,
 "older_2008_BERWW_method_not_2013_18_BERW_assumed_identical":True,
 "unweighted_PROB_only_not_national_stream_miles_estimate":True,
 "not_2008_Chla_AFDM_and_not_WWTP_effluent":True,
 "not_ecosystem_service_or_emergent_insect_mass":True,
 "QA":qa,"observed_models":results}
p=R/"NRSA0809_ORIGINAL_EPT_TAXON_MMI_TRUE_NP_PHAB_RESULTS.json"
p.write_text(json.dumps(OUT,indent=2,sort_keys=True))
print("NRSA0809_REAL_EPT_TAXON_MMI_NP_PHAB_RESULT")
print(json.dumps({"QA":{k:v for k,v in qa.items() if k!="source_ranges"},
 "responses":{name:{"gainNP":v["N_P_added_after_PHab_gain_pct"],
 "gainPHab":v["PHab_added_after_N_P_gain_pct"],
 "coefficients":v["source_original_instream_NP_cluster_associations"],
 "HUC8_bootNP":v["HUC8_fixedOOF_2000_bootstrap95_N_P_after_PHab"]} for name,v in results.items()}},sort_keys=True))
