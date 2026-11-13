#!/usr/bin/env python3
"""New genuine 2018–19 NRSA same-visit measured freshwater nutrient, periphyton
Chla, benthic EPT counted fraction, and benthic taxon-ID richness with PHab.

Sample selection SOURCE blinded through frozen source QA/PR58/PR73. Analyses
are exploratory observational co-state, not sewer effluent, mediation, adult
emergence, primary productivity or ecosystem services.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib,time,urllib.request
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm

ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui1819_triad_numeric"
ROOT.mkdir(exist_ok=True,parents=True)
PREVIOUS=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"/"AUDIT_ONLY_2018_19_PANEL.csv"
BIOLOGICAL_SRC="https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv"
BIO_SHA="4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7"
PHAB_SRC="https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv"
PHAB_SHA="69b3b56962c04492423bd7eaf4d074c6ae8d253924d683ba12eb25fd857ce063"
SOURCE_MANIFEST=pathlib.Path(os.environ["FROZEN_TRIAD_SOURCE"])
source=json.loads(SOURCE_MANIFEST.read_text())
if source["schema"]!="gesui_NRSA1819_EPT_Chla_PHab_true_native_source_gate_v1":
 raise RuntimeError("Original SOURCE-only EPT+Chla+physical study frame not frozen")
valid_uids=set(source["source_UIDs"])
if len(valid_uids)!=1221:raise RuntimeError("Expected 1221 fully source-qualified original visits")
df=pd.read_csv(PREVIOUS,dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str})
if len(df)!=1221 or set(df.UID)!=valid_uids or df.UID.duplicated().any():
 raise RuntimeError("PR73 frozen native BERW + true Chla cohort changed")

def fetch(url,expected,name):
 p=ROOT/name
 for i in range(4):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"gesui-NRSA-2018-2019-real-triad/1.0"})
   with urllib.request.urlopen(req,timeout=100) as r,p.open("wb") as f:
    while data:=r.read(1048576):f.write(data)
   if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
    raise RuntimeError("Original EPA source changed "+name)
   return p
  except Exception:
   p.unlink(missing_ok=True)
   if i==3:raise
   time.sleep(2**i)
path=fetch(PHAB_SRC,PHAB_SHA,"phab.csv")
hab=pd.read_csv(path,dtype=str,keep_default_na=False,low_memory=False)
path.unlink();hab.columns=hab.columns.str.upper()
if hab.UID.duplicated().any():raise RuntimeError("PHab has duplicate native original UID")
cols=["XPCAN","XGVEG","XFC_OHV","W1_HALL","XWIDTH","XCDENMID","XEMBED"]
if any(c not in hab for c in cols):raise RuntimeError("Original physical habitat schema changed")
hab=hab[hab.UID.isin(valid_uids)].copy().set_index("UID")
for c in cols:hab[c]=pd.to_numeric(hab[c],errors="coerce")
hab["log_width"]=np.log10(hab.XWIDTH.where(hab.XWIDTH>0))
if len(hab)!=1221:raise RuntimeError("Original PHab UID join is incomplete")
p=fetch(BIOLOGICAL_SRC,BIO_SHA,"benthic.csv")
B=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
p.unlink()
B=B[B.UID.isin(valid_uids)&B.SAMPLE_TYPE.eq("BERW")].copy()
B["count"]=pd.to_numeric(B.TOTAL300,errors="coerce")
B=B[B["count"].gt(0)&np.isfinite(B["count"])].copy()
if not set(valid_uids).issubset(set(B.UID)):raise RuntimeError("Original benthic count data incomplete")
if any(B.TAXA_ID.eq("")):raise RuntimeError("Positive native benthic count TAXA_ID missing")
duplicate_taxon=int(B.duplicated(["UID","TAXA_ID"]).sum())
if duplicate_taxon>len(B)*0.005:raise RuntimeError("Original native taxon ID duplicates exceed frozen source QC threshold")
B["is_ept"]=B.ORDER.str.strip().str.lower().isin(["ephemeroptera","plecoptera","trichoptera"])
G=B.groupby("UID")
taxa=G.TAXA_ID.nunique()
tot=G["count"].sum()
eptcount=(B["count"]*B.is_ept.astype(int)).groupby(B.UID).sum()
bio=pd.DataFrame({"native_taxonID_richness":taxa,"total_benthic_subsample_count":tot,
 "raw_EPT_fraction":(eptcount/tot)})
x=df.set_index("UID").join(hab[cols+["log_width"]]).join(bio)
if not np.allclose(x.ept_fraction.to_numpy(float),x.raw_EPT_fraction.to_numpy(float),equal_nan=True,atol=1e-10,rtol=0):
 raise RuntimeError("Independent raw EPA 2018 EPT source count disagrees with original PR73 EPT fraction")
x["log_taxon_richness"]=np.log1p(x.native_taxonID_richness)
x["log_total_subsample_count"]=np.log(x.total_benthic_subsample_count)
flag=x.chla_flag.fillna("").eq("")&x.tn_mgL_rawflag.fillna("").eq("")&x.tp_mgL_rawflag.fillna("").eq("")
core=["lat","lon","year","stream_order","season_sin","season_cos","log_cond","log_tn","log_tp",
      "log_chla","ept_fraction","log_taxon_richness","log_total_subsample_count"]
strict=x[flag].dropna(subset=core).copy()
strict=strict[(strict.huc8.notna())&(strict.huc8!="")&(strict.state.notna())&(strict.state!="")&(strict.eco9.notna())&(strict.eco9!="")]
if len(strict)!=947 or strict.huc8.nunique()!=603 or strict.state.nunique()!=47:
 raise RuntimeError(f"Original 2018/19 strict biota+chemistry+areal Chla cohort not equal n947/603/47: {len(strict)}")
# SOURCE-registered (not chosen by sign), mutually nested PHab tiers.
PHAB=[
 ("broad",["XPCAN","XGVEG","XFC_OHV","W1_HALL"],600,250,25),
 ("dense",["XPCAN","XGVEG","XFC_OHV","W1_HALL","log_width","XCDENMID","XEMBED"],430,170,20)
]
def hashed_fold(h):
 return int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5
def designs(d,phabcols,with_count):
 bcols=["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]
 if with_count:bcols=bcols+["log_total_subsample_count"]
 base=pd.concat([d[bcols].astype(float),
                 pd.get_dummies(d[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 nut=d[["log_tn","log_tp"]].astype(float)
 hab=d[phabcols].astype(float)
 return {"base":base,
         "N":pd.concat([base,nut],axis=1),
         "H":pd.concat([base,hab],axis=1),
         "NH":pd.concat([base,hab,nut],axis=1)}
def cvpred(x,Y,fold):
 a=x.to_numpy(dtype=float)
 pred=np.full(len(Y),np.nan)
 for k in range(5):
  tr=fold!=k;te=~tr
  if te.sum()<50 or tr.sum()<200:raise RuntimeError("Insufficient geographic CV holdout")
  sd=StandardScaler().fit(a[tr])
  z=sd.transform(a)
  pred[te]=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(z[tr],Y[tr]).predict(z[te])
 if not np.isfinite(pred).all():raise RuntimeError("Non-finite geographic validation")
 return pred
def coefficient(x,Y,groups,names):
 XX=x.copy()
 for c in names:
  XX[c]=(XX[c]-XX[c].mean())/XX[c].std(ddof=0)
 val=sm.OLS(Y,sm.add_constant(XX.to_numpy(dtype=float),has_constant="add")).fit(
  cov_type="cluster",cov_kwds={"groups":groups,"use_correction":True})
 return {name:{"beta_1SD_original_predictor":float(val.params[XX.columns.get_loc(name)+1]),
               "HUC8_cluster_descriptive95":list(map(float,val.conf_int()[XX.columns.get_loc(name)+1]))} for name in names}
def predict_gain(y,p0,p1):
 e0=np.mean((y-p0)**2);e1=np.mean((y-p1)**2)
 return float(100*(e0-e1)/e0)
def bootstrap(y,pred0,pred1,groups):
 ids,inv=np.unique(groups,return_inverse=True)
 e0=np.bincount(inv,weights=(y-pred0)**2,minlength=len(ids))
 e1=np.bincount(inv,weights=(y-pred1)**2,minlength=len(ids))
 rng=np.random.default_rng(20261114)
 idx=rng.integers(0,len(ids),size=(1600,len(ids)))
 a=e0[idx].sum(axis=1);b=e1[idx].sum(axis=1)
 good=a>0
 return list(map(float,np.quantile(100*(a[good]-b[good])/a[good],[.025,.5,.975])))
results={}
for name,fields,minn,minh,mins in PHAB:
 d=strict.dropna(subset=fields).sort_index().copy()
 qa={"n_original_biota_chla_quality_before_phab":947,
     "n_joined_complete_PHAB":len(d),
     "n_HUC8":d.huc8.nunique(),"n_states":d.state.nunique(),
     "n_distinct_physical_SITE_ID":d.site_id.nunique(),
     "original_habitat_vars":fields,
     "source_zero_and_missing_distinguished":True}
 if len(d)<minn or d.huc8.nunique()<minh or d.state.nunique()<mins:
  results[name]={"source_QA":qa,"status":"SOURCE_COMPLETE_NUMERIC_GATE_FAIL_NOT_AN_ECOLOGICAL_NULL"}
  continue
 fold=np.array([hashed_fold(z) for z in d.huc8])
 qa["status"]="SOURCE_COMPLETE_GATES_PASS"
 qa["XGVEG_greater_than_one_n"]=int((d.XGVEG>1).sum())
 predout={}
 groups=d.huc8.to_numpy()
 for target,extra,richen in [("ept_fraction","log_chla",False),
                             ("log_chla","ept_fraction",False),
                             ("log_taxon_richness","log_chla",True)]:
  Y=d[target].to_numpy(float)
  mats=designs(d,fields,richen)
  mats["NH+other_biota"]=mats["NH"].assign(**{extra:d[extra].to_numpy(float)})
  predictions={lab:cvpred(xx,Y,fold) for lab,xx in mats.items()}
  mse={lab:float(np.mean((Y-v)**2)) for lab,v in predictions.items()}
  improvements={
   "N_added_to_base_pct":predict_gain(Y,predictions["base"],predictions["N"]),
   "N_added_after_PHAB_pct":predict_gain(Y,predictions["H"],predictions["NH"]),
   "PHAB_added_after_N_pct":predict_gain(Y,predictions["N"],predictions["NH"]),
   "co_measured_other_biota_added_after_N_PHAB_pct":predict_gain(Y,predictions["NH"],predictions["NH+other_biota"])
  }
  slopes=coefficient(mats["NH+other_biota"],Y,groups,["log_tn","log_tp",extra])
  predout[target]={"n":len(d),"mean_response":float(np.mean(Y)),
                   "mse_HUC8_group5fold":mse,
                   "incremental_MSE_improvement_pct":improvements,
                   "clustered_descriptive_association":slopes,
                   "HUC8_fixedheldout_bootstrap95":{
                     "N_after_PHAB_pct":bootstrap(Y,predictions["H"],predictions["NH"],groups),
                     "other_biota_after_N_PHAB_pct":bootstrap(Y,predictions["NH"],predictions["NH+other_biota"],groups)},
                   "other_biota_source_predictor":extra}
 # These are raw observational ranks, not independent sample organisms.
 correlations={
  "EPT_fraction_vs_logChla_spearman":float(spearmanr(d.ept_fraction,d.log_chla).statistic),
  "taxonID_richness_vs_logChla_spearman":float(spearmanr(d.native_taxonID_richness,d.log_chla).statistic),
  "EPT_fraction_vs_taxonID_richness_spearman":float(spearmanr(d.ept_fraction,d.native_taxonID_richness).statistic)
 }
 results[name]={"source_QA":qa,"status":"OBSERVATIONAL_SAME_UID_FITS_COMPLETE",
                "observed_within_cohort_correlations":correlations,
                "responses":predout}
out={
 "schema":"gesui_NRSA1819_TRUE_SAME_UID_EPT_Chla_TAXONRICHNESS_PHAB_NP_OBSERVATIONAL_MACRO_v1",
 "original_EPA_CSV_SHA256":{"benthos":BIO_SHA,"phab":PHAB_SHA},
 "original_fixed_source_manifest_SHA256":hashlib.sha256(SOURCE_MANIFEST.read_bytes()).hexdigest(),
 "original_same_biological_EPT_quality_947_603_HUC8_47states":True,
 "new_source_numeric_comparisons_known_beforehand":True,
 "different_years_2023AFDM_not_falsely_joined":True,
 "benthic_taxonID_richness_not_species_richness":True,
 "Chla_ug_per_cm2_is_biofilm_pigment_stock_not_primary_production":True,
 "benthic_EPT_counted_fraction_not_adult_emergence_biomass":True,
 "no_measured_wastewater_origin_or_pollination_service":True,
 "posthoc_derived_associations_not_independent_study":True,
 "results":results}
outpath=ROOT/"NRSA1819_SAME_UID_BIOTIC_TRIAD_PHAB_NP_REAL_RESULTS.json"
outpath.write_text(json.dumps(out,sort_keys=True,indent=2))
print("NRSA1819_REAL_TRIAD_RESULTS_START")
print(json.dumps({"QA":{k:v["source_QA"] for k,v in results.items()},
 "results":{tier:{resp:{"n":v["n"],
   "MSE_gains":v["incremental_MSE_improvement_pct"],
   "partial":v["clustered_descriptive_association"],
   "HUC8_bootstrap":v["HUC8_fixedheldout_bootstrap95"]}
   for resp,v in d.get("responses",{}).items()}
   for tier,d in results.items()},
 "raw_rank":{k:v.get("observed_within_cohort_correlations") for k,v in results.items()}},sort_keys=True))
print("NRSA1819_REAL_TRIAD_RESULTS_END")
