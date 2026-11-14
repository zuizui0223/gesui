#!/usr/bin/env python3
"""Original EPA NRSA 2008/13/18 same-EPT-site rarefied taxonID richness at
exactly 100 organisms, comparing:
(A) PREVIOUS m200 n2209 source UID cohort (depth-only robustness) versus
(B) source-only expanded m100 original UID cohort (count-selection change).
Three original responses on each same native UID: EPT taxon richness,
original counted EPT relative fraction and nonEPT taxon richness.

Other two survey cycles train after excluding ALL original target-period HUC8,
heldout era 100% untouched. Fixed PHab and ambient N/P source models,
train-only categories and scaling, HUC8-equal weights and 2000 bootstrap.
Exploratory/noncausal; no sewage source delivery or actual ecosystem service.
"""
from __future__ import annotations
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from scipy.special import gammaln
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=TEMP/"gesui_threecycle_rarefied_taxa"
src100=json.loads(pathlib.Path(os.environ["ORIGINAL_M100_SOURCE_ONLY_QA"]).read_text())
old200=json.loads(pathlib.Path(os.environ["ORIGINAL_RARE200_MODEL_RESULT"]).read_text())
if src100["schema"]!="gesui_EPA_3cycle_NRSA_m100_vs_m200_sourceONLY_two_population_frames_v1":
 raise RuntimeError("EPA original source-only m100 cohort not frozen")
if old200["schema"]!="gesui_NR​​SA_3cycle_true_original_RARE200_EPT_taxa_vs_same_UID_share_cross-era_v1":
 raise RuntimeError("EPA original rare200 survey cycle result missing source")
count=pd.read_csv(R/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv",
 dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
frozen=pd.read_csv(TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv",
 dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
if len(frozen)!=2648 or frozen.duplicated(["source_cycle","source_UID"]).any():raise RuntimeError("Original source n2648 changed")
count["known"]=count.known.astype(str).str.lower().isin(("true","1","yes"))
count["count"]=pd.to_numeric(count["count"],errors="raise")
if count.duplicated(["cycle","UID","key"]).any() or (count["count"]<=0).any():
 raise RuntimeError("Original EPA original counting source duplicate or nonpositive taxon entry")
m=100
rows=[]
for (cycle,uid),d in count.groupby(["cycle","UID"],sort=False):
 n=d["count"].to_numpy(float)
 N=int(round(n.sum()))
 if N<m:continue
 if abs(n.sum()-N)>1e-8 or (n-n.round()).max()>1e-8:
  raise RuntimeError("Original counted TOTAL300 source N non-integer")
 known=d.known.to_numpy(bool)
 knownshare=float(n[known].sum()/N)
 if knownshare<.95-1e-10:continue
 choose=(N-n)>=m
 p=np.ones(len(n))
 if choose.any():
  q=N-n[choose]
  logp=gammaln(q+1)-gammaln(q-m+1)-gammaln(N+1)+gammaln(N-m+1)
  p[choose]=-np.expm1(np.minimum(logp,0))
 if not np.isfinite(p).all() or np.any(p<0) or np.any(p>1+1e-8):
  raise RuntimeError("Actual EPA counted source hypergeom expected taxa invalid")
 ept=d.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA")).to_numpy()
 valid=d.key.ne("").to_numpy()
 eptrich=float(p[ept&valid].sum())
 nonrich=float(p[~ept&known&valid].sum())
 upper=float(p[(ept|~known)&valid].sum())
 rows.append({"source_cycle":cycle,"source_UID":uid,
              "source_count_N":N,"source_known_order_fraction":knownshare,
              "source_rare100_ept_taxonID":eptrich,
              "source_rare100_nonEPT_taxonID":nonrich,
              "source_rare100_ept_upper_bound":upper,
              "source_fraction_from_original_count":float(n[ept].sum()/N)})
ra=pd.DataFrame(rows)
df=frozen.merge(ra,on=["source_cycle","source_UID"],how="inner",validate="one_to_one")
if df.duplicated(["source_cycle","source_UID"]).any():raise RuntimeError("Original source ID join expanded duplicated")
expected={"2008_09":(523,401,47,476),"2013_14":(1035,661,47,957),"2018_19":(845,549,47,776)}
oldcount=pd.read_csv(R/"NRSA_THREECYCLE_RARE200_NUMERIC_SOURCE_RESULT_QA.json") if False else None
# m200 originals are frozen by actual count source quality, not a newly
# chosen positive association or exposure cutoff.
m200same=df[df.source_count_N.ge(200)].copy()
if len(m200same)!=2209:raise RuntimeError("Original n2209 m200 source quality changed")
sourceQA={}
for era,(n,h,s,m200n) in expected.items():
 v=df[df.source_cycle.eq(era)]
 w=m200same[m200same.source_cycle.eq(era)]
 prev=src100["source_frames"][era]
 if [len(v),v.huc8.nunique(),v.state.nunique(),len(w)]!=[n,h,s,m200n]:
  raise RuntimeError("Source-only pre-effect n/HUC/state or m200 n mismatch "+era)
 if prev["m100_expanded_source_n"]!=n or prev["m100_expanded_HUC8"]!=h:
  raise RuntimeError("Source-only native original 100 count selected frame changed "+era)
 # Counted original EPT fraction unaffected by rarefaction.
 err=(v.source_fraction_from_original_count-v.ept_fraction).abs().max()
 if err>1e-9:raise RuntimeError("Source old original counted EPT share differed "+era)
 sourceQA[era]={
  "m100_n":n,"m100_HUC8":h,"m200_same_n":m200n,
  "new_n_after_adding_original100_to199_counted":n-m200n,
  "original_per_count_source_EPT_fraction_maxdiff":float(err),
  "source_original_rare100_EPT_taxon_ID_mean":float(v.source_rare100_ept_taxonID.mean()),
  "source_original_rare100_EPT_taxon_ID_mean_on_m200_cohort":float(w.source_rare100_ept_taxonID.mean()),
  "source_original_rare100_nonEPT_taxonID_mean":float(v.source_rare100_nonEPT_taxonID.mean()),
  "same_UID_original_rare200_source_EPT_taxonID_mean":float(old200["m200_site_cohort_sourceQA"][era]["source_expected_rare200_EPT_taxonID_mean"])}
 print("NRSA_SOURCE_RARE100_EPT_INDEPENDENT_COHORT",era,json.dumps(sourceQA[era],sort_keys=True),flush=True)
if len(df)!=2403:raise RuntimeError("Source m100 expanded total is not source-only 2403")
df["log1p_rare100_ept"]=np.log1p(df.source_rare100_ept_taxonID)
df["log1p_rare100_nonEPT"]=np.log1p(df.source_rare100_nonEPT_taxonID)
df["log1p_rare100_ept_upper"]=np.log1p(df.source_rare100_ept_upper_bound)
df=df.sort_values(["source_cycle","source_UID"]).reset_index(drop=True)
# Source only full 2648 original target cycle HUC8 used for exclusion:
# the m100 threshold never modifies the train-disjoint HUC8 decision.
YEARS=["2008_09","2013_14","2018_19"]
fullgroups={era:set(frozen.loc[frozen.source_cycle.eq(era),"huc8"]) for era in YEARS}
base=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
phab=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
NUT=["log_tn","log_tp"]
models={"B":base,"B+NP":base+NUT,"B+PHab":base+phab,"B+PHab+NP":base+phab+NUT}
resp={"EPT_rare100_taxonID_log1p":"log1p_rare100_ept",
      "EPT_counted_fraction_same_site":"ept_fraction",
      "nonEPT_rare100_taxonID_log1p":"log1p_rare100_nonEPT"}
def fit(train,test,response,features):
 prep=ColumnTransformer(
  [("numeric",StandardScaler(),features),
   ("category",OneHotEncoder(handle_unknown="ignore",sparse_output=False,dtype=float),["state","eco9"])],
  sparse_threshold=0,remainder="drop")
 X=prep.fit_transform(train);T=prep.transform(test)
 cnt=train.huc8.value_counts()
 weights=(len(train)/len(cnt))/train.huc8.map(cnt)
 out=Ridge(alpha=10,solver="svd").fit(X,train[response].to_numpy(float),
                                    sample_weight=weights.to_numpy(float)).predict(T)
 if not np.isfinite(out).all():raise RuntimeError("Original source no-leak out-of-era predictions not finite")
 return out
def calc(target,y,pred):
 u,g=np.unique(target.huc8.to_numpy(str),return_inverse=True)
 n=np.bincount(g,minlength=len(u))
 err={key:np.bincount(g,weights=(y-vals)**2,minlength=len(u))/n for key,vals in pred.items()}
 mse={key:float(val.mean()) for key,val in err.items()}
 def gain(a,b):return float(100*(mse[a]-mse[b])/mse[a])
 rng=np.random.default_rng(20261114)
 draws=rng.integers(0,len(u),size=(2000,len(u)))
 def bs(a,b):
  M=err[a][draws].mean(axis=1);N=err[b][draws].mean(axis=1)
  return list(map(float,np.quantile(100*(M-N)/M,[.025,.5,.975])))
 return {
  "source_original_UID_n":len(target),"source_original_HUC8":len(u),
  "primary_equal_HUC8_MSE":mse,
  "PHab_after_NP_gain_pct":gain("B+NP","B+PHab+NP"),
  "NP_after_PHab_gain_pct":gain("B+PHab","B+PHab+NP"),
  "HUC8_fixedOOF_2000boot95_PHab":bs("B+NP","B+PHab+NP"),
  "HUC8_fixedOOF_2000boot95_NP":bs("B+PHab","B+PHab+NP")}
results={}
for frame in ["same_2209_m200_site_UIDs_at_m100_depth","expanded_2403_m100_site_UIDs"]:
 cohort=m200same if frame.startswith("same") else df
 block={}
 for era in YEARS:
  test=cohort[cohort.source_cycle.eq(era)].reset_index(drop=True)
  train=cohort[(~cohort.source_cycle.eq(era))&
      (~cohort.huc8.isin(fullgroups[era]))].reset_index(drop=True)
  if set(train.huc8)&set(test.huc8):raise RuntimeError("Original m100 external HUC8 leak "+frame+"/"+era)
  eligibility=(len(test)>=250 and test.huc8.nunique()>=150 and test.state.nunique()>=20 and
    len(train)>=300 and train.huc8.nunique()>=200 and train.state.nunique()>=15)
  qa={"target_original_native_uid_n":len(test),"target_HUC8":test.huc8.nunique(),
      "training_source_uid_after_full_original_target_HUC8_excluded":len(train),
      "training_HUC8":train.huc8.nunique(),"training_states":train.state.nunique(),
      "HUC8_overlap_train_test":len(set(train.huc8)&set(test.huc8)),
      "source_predeclared_support_gate":bool(eligibility)}
  if not eligibility:
   block[era]={"QA":qa,"gate":"SOURCE_SUPPORT_FAIL_NO_FIT"};continue
  responses={}
  for label,endpoint in resp.items():
   predictions={key:fit(train,test,endpoint,columns) for key,columns in models.items()}
   responses[label]=calc(test,test[endpoint].to_numpy(float),predictions)
  block[era]={"QA":qa,"gate":"TRUE_EXTERNAL_ERA_WITH_ZERO_FULL_ORIGINAL_TARGET_HUC8",
              "responses":responses}
 results[frame]=block
# Compare m100 cohort fixed sites directly to earlier m200 main signed output,
# which had exactly the same target/training original SOURCE UID populations.
depth_deltas={}
for era in YEARS:
 cur=results["same_2209_m200_site_UIDs_at_m100_depth"][era]
 old=old200["external_whole_cycle_zero_original_HUC8"][era]
 if cur["gate"]!="TRUE_EXTERNAL_ERA_WITH_ZERO_FULL_ORIGINAL_TARGET_HUC8" or old["gate"]!="TRUE_OUT_OF_SURVEY_AND_ALL_ORIGINAL_TARGET_HUC8":
  depth_deltas[era]={"gate":"MISSING_ORIGINAL_SOURCE_COMPARABLE"}
  continue
 if cur["QA"]["target_original_native_uid_n"]!=old["QA"]["heldout_original_source_uid_count_m200"]:
  raise RuntimeError("Original m200–m100 same-site source target changed "+era)
 if cur["QA"]["training_source_uid_after_full_original_target_HUC8_excluded"]!=old["QA"]["train_original_after_excluding_ALL_target_period_native_HUC8"]:
  raise RuntimeError("Original m200–m100 same-site source TRAIN changed "+era)
 before=old["responses"]["EPT_rarefied_200_taxonID_log1p"]
 after=cur["responses"]["EPT_rare100_taxonID_log1p"]
 depth_deltas[era]={
  "original_same_UID_200_expected_EPT_taxonID_PHab_gain_pct":before["gain_PHab_added_AFTER_NP_pct"],
  "new_same_UID_100_expected_EPT_taxonID_PHab_gain_pct":after["PHab_after_NP_gain_pct"],
  "delta_PHab_gain_100minus200_percentage_points":after["PHab_after_NP_gain_pct"]-before["gain_PHab_added_AFTER_NP_pct"],
  "original_same_UID_200_expected_EPT_taxonID_NP_gain_pct":before["gain_NP_added_AFTER_PHab_pct"],
  "new_same_UID_100_expected_EPT_taxonID_NP_gain_pct":after["NP_after_PHab_gain_pct"],
  "delta_NP_gain_100minus200_percentage_points":after["NP_after_PHab_gain_pct"]-before["gain_NP_added_AFTER_PHab_pct"],
  "same_source_train_and_test_native_UIDs":True}
checks={}
for frame,block in results.items():
 checks[frame]=all(
  r["gate"]=="TRUE_EXTERNAL_ERA_WITH_ZERO_FULL_ORIGINAL_TARGET_HUC8" and
  r["responses"]["EPT_rare100_taxonID_log1p"]["HUC8_fixedOOF_2000boot95_PHab"][0]>0
 for r in block.values())
O={
 "schema":"gesui_NRSA3cycle_RARE100_same200vs_expanded_source_huc8_external_validated_v1",
 "old_prior_original_source_PR82_m200_n2209_run":43568477035,
 "original_m100_source_only_manifest_SHA256":hashlib.sha256(pathlib.Path(os.environ["ORIGINAL_M100_SOURCE_ONLY_QA"]).read_bytes()).hexdigest(),
 "m100_source_2403_U​​IDs":len(df),"original_m200_source_2209_U​​IDs":len(m200same),
 "original_counted_site_cohort_no_sign_tuning":True,
 "source_3cycle_each_100rarefaction_depth_exact_hypergeom":True,
 "original_3cycle_original_frozen_2648_huc8_excluded_from_train":True,
 "m100_same_sites_vs_expanded_sites_do_NOT_compare_as_pure_depth_effect":True,
 "no_synthetic_natural_species_richness_or_insect_biomass":True,
 "no_actual_municipal_WWTP_delivered_NP":True,
 "sourceQA":sourceQA,
 "source_rare100_HUC8_external_models":results,
 "strict_same_2209_UIDs_m100_vs_original_m200_depth_only":depth_deltas,
 "universal_PHab_richness_external_95_positive_per_cohort":checks}
path=R/"NRSA_THREECYCLE_RARE100_SAMESITES_VS_EXPANDED_HUC8_EXTERNAL.json"
path.write_text(json.dumps(O,indent=2,sort_keys=True))
print("ORIGINAL_EPA_RARE100_COMPARE200_EXTERNAL_SUCCESS")
print(json.dumps({
 "cohorts":{era:{"n100":s["m100_n"],"n200":s["m200_same_n"],"newUIDs":s["new_n_after_adding_original100_to199_counted"]} for era,s in sourceQA.items()},
 "same_site_depth_delta":depth_deltas,
 "universal_PHab_richness_criterion":checks,
 "m100_results":{
  frame:{era:{"gate":row["gate"],"n":row["QA"]["target_original_native_uid_n"],
   "train_n":row["QA"]["training_source_uid_after_full_original_target_HUC8_excluded"],
   "EPT_richness":{k:v for k,v in row.get("responses",{}).get("EPT_rare100_taxonID_log1p",{}).items()
     if k in ("PHab_after_NP_gain_pct","NP_after_PHab_gain_pct","HUC8_fixedOOF_2000boot95_PHab","HUC8_fixedOOF_2000boot95_NP")},
   "EPT_fraction_PHab":row.get("responses",{}).get("EPT_counted_fraction_same_site",{}).get("PHab_after_NP_gain_pct")}
   for era,row in block.items()} for frame,block in results.items()}},sort_keys=True))
