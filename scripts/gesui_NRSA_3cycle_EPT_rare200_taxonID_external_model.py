#!/usr/bin/env python3
"""EPA three real survey cycles, effort-standardized original taxon-ID richness.

Source-only QA and frozen protocol before fit: 2008 source BERWW raw total300
counts + original TAXA_ID->ORDER, 2013/18 original BERW count/ORDER. Hypergeom
expected richness at 200 individuals from the SAME original measured 2209
visits. Leave entire survey period + ALL same HUC8 out of source training.
Three responses always retained: EPT rare200 taxonID richness, EPT counted
fraction on EXACT SAME sites, and nonEPT rare200 original taxonID richness.
Only predictive observational macroecology, never sewage causality, species
richness, benthic biomass, insect emergence or realized ecosystem services.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib
import numpy as np,pandas as pd
from scipy.special import gammaln
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from scipy.stats import spearmanr

TEMP=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
R=TEMP/"gesui_threecycle_rarefied_taxa"
OLD=TEMP/"gesui_threecycle_ept_transfer"/"AUDIT_ONLY_NRSA_THREE_REAL_CYCLE_EPT_NP_PHAB.csv"
AUD=R/"AUDIT_ONLY_ORIGINAL_EPA_3CYCLE_TAXON_COUNTS.csv"
FIRST=pathlib.Path(os.environ["ORIGINAL_RARE200_SOURCE_ONLY_QA"])
orig=json.loads(FIRST.read_text())
if orig["schema"]!="gesui_EPA3cycle_SOURCEONLY_BERW_TAXAID_count200_rarefaction_population_v1" or not orig["precommitted_m200_gate_all_3cycle_pass"]:
 raise RuntimeError("Source first EPA 3cycle original m200 cohort/coverage failed")
full=pd.read_csv(OLD,dtype={"source_cycle":str,"source_UID":str,"huc8":str,"state":str,"eco9":str})
original_count=2648
if len(full)!=original_count or full.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Original previously source-locked 2648 native-UID frame missing")
counts=pd.read_csv(AUD,dtype={"cycle":str,"UID":str,"key":str,"order":str},keep_default_na=False)
if counts.duplicated(["cycle","UID","key"]).any():raise RuntimeError("Source original taxon-ID counts not unique per UID")
counts["count"]=pd.to_numeric(counts["count"],errors="raise")
if not np.isfinite(counts["count"]).all() or counts["count"].le(0).any() or (
 counts["count"]-counts["count"].round()).abs().gt(1e-8).any():
 raise RuntimeError("Original standardized TOTAL300 not positive integer")
counts["known"]=(counts.known.astype(str).str.strip().str.lower().isin(("true","1","yes")))
expected_n={"2008_09":476,"2013_14":957,"2018_19":776}
expected_state={"2008_09":46,"2013_14":47,"2018_19":47}
# No source sample-size or sign selection. "N" includes all original counted
# individual taxon records including genuinely unknown-order records.
m=200
def rarefy(T):
 N=int(round(T["count"].sum()))
 n=T["count"].to_numpy(dtype=float)
 if N<m:return None
 absent=np.zeros(len(T),dtype=float)
 mask=(N-n)>=m
 if mask.any():
  z=N-n[mask]
  logp=(gammaln(z+1)-gammaln(z-m+1)-
        gammaln(N+1)+gammaln(N-m+1))
  absent[mask]=np.exp(np.minimum(0.,logp))
 present=-np.expm1(np.log(np.maximum(absent,1e-300)))
 # Exact zero absence probability when n_i > N-m; p=1.
 present[~mask]=1.
 if not np.isfinite(present).all() or (present<0).any() or (present>1+1e-8).any():
  raise RuntimeError("Original 200-individual hypergeometric expected count invalid")
 ept=T.order.isin(("EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA")).to_numpy()
 known=T.known.to_numpy(bool)
 validkey=T.key.ne("").to_numpy()
 # unknown ORDER source stays unknown, never silently classified nonEPT.
 ept_contrib=present[ept&validkey]
 non_contrib=present[~ept&known&validkey]
 unk_contrib=present[~known&validkey]
 known_fraction=float(T.loc[known,"count"].sum()/N)
 if known_fraction<.95-1e-10:return None
 # Source original taxonID count includes nonblank taxa for "total",
 # unknown ORDER remains excluded from EPT and nonEPT estimates.
 return {
  "original_total_individual_count":N,
  "original_known_order_taxonID_individual_fraction":known_fraction,
  "rarefied_200_EPT_taxonID":float(ept_contrib.sum()),
  "rarefied_200_nonEPT_taxonID":float(non_contrib.sum()),
  "rarefied_200_total_taxonID":float(present[validkey].sum()),
  "rarefied_200_EPT_upper_if_unknown_order_all_EPT":float(ept_contrib.sum()+unk_contrib.sum()),
  "original_raw_EPT_taxonID_count":int((ept&validkey).sum()),
  "original_raw_nonEPT_taxonID_count":int((~ept&known&validkey).sum()),
  "raw_original_200_counted_fraction_concordance":float(T.loc[ept,"count"].sum()/N)}
vals=[]
for (cycle,uid),g in counts.groupby(["cycle","UID"],sort=False):
 out=rarefy(g)
 if out is not None:vals.append({"source_cycle":cycle,"source_UID":uid,**out})
rare=pd.DataFrame(vals)
if rare.duplicated(["source_cycle","source_UID"]).any():
 raise RuntimeError("Original true native UID counted unit duplicated in rarefaction")
df=full.merge(rare,on=["source_cycle","source_UID"],how="inner",validate="one_to_one")
if len(df)!=2209:raise RuntimeError("Original source-only m200 cohort changed, expected 2209, observed "+str(len(df)))
info={}
for era in expected_n:
 x=df[df.source_cycle.eq(era)]
 if len(x)!=expected_n[era] or x.state.nunique()!=expected_state[era]:
  raise RuntimeError("Original source count QA cohort differs by cycle: "+era)
 gold=orig["quality_gates_by_cycle"][era]["200individual_eligible_n_HUC8_state"]
 if [len(x),x.huc8.nunique(),x.state.nunique()]!=gold:
  raise RuntimeError("Published source-only n/HUC8/state bridge differs after hypergeom formula "+era)
 err=(x.raw_original_200_counted_fraction_concordance-x.ept_fraction).abs().max()
 if err>1e-9:raise RuntimeError("Source total300 relative EPT share changed for "+era)
 for metric in ("rarefied_200_EPT_taxonID","rarefied_200_nonEPT_taxonID","rarefied_200_total_taxonID"):
  if (x[metric]<0).any() or (~np.isfinite(x[metric])).any():
   raise RuntimeError("Expected original rare200 richness invalid "+era+"/"+metric)
 if (x.rarefied_200_EPT_taxonID>x.original_raw_EPT_taxonID_count+1e-8).any():
  raise RuntimeError("Raw source EPT taxon richness exceeded mathematically")
 info[era]={"n":len(x),"HUC8":x.huc8.nunique(),"states":x.state.nunique(),
  "source_n_counted_median":float(x.original_total_individual_count.median()),
  "source_expected_rare200_EPT_taxonID_mean":float(x.rarefied_200_EPT_taxonID.mean()),
  "source_expected_rare200_nonEPT_taxonID_mean":float(x.rarefied_200_nonEPT_taxonID.mean()),
  "source_expected_rare200_total_taxonID_mean":float(x.rarefied_200_total_taxonID.mean()),
  "EPT_fraction_vs_rare200_EPT_taxonID_Spearman":float(spearmanr(x.ept_fraction,x.rarefied_200_EPT_taxonID).statistic),
  "source_unknown_order_expected_richness_upper_gap_mean":float(
    (x.rarefied_200_EPT_upper_if_unknown_order_all_EPT-x.rarefied_200_EPT_taxonID).mean())}
df["log1p_rare200_EPT_taxonID"]=np.log1p(df.rarefied_200_EPT_taxonID)
df["log1p_rare200_nonEPT_taxonID"]=np.log1p(df.rarefied_200_nonEPT_taxonID)
df["log1p_rare200_total_taxonID"]=np.log1p(df.rarefied_200_total_taxonID)
df["log1p_rare200_EPT_upper"]=np.log1p(df.rarefied_200_EPT_upper_if_unknown_order_all_EPT)
df=df.sort_values(["source_cycle","source_UID"]).reset_index(drop=True)
df.to_csv(R/"AUDIT_ONLY_EPA_THREECYCLE_RARE200_SAME_UID.csv",index=False)
# New relative EPT fraction control on exactly the original rarefied 200
# cohort; enables effort-attribution comparison with true EPT richness.
q={
 "schema":"gesui_EPA3cycle_ORIGINAL_RARE200_EPT_nonEPT_taxonID_source_only_and_transport_population_v1",
 "original_source_only_prior_run":43568105806,
 "original_source_manifest_SHA256":hashlib.sha256(FIRST.read_bytes()).hexdigest(),
 "original_count_2648_UIDs_same_as_PR81":True,
 "source_m200_cohort_UIDs":len(df),
 "source_n_by_cycle":info,
 "sampling_effort_m200_exact_hypergeom":True,
 "original_taxon_identifier_not_species_and_not_biomass":True,
 "all_reconstructed_EPT_relative_source_shares_unchanged":True,
 "no_new_model_prediction_yet":True}
(R/"NRSA_THREECYCLE_RARE200_NUMERIC_SOURCE_RESULT_QA.json").write_text(json.dumps(q,indent=2,sort_keys=True))
print("NRSA_THREECYCLE_RARE200_NUMERIC_SOURCE_QA",json.dumps(info,sort_keys=True),flush=True)

CYCLES=["2008_09","2013_14","2018_19"]
CATS=["state","eco9"]
B=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NP=["log_tn","log_tp"]
H=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
FEATURES={"B":B,"B+NP":B+NP,"B+PHab":B+H,"B+PHab+NP":B+H+NP}
RESP={
 "EPT_rarefied_200_taxonID_log1p":"log1p_rare200_EPT_taxonID",
 "EPT_counted_fraction_same_sites":"ept_fraction",
 "nonEPT_rarefied_200_taxonID_log1p":"log1p_rare200_nonEPT_taxonID"}
def equal_weights(train):
 cnt=train.huc8.value_counts()
 n=len(train)
 w=(n/train.huc8.nunique())/train.huc8.map(cnt)
 return w.to_numpy(float)
def fit(train,test,col,features):
 prep=ColumnTransformer(transformers=[
  ("numbers",StandardScaler(),features),
  ("categories",OneHotEncoder(handle_unknown="ignore",sparse_output=False,dtype=float),CATS)
 ],sparse_threshold=0)
 A=prep.fit_transform(train);Z=prep.transform(test)
 weights=equal_weights(train)
 model=Ridge(alpha=10,solver="svd",fit_intercept=True).fit(
   A,train[col].to_numpy(float),sample_weight=weights)
 pred=model.predict(Z)
 if not np.isfinite(pred).all():raise RuntimeError("Heldout era EPT taxon source prediction not finite")
 return pred
def getgroup(target,y,p):
 u,inv=np.unique(target.huc8.to_numpy(str),return_inverse=True)
 err=np.bincount(inv,weights=(y-p)**2,minlength=len(u))/np.bincount(inv,minlength=len(u))
 return err
def score(target,y,preds):
 sq={k:getgroup(target,y,v) for k,v in preds.items()}
 mse={k:float(v.mean()) for k,v in sq.items()}
 gain=lambda b,a:float(100*(mse[b]-mse[a])/mse[b])
 def ci(b,a):
  rng=np.random.default_rng(20261114)
  idx=rng.integers(0,len(sq[b]),size=(2000,len(sq[b])))
  bas=sq[b][idx].mean(axis=1);alt=sq[a][idx].mean(axis=1)
  return list(map(float,np.quantile(100*(bas-alt)/bas,[.025,.5,.975])))
 return {"MSE_equal_HUC8":mse,
   "gain_PHab_added_AFTER_NP_pct":gain("B+NP","B+PHab+NP"),
   "gain_NP_added_AFTER_PHab_pct":gain("B+PHab","B+PHab+NP"),
   "HUC8_fixedOOF_bootstrap95_PHab_after_NP":ci("B+NP","B+PHab+NP"),
   "HUC8_fixedOOF_bootstrap95_NP_after_PHab":ci("B+PHab","B+PHab+NP"),
   "source_target_y_mean":float(np.mean(y)),
   "source_target_y_sd":float(np.std(y))}
external={};context={}
fullgroups={era:set(full.loc[full.source_cycle.eq(era),"huc8"]) for era in CYCLES}
for era in CYCLES:
 test=df[df.source_cycle.eq(era)].reset_index(drop=True)
 oth=df[~df.source_cycle.eq(era)]
 tr=oth.loc[~oth.huc8.isin(fullgroups[era])].reset_index(drop=True)
 forbidden=set(tr.huc8)&fullgroups[era]
 if forbidden:raise RuntimeError("Original full-period HUC8 leaked into rare200 training "+era)
 gate=(len(test)>=250 and test.huc8.nunique()>=150 and test.state.nunique()>=20
       and len(tr)>=300 and tr.huc8.nunique()>=200 and tr.state.nunique()>=15)
 qa={"heldout_original_source_uid_count_m200":len(test),
     "heldout_original_HUC8_m200":test.huc8.nunique(),
     "train_original_after_excluding_ALL_target_period_native_HUC8":len(tr),
     "train_original_HUC8":tr.huc8.nunique(),
     "train_states":tr.state.nunique(),
     "excluded_any_shared_original_HUC8_count":len(forbidden),
     "other_period_training_UID_attrition_due_to_original_fullHUC_exclusion":len(oth)-len(tr),
     "train_cycles":tr.source_cycle.value_counts().to_dict(),
     "pre_fixed_reduced_rarefied_cohort_gate":bool(gate)}
 if not gate:
  external[era]={"QA":qa,"gate":"SOURCE_SUPPORT_FAIL_NO_EFFECT_FIT"}
  continue
 results={}
 for label,col in RESP.items():
  pred={name:fit(tr,test,col,feat) for name,feat in FEATURES.items()}
  results[label]=score(test,test[col].to_numpy(float),pred)
 external[era]={"QA":qa,"gate":"TRUE_OUT_OF_SURVEY_AND_ALL_ORIGINAL_TARGET_HUC8",
                "responses":results}
 # Context only: within-era five whole-HUC8 disjoint folds on the same
 # original rarefied sample, not a substitute for unknown-cycle transfer.
 fold=np.array([int(hashlib.sha256(str(huc).encode()).hexdigest()[:8],16)%5 for huc in test.huc8])
 if min(np.bincount(fold,minlength=5))<25:
  context[era]={"gate":"INSUFFICIENT_EACH_CV_FOLD"}
 else:
  within={}
  for label,col in RESP.items():
   PP={k:np.full(len(test),np.nan) for k in FEATURES}
   for f in range(5):
    ktr=test.loc[fold!=f].reset_index(drop=True)
    kte=test.loc[fold==f].reset_index(drop=True)
    if set(ktr.huc8)&set(kte.huc8):raise RuntimeError("Within period HUC8 leak")
    for key,feat in FEATURES.items():
     PP[key][fold==f]=fit(ktr,kte,col,feat)
   if not all(np.isfinite(x).all() for x in PP.values()):
    raise RuntimeError("Within-cycle rare200 predictions missing")
   within[label]=score(test,test[col].to_numpy(float),PP)
  context[era]={"gate":"WHOLE_HUC8_5FOLD_SOURCE_ONLY",
                "responses":within}
# Gate universal prediction, not merely all point estimates positive.
rich_result=[v.get("responses",{}).get("EPT_rarefied_200_taxonID_log1p")
             for v in external.values()]
universality=(all(x is not None for x in rich_result) and
              all(x["HUC8_fixedOOF_bootstrap95_PHab_after_NP"][0]>0 for x in rich_result))
O={
 "schema":"gesui_NR​​SA_3cycle_true_original_RARE200_EPT_taxa_vs_same_UID_share_cross-era_v1",
 "source_original_sourceQA_run":43568105806,
 "m200_total_native_UID":len(df),
 "m200_site_cohort_sourceQA":info,
 "original_source_count_fraction_concordance":True,
 "true_EPT_taxonID_rarefaction_m200_not_species_or_drymass":True,
 "source_original_2008_BERWW_vs_later_BERW_protocol_not_identical":True,
 "source_original_2013_and_2018_periphyton_eligible_subcohorts":True,
 "no_WWTP_NP_discharge_or_realized_service":True,
 "original_precommitted_universal_PHab_3test_ci_positive":bool(universality),
 "original_precommitted_universal_claim_SUPPORTED":bool(universality),
 "external_whole_cycle_zero_original_HUC8":external,
 "within_source_cycle_HUC8_cv_context":context}
p=R/"NRSA_THREECYCLE_RARE200_EPT_SHARE_RICHNESS_EXTERNAL_TRANSPORT.json"
p.write_text(json.dumps(O,indent=2,sort_keys=True))
print("NRSA_THREECYCLE_RARE200_EPT_TAXONID_EXTERNAL_RESULT",json.dumps({
 "source_cohorts":{x:{"n":v["n"],"HUC8":v["HUC8"]} for x,v in info.items()},
 "original_bold_universal_rarefied_EPT_richness":bool(universality),
 "crossperiod":{
  era:{"gate":v["gate"],"QA":v["QA"],
    "responses":{
     name:{k:r[k] for k in
       ("gain_PHab_added_AFTER_NP_pct",
        "gain_NP_added_AFTER_PHab_pct",
        "HUC8_fixedOOF_bootstrap95_PHab_after_NP",
        "HUC8_fixedOOF_bootstrap95_NP_after_PHab")}
     for name,r in v.get("responses",{}).items()}}
  for era,v in external.items()},
 "internal":{
  era:{"gate":v["gate"],
       "responses":{
         name:{"PHab_gain":r["gain_PHab_added_AFTER_NP_pct"],
               "NP_gain":r["gain_NP_added_AFTER_PHab_pct"]}
         for name,r in v.get("responses",{}).items()}}
  for era,v in context.items()}},sort_keys=True))
