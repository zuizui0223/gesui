#!/usr/bin/env python3
"""EPA fourth cycle 2023–24 genuinely heldout survey-era + HUC8 EPT test.

This analysis is posthoc to prior 2008/2013/2018 and 2023 studies. EPA original
same-visit field nutrient, benthic counted EPT and habitat, NOT sewage-derived N/P,
actual emergent biomass, cross-ecosystem flux or randomized restoration.
"""
from __future__ import annotations
import argparse, collections, hashlib, json, pathlib
import numpy as np, pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.linear_model import Ridge
import gesui_nrsa2324_fourthcycle_source_only_20261122 as qa2023

HIST_ROOT_NAME="gesui_threecycle_ept_transfer"
OLD_EXPECT={"2008_09":(586,437,48),"2013_14":(1116,706,48),"2018_19":(946,603,47)}
CATS=["state","eco9"]
B=["lat","lon","stream_order","season_sin","season_cos","log_cond"]
NP=["log_tn","log_tp"]
PHAB=["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
ALL_COLS=["source_cycle","source_UID","site_id","huc8","state","eco9","year"]+B+NP+PHAB+["ept_fraction"]
MODEL={
 "B":B,
 "B+NP":B+NP,
 "B+PHab":B+PHAB,
 "B+NP+PHab":B+NP+PHAB,
 "B+NP+width":B+NP+["log_width"],
 "B+NP+PHab_except_width":B+NP+[c for c in PHAB if c!="log_width"],
 "B+NP+bed":B+NP+["XEMBED"],
 "B+NP+PHab_except_bed":B+NP+[c for c in PHAB if c!="XEMBED"],
 "B+NP+canopy":B+NP+["XCDENMID","XPCAN","XFC_OHV"],
 "B+NP+human":B+NP+["W1_HALL"]
}
COMPARE=[
 ("NP_after_PHab","B+PHab","B+NP+PHab"),
 ("PHab_after_NP","B+NP","B+NP+PHab"),
 ("width_alone_after_NP","B+NP","B+NP+width"),
 ("width_after_NP_other5","B+NP+PHab_except_width","B+NP+PHab"),
 ("bed_alone_after_NP","B+NP","B+NP+bed"),
 ("bed_after_NP_other5","B+NP+PHab_except_bed","B+NP+PHab"),
 ("canopy_after_NP","B+NP","B+NP+canopy"),
 ("human_after_NP","B+NP","B+NP+human")
]
def normalize_huc(raw):
 x=raw.astype(str).str.strip()
 pref=x.str.fullmatch(r"[Hh]\d{8}")
 z=x.where(~pref,x.str.slice(1))
 valid=z.str.fullmatch(r"\d{7,8}")
 if not valid.all():raise ValueError("Original 2023 HUC8 has source-unrecognized native label")
 canonical=z.str.zfill(8)
 if canonical.nunique()!=x.nunique():raise ValueError("HUC8 source id collapsed")
 return canonical
def n(x):return pd.to_numeric(x,errors="coerce")
def covariates(args,qa):
 site=pd.read_csv(args.sites,dtype=str,keep_default_na=False,low_memory=False)
 wc=pd.read_csv(args.water,dtype=str,keep_default_na=False,low_memory=False)
 ph=pd.read_csv(args.habitat,dtype=str,keep_default_na=False,low_memory=False)
 chla=pd.read_csv(args.chla,dtype=str,keep_default_na=False,low_memory=False)
 for label,frame in (("sites",site),("water",wc),("habitat",ph),("chla",chla)):
  if frame.UID.duplicated().any():raise ValueError("Duplicate original UID "+label)
 site=site.set_index("UID");wc=wc.set_index("UID");ph=ph.set_index("UID")
 chl=chla.set_index("UID")
 # Exactly SOURCE ONLY signature; the numeric panel must not change eligibility.
 source_ids=set(site.index)&set(wc.index)&set(ph.index)
 if len(source_ids)<2139:raise ValueError("2023 native source intersection shrank")
 d=site.loc[sorted(source_ids)].copy()
 w=wc.reindex(d.index);h=ph.reindex(d.index)
 dates=pd.to_datetime(d.DATE_COL,format="mixed",errors="coerce")
 tmp=pd.DataFrame(index=d.index)
 tmp["site_id"]=d.SITE_ID.str.strip()
 tmp["huc8"]=normalize_huc(d.HUC8)
 tmp["state"]=d.PSTL_CODE.str.strip()
 tmp["eco9"]=d.AG_ECO9.str.strip()
 tmp["source_cycle"]="2023_24"
 tmp["source_UID"]=d.index
 tmp["lat"]=n(d.LAT_DD83)
 tmp["lon"]=n(d.LON_DD83)
 tmp["stream_order"]=n(d.STRAH_ORD)
 tmp["year"]=n(d.YEAR)
 tmp["season_sin"]=np.sin(2*np.pi*dates.dt.dayofyear/365.25)
 tmp["season_cos"]=np.cos(2*np.pi*dates.dt.dayofyear/365.25)
 tmp["log_cond"]=np.log10(n(w.COND_RESULT).where(n(w.COND_RESULT)>0))
 tmp["log_tn"]=np.log10(qa2023.correct_units(w,"NTL","N"))
 tmp["log_tp"]=np.log10(qa2023.correct_units(w,"PTL","P"))
 for x in PHAB:
  if x=="log_width":
   q=n(h.XWIDTH);tmp[x]=np.log10(q.where(q>0))
  else:tmp[x]=n(h[x])
 tmp["site_ok"]=w.SITE_ID.eq(d.SITE_ID)&h.SITE_ID.eq(d.SITE_ID)
 tmp["visit_ok"]=w.VISIT_NO.eq(d.VISIT_NO)&h.VISIT_NO.eq(d.VISIT_NO)
 tmp["date_ok"]=(pd.to_datetime(w.DATE_COL,format="mixed",errors="coerce").to_numpy().astype("datetime64[D]")==dates.to_numpy().astype("datetime64[D]"))&(
   pd.to_datetime(h.DATE_COL,format="mixed",errors="coerce").to_numpy().astype("datetime64[D]")==dates.to_numpy().astype("datetime64[D]"))
 tmp["chem_flag_ok"]=w.NTL_NARS_FLAG.astype(str).str.strip().eq("")&w.PTL_NARS_FLAG.astype(str).str.strip().eq("")
 tmp["phab_protocol_ok"]=h.PROTOCOL.astype(str).str.strip().eq("WADEABLE")
 tmp["chla_flag_ok"]=tmp.index.isin(chl.index[chl.NARS_FLAG.str.strip().eq("")])
 if not tmp[["site_ok","visit_ok","date_ok"]].all().all():raise ValueError("Original 2023 wadeable native source ID/date mismatch")
 keep=tmp.phab_protocol_ok&tmp.chem_flag_ok
 good=tmp.loc[keep, B+NP+PHAB+["year"]].replace([np.inf,-np.inf],np.nan).notna().all(axis=1)
 candidate=tmp.loc[good.index[good]].copy()
 candidate=candidate[(candidate.site_id!="")&(candidate.state!="")&(candidate.eco9!="")]
 if not candidate.year.isin([2023,2024]).all():raise ValueError("2023 source year mismatch")
 # Response source membership ONLY; original counted taxon numbers not opened
 # until the independent source-only signed n/HUC gate has passed.
 return candidate
def make_target(args,source_qa):
 b=pd.read_csv(args.benthos,
   usecols=["UID","SAMPLE_TYPE","ORDER","TOTAL300","SITE_ID","VISIT_NO","DATE_COL"],
   dtype=str,keep_default_na=False,low_memory=False)
 b["count"]=n(b.TOTAL300)
 b=b[b.SAMPLE_TYPE.eq("BERW") & b["count"].gt(0)].copy()
 b["EPT_count"]=b["count"]*b.ORDER.str.lower().str.strip().isin(
   ("ephemeroptera","plecoptera","trichoptera")).astype(float)
 g=b.groupby("UID",sort=True).agg(count_total=("count","sum"),EPT_total=("EPT_count","sum"))
 g["ept_fraction"]=g.EPT_total/g.count_total
 d=covariates(args,source_qa)
 target=d.join(g[["ept_fraction","count_total"]],how="left")
 target=target.dropna(subset=["ept_fraction"])
 if not target.ept_fraction.between(0,1).all():raise ValueError("Original 2023 EPT counted share outside [0,1]")
 target=target.dropna(subset=B+NP+PHAB+["year"]).copy()
 target=target[(target.count_total>0)&target.site_ok&target.visit_ok&target.date_ok]
 target=target.sort_index()
 return target
def canonical_target(frame,previous_audit,kind):
 if kind=="full":t=frame
 elif kind=="periphyton":t=frame[frame.chla_flag_ok].copy()
 else:raise ValueError("Unknown source cohort")
 meta=previous_audit["2023_original_all_source_valid" if kind=="full" else "2023_original_periphyton_eligible_source_valid"]
 if (len(t)!=meta["UID_visits"] or t.huc8.nunique()!=meta["HUC8"] or
     t.site_id.nunique()!=meta["physical_SITE_ID"] or t.state.nunique()!=meta["states"]):
  raise RuntimeError("Original source-only 2023 pre-outcome membership count changed "+kind)
 if hashlib.sha256("\n".join(sorted(t.huc8.unique())).encode()).hexdigest()!=meta["original_huc8_sha256"]:
  raise RuntimeError("Original target catchment source-frozen SHA changed")
 out=t[["source_cycle","source_UID","site_id","huc8","state","eco9","year"]+B+NP+PHAB+["ept_fraction"]].copy()
 if out.duplicated(["source_cycle","source_UID"]).any():raise ValueError("Duplicate UID")
 if not np.isfinite(out[B+NP+PHAB+["ept_fraction"]].to_numpy(float)).all():raise ValueError("Invalid target numeric")
 return out.reset_index(drop=True)
def history(path):
 d=pd.read_csv(path,dtype={"source_UID":str,"huc8":str,"state":str,"eco9":str,"source_cycle":str})
 if len(d)!=2648:raise RuntimeError("Original three source-era full sample changed")
 for key,(nn,ng,ns) in OLD_EXPECT.items():
  v=d[d.source_cycle.eq(key)]
  if len(v)!=nn or v.huc8.nunique()!=ng or v.state.nunique()!=ns:
   raise RuntimeError("Original "+key+" EPA source frozen sample changed")
 if d.duplicated(["source_cycle","source_UID"]).any():raise RuntimeError("Old source duplicate UID")
 if not d.huc8.str.fullmatch(r"\d{8}").all():raise RuntimeError("Old original HUC8 normalization changed")
 return d
def equal_weights(train):
 counts=train.huc8.value_counts()
 w=(len(train)/len(counts))/train.huc8.map(counts)
 if not np.isfinite(w).all():raise ValueError("Source HUC8 weighting failure")
 return w.to_numpy(float)
def predict(train,test,columns,ycol="ept_fraction",weight_by_huc=True):
 pre=ColumnTransformer([
   ("numeric",StandardScaler(),columns),
   ("categorical",OneHotEncoder(handle_unknown="ignore",sparse_output=False),CATS)],
   sparse_threshold=0,remainder="drop")
 xa=pre.fit_transform(train);xb=pre.transform(test)
 if not np.isfinite(xa).all() or not np.isfinite(xb).all():raise RuntimeError("Source model nan")
 y=train[ycol].to_numpy(float)
 w=equal_weights(train) if weight_by_huc else None
 model=Ridge(alpha=10,solver="svd").fit(xa,y,sample_weight=w)
 out=model.predict(xb)
 if not np.isfinite(out).all():raise RuntimeError("Nonfinite source target predictions")
 return out
def evaluate(test,predictions):
 y=test.ept_fraction.to_numpy(float)
 scores={}
 errors={}
 for name,p in predictions.items():
  e=(y-p)**2
  g=pd.Series(e,index=test.huc8.to_numpy()).groupby(level=0).mean().sort_index()
  scores[name]={"MSE_equal_HUC8":float(g.mean()),"MSE_equal_UID":float(e.mean())}
  errors[name]=g
 gains={}
 for tag,a,b in COMPARE:
  x=errors[a].to_numpy();z=errors[b].to_numpy()
  seed=int(hashlib.sha256(tag.encode()).hexdigest()[:8],16)
  ix=np.random.default_rng(seed).integers(0,len(x),size=(2000,len(x)))
  ci=100*(x[ix].mean(axis=1)-z[ix].mean(axis=1))/x[ix].mean(axis=1)
  gains[tag]={
    "HUC8_equal_gain_pct":float(100*(x.mean()-z.mean())/x.mean()),
    "UID_equal_gain_pct":float(100*(scores[a]["MSE_equal_UID"]-scores[b]["MSE_equal_UID"])/scores[a]["MSE_equal_UID"]),
    "fixed_target_HUC8_bootstrap95":list(map(float,np.quantile(ci,[.025,.975])))}
 mean_huc=test.groupby("huc8").ept_fraction.mean()
 target_var=float(pd.Series((y-mean_huc.mean())**2,
   index=test.huc8.to_numpy()).groupby(level=0).mean().mean())
 return {
   "source_test_UID":len(test),"source_test_HUC8":test.huc8.nunique(),
   "EPT_mean_equal_HUC8":float(mean_huc.mean()),
   "target_HUC8_ept_variance":target_var,
   "model_errors":scores,
   "relative_to_target_huc_variance_R2":{k:1-v["MSE_equal_HUC8"]/target_var for k,v in scores.items()},
   "comparisons":gains
 }
def all_predictions(train,target,weight_by_huc=True):
 return {name:predict(train,target,columns,weight_by_huc=weight_by_huc) for name,columns in MODEL.items()}
def main():
 a=argparse.ArgumentParser()
 for k in qa2023.SHA:a.add_argument("--"+k,type=pathlib.Path,required=True)
 a.add_argument("--source-qa",type=pathlib.Path,required=True)
 a.add_argument("--old-panel",type=pathlib.Path,required=True)
 a.add_argument("--out",type=pathlib.Path,required=True)
 x=a.parse_args()
 srcq=json.loads(x.source_qa.read_text())
 if srcq.get("original_EPA_sha256")!=qa2023.SHA:raise ValueError("Original 2023 source-only SHA mismatch")
 if not srcq["2023_original_all_source_valid"]["source_go"]:
  raise RuntimeError("Target source gate failed: no numeric outcome testing")
 for name,orig_sha in qa2023.SHA.items():
  if qa2023.sha(getattr(x,name))!=orig_sha:raise ValueError("Original 2023 EPA SHA mismatch "+name)
 old=history(x.old_panel)
 target=make_target(x,srcq)
 if not srcq["2023_original_all_source_valid"]["source_go"]:
  raise RuntimeError("No fourth cycle source gate")
 # All source memberships were frozen in a different Actions run before biological
 # predictions. Do not change cohorts after seeing EPT response signs.
 tsets={k:canonical_target(target,srcq,k) for k in ("full","periphyton")}
 out={}
 for case,T in tsets.items():
  source_huc=set(T.huc8)
  other=old
  tr=other.loc[~other.huc8.isin(source_huc)].copy().reset_index(drop=True)
  ntrain,htrain=len(tr),tr.huc8.nunique()
  diag={
   "target_UID":len(T),"target_HUC8":T.huc8.nunique(),
   "historical_three_source_n":len(other),"historical_3source_union_HUC8":other.huc8.nunique(),
   "train_after_removing_every_shared_target_HUC8_UID":ntrain,
   "train_HUC8":htrain,"target_train_HUC8_overlap":len(set(tr.huc8)&source_huc),
   "train_era_counts":tr.source_cycle.value_counts().to_dict(),
   "target_eco9_unseen_in_training":sorted(set(T.eco9)-set(tr.eco9)),
   "target_states_unseen_in_training":sorted(set(T.state)-set(tr.state))}
  if diag["target_train_HUC8_overlap"]:raise RuntimeError("Catchment target leakage")
  if not (ntrain>=500 and htrain>=300):
   out[case]={"status":"TRAIN_SOURCE_STOP","source":diag};continue
  p=all_predictions(tr,T)
  metrics=evaluate(T,p)
  novel=set(T.huc8)-set(old.huc8)
  sub=T[T.huc8.isin(novel)].reset_index(drop=True)
  if len(sub)>=150 and sub.huc8.nunique()>=100:
   mask=T.huc8.isin(novel).to_numpy()
   p2={k:arr[mask] for k,arr in p.items()}
   m2=evaluate(sub,p2)
   newresult={"status":"NOVEL_TARGET_HUC8_COMPLETE","results":m2}
  else:newresult={"status":"NOVEL_TARGET_HUC8_SOURCE_STOP",
    "UID":len(sub),"HUC8":sub.huc8.nunique()}
  # Sensitivity to historical HUC weighting, same source and same 2023 target:
  pu=all_predictions(tr,T,weight_by_huc=False)
  unweighted=evaluate(T,pu)
  out[case]={"status":"EXTERNAL_FOURTH_CYCLE_COMPLETE",
    "source":diag,"target":metrics,
    "novel_target_only_HUC8":newresult,
    "source_train_visit_equal_instead_of_HUC8_equal":{
      "comparisons":{k:unweighted["comparisons"][k] for k in ("NP_after_PHab","PHab_after_NP","width_alone_after_NP","bed_alone_after_NP")}}}
  print("NRSA_FOURTH_CYCLE_ONE_RESULT",case,json.dumps(out[case],sort_keys=True),flush=True)
 report={"schema":"gesui_original_EPA_4thcycle_2023_24_EPT_NP_PHab_outofera_HUC8_holdout_v1",
   "2023_source_original_SHA":qa2023.SHA,
   "original_3cycles_frozen":OLD_EXPECT,
   "original_target_source_only_n":srcq["2023_original_all_source_valid"],
   "original_target_periphyton_n":srcq["2023_original_periphyton_eligible_source_valid"],
   "three_older_eras_used_to_predict_original_unseen_2023_24":True,
   "train_target_zero_shared_huc8_enforced":True,
   "posthoc_after_prior_original_2023_and_threecycle_signs_known":True,
   "not_historical_sewage_WWTP_dose":True,
   "not_independent_preregistered_experiment":True,
   "results":out}
 x.out.parent.mkdir(parents=True,exist_ok=True)
 x.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("NRSA_FOURTH_CYCLE_SUMMARY",json.dumps({case:{
   "status":r["status"],"test_UID":r["source"]["target_UID"],
   "test_HUC8":r["source"]["target_HUC8"],"train_UID":r["source"]["train_after_removing_every_shared_target_HUC8_UID"],
   "train_HUC8":r["source"]["train_HUC8"],
   "gains":{k:{q:float(z[q]) if isinstance(z[q],(int,float)) else z[q] for q in ("HUC8_equal_gain_pct","fixed_target_HUC8_bootstrap95")}
    for k,z in r["target"]["comparisons"].items()}}
    for case,r in out.items() if r.get("status")=="EXTERNAL_FOURTH_CYCLE_COMPLETE"},sort_keys=True),flush=True)
if __name__=="__main__":main()
