#!/usr/bin/env python3
"""Independent sklearn/statsmodels numerical QA of EPA NRSA nutrient,
periphyton and EPT observational models. Reuses source-locked *ephemeral*
analysis frames, not their regression design matrices, folds or model fits.
Original source row validation and ID join are in the prior source pipeline.
"""
import hashlib,json,os,pathlib
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
import statsmodels.api as sm

root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"
original=json.loads((root/"OBSERVED_NRSA_NUTRIENT_PERIPHYTON_EPT_TWO_ERAS_NO_WWTP_CAUSAL.json").read_text())
report={"schema":"NRSA_TWO_ERA_PERIPHYTON_EPT_SKLEARN_STATSMODELS_INDEPENDENT_NUMERIC_QA_v1",
        "original_source_sha":original["source_SHAs"],
        "not_independently_downloaded_field_data":True,
        "not_causal_or_WWTP_measured_exposure":True,"eras":{}}
base_num=["lat","lon","stream_order","year","season_sin","season_cos","log_cond"]
def fold(x):return int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5
def audit(era,strict):
 label="2018_19" if era=="2018_19" else "2023_24"
 df=pd.read_csv(root/("AUDIT_ONLY_"+label+"_PANEL.csv"),dtype={"UID":str,"huc8":str,"state":str,"eco9":str,"site_id":str},low_memory=False)
 ynames=["ept_fraction","log_chla"] if era=="2018_19" else ["log_afdm","log_chla"]
 strict_mask=df.chla_flag.fillna("").eq("")&df.tn_mgL_rawflag.fillna("").eq("")&df.tp_mgL_rawflag.fillna("").eq("")
 if era!="2018_19":strict_mask&=df.afdm_flag.fillna("").eq("")
 if not strict:
  def accepted(x):return ~x.fillna("").astype(str).str.contains("Q|NR|ND",regex=True)
  strict_mask=accepted(df.chla_flag)&accepted(df.tn_mgL_rawflag)&accepted(df.tp_mgL_rawflag)
  if era!="2018_19":strict_mask &= accepted(df.afdm_flag)
 df=df.loc[strict_mask].dropna(subset=base_num+["log_tn","log_tp"]+ynames).copy()
 df=df[(df.huc8!="")&(df.state!="")&(df.eco9!="")]
 df.reset_index(drop=True,inplace=True)
 X0=pd.concat([df[base_num].astype(float),
               pd.get_dummies(df[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 trainfold=np.array([fold(x) for x in df.huc8])
 Xdict={"baseline":X0,"TN":X0.assign(log_tn=df.log_tn.to_numpy(float)),
        "TP":X0.assign(log_tp=df.log_tp.to_numpy(float)),
        "TN+TP":X0.assign(log_tn=df.log_tn.to_numpy(float),log_tp=df.log_tp.to_numpy(float))}
 Xdict["TN+TP+chla"]=Xdict["TN+TP"].assign(log_chla=df.log_chla.to_numpy(float))
 cohort=original["years"][era]["strict" if strict else "quality_flag_sensitivity"]
 if len(df)!=cohort["source_quality_n"]:raise AssertionError("Cohort count changed")
 if df.huc8.nunique()!=cohort["n_HUC8"]:raise AssertionError("Distinct catchment count changed")
 out={"n":len(df),"n_HUC8":df.huc8.nunique(),"outcomes":{}}
 for name in ynames:
  yy=df[name].to_numpy(float)
  baseline=cohort["responses"][name]["huc8_disjoint_5fold"]
  models=list(baseline.keys())
  check={}
  for m in models:
   xx=Xdict[m].to_numpy(float)
   pred=np.empty(len(df))
   for k in range(5):
    sel=trainfold!=k;held=~sel
    scaler=StandardScaler().fit(xx[sel])
    ztrain=scaler.transform(xx[sel]);ztest=scaler.transform(xx[held])
    pred[held]=Ridge(alpha=10,fit_intercept=True,solver="svd").fit(ztrain,yy[sel]).predict(ztest)
   mse=float(np.mean((yy-pred)**2))
   expected=float(baseline[m]["MSE"])
   if abs(mse-expected)>max(1e-8,1e-6*abs(expected)):
    raise AssertionError(f"sklearn independent MSE differs {era}/{strict}/{name}/{m}: {mse} vs {expected}")
   check[m]={"sklearn_MSE":mse,"source_result_MSE":expected,
             "abs_discrepancy":abs(mse-expected)}
  # statsmodels parameter estimation independent of original pseudo-inverse
  rhs=Xdict["TN+TP+chla" if name=="ept_fraction" else "TN+TP"]
  f=sm.OLS(yy,sm.add_constant(rhs.to_numpy(float),has_constant="add"))
  robust=f.fit(cov_type="cluster",cov_kwds={"groups":df.huc8,"use_correction":True})
  slopes={}
  for nm in ["log_tn","log_tp"]+(["log_chla"] if name=="ept_fraction" else []):
   pos=rhs.columns.get_loc(nm)+1 # add_constant intercept at index 0
   b=float(robust.params[pos]*df[nm].std(ddof=0))
   primary=cohort["responses"][name]["partial_slopes"][nm]
   if abs(b-primary["beta_per_1SD"])>1e-7:
    raise AssertionError("Independent OLS per-SD slope disagree "+nm+" "+str((b,primary)))
   slopes[nm]={"statsmodels_beta_per_SD":b,"source_beta_per_SD":primary["beta_per_1SD"]}
  out["outcomes"][name]={"models":check,"slopes":slopes}
 return out
for era in ["2018_19","2023_24"]:
 report["eras"][era]={"strict":audit(era,True),"quality_flag_sensitivity":audit(era,False)}
report["conclusion"]="SUCCESS: all independent sklearn Ridge HUC8-grouped OOS MSE and statsmodels per-SD slopes match original within specified numerical tolerance; NOT independent source data or causal evidence"
(root/"NRSA_2ERA_PERIPHYTON_EPT_INDEPENDENT_NUMERIC_QA.json").write_text(json.dumps(report,indent=2,sort_keys=True))
print("NRSA_PERIPHYTON_NUMERIC_INDEPENDENT_QA_SUCCESS")
print(json.dumps({"cohorts":{era:{stage:{"n":d["n"],"huc8":d["n_HUC8"]}
                               for stage,d in report["eras"][era].items()} for era in report["eras"]},
                  "conclusion":report["conclusion"]},sort_keys=True))
