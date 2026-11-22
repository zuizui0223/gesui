#!/usr/bin/env python3
"""Exploratory original EPA NRSA2324 ambient river TP x canopy x benthos.
Prespecified posthoc source protocol in docs/PROTOCOL_20261122_*.md.
NOT wastewater-derived TP, actual insect emergence, or causal evidence.
"""
import argparse,collections,hashlib,json,pathlib
import numpy as np,pandas as pd,statsmodels.api as sm
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,OneHotEncoder
from sklearn.linear_model import Ridge

SHA={"siteinfo":"f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013",
"waterchem":"b34d3f9c4b4bdc763cb5e53a856bbab329f1646d248f223f78e9404bc4560d71",
"physicalhabitat":"8a6094bb87f8a57cf45a35a1272b137e6a9327cefd158f4bee88b025eaf9b2f4",
"benthiccount":"5290d679f7eb444e97c8d45096beb4531b450186cef190580f41f4b9430b0805",
"benthicmmi":"5ef607e8e4726ed8c3d24b9d15fe351e7300784859654b812fa595c09a4009dc"}
NUM=["lat","lon","order","year","season_sin","season_cos","log_cond","log_tn","log_width","log_tp","canopy"]
CAT=["state","eco9"]
def sha(p):
 h=hashlib.sha256()
 with open(p,"rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def one(df,name):
 if df.UID.duplicated().any():raise ValueError("Duplicate native UID "+name)
 return df.set_index("UID")
def num(x):return pd.to_numeric(x,errors="coerce")
def log(x):return np.log10(num(x).where(num(x)>0))
def chem_value(df,field,typ):
 raw=num(df[field+"_RESULT"]);u=df[field+"_UNITS"].astype(str).str.strip()
 out=np.full(len(raw),np.nan)
 m=u.isin(["mg/L","mg P/L"] if typ=="P" else ["mg/L","mg N/L"])
 out[m.to_numpy()]=raw[m].to_numpy()
 if typ=="P":
  m=u.isin(["ug/L","ug P/L"])
  out[m.to_numpy()]=raw[m].to_numpy()/1000
 return pd.Series(np.where(out>0,out,np.nan),index=df.index)
def date(df):return pd.to_datetime(df.DATE_COL,format="mixed",errors="coerce").dt.normalize()
def panel(t):
 # Original native ID / same site+visit+date gate before numeric outcomes.
 common=set.intersection(*[set(v.UID) for v in t.values()])
 if len(common)!=2139:raise ValueError("Original five-table common-UID count changed")
 keys=sorted(common);ss=one(t["siteinfo"],"siteinfo").loc[keys]
 others={k:one(t[k],"single row table").loc[keys] for k in ("waterchem","physicalhabitat","benthicmmi")}
 cnt=t["benthiccount"];cntmeta=cnt.drop_duplicates("UID").set_index("UID").loc[keys]
 check={}
 for name,f in list(others.items())+[("benthiccount",cntmeta)]:
  violations=int((ss.SITE_ID.to_numpy()!=f.SITE_ID.to_numpy()).sum()+
             (ss.VISIT_NO.to_numpy()!=f.VISIT_NO.to_numpy()).sum()+
             (date(ss).to_numpy()!=date(f).to_numpy()).sum())
  check[name]=violations
  if violations:raise ValueError("Site/visit/date mismatch "+name+": "+str(violations))
 wc=others["waterchem"];ph=others["physicalhabitat"];bm=others["benthicmmi"]
 co=cnt[cnt.UID.isin(common)].copy()
 co["n"]=num(co.TOTAL300)
 co=co[co.n.gt(0)&np.isfinite(co.n)].copy()
 co["ept"]=co.n*co.ORDER.astype(str).str.lower().str.strip().isin(
   ["ephemeroptera","plecoptera","trichoptera"]).astype(float)
 cc=co.groupby("UID")[["n","ept"]].sum()
 cc["EPT_fraction"]=cc.ept/cc.n
 df=pd.DataFrame(index=ss.index)
 for x,field in (("site","SITE_ID"),("huc8","HUC8"),("state","PSTL_CODE"),("eco9","AG_ECO9")):
  df[x]=ss[field].str.strip()
 df["method"]=bm.SAMPLE_TYPE.str.strip();df["count_method"]=cntmeta.SAMPLE_TYPE.str.strip()
 df["protocol"]=ph.PROTOCOL.str.strip()
 df["lat"]=num(ss.LAT_DD83);df["lon"]=num(ss.LON_DD83);df["order"]=num(ss.STRAH_ORD)
 day=date(ss);df["year"]=day.dt.year.astype(float)
 df["season_sin"]=np.sin(2*np.pi*day.dt.dayofyear/365.25)
 df["season_cos"]=np.cos(2*np.pi*day.dt.dayofyear/365.25)
 df["log_cond"]=log(wc.COND_RESULT)
 df["log_tn"]=np.log10(chem_value(wc,"NTL","N"))
 df["log_tp"]=np.log10(chem_value(wc,"PTL","P"))
 df["log_width"]=log(ph.XWIDTH)
 df["tp_ok"]=wc.PTL_NARS_FLAG.str.strip().eq("")
 df["tn_ok"]=wc.NTL_NARS_FLAG.str.strip().eq("")
 df["mid"]=num(ph.XCDENMID);df["bank"]=num(ph.XCDENBK)
 df["MMI_BENT"]=num(bm.MMI_BENT)
 df=df.join(cc[["n","EPT_fraction"]])
 return df,check,len(common)
def cohort(df,method,y):
 if method=="BERW":
  m=df.method.eq("BERW")&df.count_method.eq("BERW")&df.protocol.eq("WADEABLE")
  canopy="mid"
 else:
  m=df.method.eq("BETB")&df.count_method.eq("BETB")&df.protocol.eq("BOATABLE")
  canopy="bank"
 d=df[m & df.tp_ok & df.tn_ok].copy();d["canopy"]=d[canopy]
 d=d[d.canopy.between(0,100)].replace([np.inf,-np.inf],np.nan).dropna(subset=NUM+[y])
 d=d[(d.huc8!="")&(d.site!="")&(d.state!="")&(d.eco9!="")]
 if y=="EPT_fraction":d=d[(d.n>0)&d.EPT_fraction.between(0,1)]
 return d.sort_index()
def test(d,y):
 folds=np.asarray([int(hashlib.sha256(str(h).encode()).hexdigest()[:8],16)%5 for h in d.huc8])
 nh=[int(d.loc[folds==f,"huc8"].nunique()) for f in range(5)]
 if min(nh)<3:raise ValueError("Insufficient HUC8 in heldout fold")
 p0=np.full(len(d),np.nan);p1=np.full(len(d),np.nan)
 for f in range(5):
  tr=folds!=f;te=~tr
  enc=ColumnTransformer([("nums",StandardScaler(),NUM),
    ("categories",OneHotEncoder(handle_unknown="ignore",sparse_output=False),CAT)],
    sparse_threshold=0)
  xt=enc.fit_transform(d.loc[tr,NUM+CAT]);xs=enc.transform(d.loc[te,NUM+CAT])
  ti,ci=NUM.index("log_tp"),NUM.index("canopy")
  ia=(xt[:,ti]*xt[:,ci]).reshape(-1,1)
  ib=(xs[:,ti]*xs[:,ci]).reshape(-1,1)
  base=Ridge(alpha=10,solver="svd").fit(xt,y[tr])
  full=Ridge(alpha=10,solver="svd").fit(np.c_[xt,ia],y[tr])
  p0[te]=base.predict(xs);p1[te]=full.predict(np.c_[xs,ib])
 if not (np.isfinite(p0).all() and np.isfinite(p1).all()):
  raise ValueError("Not all heldout predictions finite")
 e0=(y-p0)**2;e1=(y-p1)**2
 g0=pd.Series(e0,index=d.huc8.to_numpy()).groupby(level=0).mean().sort_index()
 g1=pd.Series(e1,index=d.huc8.to_numpy()).groupby(level=0).mean().reindex(g0.index)
 A=g0.to_numpy();B=g1.to_numpy()
 ix=np.random.default_rng(20261122).integers(0,len(A),size=(2000,len(A)))
 boot=100*(A[ix].mean(axis=1)-B[ix].mean(axis=1))/A[ix].mean(axis=1)
 # Descriptive HUC8-cluster OLS, response and continuous terms standardized.
 z=d[NUM].astype(float)
 z=(z-z.mean())/z.std(ddof=0).replace(0,1)
 X=pd.concat([z.reset_index(drop=True),
    pd.get_dummies(d[CAT],drop_first=True,dtype=float).reset_index(drop=True)],axis=1)
 X["TP_x_canopy"]=z.log_tp.to_numpy()*z.canopy.to_numpy()
 X=sm.add_constant(X,has_constant="add")
 sy=float(np.std(y,ddof=0))
 if sy<=0:raise ValueError("Zero response variance")
 res=sm.OLS((y-np.mean(y))/sy,X.to_numpy(dtype=float)).fit(
    cov_type="cluster",cov_kwds={"groups":d.huc8.to_numpy(),"use_correction":True})
 i=list(X.columns).index("TP_x_canopy");ci=res.conf_int()[i]
 return {"oof_HUC8_fold_groups":nh,
  "HUC8_equal_ref_MSE":float(np.mean(A)),"HUC8_equal_full_MSE":float(np.mean(B)),
  "HUC8_equal_interaction_gain_pct":float(100*(np.mean(A)-np.mean(B))/np.mean(A)),
  "HUC8_equal_fixedOOF_bootstrap95_pct":[float(q) for q in np.quantile(boot,[.025,.975])],
  "visit_weighted_interaction_gain_pct":float(100*(np.mean(e0)-np.mean(e1))/np.mean(e0)),
  "interaction_beta_ySD_per_logTP_SD_times_canopy_SD":float(res.params[i]),
  "interaction_HUC8_cluster_descriptive95":[float(ci[0]),float(ci[1])],
  "interaction_descriptive_p":float(res.pvalues[i]),"posthoc_not_causal":True}
def main():
 a=argparse.ArgumentParser()
 for k in SHA:a.add_argument("--"+k,type=pathlib.Path,required=True)
 a.add_argument("--out",type=pathlib.Path,required=True);z=a.parse_args()
 t={}
 for k,pin in SHA.items():
  p=getattr(z,k)
  if sha(p)!=pin:raise ValueError("Official EPA source SHA mismatch "+k)
  t[k]=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
 df,qa,source_n=panel(t)
 results={}
 for method in ("BERW","BETB"):
  results[method]={}
  for response in ("MMI_BENT","EPT_fraction"):
   d=cohort(df,method,response)
   data={"status":"SOURCE_STOP","visits":len(d),"SITE_ID":d.site.nunique(),
         "HUC8":d.huc8.nunique(),"states":d.state.nunique(),
         "canopy_field":"XCDENMID" if method=="BERW" else "XCDENBK"}
   if len(d)>=250 and d.huc8.nunique()>=100:
    try:
     data.update(test(d,d[response].to_numpy(dtype=float)))
     data["status"]="EXPLORATORY_NUMERIC_COMPLETE"
    except Exception as e:
     data.update({"status":"ANALYSIS_STOP","error":type(e).__name__+": "+str(e)[:280]})
   results[method][response]=data
   print("OBSERVED_NRSA_BENTHOS_TP_CANOPY",method,response,json.dumps(data,sort_keys=True),flush=True)
 result={"schema":"gesui_NRSA2324_original_ambient_TP_canopy_benthos_posthoc_v1",
   "original_EPA_SHA256":SHA,"original_source_joint_UID":source_n,
   "original_site_visit_date_QA_mismatches":qa,"methods_not_pooled":True,
   "pre_numeric_source_min_visits":250,"pre_numeric_source_min_HUC8":100,
   "posthoc":True,"not_wastewater_TP":True,"no_causal_inference":True,
   "responses":results}
 z.out.parent.mkdir(parents=True,exist_ok=True)
 z.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("NRSA_AMBIENT_TP_CANOPY_SUMMARY",json.dumps(
  {m:{r:{"status":v["status"],"n":v["visits"],"HUC8":v["HUC8"],
  "CVgain":v.get("HUC8_equal_interaction_gain_pct"),"beta":v.get("interaction_beta_ySD_per_logTP_SD_times_canopy_SD")}
  for r,v in x.items()} for m,x in results.items()},sort_keys=True),flush=True)
if __name__=="__main__":main()
