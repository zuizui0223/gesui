#!/usr/bin/env python3
"""SOURCE-PINNED exploratory EPA NRSA matched visit measured nutrient vs.
area-unit periphyton AFDM/Chla and benthic EPT fraction, across two periods.

2018/19: original 1221 UID benthic source-frozen pre-nutrient-manifest from
PR #58; 2018/19 Chla original UID -> canonical UID via blind site/date bridge.
2023/24: 2129 original pre-effect UID intersection manifest from PR #73.
NO WWTP load attribution, primary productivity, actual biomass flux or service.
"""
import collections,csv,hashlib,json,math,os,pathlib,time,urllib.request
import numpy as np,pandas as pd
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_biomass_numeric"
ROOT.mkdir(exist_ok=True,parents=True)
EXPECTED={
 "18_sites":("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_SiteInfo.csv","1b5b2119c881290191edeb47f1d80a7711f10440e9851574cfa86fd57be770e3"),
 "18_water":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_water_chemistry_chla_-_data.csv","cb9a3349ae3708834eab533e71423d74c1778bf76aba0737f0ae1bcc8f50baa1"),
 "18_benthos":("https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_benthic_macroinvertebrate_count_-_data.csv","4e9b3a0798b25f67c7a96d4800b94f86544bbd623abaf119358e067b158324a7"),
 "18_chla":("https://www.epa.gov/system/files/other-files/2023-01/NRSA_1819_PeriChla.csv","47a3663e9cc5e99ec54e581bda09c831604c0f267153f56c56ec6558221acf0a"),
 "23_sites":("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_siteinfo.csv","f7d14cf1574d12828b3aa97bd7b08b95f22aec2115753d5b75514f8eed43d013"),
 "23_water":("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_waterchemistry.csv","b34d3f9c4b4bdc763cb5e53a856bbab329f1646d248f223f78e9404bc4560d71"),
 "23_biomass":("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonbiomass.csv","a0d1608224f65b5b576bceea2d952311894f48d9f367bd68e956a4afb441e3a5"),
 "23_chla":("https://www.epa.gov/system/files/other-files/2026-06/nrsa2324_periphytonchlorophyll-a.csv","aec44b58806481f5ee53ce66d61c46b4b9fad344d2da2cc8312fceb3e1fb192e"),
}
MANIFEST_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
manifest_path=pathlib.Path(os.environ["SOURCE_FROZEN_18_MANIFEST"])
crosswalk_path=pathlib.Path(os.environ["SOURCE_FROZEN_BIOMASS_CROSSWALK"])
if hashlib.sha256(manifest_path.read_bytes()).hexdigest()!=MANIFEST_SHA:
 raise RuntimeError("Existing original 2018 original benthic UID freeze changed")
manifest=json.loads(manifest_path.read_text())
if len(manifest["primary_berw_UID_subset"])!=1221:raise RuntimeError("Unexpected 2018 original eligible UID count")
frozen18={r["uid"] for r in manifest["primary_berw_UID_subset"]}
cross=json.loads(crosswalk_path.read_text())
for n, k, era in [("18_sites","sites","2018_19"),("18_water","water","2018_19"),
                   ("18_chla","chla","2018_19"),("23_sites","sites","2023_24"),
                   ("23_water","water","2023_24"),("23_biomass","biomass","2023_24"),
                   ("23_chla","chla","2023_24")]:
 if cross[era]["source_sha256"][k]!=EXPECTED[n][1]:
  raise RuntimeError("Original before-outcome source freeze checksum differs: "+n)
def load(name):
 url,original_sha=EXPECTED[name]
 req=urllib.request.Request(url,headers={"User-Agent":"gesui-nutrient-biomass-repro/20261113"})
 p=ROOT/(name+".csv")
 for i in range(4):
  try:
   with urllib.request.urlopen(req,timeout=90) as resp,p.open("wb") as out:
    while b:=resp.read(1048576):out.write(b)
   if hashlib.sha256(p.read_bytes()).hexdigest()!=original_sha:
    raise RuntimeError("Original EPA SOURCE SHA mismatch: "+name)
   df=pd.read_csv(p,dtype=str,keep_default_na=False,low_memory=False)
   p.unlink()
   return df
  except Exception as e:
   p.unlink(missing_ok=True)
   if i==3:raise
   time.sleep(2**i)
D={name:load(name) for name in EXPECTED}
def uniq(df,label):
 if df.UID.duplicated().any():raise RuntimeError(label+" repeated original UID, stop")
 return df.set_index("UID")
def numeric(v):
 return pd.to_numeric(v,errors="coerce")
def unit_transform(v,unit,kind):
 u=str(unit).strip()
 if kind=="N":
  if u not in ("mg N/L","mg/L"):raise RuntimeError("N nutrient unit unrecognized: "+repr(u))
  return float(v)
 if kind=="P":
  if u in ("ug/L","ug P/L"):return float(v)/1000.
  if u in ("mg/L","mg P/L"):return float(v)
  raise RuntimeError("P nutrient unit unrecognized: "+repr(u))
def chemistry(w):
 res=w.copy()
 for field,unitf,kind,target in [("NTL_RESULT","NTL_UNITS","N","tn_mgL"),
                                 ("PTL_RESULT","PTL_UNITS","P","tp_mgL")]:
  valid=numeric(res[field])
  def conv(i):
   if not np.isfinite(valid.loc[i]) or valid.loc[i]<=0:return np.nan
   return unit_transform(valid.loc[i],res.loc[i,unitf],kind)
  # Avoid silently treating missing as zero.
  res[target]=pd.Series([conv(i) for i in res.index],index=res.index)
  res[target+"_rawflag"]=res[field.replace("_RESULT","_NARS_FLAG")].fillna("")
 res["log_tn"]=np.log10(res.tn_mgL)
 res["log_tp"]=np.log10(res.tp_mgL)
 cc=numeric(res.COND_RESULT)
 res["log_cond"]=np.where(cc>0,np.log10(cc),np.nan)
 return res
def site(frame,era):
 a=frame.copy()
 cols=["UID","SITE_ID","DATE_COL","VISIT_NO","PSTL_CODE","AG_ECO9","LAT_DD83",
       "LON_DD83","HUC8","YEAR","STRAH_ORD"]
 if any(c not in a for c in cols):raise RuntimeError(era+": missing physical site column")
 if a.UID.duplicated().any():raise RuntimeError("Duplicated site UID")
 a=a.set_index("UID")
 o=pd.DataFrame(index=a.index)
 o["site_id"]=a.SITE_ID
 o["huc8"]=a.HUC8
 o["state"]=a.PSTL_CODE
 o["eco9"]=a.AG_ECO9
 for target,src in [("lat","LAT_DD83"),("lon","LON_DD83"),("stream_order","STRAH_ORD"),("year","YEAR")]:
  o[target]=numeric(a[src])
 # Original native day of year from date; no absent date silently imputed.
 date=pd.to_datetime(a.DATE_COL,errors="coerce")
 o["season_sin"]=np.sin(2*np.pi*date.dt.dayofyear/365.25)
 o["season_cos"]=np.cos(2*np.pi*date.dt.dayofyear/365.25)
 return o
def chla_frame(x,era):
 z=x.copy()
 if era=="2018_19":
  maps=cross[era]["uid_bridge_manifest"]["chla"]
  z["canonical_uid"]=z.UID.map(maps)
 else:z["canonical_uid"]=z.UID
 z=z[z.canonical_uid.notna()].copy()
 if z.canonical_uid.duplicated().any():raise RuntimeError("Chla many-to-one UID source mapping")
 z=z.set_index("canonical_uid")
 z["chla_ug_cm2"]=np.where(z.RESULT_UNITS.eq("ug/cm2"),numeric(z.RESULT),np.nan)
 z["chla_flag"]=z.NARS_FLAG
 return z[["chla_ug_cm2","chla_flag"]]
def afdm_frame(x):
 z=x.copy()
 if z.UID.duplicated().any():raise RuntimeError("2023/24 AFDM UID duplicates")
 z=z.set_index("UID")
 z["afdm_mg_cm2"]=np.where(z.RESULT_UNITS.eq("mg/cm2"),numeric(z.RESULT),np.nan)
 z["afdm_flag"]=z.NARS_FLAG
 return z[["afdm_mg_cm2","afdm_flag"]]
def ept_frame(x):
 x=x.copy()
 x=x[x.UID.isin(frozen18)&x.SAMPLE_TYPE.eq("BERW")].copy()
 x["c"]=numeric(x.TOTAL300)
 x=x[x.c.gt(0)&np.isfinite(x.c)]
 x["ept"]=x.ORDER.str.strip().str.lower().isin(["ephemeroptera","plecoptera","trichoptera"])
 x["ept_c"]=x.c*x.ept.astype(int)
 g=x.groupby("UID",sort=True)
 y=pd.DataFrame({"count_total":g.c.sum(),"count_ept":g.ept_c.sum(),"n_taxon_rows":g.size()})
 y["ept_fraction"]=y.count_ept/y.count_total
 if len(y)!=1221:raise RuntimeError("Previously frozen pre-effect EPT sample changed")
 return y
def prepare(era):
 tag="18" if era=="2018_19" else "23"
 a=site(D[tag+"_sites"],era)
 w=chemistry(uniq(D[tag+"_water"],tag+"_water"))
 c=chla_frame(D[tag+"_chla"],era)
 panel=a.join(w[["log_tn","log_tp","log_cond","tn_mgL","tp_mgL","tn_mgL_rawflag","tp_mgL_rawflag"]],how="inner")
 panel=panel.join(c,how="inner")
 if tag=="18":
  panel=panel.join(ept_frame(D["18_benthos"]),how="inner")
  panel=panel[panel.index.isin(frozen18)].copy()
 else:
  panel=panel.join(afdm_frame(D["23_biomass"]),how="inner")
  eligible_source=set(cross[era]["uid_bridge_manifest"]["biomass"].values())
  eligible_source &=set(cross[era]["uid_bridge_manifest"]["chla"].values())
  eligible_source &=set(cross[era]["uid_bridge_manifest"]["water"].values())
  panel=panel[panel.index.isin(eligible_source)].copy()
 panel["log_chla"]=np.where(panel.chla_ug_cm2>0,np.log10(panel.chla_ug_cm2),np.nan)
 if tag=="23":panel["log_afdm"]=np.where(panel.afdm_mg_cm2>0,np.log10(panel.afdm_mg_cm2),np.nan)
 return panel
def folds(eligible):
 return np.array([int(hashlib.sha256(str(x).encode()).hexdigest()[:8],16)%5 for x in eligible.huc8])
def matrices(df,targets):
 M=pd.concat([df[["lat","lon","stream_order","year","season_sin","season_cos","log_cond"]].astype(float),
              pd.get_dummies(df[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 M.insert(0,"intercept",1.)
 x0=M.to_numpy(dtype=float)
 xs={"baseline":x0,"TP":np.column_stack([x0,df.log_tp.to_numpy(float)]),
     "TN":np.column_stack([x0,df.log_tn.to_numpy(float)]),
     "TN+TP":np.column_stack([x0,df.log_tn.to_numpy(float),df.log_tp.to_numpy(float)])}
 if "ept_fraction" in targets:
  xs["TN+TP+chla"]=np.column_stack([xs["TN+TP"],df.log_chla.to_numpy(float)])
 return xs
def ridge(X,Y,fold):
 outs=np.empty(len(Y))
 for f in sorted(set(fold)):
  train=fold!=f;test=~train
  A=X.astype(float).copy()
  mu=A[train,1:].mean(0);sd=A[train,1:].std(0)
  sd=np.where(sd>1e-10,sd,1)
  A[:,1:]=(A[:,1:]-mu)/sd
  Z=A[train];T=A[test]
  reg=np.eye(Z.shape[1])*10;reg[0,0]=0
  outs[test]=T@np.linalg.solve(Z.T@Z+reg,Z.T@Y[train])
 return outs
def ols_cluster(X,y,blocks,field_labels):
 # Per +1 SD raw logTN/logTP/chla, with fixed original regional dummies.
 p=X.copy().astype(float)
 names=field_labels
 for index in names:
  sd=np.std(p[:,index])
  if sd>1e-8:p[:,index]=(p[:,index]-np.mean(p[:,index]))/sd
 b=np.linalg.lstsq(p,y,rcond=None)[0];u=y-p@b
 meat=np.zeros((p.shape[1],p.shape[1]))
 for g in set(blocks):
  k=(blocks==g);v=p[k].T@u[k];meat+=np.outer(v,v)
 inv=np.linalg.pinv(p.T@p)
 v=inv@meat@inv
 G=len(set(blocks));N=len(y);K=np.linalg.matrix_rank(p)
 if G>1 and N>K:v*=G/(G-1)*(N-1)/(N-K)
 return {str(index):{"beta_per_1SD":float(b[index]),
     "HUC8_cluster_descriptive95":[float(b[index]-1.96*np.sqrt(max(v[index,index],0))),
                                    float(b[index]+1.96*np.sqrt(max(v[index,index],0)))]}
          for index in names}
def fit_one(era,p,strict=True):
 targets=["ept_fraction","log_chla"] if era=="2018_19" else ["log_afdm","log_chla"]
 columns=["log_tn","log_tp","log_cond","lat","lon","stream_order","year","season_sin","season_cos"]+targets
 qa={"joint_IDs":len(p)}
 if strict:
  keep=p.chla_flag.eq("") & p.tn_mgL_rawflag.eq("") & p.tp_mgL_rawflag.eq("")
  if era=="2023_24":keep=keep&p.afdm_flag.eq("")
 else:
  allowed=lambda z:~z.astype(str).str.contains("Q|NR|ND",regex=True)
  keep=allowed(p.chla_flag)&allowed(p.tn_mgL_rawflag)&allowed(p.tp_mgL_rawflag)
  if era=="2023_24":keep=keep&allowed(p.afdm_flag)
 p=p[keep].dropna(subset=columns).copy()
 p=p[(p.huc8!="")&(p.state!="")&(p.eco9!="")]
 if era=="2018_19":minimum=400
 else:minimum=850
 qa.update({"source_quality_n":len(p),"n_HUC8":p.huc8.nunique(),
             "n_states":p.state.nunique(),"n_sites":p.site_id.nunique(),
             "n_eco9":p.eco9.nunique(),"flags_strict":strict,
             "negative_or_zero_afdm_excluded":True,"no_numeric_missing_imputation":True})
 if len(p)<minimum or p.huc8.nunique()<100 or p.state.nunique()<20:
  qa["result"]="SOURCE/NUMERIC SUPPORT GATE FAIL - NO fitted coefficients"
  return qa
 X=matrices(p,targets);fold=folds(p);blocks=p.huc8.astype(str).to_numpy()
 outputs={}
 for ykey in targets:
  Y=p[ykey].to_numpy(dtype=float)
  base=float(np.mean((Y-ridge(X["baseline"],Y,fold))**2))
  models={};modelnames=["baseline","TN","TP","TN+TP"]
  if ykey=="ept_fraction":modelnames+=["TN+TP+chla"]
  for name in modelnames:
   pred=ridge(X[name],Y,fold)
   err=float(np.mean((Y-pred)**2))
   models[name]={"MSE":err,"heldout_MSE_gain_vs_baseline_pct":100*(base-err)/base}
  xfull=X["TN+TP"+("+chla" if ykey=="ept_fraction" else "")]
  ids=list(range(xfull.shape[1]-3,xfull.shape[1])) if ykey=="ept_fraction" else list(range(xfull.shape[1]-2,xfull.shape[1]))
  slopevals=ols_cluster(xfull,Y,blocks,ids)
  slopenames=["log_tn","log_tp","log_chla"] if ykey=="ept_fraction" else ["log_tn","log_tp"]
  outputs[ykey]={"n":len(p),"observed_mean":float(Y.mean()),"observed_SD":float(Y.std()),
                 "huc8_disjoint_5fold":models,
                 "partial_slopes":{name:slopevals[str(index)] for name,index in zip(slopenames,ids)}}
 qa["result"]="COMPUTED - OBSERVATIONAL ONLY";qa["responses"]=outputs
 return qa
p18=prepare("2018_19");p23=prepare("2023_24")
# Ephemeral per-visit analysis frames for independently implemented numerical QA.
# Neither rows nor site coordinates are uploaded as GitHub Actions artifacts.
p18.to_csv(ROOT/"AUDIT_ONLY_2018_19_PANEL.csv",index_label="UID")
p23.to_csv(ROOT/"AUDIT_ONLY_2023_24_PANEL.csv",index_label="UID")
out={
 "schema":"gesui_NRSA_2018_19_2023_24_SAME_UID_NUTRIENT_BIOMASS_EPT_EXPLORATORY_v1",
 "source_SHAs":{n:s for n,(u,s) in EXPECTED.items()},
 "2018_benthos_pre_effect_original_uid_sha":MANIFEST_SHA,
 "2018_chla_bridge_from_pre_effect_source":True,
 "2023_original_source_before_bio_effect_uid_bridge":True,
 "not_WWTP_sewage_exposure":True,
 "not_ecosystem_service_or_primary_production":True,
 "not_causal":True,
 "not_temporal_change_or_repeated_site":True,
 "years":{"2018_19":{"strict":fit_one("2018_19",p18,True),
                     "quality_flag_sensitivity":fit_one("2018_19",p18,False)},
          "2023_24":{"strict":fit_one("2023_24",p23,True),
                     "quality_flag_sensitivity":fit_one("2023_24",p23,False)}},
 "avoid_interpretation":["2018 AFDM mg/L extracted volume is NOT comparable to 2023 AFDM mg/cm2;",
                         "2018 EPT=standardized counted larval composition not biomass/abundance;",
                         "Chlorophyll-a ug/cm2 is a standing concentration, not instantaneous productivity;",
                         "2023/24 original benthic count publisher link 404; do not fabricate matching EPT;",
                         "HUC8 grouped forecast and HUC8 cluster intervals do NOT remove confounding."]}
t=ROOT/"OBSERVED_NRSA_NUTRIENT_PERIPHYTON_EPT_TWO_ERAS_NO_WWTP_CAUSAL.json"
t.write_text(json.dumps(out,sort_keys=True,indent=2))
print("NRSA_BIOMASS_EFFECT_RESULT_START")
print(json.dumps(out,sort_keys=True,default=str))
print("NRSA_BIOMASS_EFFECT_RESULT_END")
