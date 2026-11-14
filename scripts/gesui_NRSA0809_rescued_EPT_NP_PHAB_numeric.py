#!/usr/bin/env python3
"""Actual EPA 2008/09 source-rescued wadeable BERWW EPT counted fraction
from original 300-organism counts × authentic TAXA_ID/ORDER. Compare native
measured ambient stream TN/TP, physical habitat and original EPT taxa/MMI.
NOT sewage N/P load, 2008 periphyton pigment, emergent insect biomass,
or causal ecosystem services. Old published-only EPT_PIND n286 FAIL preserved.
"""
from __future__ import annotations
import collections,hashlib,json,os,pathlib,re,time,urllib.request
import numpy as np,pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge
import statsmodels.api as sm
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui0809_reconstructed_EPT"
ROOT.mkdir(parents=True,exist_ok=True)
first=json.loads(pathlib.Path(os.environ["ORIGINAL_2008_SOURCE"]).read_text())
taxqa=json.loads(pathlib.Path(os.environ["ORIGINAL_TAXON_QA"]).read_text())
fidelity=json.loads(pathlib.Path(os.environ["ORIGINAL_EPT_FIDELITY"]).read_text())
if first["schema"]!="gesui_NRSA0809_AUTHENTIC_EPA_WATER_BENTHOS_PHAB_MMI_SOURCE_ONLY_v1":
 raise RuntimeError("Original 2008 five-source schema changed")
if taxqa["schema"]!="gesui_2008_original_source_TAXA_ID_ORDER_BERWW_TOTAL300_bridge_v1":
 raise RuntimeError("Original taxonomy gate not source frozen")
if fidelity["schema"]!="gesui_EPA0809_original_BERWW_TOTAL300_EPT_ORDER_with_legacy_metric_fidelity_v1":
 raise RuntimeError("EPA published EPT source fidelity not independently audited")
if not fidelity["true_predeclared_fidelity_gate_pass"]:raise RuntimeError("EPA original EPT raw/published discrepancy gate failed")
allowed=set(fidelity["eligible_native_UIDs_for_reconstructed_share"])
if len(allowed)!=1232:raise RuntimeError("Independent 2008 source-native raw taxon EPT candidate n1232 changed")
META=first["original_data_sources"]
src={
 "siteinfo":META["siteinfo"]["url"],
 "water":META["water"]["url"],
 "phab":META["phab"]["url"],
 "bentcount":META["bentcount"]["url"],
 "bentmmi":META["bentmmi"]["url"],
 "taxa":"https://www.epa.gov/sites/default/files/2016-06/nrsa_0809_benttaxa.csv"
}
def get(k):
 p=ROOT/(k+".csv")
 for i in range(3):
  try:
   req=urllib.request.Request(src[k],headers={"User-Agent":"gesui-public-EPA2008-raw-EPT-count-NP-v1"})
   with urllib.request.urlopen(req,timeout=120) as r,p.open("wb") as out:
    while z:=r.read(1<<20):out.write(z)
   h=hashlib.sha256(p.read_bytes()).hexdigest()
   old=taxqa["original_EPA_ID_order_mapping_QA"]["taxonomy_csv_sha256"] if k=="taxa" else META[k]["source_SHA256"]
   if h!=old:raise RuntimeError("Original source SHA changed "+k)
   break
  except Exception:
   p.unlink(missing_ok=True)
   if i==2:raise
   time.sleep(2**i)
 df=pd.read_csv(p,dtype=str,keep_default_na=False,encoding="latin1",low_memory=False)
 p.unlink()
 df.columns=df.columns.str.upper().str.strip()
 return df
T=get("taxa");C=get("bentcount");S=get("siteinfo");W=get("water");P=get("phab");M=get("bentmmi")
# Original field taxonomy mapping exactly one ORDER per official taxon ID.
T["key"]=T.TAXA_ID.str.strip().str.upper()
T["ord"]=T.ORDER.str.strip().str.upper()
tx=T[T.key.ne("")].drop_duplicates(["key","ord"])
if (tx.groupby("key").ord.nunique()>1).any():raise RuntimeError("Source 2008 taxonomy contradicts own source")
byorder=tx.drop_duplicates("key").set_index("key").ord
C=C[C.UID.isin(allowed)&C.SAMPLE_TYPE.str.strip().str.upper().eq("BERWW")&
    C.SAMPLE_CAT.str.strip().str.upper().eq("P")].copy()
C["N"]=pd.to_numeric(C.TOTAL300,errors="coerce")
C=C[C.N.gt(0)].copy()
C["key"]=C.TAXA_ID.str.strip().str.upper()
C["order"]=C.key.map(byorder).fillna("")
C["EPT"]=C.order.isin(["EPHEMEROPTERA","PLECOPTERA","TRICHOPTERA"])
C["known"]=C.order.ne("")
G=C.groupby("UID")
n=G.N.sum()
nE=(C.N*C.EPT.astype(float)).groupby(C.UID).sum()
nK=(C.N*C.known.astype(float)).groupby(C.UID).sum()
if set(n.index)!=allowed:raise RuntimeError("Raw original EPT numerator count missing some pre-audited eligible UIDs")
ratio=(nE/n).reindex(sorted(allowed))
known_ratio=(nK/n).reindex(sorted(allowed))
if (known_ratio<.95-1e-10).any():raise RuntimeError("Original 2008 order mapping positive count coverage less than source-frozen 95%")
Traw=S[S.UID.isin(allowed)].set_index("UID").sort_index()
Ph=P[P.UID.isin(allowed)].set_index("UID").sort_index()
Mi=M[M.UID.isin(allowed)].set_index("UID").sort_index()
Wa=W[W.UID.isin(allowed)&W.SAMPLE_TYPE.str.strip().str.upper().eq("CHEMW")&
    W.SAMPLE_CAT.str.strip().str.upper().eq("P")].set_index("UID").sort_index()
for name,x in (("SiteInfo",Traw),("PHab",Ph),("BENT MMI",Mi),("CHEMW",Wa)):
 if len(x)!=1232 or x.index.duplicated().any():raise RuntimeError("Source original 2008 same native UID from 1232 is not 1:1 "+name)
for name,x in (("PHab",Ph),("BENT MMI",Mi),("CHEMW",Wa)):
 if not x.SITE_ID.equals(Traw.SITE_ID):raise RuntimeError("Source native UID but wrong SITE_ID "+name)
if not Traw.SITE_CLASS.str.upper().eq("PROB").all():raise RuntimeError("Original 2008 probability site eligibility not met")
if not Mi.BENT_SAMPLE_TYPE.str.upper().eq("BERWW").all():
 raise RuntimeError("Original BERWW raw fraction doesn't match native published wadeable MMI")
for name,x in (("PHab",Ph),("BENT MMI",Mi),("CHEMW",Wa)):
 dat=pd.to_datetime(Traw.DATE_COL,errors="coerce",format="mixed")
 compare=pd.to_datetime(x.DATE_COL,errors="coerce",format="mixed")
 discrepancy=dat.notna()&compare.notna()&dat.dt.normalize().ne(compare.dt.normalize())
 if discrepancy.any():raise RuntimeError("Original native visit date mismatch "+name)
def N(X,k):return pd.to_numeric(X[k],errors="coerce")
# True source ordinal categories, NOT original missingness.
order_labels={"1ST":1,"2ND":2,"3RD":3,"4TH":4,"5TH":5,"6TH":6,"7TH":7}
original_ordinal=Traw.STRAHLERORDER.str.strip().str.upper()
actualorder=original_ordinal.map(order_labels)
if int(actualorder.notna().sum())<1200:raise RuntimeError("Actual Strahler order labels source unsupported")
if not set(original_ordinal.unique())<=set(order_labels)|{"8TH+"}:
 raise RuntimeError("Original 2008 uncatalogued Strahler source labels: "+repr(original_ordinal.unique()[:10]))
D=pd.DataFrame(index=Traw.index)
D["state"]=Traw.STATE.str.strip()
D["eco9"]=Traw.AGGR_ECO9_2015.str.strip()
D["huc8"]=Traw.HUC8.str.strip()
D["site_id"]=Traw.SITE_ID.str.strip()
D["lat"]=N(Traw,"LAT_DD83");D["lon"]=N(Traw,"LON_DD83")
D["year"]=N(Traw,"YEAR");D["stream_order"]=actualorder
dates=pd.to_datetime(Traw.DATE_COL,errors="coerce",format="mixed")
D["season_sin"]=np.sin(2*np.pi*dates.dt.dayofyear/365.2425)
D["season_cos"]=np.cos(2*np.pi*dates.dt.dayofyear/365.2425)
D["log_cond"]=np.log10(N(Wa,"COND").where(lambda a:a>0))
D["log_tn"]=np.log10((N(Wa,"NTL")/1000).where(lambda a:a>0))
D["log_tp"]=np.log10((N(Wa,"PTL")/1000).where(lambda a:a>0))
D["strict_chem"]=Wa.NTL_ALERT.str.strip().eq("")&Wa.PTL_ALERT.str.strip().eq("")&Wa.COND_ALERT.str.strip().eq("")
for v in ("XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","XWIDTH"):D[v]=N(Ph,v)
D["log_width"]=np.log10(D.XWIDTH.where(D.XWIDTH>0))
D["ept_share_raw_count_lower_bound"]=ratio
D["ept_share_source_upper_bound"]=((nE+(n-nK))/n).reindex(D.index)
D["source_order_mapped_fraction"]=known_ratio
D["ept_native_published_taxon_count"]=N(Mi,"EPT_NTAX")
D["log1p_ept_taxa"]=np.log1p(D.ept_native_published_taxon_count.where(D.ept_native_published_taxon_count>=0))
D["benthic_mmi_score"]=N(Mi,"MMI_BENT")
D["mmi_has_observed_berww"]=Mi.SAMPLED_BENT.str.upper().str.strip().eq("YES")
if (D.ept_share_raw_count_lower_bound.dropna()<0).any() or (D.ept_share_raw_count_lower_bound.dropna()>1).any():
 raise RuntimeError("Original source EPT counted fraction invalid")
if D.ept_native_published_taxon_count.dropna().lt(0).any():raise RuntimeError("Native EPT_NTAX cannot be negative")
if D.benthic_mmi_score.dropna().lt(0).any() or D.benthic_mmi_score.dropna().gt(100).any():
 raise RuntimeError("Source original EPA benthic MMI score outside 0..100")
all_phys=["lat","lon","year","stream_order","season_sin","season_cos","log_cond",
          "log_tn","log_tp","XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]
df=D[D.strict_chem&D.mmi_has_observed_berww].copy()
pop={"original_source_count_reconstructed_native_UID":1232,
 "eligible_clean_unalerted_chem":int(D.strict_chem.sum()),
 "eligible_published_mmi_bent_with_berww":int(D.mmi_has_observed_berww.sum()),
 "numeric_common_baseline_and_phys":int(df.dropna(subset=all_phys).shape[0]),
 "native_ordinal_8th_plus_unresolved":int(original_ordinal.eq("8TH+").sum()),
 "source_ept_count_ratio_known_order_mean_pct":float(known_ratio.mean()*100),
 "sample_2008_2009_site_counts":{str(k):int(v) for k,v in D.year.value_counts().items()}}
results={}
for target in ("ept_share_raw_count_lower_bound","log1p_ept_taxa","benthic_mmi_score"):
 sample=df.dropna(subset=all_phys+[target]).copy()
 sample=sample[(sample.huc8!="")&(sample.state!="")&(sample.eco9!="")&sample.year.isin([2008,2009])]
 if target=="ept_share_raw_count_lower_bound":
  sample.to_csv(ROOT/"AUDIT_ONLY_2008_RAW_EPT_CHEMW_PHAB_SOURCE_PANEL.csv",index_label="UID")
 qa={"n":len(sample),"HUC8":sample.huc8.nunique(),"states":sample.state.nunique(),
     "source_physical_SITE_ID":sample.site_id.nunique()}
 if len(sample)<300 or sample.huc8.nunique()<150 or sample.state.nunique()<20:
  results[target]={"source_QA":qa,"status":"SOURCE_QUALITY_NUMERIC_GATE_FAIL_NO_FALSE_NULL"}
  continue
 X0=pd.concat([sample[["lat","lon","year","stream_order","season_sin","season_cos","log_cond"]].astype(float),
               pd.get_dummies(sample[["state","eco9"]],drop_first=True,dtype=float)],axis=1)
 Xp=sample[["XEMBED","XCDENMID","XPCAN","XFC_OHV","W1_HALL","log_width"]].astype(float)
 Xn=sample[["log_tn","log_tp"]].astype(float)
 mats={"geo":X0,
       "geo+N_P":pd.concat([X0,Xn],axis=1),
       "geo+PHab":pd.concat([X0,Xp],axis=1),
       "geo+PHab+N_P":pd.concat([X0,Xp,Xn],axis=1)}
 y=sample[target].to_numpy(float)
 groups=sample.huc8.to_numpy()
 fold=np.asarray([int(hashlib.sha256(str(g).encode()).hexdigest()[:8],16)%5 for g in groups])
 if min(np.bincount(fold,minlength=5))<20:raise RuntimeError("Out of HUC8 third-cycle 2008 fold support")
 def CV(mat):
  a=mat.to_numpy(float);p=np.full(len(y),np.nan)
  for k in range(5):
   tr=fold!=k;te=~tr
   sc=StandardScaler().fit(a[tr]);z=sc.transform(a)
   p[te]=Ridge(alpha=10,solver="svd").fit(z[tr],y[tr]).predict(z[te])
  if not np.isfinite(p).all():raise RuntimeError("2008 original source missing HUC8 heldout prediction")
  return p
 preds={k:CV(v) for k,v in mats.items()}
 mse={k:float(np.mean((y-p)**2)) for k,p in preds.items()}
 gain=lambda a,b:float(100*(mse[a]-mse[b])/mse[a])
 X=mats["geo+PHab+N_P"].copy()
 for col in ("log_tn","log_tp"):
  X[col]=(X[col]-X[col].mean())/X[col].std(ddof=0)
 fit=sm.OLS(y,sm.add_constant(X.to_numpy(float),has_constant="add")).fit(
   cov_type="cluster",cov_kwds={"groups":groups,"use_correction":True})
 effect={col:{"beta_1SD_log10_observed_instream":float(fit.params[X.columns.get_loc(col)+1]),
          "HUC8_cluster_descriptive95":list(map(float,fit.conf_int()[X.columns.get_loc(col)+1]))}
         for col in ("log_tn","log_tp")}
 ids,inv=np.unique(groups,return_inverse=True)
 e0=np.bincount(inv,weights=(y-preds["geo+PHab"])**2,minlength=len(ids))
 e1=np.bincount(inv,weights=(y-preds["geo+PHab+N_P"])**2,minlength=len(ids))
 rng=np.random.default_rng(20261114)
 samples=rng.integers(0,len(ids),size=(2000,len(ids)))
 den=e0[samples].sum(1);alt=e1[samples].sum(1)
 boot=list(map(float,np.quantile(100*(den-alt)/den,[.025,.5,.975])))
 results[target]={"source_QA":qa,"status":"REAL_ORIGINAL_EPA_NUMERIC_FIT",
   "source_y_mean":float(y.mean()),"source_y_min":float(y.min()),"source_y_max":float(y.max()),
   "original_HUC8_spatial5fold_MSE":mse,
   "extra_N_P_after_PHab_gain_pct":gain("geo+PHab","geo+PHab+N_P"),
   "extra_PHab_after_N_P_gain_pct":gain("geo+N_P","geo+PHab+N_P"),
   "HUC8_bootstrap2000_fixed_OOF_extra_NP_descriptive95":boot,
   "fully_adjusted_nutrient_associations":effect}
out={"schema":"gesui_NRSA0809_original_rescued_EPT_count_share_EPT_taxon_MMI_NP_PHAB_2008_v1",
 "source original EPT total300 validated":"EPA original published EPT_PIND 279 same visits median error 0.003 pp",
 "source_fidelity_manifest_SHA256":hashlib.sha256(pathlib.Path(os.environ["ORIGINAL_EPT_FIDELITY"]).read_bytes()).hexdigest(),
 "source_original_taxon_SHA256":taxqa["original_EPA_ID_order_mapping_QA"]["taxonomy_csv_sha256"],
 "original_taxon_order_no_contradictions":True,
 "source_relative_EPT_share_lower_bound_not_absolute_insects":True,
 "not_2008_periphyton_chlorophyll_or_AFDM":True,"not_measured_WWTP_NP_delivery":True,
 "no_realized_services_or_causation":True,
 "sample_source_availability":pop,"responses":results}
p=ROOT/"NRSA0809_RECONSTRUCTED_EPT_REAL_AMBIENT_NP_PHab_NUMERIC_RESULTS.json"
p.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA0809_RECONSTRUCTED_EPT_CHEMW_PHAB_TRUE_ORIGINAL_RESULT")
print(json.dumps({"SOURCE":pop,
 "results":{target:{"QA":m["source_QA"],
  "status":m["status"],
  "extra_NP":m.get("extra_N_P_after_PHab_gain_pct"),
  "extra_PHab":m.get("extra_PHab_after_N_P_gain_pct"),
  "nutrient_coefficients":m.get("fully_adjusted_nutrient_associations"),
  "NP_boot95":m.get("HUC8_bootstrap2000_fixed_OOF_extra_NP_descriptive95")}
  for target,m in results.items()}},sort_keys=True))
