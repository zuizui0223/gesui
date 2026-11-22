#!/usr/bin/env python3
"""Independent NumPy weighted-normal-equation audit of original 4country EPT.
Does not use the original country-heldout or covariance decomposition module.
Only reuses the publisher OOXML read-only parser; all statistical arithmetic
is independently implemented from the source normalized taxon/sediment rows.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np
import gesui_external_fourcountry_xlsx_readonly_ooxml as oo

SOURCE_SHA="a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48"
SOURCE_MD5="1e4d6e37b73c663406e5567dd79d9452"
SHEETS={
 "Australia":("Australia taxa","Australia sediment",3099),
 "Brazil":("Brazil taxa","Brazil sediment",232),
 "New Zealand":("NZ taxa","Nz sediment",342),
 "UK":("UK taxa","UK Sediment",2818)}
def sha(path,algo):
 h=hashlib.new(algo)
 with path.open("rb") as f:
  for b in iter(lambda:f.read(1048576),b""):h.update(b)
 return h.hexdigest()
def equal(a,b,label):
 if not np.allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),
                    rtol=2e-9,atol=2e-10):
  raise RuntimeError("Independent publisher source/covariance QA FAILED "+label+" "+str((a,b)))
def read_publisher(path,taxa):
 outputs={}
 for country,(taxsheet,sedsheet,expected) in SHEETS.items():
  head,values=oo.read_all(path,taxsheet)
  header,sed=oo.read_all(path,sedsheet)
  if len(values)!=expected or len(sed)!=expected or head[0]!="ANALYSIS_ID":
   raise RuntimeError("Original publisher country sample grain changed")
  names=head[1:]
  if not all(name in taxa for name in names):raise RuntimeError("Taxon header not in source-only taxonomy")
  s_by_id={oo.normalize_id(r[0]):float(oo.as_float(r[1])) for r in sed}
  if len(s_by_id)!=expected:raise RuntimeError("Original source sediment duplicate IDs")
  ye=[];upper=[];s=[];ids=set()
  iEPT=[j for j,t in enumerate(names) if taxa[t]=="EPT"]
  iUNK=[j for j,t in enumerate(names) if taxa[t]=="UNKNOWN"]
  for row in values:
   uid=oo.normalize_id(row[0])
   if uid in ids or uid not in s_by_id:raise RuntimeError("Duplicate or missing publisher paired ID")
   ids.add(uid)
   v=np.array([oo.as_float(z) for z in row[1:]],float)
   if not np.isfinite(v).all() or np.any(v<0) or abs(v.sum()-1)>1e-9:
    raise RuntimeError("Publisher original row not a normalized nonnegative taxon vector")
   yy=float(v[iEPT].sum());unk=float(v[iUNK].sum())
   ye.append(yy);upper.append(yy+unk);s.append(s_by_id[uid]/100.)
  if not set(ids)==set(s_by_id):raise RuntimeError("Original source taxon/sediment ID mismatch")
  outputs[country]={"s":np.asarray(s,float),"confirmed_EPT_lower":np.asarray(ye,float),
    "all_unknown_as_EPT_upper":np.asarray(upper,float)}
 return outputs
def country_moments(group,ycol):
 out={}
 for name,d in group.items():
  s=d["s"];y=d[ycol]
  s0=float(s.mean());y0=float(y.mean())
  z=s-s0;q=y-y0
  out[name]={"n":len(y),"mu_s":s0,"mu_y":y0,
   "within_cov":float(np.dot(z,q)/len(z)),
   "within_var":float(np.dot(z,z)/len(z))}
 return out
def matrix_fit(train,ycol,alpha,kind,visit_equal=False):
 mom=country_moments(train,ycol)
 labels=sorted(train)
 sx=np.concatenate([train[c]["s"] for c in labels])
 yy=np.concatenate([train[c][ycol] for c in labels])
 if visit_equal:
  w=np.full(len(sx),3./len(sx))
 else:
  w=np.concatenate([np.full(len(train[c]["s"]),1./len(train[c]["s"])) for c in labels])
 assert abs(w.sum()-3)<1e-9
 mu_s=float(np.dot(w,sx)/w.sum());mu_y=float(np.dot(w,yy)/w.sum())
 if kind=="pooled":
  # Explicit 2x2 weighted penalized normal equation with UNPENALIZED intercept.
  X=np.column_stack((np.ones(len(sx)),sx))
  lhs=X.T@(X*w[:,None])+np.diag([0.,alpha])
  rhs=X.T@(w*yy)
  intercept,beta=np.linalg.solve(lhs,rhs)
  equal(intercept,mu_y-beta*mu_s,"source intercept unpenalized")
 elif kind=="within":
  # Subtract each training country's own source means BEFORE calculating beta.
  # The final intercept still predicts the training country-EQUAL mean.
  sg=np.concatenate([train[c]["s"]-mom[c]["mu_s"] for c in labels])
  yg=np.concatenate([train[c][ycol]-mom[c]["mu_y"] for c in labels])
  beta=float((np.dot(w,sg*yg))/(np.dot(w,sg*sg)+alpha))
  intercept=mu_y-beta*mu_s
 else:raise ValueError("Unexpected source model")
 return {"intercept":float(intercept),"beta":float(beta),"mu_s":mu_s,
    "mu_y":mu_y,"source_by_country":mom}
def expected(target,model,ycol):
 y=target[ycol];s=target["s"]
 p0=np.full(len(y),model["mu_y"])
 p1=model["intercept"]+model["beta"]*s
 ref=float(np.mean((y-p0)**2));fit=float(np.mean((y-p1)**2))
 cov=float(np.mean((s-s.mean())*(y-y.mean())))
 var_s=float(np.mean((s-s.mean())**2))
 ds=float(s.mean()-model["mu_s"]);dy=float(y.mean()-model["mu_y"])
 b=model["beta"]
 within=2*b*cov;between=2*b*ds*dy
 penalty=-b*b*(var_s+ds*ds)
 equal(within+between+penalty,ref-fit,"independent target MSE identity")
 return {"ref":ref,"full":fit,"MSE_gain_pct":100*(ref-fit)/ref,
    "target_cov":cov,"source_target_mean_s_diff":ds,"source_target_mean_y_diff":dy,
    "within_gain_component_pct":100*within/ref,
    "offset_gain_component_pct":100*between/ref,
    "penalty_gain_component_pct":100*penalty/ref}
def validate_one(expected_signed,actual,original_country,model_kind):
 equal(actual["MSE_gain_pct"],expected_signed["MSE_gain_pct"],
       original_country+" "+model_kind+" heldout MSE gain")
 equal(actual["ref"],expected_signed["MSE0"],original_country+" baseline MSE")
 equal(actual["full"],expected_signed["MSE1"],original_country+" full MSE")
 for k,field in (("within_gain_component_pct","within_target_covariance"),
                 ("offset_gain_component_pct","mean_domain_offset"),
                 ("penalty_gain_component_pct","slope_variance_penalty")):
  equal(actual[k],expected_signed["MSE_gain_components_pctpoints"][field],
        original_country+" "+model_kind+" "+k)
def main():
 p=argparse.ArgumentParser()
 for key in ("xlsx","source_taxonomy","original_decomposition","out"):
  p.add_argument("--"+key,required=True,type=pathlib.Path)
 a=p.parse_args()
 orig=json.loads(a.original_decomposition.read_text())
 tax=json.loads(a.source_taxonomy.read_text())
 if orig["schema"]!="gesui_external4country_EPT_country_within_between_source_gain_decomposition_v1":
  raise RuntimeError("Signed original 4country decomposition schema changed")
 if orig["original_publisher_Figshare_SHA256"]!=SOURCE_SHA or tax["Figshare_original_XLSX_SHA256"]!=SOURCE_SHA:
  raise RuntimeError("Signed taxonomy or source Figshare SHA mismatch")
 if orig["original_hierarchy_taxonomy_sha256"]!=sha(a.source_taxonomy,"sha256"):
  raise RuntimeError("Signed original taxonomy source digest mismatch")
 if sha(a.xlsx,"sha256")!=SOURCE_SHA or sha(a.xlsx,"md5")!=SOURCE_MD5:
  raise RuntimeError("Original publisher XLSX source SHA/MD5 changed")
 classification={r["original_label"]:r["state"] for r in tax["taxon_label_source_lookup"]}
 if len(classification)!=302:raise RuntimeError("Original taxonomy label census changed")
 data=read_publisher(a.xlsx,classification)
 checks={}
 for bound in ("confirmed_EPT_lower","all_unknown_as_EPT_upper"):
  checks[bound]={}
  for domain,df in (
      ("full_original_0to100_visual_fines_domain",data),
      ("common_original_0to50_visual_fines_domain",
       {k:{f:v[old["s"]<=0.5] for f,v in old.items()} for k,old in data.items()})):
   signed=orig["results"][bound][domain]
   checks[bound][domain]={}
   for alpha_str,country in signed.items():
    alpha=float(alpha_str)
    checks[bound][domain][alpha_str]={}
    for name in sorted(df):
     train={k:v for k,v in df.items() if k!=name}
     if len(train)!=3:raise RuntimeError("Country-heldout source data split changed")
     src_p=matrix_fit(train,bound,alpha,"pooled")
     src_w=matrix_fit(train,bound,alpha,"within")
     src_v=matrix_fit(train,bound,alpha,"pooled",visit_equal=True)
     result_p=expected(df[name],src_p,bound)
     result_w=expected(df[name],src_w,bound)
     result_v=expected(df[name],src_v,bound)
     record=country[name]
     validate_one(record["heldout_pooled_country_equal"],result_p,name,"pool")
     validate_one(record["heldout_within_country_train_only"],result_w,name,"within")
     equal(src_p["beta"],record["train_pooled_beta"],name+" pooled beta")
     equal(src_w["beta"],record["train_within_country_only_beta"],name+" within beta")
     if alpha==10:
      validate_one(record["heldout_visit_equal_train_only"],result_v,name,"visit")
      equal(src_v["beta"],record["train_visit_equal_beta"],name+" visit equal beta")
     # Independent source within/between covariance identity using no statsmodels/sklearn.
     moms=src_p["source_by_country"]
     ks=sorted(moms)
     w_c=float(sum(moms[k]["within_cov"] for k in ks)/3)
     b_c=float(sum((moms[k]["mu_s"]-src_p["mu_s"])*(moms[k]["mu_y"]-src_p["mu_y"]) for k in ks)/3)
     mem=record["original_training_source_covariance_decomposition"]
     equal(w_c,mem["train_within_covariance_EPT_sediment"],name+" source within covariance")
     equal(b_c,mem["train_between_country_mean_covariance"],name+" source between covariance")
     equal(result_p["target_cov"],record["heldout_pooled_country_equal"]["target_within_country_covariance"],
           name+" target within covariance")
     checks[bound][domain][alpha_str][name]={
       "all_signed_model_and_exact_identity_reproduced":True,
       "alpha":alpha,
       "heldout_n":len(df[name]["s"]),
       "pooled_beta":src_p["beta"],
       "within_beta":src_w["beta"],
       "pooled_MSE_gain_pct":result_p["MSE_gain_pct"],
       "within_only_MSE_gain_pct":result_w["MSE_gain_pct"]}
 output={
  "schema":"gesui_external_fourcountry_within_between_heldout_manual_weighted_WLS_independent_QA_v1",
  "publisher_original_XLSX_SHA256":SOURCE_SHA,
  "original_signed_decomposition_sha256":sha(a.original_decomposition,"sha256"),
  "original_hierarchy_taxonomy_sha256":sha(a.source_taxonomy,"sha256"),
  "original_sample_count":6491,"all_country_domains":4,
  "all_original_source_country_alpha_domain_model_gains_and_covariance_decompositions_reproduced":True,
  "manual_2x2_weighted_normal_equations_no_original_model_helper":True,
  "no_site_independent_inference":True,
  "not_original_EPA_HUC8_history_test":True,
  "result_checks":checks}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("FOURCOUNTRY_WITHIN_BETWEEN_INDEPENDENT_MANUAL_QA",json.dumps({
  "all":True,"source_n":6491,
  "full_country_alpha10":{country:{
    "pooled_gain":q["pooled_MSE_gain_pct"],"within_only_gain":q["within_only_MSE_gain_pct"],
    "pooled_beta":q["pooled_beta"],"within_beta":q["within_beta"]}
   for country,q in checks["confirmed_EPT_lower"]["full_original_0to100_visual_fines_domain"]["10"].items()},
  "full_country_alpha0":{country:{
    "pooled_gain":q["pooled_MSE_gain_pct"],"within_only_gain":q["within_only_MSE_gain_pct"]}
   for country,q in checks["confirmed_EPT_lower"]["full_original_0to100_visual_fines_domain"]["0"].items()}
  },sort_keys=True),flush=True)
if __name__=="__main__":main()
