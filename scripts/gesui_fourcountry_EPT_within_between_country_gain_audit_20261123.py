#!/usr/bin/env python3
"""Decompose original 4-country fine-sediment EPT heldout gain into
within-country target covariance and between-country mean-offset correction.

No new taxonomy selection; use exact original Figshare 6,491 normalized
relative EPT rows and source-signed GBIF hierarchy EPT labels. Each heldout
country is never used to fit its slope or intercept. Descriptive ONLY.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib
import numpy as np
import gesui_external_fourcountry_fixed_country_holdout_EPT_20261122 as original
from gesui_external_fourcountry_source_grain_20261122 import ORIGINAL_SHA,ORIGINAL_MD5,digest,SHEETS

COUNTRIES=tuple(sorted(SHEETS))
ALPHAS=(0.,1.,10.,100.)
def near(a,b,label,rtol=1e-10,atol=1e-11):
 if not np.allclose(a,b,rtol=rtol,atol=atol):
  raise RuntimeError("Signed exact EPA? External publisher source replay mismatch "+label+": "+str((a,b)))
def cov(x,y):
 return float(np.mean((x-np.mean(x))*(y-np.mean(y))))
def var(x):
 return float(np.mean((x-np.mean(x))**2))
def source_moments(train,field):
 countries=sorted(train)
 c={}
 for k in countries:
  s=train[k]["sediment_fraction"];y=train[k][field]
  if len(s)==0 or not np.isfinite(s).all() or not np.isfinite(y).all():
   raise RuntimeError("Original country source missing finite sample")
  c[k]={"original_source_samples":len(s),"mean_s":float(s.mean()),
   "mean_y":float(y.mean()),"within_cov":cov(s,y),
   "within_var_s":var(s),"within_var_y":var(y)}
 mu_s=float(np.mean([r["mean_s"] for r in c.values()]))
 mu_y=float(np.mean([r["mean_y"] for r in c.values()]))
 w_cov=float(np.mean([r["within_cov"] for r in c.values()]))
 w_var=float(np.mean([r["within_var_s"] for r in c.values()]))
 b_cov=float(np.mean([(r["mean_s"]-mu_s)*(r["mean_y"]-mu_y) for r in c.values()]))
 b_var=float(np.mean([(r["mean_s"]-mu_s)**2 for r in c.values()]))
 return {"train_countries":c,"n_train_countries":3,"train_source_mu_s":mu_s,
  "train_source_mu_y":mu_y,"train_within_covariance_EPT_sediment":w_cov,
  "train_between_country_mean_covariance":b_cov,
  "train_within_sediment_variance":w_var,
  "train_between_country_mean_sediment_variance":b_var,
  "train_pooled_total_covariance":w_cov+b_cov,
  "train_pooled_total_sediment_variance":w_var+b_var}
def slope(moments,alpha,kind):
 num=3*(moments["train_within_covariance_EPT_sediment"]+
       (moments["train_between_country_mean_covariance"] if kind=="pooled" else 0.))
 den=3*(moments["train_within_sediment_variance"]+
       (moments["train_between_country_mean_sediment_variance"] if kind=="pooled" else 0.))+alpha
 if not den>0:raise RuntimeError("No source sediment gradient for source-only fixed slope")
 return float(num/den)
def visit_equal_fit(train,field,alpha):
 S=np.concatenate([train[c]["sediment_fraction"] for c in sorted(train)])
 Y=np.concatenate([train[c][field] for c in sorted(train)])
 mu_s=float(S.mean());mu_y=float(Y.mean())
 # Total source weight 3 to match country-balanced source precision exactly.
 beta=3*cov(S,Y)/(3*var(S)+alpha)
 return {"train_source_mu_s":mu_s,"train_source_mu_y":mu_y,
   "slope":float(beta),"train_visit_weight_equal_total_weight_3":True}
def model_predictions(test,mu_s,mu_y,beta,field):
 s=test["sediment_fraction"];y=test[field]
 p0=np.full(len(y),mu_y,dtype=float)
 p1=mu_y+beta*(s-mu_s)
 e0=float(np.mean((y-p0)**2))
 e1=float(np.mean((y-p1)**2))
 if e0<=0 or not np.isfinite(e0+e1):raise RuntimeError("Source target baseline variance absent")
 delta_s=float(s.mean()-mu_s);delta_y=float(y.mean()-mu_y)
 cov_t=cov(s,y);var_t=var(s)
 within=2*beta*cov_t
 between=2*beta*delta_s*delta_y
 penalty=-beta*beta*(var_t+delta_s*delta_s)
 improvement=e0-e1
 near(within+between+penalty,improvement,"Exact per-country MSE gain decomposition",rtol=1e-9,atol=1e-13)
 return {"heldout_country_sample_n":len(y),"target_source_mean_s":float(s.mean()),
  "target_source_mean_y":float(y.mean()),
  "target_within_country_covariance":cov_t,
  "target_within_sediment_variance":var_t,
  "target_descriptive_within_country_slope":float(cov_t/var_t) if var_t>0 else None,
  "target_minus_train_mean_s":delta_s,"target_minus_train_mean_y":delta_y,
  "MSE0":e0,"MSE1":e1,
  "MSE_gain_pct":float(100*improvement/e0),
  "MSE_gain_components_pctpoints":{
    "within_target_covariance":float(100*within/e0),
    "mean_domain_offset":float(100*between/e0),
    "slope_variance_penalty":float(100*penalty/e0)},
  "source_variance_formula_exact":True}
def revalidate_original(d,source_signed):
 if source_signed.get("schema")!="gesui_external_original_fourcountry_relative_EPT_country_heldout_DESCRIPTIVE_v1":
  raise RuntimeError("Original publisher fixed 4-country model schema changed")
 if source_signed["original_loughborough_figshare_SHA256"]!=ORIGINAL_SHA:
  raise RuntimeError("Original publisher source SHA mismatch")
 replay={}
 for bound,field in (("confirmed_EPT_lower","known_EPT_fraction"),
                     ("all_unknown_as_EPT_upper","unknown_bound_EPT_upper_fraction")):
  replay[bound]={}
  for domain,groups in (("full_original_0to100_visual_fines_domain",d),
                        ("common_original_0to50_visual_fines_domain",original.common_domain(d))):
   saved=source_signed["results"][bound][domain]
   fresh=original.compare(groups,field)
   for country in COUNTRIES:
    a=fresh[country];b=saved[country]
    for f in ("heldout_original_samples","heldout_MSE_train_country_mean_baseline",
              "heldout_MSE_add_linear_ridge10_fines","heldout_MSE_gain_pct",
              "observed_target_mean_relative_EPT","observed_target_mean_visual_fines_fraction"):
     near(a[f],b[f],"prior original model "+bound+"/"+domain+"/"+country+"/"+f)
    for f in ("train_3country_equal_country_baseline_mean_EPT_fraction",
              "train_3country_equal_country_mean_fines_fraction",
              "Ridge10_slope_EPT_fraction_per_fines_0to1"):
     near(a["model"][f],b["model"][f],"prior original model "+bound+"/"+domain+"/"+country+"/"+f)
   replay[bound][domain]={"countries_exact_numeric_reproduced":4,
     "all_source_sample_counts_and_original_MSE_slopes_reproduced":True}
 return replay
def fit_all(d,field,alpha,domain,include_offset=True):
 out={}
 for target_country in COUNTRIES:
  train={c:d[c] for c in COUNTRIES if c!=target_country}
  mom=source_moments(train,field)
  beta_pool=slope(mom,alpha,"pooled")
  beta_within=slope(mom,alpha,"within")
  trainv=visit_equal_fit(train,field,alpha)
  target=d[target_country]
  p=model_predictions(target,mom["train_source_mu_s"],mom["train_source_mu_y"],beta_pool,field)
  w=model_predictions(target,mom["train_source_mu_s"],mom["train_source_mu_y"],beta_within,field)
  vs=model_predictions(target,trainv["train_source_mu_s"],trainv["train_source_mu_y"],trainv["slope"],field)
  if alpha==10.:
   # Strict cross-check against the original signed M_pool model.
   originalfit=original.fit_three(train,field)
   near(beta_pool,originalfit["Ridge10_slope_EPT_fraction_per_fines_0to1"],
        "alpha10 original pooled exact source beta "+target_country)
   near(mom["train_source_mu_s"],originalfit["train_3country_equal_country_mean_fines_fraction"],"original mean-s")
   near(mom["train_source_mu_y"],originalfit["train_3country_equal_country_baseline_mean_EPT_fraction"],"original mean-y")
  out[target_country]={
   "heldout_source_country":target_country,
   "source_sample_n":len(target[field]),
   "alpha":alpha,
   "train_pooled_beta":beta_pool,
   "train_within_country_only_beta":beta_within,
   "train_visit_equal_beta":trainv["slope"],
   "original_training_source_covariance_decomposition":mom if include_offset else None,
   "heldout_pooled_country_equal":p,
   "heldout_within_country_train_only":w,
   "heldout_visit_equal_train_only":vs if alpha==10. else None,
   "predeclared_within_country_signal_direction_agrees":(
      beta_within*float(p["target_within_country_covariance"])>0),
   "original_pool_gain_mean_offset_abs_gt_within_abs":(
      abs(p["MSE_gain_components_pctpoints"]["mean_domain_offset"])>
      abs(p["MSE_gain_components_pctpoints"]["within_target_covariance"])),
   "no_independent_site_id":True}
 return out
def main():
 ap=argparse.ArgumentParser()
 for k in ("xlsx","source_gate","source_taxonomy","original_result","out"):
  ap.add_argument("--"+k,required=True,type=pathlib.Path)
 args=ap.parse_args()
 signed=json.loads(args.source_gate.read_text())
 tax=json.loads(args.source_taxonomy.read_text())
 old=json.loads(args.original_result.read_text())
 if signed["country_heldout_descriptive_GO"] is not True or signed["all_four_country_unknown_taxa_mass_GO"] is not True:
  raise RuntimeError("Source-only original taxon uncertainty GO failed; stop analytical model")
 if tax["schema"]!="gesui_original_McKenzie_4country_taxonomy_hierarchy_source_amendment_v2":
  raise RuntimeError("Original hierarchy source taxonomy no longer valid")
 if signed["taxonomy_original_source_json_sha256"]!=hashlib.sha256(args.source_taxonomy.read_bytes()).hexdigest():
  raise RuntimeError("Signed source-only taxonomy digest mismatch")
 if old["source_hierarchy_taxonomy_sha256"]!=hashlib.sha256(args.source_taxonomy.read_bytes()).hexdigest():
  raise RuntimeError("Original old PR source taxonomy drift")
 if original.digest(args.xlsx,"sha256")!=ORIGINAL_SHA or original.digest(args.xlsx,"md5")!=ORIGINAL_MD5:
  raise RuntimeError("Original Loughborough publisher exact file bytes changed")
 names={x["original_label"]:x for x in tax["taxon_label_source_lookup"]}
 if len(names)!=302:raise RuntimeError("Source original header names changed")
 all_source=original.read_dataset(args.xlsx,names)
 if sum(len(y["known_EPT_fraction"]) for y in all_source.values())!=6491:
  raise RuntimeError("Publisher original 6491 paired source samples changed")
 for country,rows in all_source.items():
  unknown=float(np.mean(rows["unknown_bound_EPT_upper_fraction"]-rows["known_EPT_fraction"]))
  near(unknown,signed["country_profiles"][country]["mean_unknown_taxon_fraction"],
       "Prior source-only unknown-mass gate")
 replays=revalidate_original(all_source,old)
 blocks={}
 for bound,field in (("confirmed_EPT_lower","known_EPT_fraction"),
                     ("all_unknown_as_EPT_upper","unknown_bound_EPT_upper_fraction")):
  blocks[bound]={}
  for domain,d in (("full_original_0to100_visual_fines_domain",all_source),
                   ("common_original_0to50_visual_fines_domain",original.common_domain(all_source))):
   alpha_results={}
   for alpha in ALPHAS if bound=="confirmed_EPT_lower" else (10.,):
    alpha_results[str(int(alpha))]=fit_all(d,field,alpha,domain)
   blocks[bound][domain]=alpha_results
 result={
  "schema":"gesui_external4country_EPT_country_within_between_source_gain_decomposition_v1",
  "original_publisher_Figshare_SHA256":ORIGINAL_SHA,
  "original_hierarchy_taxonomy_sha256":hashlib.sha256(args.source_taxonomy.read_bytes()).hexdigest(),
  "original_source_unknown_gate_sha256":hashlib.sha256(args.source_gate.read_bytes()).hexdigest(),
  "original_prior_countryholdout_sha256":hashlib.sha256(args.original_result.read_bytes()).hexdigest(),
  "all_original_16_signed_countryholdout_models_exactly_replayed":True,
  "original_replay_checks":replays,
  "alpha_grid_frozen":[0,1,10,100],
  "country_equal_training_weights_sum_3":True,
  "visit_equal_sensitivity_weights_sum_3":True,
  "fixed_beta_decomposition_formula":"gain=2beta Cov_T(s,y)+2beta Delta_mean_s Delta_mean_y - beta^2(Var_T(s)+Delta_s^2)",
  "no_sample_site_independence_or_uncertainty_inferred":True,
  "model_result_posthoc":True,
  "no_WWTP_dose_or_epa_HUC8_status_transfer":True,
  "results":blocks}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 quick={}
 for domain in ("full_original_0to100_visual_fines_domain","common_original_0to50_visual_fines_domain"):
  original10=blocks["confirmed_EPT_lower"][domain]["10"]
  quick[domain]={country:{
    "source_n":r["source_sample_n"],
    "train_pooled_beta":r["train_pooled_beta"],
    "train_within_only_beta":r["train_within_country_only_beta"],
    "original_pooled_MSEgain_pct":r["heldout_pooled_country_equal"]["MSE_gain_pct"],
    "train_within_only_MSEgain_pct":r["heldout_within_country_train_only"]["MSE_gain_pct"],
    "train_visit_equal_MSEgain_pct":r["heldout_visit_equal_train_only"]["MSE_gain_pct"],
    "within_target_source_component_pctpoints":r["heldout_pooled_country_equal"]["MSE_gain_components_pctpoints"]["within_target_covariance"],
    "country_mean_offset_component_pctpoints":r["heldout_pooled_country_equal"]["MSE_gain_components_pctpoints"]["mean_domain_offset"],
    "slope_penalty_component_pctpoints":r["heldout_pooled_country_equal"]["MSE_gain_components_pctpoints"]["slope_variance_penalty"],
    "train_source_within_country_cov":r["original_training_source_covariance_decomposition"]["train_within_covariance_EPT_sediment"],
    "train_source_between_country_cov":r["original_training_source_covariance_decomposition"]["train_between_country_mean_covariance"],
    "target_within_country_cov":r["heldout_pooled_country_equal"]["target_within_country_covariance"],
    "within_direction_agrees":r["predeclared_within_country_signal_direction_agrees"],
    "mean_offset_abs_gt_within":r["original_pool_gain_mean_offset_abs_gt_within_abs"]
  } for country,r in original10.items()}
 print("FOURCOUNTRY_EPT_WITHIN_BETWEEN_EXACT_GAINS",json.dumps(quick,sort_keys=True),flush=True)
 print("FOURCOUNTRY_EPT_ALPHA_GRID",json.dumps({domain:{
   a:{name:{
     "original_pooled_gain":r["heldout_pooled_country_equal"]["MSE_gain_pct"],
     "train_within_only_gain":r["heldout_within_country_train_only"]["MSE_gain_pct"]}
     for name,r in countries.items()}
   for a,countries in blocks["confirmed_EPT_lower"][domain].items()}
  for domain in quick},sort_keys=True),flush=True)
if __name__=="__main__":main()
