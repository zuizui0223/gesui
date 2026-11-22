#!/usr/bin/env python3
"""Fixed descriptive original 4-country fine-sediment vs EPT share OOD test.

Only execute if independently verified source taxonomy quality all 4 countries
passes original strict UNKNOWN<=2% and >5% UNKNOWN fraction in <=10% samples.
No site IDs, no independent sample-level p/CI, not EPA XEMBED/source history.
"""
from __future__ import annotations
import argparse,hashlib,json,math,pathlib
import numpy as np
import gesui_external_fourcountry_xlsx_readonly_ooxml as oo
from gesui_external_fourcountry_source_grain_20261122 import SHEETS,ORIGINAL_SHA,ORIGINAL_MD5,digest

ALPHA=10.
def read_dataset(path,taxon_map):
 out={}
 for country,(tx,sed,expected) in SHEETS.items():
  head,rows=oo.read_all(path,tx);header,env=oo.read_all(path,sed)
  if len(rows)!=expected or len(env)!=expected:raise RuntimeError("Publisher original source sample count changed")
  taxa=head[1:]
  if any(t not in taxon_map for t in taxa):raise RuntimeError("Publisher taxon source names not classified")
  ii_EPT=[i for i,t in enumerate(taxa) if taxon_map[t]["state"]=="EPT"]
  ii_UNK=[i for i,t in enumerate(taxa) if taxon_map[t]["state"]=="UNKNOWN"]
  original_fines={oo.normalize_id(row[0]):oo.as_float(row[1]) for row in env}
  if len(original_fines)!=expected or len(set(oo.normalize_id(row[0]) for row in rows))!=expected:
   raise RuntimeError("Original ANALYSIS_ID duplicate or missing sample")
  y=[];yupper=[];s=[];missing=0
  for row in rows:
   vals=np.array([oo.as_float(q) for q in row[1:]],dtype=float)
   if len(vals)!=len(taxa) or not np.isfinite(vals).all() or np.any(vals<0):
    raise RuntimeError("Original normalized publisher abundance row invalid")
   if abs(float(vals.sum())-1)>1e-9:raise RuntimeError("Publisher original row not normalized")
   t=oo.normalize_id(row[0])
   if t not in original_fines:raise RuntimeError("Unpaired original sediment ANALYSIS_ID")
   fines=float(original_fines[t])
   if not 0<=fines<=100:raise RuntimeError("Publisher visual fines outside 0–100")
   y.append(float(vals[ii_EPT].sum()))
   yupper.append(float(vals[ii_EPT].sum()+vals[ii_UNK].sum()))
   s.append(fines/100)
  out[country]={"sediment_fraction":np.array(s),
    "known_EPT_fraction":np.array(y),
    "unknown_bound_EPT_upper_fraction":np.array(yupper)}
 return out
def fit_three(train,yfield):
 s=np.concatenate([train[c]["sediment_fraction"] for c in train])
 y=np.concatenate([train[c][yfield] for c in train])
 w=np.concatenate([np.full(len(train[c][yfield]),1/len(train[c][yfield])) for c in train])
 if not (np.isfinite(s).all() and np.isfinite(y).all() and np.isfinite(w).all()):
  raise RuntimeError("Invalid original measurement source")
 if not np.isclose(w.sum(),len(train)):raise RuntimeError("Unbalanced train country normalization")
 mus=float(np.sum(w*s)/np.sum(w))
 muy=float(np.sum(w*y)/np.sum(w))
 numerator=float(np.sum(w*(s-mus)*(y-muy)))
 denominator=float(np.sum(w*(s-mus)**2))
 coef=numerator/(denominator+ALPHA)
 return {"train_3country_equal_country_baseline_mean_EPT_fraction":muy,
   "train_3country_equal_country_mean_fines_fraction":mus,
   "Ridge10_slope_EPT_fraction_per_fines_0to1":coef,
   "Ridge10_positive_nonnegative_penalty":True}
def predict(target,model):
 s=target["sediment_fraction"]
 mean=model["train_3country_equal_country_baseline_mean_EPT_fraction"]
 ms=model["train_3country_equal_country_mean_fines_fraction"]
 beta=model["Ridge10_slope_EPT_fraction_per_fines_0to1"]
 return np.full(len(s),mean),np.asarray(mean+beta*(s-ms),float)
def compare(data,yfield):
 out={}
 for country in sorted(data):
  train={c:data[c] for c in sorted(data) if c!=country}
  test=data[country]
  if set(train)=={country} or len(train)!=3:raise RuntimeError("Country holdout group wrong")
  if len(test[yfield])<1 or any(len(train[k][yfield])<1 for k in train):
   raise RuntimeError("Fixed original 0–50% visual fine cover support absent in heldout/train country")
  model=fit_three(train,yfield)
  p0,p1=predict(test,model)
  y=test[yfield]
  e0=float(np.mean((y-p0)**2));e1=float(np.mean((y-p1)**2))
  if not(e0>0 and np.isfinite(e1)):
   raise RuntimeError("Invalid heldout source domain baseline MSE")
  out[country]={"heldout_country":country,"heldout_original_samples":len(y),
    "train_countries_original_samples":{c:len(train[c][yfield]) for c in train},
    "observed_target_mean_relative_EPT":float(y.mean()),
    "observed_target_mean_visual_fines_fraction":float(test["sediment_fraction"].mean()),
    "model":model,"heldout_MSE_train_country_mean_baseline":e0,
    "heldout_MSE_add_linear_ridge10_fines":e1,
    "heldout_MSE_gain_pct":float(100*(e0-e1)/e0)}
 return out
def common_domain(original):
 return {c:{k:v[mask] for k,v in d.items()}
   for c,d in original.items() for mask in [d["sediment_fraction"]<=.5]}
def main():
 p=argparse.ArgumentParser()
 for k in ("xlsx","source_gate","source_taxonomy","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 signed=json.loads(a.source_gate.read_text())
 tax=json.loads(a.source_taxonomy.read_text())
 if signed.get("schema")!="gesui_external_fourcountry_EPT_taxon_unknown_mass_source_gate_v1":
  raise RuntimeError("Independent source unknown-mass gate missing")
 if signed["country_heldout_descriptive_GO"] is not True or signed["all_four_country_unknown_taxa_mass_GO"] is not True:
  output={"schema":"gesui_4country_fixed_relative_EPT_crosscountry_GO_STOP",
    "status":"SOURCE_STOP_UNKNOWN_TAXON_MASS","original_threshold_NOT_RELAXED":True,
    "numeric_EPT_sediment_models_fitted":False}
  a.out.parent.mkdir(parents=True,exist_ok=True)
  a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
  print("FOURCOUNTRY_SOURCE_STOP_NOT_FITTED",json.dumps(output,sort_keys=True))
  return
 if signed["hierarchy_amendment_to_original_GBIF_taxonomy"] is not True:
  raise RuntimeError("Original source was not authority-repaired before run")
 if signed["publisher_original_XLSX_SHA256"]!=ORIGINAL_SHA or tax["Figshare_original_XLSX_SHA256"]!=ORIGINAL_SHA:
  raise RuntimeError("Original publisher SHA disagreement")
 if hashlib.sha256(a.source_taxonomy.read_bytes()).hexdigest()!=signed["taxonomy_original_source_json_sha256"]:
  raise RuntimeError("Source gate taxonomy content digest changed")
 if tax["schema"]!="gesui_original_McKenzie_4country_taxonomy_hierarchy_source_amendment_v2":
  raise RuntimeError("Source taxonomy authority v2 not used")
 if digest(a.xlsx,"sha256")!=ORIGINAL_SHA or digest(a.xlsx,"md5")!=ORIGINAL_MD5:
  raise RuntimeError("Downloaded original publisher XLSX SHA/MD5 mismatch")
 mapping={x["original_label"]:x for x in tax["taxon_label_source_lookup"]}
 if len(mapping)!=302:raise RuntimeError("Original 302 taxa count changed")
 d=read_dataset(a.xlsx,mapping)
 for country,z in d.items():
  unknown=np.mean(z["unknown_bound_EPT_upper_fraction"]-z["known_EPT_fraction"])
  allowed=signed["country_profiles"][country]["mean_unknown_taxon_fraction"]
  if not np.isclose(unknown,allowed,rtol=0,atol=1e-12):
   raise RuntimeError("Original source gate unknown fraction failed numerical reinspection")
  if len(z["known_EPT_fraction"])!=SHEETS[country][2]:
   raise RuntimeError("Country source count differs")
 results={}
 for spec,field in (("confirmed_EPT_lower","known_EPT_fraction"),
                    ("all_unknown_as_EPT_upper","unknown_bound_EPT_upper_fraction")):
  results[spec]={"full_original_0to100_visual_fines_domain":compare(d,field),
   "common_original_0to50_visual_fines_domain":compare(common_domain(d),field)}
 out={"schema":"gesui_external_original_fourcountry_relative_EPT_country_heldout_DESCRIPTIVE_v1",
  "original_loughborough_figshare_SHA256":ORIGINAL_SHA,
  "source_hierarchy_taxonomy_sha256":hashlib.sha256(a.source_taxonomy.read_bytes()).hexdigest(),
  "source_unknown_gate_digest":hashlib.sha256(a.source_gate.read_bytes()).hexdigest(),
  "source_unknown_gate_GO":True,
  "original_source_total_samples":sum(len(v["known_EPT_fraction"]) for v in d.values()),
  "original_countries":list(sorted(d)),
  "fine_sediment_concept":"publisher visual percent channel-bed covered by particles <2mm, NOT EPA XEMBED",
  "fitted_already_published_ecological_question":True,
  "no_independent_site_ID_available":True,
  "no_independent_site_based_CI_or_pvalues":True,
  "not_EPA_old_source_history_interaction_replication":True,
  "not_sewage_exposure_or_land_aquatic_food_flux":True,
  "model_config":{"train_country_weight":1,"target_country_holdout":True,
    "Ridge_alpha":ALPHA,"linear_predictor":"visual fine cover fraction 0-1",
    "intercept_unpenalized":True},
  "results":results}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
 print("FOURCOUNTRY_ORIGINAL_DESCRIPTIVE_COUNTRY_HELDOUT",json.dumps(
  {spec:{domain:{country:{
    "n":res["heldout_original_samples"],"gain_MSE_pct":res["heldout_MSE_gain_pct"],
    "train_beta":res["model"]["Ridge10_slope_EPT_fraction_per_fines_0to1"],
    "mean_target_EPT":res["observed_target_mean_relative_EPT"]}
    for country,res in blocks.items()} for domain,blocks in domains.items()}
   for spec,domains in results.items()},sort_keys=True))
if __name__=="__main__":main()
