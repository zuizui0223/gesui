#!/usr/bin/env python3
"""Exact scientific identity QA; allow only insignificant floating serialization noise.

Original HUC8 pair digest, sample group membership, regions, all integer gates
and source identity MUST match exactly; floating source summaries may vary at
last binary representation bits, no outcome signs or findings can change.
"""
import argparse,json,math,pathlib
def compare(left,right,root="root",errors=None):
 if errors is None:errors=[]
 if type(left)!=type(right):
  errors.append((root,"TYPE",type(left).__name__,type(right).__name__))
 elif isinstance(left,dict):
  if set(left)!=set(right):
   errors.append((root,"KEYS",sorted(set(left)^set(right))))
  else:
   for k in sorted(left):compare(left[k],right[k],root+"."+str(k),errors)
 elif isinstance(left,list):
  if len(left)!=len(right):errors.append((root,"LENGTH",len(left),len(right)))
  else:
   for i,(x,y) in enumerate(zip(left,right)):compare(x,y,root+"["+str(i)+"]",errors)
 elif isinstance(left,float):
  if not math.isclose(left,right,rel_tol=1e-12,abs_tol=1e-12):
   errors.append((root,"FLOAT",left,right))
 elif left!=right:errors.append((root,"EXACT",left,right))
 return errors
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--original",type=pathlib.Path,required=True)
 p.add_argument("--rerun",type=pathlib.Path,required=True)
 a=p.parse_args()
 x=json.loads(a.original.read_text());y=json.loads(a.rerun.read_text())
 expected="gesui_NRSA2023_original_source_only_HUC8_novelty_env_domain_v1"
 if x.get("schema")!=expected or y.get("schema")!=expected:
  raise RuntimeError("Unexpected signed source-only contract schema")
 for key in ("source_2023_groups","historical_training_HUC8",
   "region_counts","within_ecoregion_source_only_match"):
  if key not in x or key not in y:raise RuntimeError("Missing original source metadata")
 pair_key="source_pair_membership_and_distance_SHA256"
 assert x["within_ecoregion_source_only_match"][pair_key]==y["within_ecoregion_source_only_match"][pair_key]
 assert x["within_ecoregion_source_only_match"]["source_go"] is True
 errs=compare(x,y)
 if errs:raise RuntimeError("Original source identity fails: "+str(errs[:15]))
 print("NRSA_NOVEL_HUC_SOURCE_SNAPSHOT_QA_SUCCESS",
   json.dumps({"HUC8":x["source_2023_HUC8"],
    "fixed_matched_HUC8_pairs":x["within_ecoregion_source_only_match"]["matched_pair_count"],
    "source_pair_SHA256":x["within_ecoregion_source_only_match"][pair_key],
    "floating_source_summary_rtol":1e-12,
    "all_integer_and_string_source_fields_exact":True},sort_keys=True))
if __name__=="__main__":main()
