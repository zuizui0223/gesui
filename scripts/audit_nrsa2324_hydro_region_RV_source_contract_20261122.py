#!/usr/bin/env python3
"""Independent cross-artifact aggregate-only QA for hydrology, region and RV overlap.

No field-level outcome rows. Falsify inflated 'universal width sign reversal'.
"""
import argparse,hashlib,json,math,pathlib
SPECS={"BERW":("MMI_BENT","EPT_fraction"),"BETB":("MMI_BENT","EPT_fraction")}
SUBSETS=("RV_source_only","RV_common_flow","RV_common_flow_and_width")
def load(p):
 return json.loads(pathlib.Path(p).read_text())
def assert_equal(a,b,label):
 if a!=b:raise RuntimeError(label+": disagreement")
def interval_excludes_zero(ci):
 return ci[0]>0 or ci[1]<0
def main():
 p=argparse.ArgumentParser()
 for k in ("hydro","ecoregion","overlap_source","overlap_numeric","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 args=p.parse_args()
 h=load(args.hydro);e=load(args.ecoregion);s=load(args.overlap_source);v=load(args.overlap_numeric)
 assert_equal(h["original_joint_UID"],2139,"national five-table original UID count")
 for key,x in (("e",e["origin_5table_UID"]),("s",s["common_UID_all_original"]),("v",v["same_original_UID"])):
  assert_equal(x,2139,key+" original national UID count")
 for k in ("siteinfo","physicalhabitat"):
  assert_equal(h["original_sha"][k],s["source_SHA"][k],k+" source original SHA")
 for x in (e["EPA_original_SHA"],v["original_EPA_SHA"]):
  assert_equal(x,h["original_sha"],"five original EPA signed table SHAs")
 sourcefile_sha=hashlib.sha256(args.overlap_source.read_bytes()).hexdigest()
 assert_equal(v["source_only_artifact_SHA256"],sourcefile_sha,"RV outcome uses exact prior source cutoffs")
 for k in ("fixed_common_log10_flow_bounds","fixed_common_log10_width_bounds"):
  assert_equal(s[k],v["source_only_preoutcome_RV_cutoffs"][k],"frozen ecological domain")
 qa={"schema":"gesui_NRSA2023_24_width_sign_cross_workflow_independent_QA_v1",
  "original_source_five_table_UIDs":2139,
  "source_sha256":h["original_sha"],
  "primary_hydrological_covariate":"modeled NHDPlusHR mean annual flow, not measured visit discharge",
  "cohort":{},"source_QA_PASSED":True,
  "observational_exploratory_only":True}
 for method, endpoints in SPECS.items():
  qa["cohort"][method]={}
  for response in endpoints:
   z=h["results"][method][response]
   r=e["results"][method][response]
   sub=v["results"][method][response]["fixed_source_overlap_subsets"]
   if z.get("status")!="CORE_NUMERIC_COMPLETE":raise RuntimeError("National hydrology analysis incomplete")
   if not r.get("huc8_previous_ref_reproduced"):raise RuntimeError("Other ecoregion model failed original reference")
   if r["cohort"]["UID"]!=z["source_support"]["all_hydro_complete"]:
    raise RuntimeError("HUC8 and heldout ecoregion cohorts differ")
   if r["cohort"]["HUC8"]!=z["source_support"]["all_hydro_HUC8"]:
    raise RuntimeError("Region source original watershed cohort mismatch")
   region=r["region_extrapolation"]
   if region.get("status")!="COMPLETE_NINE_REGIONS" or region["heldout_ecoregions"]!=9:
    raise RuntimeError("Not all nine external ecological regions were evaluated")
   if not region["zero_shared_HUC8_between_each_fold_train_and_test"]:
    raise RuntimeError("Geographic source leakage")
   prev_n=z["source_support"]["all_hydro_complete"]
   prev_h=z["source_support"]["all_hydro_HUC8"]
   rv_summary={}
   for tag in SUBSETS:
    x=sub[tag]
    original=s["methods"][method]["survey_subsets"][tag]
    assert_equal(x["source_only"],original,method+" "+tag+" exact source-only gate")
    if x["observed_UID"]>original["visits"] or x["independent_HUC8"]>original["HUC8"]:
     raise RuntimeError("Biological cohort larger than its prior source-only opportunity")
    if x["observed_UID"]>prev_n or x["independent_HUC8"]>prev_h:
     raise RuntimeError("RV complete cases exceed national source cohort")
    prev_n=x["observed_UID"];prev_h=x["independent_HUC8"]
    ci=x.get("adjusted_width",{}).get("cluster_HUC8_CI95")
    pred=x.get("HUC8_prediction",{})
    rv_summary[tag]={"UID":x["observed_UID"],"HUC8":x["independent_HUC8"],
      "status":x["status"],
      "width_standardized_beta":x.get("adjusted_width",{}).get("beta_width_ySD_per_1SD_logwidth"),
      "width_beta_cluster95":ci,
      "width_beta_direction_confirmed":interval_excludes_zero(ci) if ci else None,
      "HUC8_OOF_MSE_gain_pct":pred.get("HUC8_equal_MSE_gain_pct"),
      "HUC8_fixedOOF_boot95":pred.get("HUC8_fixedOOF_bootstrap95_pct")}
   national=z["core"]["comparisons"]["width_after_flow_and_area"]
   regional=region["region_bootstrap_fixedOOF_95_pct"]
   qa["cohort"][method][response]={
     "national_UID":r["cohort"]["UID"],"national_HUC8":r["cohort"]["HUC8"],
     "national_HUC8_OOF_width_gain_pct":national["HUC8_equal_MSE_gain_pct"],
     "national_HUC8_fixedOOF_boot95":national["HUC8_fixedOOF_bootstrap95_pct"],
     "nine_ecoregion_equal_MSE_width_gain_pct":region["region_equal_MSE_gain_pct"],
     "nine_ecoregion_boot95":regional,
     "nine_ecoregions_improve":region["n_regions_positive"],
     "national_adjusted_width_beta":r["width_descriptive"].get("standardized_width_beta"),
     "original_RV_common_support_results":rv_summary}
  # No inference about comparing coefficients across incompatible survey protocols.
 # A universal within-flow sign reversal would require non-overlapping opposite signed
 # adjusted effects in each method at the SAME locked RV flow+width support.
 reversal={}
 for endpoint in ("MMI_BENT","EPT_fraction"):
  wa=qa["cohort"]["BERW"][endpoint]["original_RV_common_support_results"]["RV_common_flow_and_width"]
  bo=qa["cohort"]["BETB"][endpoint]["original_RV_common_support_results"]["RV_common_flow_and_width"]
  cwa=wa["width_beta_cluster95"];cbo=bo["width_beta_cluster95"]
  supported=bool(cwa and cbo and
    cwa[0]>0 and cbo[1]<0 and
    wa["HUC8_OOF_MSE_gain_pct"]>0 and bo["HUC8_OOF_MSE_gain_pct"]>0)
  reversal[endpoint]={"opposite_direction_same_source_domain_supported":supported,
    "wadeable_interval":cwa,"boatable_interval":cbo}
 qa["fixed_source_common_domain_opposite_width_slope_claim"]=reversal
 qa["universal_direction_flip_GO"]=all(x["opposite_direction_same_source_domain_supported"] for x in reversal.values())
 qa["eco9_generalization_claim_per_method"]={
  method:all(qa["cohort"][method][endpoint]["nine_ecoregion_boot95"][0]>0
             for endpoint in ("MMI_BENT","EPT_fraction"))
  for method in ("BERW","BETB")}
 qa["source_conclusion"]="National width associations and full-ecoregion transport do NOT identify a universal biological sign reversal; common-flow/width source-domain falsifier must be reported."
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(qa,indent=2,sort_keys=True)+"\n")
 print("WIDTH_ECO9_RV_INDEPENDENT_QA",json.dumps({"source_QA":qa["source_QA_PASSED"],
    "universal_reversal_GO":qa["universal_direction_flip_GO"],
    "regional_transfer":qa["eco9_generalization_claim_per_method"],
    "matched_RV_sign_test":reversal},sort_keys=True),flush=True)
if __name__=="__main__":main()
