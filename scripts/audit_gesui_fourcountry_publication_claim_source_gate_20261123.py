#!/usr/bin/env python3
"""Inspectably machine-enforced publication-claim evidence boundary.

Only signed previous ORIGINAL publisher source-design, prior 4country fit and
independent manual WLS QA aggregate JSON. NO raw biological sources or models
opened here; separates descriptive model GO from mechanism/inference STOP.
"""
from __future__ import annotations
import argparse,hashlib,json,pathlib

SOURCE_SHA="a02573bef0e11ea19c9465a78ebecc130fb9c39d76d35a79863c2efd5debcd48"
DOC_SHA="20dd6fc28888e6e28efeb5a9ce65e8e32fc01a3cf55d9fb6251868c556fff1bf"
N={"Australia":3099,"Brazil":232,"New Zealand":342,"UK":2818}
def digest(p):
 return hashlib.sha256(p.read_bytes()).hexdigest()
def require(ok,detail):
 if not ok:raise RuntimeError("Source-locked claim contract inconsistency: "+detail)
def main():
 p=argparse.ArgumentParser()
 for k in ("original_design","original_alpha10","original_decomp","independent_math","original_taxonomy","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 d=json.loads(a.original_design.read_text())
 old=json.loads(a.original_alpha10.read_text())
 audit=json.loads(a.original_decomp.read_text())
 ind=json.loads(a.independent_math.read_text())
 tax=json.loads(a.original_taxonomy.read_text())
 require(d.get("schema")=="gesui_fourcountry_original_source_design_sampling_EXACT_RELEVANT_SIGNALS_v1","Author original design schema")
 require(d["source_original_docx_SHA256"]==DOC_SHA,"Publisher original DOCX SHA")
 require(d["source_original_xlsx_SHA256"]==SOURCE_SHA,"Publisher original XLSX SHA")
 require(audit.get("schema")=="gesui_external4country_EPT_country_within_between_source_gain_decomposition_v1","Original fourcountry numerical decomposition")
 require(ind.get("schema")=="gesui_external_fourcountry_within_between_heldout_manual_weighted_WLS_independent_QA_v1","Independent source WLS QA schema")
 require(old.get("schema")=="gesui_external_original_fourcountry_relative_EPT_country_heldout_DESCRIPTIVE_v1","Original fourcountry country holdout schema")
 require(tax.get("schema")=="gesui_external_fourcountry_EPT_taxon_unknown_mass_source_gate_v1","Original repaired EPT taxonomy source gate")
 require(audit["original_publisher_Figshare_SHA256"]==SOURCE_SHA and old["original_loughborough_figshare_SHA256"]==SOURCE_SHA and ind["publisher_original_XLSX_SHA256"]==SOURCE_SHA,"All fourcountry source original SHA")
 require(ind["original_signed_decomposition_sha256"]==digest(a.original_decomp),"Signed alpha grid verified by independent manual WLS")
 require(ind["all_original_source_country_alpha_domain_model_gains_and_covariance_decompositions_reproduced"] is True,"Independent math has not passed")
 require(audit["all_original_16_signed_countryholdout_models_exactly_replayed"] is True,"Reproduction of original signed alpha10 incomplete")
 require(tax["all_four_country_unknown_taxa_mass_GO"] and tax["country_heldout_descriptive_GO"],"Original EPT source taxon classification STOP")
 require(audit["original_source_unknown_gate_sha256"]==digest(a.original_taxonomy),"Original fixed source taxonomy gate digest mismatch")
 flags=d["author_documented_source_filters"]
 for key in ("UK_environment_agency_sewage_discharge_downstream_exclusion","UK_academic_source_both_spring_and_autumn_sampling",
    "Australia_reference_site_AUS_RIVAS_filter","NZ_multiple_different_academic_studies","NZ_experimental_before_only_filter"):
  require(flags[key]["verified"] is True,"Missing original-author factual sample selection "+key)
 require(d["no_independent_physical_site_ID_in_original_publisher_XLSX"] is True,"Independent physical site ID contradiction")
 require(d["new_zealand_publisher_original_table_n"]==352 and d["new_zealand_XLSX_sample_n"]==342,"NZ methods vs source released count attrition changed")
 nz=d["new_zealand_XLSX_first_token_ID_patterns"]
 require(sum(z["source_observation_rows"] for z in nz.values())==342,"NZ source original samples")
 prefixes=sum(z["source_prefix_first_two_tokens_unique"] for z in nz.values())
 rows_repeated=sum(z["source_rows_sharing_repeated_two_token_prefix"] for z in nz.values())
 require(prefixes==178 and rows_repeated==272,"Corrected string prefix reproducibility")
 full=audit["results"]["confirmed_EPT_lower"]["full_original_0to100_visual_fines_domain"]
 sign={}
 for country,n in N.items():
  a0=full["0"][country]
  a10=full["10"][country]
  original=old["results"]["confirmed_EPT_lower"]["full_original_0to100_visual_fines_domain"][country]
  require(a0["source_sample_n"]==n and a10["source_sample_n"]==n and original["heldout_original_samples"]==n,"Original 4country n")
  g0=a0["heldout_pooled_country_equal"]["MSE_gain_pct"]
  g10=a10["heldout_pooled_country_equal"]["MSE_gain_pct"]
  require(g0>g10>0,"Fixed 4country descriptive gains not positive at alpha0 and10 "+country)
  require(abs(g10-original["heldout_MSE_gain_pct"])<1e-10,"Original fixed alpha10 sign/magnitude changed")
  sign[country]={"publisher_original_samples":n,
    "whole_country_heldout_original_alpha10_gain_pct":g10,
    "whole_country_heldout_posthoc_alpha0_gain_pct":g0,
    "source_specific_model_shrinkage_caution":True}
 country_bias={
    c:r["heldout_pooled_country_equal"]["MSE_gain_components_pctpoints"]
    for c,r in full["10"].items()}
 require(country_bias["Brazil"]["mean_domain_offset"]>country_bias["Brazil"]["within_target_covariance"],"Brazil math source not mean-offset dominated")
 decisions={
  "original_publisher_authoritative_provenance":{
    "verdict":"GO","why":"Original Loughborough 6491 normalized paired sample source and author DOCX/XLSX SHA verified"},
  "country_holdout_conditional_descriptive_MSE":{
    "verdict":"GO_DESCRIPTIVE_ONLY","why":"Four whole-country heldout source responses and model outcomes independently manually replayed; source country weight shrinkage disclosed"},
  "slope_within_country_direction":{
    "verdict":"GO_DESCRIPTIVE_ONLY","why":"Target original within-country covariance observed negative across all four; original paper already established regional EPT sediment responses"},
  "claim_fourcountry_alpha10_small_signal_is_ecologically_weak":{
    "verdict":"NO_GO_MISLEADING","why":"Original training equal country weight totals only 3 while penalizing fine-cover slope alpha10; alpha0 gains 25-34%, model sensitivity"},
  "claim_6491_independent_physical_streams":{
    "verdict":"SOURCE_STOP","why":"Publisher XLSX ANALYSIS_ID is sample join, not proven SITE_ID, UK spring/autumn at sites; NZ code repeats and source mixing"},
  "claim_site_independent_CI_pvalue":{
    "verdict":"SOURCE_STOP","why":"No valid physical site clustering metadata or observation inclusion probabilities"},
  "claim_unselected_national_river_population":{
    "verdict":"SOURCE_STOP","why":"UK EA screens out WWTP-downstream sites, Australia filtered AusRivAS/reference condition, NZ before-only source selected"},
  "claim_sewage_treatment_plant_causal_effect":{
    "verdict":"SOURCE_STOP","why":"No delivered WWTP exposure dose and UK EA source specifically excludes sampled WWTP downstream sites"},
  "claim_EPA_XEMBED_historical_HUC_status_replication":{
    "verdict":"SOURCE_STOP","why":"Fine<2mm visual cover differs from EPA embeddedness protocol, no same older EPA HUC8 historical source membership"},
  "claim_aquatic_to_terrestrial_transfer_or_services":{
    "verdict":"SOURCE_STOP","why":"No quantified subsidy, emergent insect biomass, riparian consumers, host plants, terrestrial butterflies or measured water-to-land flux"},
  "claim_new_global_EPT_fine_sediment_principle":{
    "verdict":"NOT_NOVEL","why":"McKenzie et al GCB 30:e17084 already documented EPT indicators and regional fine sediment response heterogeneity"}
 }
 output={"schema":"gesui_fourcountry_actual_source_DESCRIPTIVE_VS_INFERENCE_publication_claim_gate_v1",
  "original_DOCX_SHA256":DOC_SHA,"original_XLSX_SHA256":SOURCE_SHA,
  "source_design_receipt_sha256":digest(a.original_design),
  "source_taxonomy_gate_receipt_sha256":digest(a.original_taxonomy),
  "original_country_holdout_receipt_sha256":digest(a.original_alpha10),
  "independent_alpha_sensitivity_receipt_sha256":digest(a.original_decomp),
  "independent_math_receipt_sha256":digest(a.independent_math),
  "total_paired_samples":sum(N.values()),"verified_independent_physical_sites":None,
  "NZ_string_prefix_code_groups_not_sites":prefixes,
  "NZ_rows_in_repeated_string_prefix_groups_not_independent_site_estimate":rows_repeated,
  "NZ_method_source_samples":352,"NZ_released_rows":342,
  "UK_EA_WWTP_downstream_source_excluded":True,
  "published_McKenzie_GCB_DOI":"10.1111/gcb.17084",
  "published_Loughborough_data_DOI":"10.17028/rd.lboro.24615570.v1",
  "country_gains":sign,
  "original_alpha10_predictive_MSE_contributions":country_bias,
  "publication_claim_boundaries":decisions,
  "raw_copyrighted_publisher_original_documents_exported":False,
  "no_new_ecological_model_fitted":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("FOURCOUNTRY_PUBLISHABLE_CLAIM_BOUNDARY_QA",json.dumps({
  "all_source_receipts_verified":True,
  "country_source_samples":sum(N.values()),
  "NZ_corrected_prefix_codes":prefixes,
  "NZ_code_recurring_rows":rows_repeated,
  "UK_EA_WWTP_exclusion":True,
  "original_country_gains":sign,
  "verdicts":{k:v["verdict"] for k,v in decisions.items()}},sort_keys=True),flush=True)
if __name__=="__main__":main()
