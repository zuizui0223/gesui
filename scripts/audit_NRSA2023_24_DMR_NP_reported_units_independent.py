#!/usr/bin/env python3
"""Independent source-only historical EPA DMR numerical TYPE/UNIT quality audit.

The true year-by-year observed N/P presence is NOT converted into mass or
sewage-related biodiversity coefficients. This script contains no nutrient
measurement magnitudes, biological outcomes, or derived per-site doses.
"""
import argparse,collections,json,pathlib
p=argparse.ArgumentParser()
p.add_argument("--input",type=pathlib.Path,required=True)
p.add_argument("--output",type=pathlib.Path,required=True)
a=p.parse_args()
x=json.loads(a.input.read_text())
if x["schema"]!="gesui_NRSA2023_24_EPA_FY2023_25_TRUE_ORIGINAL_UPSTREAM_EXACT_DMR_TN_TP_SOURCE_REPORT_GATE_v1":
 raise RuntimeError("Unexpected original three-FY report schema")
if sorted(map(int,x["all_fiscal_archives"]))!=[2023,2024,2025]:
 raise RuntimeError("Original source archives incomplete")
profile=collections.defaultdict(collections.Counter)
for r in x["observed_DMR_type_unit_profiles"]:
 yy=r["year"];nm=r["parameter"]
 v=str(r["value_type"]).strip().upper();u=str(r["original_unit"]).strip()
 cat="Q_quantity" if v.startswith("Q") else ("C_concentration" if v.startswith("C") else "OTHER")
 profile[(yy,nm,cat)][(v,u)]+=r["reported_rows"]
detail=[]
for (year,param,typ),distribution in sorted(profile.items()):
 detail.append({"year":year,"param":param,"class":typ,"n_original_reporting_rows":sum(distribution.values()),
                "original_unit_types":[{"dmr_value_type":v,"native_unit":u,"rows":cnt}
                                       for (v,u),cnt in distribution.most_common()]})
yearly={}
for year in ("2023","2024"):
 z=x["calendar_source_years"][year]
 yearly[year]={"pre_source_UT_site_visits":z["source_frozen_UT_site_UIDs"]}
 for typ in ("quantity","concentration"):
  c=z["by_report_type"][typ]
  yearly[year][typ]={"at_least_one_upstream_reported_NP_same_EXO_sites":c["source_associated_NRSA_site_UIDs_both_N_P"],
                     "positive_both_NP_exact_outfall_features":c["same_exact_source_feature_positive_both_N_and_P"],
                     "original_gate_met":c["meets_min40source_visits_25features_10HUC8_5states"]}
  if c["meets_min40source_visits_25features_10HUC8_5states"]:
   raise RuntimeError("Original gate changed; do not retain STOP conclusion")
out={"schema":"NRSA2023_24_DMR_ORIGINAL_NP_UNITS_AND_NONIDENTIFIABILITY_QA",
     "all_original_fiscal_archives_present":True,"original_positive_reporting_units":detail,"cohort":yearly,
     "positive_DMR_quantity_or_concentration_is_NOT_assumed_true_annual_mass_or_delivery":True,
     "unreported_DMR_not_environmental_zero":True,
     "original_predeclared_both_NP_gate_fails_both_years":True,
     "no_biodiversity_biomass_ecosystem_service_effects_estimated":True}
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(out,sort_keys=True,indent=2))
print("INDEPENDENT_DMR_SOURCE_TYPES_SUCCESS",json.dumps({"years":yearly,"categories":detail},sort_keys=True))
