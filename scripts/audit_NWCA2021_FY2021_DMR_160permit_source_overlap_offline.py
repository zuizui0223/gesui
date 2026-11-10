#!/usr/bin/env python3
"""Offline independent gate arithmetic for original 160 FY2021 NPDES permits.

Original national EPA FY2021 raw archive SHA was source-verified by earlier
~576MB streaming workflow. This independent replay checks only ORIGINAL
pre-frozen 128 wetland sites, registered NPDES admin-permit hashes and nutrient
monitoring classification. No H_ALL, NTL/PTL wetland outcomes, no inference.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
GEOM=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
DMR=D/"EPA_NWCA2021_FY2021_REAL_NPDES_DMR_160MUNICIPAL_PERMIT_TNTP_MONITOR_SOURCE_GATE.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def audit(root):
    o=json.loads((root/GEOM).read_text())
    v=json.loads((root/DMR).read_text())
    if v["original_2026_NPDES_registered_outfall_568_NWCA_source_SHA256"]!=sha(root/GEOM):
        raise ValueError("Original frozen nearNPDES WWTP location source not matched to DMR")
    if not v["this_stage_NEVER_accessed_2021_wetland_Shannon_FQAI_or_source_TNTP_values"]:
        raise ValueError("Source DMR stage accessed posthoc biological outcomes")
    if o["source_precommitted_GEO_eligibility_gate"]["permitted_2021_outfall_distance_association_source_support_PASS"] is not True:
        raise ValueError("Original source near/far qualification changed")
    near=[x for x in o["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
          if x["outfall_rows_within_10km"]>0]
    permits={z for site in near for z in site["independent_10km_NPDES_permits_SHA256"]}
    if len(near)!=128 or len(permits)!=160:
        raise ValueError("Original firstvisit wetland 128 and administrative WWTP permit IDs 160 changed")
    original=v["every_160_original_nearby_candidate_administrative_sewage_permit_DMR_availability_QA"]
    if set(original)!=permits:raise ValueError("Source DMR original 160 admin permit status incomplete")
    def subset(key):
        return {pid for pid,r in original.items() if r.get(key,0)>0}
    anynum=subset("original_DMR_ANY_parameter_numeric_n")
    total_n=subset("NUMERIC_TN_00600")
    total_p=subset("NUMERIC_TP_00665")
    combos={
      "any_numeric_DMR_FY2021":anynum,
      "total_N_00600_numeric_FY2021":total_n,
      "total_P_00665_numeric_FY2021":total_p,
      "total_N_OR_P_numeric_FY2021":total_n|total_p,
      "total_N_AND_P_numeric_FY2021":total_n&total_p
    }
    stats={}
    for name,pids in combos.items():
        selected=[r for r in near if pids.intersection(r["independent_10km_NPDES_permits_SHA256"])]
        counts={
          "original_NWCA_568_source_near10_sites_with_geographically_associated_PERMIT_reporting_n":len(selected),
          "original_reg_sewage_NPDES_permits_with_FY2021_monitoring_data":len(pids),
          "original_distinct_US_states_with_monitor_and_near_site":len({r["original_state"] for r in selected}),
          "original_ecoregions_with_monitor_and_near_site":len({r["original_ecoregion"] for r in selected})
        }
        original_result=v["original_160_registered_sewage_source_permit_count_summary"][name]
        for k,value in counts.items():
            if value!=original_result[k]:
                raise ValueError("Original official DMR source overlap counter not reproduced: "+name+"."+k)
        if sorted(r["NWCA2021_source_site_visit_SHA256"] for r in selected)!=sorted(original_result["original_near_source_2021_wetland_site_hashes_only"]):
            raise ValueError("NPDES-to-original-wetland site SHA crosswalk not faithfully reproduced")
        stats[name]=counts
    a=stats["total_N_OR_P_numeric_FY2021"]
    gate=v["original_160_registered_near_outfall_sewage_permits_FY2021_anyNutrient_model_support"]
    passed=(a["original_NWCA_568_source_near10_sites_with_geographically_associated_PERMIT_reporting_n"]>=40 and
            a["original_reg_sewage_NPDES_permits_with_FY2021_monitoring_data"]>=15 and
            a["original_distinct_US_states_with_monitor_and_near_site"]>=8)
    if passed!=gate["2021_FY_Nutrient_MONITORING_SOURCE_GATE_PASS"]:
        raise ValueError("Source-only objective DMR evidence cutoff changed")
    if v["FY2021_national_original_DMR_source_row_count_QA"]["original_national_DMR_csv_rows"]<100_000:
        raise ValueError("Incomplete original USA EPA FY2021 national DMR source")
    return {
      "schema":"gesui_EPA_NWCA2021_128source_wetlands_160permitted_WWTP_FY2021_DMR_offline_reproduction_v01",
      "original_568_NWCA2021_ECHO2026_WWTP_GEO_source_SHA256":sha(root/GEOM),
      "original_FY2021_ECHO_national_DMR_rawStream_result_SHA256":sha(root/DMR),
      "source_national_FY2021_EPA_DMR_ZIP_SHA256":v["original_EPA_fiscalyear2021_national_NPDES_DMR_ZIP_SHA256"],
      "all_original_160_WWTP_admin_source_permits_present":True,
      "all_original_128_near_WWTP_wetland_source_sites_reproduced":True,
      "all_5_numeric_DMR_NP_source_monitoring_groups_recomputed_exactly":True,
      "original_source_national_DMR_fiscal2021_row_count":v["FY2021_national_original_DMR_source_row_count_QA"]["original_national_DMR_csv_rows"],
      "independently_reconstructed_source_perm_to_site_monitoring":stats,
      "precommitted_source_only_monitoring_gate":passed,
      "fiscal2021_measured_any_under_same_administrative_permit_is_not_same_nearby_outfall_or_receiving_wetland":True,
      "2021_NWCA_plant_and_water_values_not_opened_in_DMR_work":True,
      "NO_SEWAGE_CAUSAL_ECOLOGICAL_EFFECT":True
    }

def main():
    import datetime
    a=argparse.ArgumentParser()
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=audit(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(o,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in o.items() if k in (
        "source_national_FY2021_EPA_DMR_ZIP_SHA256",
        "original_source_national_DMR_fiscal2021_row_count",
        "independently_reconstructed_source_perm_to_site_monitoring",
        "precommitted_source_only_monitoring_gate")},indent=2))

if __name__=="__main__":
    main()
