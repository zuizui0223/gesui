#!/usr/bin/env python3
"""Join exact FY2021+FY2022 original EPA DMR permit reporting support in CY2021.

No NWCA H_ALL/FQAI/2021 wetland TN/TP response read. Permit-only nutrient
monitoring is NOT physical outfall-to-wetland exposure. Exact feature
geography/hydrology must be separately verified after this stage.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
PERM=D/"EPA_NWCA2021_160_SEWAGE_PERMITS_BEFORE_DMR2021_SOURCE.json"
GEO=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(root,first,second):
    permit=json.loads((root/PERM).read_text())
    geo=json.loads((root/GEO).read_text())
    docs=[json.loads(x.read_text()) for x in (first,second)]
    if sorted(x["source_fiscal_year"] for x in docs)!=[2021,2022]:
        raise ValueError("Both fiscal DMR years necessary for calendar 2021")
    if len({x["source_original_160permit_preDMR_manifest_SHA256"] for x in docs})!=1 or docs[0]["source_original_160permit_preDMR_manifest_SHA256"]!=sha(root/PERM):
        raise ValueError("Original 160 source-manifest changed between DMR datasets")
    if permit["original_EPA_NWCA2021_568source_2026_NPDES_geo_SHA256"]!=sha(root/GEO):
        raise ValueError("Pre-outcome EPA 2026 NPDES geographical source changed")
    ids=set(permit["exact_source_frozen_near10_registered_municipal_permit_ID_SHA256"])
    if len(ids)!=160:raise ValueError("Original national 160 near-source permits not frozen")
    # Exact original row hash can repeat in consecutive fiscal-year releases.
    seen=set();unique=[];duplicate=Counter()
    for doc in docs:
        for r in doc["original_2021_full_source_selected_DMR_rows_no_raw_reported_values"]:
            k=(r["source_permit_ID_SHA256"],r["source_permit_feature_ID_SHA256"],
               r["original_2021_DMR_VALUE_ID_SHA256"],r["original_monitoring_period_end_date"],
               r["original_parameter_code"],r["original_DMR_VALUE_TYPE_CODE"])
            if k in seen:
                duplicate["same_DMR_key_seen_twice"]+=1
                continue
            seen.add(k)
            if r["source_permit_ID_SHA256"] not in ids:
                raise ValueError("ECHO 2021 DMR result introduced unselected NPDES permits")
            unique.append(r)
    group=defaultdict(lambda:{"N":set(),"P":set(),"flow":set(),"outfalls":set(),"NP_pair":set()})
    both_feature=defaultdict(lambda:defaultdict(set))
    for r in unique:
        if not r["strict_external_outfall_effluent_numeric_uncensored"]:continue
        code=r["original_parameter_code"]
        ty={"00600":"N","00665":"P","50050":"flow"}[code]
        permit_id=r["source_permit_ID_SHA256"]
        feature=r["source_permit_feature_ID_SHA256"]
        group[permit_id][ty].add((feature,r["original_monitoring_period_end_date"]))
        group[permit_id]["outfalls"].add(feature)
        both_feature[(permit_id,feature)][ty].add(r["original_monitoring_period_end_date"])
    summaries={}
    hasN={k for k,v in group.items() if v["N"]}
    hasP={k for k,v in group.items() if v["P"]}
    hasFlow={k for k,v in group.items() if v["flow"]}
    bothPermit=hasN&hasP
    exact_same_featureNP={k for k,v in both_feature.items() if v["N"] and v["P"]}
    same_feature_same_periodNP={k for k,v in both_feature.items() if v["N"]&v["P"]}
    records=permit["original_site_near10_permits_2021outfallPlausible_from_2026_registry"]
    sites={x["NWCA_site_visit_SHA256"]:x for x in records}
    if len(sites)!=128:raise ValueError("Original 128 source-near NWCA sites changed")
    geo_sites=geo["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    near=[x for x in geo_sites if x["outfall_rows_within_10km"]>0]
    far=[x for x in geo_sites if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or x["nearest_outfall_geodesic_km_within25km_else_null"]>=20]
    strata={x["original_state"]+"|"+x["original_ecoregion"] for x in near}&{
           x["original_state"]+"|"+x["original_ecoregion"] for x in far}
    matched_near={x["NWCA2021_source_site_visit_SHA256"] for x in near if x["original_state"]+"|"+x["original_ecoregion"] in strata}
    if len(matched_near)!=103:
        raise ValueError("Original 240 site source geography changed before DMR")
    def sitecount(permset):
        total=[s for s in records if set(s["permitted_sewage_NPDES_ID_SHA256"])&permset]
        return {
            "all_original_128_near_source_sites_having_any_reported_permit":len(total),
            "source_matched_state_ecoregion_103_near_sites_having_any_reported_permit":
                sum(s["NWCA_site_visit_SHA256"] in matched_near for s in total),
            "distinct_original_state_codes":len({s["state"] for s in total}),
        }
    out={
       "schema":"gesui_NWCA2021_160_originalPOTW_permit_FY2021_FY2022_calendar2021_DMR_nutrient_SOURCE_FEASIBILITY_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "source_preDMR_160_permits_SHA256":sha(root/PERM),
       "source_geo_568_NWCA_registered_outfall_SHA256":sha(root/GEO),
       "source_FY2021_DMR_fiscal_JSON_SHA256":sha(first) if docs[0]["source_fiscal_year"]==2021 else sha(second),
       "source_FY2022_DMR_fiscal_JSON_SHA256":sha(second) if docs[1]["source_fiscal_year"]==2022 else sha(first),
       "source_original_EPA_national_DMRS_ZIP_SHA256":{str(x["source_fiscal_year"]):x["source_original_FY_DMR_ZIP_SHA256"] for x in docs},
       "source_original_wetlands_with_near10_source_n":128,
       "source_2021_municipal_NPDES_near_permit_population_n":160,
       "source_original_state_eco_plant_comparison_near_sites_n":103,
       "source_2021_TN_00600_and_TP_00665_DMR_count_rows_after_year_FY_dedup":len(unique),
       "source_original_Monthly_FY_duplicate_rows":dict(duplicate),
       "source_permit_level_numeric_final_effluent_counts":{
          "160_original_near_permits_with_numeric_TN00600":len(hasN),
          "160_original_near_permits_with_numeric_TP00665":len(hasP),
          "160_original_near_permits_with_both_NP":len(bothPermit),
          "160_original_near_permits_with_numeric_flow50050":len(hasFlow),
          "permit_feature_pairs_with_any_numeric_TN":sum(bool(v["N"]) for v in both_feature.values()),
          "permit_feature_pairs_with_any_numeric_TP":sum(bool(v["P"]) for v in both_feature.values()),
          "permit_feature_pairs_with_numeric_TN_and_TP":len(exact_same_featureNP),
          "permit_feature_pairs_with_TNandTP_same_monitoring_period_end_date":len(same_feature_same_periodNP)
       },
       "source_wetland_candidates_with_nearby_PERMIT_reporting":{
          "TN_permit_only":sitecount(hasN),
          "TP_permit_only":sitecount(hasP),
          "bothNP_permit_only":sitecount(bothPermit),
          "anyNP_permit_only":sitecount(hasN|hasP),
          "any_reported_Flow_permit_only":sitecount(hasFlow)
       },
       "source_feature_specific_DMR_reported_permits_and_months_no_numeric_values":[
          {"NPDES_permit_SHA256":k,
           "numeric_TN_2021_feature_date_n":len(v["N"]),
           "numeric_TP_2021_feature_date_n":len(v["P"]),
           "numeric_flow_2021_feature_date_n":len(v["flow"]),
           "same_exact_DMR_outfall_feature_n":len(v["outfalls"])}
          for k,v in sorted(group.items())],
       "DMR_permit_match_does_not_verify_2026_nearby_OUTFALL_NUMBER_or_2021_WWTP_to_wetland_flow":True,
       "future_exact_2021_DMR_outfall_feature_GEOMETRY_and_NHDPLUS_requires_separate_stage":True,
       "actual_2021_nutrient_loads_not_calculated_without_common_unit_flow_duration_and_hydrology":True,
       "source_no_plant_diversity_HALL_FQAI_or_wetland_TNTP_values_read":True
    }
    return out

def main():
 p=argparse.ArgumentParser()
 p.add_argument("--fy2021",type=Path,required=True);p.add_argument("--fy2022",type=Path,required=True)
 p.add_argument("--out",type=Path,required=True)
 a=p.parse_args();root=Path(__file__).resolve().parents[1]
 output=run(root,a.fy2021,a.fy2022)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
 print(json.dumps({"permitted_SOURCE_near":output["source_2021_municipal_NPDES_near_permit_population_n"],
      "measured_effluent":output["source_permit_level_numeric_final_effluent_counts"],
      "site_examples_NO_BIOLOGY":output["source_wetland_candidates_with_nearby_PERMIT_reporting"]},indent=2))
if __name__=="__main__":main()
