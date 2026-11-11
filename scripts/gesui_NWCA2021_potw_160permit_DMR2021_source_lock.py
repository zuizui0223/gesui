#!/usr/bin/env python3
"""Source-only lock 160 municipal POTW registered permit hashes near 128 NWCA sites.

Do not read original NWCA vegetation/nutrient responses, DMR monitoring data,
or alter PR #49 original 2026 ECHO outfall source definition. DMR source scope
must be frozen before national FY2021/FY2022 original DMR download/analysis.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

PATH=Path("results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json")
SOURCE_EXPECT="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
def sha(b):return hashlib.sha256(b).hexdigest()

def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
 args=p.parse_args();root=Path(__file__).resolve().parents[1]
 file=root/PATH;original=json.loads(file.read_text())
 if original["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]!=SOURCE_EXPECT:
  raise ValueError("Registered municipal outfall original zip changed")
 if not original["source_precommitted_GEO_eligibility_gate"]["permitted_2021_outfall_distance_association_source_support_PASS"]:
  raise ValueError("Original 568 NWCA site/2026 permit source eligibility failed")
 if not original["no_original_NWCA_H_ALL_FQAI_or_TNTP_outcomes_accessed_in_source_stage"]:
  raise ValueError("Original 2026 NPDES source selected after plant biology")
 records=original["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
 if len(records)!=568 or len({r["NWCA2021_source_site_visit_SHA256"] for r in records})!=568:
  raise ValueError("Original NWCA source 568 site IDs mutated")
 near=[r for r in records if r["outfall_rows_within_10km"]>0]
 if len(near)!=128 or len({r["original_state"] for r in near})!=32:
  raise ValueError("Original municipal ECHO strict near site 128/32 changed")
 permits=sorted({perm for r in near for perm in r["independent_10km_NPDES_permits_SHA256"]})
 if len(permits)!=160 or any(len(k)!=64 for k in permits):
  raise ValueError("Original near10 matched permit population changed from 160")
 groups=defaultdict(list)
 for r in near:
  groups[r["original_state"]].append(r["NWCA2021_source_site_visit_SHA256"])
 record={
  "schema":"gesui_US_NWCA2021_original128near_160municipalNPDESpermitted_2021DMR_PREOUTCOME_source_v01",
  "UTC":datetime.now(timezone.utc).isoformat(),
  "original_EPA_NWCA2021_568source_2026_NPDES_geo_SHA256":sha(file.read_bytes()),
  "original_2026_ECHO_NPDES_source_zip_SHA256":SOURCE_EXPECT,
  "original_568_wetland_survey_firstvisit_source_sites":568,
  "original_2021_study_potential_near10_sites":128,
  "original_2021_study_potential_near10_states":32,
  "original_registered_municipal_NPDES_permit_ids_SHA256_COUNT":len(permits),
  "exact_source_frozen_near10_registered_municipal_permit_ID_SHA256":permits,
  "original_site_near10_permits_2021outfallPlausible_from_2026_registry":[
   {"NWCA_site_visit_SHA256":r["NWCA2021_source_site_visit_SHA256"],
    "state":r["original_state"],"ecoregion":r["original_ecoregion"],
    "permitted_sewage_NPDES_ID_SHA256":r["independent_10km_NPDES_permits_SHA256"],
    "registered_potential_outfall_rows_10km":r["outfall_rows_within_10km"]}
   for r in near],
  "original_source_study_2026_registered_permitted_feature_number_not_preserved_in_previous_geo_json":True,
  "2026_weekly_permit_date_plausible_is_NOT_proof_2021_actual_discharge":True,
  "future_FY2021_AND_FY2022_source_DMR_fiscal_year_required_for_calendar2021":True,
  "DMR_calendar2021_parameters":{"00600":"total nitrogen as N","00665":"total phosphorus as P","50050":"effluent flow, separate optional"},
  "precommitted_DMR_monitored_locations":["1","2","EG"],
  "precommitted_DMR_external_outfall_type":"EXO",
  "precommitted_year_window":["2021-01-01","2021-12-31"],
  "NO_butterfly_or_wetland_plant_HALL_or_measured_2021_wetland_TNTP_accessed":True,
  "NO_FY2021_FY2022_DMR_measurements_accessed_in_this_source_stage":True
 }
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(record,sort_keys=True,indent=2)+"\n")
 print(json.dumps({"near_sites":len(near),"distinct_permits":len(permits),"states":len(groups),"SHA":record["original_EPA_NWCA2021_568source_2026_NPDES_geo_SHA256"]},indent=2))
if __name__=="__main__":main()
