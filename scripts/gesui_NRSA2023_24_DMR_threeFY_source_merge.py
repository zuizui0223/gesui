#!/usr/bin/env python3
"""Source-only combine all three EPA fiscal DMR original extracts for 2023/24.

Pre-existing NLDI 195 unique NRSA visits/657 exact outfalls selected before
DMR outcomes. Distinct permitting features, HUC8, states, year and DMR types.
NO absent reporting interpreted as zero; no ecological outcome regressions.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib
import pandas as pd

FILES={2023:"NRSA_POTW_DMR_SOURCE_FY2023.json",2024:"NRSA_POTW_DMR_SOURCE_FY2024.json",2025:"NRSA_POTW_DMR_SOURCE_FY2025.json"}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def main():
 a=argparse.ArgumentParser()
 a.add_argument("--fiscal-dir",type=pathlib.Path,required=True)
 a.add_argument("--source-qa",type=pathlib.Path,required=True)
 a.add_argument("--geo-source",type=pathlib.Path,required=True)
 a.add_argument("--out",type=pathlib.Path,required=True)
 q=a.parse_args()
 source=load(q.source_qa)
 geo=load(q.geo_source)
 if source["original_geocandidate_manifest_SHA256"]!=sha(q.geo_source):
  raise RuntimeError("Historical pre-outcome geographic source SHA differs")
 u=source["eras"]["2023_24"]
 upstream=u["verified_source_only_upstream_candidate_pairs"]
 unique_uids={x["uid"] for x in upstream}
 unique_features={(x["permit"].upper(),x["feature"].upper().zfill(3) if x["feature"].isdigit() and len(x["feature"])<=3 else x["feature"].upper()) for x in upstream}
 if not (len(unique_uids)==195 and len(unique_features)==657):
  raise RuntimeError("Pre-DMR source population 195 UIDs/657 outfalls is not intact")
 orig={x["uid"]:x for x in geo["eras"]["2023_24"]["source_records"]}
 mapping=collections.defaultdict(set)
 site_year={}
 meta={}
 for rec in upstream:
  uid=rec["uid"]
  if uid not in orig:raise RuntimeError("Frozen original source UID no longer matches")
  dt=pd.to_datetime(orig[uid]["date_col"],errors="coerce")
  if pd.isna(dt) or int(dt.year) not in (2023,2024):
   raise RuntimeError("NRSA original site has invalid survey year")
  site_year[uid]=int(dt.year)
  feature=rec["feature"].upper().zfill(3) if rec["feature"].isdigit() and len(rec["feature"])<=3 else rec["feature"].upper()
  mapping[(int(dt.year),rec["permit"].upper(),feature)].add(uid)
  meta[uid]={"year":int(dt.year),"state":rec["state"],"huc8":rec["huc8"]}

 values={}
 archives={}
 official_report_years=collections.Counter()
 code_units=collections.Counter()
 for fy,fname in FILES.items():
  files=list(q.fiscal_dir.rglob(fname))
  if len(files)!=1:raise RuntimeError(f"Missing or duplicate official fiscal original source file for FY{fy}")
  entry=load(files[0])
  if entry["schema"]!="gesui_NRSA2023_24_EXACT_UPSTREAM_POTW_OFFICIAL_DMR_SOURCE_ONLY_FY_v1":
   raise RuntimeError("Original DMR fiscal output schema changed")
  if entry["fiscal_year"]!=fy:raise RuntimeError("Original archive and report fiscal label disagree")
  if entry["source_NLDI_independent_6shard_manifest_SHA256"]!=sha(q.source_qa):
   raise RuntimeError("Different original pre-DMR source census")
  if entry["source_ECHO_historical_geometry_manifest_SHA256"]!=sha(q.geo_source):
   raise RuntimeError("Different original ECHO registered source across fiscal scans")
  archives[fy]={"original_zip_sha256":entry["source_archive_sha256"],"source_counts":entry["source_counts"],
                "qualifying_source_feature_reported_rows":len(entry["exact_outfall_year_param_observed_categories"])}
  for x in entry["year_param_unit_value_type_summary"]:
   code_units[(x["calendar_year"],x["parameter"],x["value_type_code"],x["dmr_unit"])]+=x["rows"]
  for x in entry["exact_outfall_year_param_observed_categories"]:
   k=(x["calendar_year"],x["permit"],x["feature"],x["parameter"],x["class"])
   if k[:3] not in mapping:raise RuntimeError("DMR reports source feature not in original NRSA upstream source+year")
   if k in values:
    raise RuntimeError("Same feature/year/param/class occurred across fiscal archives; independent dedup required")
   values[k]=(int(x["reported_numeric_rows"]),int(x["strictly_positive_numeric_rows"]))
   official_report_years[x["calendar_year"]]+=int(x["reported_numeric_rows"])
 report={}
 for year in (2023,2024):
  cohort={uid for uid,y in site_year.items() if y==year}
  by_class={}
  for category in ("quantity","concentration"):
   by_param={}
   for param in ("total_N","total_P"):
    observed={(p,f) for (y,p,f,z,cl),v in values.items() if y==year and z==param and cl==category and v[0]>0}
    positive={(p,f) for (y,p,f,z,cl),v in values.items() if y==year and z==param and cl==category and v[1]>0}
    site_any={uid for (y,p,f),uids in mapping.items() if y==year and (p,f) in observed for uid in uids}
    site_pos={uid for (y,p,f),uids in mapping.items() if y==year and (p,f) in positive for uid in uids}
    by_param[param]={"reported_source_features":len(observed),"strict_positive_source_features":len(positive),
                     "site_UIDs_with_at_least_one_reported_source":len(site_any),
                     "site_UIDs_with_at_least_one_positive_source":len(site_pos)}
   bN={(p,f) for (y,p,f,z,cl),v in values.items() if y==year and z=="total_N" and cl==category and v[1]>0}
   bP={(p,f) for (y,p,f,z,cl),v in values.items() if y==year and z=="total_P" and cl==category and v[1]>0}
   both=bN&bP
   s2={uid for (y,p,f),uids in mapping.items() if y==year and (p,f) in both for uid in uids}
   by_class[category]={"by_parameter":by_param,
       "same_exact_source_feature_positive_both_N_and_P":len(both),
       "source_associated_NR​​SA_site_UIDs_both_N_P":len(s2),
       "states_with_both_N_P":len({meta[uid]["state"] for uid in s2}),
       "HUC8_with_both_N_P":len({meta[uid]["huc8"] for uid in s2}),
       "meets_min40source_visits_25features_10HUC8_5states":(
            len(s2)>=40 and len(both)>=25 and len({meta[uid]["huc8"] for uid in s2})>=10 and len({meta[uid]["state"] for uid in s2})>=5)}
  report[year]={"source_frozen_UT_site_UIDs":len(cohort),
                "source_frozen_unique_UT_features":len({(p,f) for (y,p,f) in mapping if y==year}),
                "by_report_type":by_class}
 result={
 "schema":"gesui_NRSA2023_24_EPA_FY2023_25_TRUE_ORIGINAL_UPSTREAM_EXACT_DMR_TN_TP_SOURCE_REPORT_GATE_v1",
 "original_source_geo_sha256":sha(q.geo_source),
 "original_source_nldi_6shard_sha256":sha(q.source_qa),
 "all_fiscal_archives":archives,"calendar_source_years":report,
 "observed_DMR_type_unit_profiles":[{"year":k[0],"parameter":k[1],"value_type":k[2],"original_unit":k[3],"reported_rows":v}
                                   for k,v in sorted(code_units.items())],
 "verified_2023_24_upstream_registered_source_UIDs":len(unique_uids),
 "verified_distinct_registered_outfalls":len(unique_features),
 "true_treated_versus_no_treatment_control_sample_NOT_identified":True,
 "true_reported_annual_effluent_mass_NOT_estimated":True,
 "not_actual_NP_delivered_at_NR​​SA_visited_reach":True,
 "no_biodiversity_or_biomass_source_specific_effect_regression":True,
 "unreported_NP_is_NOT_zero":True,
 "prohibition_on_causal_mediation_from_source_only_data":True
 }
 q.out.parent.mkdir(parents=True,exist_ok=True)
 q.out.write_text(json.dumps(result,sort_keys=True,indent=2))
 print("NRSA_DMR2023_24_SOURCE_GATE_RESULTS_START")
 print(json.dumps({"year":report,"archive":archives},sort_keys=True))
 print("NRSA_DMR2023_24_SOURCE_GATE_RESULTS_END")
if __name__=="__main__":main()
