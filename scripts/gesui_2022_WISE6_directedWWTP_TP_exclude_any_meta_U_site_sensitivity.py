#!/usr/bin/env python3
"""POSTHOC QUALITY sensitivity: remove entire source WISE6 2022 TP sites that
contain ANY original metadata_observationStatus='U', then refit the exact
171-site original hydrologically directed WWTP source vs downstream placebo
2degree blocked CV.

A prior EEA2026 full 102,569,057-record source status audit already fixed the
9 original U rows on 9 DE stations; no new raw WQ or butterfly data accessed.
The original model and every source/TP sample remain frozen unchanged.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import gesui_2022_WISE6_realTP_directed_WWTP_source_vs_downstream_CV as original

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-status",type=Path,required=True)
    p.add_argument("--original-directed-source",type=Path,required=True)
    p.add_argument("--original-realTP-numeric",type=Path,required=True)
    p.add_argument("--original-WorldCover",type=Path,required=True)
    p.add_argument("--original-171site-CV",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    status=json.loads(a.original_status.read_text())
    source=json.loads(a.original_directed_source.read_text())
    tp=json.loads(a.original_realTP_numeric.read_text())
    wc=json.loads(a.original_WorldCover.read_text())
    first=json.loads(a.original_171site_CV.read_text())
    if (not source["no_numeric_TP_opened_source_gate_pass"] or
        not first["source_and_numeric_support_regression_gate_pass"]):
        raise ValueError("Original 2022 source-only/numeric TP quality gate not passed")
    if first["source_only_directed_WWTP_source_gate_SHA256"]!=sha(a.original_directed_source):
        raise ValueError("Original source SHA changed after successful initial CV")
    if first["2022_WISE6_original_numeric_TP_site_year_source_SHA256"]!=sha(a.original_realTP_numeric):
        raise ValueError("Original 2022 TP numeric source mutated")
    if status["original_2022_TP_source187_sites"]!=187 or status["original_2022_TP_source_samples_flagged"]!=2471:
        raise ValueError("Original EEA 2022 187 station status QA source population changed")
    U={x["site_original_SHA256"] for x in status["original_site_quality_flags"]
       if (x["flag_values_original"].get("metadata_observationStatus",{}).get("U",0)>0)}
    if len(U)!=9:raise ValueError("Original 9 metadata U sites unexpectedly changed")
    counts=status["observation_flag_field_values_overall"]
    if counts["metadata_observationStatus"]!={"A":2462,"U":9}:
        raise ValueError("Original WISE6 metadata code distribution changed")
    if counts["resultQualityObservedValueBelowLOQ"]!={"0":2418,"1":53}:
        raise ValueError("Original LOQ code distribution changed")
    seasource={x["original_station_siteKey_scheme_country_SHA256"]:x
        for x in tp["original_source_187_sites_all_numeric_and_censor_QC"]}
    worlds={(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))):x
        for x in wc["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]}
    if len(seasource)!=187:raise ValueError("Numeric original TP source did not contain 187 sites")
    records=[]
    for r in source["original_187_WISE6_TP_stations_with_public_source_reach_exposures"]:
        idhash=r["original_WISE6_country_scheme_site_SHA256"]
        q=seasource[idhash]
        if not q["numeric_site_support_ge4days_ge3months_LOQfraction_lehalf"] or r["HYRIV_ID"] is None:
            continue
        flow=r["HydroRIVERS_original_modeled_DIS_AV_CMS"]
        catch=r["HydroRIVERS_original_upstream_basin_UPLAND_SKM"]
        if not isinstance(flow,(int,float)) or flow<=0 or not isinstance(catch,(int,float)) or catch<=0:continue
        wcrow=worlds[(r["countryCode"],int(float(r["country_1deg_lat_floor"])),int(float(r["country_1deg_lon_floor"])))]
        z=dict(r)
        z["TP_original_2022_annual_site_median_mg_P_per_L"]=q["original_EEA2022_TP_station_median_of_daily_medians_mg_P_per_L"]
        z["WC2021_512_built_land_fraction"]=wcrow["WC2021_512_built_land_fraction"]
        z["WC2021_512_crop_land_fraction"]=wcrow["WC2021_512_crop_land_fraction"]
        z["WC2021_512_country_terrestrial_area_km2"]=wcrow["WC2021_512_country_terrestrial_area_km2"]
        records.append(z)
    if len(records)!=171:raise ValueError("Exact old 171 source-only + original numeric TP input not recovered")
    original_ref=first["original_2022_river_2021_land_2022_EEA_site_source_model_results"]
    selected=[
        "strict_upstream_reported_P_mass_t_year",
        "strict_downstream_reported_P_mass_placebo",
        "source_reported_P_over_modeled_receptor_flow_mgL",
        "strict_upstream_qualified_WWTP_count",
        "strict_upstream_reported_N_mass_t_year",
    ]
    validated={}
    for name in selected:
        original_recomputed=original.blocked_cv(records,name)
        oldgain=original_ref[name]["relative_source_added_MSE_improvement_pct"]
        if abs(original_recomputed["relative_source_added_MSE_improvement_pct"]-oldgain)>1e-10:
            raise ValueError("Original geography/TP regression no longer reproduced: "+name)
        validated[name]=oldgain
    clean=[r for r in records if r["original_WISE6_country_scheme_site_SHA256"] not in U]
    if len(clean)!=162:raise ValueError("Original exactly nine U sites not excluded from 171 TP cohort")
    positive=[r for r in clean if r["strict_upstream_qualified_WWTP_count_50km"]>0]
    negative=[r for r in clean if r["strict_upstream_qualified_WWTP_count_50km"]==0]
    c=Counter(x["countryCode"] for x in positive)
    regions=set(r["countryCode"]+"|"+str(r["original_country2degree_block"]) for r in positive)
    reaches=set(r["HYRIV_ID"] for r in positive)
    if not(len(positive)>=20 and len(negative)>=20 and c["DE"]>=5 and c["AT"]>=5
           and len(regions)>=5 and len(reaches)>=15):
        raise ValueError("Conservative raw status U exclusion invalidates original station/source eligibility gate")
    results={}
    for name in selected:
        v=original.blocked_cv(clean,name,include_predictions=True)
        results[name]=v
        print("GESUI_WISE6_ORIGINAL_TP_2022_EXCLUDE_ANY_U_SITE_GEOGRAPHIC_SOURCE_NEGATIVE_CONTROL",
          json.dumps({"source":name,"n":len(clean),
            "original_gain":validated[name],
            "new_after_exclusion_gain":v["relative_source_added_MSE_improvement_pct"],
            "country":v["country_specific"]},sort_keys=True),flush=True)
    out={
      "schema":"gesui_2022_WISE6_TP_directed_WWTP_original_meta_observation_U_9station_exclusion_POSTHOC_QC_v01",
      "date_utc":datetime.now(timezone.utc).isoformat(),
      "explicitly_POSTHOC_status_exclusion_after_original_171site_model":True,
      "original_EEA2026_WISE6_187site_raw_status_QC_SHA256":sha(a.original_status),
      "original_HydroRIVERS_directed_WWTP_sourceonly_SHA256":sha(a.original_directed_source),
      "original_2022_numeric_WISE6_TP_QC_SHA256":sha(a.original_realTP_numeric),
      "original_171station_CV_source_SHA256":sha(a.original_171site_CV),
      "original_187_sites_total":187,
      "original_171_sites_with_numeric_QC_and_500m_directed_reach":len(records),
      "original_meta_U_rows_n":9,
      "original_any_U_station_SHA256_removed":sorted(U),
      "source_only_stations_remaining_no_any_U_original_flag":len(clean),
      "numeric_and_original_source_support_remaining":{
         "source_positive":len(positive),"registered_source_zero":len(negative),
         "positive_by_country":dict(c),"positive_unique_reaches":len(reaches),
         "positive_country2degree_blocks":len(regions)},
      "original_all_171_recomputed_source_gain_percent":validated,
      "conservative_drop_entire_site_with_any_original_meta_U_2degree_holdout_CV":results,
      "metadata_U_definition_from_secondary_dictionary_requires_EEA_original_confirmation":True,
      "cannot_infer_WWTP_causes_or_does_not_cause_river_total_P":True,
      "no_2022_Lycaenidae_photo_or_wildlife_abundance_accessed":True,
      "original_WISE6_values_NOT_REEXTRACTED_or_changed":True
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_ORIGINAL_2022_WISE6_FLAG_U_CONSERVATIVE_ALL_SITE_REMOVAL_COMPLETED",
      json.dumps({"original_n":len(records),"without_U":len(clean),
        "primary_upstreamP_gain":results["strict_upstream_reported_P_mass_t_year"]["relative_source_added_MSE_improvement_pct"],
        "downstream_placebo_gain":results["strict_downstream_reported_P_mass_placebo"]["relative_source_added_MSE_improvement_pct"]},sort_keys=True),flush=True)

if __name__=="__main__":main()
