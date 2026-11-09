#!/usr/bin/env python3
"""Freeze a REAL, outcome-blind directed WWTP→river 2022 TP source panel.

Join original WISE6 2022 site+source graph (187 sites) to already
source-frozen *numeric-quality eligibility flags only*, not actual measured
TP values, then original ESA WorldCover2021 country/built/crop covariates.

This is POSTHOC source-to-receptor physics after 2022 photo association was
seen. Selection never depends on TP magnitude. Numeric TP site value must be
opened only after this exact hash-frozen source selection file is committed.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

MIN_POS=20
MIN_ZERO=20
MIN_COUNTRY_POS=5
MIN_SPATIAL_POS_BLOCKS=5
MIN_POS_REACHES=15
EXPECTED_ORIGINAL_SOURCE_SITES=187
EXPECTED_NUMERIC_QUALITY_SITES=184
EXPECTED_2022_PUBLIC_WISE_ZIP_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
EXPECTED_HYDRO_SOURCE_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for x in iter(lambda:f.read(1048576),b""):h.update(x)
    return h.hexdigest()

def key(c,lat,lon):return (str(c),int(float(lat)),int(float(lon)))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--directed-source",required=True,type=Path)
    p.add_argument("--numeric-site-qc",required=True,type=Path)
    p.add_argument("--worldcover-original",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    q=json.loads(a.directed_source.read_text())
    n=json.loads(a.numeric_site_qc.read_text())
    world=json.loads(a.worldcover_original.read_text())
    if not q["no_numeric_TP_opened_source_gate_pass"] or q["hydrorivers_thirdparty_EU_v1_GeoParquet_SHA256"]!=EXPECTED_HYDRO_SOURCE_SHA:
        raise ValueError("Original directed/source-only 187-station gate not verified")
    if not q["original_2022_WISE6_numeric_TP_concentration_not_loaded_or_queried"]:
        raise ValueError("Original river graph source panel already outcome exposed")
    if not n["numeric_P_water_TP_to_Lycaenidae_photo_geographic_gate_pass"]:
        raise ValueError("Original numerical QUALITY-not-magnitude station gate not verified")
    if n["original_2022_EEA2026_WISE6_DisaggregatedData_ZIP_verified"]["original_archive_sha256"]!=EXPECTED_2022_PUBLIC_WISE_ZIP_SHA:
        raise ValueError("Original WISE6 station source version changed")
    if world["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED":
        raise ValueError("Independent 2021 ESA WorldCover country land cover source not frozen")
    available=world["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    landscapes={key(z["countryCode"],z["lat_floor"],z["lon_floor"]):z for z in available}
    if len(landscapes)!=43:raise ValueError("Original 2021 WorldCover 43 source cells changed")
    quality={z["original_station_siteKey_scheme_country_SHA256"]:z
             for z in n["original_source_187_sites_all_numeric_and_censor_QC"]}
    if len(quality)!=EXPECTED_ORIGINAL_SOURCE_SITES:
        raise ValueError("Original 187 WISE6 site numeric quality flag identifiers changed")
    if sum(bool(z["numeric_site_support_ge4days_ge3months_LOQfraction_lehalf"]) for z in quality.values())!=EXPECTED_NUMERIC_QUALITY_SITES:
        raise ValueError("Original 184 outcome-eligible site-count changed")
    points=q["original_187_WISE6_TP_stations_with_public_source_reach_exposures"]
    if len(points)!=EXPECTED_ORIGINAL_SOURCE_SITES:raise ValueError("Original station directed source frame not 187")
    record=[];counters=Counter();reach=defaultdict(list)
    for s in points:
        skey=s["original_WISE6_country_scheme_site_SHA256"]
        qc=quality.get(skey)
        if qc is None:raise ValueError("Site hash missing original numeric QC")
        if not qc["numeric_site_support_ge4days_ge3months_LOQfraction_lehalf"]:
            counters["original_numeric_site_LOQ_or_date_quality_ineligible"]+=1
            continue
        if s["HYRIV_ID"] is None:
            counters["original_site_not_undambiguously_snapped_500m"]+=1
            continue
        w=landscapes.get(key(s["countryCode"],s["country_1deg_lat_floor"],s["country_1deg_lon_floor"]))
        if w is None:raise ValueError("Originally preselected 187 water sites not in original 2021 WorldCover source cells")
        def finite_positive(field):
            v=s.get(field)
            return isinstance(v,(int,float)) and math.isfinite(v) and v>0
        if not all(finite_positive(field) for field in ("HydroRIVERS_original_modeled_DIS_AV_CMS",
                                                       "HydroRIVERS_original_upstream_basin_UPLAND_SKM")):
            counters["original_source_reach_flow_or_area_unusable"]+=1
            continue
        upcount=int(s["strict_upstream_qualified_WWTP_count_50km"])
        pMass=float(s["strict_upstream_qualified_WWTP_reported_P_native_t_year"])
        nMass=float(s["strict_upstream_qualified_WWTP_reported_N_native_t_year"])
        dp=float(s["strict_downstream_50km_WWTP_placebo_P_native_t_year"])
        if min(pMass,nMass,dp)<0 or not all(math.isfinite(x) for x in (pMass,nMass,dp)):
            raise ValueError("Original public 2022 source reported annual mass nonfinite or negative")
        if (upcount>0)!=(pMass>0):
            raise ValueError("Outfall source count does not correspond to original positive P report")
        rec={
          "original_country_scheme_site_SHA256":skey,
          "countryCode":s["countryCode"],
          "lat":float(s["public_monitored_WGS84_lat"]),
          "lon":float(s["public_monitored_WGS84_lon"]),
          "original_country_1deg_lat_floor":int(float(s["country_1deg_lat_floor"])),
          "original_country_1deg_lon_floor":int(float(s["country_1deg_lon_floor"])),
          "original_country_2degree_block":str(s["original_country2degree_block"]),
          "original_HydroRIVERS_HYRIV_ID":int(s["HYRIV_ID"]),
          "original_site_reach_snapping_distance_m":float(s["original_receptor_snap_m"]),
          "original_HydroRIVERS_modeled_flow_cms":float(s["HydroRIVERS_original_modeled_DIS_AV_CMS"]),
          "original_HydroRIVERS_modelled_upstream_catchment_km2":float(s["HydroRIVERS_original_upstream_basin_UPLAND_SKM"]),
          "original_EEA2022_upstream_qualified_source_n":upcount,
          "original_EEA2022_upstream_registered_P_t_per_year":pMass,
          "original_EEA2022_upstream_registered_N_t_per_year":nMass,
          "original_EEA2022_downstream_placebo_source_n":int(s["strict_downstream_50km_WWTP_placebo_count"]),
          "original_EEA2022_downstream_placebo_P_t_per_year":dp,
          "original_load_per_modeled_flow_proxy_mg_P_L":s["source_apparent_P_load_over_modeled_avg_flow_mg_L"],
          "country_ESAWC2021_built_fraction":float(w["WC2021_512_built_land_fraction"]),
          "country_ESAWC2021_crop_fraction":float(w["WC2021_512_crop_land_fraction"]),
          "country_ESAWC2021_terrestrial_area_km2":float(w["WC2021_512_country_terrestrial_area_km2"]),
          "source_original_2022_TP_qualifying_numeric_site_QC_without_observing_value":True,
          "no_original_river_TP_concentration_accessed_in_this_script":True,
        }
        record.append(rec);reach[int(s["HYRIV_ID"])].append(rec)
        counters["eligible_station_with_reach_and_quality"]+=1
        counters["upstream_source_positive" if upcount else "no_registered_qualifying_upstream_source"]+=1
        counters["country_"+s["countryCode"]+"_upstream" if upcount else "country_"+s["countryCode"]+"_sourcezero"]+=1
    if len(record)!=len({r["original_country_scheme_site_SHA256"] for r in record}):
        raise ValueError("Repeated station source hashes")
    for rid,vals in reach.items():
        if len(vals)<2:continue
        base=vals[0]
        for v in vals[1:]:
            if v["original_EEA2022_upstream_qualified_source_n"]!=base["original_EEA2022_upstream_qualified_source_n"] or abs(v["original_EEA2022_upstream_registered_P_t_per_year"]-base["original_EEA2022_upstream_registered_P_t_per_year"])>1e-10:
                raise ValueError("Same original reach assigned contradictory source exposure")
    positive=[r for r in record if r["original_EEA2022_upstream_qualified_source_n"]>0]
    zero=[r for r in record if r["original_EEA2022_upstream_qualified_source_n"]==0]
    positive_reaches={r["original_HydroRIVERS_HYRIV_ID"] for r in positive}
    positive_blocks={r["countryCode"]+"|"+r["original_country_2degree_block"] for r in positive}
    gate=len(positive)>=MIN_POS and len(zero)>=MIN_ZERO and counters["country_DE_upstream"]>=MIN_COUNTRY_POS and counters["country_AT_upstream"]>=MIN_COUNTRY_POS and len(positive_reaches)>=MIN_POS_REACHES and len(positive_blocks)>=MIN_SPATIAL_POS_BLOCKS
    summary={
      "schema":"gesui_WISE6_2022_directed_EEA_WWTP_source_170_TP_QC_sites_preOutcome_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "source_directed_187_WISE_station_QC_sha256":sha(a.directed_source),
      "source_WISE6_187_numeric_quality_flags_file_sha256":sha(a.numeric_site_qc),
      "original_WorldCover_2021_43_country_cells_SHA256":sha(a.worldcover_original),
      "original_full_EEA2026_WISE6_2022_zip_sha256":EXPECTED_2022_PUBLIC_WISE_ZIP_SHA,
      "source_THIRD_PARTY_HydroRIVERS_converted_2022_sha256":EXPECTED_HYDRO_SOURCE_SHA,
      "original_2022_TP_source_site_QC_exclusion_counts":dict(counters),
      "country_source_count_after_2022_water_TP_numeric_QC":{
        "DE_upstream":counters["country_DE_upstream"],"DE_sourcezero":counters["country_DE_sourcezero"],
        "AT_upstream":counters["country_AT_upstream"],"AT_sourcezero":counters["country_AT_sourcezero"]},
      "qualifying_2022_TP_stations_snapped_with_modeled_flow_and_built_crop":len(record),
      "positive_strictly_upstream_source_stations":len(positive),
      "sourcezero_original_registered_2022_qualifying_WWTP_stations":len(zero),
      "unique_original_directed_river_reaches_with_sites":len(reach),
      "positive_original_source_reaches_n":len(positive_reaches),
      "positive_source_original_country2degree_blocks":len(positive_blocks),
      "strict_original_source_support_gate_20positive_20zero_5percountry_15reaches_5blocks":gate,
      "station_source_panel_selection_never_used_actual_2022_river_TP_numeric_concentration":True,
      "station_source_panel_187_full_QC_including_unsnapped_and_ineligible_is_original_separate_result":True,
      "third_party_converted_HydroRIVERS_not_official_topology_verified":True,
      "NO_2022_TP_VALUE_OBSERVED_IN_THIS_SCRIPT":True,
      "all_qc_stations_with_original_directed_source_and_2021_WorldCover":record,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_2022_WISE6_DIRECTED_WWTP_ORIGINAL_SOURCE_SAMPLE_FROZEN_BEFORE_TP_MODEL",
      json.dumps({"sites":len(record),"sourcePositive":len(positive),"sourceZero":len(zero),
          "country":summary["country_source_count_after_2022_water_TP_numeric_QC"],
          "uniqueReaches":len(reach),"sourcePositiveReaches":len(positive_reaches),
          "positiveBlocks":len(positive_blocks),"pass":gate,
          "no_numeric_concentration_loaded":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
