#!/usr/bin/env python3
"""Official HydroRIVERS v1 NEXT_DOWN source-only screening for real 2017/18 EA TP.

1. Original EIDC ChemPop biodiversity source-qualified BIO_SITE_ID (no BMWP
   scores opened), year-specific 2017/2018 monitoring dates.
2. Original Environment Agency final treated Water Company (SA) total
   phosphorus as P actual observations, same year, numeric mg/L, at native
   samplingPoint.notation; no fabricated outlet or annual mass.
3. Original HydroSHEDS HydroRIVERS v1 Europe + Middle East (68MB zip) with
   true directed NEXT_DOWN geometry, snapped to original site/sampling
   point BNG coordinates within 500m.

The WQE SA sampling POINT may lie within an STW site rather than precisely
at the registered physical outlet. Hence even if HydroRIVERS reach is
upstream, this is a directionally SCREENED candidate, NOT direct proof of
hydraulic sewage nutrient delivery. Missing smaller headwater river reaches
must be reported as missing routing support, never a "zero discharge" sample.

Do not rehost original station names, IDs, geometries or chemical responses,
only aggregate QA counts, original ZIP/member SHA and STOP conditions.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,shutil,sys,tempfile,urllib.request,zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as nat
import ea_WQE_national_SA_actual2017_18_TP_chemistry_bio_screen as hist
import ea_SA_TP_to_ChemPop_location_sourcegate as prior

HYDRO_URL="https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip"
HYDRO_SHA="500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5"
EXPECTED_REACH_EU=938544
MAX_ARCHIVE_BYTES=130*1024*1024
WGS84_BBOX=(-9.0,49.0,4.0,58.0)  # downstream up to 50km from England sites
NETWORK_MAX_KM=50.0
SNAP_MAX_M=500.0
NEAR_TIE_THRESHOLD_M=1.0
N_MIN_DOWNSTREAM_BIO_SITE=50
N_MIN_SITEYEAR=100
N_MIN_HYDROBAS_L12=5
STW_EXPECTED={2017:67,2018:71}
BIO_EXPECTED={2017:389,2018:366}


def get_rivers(base):
    f=base/"eu_original.zip"
    sha=hashlib.sha256();total=0
    req=urllib.request.Request(HYDRO_URL,headers={
        "User-Agent":"gesui-source-only-HydroRIVERS-v1-England-physical-flow-screen"})
    with urllib.request.urlopen(req,timeout=140) as resp,f.open("wb") as dest:
        if resp.status!=200:raise RuntimeError("HydroRIVERS official source HTTP != 200")
        while True:
            block=resp.read(1048576)
            if not block:break
            total+=len(block)
            if total>MAX_ARCHIVE_BYTES:raise RuntimeError("HydroRIVERS original source exceeds cap")
            sha.update(block);dest.write(block)
    if sha.hexdigest()!=HYDRO_SHA:
        raise RuntimeError("HydroRIVERS official original ZIP hash changed from source-locked run")
    with zipfile.ZipFile(f) as z:
        names=z.namelist()
        members=[n for n in names if n.lower().endswith(".shp")]
        if len(members)!=1:raise RuntimeError("HydroRIVERS EU source has ambiguous shapefile members")
        target=base/"hydrorivers"
        target.mkdir()
        source_name=Path(members[0]).stem
        for name in names:
            if (Path(name).stem!=source_name or
                Path(name).suffix.lower() not in {".shp",".shx",".dbf",".prj",".cpg"}):
                continue
            out=target/Path(name).name
            with z.open(name) as src,out.open("wb") as dst:
                shutil.copyfileobj(src,dst)
        shp=target/(source_name+".shp")
    g=gpd.read_file(str(shp),bbox=WGS84_BBOX,engine="pyogrio")
    if g.crs is None or g.crs.to_epsg()!=4326:
        raise RuntimeError("Original HydroRIVERS source WGS84 CRS not independently recognized")
    fields={"HYRIV_ID","NEXT_DOWN","LENGTH_KM","MAIN_RIV","HYBAS_L12"}
    if not fields.issubset(g.columns):
        raise RuntimeError("Directed source network required fields missing")
    if g.HYRIV_ID.duplicated().any():raise RuntimeError("HydroRIVERS source reach IDs duplicated")
    if len(g)<3000:raise RuntimeError("Official UK river network bbox source unexpectedly sparse")
    g=g.to_crs("EPSG:27700")
    g["HYRIV_ID"]=pd.to_numeric(g.HYRIV_ID,errors="raise").astype("int64")
    g["NEXT_DOWN"]=pd.to_numeric(g.NEXT_DOWN,errors="coerce").fillna(0).astype("int64")
    g["LENGTH_KM"]=pd.to_numeric(g.LENGTH_KM,errors="coerce")
    if g.LENGTH_KM.isna().any() or (g.LENGTH_KM<=0).any():
        raise RuntimeError("HydroRIVERS original reach length null/nonpositive")
    graph={
       int(r.HYRIV_ID):(int(r.NEXT_DOWN),float(r.LENGTH_KM),str(r.MAIN_RIV),str(r.HYBAS_L12))
       for r in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM","MAIN_RIV","HYBAS_L12"]].itertuples(index=False)
    }
    return g,graph,{
      "original_official_HydroRIVERS_EU_ZIP_SHA256":HYDRO_SHA,
      "network_reaches_in_england_padded_WGS84_bbox":len(g),
      "original_BNG_transformed_from_official_WGS84_geometry":True,
      "original_next_down_present":True,
      "river_lines_omit_small_headwaters_below_10km2_or_0_1_cms":True,
      "original_river_v1_release_not_modern_time_variable_drainage":True,
    }


def snap(points,rivers,kind):
    # Input points are native EA/EIDC BNG source coordinates; do not allow
    # arbitrary lon/lat or broken zero values.
    keys=sorted(points)
    frame=gpd.GeoDataFrame(
       {"point_key":keys},
       geometry=[Point(points[k]) for k in keys],
       crs="EPSG:27700")
    target=rivers[["HYRIV_ID","geometry"]]
    joined=gpd.sjoin_nearest(frame,target,how="left",
                            max_distance=SNAP_MAX_M,distance_col="distance_m")
    records={}
    ambiguous=0
    far=0
    for key,rows in joined.groupby("point_key",sort=False):
        rows=rows[rows.HYRIV_ID.notna()]
        if len(rows)==0:
            far+=1;continue
        mi=float(rows.distance_m.min())
        ties=rows[rows.distance_m-mi<=NEAR_TIE_THRESHOLD_M]
        if ties.HYRIV_ID.nunique()!=1:
            ambiguous+=1;continue
        rr=ties.sort_values(["HYRIV_ID"],kind="stable").iloc[0]
        records[key]={"reach":int(rr.HYRIV_ID),"snap_distance_m":mi}
    out={"point_type":kind,"input_points_n":len(points),
         "snapped_to_unique_nearest_reach_within_500m_n":len(records),
         "no_river_within_500m_n":far,
         "ambiguous_confluence_tied_reach_n":ambiguous,
         "median_snap_distance_m":float(np.median([v["snap_distance_m"] for v in records.values()])) if records else None,
         "all_raw_site_ids_and_coordinates_omitted":True}
    if len(records)+far+ambiguous!=len(points):
        raise RuntimeError("Snapping source disjoint QA counts do not add to input")
    return records,out


def downstream_trace(start,graph,maxkm=NETWORK_MAX_KM):
    result={}
    seen={start}
    k=int(start)
    dist=0.0
    for step in range(2000):
        if k not in graph:return result,"left_bbox_or_network"
        nxt,length,*_=graph[k]
        dist+=length
        if nxt<=0 or nxt in seen:return result,"terminal_or_cycle"
        if dist>maxkm:return result,"over_50km"
        if nxt not in graph:return result,"downstream_reach_not_available_in_bbox"
        result[nxt]=dist
        seen.add(nxt)
        k=nxt
    return result,"path_too_long"


def by_year_route(yr,bio,stw_source,survey,biosnap,stwsnap,graph):
    qualified_sites=sorted(survey[yr])
    station_codes=sorted(stw_source[yr])
    bio_unique={k:biosnap[k] for k in qualified_sites if k in biosnap}
    stw_unique={k:stwsnap[k] for k in station_codes if k in stwsnap}
    counts=collections.Counter()
    delivered_bio=set()
    reverse_bio=set()
    ambiguous_same_reach=set()
    delivered_source=set()
    reachable_HYBAS_L12=set()
    parent_main_riv=set()
    directed_network_distances=[]
    per_bio_source=collections.Counter()
    src_paths={}
    for source_id,r in stw_unique.items():
        src_paths[source_id]=downstream_trace(r["reach"],graph)
    bio_paths={}
    for bio_id,r in bio_unique.items():
        bio_paths[bio_id]=downstream_trace(r["reach"],graph)
    for b_id,b in bio_unique.items():
        b_reach=b["reach"]
        for s_id,s in stw_unique.items():
            s_reach=s["reach"]
            counts["all_possible_pairs_with_both_reaches_snapped"]+=1
            if s_reach==b_reach:
                counts["same_HydroRIVERS_reach_ORDER_UNKNOWN"]+=1
                ambiguous_same_reach.add(b_id);continue
            ds,why=src_paths[s_id]
            if b_reach in ds:
                counts["source_reach_upstream_reaching_BIO_reach_within50km"]+=1
                delivered_bio.add(b_id)
                delivered_source.add(s_id)
                directed_network_distances.append(ds[b_reach])
                per_bio_source[b_id]+=1
                reachable_HYBAS_L12.add(graph[b_reach][3])
                parent_main_riv.add(graph[b_reach][2])
                continue
            rb,why=bio_paths[b_id]
            if s_reach in rb:
                counts["source_reach_DOWNSTREAM_of_BIO_placebo_within50km"]+=1
                reverse_bio.add(b_id)
            else:
                counts["different_network_branch_or_longer_than50km_or_missing"]+=1
    return {
      "year":yr,
      "n_source_BIO_sites_eligible":len(qualified_sites),
      "n_stw_tp_final_SA_native_source_ids_observed_in_year":len(station_codes),
      "n_source_BIO_nearest_river_snap_pass":len(bio_unique),
      "n_source_STW_SA_nearest_river_snap_pass":len(stw_unique),
      "pair_relation_counts":dict(counts),
      "n_unique_BIO_sites_on_potentially_downstream_reach_from_SA":len(delivered_bio),
      "n_unique_STW_SA_source_monitoring_sites_on_potential_upstream_reach":len(delivered_source),
      "n_unique_BIO_sites_with_STW_DOWNSTREAM_reverse_placebo":len(reverse_bio),
      "n_unique_BIO_sites_with_same_reach_source_order_unknown":len(ambiguous_same_reach),
      "n_true_HYBAS_L12_parent_blocks_among_routed_BIO":len(reachable_HYBAS_L12),
      "n_MAIN_RIV_network_components_among_routed_BIO":len(parent_main_riv),
      "network_upstream_path_km_median":float(np.median(directed_network_distances)) if directed_network_distances else None,
      "n_sites_with_two_or_more_network_upstream_STW_sources":sum(v>=2 for v in per_bio_source.values()),
      "n_bio_sites_with_no_historical_STW_SA_near_source_snapped_river_or_directed_path_NOT_zero_exposure":
            len(qualified_sites)-len(delivered_bio),
      "river_network_is_HydroRIVERS_approximate_NOT_verification_of_physical_outlet":True,
      "same_reach_distance_order_unknown_excluded":True,
      "original_source_STW_current_registry_2km_candidate_selection_bias":True,
    },delivered_bio


def main(a):
    out={"schema":"GESUI_20261120_TRUE_EA_STW_TP_HYDRORIVERS_ENGLAND_UPSTREAM_BIO_SOURCE_GATE_V1",
         "posthoc_source_only_no_BMWP_response_or_effect_estimate":True,
         "original_HydroRIVERS_v1_missing_small_streams_and_samplepoint_not_outfall":True,
         "no_WWTP_annual_NP_mass_or_causal_delivered_dose":True,
         "source_candidate_SA_original_frame_is_only_2km_proximity_to_2017_18_eligible_BIO":True,
         "status":"SOURCE_STOP"}
    try:
        registry,regqa=nat.collect()
        bio,year_bio,bioqa=prior.chempop_geo()
        names,candqa=hist.source_candidates(registry,year_bio,bio)
        measured,measqa=hist.candidate_history(names)
        if ({yr:len(measured[yr]) for yr in (2017,2018)}!=STW_EXPECTED):
            raise RuntimeError("Original 67/71 year-specific SA measured TP source frame changed")
        with tempfile.TemporaryDirectory(prefix="gesui_uk_directed_source_") as tmp:
            rivers,graph,sourceqa=get_rivers(Path(tmp))
            subset_bio={k:bio[k] for k in year_bio[2017]|year_bio[2018]}
            subset_stw={k:registry[k] for k in set(measured[2017])|set(measured[2018])}
            biosnap,qa_biosnap=snap(subset_bio,rivers,"ChemPop native BIO site")
            stwsnap,qa_stwsnap=snap(subset_stw,rivers,"EA SA water-company treated TP sampling point")
            source_by_year={}
            route_site_year={}
            for yr in (2017,2018):
                detail,eligible=by_year_route(
                     yr,bio,measured,year_bio,biosnap,stwsnap,graph)
                source_by_year[str(yr)]=detail
                route_site_year[str(yr)]=eligible
        union=set().union(*route_site_year.values())
        total_siteyear=sum(len(v) for v in route_site_year.values())
        n_parent_blocks=sum(x["n_true_HYBAS_L12_parent_blocks_among_routed_BIO"] for x in source_by_year.values())
        out.update({
            "official_EA_WQE_STW_SA_original_source_qa":{
                "national_SA_current_registered_n":regqa["complete_source_unique_native_SA_sampling_point_count"],
                "actual_2017_2018_STW_TP_n":measqa["historically_measured_STW_sampling_points_by_year"],
                "actual_2017_2018_STW_same_native_STW_point_n":measqa["historically_measured_STW_sampling_points_both_years"],
                "source_2km_candidate_STW_n":candqa["candidate_current_SA_sites_within_2km_source_screen"],
            },
            "official_HydroRIVERS_network_source_qa":sourceqa,
            "source_original_BIO_snap_qa":qa_biosnap,
            "source_original_STW_sampling_point_snap_qa":qa_stwsnap,
            "source_original_each_year_network_relation":source_by_year,
            "potentially_downstream_site_years_total":total_siteyear,
            "potentially_downstream_distinct_BIO_sites_union":len(union),
            "network_source_gates":{
                "original_source_2017_18_STW_and_BIO_source_year_match":True,
                "original_HydroRIVERS_directed_NEXT_DOWN_used_not_radius":True,
                "source_origin_is_TP_effluent_SAMPLING_POINT_not_verified_permit_outfall":True,
                "benthic_sites_with_potential_upstream_STW_georeach_at_least50":
                    len(union)>=N_MIN_DOWNSTREAM_BIO_SITE,
                "potentially_upstream_linked_siteyears_at_least100":
                    total_siteyear>=N_MIN_SITEYEAR,
                "at_least5_independent_hydrological_source_groups_NOT_PROVEN":
                    False,
                "original_physical_STW_outfall_location_exactly_verified":False,
                "mass_kg_per_year_or_flow_verified":False,
                "source_directed_candidate_is_NOT_confirmed_pollution_delivery":True,
                "ecological_regression_authorized":False,
            },
            "status":"OFFICIAL_2017_18_STW_SAMPLEPOINT_HYDRORIVERS_DIRECTED_SCREEN_NOT_CAUSAL"
        })
    except Exception as e:
        out.update({"status":"SOURCE_HYDROLOGY_OR_SITE_GATE_STOP",
                    "error":type(e).__name__+": "+str(e)[:780]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EA_STW_TP_HYDRORIVERS_DIRECTED_SCREEN="+json.dumps({
        "status":out["status"],"error":out.get("error"),
        "network_source":out.get("official_HydroRIVERS_network_source_qa"),
        "biological_source":out.get("source_original_BIO_snap_qa"),
        "STW_source":out.get("source_original_STW_sampling_point_snap_qa"),
        "source_years":out.get("source_original_each_year_network_relation"),
        "distinct_bio_total":out.get("potentially_downstream_distinct_BIO_sites_union"),
        "siteyears_total":out.get("potentially_downstream_site_years_total"),
        "gates":out.get("network_source_gates"),
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
