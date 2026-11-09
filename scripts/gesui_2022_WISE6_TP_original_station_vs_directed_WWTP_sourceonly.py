#!/usr/bin/env python3
"""SOURCE ONLY: original 187 WISE6 real 2022 river TP sites vs directed EU HydroRIVERS
and EEA 2022 both measured-positive WWTP discharges.

Never opens the numeric river TP site-year medians, 2022 Lycaenidae photos,
or any WISE6 water numeric concentration field. Reuses SHA-pinned original
gesui 2022 EEA source loader and HydroRIVERS v1 converted geometry.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,geopandas as gpd

OLD_SHA="9c07ce9de7f5104b47052c7598a943df405a4b83"
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
SPATIAL_SHA="afd64e1d5582492c09453f5a05aae0b1a9d9a38c8814778afdb5c07e43a1c974"
BBOX=(-12,44,26,64)
LIMIT_M=500.
NEAREST_TIE_M=2.
MAX_GRAPH_KM=50.
EXPECTED_SOURCE_FACILITIES=5492
EXPECTED_SOURCE_DCP=24416
EXPECTED_EU_QUALIFY_FW=4057
EXPECTED_SOURCE_SITE=187

def sha_file(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
    return h.hexdigest()

def original_blob(p):
    content=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(content)).encode()+b"\0"+content).hexdigest()

def sitehash(country,scheme,site):
    return hashlib.sha256((country+"|"+scheme+"|"+site).encode()).hexdigest()

def original_point_site_catalog(meta,positions):
    expected={}
    for cell in meta["original_42_country_1degree_photo_cells_and_WISE6_TP_stations"]:
        for station in cell["original_non_numeric_WISE6_2022_station_sampling_support"]:
            key=station["monitoring_site_ID_with_scheme_sha256"]
            if key in expected:raise ValueError("Duplicate original WISE source site")
            expected[key]={"country":cell["country"],"lat_floor":cell["1deg_lat_floor"],
                           "lon_floor":cell["1deg_lon_floor"],
                           "original_country2deg_block":cell["original_2deg_spatial_validation_block"],
                           "2022_original_sample_days":station["original_2022_TP_distinct_sample_dates"],
                           "2022_original_sample_months":station["original_2022_TP_distinct_sample_months"]}
    if len(expected)!=EXPECTED_SOURCE_SITE:
        raise ValueError("Original WISE6 2022 real TP station source count drift")
    seen={}
    duplicates=Counter()
    for key,position in positions.items():
        country,sid,scheme=key
        h=sitehash(country,scheme,sid)
        if h not in expected:continue
        if position is None:
            duplicates["source_original_station_ambiguous_geolocation"]+=1
            continue
        lat,lon=position
        if h in seen:
            duplicates["duplicate_orig_key_in_source_positions"]+=1
            continue
        if not (44<=lat<=64 and -12<=lon<=26):
            duplicates["source_original_station_outside_eu_study_bbox"]+=1
            continue
        orig=expected[h]
        if country!=orig["country"] or (math.floor(lat),math.floor(lon))!=(int(orig["lat_floor"]),int(orig["lon_floor"])):
            raise ValueError("Original WISE6 2022 country 1degree source cell shifted")
        seen[h]={**orig,"original_WISE6_country_scheme_site_SHA256":h,
                 "source_original_public_WGS84_lat":lat,
                 "source_original_public_WGS84_lon":lon,
                 "source_row_id":len(seen)}
    return seen,{"original_187_station_hashes":len(expected),
                 "geocoded_unambiguous_public_source_stations":len(seen),
                 "unmatched_original_station_hashes":len(expected)-len(seen),
                 "geocoding_issues":dict(duplicates)}

def snap_unambiguous(df,river,tag):
    if len(df)==0:return pd.DataFrame(),{"no_input_points":True}
    left=gpd.GeoDataFrame(df.copy(),
        geometry=gpd.points_from_xy(df["lon"],df["lat"]),crs="EPSG:4326").to_crs("EPSG:3035")
    right=river[["HYRIV_ID","geometry"]].to_crs("EPSG:3035")
    matches=gpd.sjoin_nearest(left,right,how="left",max_distance=LIMIT_M,distance_col="snapped_m")
    selected=[];stats=Counter()
    for key,rows in matches.groupby("source_row_id",sort=False):
        stats["all_original_points_attempted"]+=1
        valid=rows[rows["HYRIV_ID"].notna()].sort_values(["snapped_m","HYRIV_ID"])
        if not len(valid):
            stats["points_not_snapped_within_500m"]+=1
            continue
        min_dist=float(valid["snapped_m"].iloc[0])
        # sjoin_nearest returns only exact distance ties, NOT all reaches
        # within 2 metres of the best. Query all original projected lines
        # intersecting a radius (best_distance + 2m) around the point
        # before rejecting competing near-confluence assignments.
        pt=valid.geometry.iloc[0]
        candidate_offsets=right.sindex.query(
            pt.buffer(min_dist+NEAREST_TIE_M),predicate="intersects")
        competing=set()
        for offset in candidate_offsets:
            feature=right.iloc[int(offset)]
            if float(feature.geometry.distance(pt))<=min_dist+NEAREST_TIE_M+1e-7:
                competing.add(int(feature.HYRIV_ID))
        if len(competing)>1:
            stats["ambiguous_competing_nearest_reaches_within2m"]+=1
            continue
        record=rows.iloc[0].to_dict()
        record["HYRIV_ID"]=int(valid["HYRIV_ID"].iloc[0])
        record["snapped_m"]=min_dist
        record.pop("geometry",None)
        record.pop("index_right",None)
        selected.append(record)
        stats["unambiguous_snap_within500m"]+=1
    if len(df)!=stats["all_original_points_attempted"]:
        raise ValueError("Unexpected original 2022 point group expansion/loss during snap")
    if selected:
        stats["median_source_snapped_m"]=float(np.median([r["snapped_m"] for r in selected]))
        stats["max_source_snapped_m"]=float(max(r["snapped_m"] for r in selected))
    return pd.DataFrame(selected),dict(stats)

def river_topology(rivers):
    needed=["HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","UPLAND_SKM","geometry"]
    if any(c not in rivers.columns for c in needed):
        raise ValueError("Original 2022 HydroRIVERS source lacks directed topology or dilution fields")
    river={}
    for r in rivers.itertuples(index=False):
        try:
            riv=int(r.HYRIV_ID)
            down=int(r.NEXT_DOWN) if pd.notna(r.NEXT_DOWN) else 0
            length=float(r.LENGTH_KM)
            flow=float(r.DIS_AV_CMS) if pd.notna(r.DIS_AV_CMS) else None
            basin=float(r.UPLAND_SKM) if pd.notna(r.UPLAND_SKM) else None
        except (TypeError,ValueError,OverflowError):
            continue
        if not math.isfinite(length) or length<0:continue
        if riv in river:raise ValueError("Multiple conflicting directed HydroRIVERS reach IDs")
        river[riv]=(down,length,flow,basin)
    return river

def load_dilution(tpy,cms):
    if not isinstance(cms,(int,float)) or not math.isfinite(cms) or cms<=0:return None
    return tpy*1e9/(cms*(365.25*86400)*1000)

def accumulate_strict_upstream(sources,river):
    # Every source-outfall's position within reach is unknown. Source can
    # only be assigned to strictly DOWNSTREAM reaches, not its own.
    out=defaultdict(lambda:[0,0.,0.,set()])
    original_source_by_reach=defaultdict(list)
    qc=Counter()
    for r in sources.itertuples(index=False):
        origin=int(r.HYRIV_ID)
        if origin not in river:continue
        n=float(r.uwwNDischargeMeasured);p=float(r.uwwPDischargeMeasured)
        if not (math.isfinite(n) and math.isfinite(p) and n>0 and p>0):continue
        original_source_by_reach[origin].append((n,p))
        seen={origin}
        cursor=river[origin][0]
        dist=river[origin][1]
        while cursor in river and cursor>0 and cursor not in seen and dist<=MAX_GRAPH_KM:
            out[cursor][0]+=1
            out[cursor][1]+=n
            out[cursor][2]+=p
            out[cursor][3].add(origin)
            qc["strict_source_downstream_reach_assignments"]+=1
            seen.add(cursor)
            down,length,_,_=river[cursor]
            dist+=length
            cursor=down
        if cursor in seen:qc["source_reach_path_cycle_detected"]+=1
    qc["distinct_original_source_reaches_with_eligible_WWTPs"]=len(original_source_by_reach)
    qc["receptor_reaches_with_strict_upstream_source"]=len(out)
    return out,original_source_by_reach,dict(qc)

def downstream_placebo(reach,topology,byreach):
    if reach not in topology:return [0,0.,0.]
    seen={reach}
    cursor=topology[reach][0];dist=topology[reach][1]
    n=0.;p=0.;cnt=0
    while cursor in topology and cursor>0 and cursor not in seen and dist<=MAX_GRAPH_KM:
        for massN,massP in byreach.get(cursor,[]):
            cnt+=1;n+=massN;p+=massP
        seen.add(cursor);down,length,_,_=topology[cursor]
        dist+=length;cursor=down
    return [cnt,n,p]

def main():
    par=argparse.ArgumentParser()
    par.add_argument("--TP-sourceonly",type=Path,required=True)
    par.add_argument("--spatial-zip",type=Path,required=True)
    par.add_argument("--river-parquet",type=Path,required=True)
    par.add_argument("--old-source-script",type=Path,required=True)
    par.add_argument("--original-2022-API-source-QC",type=Path,required=True)
    par.add_argument("--out",type=Path,required=True)
    a=par.parse_args()
    if sha_file(a.river_parquet)!=RIVER_SHA:raise ValueError("Different third party HydroRIVERS converted v1 Europe source")
    if sha_file(a.spatial_zip)!=SPATIAL_SHA:raise ValueError("Different original public EEA2026 WISE6 spatial table source")
    expectedsha=OLD_SHA
    if original_blob(a.old_source_script)!=expectedsha:raise ValueError("Original prior 2022 EEA UWWTD directed helper script blob drift")
    sys.path.insert(0,str(a.old_source_script.parent))
    import gesui_eea_nutrient_load_river_crosswalk as prior
    import gesui_eea2026_wise6_original_2022_TP_annual_station_coverage as WISE
    source_meta=json.loads(a.TP_sourceonly.read_text())
    if not source_meta["source_preoutcome_geography_repeated_site_gate_20cells_DE10_AT6_blocks10"]:
        raise ValueError("Original WISE 187 site geography pre-outcome support failed")
    srcpos,geoqc=WISE.spatial_source(a.spatial_zip)
    sites,sq=original_point_site_catalog(source_meta,srcpos)
    ww,prec=prior.download_layer(prior.PLANT,prior.BOTH_WHERE,prior.PLANT_FIELDS,EXPECTED_SOURCE_FACILITIES)
    dp,drec=prior.download_layer(prior.OUTFALL,"1=1",prior.OUT_FIELDS,EXPECTED_SOURCE_DCP)
    original2022=json.loads(a.original_2022_API_source_QC.read_text())
    if prec+drec!=original2022["source_ArcGIS_all_page_sha256"]:
        raise ValueError("Original source-only EEA2022 WWTP API pages changed from original independently archived pre-butterfly PR9")
    WW=prior.normalize_key(ww);DP=prior.normalize_key(dp)
    WW=WW[WW.compound_join_key.notna()].copy()
    DP=DP[DP.compound_join_key.notna()].copy()
    wcount=WW.compound_join_key.value_counts()
    dcount=DP.compound_join_key.value_counts()
    WW=WW[WW.compound_join_key.map(wcount).eq(1)].copy()
    WW["original_discharge_per_plant"]=WW.compound_join_key.map(dcount).fillna(0).astype(int)
    unique=WW[WW.original_discharge_per_plant.eq(1)].merge(
        DP.drop_duplicates("compound_join_key"),on="compound_join_key",how="inner",
        suffixes=("_WWTP","_DCP"),validate="one_to_one")
    for key in ("dcpLatitude","dcpLongitude","uwwNDischargeMeasured","uwwPDischargeMeasured"):
        unique[key]=pd.to_numeric(unique[key],errors="coerce")
    unique=unique[
       unique.dcpLatitude.between(BBOX[1],BBOX[3])&
       unique.dcpLongitude.between(BBOX[0],BBOX[2])].copy()
    unique=unique[unique.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")].copy()
    if len(unique)!=EXPECTED_EU_QUALIFY_FW:
        raise ValueError(f"Original EEA 2022 qualifying FW source number changed {len(unique)} != {EXPECTED_EU_QUALIFY_FW}")
    g=gpd.read_parquet(a.river_parquet).to_crs("EPSG:4326")
    g=g.rename(columns={q:q.upper() for q in g.columns if q!=g.geometry.name})
    g=g.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    for col in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","UPLAND_SKM"):
        g[col]=pd.to_numeric(g[col],errors="coerce")
    if g.HYRIV_ID.duplicated().any():raise ValueError("Converted HydroRIVERS EU source contains duplicate reach IDs")
    g=g[g.HYRIV_ID.notna()].copy()
    river=river_topology(g)
    sourceframe=unique.copy()
    sourceframe["source_row_id"]=np.arange(len(sourceframe))
    sourceframe["lat"]=sourceframe.dcpLatitude
    sourceframe["lon"]=sourceframe.dcpLongitude
    snaps,srcqc=snap_unambiguous(sourceframe,g,"WWTP")
    if len(snaps):
        snaps=snaps[snaps.HYRIV_ID.isin(river)].copy()
    accum,byreach,graphqc=accumulate_strict_upstream(snaps,river)
    site_df=pd.DataFrame(list(sites.values()))
    site_snapped,siteqc=snap_unambiguous(site_df.assign(
            lat=site_df.source_original_public_WGS84_lat,
            lon=site_df.source_original_public_WGS84_lon),g,"WISE6")
    station_snap_map={r.original_WISE6_country_scheme_site_SHA256:r
        for r in site_snapped.itertuples(index=False)}
    records=[];support=Counter()
    for idhash,site in sites.items():
        r=station_snap_map.get(idhash)
        rid=int(r.HYRIV_ID) if r is not None else None
        stats=accum.get(rid,[0,0.,0.,set()]) if rid in river else [0,0.,0.,set()]
        down=downstream_placebo(rid,river,byreach) if rid in river else [0,0.,0.]
        flow=river[rid][2] if rid in river else None
        basin=river[rid][3] if rid in river else None
        stat={
          "original_WISE6_country_scheme_site_SHA256":idhash,
          "countryCode":site["country"],"country_1deg_lat_floor":site["lat_floor"],
          "country_1deg_lon_floor":site["lon_floor"],
          "original_country2degree_block":site["original_country2deg_block"],
          "source_WISE6_TP_2022_distinct_sample_days":site["2022_original_sample_days"],
          "source_WISE6_TP_2022_distinct_sample_months":site["2022_original_sample_months"],
          "public_monitored_WGS84_lat":site["source_original_public_WGS84_lat"],
          "public_monitored_WGS84_lon":site["source_original_public_WGS84_lon"],
          "HYRIV_ID":rid,
          "original_receptor_snap_m":float(r.snapped_m) if r is not None else None,
          "strict_upstream_qualified_WWTP_count_50km":int(stats[0]) if rid is not None else None,
          "strict_upstream_qualified_WWTP_reported_P_native_t_year":float(stats[2]) if rid is not None else None,
          "strict_upstream_qualified_WWTP_reported_N_native_t_year":float(stats[1]) if rid is not None else None,
          "distinct_upstream_qualified_source_reaches":len(stats[3]) if rid is not None else None,
          "strict_downstream_50km_WWTP_placebo_count":int(down[0]) if rid is not None else None,
          "strict_downstream_50km_WWTP_placebo_P_native_t_year":float(down[2]) if rid is not None else None,
          "HydroRIVERS_original_modeled_DIS_AV_CMS":float(flow) if flow is not None and math.isfinite(flow) else None,
          "HydroRIVERS_original_upstream_basin_UPLAND_SKM":float(basin) if basin is not None and math.isfinite(basin) else None,
          "source_apparent_P_load_over_modeled_avg_flow_mg_L":
            load_dilution(float(stats[2]),float(flow)) if rid is not None and flow is not None else None,
          "source_same_reach_omitted_unknown_relative_withinreach_order":True,
          "original_2022_TP_numeric_values_not_read_for_source_matching":True,
        }
        records.append(stat)
        support["original_2022_source_stations"]+=1
        if rid is None:support["station_not_unambiguously_snapped_within500m"]+=1
        else:
            support["station_unambiguously_snap_within500m"]+=1
            if stat["strict_upstream_qualified_WWTP_count_50km"]>0:
                support["strict_upstream_exposed_"+site["country"]]+=1
                support["strict_upstream_exposed"]+=1
            else:support["registered_qualifying_upstream_source_zero"]+=1
            if stat["strict_downstream_50km_WWTP_placebo_count"]>0:
                support["downstream_placebo_source_positive"]+=1
    positives=[x for x in records if x["strict_upstream_qualified_WWTP_count_50km"] is not None and x["strict_upstream_qualified_WWTP_count_50km"]>0]
    zero=[x for x in records if x["strict_upstream_qualified_WWTP_count_50km"]==0]
    distinct_reaches=len({x["HYRIV_ID"] for x in positives})
    posblocks=len({x["countryCode"]+"|"+str(x["original_country2degree_block"]) for x in positives})
    # Explicit source-geometry-only conservative gate. Numeric-eligible subset
    # will separately need to retain this support before any concentration fit.
    sourcegate=(len(positives)>=20 and len(zero)>=20 and
       support["strict_upstream_exposed_DE"]>=5 and support["strict_upstream_exposed_AT"]>=5 and
       distinct_reaches>=15 and posblocks>=5)
    report={
        "schema":"gesui_2022_EEA_WISE6_original_TP_187_stations_HydroRIVERS_directed_WWTP2022_SOURCE_ONLY_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2022_WISE6_metadata_QC_SHA256":sha_file(a.TP_sourceonly),
        "original_2026_WISE6_spatial_source_SHA256":SPATIAL_SHA,
        "original_original_directed_EEA2022_source_script_git_blob_SHA1":expectedsha,
        "hydrorivers_thirdparty_EU_v1_GeoParquet_SHA256":RIVER_SHA,
        "thirdparty_converted_HydroRIVERS_mirror_NOT_independently_verified_against_official":True,
        "original_WWTP_reported_both_positive_facility_rows":EXPECTED_SOURCE_FACILITIES,
        "original_EEA_DCP_source_rows":EXPECTED_SOURCE_DCP,
        "original_direct_FW_unique_outfall_before_500m_snap":len(unique),
        "original_source_API_page_response_SHA_receipts":prec+drec,
        "original_prewater_2022_EEA_API_source_QC_sha256":sha_file(a.original_2022_API_source_QC),
        "original_2022_EEA_API_source_pages_match_prior_PR9_original":True,
        "hydro_river_source_segment_count_in_EU_extent":len(river),
        "original_source_direct_FW_snap_500m_QC":srcqc,
        "original_source_directed_network_accumulation_QC":graphqc,
        "original_2022_WISE_public_station_geocoding_QC":geoqc,
        "original_2022_WISE_public_station_snap_500m_QC":siteqc,
        "original_station_source_support":dict(support),
        "strict_source_positive_reaches_distinct":distinct_reaches,
        "country2degree_blocks_with_strict_upstream_source":posblocks,
        "no_numeric_TP_opened_source_gate_pass":sourcegate,
        "original_187_WISE6_TP_stations_with_public_source_reach_exposures":records,
        "WWTP_same_reach_excluded_without_within_reach_order":True,
        "source_zero_means_no_qualified_REGISTERED_original_WWTP_upstream":True,
        "source_load_native_tonnes_year_as_reported_not_actual_2022_inriver_P":True,
        "original_2022_WISE6_numeric_TP_concentration_not_loaded_or_queried":True,
        "no_Lycaenidae_photo_individual_or_behavior_data_used":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_EEA2022_WISE6_TP_187_SOURCE_ONLY_DIRECTED_RIVER_WWTP_GATE",
      json.dumps({"total":len(records),"snapped":len(records)-support["station_not_unambiguously_snapped_within500m"],
             "upstream_positive":len(positives),"source_zero":len(zero),
             "country_source_positive":{c:support[f"strict_upstream_exposed_{c}"] for c in ("DE","AT")},
             "source_positive_unique_reaches":distinct_reaches,"positive_country2deg_blocks":posblocks,
             "source_gate":sourcegate,"original_fw_source_snap_QC":srcqc,"original_station_snap_QC":siteqc},
          sort_keys=True),flush=True)
if __name__=="__main__":main()
