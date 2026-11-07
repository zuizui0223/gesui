#!/usr/bin/env python3
"""EEA UWWTD 2022 primary reported-measured TN/TP WWTP LOAD source geometry.

Read EEA official TreatmentPlants and DischargePoints ArcGIS REST pages; only
both-positive *measured* N/P annual reported load facilities; join unique
discharge points via normalized reporter+plant ID; snap to Europe HydroRIVERS
v1-derived SHA-locked network within 500m; accumulate hydrologically upstream
facility load to existing Peters river-station IDs within 50km.

All external source observations, source units, code flags and joins are
audited. No river TN/TP outcome and no GBIF plants are used here.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
from datetime import datetime,timezone
import hashlib,json,math,re,time
from pathlib import Path
import geopandas as gpd
import numpy as np,pandas as pd
import requests

BASE="https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD"
PLANT=BASE+"/TreatmentPlants_2022/MapServer/0"
OUTFALL=BASE+"/DischargePoints_2022/MapServer/0"
BOTH_WHERE="uwwNDischargeMeasured > 0 AND uwwPDischargeMeasured > 0"
PLANT_FIELDS=("OBJECTID","repCode","uwwCode","uwwLatitude","uwwLongitude",
              "uwwNDischargeMeasured","uwwPDischargeMeasured",
              "uwwNRemoval","uwwPRemoval","uwwState","uwwWasteWaterTreated")
OUT_FIELDS=("OBJECTID","repCode","uwwCode","dcpCode","dcpLatitude","dcpLongitude",
            "dcpWaterBodyType","dcpState")
RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
SOURCE_FEATURE_SHA="47fb5bb92726602170a98473d02d2c810f98e49a51ef867d017190a85ae2f10b"
STATION_LIMIT_M=500
SOURCE_LIMIT_M=500
MAX_GRAPH_KM=50
BBOX=(-12,44,26,64)
PAGE=1000

def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def api(url,params):
    last=None
    for i in range(4):
        try:
            r=requests.get(url,params=params,timeout=(12,65),
              headers={"User-Agent":"gesui-EEA-measured-TNTP-river-source/0.1",
                       "Accept":"application/json"})
            r.raise_for_status()
            data=r.json()
            if "error" in data:raise ValueError(str(data["error"])[:400])
            return data,hashlib.sha256(r.content).hexdigest()
        except (requests.RequestException,ValueError) as e:
            last=e
            if i<3:time.sleep((i+1)*3)
    raise RuntimeError("Official EEA REST request failed: "+str(last))

def download_layer(base,where,fields,expected):
    c,_=api(base+"/query",{"f":"json","where":where,"returnCountOnly":"true"})
    n=int(c["count"])
    if expected is not None and n!=expected:
        raise ValueError(f"Official source count drift {base}: {n} vs {expected}")
    data=[]
    receipts=[]
    for offset in range(0,n,PAGE):
        d,digest=api(base+"/query",{"f":"json","where":where,
             "returnGeometry":"false","outFields":",".join(fields),
             "orderByFields":"OBJECTID ASC",
             "resultOffset":offset,"resultRecordCount":PAGE})
        part=[x["attributes"] for x in d.get("features",[])]
        if not part or len(part)>PAGE:
            raise ValueError("EEA returned empty/truncated API page at "+str(offset))
        data.extend(part)
        receipts.append({"offset":offset,"rows":len(part),"api_body_sha256":digest})
    if len(data)!=n:raise ValueError(f"EEA missing source pages {len(data)}/{n}")
    df=pd.DataFrame(data)
    if set(fields)-set(df.columns):raise ValueError("Official layer missing fields")
    if df.OBJECTID.duplicated().any():raise ValueError("Layer OBJECTID duplicate across pages")
    return df,receipts

def normalize_country_reporter(reporter):
    value=str(reporter or "").strip().upper()
    m=re.match(r"^([A-Z]{2})[_-]",value)
    return m.group(1) if m else ""

def normalize_key(frame):
    frame=frame.copy()
    # EEA repCode isn't standardized to a two-letter prefix in every country.
    # The original crosswalk accidentally dropped 3,111 of 5,492 eligible
    # facilities solely because it enforced an optional country-code regex.
    # Use the FULL official reporting-cycle code with the official facility
    # uwwCode, exactly as present in *both* EEA public layers.
    # This source-key repair uses no observed TN/TP receptor chemistry outcomes.
    frame["country_reported"]=frame.repCode.map(normalize_country_reporter)
    frame["plant_key"]=frame.uwwCode.astype("string").str.strip()
    frame["report_code_norm"]=frame.repCode.astype("string").str.strip().str.upper()
    frame["compound_join_key"]=frame["report_code_norm"].fillna("")+"|"+frame["plant_key"].fillna("")
    frame.loc[
        frame.report_code_norm.isna() | frame.report_code_norm.eq("") |
        frame.plant_key.isna() | frame.plant_key.eq(""),
        "compound_join_key"
    ]=None
    return frame

def gdf_points(df,lat,lon):
    return gpd.GeoDataFrame(df.copy(),
       geometry=gpd.points_from_xy(df[lon],df[lat]),crs="EPSG:4326")

def snap(points,rivers,limit_m):
    left=points.to_crs("EPSG:3035")
    right=rivers[["HYRIV_ID","geometry"]].to_crs("EPSG:3035")
    assigned=gpd.sjoin_nearest(left,right,how="left",
                max_distance=limit_m,distance_col="snap_distance_m")
    assigned=assigned.sort_values(["source_row_id","snap_distance_m","HYRIV_ID"],
                        na_position="last").drop_duplicates("source_row_id")
    return assigned

def accum_upstream(source_rows,river):
    totals=defaultdict(lambda:[0,0.,0.,0])
    byreach=defaultdict(list)
    source_events=0
    for row in source_rows.itertuples(index=False):
        orig=int(row.HYRIV_ID)
        n=float(row.uwwNDischargeMeasured)
        p=float(row.uwwPDischargeMeasured)
        if orig not in river or not math.isfinite(n) or not math.isfinite(p) or n<=0 or p<=0:
            continue
        byreach[orig].append((n,p))
        # outfall is on originating reach with unknown ordering; do not assign there
        seen={orig}
        cursor=river[orig][0]
        dist=float(river[orig][1])
        while cursor in river and cursor>0 and cursor not in seen and dist<=MAX_GRAPH_KM:
            totals[cursor][0]+=1
            totals[cursor][1]+=n
            totals[cursor][2]+=p
            totals[cursor][3]+=1
            source_events+=1
            seen.add(cursor)
            seg=river[cursor]
            dist+=max(0.,float(seg[1]))
            cursor=seg[0]
    return totals,byreach,source_events

def downstream_placebo(reach,river,byreach):
    seen={reach}
    current=river[reach][0]
    distance=float(river[reach][1])
    cnt=0
    n=0.
    p=0.
    while current in river and current>0 and current not in seen and distance<=MAX_GRAPH_KM:
        for v in byreach.get(current,[]):
            cnt+=1;n+=v[0];p+=v[1]
        seen.add(current);distance+=max(0.,float(river[current][1]))
        current=river[current][0]
    return cnt,n,p

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--river-parquet",type=Path,required=True)
    arg.add_argument("--original-monitor-source",type=Path,required=True)
    arg.add_argument("--out",type=Path,required=True)
    a=arg.parse_args()
    if sha256(a.river_parquet)!=RIVER_SHA:raise ValueError("Changed HydroRIVERS EU source SHA")
    if sha256(a.original_monitor_source)!=SOURCE_FEATURE_SHA:
        raise ValueError("Original 2026 Peters station source-only selection SHA mismatch")
    WW,prec=download_layer(PLANT,BOTH_WHERE,PLANT_FIELDS,5492)
    DP,drec=download_layer(OUTFALL,"1=1",OUT_FIELDS,24416)
    WW=normalize_key(WW)
    DP=normalize_key(DP)
    quality={}
    quality["unrecognized_reporter_prefix_wwtp"]=int(WW.country_reported.eq("").sum())
    quality["unrecognized_reporter_prefix_dcp"]=int(DP.country_reported.eq("").sum())
    quality["wwtp_missing_compound_key"]=int(WW.compound_join_key.isna().sum())
    quality["dcp_missing_compound_key"]=int(DP.compound_join_key.isna().sum())
    WW=WW.dropna(subset=["compound_join_key"]).copy()
    DP=DP.dropna(subset=["compound_join_key"]).copy()
    wcount=WW.compound_join_key.value_counts()
    WW=WW[WW.compound_join_key.map(wcount).eq(1)].copy()
    counts=DP.compound_join_key.value_counts()
    WW["discharge_points_per_plant"]=WW.compound_join_key.map(counts).fillna(0).astype(int)
    quality["initial_both_measured_n"]=5492
    quality["unique_country_and_plant_key"]=len(WW)
    quality["matched_any_dcp"]=int(WW.discharge_points_per_plant.gt(0).sum())
    quality["unambiguous_one_dcp"]=int(WW.discharge_points_per_plant.eq(1).sum())
    quality["multiple_dcp_cannot_allocate_total_mass"]=int(WW.discharge_points_per_plant.gt(1).sum())
    quality["no_discharge_point"]=int(WW.discharge_points_per_plant.eq(0).sum())
    unique=WW[WW.discharge_points_per_plant.eq(1)].merge(
        DP.drop_duplicates("compound_join_key"),on="compound_join_key",how="inner",
        suffixes=("_treatment","_dcp"),validate="one_to_one")
    for k in ("dcpLatitude","dcpLongitude","uwwPDischargeMeasured","uwwNDischargeMeasured"):
        unique[k]=pd.to_numeric(unique[k],errors="coerce")
    unique=unique[unique.dcpLatitude.between(BBOX[1],BBOX[3])&
                  unique.dcpLongitude.between(BBOX[0],BBOX[2])].copy()
    quality["unique_dcp_in_europe_study_bbox"]=len(unique)
    quality["dcp_waterbody_code_distribution"]={
        str(k):int(v) for k,v in unique.dcpWaterBodyType.fillna("MISSING").value_counts().head(20).items()}
    g=gpd.read_parquet(a.river_parquet).to_crs("EPSG:4326")
    g=g.rename(columns={x:x.upper() for x in g.columns if x!=g.geometry.name})
    if set(("HYRIV_ID","NEXT_DOWN","LENGTH_KM"))-set(g.columns):
        raise ValueError("Source river schema changed")
    g=g.cx[BBOX[0]:BBOX[2],BBOX[1]:BBOX[3]].copy()
    for k in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM"):
        g[k]=pd.to_numeric(g[k],errors="coerce")
    g=g[g.HYRIV_ID.notna()].drop_duplicates("HYRIV_ID").copy()
    river={int(r.HYRIV_ID):(int(r.NEXT_DOWN) if pd.notna(r.NEXT_DOWN) else 0,float(r.LENGTH_KM))
           for r in g[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    unique["source_row_id"]=np.arange(len(unique),dtype=int)
    zipped=snap(gdf_points(unique,"dcpLatitude","dcpLongitude"),g,SOURCE_LIMIT_M)
    zipped=zipped[zipped.HYRIV_ID.notna()].copy()
    zipped["HYRIV_ID"]=zipped.HYRIV_ID.astype(int)
    zipped=zipped[zipped.HYRIV_ID.isin(river)].copy()
    quality["dcp_snapped_to_river_within_500m"]=len(zipped)
    quality["dcp_snap_m_median"]=float(zipped.snap_distance_m.median()) if len(zipped) else None
    quality["dcp_snap_m_max"]=float(zipped.snap_distance_m.max()) if len(zipped) else None
    totals,byreach,links=accum_upstream(zipped,river)
    quality["source_downstream_50km_aggregate_connections"]=links
    sites=pd.read_csv(a.original_monitor_source,dtype={"WQ.ID":str},low_memory=False)
    if set(sites.nutrient)!={"TP","TN"} or sites[["nutrient","WQ.ID"]].duplicated().any():
        raise ValueError("Original river monitoring station source table changed")
    stations=sites.copy()
    stations["eea2022_upstream_wwtp_both_positive_count"]=stations.HYRIV_ID.map(
        lambda k: totals.get(int(k),[0,0.,0.,0])[0])
    stations["eea2022_upstream_N_measured_load_nativeunit"]=stations.HYRIV_ID.map(
        lambda k: totals.get(int(k),[0,0.,0.,0])[1])
    stations["eea2022_upstream_P_measured_load_nativeunit"]=stations.HYRIV_ID.map(
        lambda k: totals.get(int(k),[0,0.,0.,0])[2])
    d=[downstream_placebo(int(k),river,byreach) for k in stations.HYRIV_ID]
    stations["eea2022_downstream_wwtp_count_placebo"]=[x[0] for x in d]
    stations["eea2022_downstream_N_measured_load_nativeunit_placebo"]=[x[1] for x in d]
    stations["eea2022_downstream_P_measured_load_nativeunit_placebo"]=[x[2] for x in d]
    stats={}
    for par,part in stations.groupby("nutrient"):
        stats[par]={
         "fixed_2011_2024_eligible_station_count":len(part),
         "source_eea_measured_NP_any_upstream_50km":int((part.eea2022_upstream_wwtp_both_positive_count>0).sum()),
         "source_eea_measured_NP_any_downstream_50km":int((part.eea2022_downstream_wwtp_count_placebo>0).sum()),
         "source_eea_published_both_load_positive_stations":int((part.eea2022_upstream_N_measured_load_nativeunit.gt(0)&
                  part.eea2022_upstream_P_measured_load_nativeunit.gt(0)).sum()),
         "spatial_3degree_blocks_with_upstream_source":int(part.loc[part.eea2022_upstream_wwtp_both_positive_count.gt(0),"spatial_3deg_block"].nunique()),
         "source_wwtp_measured_load_N_total_at_receptors_unweighted_nativeunit":
             float(part.eea2022_upstream_N_measured_load_nativeunit.sum()),
         "source_wwtp_measured_load_P_total_at_receptors_unweighted_nativeunit":
             float(part.eea2022_upstream_P_measured_load_nativeunit.sum()),
         "source_covariates_crop_pop_available_with_upstream_source":int(
             (part.eea2022_upstream_wwtp_both_positive_count.gt(0)&part[["crop2","pop_dens"]].notna().all(axis=1)).sum()),
        }
        print("GESUI_EEA_MEASURED_SOURCE_RIVER_STATION_SUPPORT",
              json.dumps({par:stats[par]},sort_keys=True),flush=True)
    a.out.mkdir(parents=True,exist_ok=True)
    out_cols=["nutrient","WQ.ID","WQ.lat","WQ.lon","HYRIV_ID","spatial_3deg_block",
      "DIS_AV_CMS","UPLAND_SKM","start","end","sum","crop2","pop_dens","N_deposition",
      "eea2022_upstream_wwtp_both_positive_count",
      "eea2022_upstream_N_measured_load_nativeunit","eea2022_upstream_P_measured_load_nativeunit",
      "eea2022_downstream_wwtp_count_placebo",
      "eea2022_downstream_N_measured_load_nativeunit_placebo",
      "eea2022_downstream_P_measured_load_nativeunit_placebo"]
    stations[out_cols].to_csv(a.out/"original_Peters_station_EEA_measured_load_source_only.csv",
                               index=False,float_format="%.10g")
    rec={
      "schema":"gesui_EEA_UWWTD2022_measured_NP_2022_to_Peters_river_monitor_source_v01",
      "recorded_at_utc":datetime.now(timezone.utc).isoformat(),
      "source_plants_API_count_expected":5492,
      "source_outfalls_API_count_expected":24416,
      "source_ArcGIS_original_pages_NP":[len(prec),len(drec)],
      "source_ArcGIS_all_page_sha256":prec+drec,
      "official_source_base":BASE,
      "source_2022_reference_year":True,
      "source_reported_measured_mass_units_need_dictionary_crosscheck":True,
      "hydrorivers_GeoParquet_SHA256":RIVER_SHA,
      "original_monitor_station_source_SHA256":SOURCE_FEATURE_SHA,
      "UWWTD_discharge_point_source_quality":quality,
      "monitoring_station_source_support":stats,
      "no_Peters_measurement_concentrations_accessed":True,
      "no_plants_or_GBIF_records_accessed":True,
      "estimated_and_calculated_UWWTD_NP_loads_excluded":True,
      "ambiguous_multiple_outlets_excluded":True,
      "same_river_reach_source_excluded_due_to_unknown_inreach_order":True,
      "upstream_path_max_km":MAX_GRAPH_KM,
      "negative_control_downstream_loads":True,
      "caveats":[
        "Reported EEA facility load is a 2022 value and does not reconstruct 2011-2024",
        "WWTP annual load requires data dictionary unit check; source nativeunit stored without conversion",
        "Only facilities with positive *both* reported measured N and P are included; source selection undercoverage",
        "Multiple outlets excluded because single facility mass cannot be uniquely allocated to one discharge point",
        "Third-party converted HydroRIVERS Europe network, not independently official-original topology",
        "2022 receptor nutrient measurements are NOT yet acquired in this source gate",
        "Nearby agricultural or other point nutrient sources not included"
      ]}
    (a.out/"EEA_MEASURED_LOAD_SOURCE_QC.json").write_text(json.dumps(rec,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA2022_REPORTED_SOURCE_GATE_FINAL",json.dumps({
       "joined_unique_europe_outfalls":quality["unique_dcp_in_europe_study_bbox"],
       "river_snapped_reported_NP_facilities":quality["dcp_snapped_to_river_within_500m"],
       "monitoring_source_support":stats,"source_only":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
