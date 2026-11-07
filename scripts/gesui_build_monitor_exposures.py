#!/usr/bin/env python3
"""Independent source-only EU station-by-river wastewater exposure panel.

Does NOT read measured TP/TN concentration values and DOES NOT read plant data.
Source: pinned HydroWASTE 58,502, v1-derived HydroRIVERS GeoParquet,
Peters2026 Figshare v5 measured-site metadata and independent land cover.
River topology is directed. Same-reach WWTPs excluded from upstream exposure.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlencode
import argparse,hashlib,io,json,math,zipfile
import geopandas as gpd
import numpy as np,pandas as pd,requests
SOURCE_RIVER_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
SOURCE_WWTPS_SHA="0506f873f6541c5f324478894a9731290d722bce"
DOI="10.25400/lincolnuninz.30546371.v5"
ARTICLE="https://api.figshare.com/v2/articles/30546371"
OUTER=(-12,44,26,64)
INTERIOR=(-11,45,25,63)
YEARS=("2011-01-01","2024-12-31")
HALO_KM=50
EARTH=6371.0088
MAX_SNAP_METERS=500
COVARS=("crop2","pasture2","pop_dens","N_deposition","olsenP",
        "area","slope_mean","forest","WWP","BIO12","Q_log")
BANNED={"closed","decommissioned","non-operational","projected",
        "proposed","under construction"}
def get(url,limit=32_000_000):
    r=requests.get(url,timeout=(10,100),stream=True,
        headers={"User-Agent":"gesui-EU-river-nutrient-study-source-only/0.1"})
    r.raise_for_status();b=bytearray()
    for x in r.iter_content(131072):
        b.extend(x)
        if len(b)>limit:raise ValueError("Source exceeds bounded file size "+url)
    r.close();return bytes(b)
def sha256(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
    return h.hexdigest()
def gitblob(p):
    dat=Path(p).read_bytes()
    return hashlib.sha1(b"blob "+str(len(dat)).encode()+b"\0"+dat).hexdigest()
def hull(x,y,low,high):
    return x.between(low[0],high[0])&y.between(low[1],high[1])
def collect_river_source(geoparquet):
    if sha256(geoparquet)!=SOURCE_RIVER_SHA:raise ValueError("HydroRIVERS GeoParquet digest unexpected")
    g=gpd.read_parquet(geoparquet).to_crs("EPSG:4326")
    g=g.rename(columns={x:x.upper() for x in g.columns if x!=g.geometry.name})
    cols=("HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","HYBAS_L12","UPLAND_SKM")
    if set(cols)-set(g.columns):raise ValueError("HydroRIVERS missing expected field")
    g=g.cx[OUTER[0]:OUTER[2],OUTER[1]:OUTER[3]].copy()
    for x in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM","DIS_AV_CMS","HYBAS_L12","UPLAND_SKM"):
        g[x]=pd.to_numeric(g[x],errors="coerce")
    g=g[g.HYRIV_ID.notna()].drop_duplicates("HYRIV_ID").copy()
    river={int(r.HYRIV_ID):{"down":int(r.NEXT_DOWN) if pd.notna(r.NEXT_DOWN) else 0,
       "length_km":float(r.LENGTH_KM),"q_m3sec":float(r.DIS_AV_CMS),
       "hybas_l12":int(r.HYBAS_L12) if pd.notna(r.HYBAS_L12) else None,
       "upland_skm":float(r.UPLAND_SKM)} for r in g.itertuples(index=False)}
    return g,river
def calc_upstream_source(ww,river):
    src=dict()
    counts=defaultdict(int);flows=defaultdict(float);opflows=defaultdict(float)
    starts=Counter();linked=0;traversals=0
    for r in ww.itertuples(index=False):
        k=int(r.HYRIV_ID)
        if k not in river:continue
        linked+=1
        q=float(r.WASTE_DIS)
        starts[k]+=1
        src[k]=src.get(k,0.)+q
        cursor=river[k]["down"]
        km=float(river[k]["length_km"])
        visited={k}
        while cursor>0 and cursor in river and cursor not in visited and km<=HALO_KM:
            counts[cursor]+=1
            flows[cursor]+=q
            if str(r.STATUS).lower()=="operational":
                opflows[cursor]+=q
            traversals+=1;visited.add(cursor)
            seg=river[cursor]
            km+=max(0,float(seg["length_km"]))
            cursor=seg["down"]
    return flows,counts,opflows,src,{"WWTP_with_matching_river_id":linked,"directed_50km_WWTP_segment_exposure_links":traversals}
def source_downstream(k,river,wwtp_same):
    n=0;q=0.;cursor=river[k]["down"];km=float(river[k]["length_km"]);seen={k}
    while cursor>0 and cursor in river and cursor not in seen and km<=HALO_KM:
        n+=int(wwtp_same.get(cursor,0)>0)
        q+=float(wwtp_same.get(cursor,0))
        seen.add(cursor);km+=max(0.,river[cursor]["length_km"]);cursor=river[cursor]["down"]
    return q,n
def load_peters():
    article=json.loads(get(ARTICLE,500_000))
    if article["doi"]!=DOI:raise ValueError("Peters source DOI version changed")
    idx={f["name"]:f for f in article["files"]}
    raw={}
    receipts=[]
    for p in ("TN","TP"):
        for suffix in ("2_meta.csv","2_years.csv","4_land.csv","2_-_Harmonised.zip"):
            name=p+"."+suffix
            f=idx[name]
            blob=get(f["download_url"],32_000_000)
            md5=hashlib.md5(blob).hexdigest()
            if md5!=f["supplied_md5"]:raise ValueError("Peters source v5 MD5 "+name)
            if name.endswith(".zip"):
                z=zipfile.ZipFile(io.BytesIO(blob))
                raw[name]=set(Path(x).stem for x in z.namelist() if x.endswith(".csv"))
            else:
                raw[name]=pd.read_csv(io.BytesIO(blob),low_memory=False,dtype={"WQ.ID":str,"WQ_ID":str})
            receipts.append({"name":name,"id":f["id"],"md5":md5,"bytes":len(blob)})
    return raw,receipts
def station_source(p,files):
    df=files[p+".2_meta.csv"].copy()
    df["WQ.lat"]=pd.to_numeric(df["WQ.lat"],errors="coerce")
    df["WQ.lon"]=pd.to_numeric(df["WQ.lon"],errors="coerce")
    df["sum"]=pd.to_numeric(df["sum"],errors="coerce")
    df["WQ.available"]=pd.to_numeric(df["WQ.available"],errors="coerce")
    df["start"]=pd.to_datetime(df["start"],errors="coerce",utc=True)
    df["end"]=pd.to_datetime(df["end"],errors="coerce",utc=True)
    years=files[p+".2_years.csv"].copy()
    years["start"]=pd.to_datetime(years["start"],errors="coerce",utc=True)
    years["end"]=pd.to_datetime(years["end"],errors="coerce",utc=True)
    ids=set(years.loc[(years.start<=pd.Timestamp(YEARS[1],tz="UTC"))&
                  (years.end>=pd.Timestamp(YEARS[0],tz="UTC")),"WQ.ID"])
    filt=df["WQ.lat"].between(INTERIOR[1],INTERIOR[3]) &df["WQ.lon"].between(INTERIOR[0],INTERIOR[2])
    filt &= df["sum"].ge(12) &df["WQ.available"].eq(1)
    filt &= (df.start<=pd.Timestamp(YEARS[1],tz="UTC")) &(df.end>=pd.Timestamp(YEARS[0],tz="UTC"))
    filt &= df["WQ.ID"].isin(ids) &df["WQ.ID"].isin(files[p+".2_-_Harmonised.zip"])
    df=df.loc[filt].drop_duplicates("WQ.ID").copy()
    land=files[p+".4_land.csv"].copy()
    if "WQ_ID" not in land.columns:raise ValueError("land source lack WQ_ID")
    if land.WQ_ID.duplicated().any():raise ValueError("source land table duplicate WQ id")
    # Critical: drop C_mean, C_median, L, Y response values from the source-provided model table;
    # use ONLY independent land-use, human and physical covariates below.
    vars=[k for k in COVARS if k in land.columns]
    land=land[["WQ_ID"]+vars].copy()
    land=land.rename(columns={"WQ_ID":"WQ.ID","WWP":"WWP_land"})
    df=df.merge(land,how="left",on="WQ.ID",validate="one_to_one")
    return df
def assign_stations(stations,g,river):
    point=gpd.GeoDataFrame(stations[["WQ.ID"]].copy(),
               geometry=gpd.points_from_xy(stations["WQ.lon"],stations["WQ.lat"]),
               crs="EPSG:4326")
    snapped=gpd.sjoin_nearest(point.to_crs("EPSG:3035"),
               g[["HYRIV_ID","geometry"]].to_crs("EPSG:3035"),
               how="left",max_distance=MAX_SNAP_METERS,distance_col="snap_m")
    snapped=snapped.sort_values(["WQ.ID","snap_m","HYRIV_ID"],na_position="last").drop_duplicates("WQ.ID")
    station=stations.merge(snapped[["WQ.ID","HYRIV_ID","snap_m"]],on="WQ.ID",how="left",validate="one_to_one")
    station=station[station.HYRIV_ID.notna()].copy()
    station["HYRIV_ID"]=station.HYRIV_ID.astype("int64")
    station=station[station.HYRIV_ID.isin(river)].copy()
    return station
def run(args):
    if gitblob(args.hydrowaste)!=SOURCE_WWTPS_SHA:raise ValueError("HydroWASTE blob mismatch")
    ww=pd.read_csv(args.hydrowaste,encoding="latin1",low_memory=False)
    if len(ww)!=58502 or ww.WASTE_ID.nunique()!=58502:raise ValueError("HydroWASTE original size invalid")
    for col in ("HYRIV_ID","QUAL_LOC","QUAL_WASTE","WASTE_DIS"):
        ww[col]=pd.to_numeric(ww[col],errors="coerce")
    src=ww[(~ww.STATUS.astype(str).str.lower().isin(BANNED))&
            ww.QUAL_LOC.le(2)&ww.QUAL_WASTE.le(3)&ww.WASTE_DIS.gt(0)&
            ww.HYRIV_ID.notna()].copy()
    src["HYRIV_ID"]=src.HYRIV_ID.astype("int64")
    g,river=collect_river_source(args.rivers)
    src=src[src.HYRIV_ID.isin(river)].copy()
    flow,counts,opflow,same,geom=calc_upstream_source(src,river)
    source_files,receipts=load_peters()
    tables=[]
    stats={}
    for nutrient in ("TN","TP"):
        df=station_source(nutrient,source_files)
        meta_n=len(df)
        df=assign_stations(df,g,river)
        for k,name in (("upstream_50km_ww_m3day",flow),
                       ("upstream_50km_ww_count",counts),
                       ("operational_reported_upstream_ww_m3day",opflow),
                       ("same_reach_WWTP_m3day",same)):
            df[k]=df.HYRIV_ID.map(name).fillna(0)
        for k in ("DIS_AV_CMS","UPLAND_SKM","HYBAS_L12"):
            f={"DIS_AV_CMS":"q_m3sec","UPLAND_SKM":"upland_skm","HYBAS_L12":"hybas_l12"}[k]
            df[k]=df.HYRIV_ID.map(lambda r:river[int(r)][f])
        downstream=[source_downstream(int(x),river,same) for x in df.HYRIV_ID]
        df["downstream_50km_ww_m3day"]=[q for q,n in downstream]
        df["downstream_50km_ww_reach_count"]=[n for q,n in downstream]
        df=df[np.isfinite(df.DIS_AV_CMS)&df.DIS_AV_CMS.ge(0.5)].copy()
        df["upstream_flow_fraction_potential"]=df.upstream_50km_ww_m3day /(
            df.upstream_50km_ww_m3day+df.DIS_AV_CMS*86400)
        df["downstream_flow_fraction_placebo"]=df.downstream_50km_ww_m3day /(
            df.downstream_50km_ww_m3day+df.DIS_AV_CMS*86400)
        df["nutrient"]=nutrient
        df["lat3_block"]=np.floor(df["WQ.lat"]/3).astype(int)
        df["lon3_block"]=np.floor(df["WQ.lon"]/3).astype(int)
        df["spatial_3deg_block"]=df.lat3_block.astype(str)+"_"+df.lon3_block.astype(str)
        land_present=df[["crop2","pop_dens"]].notna().all(axis=1)
        stats[nutrient]={
           "original_monitoring_metadata_rows":int(len(source_files[nutrient+".2_meta.csv"])),
           "source_time_and_licence_filtered_station_ids":int(meta_n),
           "valid_river_snap_and_q_ge_0_5":int(len(df)),
           "positive_upstream_ww50km":int(df.upstream_50km_ww_count.gt(0).sum()),
           "positive_downstream_ww50km":int(df.downstream_50km_ww_m3day.gt(0).sum()),
           "independent_crop_pop_source_join_n":int(land_present.sum()),
           "both_up_and_down_ww_positive":int((df.upstream_50km_ww_count.gt(0)&df.downstream_50km_ww_m3day.gt(0)).sum()),
           "independent_3deg_blocks":int(df.spatial_3deg_block.nunique()),
           "unique_HYBAS_L12":int(df.HYBAS_L12.nunique()),
           "upstream_flow_fraction_percentile_p50":float(df.upstream_flow_fraction_potential.median()) if len(df) else None,
           "upstream_flow_fraction_percentile_p90":float(df.upstream_flow_fraction_potential.quantile(.9)) if len(df) else None,
           "dropped_for_missing_independent_land_covariates":int((~land_present).sum())
        }
        tables.append(df)
        print("GESUI_NUTRIENT_MACRO_SOURCE_SUPPORT",json.dumps({nutrient:stats[nutrient]},sort_keys=True),flush=True)
    frames=pd.concat(tables,ignore_index=True)
    args.out.mkdir(parents=True,exist_ok=True)
    # Store source-only station exposure table in ephemeral artifact; never includes measured concentration
    # and never re-vendors the full Peters original observation records.
    columns=["nutrient","WQ.ID","WQ.lat","WQ.lon","HYRIV_ID","HYBAS_L12",
        "DIS_AV_CMS","UPLAND_SKM","spatial_3deg_block","start","end","sum","WQ.available",
        "upstream_50km_ww_m3day","upstream_50km_ww_count","downstream_50km_ww_m3day",
        "downstream_50km_ww_reach_count","same_reach_WWTP_m3day",
        "operational_reported_upstream_ww_m3day","upstream_flow_fraction_potential",
        "downstream_flow_fraction_placebo","snap_m"]+[c for c in COVARS if c in frames.columns]
    cols=list(dict.fromkeys(columns))
    frames[cols].to_csv(args.out/"nutrient_monitor_station_source_features.csv",index=False,float_format="%.9g")
    summary={
       "schema":"gesui_monitored_river_wwtp_upstream_source_features_v01",
       "status":"SOURCE_ONLY_NO_NUTRIENT_VALUES_ANALYSED",
       "source_git_blob_HydroWASTE":SOURCE_WWTPS_SHA,
       "source_hash_HydroRIVERS_Europe":SOURCE_RIVER_SHA,
       "Peters_figshare_doi":DOI,
       "Peters_source_original_md5_files":receipts,
       "HydroRIVERS_official_original_independently_verified":False,
       "source_date_utc":datetime.now(timezone.utc).isoformat(),
       "hydrorivers_reaches_within_outer_BBOX":len(river),
       "wwtp_source_active_or_not_reported_hq_in_BBOX":int(len(src)),
       "wwtp_upstream_network_50km":geom,
       "sampling_interior_bbox":INTERIOR,
       "source_network_outer_bbox":OUTER,
       "max_snap_m":MAX_SNAP_METERS,
       "upstream_50km_excludes_same_source_reach":True,
       "source_50km_downstream_WWTP_placebo":True,
       "measured_nutrient_values_accessed":False,
       "published_land_covariates_source_subset_only":True,
       "excluded_peters_source_columns":["C_mean","C_median","L","Y"],
       "source_support_by_nutrient":stats,
       "caveats":[
         "Peters .4_land source was assembled from model-input sites; land-use-complete subset may not represent all monitoring stations",
         "HydroWASTE facility operational status often 'Not Reported'; source WWTP volumes are estimates, not measured effluent nitrogen",
         "Mean river discharge does not capture historical or seasonal mixing",
         "Upstream and downstream river-graph exposure is not a unique wastewater signature; all municipal/agricultural urban confounds remain",
         "Some reaches truncated at BBOX edges; interior selection and residual edge diagnostics required",
         "No measured river nutrient outcome or plant community outcome inspected by this source gate",
         "A positive source geometry association is not N-load causation"
       ]}
    (args.out/"SOURCE_FEATURES_QC.json").write_text(json.dumps(summary,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_NUTRIENT_MACRO_SOURCE_GATE",json.dumps({k:v for k,v in summary.items()
        if k not in ("Peters_source_original_md5_files","caveats")},sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--hydrowaste",required=True,type=Path)
    p.add_argument("--rivers",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    run(p.parse_args())
