#!/usr/bin/env python3
"""Source-only 2022 Austria67 hydrologically UPSTREAM 50km actual river
network corridor opportunity from fixed original HydroRIVERS v1.

This is NOT full contributing-catchment landcover, and does NOT access
WWTP P/N, river TP, benthic EQR numerical values or 2022 NDVI.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math
from pathlib import Path
import geopandas as gpd
import numpy as np
from shapely.geometry import Point

ORIGINAL_HYDRO_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
ORIGINAL_EEA_BQE_GITBLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
BQE=Path("results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv")
REGIONS={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
EXCLUDED="ATFW54110117"
SCHEMA="GESUI_20261117_AT67_HYDRO_RIVERS_ACTUAL_UPSTREAM50KM_CORRIDOR_SOURCE_TOPOLOGY_V1"
MAX_UPSTREAM_GRAPH_DISTANCE_KM=50
MAX_VISITED_REACHES_PER_SITE=5000

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gitblob(p):
    b=p.read_bytes();return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()

def water_sites():
    if gitblob(BQE)!=ORIGINAL_EEA_BQE_GITBLOB:raise RuntimeError("Original EEA source changed")
    with BQE.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    if len(rows)!=2475:raise RuntimeError("Original BQE 2022 n2475 drift")
    out=[]
    for row in rows:
        if row["countryCode"]!="AT":continue
        pk=row["site_key"].split("|")[1]
        if pk==EXCLUDED:continue
        out.append({"source_PKINT":pk,"hydrorivers_source_reach":int(row["HYRIV_ID"]),
            "source_country":"AT","source_region":row["region_3deg_block"],
            "source_lat":float(row["lat"]),"source_lon":float(row["lon"])})
    if len(out)!=67 or collections.Counter(x["source_region"] for x in out)!=REGIONS:
        raise RuntimeError("Original selected 67 site source regional support changed")
    return out

def original_net(path):
    if sha(path)!=ORIGINAL_HYDRO_SHA:raise RuntimeError("Original immutable Europe HydroRIVERS source changed")
    df=gpd.read_parquet(path)
    if df.crs is None:raise RuntimeError("Original HydroRIVERS lacks source geospatial CRS")
    df=df.rename(columns={k:k.upper() for k in df.columns if k!=df.geometry.name})
    req={"HYRIV_ID","NEXT_DOWN","LENGTH_KM"}
    if not req.issubset(df.columns):
        raise RuntimeError("Native HydroRIVERS lacks expected original reach length/topology "+str(sorted(set(req)-set(df.columns))))
    ids=df["HYRIV_ID"].astype(np.int64).to_numpy()
    if len(ids)!=len(set(ids)):
        raise RuntimeError("Native original river ID key not unique")
    down={}
    ups=collections.defaultdict(list)
    length={}
    for ident,nxt,ln in zip(ids,df["NEXT_DOWN"].to_numpy(),df["LENGTH_KM"].to_numpy()):
        ident=int(ident)
        try:v=int(nxt)
        except(TypeError,ValueError):v=-1
        try:l=float(ln)
        except(TypeError,ValueError):l=math.nan
        if not math.isfinite(l) or l<0:
            raise RuntimeError("Actual HydroRIVERS length invalid; source must STOP")
        length[ident]=max(l,0.001)
        down[ident]=v
        if v>0:ups[v].append(ident)
    return df,down,ups,length

def upstream(uid,up_map,length):
    if uid not in length:return None
    todo=[(x,0.,(uid,)) for x in up_map.get(uid,())]
    reached=[]
    seen=set()
    while todo:
        rid,dist,anc=todo.pop()
        if rid in seen:continue
        seen.add(rid)
        if len(seen)>MAX_VISITED_REACHES_PER_SITE:
            return {"source_STOP":True,"source_reason":"upstream_reaches_more_than5000"}
        if rid in anc:
            return {"source_STOP":True,"source_reason":"source_upstream_topology_cycle"}
        if rid not in length:
            return {"source_STOP":True,"source_reason":"original_upstream_source_reach_missing"}
        nkm=min(float(length[rid]),MAX_UPSTREAM_GRAPH_DISTANCE_KM-dist)
        if nkm<=0:continue
        reached.append((rid,float(dist),float(nkm)))
        if dist+length[rid]<MAX_UPSTREAM_GRAPH_DISTANCE_KM:
            for other in up_map.get(rid,()):
                todo.append((other,dist+length[rid],anc+(rid,)))
    return {"source_STOP":False,
        "source_strict_upstream_branching_reach_count":len(reached),
        "source_network_upstream_km_total_summed_all_branches":float(sum(z[2] for z in reached)),
        "source_full_reaches":sum(abs(z[2]-length[z[0]])<1e-6 for z in reached),
        "source_terminal_partials":sum(z[2]<length[z[0]] for z in reached)}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    original=water_sites()
    df,down,up,length=original_net(a.rivers)
    result=[];fail=[]
    for site in original:
        rid=site["hydrorivers_source_reach"]
        info=upstream(rid,up,length)
        if info is None:info={"source_STOP":True,"source_reason":"original_monitor_reach_ID_not_in_HydroRIVERS"}
        q={**site,**info}
        result.append(q)
        if info["source_STOP"]:fail.append(q)
    valid=[q for q in result if not q["source_STOP"] and
          q["source_strict_upstream_branching_reach_count"]>=1 and
          q["source_network_upstream_km_total_summed_all_branches"]>=10]
    reg=collections.Counter(x["source_region"] for x in valid)
    gate=len(valid)>=60 and len(reg)==5 and min(reg.values())>=5
    out={"schema":SCHEMA,
      "source_2022_actual_river_benthic_stations_67":len(original),
      "source_original_HydroRIVERS_SHA256":ORIGINAL_HYDRO_SHA,
      "source_all_original_EU_HydroRIVERS_n":len(df),
      "source_EEA_original_BQE_2022_gitblob":ORIGINAL_EEA_BQE_GITBLOB,
      "actual_unmodified_2022_station_native_HYRIV_reach_IDs":len(set(s["hydrorivers_source_reach"] for s in original)),
      "river_network_distance_per_branch_km":MAX_UPSTREAM_GRAPH_DISTANCE_KM,
      "source_outcome_bearing_river_reach_itself_EXCLUDED_from_upstream_landcover_sampling":True,
      "full_contributing_catchment_polygon_NOT_OBSERVED":True,
      "source_upstream_river_CORRIDOR_not_watershed_landcover":True,
      "source_regions_original":REGIONS,
      "source_current_station_2022_upstream_50km_topology_coverage_by_region":dict(reg),
      "source_valid_for_river_corridor_landuse_qa":len(valid),
      "source_nonmatching_or_missing_upstream":fail,
      "source_all_original67_station_geometry_and_branch_counts":result,
      "source_geometric_gate_PASS":bool(gate),
      "status":"PASS_ORIGINAL_2022_AT67_50KM_RIVER_NETWORK_GEOGRAPHY" if gate else "STOP_ORIGINAL_RIVER_NETWORK_UPSTREAM_TOPOLOGY_OPPORTUNITY",
      "source_no_landcover_crops_built_points_or_2022_TP_EQR_NP_numbers_read":True,
      "source_still_no_2022_sentinel_landgreenness_5of10_prelocked_STOP":True}
    (a.output/"AT67_ORIGINAL_HYDRORIVERS_50KM_UPSTREAM_ONLY_NETWORK_QA.json").write_text(
      json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_AT67_UPSTREAM_NETWORK_ONLY_SOURCE_GATE="+json.dumps({
      "status":out["status"],
      "true_original_Austria67":len(original),
      "valid":len(valid),"original_geographic_regions":dict(reg),
      "original_source_hydroriver_reaches":len(df),
      "invalid":len(fail),
      "network_total_km_min":min((x["source_network_upstream_km_total_summed_all_branches"] for x in valid),default=None),
      "network_total_km_max":max((x["source_network_upstream_km_total_summed_all_branches"] for x in valid),default=None),
      "true_watershed_polygon_not_tested":True},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--rivers",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
