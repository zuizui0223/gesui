#!/usr/bin/env python3
"""USGS NHDPlusV2 NLDI stream position/river-snap SOURCE-only on 71 wetlands.

Exactly 71 source-freeze NWCA2021 wetland sites and 86 registered POTW
permitted feature points with original 2021 numeric TN/TP effluent reports.
Source original ECHO 2026-11-08 ZIP SHA, NWCA four official CSVs and original
DMR feature dates are immutable. Do NOT read NWCA H_ALL/FQAI/NTL/PTL results.
"""
from __future__ import annotations
import argparse,hashlib,json,math,urllib.error,urllib.parse,urllib.request,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from gesui_NWCA2021_568source_EPA_NPDES_registered_sewage_outfall_geoONLY import (
    original_nwca_source_sites,ndpes_outfall_source)
D=Path("results/us_NWCA2021_plant_TNTP_20261110")
GEO=D/"EPA_NWCA2021_128WETLAND_160POTW_EXACT_NEAR_OUTFALL_FEATURE_SOURCE.json"
DMR=D/"EPA_NWCA2021_ORIGINAL_160POTW_2021_DMR_TNTP_SOURCE_FEASIBILITY.json"
GATE=D/"EPA_NWCA2021_EXACT_344POTW_OUTFALL_2021_TNTP_DMR_NHDPLUS_SOURCE_GATE.json"
BASE="https://api.water.usgs.gov/nldi/linked-data"
HEAD={"Accept":"application/json","User-Agent":"gesui-NWCA2021-2021POTW-NLDI-source-feasibility-v1 (ecological research)"}
SHARDS=4
SLEEP=.75
def sha(b):return hashlib.sha256(b).hexdigest()
def sha_file(p):return sha(p.read_bytes())
def hashkey(s):return sha(str(s).encode())
def haversine_m(lat1,lon1,lat2,lon2):
    a,b,c,d=map(math.radians,(lat1,lon1,lat2,lon2))
    h=math.sin((c-a)/2)**2+math.cos(a)*math.cos(c)*math.sin((d-b)/2)**2
    return 12742017.6*math.asin(min(1,math.sqrt(max(0,h))))
def api(url):
    """Returns JSON, SHA, HTTP status, or explicit API_ERROR, never fabricated no-river."""
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers=HEAD)
            with urllib.request.urlopen(req,timeout=48) as r:
                raw=r.read()
                val=json.loads(raw)
                if not isinstance(val,dict):raise ValueError("USGS NLDI invalid JSON root")
                return val,sha(raw),200,None
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,ValueError) as e:
            if isinstance(e,urllib.error.HTTPError) and e.code in (400,404,422):
                return None,None,e.code,str(e)[:160]
            if attempt==3:return None,None,None,(type(e).__name__+":"+str(e))[:160]
            time.sleep(2*(attempt+1))
    return None,None,None,"NLDI_REQUEST_UNRESOLVED"
def features(v):return v.get("features",[]) if isinstance(v,dict) and isinstance(v.get("features"),list) else []
def get_comid(v):
    for f in features(v):
        p=f.get("properties") or {}
        if str(p.get("comid","")).strip().isdigit():
            return str(p["comid"])
    return None
def source_points(root,nwca_dir,echo):
    geo=json.loads((root/GEO).read_text())
    dmr=json.loads((root/DMR).read_text())
    gate=json.loads((root/GATE).read_text())
    if gate["original_pre2021_DMR_exact344geo_features_SHA256"]!=sha_file(root/GEO) or gate["original_FY2021andFY2022_source_DMR_join_SHA256"]!=sha_file(root/DMR):
        raise ValueError("Original 2021 DMR actual source not frozen for NLDI")
    if not gate["source_precommitted_2021_nutrient_monitoring_feature_opportunity_gate"]["eligible_to_proceed_to_separate_NHDPlus_downstream_flow_route"]:
        raise ValueError("Original 2021 DMR source gate did not permit hydrology")
    pairs=geo["original_source_case_site_to_exact_permit_feature_sha256_and_distance"]
    measured={v["NPDES_exact_permit_feature_SHA256"] for v in dmr["source_original_2021_exact_permit_feature_DMR_reporting_summary"] if
              v["reported_original_numeric_2021_TN00600_period_end_dates"] or
              v["reported_original_numeric_2021_TP00665_period_end_dates"]}
    eligible=[p for p in pairs if p["NPDES_exact_permit_feature_SHA256"] in measured]
    wetlands={p["NWCA2021_site_visit_SHA256"] for p in eligible}
    outfallids={p["NPDES_exact_permit_feature_SHA256"] for p in eligible}
    if not(len(wetlands)==71 and len(outfallids)==86):
        raise ValueError("Original 2021 actual-feature DMR source differs 71 wetlands /86 outfalls")
    sites,orig=original_nwca_source_sites(root,nwca_dir)
    wpoints={}
    for r in sites.itertuples(index=False):
        pid=r.source_site_visit_SHA256
        if pid in wetlands:
            wpoints[pid]={"original_hashed_point_id":pid,"kind":"NWCA_wetland","lat":float(r.LAT_ANALYS),
                          "lon":float(r.LON_ANALYS),"original_state":str(r.PSTL_CODE),
                          "NPDES_permit_ids_sha256":[]}
    if len(wpoints)!=71:raise ValueError("Original fixed 71 wetland point coordinates not recovered")
    facilities,fc=ndpes_outfall_source(echo)
    facility=facilities.loc[facilities.PERMITTED_FEATURE_LATLON & facilities.PERMIT_DATES_PLAUSIBLE_2021].copy()
    outfall_coords=defaultdict(set)
    permits_by_outfall=defaultdict(set)
    for r in facility.itertuples(index=False):
        permit=str(r.PERMIT_ID).strip().upper()
        feature=str(r.FEATURE_ID).strip().upper()
        key=hashkey(permit+"|"+feature)
        if key not in outfallids:continue
        outfall_coords[key].add((round(float(r.LAT),9),round(float(r.LON),9)))
        permits_by_outfall[key].add(hashkey(permit))
    if set(outfall_coords)!=outfallids:
        raise ValueError("Frozen 86 registered NPDES POTW feature coordinates cannot all be reconstituted")
    if any(len(v)!=1 for v in outfall_coords.values()):
        # Coordinates for same permit+feature may be duplicated in weekly source.
        # Never choose favorable one from multiple original positions.
        raise ValueError("Current 2026 original ECHO has ambiguous lat/lon for an otherwise fixed feature ID")
    fpoints={}
    for key,v in outfall_coords.items():
        lat,lon=next(iter(v))
        fpoints[key]={"original_hashed_point_id":key,"kind":"POTW_outfall",
                      "lat":lat,"lon":lon,"original_state":None,
                      "NPDES_permit_ids_sha256":sorted(permits_by_outfall[key])}
    return sorted(list(wpoints.values())+list(fpoints.values()),key=lambda x:x["original_hashed_point_id"]),eligible,gate

def index_one(p):
    lat,lon=p["lat"],p["lon"]
    query=urllib.parse.urlencode({"coords":f"POINT({lon:.9f} {lat:.9f})","f":"json"})
    urlpos=BASE+"/comid/position?"+query
    urlhydro=BASE+"/hydrolocation?"+query
    a,aa,ast,ae=api(urlpos)
    time.sleep(SLEEP)
    b,bb,bst,be=api(urlhydro)
    comid_position=get_comid(a)
    hydroid=None
    measured=None
    measure=None
    for feature in features(b):
        props=feature.get("properties") or {}
        geom=feature.get("geometry") or {}
        if props.get("type")!="hydrolocation" or geom.get("type")!="Point":continue
        coords=geom.get("coordinates",[])
        if not isinstance(coords,list) or len(coords)<2:continue
        try:
            x,y=float(coords[0]),float(coords[1])
            hydroid=str(int(props["comid"]))
            measured=haversine_m(lat,lon,y,x)
            measure=props.get("measure")
        except (ValueError,TypeError,KeyError,OverflowError):
            pass
        break
    strict=bool(comid_position is not None and hydroid is not None and
                comid_position==hydroid and measured is not None and measured<=200.)
    sens=bool(hydroid is not None and measured is not None and measured<=500.)
    return {
       "source_hashed_point_id":p["original_hashed_point_id"],
       "original_source_type":p["kind"],
       "original_registered_NPDES_permit_sha256":p["NPDES_permit_ids_sha256"],
       "original_NWCA_state_if_site":p["original_state"],
       "official_USGS_NLDI_catchment_comid":comid_position,
       "official_USGS_NLDI_hydrolocation_comid":hydroid,
       "original_distance_to_hydrolocation_m":round(measured,2) if measured is not None else None,
       "original_network_reach_measure_from_NLDI":measure,
       "USGS_NLDI_catchment_same_as_hydrolocation":comid_position is not None and comid_position==hydroid,
       "source_strict_200m_DIRECT_stream_candidate":strict,
       "source_500m_modelled_sensitivity_only":sens,
       "official_NLDI_position_SHA256":aa,
       "official_NLDI_hydrolocation_SHA256":bb,
       "NLDI_position_HTTP_status":ast,
       "NLDI_hydrolocation_HTTP_status":bst,
       "NLDI_api_errors":[e for e in (ae,be) if e],
       "NLDI_pos_and_hydrolocation_query_URLs_NOT_PUBLISHED_because_include_sensitive_site_coordinates":True,
       "hydrolocation_beyond200m_is_possible_raindrop_MODEL_not_direct_wetland_to_watercourse":True,
    }

def run(root,nwca_dir,echo,shard):
    points,pairs,gate=source_points(root,nwca_dir,echo)
    if len(points)!=157 or len(pairs)<86:raise ValueError("Expected 71 wetland and 86 POTW feature source input count")
    selected=[r for i,r in enumerate(points) if i%SHARDS==shard]
    if len(selected) not in (39,40):raise ValueError("Original 157 source points not split into 4 deterministic shards")
    res=[]
    for i,p in enumerate(selected,1):
        result=index_one(p)
        res.append(result)
        time.sleep(SLEEP)
        if i%10==0 or i==len(selected):
            print("USGS_NLDI_SHARD",shard,"done",i,"/",len(selected),
                  "strict200m",sum(r["source_strict_200m_DIRECT_stream_candidate"] for r in res),
                  "apiErrors",sum(bool(r["NLDI_api_errors"]) for r in res),flush=True)
    return {
      "schema":"gesui_NWCA2021_71wetland_86POTW_actualDMR_USGS_NLDI_comid200m_fourShard_sourceonly_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "source_2026_original_POTW_ECHO_zip_SHA256":sha_file(echo),
      "source_original_exact344_geo_SHA256":sha_file(root/GEO),
      "source_original_2021_DMR160_FY2021FY2022_join_SHA256":sha_file(root/DMR),
      "source_original_exact_2021_reported_NP_source_gate_SHA256":sha_file(root/GATE),
      "source_original_71wetlands_and86POTW_exposure_candidates":True,
      "USGS_NLDI_API_BASE":BASE,
      "original_source_shard_id":shard,
      "total_shards":SHARDS,
      "original_shard_points_n":len(res),
      "source_200m_direct_snap_candidate_n":sum(r["source_strict_200m_DIRECT_stream_candidate"] for r in res),
      "source_500m_sensitivity_modelled_n":sum(r["source_500m_modelled_sensitivity_only"] for r in res),
      "source_original_USGS_point_source_replies":res,
      "all_NWCA2021_plant_H_ALL_FQAI_or_2021_nutrient_measurements_NOT_READ":True,
      "all_2021_actual_DMR_numeric_NP_chemistry_VALUES_not_read_this_stage":True,
      "NHDPlusV2_CATCHMENT_is_not_proof_wetland_in_receiving_river":True,
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nwca-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--shard",type=int,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if a.shard not in range(SHARDS):raise ValueError("Only 4 fixed source shards allowed")
    root=Path(__file__).resolve().parents[1]
    output=run(root,a.nwca_dir,a.npdes_zip,a.shard)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"shard":a.shard,"source_n":output["original_shard_points_n"],
                      "strict_200m":output["source_200m_direct_snap_candidate_n"],
                      "model_500m":output["source_500m_sensitivity_modelled_n"]},indent=2))
if __name__=="__main__":main()
