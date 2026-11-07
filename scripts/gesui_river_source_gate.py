#!/usr/bin/env python3
"""Real source gate: HydroWASTE v1 + European HydroRIVERS directed reaches.

Source-only study design: NO GBIF or biological response is read. An outfall
connected to a river reach is NOT measured N exposure and not an organism corridor.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import zipfile

import geopandas as gpd
import numpy as np
import pandas as pd

SOURCE_HYDROWASTE_BLOB = "0506f873f6541c5f324478894a9731290d722bce"
R_EARTH = 6371.0088
COUNTRIES = ("DEU", "FRA", "GBR", "NLD", "BEL", "AUT", "CHE", "CZE", "POL", "DNK", "SWE")
BBOX = (-12, 44, 26, 64)
BANNED_STATUS = {"closed", "decommissioned", "non-operational", "projected",
                 "proposed", "under construction"}


def git_blob_sha(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(b"blob " + str(len(data)).encode("ascii") + b"\x00" + data).hexdigest()


def haversine(a_lat, a_lon, b_lat, b_lon):
    p1, p2 = math.radians(float(a_lat)), math.radians(float(b_lat))
    dlat = p2 - p1
    dlon = math.radians(float(b_lon) - float(a_lon))
    z = math.sin(dlat / 2)**2 + math.cos(p1)*math.cos(p2)*math.sin(dlon / 2)**2
    return 2*R_EARTH*math.asin(math.sqrt(min(1., max(0., z))))


def fraction_treated(wwtp_m3day, river_m3sec):
    if not np.isfinite(wwtp_m3day) or not np.isfinite(river_m3sec):
        return float("nan")
    if river_m3sec <= 0 or wwtp_m3day < 0:
        return float("nan")
    q = river_m3sec * 86400
    return float(wwtp_m3day / (q + wwtp_m3day))


def extract_shp(source_zip: Path, target: Path) -> Path:
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source_zip) as arc:
        members = arc.namelist()
        if not any(x.lower().endswith(".shp") for x in members):
            raise ValueError("HydroRIVERS source archive has no shapefile")
        arc.extractall(target)
    found = list(target.rglob("*.shp"))
    if len(found) != 1:
        raise ValueError(f"Expected one HydroRIVERS shapefile, found {len(found)}")
    return found[0]


def read_river_subset(shp: Path):
    fields = ["HYRIV_ID", "NEXT_DOWN", "DIS_AV_CMS", "LENGTH_KM", "UPLAND_SKM",
              "HYBAS_L12", "DIST_DN_KM"]
    gdf = gpd.read_file(shp, bbox=BBOX, engine="pyogrio",
                        columns=fields, use_arrow=False)
    if gdf.crs is None:
        raise ValueError("HydroRIVERS has no CRS")
    gdf = gdf.to_crs("EPSG:4326")
    missing = set(fields) - set(gdf.columns)
    if missing:
        raise ValueError(f"HydroRIVERS missing {missing}")
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty]
    gdf["HYRIV_ID"] = pd.to_numeric(gdf.HYRIV_ID, errors="coerce")
    gdf["NEXT_DOWN"] = pd.to_numeric(gdf.NEXT_DOWN, errors="coerce")
    gdf = gdf[gdf.HYRIV_ID.notna()].drop_duplicates("HYRIV_ID")
    rid = {int(r.HYRIV_ID): r for r in gdf.itertuples(index=False)}
    upstream = defaultdict(list)
    for r in gdf.itertuples(index=False):
        if pd.notna(r.NEXT_DOWN) and int(r.NEXT_DOWN) > 0:
            upstream[int(r.NEXT_DOWN)].append(int(r.HYRIV_ID))
    return rid, upstream, len(gdf)


def choose_primary_upstream(reach_id: int, upstream: dict, river: dict):
    opts = [k for k in upstream.get(reach_id, []) if k in river]
    if not opts:
        return None
    def rank(k):
        r = river[k]
        u = float(r.UPLAND_SKM) if pd.notna(r.UPLAND_SKM) else -1.
        q = float(r.DIS_AV_CMS) if pd.notna(r.DIS_AV_CMS) else -1.
        return (u, q, -k)
    return max(opts, key=rank)


def river_chain(waste_rid: int, upstream: dict, river: dict):
    if waste_rid not in river:
        return None
    r = river[waste_rid]
    if pd.isna(r.NEXT_DOWN) or int(r.NEXT_DOWN) not in river:
        return None
    d1 = int(r.NEXT_DOWN)
    u1 = choose_primary_upstream(waste_rid, upstream, river)
    if u1 is None:
        return None
    u2 = choose_primary_upstream(u1, upstream, river)
    if u2 is None or len({u2,u1,waste_rid,d1}) < 4:
        return None
    if (int(river[u2].NEXT_DOWN) != u1
        or int(river[u1].NEXT_DOWN) != waste_rid
        or int(river[waste_rid].NEXT_DOWN) != d1):
        raise AssertionError("Broken directed upstream/downstream source chain")
    return (u2, u1, waste_rid, d1)


def center(r):
    pt = r.geometry.interpolate(0.5, normalized=True)
    return float(pt.y), float(pt.x)


def cohort_candidates(wwtp: pd.DataFrame, river: dict, upstream: dict):
    candidates = []
    chains_valid = 0
    for w in wwtp.itertuples(index=False):
        if pd.isna(w.HYRIV_ID):
            continue
        chain = river_chain(int(w.HYRIV_ID), upstream, river)
        if chain is None:
            continue
        chains_valid += 1
        u2,u1,center_id,d1 = chain
        sites = [center(river[k]) for k in (u2,u1,d1)]
        # 1-km GBIF geometry buffers: centers must be >2.2 km apart
        sep = [haversine(*sites[a], *sites[b]) for a,b in ((0,1),(1,2),(0,2))]
        if min(sep) < 2.2 or max(sep) > 35.0:
            continue
        fraction = fraction_treated(float(w.WASTE_DIS), float(w.RIVER_DIS))
        if not math.isfinite(fraction):
            continue
        candidates.append({
            "wwtp_id": int(w.WASTE_ID), "country":str(w.CNTRY_ISO),
            "river_id":int(center_id), "hydrobasin_l12":str(river[center_id].HYBAS_L12),
            "source_status":str(w.STATUS), "location_quality":int(w.QUAL_LOC),
            "waste_quality":int(w.QUAL_WASTE),
            "outfall_lat":float(w.LAT_OUT), "outfall_lon":float(w.LON_OUT),
            "wwtp_lat":float(w.LAT_WWTP), "wwtp_lon":float(w.LON_WWTP),
            "estimated_ww_m3_day":float(w.WASTE_DIS),
            "estimated_river_m3_sec":float(w.RIVER_DIS),
            "effluent_fraction_proxy":fraction,
            "river_hydraulic_chain":"U2>U1>OUTFALL_REACH>D1",
            "u2_river_id":int(u2), "u1_river_id":int(u1), "d1_river_id":int(d1),
            "u2_lat":sites[0][0], "u2_lon":sites[0][1],
            "u1_lat":sites[1][0], "u1_lon":sites[1][1],
            "d1_lat":sites[2][0], "d1_lon":sites[2][1],
            "distance_u2_u1_km":sep[0], "distance_u1_d1_km":sep[1],
            "distance_u2_d1_km":sep[2],
            "upstream_junction_count":len(upstream.get(center_id,[])),
        })
    return candidates, chains_valid


def select_preoutcome(candidates, max_per_country_stratum=2, min_separation_km=45):
    """Deterministic, source-only high/low fraction strata. No plant records."""
    chosen=[]
    for country in COUNTRIES:
        rows=sorted((x for x in candidates if x["country"]==country),
                    key=lambda x:(x["effluent_fraction_proxy"],x["wwtp_id"]))
        if len(rows) < 8:
            continue
        n=len(rows)
        strata={"low":rows[:max(2,n//4)],"high":rows[-max(2,n//4):]}
        for name in ("high","low"):
            eligible = list(reversed(strata[name])) if name=="high" else list(strata[name])
            collected=0
            for r in eligible:
                if collected>=max_per_country_stratum:
                    break
                if any(haversine(r["outfall_lat"],r["outfall_lon"],
                                 q["outfall_lat"],q["outfall_lon"]) < min_separation_km
                       for q in chosen):
                    continue
                chosen.append(dict(r, exposure_stratum=name))
                collected+=1
    return sorted(chosen,key=lambda x:(x["country"],x["exposure_stratum"],x["wwtp_id"]))


def summarize_source(df: pd.DataFrame, eligible: pd.DataFrame, rivers_n: int,
                     candidate_n: int, chain_valid_n: int, selected: list, zip_sha: str):
    return {
        "schema":"gesui_HydroWASTE_HydroRIVERS_directed_river_source_gate_v0.1",
        "status":"ACTUAL_PUBLIC_SOURCE_SPATIAL_GATE__NO_PLANT_OUTCOME_QUERIED",
        "sources":{
            "hydrowaste_csv_git_blob": SOURCE_HYDROWASTE_BLOB,
            "hydrowaste_original_facilities":int(len(df)),
            "hydrorivers_europe_shp_zip_sha256":zip_sha,
            "hydrorivers_europe_bbox":BBOX,
            "river_reaches_in_bbox":rivers_n,
        },
        "quality":{
            "hydrowaste_facilities_with_linked_hyriv":int(df.HYRIV_ID.notna().sum()),
            "hydrowaste_status_counts":{str(k):int(v) for k,v in df.STATUS.value_counts().to_dict().items()},
            "eligible_facilities_in_western_and_central_europe":int(len(eligible)),
            "topology_chains":chain_valid_n,"sites_meeting_nonoverlap_gate":candidate_n,
            "selected_outfall_units":len(selected),
            "selected_countries":dict(Counter(x["country"] for x in selected)),
            "selected_exposure_strata":dict(Counter(x["exposure_stratum"] for x in selected)),
        },
        "interpretation":"Exposure is potential effluent water fraction at reported/estimated discharge, not measured nitrogen, ecological response, or biologically realized dispersal.",
        "future_biological_gate":"All-Plantae GBIF species facets require independently pre-registered habitat/record-effort checks; no plant data queried by this source gate.",
    }


def run(args):
    if git_blob_sha(args.hydrowaste) != SOURCE_HYDROWASTE_BLOB:
        raise ValueError("HydroWASTE source Git blob mismatch")
    ww=pd.read_csv(args.hydrowaste,encoding="latin1",low_memory=False)
    if len(ww)!=58502 or ww.WASTE_ID.nunique()!=58502:
        raise ValueError("HydroWASTE original facility count / ID mismatch")
    required=["WASTE_ID","HYRIV_ID","LAT_OUT","LON_OUT","LAT_WWTP","LON_WWTP",
              "STATUS","CNTRY_ISO","QUAL_LOC","QUAL_WASTE","WASTE_DIS",
              "RIVER_DIS","COAST_10KM"]
    if set(required)-set(ww.columns):
        raise ValueError(f"HydroWASTE missing: {set(required)-set(ww.columns)}")
    for field in ("HYRIV_ID","QUAL_LOC","QUAL_WASTE","WASTE_DIS","RIVER_DIS","COAST_10KM"):
        ww[field]=pd.to_numeric(ww[field],errors="coerce")
    eligible=ww[
        ww.CNTRY_ISO.isin(COUNTRIES)
        & (~ww.STATUS.astype(str).str.lower().isin(BANNED_STATUS))
        & (ww.QUAL_LOC<=2)
        & (ww.QUAL_WASTE<=3)
        & (ww.WASTE_DIS>0)
        & (ww.RIVER_DIS>=0.5)
        & ww.HYRIV_ID.notna()
        & (ww.COAST_10KM==0)
        & ww.LAT_OUT.between(BBOX[1],BBOX[3])
        & ww.LON_OUT.between(BBOX[0],BBOX[2])
    ].copy()
    zip_sha=hashlib.sha256(args.hydrorivers_zip.read_bytes()).hexdigest()
    shp=extract_shp(args.hydrorivers_zip,args.work_dir)
    river,upstream,rivers_n=read_river_subset(shp)
    cand,valid_n=cohort_candidates(eligible,river,upstream)
    selected=select_preoutcome(cand)
    summary=summarize_source(ww,eligible,rivers_n,len(cand),valid_n,selected,zip_sha)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    (args.out_dir/"source_qc.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n",encoding="utf8")
    (args.out_dir/"candidate_site_design.json").write_text(
        json.dumps({"metadata":{"sample_selected_without_GBIF_response":True,
                                "outfall_count":len(selected),
                                "site_radius_km":1.0,
                                "site_roles":["u2","u1","d1"],
                                "status":"HYPOTHESIS_BLIND_SOURCE_SELECTION"},
                    "selected":selected},indent=2,sort_keys=True)+"\n",encoding="utf8")
    print("GESUI_RIVER_SOURCE_GATE",json.dumps(summary,sort_keys=True),flush=True)
    if not selected:
        raise ValueError("No source-matched directed outfalls; cannot query plants")
    print("GESUI_PREOUTCOME_SELECTION_COUNT",len(selected),flush=True)
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--hydrowaste",type=Path,required=True)
    p.add_argument("--hydrorivers-zip",type=Path,required=True)
    p.add_argument("--work-dir",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    run(p.parse_args())

if __name__=="__main__":
    main()
