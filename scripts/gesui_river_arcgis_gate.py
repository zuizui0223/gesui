#!/usr/bin/env python3
"""Provisional HydroRIVERS global ArcGIS mirror for official download HTTP 403.
This mirror is not an authenticated official HydroRIVERS release.
"""
from __future__ import annotations
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from types import SimpleNamespace
import argparse
import hashlib
import json
import math
import time

import numpy as np
import pandas as pd
from shapely.geometry import shape
from gesui_river_source_gate import (
    BBOX, BANNED_STATUS, COUNTRIES, SOURCE_HYDROWASTE_BLOB, git_blob_sha,
    cohort_candidates, fraction_treated, select_preoutcome, haversine,
)

SERVICE = "https://maps.fsc.org/server/rest/services/Hosted/Optimized_Hyrdo/FeatureServer/0/query"
FIELDS = "hyriv_id,next_down,upland_skm,dis_av_cms,length_km,dist_dn_km,hybas_l12"


def chunks(ids, n=35):
    a=sorted(set(int(x) for x in ids if x is not None and int(x)>0))
    for i in range(0,len(a),n):
        yield a[i:i+n]


def request(where, receipts):
    url=SERVICE+"?"+urlencode({
        "f":"geojson","where":where,"outFields":FIELDS,"returnGeometry":"true",
        "outSR":"4326","resultRecordCount":"2000"
    })
    err=None
    for i in range(4):
        try:
            req=Request(url,headers={"User-Agent":"gesui-ecology-source-probe/0.1",
                                     "Accept":"application/json"})
            with urlopen(req,timeout=70) as resp:
                raw=resp.read()
            data=json.loads(raw)
            if "error" in data:
                raise ValueError("ArcGIS response error "+str(data["error"])[:200])
            if "features" not in data or data.get("exceededTransferLimit"):
                raise ValueError("ArcGIS response missing features or truncated")
            receipts.append({"where":where,"n":len(data["features"]),
                             "sha256":hashlib.sha256(raw).hexdigest(),"bytes":len(raw)})
            return data["features"]
        except (HTTPError,URLError,TimeoutError,ValueError) as e:
            err=e
            time.sleep(2*(i+1))
    raise RuntimeError("Mirror request failed: "+str(err))


def query(ids, selector, receipts):
    ret=[]
    for group in chunks(ids):
        ret.extend(request(selector+" IN ("+",".join(str(x) for x in group)+")",receipts))
    return ret


def read_float(x):
    try:return float(x)
    except (TypeError,ValueError):return float("nan")


def ingest(features, river, up):
    for item in features:
        props={k.lower():v for k,v in (item.get("properties") or {}).items()}
        if not props.get("hyriv_id") or not item.get("geometry"):
            continue
        rid=int(props["hyriv_id"])
        down=int(props.get("next_down") or 0)
        obj=SimpleNamespace(
            HYRIV_ID=rid,NEXT_DOWN=down,
            UPLAND_SKM=read_float(props.get("upland_skm")),
            DIS_AV_CMS=read_float(props.get("dis_av_cms")),
            LENGTH_KM=read_float(props.get("length_km")),
            DIST_DN_KM=read_float(props.get("dist_dn_km")),
            HYBAS_L12=str(props.get("hybas_l12") or ""),
            geometry=shape(item["geometry"]),
        )
        if rid not in river:
            river[rid]=obj
            if down>0:
                up[down].append(rid)
        elif river[rid].NEXT_DOWN!=down:
            raise ValueError("Mirror returned conflicting river topology")


def preselect(rows, per_country_stratum=20):
    res=[]
    for country in COUNTRIES:
        x=rows.loc[rows.CNTRY_ISO==country].copy()
        if len(x)<8: continue
        x["water_fraction"]=x.apply(lambda r:fraction_treated(
            float(r.WASTE_DIS),float(r.RIVER_DIS)),axis=1)
        x=x[np.isfinite(x.water_fraction)].sort_values(["water_fraction","WASTE_ID"])
        if len(x)<8:continue
        q=max(2,len(x)//4)
        for t in (x.iloc[-q:].iloc[::-1],x.iloc[:q]):
            n=0
            for row in t.itertuples(index=False):
                if n>=per_country_stratum:break
                if any(haversine(float(row.LAT_OUT),float(row.LON_OUT),
                                 float(r.LAT_OUT),float(r.LON_OUT))<12 for r in res):
                    continue
                res.append(row)
                n+=1
    return res


def run(hydro,out_dir):
    if git_blob_sha(hydro)!=SOURCE_HYDROWASTE_BLOB:
        raise ValueError("HydroWASTE raw source blob has changed")
    ww=pd.read_csv(hydro,encoding="latin1",low_memory=False)
    if len(ww)!=58502 or ww.WASTE_ID.nunique()!=58502:
        raise ValueError("Original HydroWASTE facility count invalid")
    for c in ("QUAL_LOC","QUAL_WASTE","WASTE_DIS","RIVER_DIS","COAST_10KM","HYRIV_ID"):
        ww[c]=pd.to_numeric(ww[c],errors="coerce")
    eligible=ww[
        ww.CNTRY_ISO.isin(COUNTRIES)
        & ~ww.STATUS.astype(str).str.lower().isin(BANNED_STATUS)
        & (ww.QUAL_LOC<=2) & (ww.QUAL_WASTE<=3)
        & (ww.WASTE_DIS>0) & (ww.RIVER_DIS>=0.5)
        & (ww.COAST_10KM==0) & ww.HYRIV_ID.notna()
        & ww.LAT_OUT.between(BBOX[1],BBOX[3])
        & ww.LON_OUT.between(BBOX[0],BBOX[2])
    ].copy()
    candidate=preselect(eligible)
    if not candidate: raise ValueError("No eligible preselected WWTP")
    river={}
    up=defaultdict(list)
    receipts=[]
    seed=[int(x.HYRIV_ID) for x in candidate]
    ingest(query(seed,"hyriv_id",receipts),river,up)
    matched=[k for k in seed if k in river]
    if not matched:raise ValueError("Mirror has no HydroWASTE matching reaches")
    downstream=[river[k].NEXT_DOWN for k in matched]
    ingest(query(downstream,"hyriv_id",receipts),river,up)
    ingest(query(matched,"next_down",receipts),river,up)
    first_up=[x for k in matched for x in up.get(k,[])]
    ingest(query(first_up,"next_down",receipts),river,up)
    df=pd.DataFrame([r._asdict() for r in candidate])
    full,topological=cohort_candidates(df,river,up)
    selected=select_preoutcome(full,max_per_country_stratum=2,min_separation_km=45)
    out_dir.mkdir(parents=True,exist_ok=True)
    qc={
        "schema":"gesui_pinned_HydroWASTE_global_HydroRIVERS_mirror_source_gate_v01",
        "source_status":"SOURCE_GEOMETRY_ONLY_NO_PLANTS",
        "original_hydrowaste_records":len(ww),"hydrowaste_git_blob":SOURCE_HYDROWASTE_BLOB,
        "eligible_europe_wwtp":len(eligible),"preselected_facilities":len(candidate),
        "mirror_url":SERVICE,"mirror_official_version_verified":False,
        "query_receipts":receipts,"mirror_reaches":len(river),
        "mirrored_outfall_reaches":len(matched),
        "valid_four_reach_chains":topological,
        "nonoverlap_candidate_chains":len(full),
        "selected_outfalls":len(selected),
        "selected_countries":dict(Counter(x["country"] for x in selected)),
        "selected_strata":dict(Counter(x["exposure_stratum"] for x in selected)),
        "warnings":["Mirror external and not original SHA-verified official HydroRIVERS",
                    "Flow fraction potential, not N/P mass or dose",
                    "No plant records queried and no causal analysis",
                    "Other upstream wastewater facilities are not yet mapped"],
    }
    (out_dir/"source_qc.json").write_text(json.dumps(qc,indent=2,sort_keys=True)+"\n")
    (out_dir/"candidate_site_design.json").write_text(json.dumps({
        "metadata":{"source_preselected_before_plant_outcomes":True,"site_buffer_km":1,
                    "sampling_sites":["u2","u1","d1"]},
        "selected":selected},indent=2,sort_keys=True)+"\n")
    print("GESUI_MIRRORED_RIVER_SOURCE_GATE",
          json.dumps({k:v for k,v in qc.items() if k!="query_receipts"},sort_keys=True),flush=True)
    if not selected:raise ValueError("No geometry-certified sites; do not query plants")
    return qc


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--hydrowaste",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    run(a.hydrowaste,a.out_dir)

if __name__=="__main__":
    main()
