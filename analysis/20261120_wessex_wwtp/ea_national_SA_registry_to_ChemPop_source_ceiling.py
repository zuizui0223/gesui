#!/usr/bin/env python3
"""National current EA WQE 5722 Water Company treated-effluent SA sampling points.

Completes source-type-filtered 250/page POST /data/sampling-point (asWKT BNG)
and screens distance to authentic ChemPop native BIO_SITE_ID in 2017/18
quality-eligible source-only subset.

CRITICAL: Current registered SA does NOT mean observed TP in year 2017/18.
No source-to-receptor directed connectivity, outfall geometry, mass, exposure,
biological effect coefficient or 2018 facility activity can be inferred here.
No original site names/IDs/coordinates or individual biology rows are rehosted.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,re,sys,time,urllib.error,urllib.parse,urllib.request
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_SA_TP_to_ChemPop_location_sourcegate as sourcegeo

URL="https://environment.data.gov.uk/water-quality/data/sampling-point"
HEADERS={"User-Agent":"gesui-authentic-EA-WQE-nationwide-SA-source-ceiling-20261120",
         "API-Version":"1","Accept":"application/ld+json","Content-Type":"application/json"}
SOURCE_TYPE="SA"
SOURCE_N_EXPECTED=5722
LIMIT=250
MAX_SOURCE_PAGES=40
MAX_REPLY_BYTES=1400000
REQUIRED_VALID_GEOMETRY_FRACTION=.95
SOURCE_SCREEN_THRESHOLDS=(500,1000,2000,5000,10000)
BNG_WKT_RE=re.compile(r"POINT\s*(?:Z\s*)?\(\s*([-+]?\d+(?:\.\d+)?)\s+([-+]?\d+(?:\.\d+)?)\s*\)",re.I)


def get_page(skip):
    query=urllib.parse.urlencode({"samplingPointType":SOURCE_TYPE,"limit":LIMIT,"skip":skip})
    req=urllib.request.Request(URL+"?"+query,method="POST",data=b"",headers=HEADERS)
    out={"skip":skip,"http":None,"source_sha256":None,"n_source_members":0}
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req,timeout=60) as r:
                out["http"]=int(r.status)
                if out["http"]!=200:raise RuntimeError("Official EA SA registry not HTTP200")
                data=r.read(MAX_REPLY_BYTES+1)
            break
        except urllib.error.HTTPError as exc:
            if exc.code in (502,503,504) and attempt<2:
                time.sleep(1.5*(attempt+1))
                continue
            raise RuntimeError("Official EA SA source registry page HTTP "+str(exc.code)+" at skip "+str(skip)) from exc
    if len(data)>MAX_REPLY_BYTES:
        raise RuntimeError("Official original filtered JSON-LD page exceeded source-only limit")
    out["n_bytes"]=len(data)
    out["source_sha256"]=hashlib.sha256(data).hexdigest()
    doc=json.loads(data)
    members=doc.get("member")
    if not isinstance(members,list):raise RuntimeError("Missing official original JSON-LD member array")
    out["n_source_members"]=len(members)
    out["totalItems"]=doc.get("totalItems")
    return members,out


def as_bng(raw):
    """Read original WKT coordinates. Never assume WGS84 degrees are metres."""
    if not raw:return None,"source_geometry_missing"
    s=str(raw)
    if ("EPSG/0/4326" in s or "CRS84" in s or
        "EPSG:4326" in s):
        return None,"source_explicit_WGS84_requires_proven_conversion"
    m=BNG_WKT_RE.search(s)
    if not m:return None,"source_WKT_unrecognized"
    e,n=float(m.group(1)),float(m.group(2))
    if not math.isfinite(e) or not math.isfinite(n):
        return None,"source_WKT_nonfinite"
    if not sourcegeo.bng_source_valid(e,n):
        return None,"source_WKT_not_plausible_BNG_easting_northing"
    return (e,n),"original_WKT_native_meter_point"


def collect():
    points={}
    status=collections.Counter()
    geofail=collections.Counter()
    region=collections.Counter()
    records=[]
    total_expected=None
    member_total=0
    complete=False
    for page in range(MAX_SOURCE_PAGES):
        mem,receipt=get_page(page*LIMIT)
        records.append(receipt)
        current=receipt["totalItems"]
        if current!=SOURCE_N_EXPECTED:
            raise RuntimeError("Source totalItems drifted from independently verified 5722 SA points")
        if total_expected is None:total_expected=current
        if current!=total_expected:raise RuntimeError("Source count changed mid-paging")
        for v in mem:
            typ=v.get("samplingPointType") or {}
            code=typ.get("notation") if isinstance(typ,dict) else str(typ)
            if str(code)!=SOURCE_TYPE:
                raise RuntimeError("Official API silently returned non-SA source in SA-only query")
            notation=str(v.get("notation","")).strip()
            if not notation:
                raise RuntimeError("Official SA sampling-point notation missing")
            if notation in points:
                raise RuntimeError("Original nationwide SA notation duplicate across source pages")
            shape=v.get("geometry") or {}
            if not isinstance(shape,dict):shape={}
            coord,why=as_bng(shape.get("asWKT"))
            if coord is None:geofail[why]+=1
            stat=v.get("samplingPointStatus") or {}
            status[str(stat.get("notation") if isinstance(stat,dict) else stat)]+=1
            reg=v.get("region") or {}
            region[str(reg.get("notation") if isinstance(reg,dict) else reg)]+=1
            points[notation]=coord
            member_total+=1
        if len(mem)<LIMIT:
            complete=True;break
    if not complete or member_total!=SOURCE_N_EXPECTED or len(points)!=SOURCE_N_EXPECTED:
        raise RuntimeError("National official SA registry coverage incomplete: "+str((member_total,len(points),complete)))
    valid={s:coords for s,coords in points.items() if coords is not None}
    frac=len(valid)/len(points)
    qa={
       "official_data_nationwide_SA_totalItems":int(total_expected),
       "complete_source_unique_native_SA_sampling_point_count":len(points),
       "source_pages_used":len(records),
       "source_page_SHA256":[x["source_sha256"] for x in records],
       "source_type":"SA Water Company FINAL/TREATED effluent monitoring point",
       "source_geo_valid_BNG_point_n":len(valid),
       "source_geo_valid_fraction":frac,
       "source_missing_unparsable_geometry_counts":dict(geofail),
       "original_registry_status_notations":dict(status),
       "original_region_notation_counts":dict(region),
       "source_all_station_code_SHA256":hashlib.sha256(
            ("\n".join(sorted(points))).encode()).hexdigest(),
       "source_valid_BNG_station_code_SHA256":hashlib.sha256(
            ("\n".join(sorted(valid))).encode()).hexdigest(),
       "published_registry_is_2026_current_not_historical_2017_2018_sampling":True,
       "no_permit_or_outfall_massing_or_routing":True
    }
    if frac<REQUIRED_VALID_GEOMETRY_FRACTION:
        raise RuntimeError("National source BNG coverage <95%; SOURCE STOP; QA="+json.dumps(qa,sort_keys=True))
    return valid,qa


def compare_bio(stw,biosites):
    if not stw or not biosites:
        return {"source_gate_status":"NO_BIO_OR_SA_SOURCE","n_biosites":len(biosites),"n_SA":len(stw)}
    ids=sorted(biosites)
    keys=sorted(stw)
    coords=np.asarray([biosites[k] for k in ids],float)
    src=np.asarray([stw[k] for k in keys],float)
    mindist=np.full(len(ids),np.inf)
    within_site={th:np.zeros(len(ids),bool) for th in SOURCE_SCREEN_THRESHOLDS}
    within_source={th:np.zeros(len(keys),bool) for th in SOURCE_SCREEN_THRESHOLDS}
    paircount=collections.Counter()
    # Streaming original site geometry in blocks avoids >100MiB distance matrix.
    for lo in range(0,len(ids),200):
        hi=min(len(ids),lo+200)
        delta=coords[lo:hi,None,:]-src[None,:,:]
        distances=np.sqrt(np.sum(delta*delta,axis=2))
        mindist[lo:hi]=distances.min(axis=1)
        for th in SOURCE_SCREEN_THRESHOLDS:
            mask=distances<=th
            within_site[th][lo:hi]=mask.any(axis=1)
            within_source[th] |= mask.any(axis=0)
            paircount[th]+=int(mask.sum())
    return {
        "n_BIO_sites":len(ids),
        "n_WQE_current_registry_SA_sites":len(keys),
        "nearest_current_SA_point_BNG_distance_m_p10_p50_p90":[float(x) for x in np.quantile(mindist,[.1,.5,.9])],
        "proximity":{"%sm"%th:{
          "n_original_BIO_site_with_any_current_registered_SA_point_within":int(within_site[th].sum()),
          "n_SA_point_with_any_BIO_site_within":int(within_source[th].sum()),
          "candidate_geographic_point_pairs":int(paircount[th])
        } for th in SOURCE_SCREEN_THRESHOLDS},
        "not_historical_STW_operation_or_downstream_receiving_link":True
    }


def main(a):
    result={
        "schema":"GESUI_20261120_NATIONWIDE_SA_WQE_CURRENT_REGISTRY_CHEMPOP_BNG_SOURCE_CEILING_V1",
        "source_only_no_biological_outcome_fitted":True,
        "current_SA_samplingpoint_is_NOT_observed_2017_18_TP_concentration":True,
        "does_not_establish_STW_flow_direction_to_BIO":True,
        "does_not_calculate_STW_load_kg_per_year":True,
        "status":"SOURCE_STOP",
    }
    try:
        national,regqa=collect()
        bio,year_elig,bioqa=sourcegeo.chempop_geo()
        result["official_national_SA_current_registry"]=regqa
        result["original_ChemPop_BNG_site_source"]=bioqa
        result["all_original_BIO_sites_current_registry_screen"]=compare_bio(national,bio)
        result["2017_2018_original_biology_years_current_registry_screen"]={
           str(yr):compare_bio(national,{k:bio[k] for k in year_elig[yr]}) for yr in (2017,2018)}
        union=year_elig[2017]|year_elig[2018]
        result["any_original_BIO_biology_in_2017_18_current_registry_screen"]=compare_bio(
            national,{k:bio[k] for k in union})
        result["source_gates"]={
            "complete_national_sampling_point_SA_registry_verified":True,
            "geometry_95percent_quality_GATE_PASS":True,
            "historic_tp_measured_for_each_current_source_SA_point":False,
            "true_year_specifically_active_STW_and_BIO_source_record_concordance":False,
            "directed_original_river_STW_to_BIO_connection_verified":False,
            "annual_STW_N_or_P_mass_verified":False,
            "source_is_geographical_opportunity_ceiling_only":True,
        }
        result["status"]="NATIONAL_CURRENT_SA_REGISTRY_CHEMPOP_BNG_SOURCE_SCREEN_COMPLETED"
    except Exception as e:
        result.update({"status":"SOURCE_REGISTRY_OR_BNG_GATE_STOP",
                       "error":type(e).__name__+": "+str(e)[:650]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+"\n")
    summary={
       "status":result["status"],"error":result.get("error"),
       "SA_registry":{k:v for k,v in result.get("official_national_SA_current_registry",{}).items()
                       if k not in ("source_page_SHA256",)},
       "original_BIO":result.get("original_ChemPop_BNG_site_source"),
       "nationwide_all_BIO":result.get("all_original_BIO_sites_current_registry_screen"),
       "bio_2017_18":result.get("2017_2018_original_biology_years_current_registry_screen"),
       "bio_either_2017_18":result.get("any_original_BIO_biology_in_2017_18_current_registry_screen"),
       "source_gates":result.get("source_gates")
    }
    print("GESUI_NATIONAL_WQE_SA_CHEMPOP_SOURCE_CEILING="+json.dumps(
        summary,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
