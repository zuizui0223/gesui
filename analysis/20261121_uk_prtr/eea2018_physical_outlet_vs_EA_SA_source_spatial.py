#!/usr/bin/env python3
"""Source-only geographic crosswalk 2018 UK EEA UWWTD historic outlets to EA SA.

2018 archived EEA UWWTD has ACTUAL reported physical discharge points:
UK5 1799 urban WWTP plant rows; 1809 reported discharge point rows,
1808 plausible UK WGS84, 1498 FW freshwater outlet classifications,
1798 native plant IDs with >=1 outfall, 1787 with exactly 1 outlet.
Original UK source positive M-assessed N/P loads in original UWWTD table:
ZERO (do not replace with PRTR facility reports or WQE TP mg/L as load).

Join original 2018 UK UWWTD native repCode/uwwCode/dcpCode to its physical
outlet, then source-only nearest geometry screening vs *modern* 2026 EA
WQE water company FINAL treated SA sampling-points. No raw names, IDs,
coordinates or mass values are published. Proximity is NOT verified
facility identity. SA water company monitoring point is not the
reported physical outlet; original 2018 location is not confirmed active
in 2017 or 2026. Do NOT fit biodiversity outcome or P-load model.
"""
from __future__ import annotations

import argparse,collections,csv,hashlib,io,json,math,sys,tempfile,urllib.request,zipfile
from pathlib import Path
import numpy as np
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as sa
sys.path.insert(0,str(ROOT/"analysis"/"20261121_uk_prtr"))
import eea_UWWTD2018_UK_outlet_source_gate as eea

EXPECTED_UK_REPORTER="UK5"
EXPECTED_UK_PLANTS=1799
EXPECTED_UK_OUTLETS=1809
EXPECTED_FRESHWATER_FW=1498
EXPECTED_N_EA_SA=5722
THRESHOLDS_M=(100,250,500,1000,2000)
BBOX=(49,61,-12,4)

def exact_source():
    original={}
    with tempfile.TemporaryDirectory(prefix="gesui_UWWTD2018_geosource_") as temp:
        archive=Path(temp)/"EEA_original.zip"
        req=urllib.request.Request(eea.URL,headers={"User-Agent":"gesui-original-UWWTD2018-UK-outlet-to-SA-source-only"})
        sha=hashlib.sha256();sz=0
        with urllib.request.urlopen(req,timeout=100) as r,archive.open("wb") as f:
            if r.status!=200:raise RuntimeError("Official EEA UWWTD source not HTTP200")
            while True:
                data=r.read(131072)
                if not data:break
                sz+=len(data)
                if sz>eea.MAX_BYTES:raise RuntimeError("Official EEA UK original archive exceeds bound")
                sha.update(data);f.write(data)
        with zipfile.ZipFile(archive) as z:
            for name,digest in eea.EXPECT.items():
                b=z.read(eea.BASE+name)
                if hashlib.sha256(b).hexdigest()!=digest:
                    raise RuntimeError("Original EEA UK 2018 source member SHA not independently frozen")
                original[name]=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig","replace"))))
    periods=original["ReportPeriod.csv"]
    uk=[r for r in periods if r.get("repCode","").strip()==EXPECTED_UK_REPORTER
        and r.get("repReportedPeriod","").strip()=="2018"
        and "2018" in r.get("repSituationAt","")]
    if len(uk)!=1:
        raise RuntimeError("Original UK5 2018 source reporter not uniquely verified")
    plants=[r for r in original["UWWTPS.csv"] if r.get("repCode","").strip()==EXPECTED_UK_REPORTER]
    outlets=[r for r in original["DischargePoints.csv"] if r.get("repCode","").strip()==EXPECTED_UK_REPORTER]
    if len(plants)!=EXPECTED_UK_PLANTS or len(outlets)!=EXPECTED_UK_OUTLETS:
        raise RuntimeError("Original EEA 2018 1799 UK plants/1809 discharge points unexpectedly changed")
    ids=[r.get("uwwCode","").strip() for r in plants]
    if not all(ids) or len(set(ids))!=len(ids):
        raise RuntimeError("Original reported UK UWWTD plant native ID invalid/duplicate")
    outf=collections.Counter(r.get("uwwCode","").strip() for r in outlets)
    if sum(outf[k]==1 for k in ids)!=1787:
        raise RuntimeError("Original UK UWWTD one-outlet plant source drift")
    fw=[r for r in outlets if r.get("dcpWaterBodyType","").strip().upper()=="FW"]
    if len(fw)!=EXPECTED_FRESHWATER_FW:
        raise RuntimeError("UK 2018 1498 UWWTD freshwater outlet source count drift")
    # The full 1498 include plants having 2+ outlets or not source-uniquely
    # assignable. Report them as candidates, not independent plant exposure.
    coords={}
    bad=0;dupes=0
    for r in fw:
        lat,lon=eea.finite(r.get("dcpLatitude")),eea.finite(r.get("dcpLongitude"))
        if lat is None or lon is None or not(49<=lat<=61 and -12<=lon<=4):
            bad+=1;continue
        native=(r["repCode"].strip(),r["uwwCode"].strip(),r["dcpCode"].strip())
        if any(not v for v in native):
            raise RuntimeError("Original UK FW outlet lacks required native source code")
        if native in coords:
            dupes+=1
            if coords[native]!=(lon,lat):
                raise RuntimeError("Source outlet native ID duplicated with conflicting WGS84 coordinates")
        else:coords[native]=(lon,lat)
    result={
      "official_2018_EEA_archive_sha256":sha.hexdigest(),
      "official_2018_EEA_archive_size":sz,
      "exact_three_2018_source_member_SHA256":eea.EXPECT,
      "UK_reporter_source":"UK5, original reporting year 2018",
      "original_2018_UK_UWWTD_plant_rows":len(plants),
      "original_2018_UK_all_outlet_rows":len(outlets),
      "original_2018_UK_FRESHWATER_outlet_rows":len(fw),
      "original_2018_UK_FW_valid_distinct_native_outlet_geographic_points":len(coords),
      "original_2018_UK_FW_outlet_coordinate_source_missing":bad,
      "original_2018_UK_FW_native_outlet_ID_duplicates":dupes,
      "original_2018_UK_single_discharge_point_plants":sum(outf[k]==1 for k in ids),
      "original_2018_UK_multiple_discharge_point_plants":sum(outf[k]>1 for k in ids),
      "origin_is_2018_reported_discharge_point_not_STW_sample":True,
      "not_measured_annual_P_or_N_source":True,
    }
    return coords,outf,result

def distance_screen(uwwtd,sa_points,one_outlet):
    uk_ids=sorted(uwwtd)
    sa_ids=sorted(sa_points)
    b=np.asarray([sa_points[k] for k in sa_ids],float)
    ut=np.asarray([uwwtd[k] for k in uk_ids],float)
    if not len(uk_ids) or not len(sa_ids):
        raise RuntimeError("Source original outlet and SA geometries empty")
    # 1500x5700 distances <1e7 float values; blockwise to bound memory.
    counts={r:{"n_outlet_with_at_least_one_SA":0,
                "n_outlet_with_EXACTLY_ONE_SA":0,
                "n_outlet_with_2plus_SA_ambiguous":0,
                "n_pairs":0,"n_single_discharge_plant_outlet_with_one_SA":0,
                "SA_source_ids":set()}
            for r in THRESHOLDS_M}
    nearest=[]
    for start in range(0,len(uk_ids),80):
        end=min(len(uk_ids),start+80)
        dist=np.linalg.norm(ut[start:end,None,:]-b[None,:,:],axis=2)
        nearest.extend(dist.min(axis=1).tolist())
        for radius,c in counts.items():
            eligible=dist<=radius
            c["n_outlet_with_at_least_one_SA"]+=int(eligible.any(axis=1).sum())
            c["n_outlet_with_EXACTLY_ONE_SA"]+=int((eligible.sum(axis=1)==1).sum())
            c["n_outlet_with_2plus_SA_ambiguous"]+=int((eligible.sum(axis=1)>1).sum())
            c["n_pairs"]+=int(eligible.sum())
            for i,row in enumerate(eligible):
                if int(row.sum())==1 and one_outlet[uk_ids[start+i][1]]==1:
                    c["n_single_discharge_plant_outlet_with_one_SA"]+=1
                for j in np.flatnonzero(row):
                    c["SA_source_ids"].add(sa_ids[int(j)])
    out={}
    for radius,c in counts.items():
        out[str(radius)+"m"]={
            k:v for k,v in c.items() if k!="SA_source_ids"
        }
        out[str(radius)+"m"]["unique_current_SA_points_with_2018_outlet_geographic_candidate"]=len(c["SA_source_ids"])
    return {"n_historic_2018_UK_FW_reported_outlet_points":len(uk_ids),
            "n_current_EA_WQE_WaterCompany_SA_points":len(sa_ids),
            "nearest_SA_distance_m_p10_median_p90":[float(x) for x in np.quantile(nearest,[.1,.5,.9])],
            "original_outlet_to_2026_sampling_point_candidate_thresholds":out,
            "historical_2018_outlet_to_current_2026_SA_is_NOT_proven_same_facility_or_year_activity":True}

def main(a):
    result={"schema":"GESUI_20261121_UWWTD2018_UK_HISTORICAL_DISCHARGEPOINT_VS_CURRENT_WQE_SA_SOURCE_CROSSWALK_V1",
            "source_only":True,"status":"SOURCE_STOP",
            "no_facility_PN_kg_year_or_macroinvertebrate_BMWP_fit":True,
            "no_original_native_ID_or_coordinate_republished":True,
            "not_verified_permit_to_EA_SA_point_same_facility":True}
    try:
        points,origin_counts,qa=exact_source()
        current,currentqa=sa.collect()
        if len(current)!=EXPECTED_N_EA_SA:
            raise RuntimeError("Current EA 5722 strict SA source registry changed")
        transform=Transformer.from_crs("EPSG:4326","EPSG:27700",always_xy=True)
        bng={}
        for key,(lon,lat) in points.items():
            e,n=transform.transform(lon,lat)
            if not(0<e<800000 and 0<n<1300000):
                raise RuntimeError("Valid EEA UK discharge point transformed outside UK BNG")
            bng[key]=(e,n)
        result.update({
          "status":"EEA_UWWTD2018_UK_FW_OUTLETS_VS_EA_WQE_SA_GEOGRAPHIC_SOURCE_VERIFIED",
          "original_2018_UWWTD_UK_official_source":qa,
          "official_EA_2026_strict_SA_registry_n":len(current),
          "official_EA_SA_source_SHA256":currentqa["source_all_station_code_SHA256"],
          "UK_2018_UWWTD_FW_outlet_to_current_EA_SA_geography":distance_screen(
              bng,current,origin_counts),
          "source_gates":{
            "2018_original_UK_UWWTD_reported_physical_discharge_locations_exist":len(points)>0,
            "2018_UK_UWWTD_outlet_HAS_verified_EA_SA_native_identity":False,
            "2018_UK_UWWTD_UWWTP_reports_has_joint_measured_positive_PN":False,
            "historical_outlet_source_location_not_registered_monitoring_station_sample":True,
            "source_PN_annual_mass_kg_from_EA_WQE_TP_mgL_NOT_authorized":True,
            "causal_ChemPop_benthic_response_regression_NOT_authorized":True,
        }})
    except Exception as exc:
        result.update({"status":"SOURCE_UWWTD_OUTLET_SA_CROSSWALK_STOP",
                       "error":type(exc).__name__+": "+str(exc)[:650]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_UWWTD2018_UK_FW_OUTLET_VS_WQE_SA_GEO="+json.dumps({
        "status":result["status"],"error":result.get("error"),
        "2018_UK_source":result.get("original_2018_UWWTD_UK_official_source"),
        "outlet_SA_geometry":result.get("UK_2018_UWWTD_FW_outlet_to_current_EA_SA_geography"),
        "gates":result.get("source_gates")
    },sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
