#!/usr/bin/env python3
"""Source-only original PRTR urban STW facility vs EA WQE final SA station BNG.

The UK PRTR 2017/18 original TRUE ReportingYear (not CSV filename)
contains annual water release total N, some total P at large urban WWTPs,
but has NO same-year same-facility paired P/N in source QA. Check if an
original England Environment Agency regulated PRTR FACILITY POINT is near
one or more official 2026-current WQE SA Water Company treated-effluent
SAMPLING POINTS. This is a GEOGRAPHIC identity candidate only; neither
point establishes a physically permitted OUTFALL location or river
chemistry delivered downstream. NOT a source-to-BMWP response model.

Official PRTR releases = FACILITY/PERMIT mass kg/year, mixed methods
M/C/E and reporting thresholds. UK PRTR covers all nations; source
subcohort requires original CompetentAuthorityPartyName=="EA", 5.(f)
and NACE37.00, original WATER and original TOTAL - NITROGEN/PHOSPHORUS.
Native NationalID is a facility identifier, and NationalID from
different PRTR wrappers may be updated; no constructed facility matching
is used. EA WQE SA 2026 registry is NOT actual 2017/18 observations.

No raw coordinates, names, native IDs or facility mass quantities emitted.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,sys,tempfile
from pathlib import Path
import numpy as np
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261121_uk_prtr"))
import source_only_2017_2018_prtr_schema as original
import census_true_year_WWTW_PN_2017_2018 as prtr
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as sa

EXPECTED_SHA={
 "2017":"4a8dad650ea37d4a654aa06d2fc3c9c3003062d2150ce4c40009f0d790858820",
 "2018":"2eaf3ba3f679e5de8b56e3ddc0723c9baa179909e700751fbfb12f8a619a57b5",
 "2019":"6dc05a80d7149f0dab959b6605432a6133833a869899f7b9b85b8071666d507a",
}
EXPECTED_EA_STRICT={
 "2017":{"TOTAL - NITROGEN":103,"TOTAL - PHOSPHORUS":0},
 "2018":{"TOTAL - NITROGEN":98,"TOTAL - PHOSPHORUS":8}
}
D_THRESHOLDS_M=(100,250,500,1000,2000)
MAX_SOURCE_POINTS=5722

def matches():
    original_rows=collections.defaultdict(list)
    source_reports={}
    with tempfile.TemporaryDirectory(prefix="gesui_prtr_facility_SA_") as temp:
        for label,url in prtr.FILES.items():
            path=Path(temp)/("prtr_"+label+".csv")
            receipt=original.download(url,path)
            if receipt["status"]!="OFFICIAL_ORIGINAL_CSV_BYTES_OBTAINED":
                raise RuntimeError("Official PRTR wrapper unavailable "+label)
            if receipt["source_sha256"]!=EXPECTED_SHA[label]:
                raise RuntimeError("Original PRTR CSV SHA changed "+label)
            source_reports[label]={"SHA256":receipt["source_sha256"],
                                   "bytes":receipt["source_bytes"]}
            for row in prtr.readcsv(path):
                yr=row["ReportingYear"].strip()
                if (yr not in prtr.VALID_YEARS
                    or row["MediumCode"].strip()!=prtr.MEDIUM
                    or row["PollutantCode"].strip() not in prtr.TARGET_POLLUTANTS
                    or row["AnnexIActivityCode"].strip()!=prtr.WWTP_ACTIVITY
                    or row["NACEMainEconomicActivityCode"].strip()!=prtr.WWTP_NACE
                    or row["CompetentAuthorityPartyName"].strip()!="EA"):
                    continue
                val=prtr.numeric(row["TotalQuantity"])
                if val is None or val<=0:continue
                natid=row["NationalID"].strip()
                if not natid:raise RuntimeError("Strict original PRTR selected NationalID empty")
                lat=prtr.numeric(row["LatitudeMeasure"])
                lon=prtr.numeric(row["LongitudeMeasure"])
                original_rows[(yr,row["PollutantCode"].strip(),natid)].append(
                    {"lat":lat,"lon":lon,"method":row["MethodBasisCode"].strip(),
                     "TotalQuantity":val,"source_file":label})
    # A single facility may have multiple original designations: NEVER
    # integrate quantities. A facility POINT can be admitted only when
    # the reported original coordinates are consistent across release rows.
    points={}
    sourceqa=collections.Counter()
    basis={}
    for key,observations in original_rows.items():
        yr,sub,facid=key
        coords=set((v["lat"],v["lon"]) for v in observations)
        if len(coords)!=1:
            sourceqa["facility_year_pollutant_ambiguous_original_coords"]+=1
            continue
        lat,lon=next(iter(coords))
        if lat is None or lon is None or not (49<=lat<=57 and -9<=lon<=4):
            sourceqa["missing_or_invalid_original_UK_facility_WGS84"]+=1
            continue
        kinds=set((v["method"],round(v["TotalQuantity"],8)) for v in observations)
        if len(kinds)>1:
            sourceqa["same_facility_pollutant_year_multiple_amount_or_method_reports"]+=1
        points[key]=(lon,lat)
        basis[key]=tuple(sorted(set(v["method"] for v in observations)))
    for yr in prtr.VALID_YEARS:
        for p in prtr.TARGET_POLLUTANTS:
            n=len({nid for yy,pp,nid in points if yy==yr and pp==p})
            if n!=EXPECTED_EA_STRICT[yr][p]:
                raise RuntimeError(f"Historical EA genuine {yr} {p} original source {n} differs from earlier independently audited {EXPECTED_EA_STRICT[yr][p]}")
    return points,basis,{
        "original_strict_facility_year_pollutant_n":len(points),
        "source_duplicate_or_invalid_geometry_QA":dict(sourceqa),
        "original_2017_18_2019_official_files":source_reports,
        "source_facility_conforms_5f_NACE37_and_EA":True,
    }

def score_geo(easting_by_fac,sa_points):
    wqe_keys=sorted(sa_points)
    if not wqe_keys:raise RuntimeError("EA official native SA points unavailable")
    s=np.asarray([sa_points[k] for k in wqe_keys],float)
    out={}
    for year in ("2017","2018"):
        peryear={}
        for pol in prtr.TARGET_POLLUTANTS:
            pts={k:v for k,v in easting_by_fac.items() if k[0]==year and k[1]==pol}
            n=len(pts)
            frame=np.asarray([pts[k] for k in sorted(pts)],float).reshape((-1,2))
            if n:
                dist=np.linalg.norm(frame[:,None,:]-s[None,:,:],axis=2)
                nearest=dist.min(axis=1)
            else:
                dist=np.empty((0,len(s)),float)
                nearest=np.empty(0,float)
            thresholds={}
            for radius in D_THRESHOLDS_M:
                m=dist<=radius
                thresholds[str(radius)+"m"]={
                    "source_PRTR_facility_points_with_geographically_near_SA":int(m.any(axis=1).sum()),
                    "source_PRTR_facility_points_with_EXACTLY_ONE_geographically_near_SA":int((m.sum(axis=1)==1).sum()),
                    "source_PRTR_facility_points_with_2_or_more_SA_points_ambiguous":int((m.sum(axis=1)>=2).sum()),
                    "EA_current_SA_sampling_points_with_any_PRTR_facility_point_within":int(m.any(axis=0).sum()),
                    "gross_geometric_facility_sampling_point_pairs":int(m.sum()),
                }
            peryear[pol]={
                "original_true_reporting_year":year,
                "n_strict_urban_WWTW_EA_facility_PN_water_mass_reports_with_geography":n,
                "PRTR_facility_to_current_WQE_SA_nearest_distance_median":float(np.median(nearest)) if n else None,
                "nearest_distance_m_p10_p90":[float(v) for v in np.quantile(nearest,[.1,.9])] if n else None,
                "distance_source_candidates_not_confirmed_same_permit_or_outfall":True,
                "thresholds":thresholds,
            }
        out[year]=peryear
    return out

def main(args):
    result={
        "schema":"GESUI_20261121_PRTR_URBAN_WWTW_FACILITY_VS_WQE_TREATED_SA_GEOSOURCE_V1",
        "source_only_no_ChemPop_BMWP_or_biological_slopes":True,
        "PRTR_facility_latlon_not_physical_permitted_river_outfall":True,
        "EA_SA_2026_monitoringpoint_not_2017_18_activity_or_outfall":True,
        "PRTR_reported_annual_kgP_kgN_source_not_routed_mass":True,
        "no_republication_of_exact_facility_mass_coordinates_or_source_native_ids":True,
        "status":"SOURCE_STOP",
    }
    try:
        original_pts,basis,prtrqa=matches()
        sa_points,saqa=sa.collect()
        if len(sa_points)!=MAX_SOURCE_POINTS:raise RuntimeError("Original WQE SA current 5722 native registry count changed")
        transformer=Transformer.from_crs("EPSG:4326","EPSG:27700",always_xy=True)
        bng={}
        for key,(lon,lat) in original_pts.items():
            east,north=transformer.transform(lon,lat)
            if not (0<east<800000 and 0<north<1300000):
                raise RuntimeError("PRTR EA WGS84 point source failed original BNG projection")
            bng[key]=(east,north)
        data=score_geo(bng,sa_points)
        byyear_identity={
            yr:{k[2] for k in original_pts if k[0]==yr and k[1]=="TOTAL - NITROGEN"}
            for yr in ("2017","2018")
        }
        result.update({
            "status":"PRTR_TRUE_YEAR_FACILITY_TO_CURRENT_WQE_SA_GEOGRAPHIC_SOURCE_CANDIDATES_COMPLETED",
            "UKPRTR_original_source_validation":prtrqa,
            "EA_current_SA_source_registry_n":len(sa_points),
            "EA_current_SA_source_registry_SHA256":saqa["source_all_station_code_SHA256"],
            "n_original_EA_urban_STW_totalN_facility_IDs_repeated_both_2017_2018":
                len(byyear_identity["2017"]&byyear_identity["2018"]),
            "n_original_EA_urban_STW_totalN_facility_IDs_either_2017_2018":
                len(byyear_identity["2017"]|byyear_identity["2018"]),
            "geographical_opportunity_by_true_source_year":data,
            "source_gates":{
                "strict_original_PRTR_EA_5f_NACE37_and_report_year_and_PN_verified":True,
                "current_WQE_5722_SA_registry_geometry_verified":True,
                "official_sampling_point_to_PRTR_facility_NationalID_crosswalk_verified":False,
                "physical_permitted_outfall_location_verified":False,
                "actual_2017_18_EA_WQE_TP_sampling_records_for_matched_source_PN":False,
                "true_upstream_WWTW_annual_mass_to_ChemPop_BIO_delivered":False,
                "same_year_PN_mass_at_same_strict_urban_WWTW_absent":True,
                "cannot_authorize_causal_BMWP_response_model":True
            }
        })
    except Exception as e:
        result.update({"status":"SOURCE_PRTR_FACILITY_TO_SA_GEO_GATE_STOP",
                       "source_error":type(e).__name__+": "+str(e)[:680]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_PRTR_2017_18_FACILITY_VS_SA_GEO="+json.dumps({
       "status":result["status"],"error":result.get("source_error"),
       "PRTR_source":result.get("UKPRTR_original_source_validation"),
       "N_facility_repeated_both_years":result.get("n_original_EA_urban_STW_totalN_facility_IDs_repeated_both_2017_2018"),
       "source_year_geographical_matches":result.get("geographical_opportunity_by_true_source_year"),
       "source_gates":result.get("source_gates")
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
