#!/usr/bin/env python3
"""Original 2018 EEA UK freshwater outfall -> EA final TP sampling -> BIO river.

SOURCE-ONLY original historical UK physical UWWTD discharge point joined
geometrically (not by a proven shared licence ID!) with current EA WQE
Water Company FINAL/TREATED SA sampling point: unique 250 m mutual match
and unique original UK UWWTD outlet per plant. Then actual 2018 SA
Phosphorus Total as P 0348 positive numeric mg/L exists. Use the ORIGINAL
historical EEA physical reported outfall coordinates, *not* SA coordinates,
to snap to official HydroRIVERS NEXT_DOWN and identify 2018 ChemPop
source-qualified BIO samples on downstream reaches <=50km, reversed
placebo, and same-reach unknown. Compare sample day pre90 and after90
with actual original 2018 BIO_SAMPLE_DATE.

This does NOT establish a matched permit licence ID, true WWTP pollutant
delivery, site-level P/N source mass kg/year, treatment intervention,
ecological selection, or BMI causality. No raw original IDs, locations,
TP concentration values, facility names or BIO responses rehosted.
"""
from __future__ import annotations
import argparse,collections,datetime as dt,hashlib,json,sys,tempfile
from pathlib import Path
import numpy as np
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261121_uk_prtr"))
import eea2018_physical_outlet_vs_EA_SA_source_spatial as earlier
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as registry
import ea_WQE_national_SA_actual2017_18_TP_chemistry_bio_screen as actual
import ea_SA_TP_to_ChemPop_location_sourcegate as bio_source
import ea_SA_TP_ChemPop_HydroRIVERS_directed_route_gate as flow
import ea_SA_TP_upstream_before_ChemPop_bio_dates as sample_date

ORIGINAL_UWWTD_2018_FW=1498
ORIGINAL_CURRENT_SA=5722
PRIMARY_MUTUAL_RADIUS_METERS=250
SNAP_RIVER_M=500
RIVER_DISTANCE_KM=50
SOURCE_YEAR=2018
BIO_YEAR_SITE_N=366
MIN_SOURCE_SITES=50
MIN_SOURCE_SITEYEARS=100
PRE_DAYS=90
POST_DAYS=90


def select_current_SA_2018_UWWTD_historic_outlets(uwwtd,outcount,sa):
    ids=sorted(uwwtd)
    registry_ids=sorted(sa)
    transf=Transformer.from_crs("EPSG:4326","EPSG:27700",always_xy=True)
    xy=[]
    for k in ids:
        lon,lat=uwwtd[k]
        e,n=transf.transform(lon,lat)
        if not(0<e<800000 and 0<n<1300000):
            raise RuntimeError("Original UWWTD 2018 FW outlet cannot convert to BNG")
        xy.append((e,n))
    a=np.asarray(xy,float)
    b=np.asarray([sa[k] for k in registry_ids],float)
    d=np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)
    if len(ids)!=ORIGINAL_UWWTD_2018_FW or len(registry_ids)!=ORIGINAL_CURRENT_SA:
        raise RuntimeError("2018 historical UWWTD outlet or modern SA source population changed")
    summaries={}
    matched_main={}
    for radius in (100,PRIMARY_MUTUAL_RADIUS_METERS,500):
        mask=d<=radius
        n_per_outlet=mask.sum(axis=1)
        n_per_sa=mask.sum(axis=0)
        matches={}
        for i,cnt in enumerate(n_per_outlet):
            if cnt!=1:continue
            native=ids[i]
            if outcount[native[1]]!=1:continue
            sa_index=int(np.flatnonzero(mask[i])[0])
            if n_per_sa[sa_index]!=1:continue
            code=registry_ids[sa_index]
            matches[code]={"outlet_id":native,
                "outlet_coord":tuple(xy[i]),
                "original2018_FW_outlet_to_SA_distance_m":float(d[i,sa_index])}
        if len(matches)!=len(set(m["outlet_id"] for m in matches.values())):
            raise RuntimeError("EEA UWWTD outlet-to-WQE SA mutual uniqueness broken")
        summaries[str(radius)+"m"]={
            "historical2018_UWWTD_FW_outlets_with_a_SA_candidate":int(mask.any(axis=1).sum()),
            "historical2018_UWWTD_FW_outlets_with_unique_SA_ignoring_reciprocal":int((n_per_outlet==1).sum()),
            "source_geographically_UNIQUE_mutual_SA_to_outlet_and_single_outlet_plant":len(matches),
            "official_current_EA_SA_codes_geographically_mutual":len(matches),
            "still_not_original_SA_to_UWWTD_permit_ID_verified":True}
        if radius==PRIMARY_MUTUAL_RADIUS_METERS:
            matched_main=matches
    qa={"original_UWWTD_2018_reported_freshwater_outlet_n":len(ids),
        "original_current_WQE_WaterCompany_final_treated_SA_n":len(registry_ids),
        "main_matching_tolerance_m":PRIMARY_MUTUAL_RADIUS_METERS,
        "source_direction_uses_original_HISTORIC_UWWTD_OUTLET_not_SA_sample_coord":True,
        "mutual_unique_original_source_geographic_candidates_by_radius":summaries,
        "main_WQE_SA_station_source_hash":hashlib.sha256(
            "\n".join(sorted(matched_main)).encode()).hexdigest(),
        "2026_SA_current_registration_not_historic2018_TP_monitoring":True}
    return matched_main,qa

def bio_routing(actual_2018,matching,bio_year_sites,bio,rivers,graph):
    source_geopoints={k:v["outlet_coord"] for k,v in matching.items()
        if k in actual_2018}
    b_points={k:bio[k] for k in bio_year_sites}
    source_snap,qa_stw=flow.snap(source_geopoints,rivers,
        "historical 2018 UK UWWTD actual freshwater outlet")
    bio_snap,qa_bio=flow.snap(b_points,rivers,
        "original 2018 ChemPop source-eligible BIO site")
    routes=collections.defaultdict(set)
    reverse=collections.defaultdict(set)
    same=collections.defaultdict(set)
    source_to_bio=collections.defaultdict(set)
    group=collections.Counter()
    upstream_distance=[]
    for stw,r in source_snap.items():
        dst,_=flow.downstream_trace(r["reach"],graph,maxkm=RIVER_DISTANCE_KM)
        for bid,br in bio_snap.items():
            if r["reach"]==br["reach"]:
                same[bid].add(stw)
                group["same_reach_unknown_pairs"]+=1
                continue
            if br["reach"] in dst:
                routes[bid].add(stw)
                source_to_bio[stw].add(bid)
                group["upstream_direction_pairs"]+=1
                upstream_distance.append(dst[br["reach"]])
                continue
            back,_=flow.downstream_trace(br["reach"],graph,maxkm=RIVER_DISTANCE_KM)
            if r["reach"] in back:
                reverse[bid].add(stw)
                group["reverse_STW_downstream_from_BIO_pairs"]+=1
    out={
        "n_true_2018_same_year_measured_SA_with_unique_HISTORIC_EEA_discharge_point_geo_match":
            len(source_geopoints),
        "n_2018_reported_outlet_matched_actual_SA_TP_sampling_site_and_snapped_to_river500m":
            len(source_snap),
        "n_2018_original_ChemPop_bio_sites":len(b_points),
        "n_2018_original_ChemPop_bio_sites_snap_nearest_river500m":len(bio_snap),
        "unique_BIO_sites_with_directional_upstream_historical_UWWTD_point_and_actual_EA_TP":
            len(routes),
        "unique_BIO_sites_with_STW_source_only_DOWNSTREAM_placebo":len(reverse),
        "unique_BIO_sites_same_reach_source_order_unknown":len(same),
        "unique_UWWTD_historical_physical_discharge_points_on_routed_BIO_upstream_reaches":
            len(source_to_bio),
        "true_upstream_path_distance_km_median":float(np.median(upstream_distance))
            if upstream_distance else None,
        "pair_count_source_qa":dict(group),
        "original_UWWTD_reported_outlet_coordinates_not_verification_of_sewage_chemical_arrival":True,
        "has_no_measured_source_2018_UWWTD_water_release_PN_kg":True
    }
    return routes,reverse,same,out,qa_stw,qa_bio

def samples_before_after(routes,measured):
    original_dates,date_source=sample_date.origin_bio_sample_dates()
    if len([a for a,b in original_dates if a==SOURCE_YEAR])!=BIO_YEAR_SITE_N:
        raise RuntimeError("Original 2018 ChemPop native sample site-date source support drift")
    out={
        "no_same_day_assumed_previous":True,
        "time_window_pre90_max_days":PRE_DAYS,
        "time_window_post90_max_days":POST_DAYS,
        "2018_source_actual_BIO_SAMPLE_DATE_provenance":date_source,
        "source_day_counts":{"BIO_sites_with_at_least_one_pre90_TP":0,
          "BIO_sites_with_all_actual_days_pre90_TP":0,
          "BIO_sites_with_at_least_one_post90_TP":0,
          "BIO_sites_same_day_TP_unknown_intraday_order":0,
          "BIO_site_dates_total":0,"BIO_site_dates_pre90":0,
          "BIO_site_dates_post90":0,
          "BIO_site_dates_same_day":0},
    }
    c=out["source_day_counts"]
    for bid,stw in routes.items():
        bd=sorted(original_dates[(SOURCE_YEAR,bid)])
        if not bd:raise RuntimeError("True 2018 BIO source dates disappeared")
        obs=set()
        for s in stw:
            for d in measured[s]["valid_days"]:
                day=dt.date.fromisoformat(d)
                if day.year==SOURCE_YEAR:obs.add(day)
        if not obs:raise RuntimeError("Actually monitored source STW TP has no 2018 date")
        ahead=[any(0<(b-t).days<=PRE_DAYS for t in obs) for b in bd]
        after=[any(0<(t-b).days<=POST_DAYS for t in obs) for b in bd]
        same=[b in obs for b in bd]
        c["BIO_sites_with_at_least_one_pre90_TP"]+=int(any(ahead))
        c["BIO_sites_with_all_actual_days_pre90_TP"]+=int(all(ahead))
        c["BIO_sites_with_at_least_one_post90_TP"]+=int(any(after))
        c["BIO_sites_same_day_TP_unknown_intraday_order"]+=int(any(same))
        c["BIO_site_dates_total"]+=len(bd)
        c["BIO_site_dates_pre90"]+=sum(ahead)
        c["BIO_site_dates_post90"]+=sum(after)
        c["BIO_site_dates_same_day"]+=sum(same)
    return out

def main(args):
    out={"schema":"GESUI_20261121_HISTORIC_UWWTD2018_REPORTED_OUTLET_TO_TRUE_EA_TP_CHEMPOP_UPSTREAM_SOURCE_ONLY_V1",
         "source_only_no_BMWP_or_causal_PN_flux":True,
         "UK_UWWTD_2018_plant_measured_both_NP_positive_zero":True,
         "source_outlet_to_EA_WQE_SA_native_permit_ID_not_verified":True,
         "historical_outlet_instead_of_facility_or_SA_geocentre":True,
         "status":"SOURCE_STOP"}
    try:
        outlets,plantcounts,eeaqa=earlier.exact_source()
        registered,saqa=registry.collect()
        matched,matchqa=select_current_SA_2018_UWWTD_historic_outlets(
             outlets,plantcounts,registered)
        if not matched:
            raise RuntimeError("Original official historic outfall-to-current SA mutual geosource contains zero")
        actual_samp,actualqa=actual.candidate_history(sorted(matched))
        sources=actual_samp[SOURCE_YEAR]
        bio,bioyear,bioqa=bio_source.chempop_geo()
        if len(bioyear[SOURCE_YEAR])!=BIO_YEAR_SITE_N:
            raise RuntimeError("Original source-qualified 2018 BIO population no longer 366")
        with tempfile.TemporaryDirectory(prefix="gesui_eea2018_originalFW_routes_") as folder:
            rivers,graph,riversqa=flow.get_rivers(Path(folder))
            positive,reverse,same,geo,STWQA,BIOQA=bio_routing(
               sources,matched,bioyear[SOURCE_YEAR],bio,rivers,graph)
        temporal=samples_before_after(positive,sources)
        out.update({
            "status":"SOURCE_2018_HISTORICAL_UWWTD_FW_OUTLET_TO_TRUE_TP_SA_AND_BIO_DIRECTED_GATES_AUDITED",
            "original_EEA2018_UKFW_source_QA":eeaqa,
            "historical_outlet_to_current_EA_WQE_SA_source_geography":matchqa,
            "official_EA_SA_2018_actual_TP_source_qa":{
                k:v for k,v in actualqa.items() if k not in ("source_page_SHA256",)},
            "source_original_HydroRIVERS_v1":riversqa,
            "original_historic_UWWTD_outlet_river_snap":STWQA,
            "original_2018_BIO_river_snap":BIOQA,
            "original_2018_TRUE_UWWTD_outlet_EA_TP_BIO_reach_upstream_vs_reverse":geo,
            "genuine_original_2018_SA_TP_sample_before_BIO_sample_date":temporal,
            "source_gates":{
                "main_250m_mutual_geographic_unique_source_match_found":True,
                "actual_2018_final_EA_SA_totalP_monitoring_confirmed":len(sources)>0,
                "2018_bio_historical_dischargepoint_NETWORK_upstream_site_count_ge50":
                    geo["unique_BIO_sites_with_directional_upstream_historical_UWWTD_point_and_actual_EA_TP"]>=MIN_SOURCE_SITES,
                "2018_bio_historical_dischargepoint_NETWORK_upstream_site_years_ge100":
                    geo["unique_BIO_sites_with_directional_upstream_historical_UWWTD_point_and_actual_EA_TP"]>=MIN_SOURCE_SITEYEARS,
                "permit_original_same_facility_identity_confirmed":False,
                "annual_NP_mass_kg_discharge_confirmed":False,
                "source_to_BMWP_ecological_causal_regression_authorized":False,
            }
        })
    except Exception as exc:
        out.update({"status":"SOURCE_HISTORICAL_OUTLET_OR_NETWORK_GATE_STOP",
                    "source_error":type(exc).__name__+": "+str(exc)[:680]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_EEA2018_UK_HISTORIC_FW_OUTLET_EA_REAL_TP_TO_BIO="+json.dumps({
       "status":out["status"],"error":out.get("source_error"),
       "source_matching":out.get("historical_outlet_to_current_EA_WQE_SA_source_geography"),
       "source_chemistry":out.get("official_EA_SA_2018_actual_TP_source_qa"),
       "river_directional":out.get("original_2018_TRUE_UWWTD_outlet_EA_TP_BIO_reach_upstream_vs_reverse"),
       "source_sample_dates":out.get("genuine_original_2018_SA_TP_sample_before_BIO_sample_date",{}).get("source_day_counts"),
       "gates":out.get("source_gates"),
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
