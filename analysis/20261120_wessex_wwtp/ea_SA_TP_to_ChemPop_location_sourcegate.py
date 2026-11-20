#!/usr/bin/env python3
"""EA Water Company treated effluent TP source locations vs ChemPop biological sites.

Strict SOURCE ONLY: 2017-18 official EA WQE 0348 SA treated effluent point
BNG easting/northing, official EIDC ChemPop v2 original BIO site matched BNG
geography and 2017/18 source-eligible biological sample site-years.

Distances are straight-line geographic OPPORTUNITY SCREENS, NOT directed
upstream STW exposure; no downstream reach matching, effluent P/N mass or
biotic response coefficients are calculated. No individual source points,
STW names/coordinates or biological samples are republished.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,sys,tempfile,zipfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261119_europe_benthic_taxa"))
import chempop_v2_source_census as csource
import chempop_v2_bmwp_filter as cfilter
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_wqe_2017_2018_SA_tp_siteyear_census as ea

TARGET_YEARS=(2017,2018)
N_EXPECTED_SA={2017:46,2018:49}
N_EXPECTED_SHARED=46
EXPECTED_SOURCE_CHEMPOP_SAMPLES=7836
EXPECTED_SOURCE_CHEMPOP_SITES=707
EXPECTED_BMWP_SCORE_SHA="4f9be755ce6d24e645b46705dbec9c0367bc71a9ce6b335d3bf2657d8c10ac5a"
SCREENS_METERS=(500,1000,2000,5000,10000)
MAX_COORD_DRIFT_METERS=100.0
MIN_CANDIDATE_SITES_FOR_NETWORK_EXPLORATION=10


def number(v):
    try:
        x=float(str(v).strip())
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None

def bng_source_valid(e,n):
    return (e is not None and n is not None and
            0 < e < 800000 and 0 < n < 1300000)

def ea_geo():
    year={}
    all_station=collections.defaultdict(list)
    all_dates=collections.defaultdict(set)
    count=collections.Counter()
    site_uwwtd=collections.defaultdict(set)
    yearly_receipts={}
    for yr in TARGET_YEARS:
        page_rows=[]
        receipts=[]
        completed=False
        for j in range(ea.MAX_PAGES_PER_YEAR):
            rows,receipt=ea.fetch(yr,j*ea.LIMIT)
            receipts.append(receipt)
            if rows is None:raise RuntimeError("Official EA WQE page source STOP")
            page_rows.extend(rows)
            if len(rows)<ea.LIMIT:
                completed=True;break
        if not completed:raise RuntimeError("Official EA WQE source page cap hit")
        scan=collections.defaultdict(list)
        for q in page_rows:
            if q["determinand.notation"].strip()!=ea.ANALYTE_TOTAL_P_AS_P:
                raise RuntimeError("Official source returned wrong constituent")
            if q["samplingPoint.samplingPointType"].strip()!=ea.SA:
                count["source_non_SA_other_type_rows"]+=1
                continue
            if q["unit"].strip()!="MILLIGRAM PER LITRE":
                count["source_SA_non_mgL_rows"]+=1
                continue
            if q["sampleMaterialType"].strip()!="FINAL SEWAGE EFFLUENT":
                count["source_SA_unverified_material_rows"]+=1
                continue
            key=q["samplingPoint.notation"].strip()
            e=number(q.get("samplingPoint.easting"))
            n=number(q.get("samplingPoint.northing"))
            if not key or not bng_source_valid(e,n):
                count["SA_missing_source_BNG_identity"]+=1
                continue
            val=number(q["result"])
            if val is None or val<0:
                count["SA_invalid_or_censored_TP"]+=1
                continue
            rawdate=q["phenomenonTime"].strip()[:10]
            if not rawdate.startswith(str(yr)):
                raise RuntimeError("Original source date-year mismatch")
            scan[key].append((e,n))
            all_station[key].append((e,n))
            all_dates[(yr,key)].add(rawdate)
            if "UWWTD" in q["samplingPurpose"].upper():
                site_uwwtd[yr].add(key)
            count["observed_SA_TP_mgL_"+str(yr)]+=1
        yearly_receipts[str(yr)]={
            "source_pages":len(receipts),
            "all_source_observation_rows":len(page_rows),
            "source_page_sha256":[r["sha256"] for r in receipts],
            "source_SA_eligible_unique_station_n":len(scan),
            "source_SA_eligible_unique_station_dates":sum(len(all_dates[(yr,k)]) for k in scan),
            "source_SA_UWWTD_station_n":len(site_uwwtd[yr])}
        year[yr]=scan
        if len(scan)!=N_EXPECTED_SA[yr]:
            raise RuntimeError("Original source initial official native SA station count drift from audited "+str(yr))
    common=set(year[2017])&set(year[2018])
    if len(common)!=N_EXPECTED_SHARED:
        raise RuntimeError("Original EA SA repeated stations cross-year source count changed")
    coords={}
    drift=[]
    for station,pts in all_station.items():
        px=np.asarray(pts,float)
        center=px.mean(axis=0)
        if not bng_source_valid(*center):
            raise RuntimeError("Original EA SA station BNG centroid invalid")
        deviations=np.linalg.norm(px-center[None,:],axis=1)
        drift.append(float(deviations.max()))
        if deviations.max()>MAX_COORD_DRIFT_METERS:
            raise RuntimeError("Original EA native sampling point moved >100m")
        coords[station]=tuple(center)
    for yr in TARGET_YEARS:
        year[yr]={k:coords[k] for k in year[yr]}
    qa={
        "official_WQE_source_point_type":ea.SA,
        "source_determinand":ea.ANALYTE_TOTAL_P_AS_P,
        "source_unit":"mg/L",
        "actual_source_sample_material":"FINAL SEWAGE EFFLUENT",
        "source_page_coverage":yearly_receipts,
        "source_station_same_year_2017":len(year[2017]),
        "source_station_same_year_2018":len(year[2018]),
        "source_station_both_years":len(common),
        "station_geometry_max_within_station_drift_m":max(drift) if drift else None,
        "station_geometry_max_drift_gate_m":MAX_COORD_DRIFT_METERS,
        "candidate_STW_total_unique":len(coords),
        "SA_STW_native_point_codes_SHA256":hashlib.sha256(
           ("\n".join(sorted(coords))).encode()).hexdigest(),
        "no_outfall_to_river_topology_or_source_mass_inferred":True,
    }
    return year,site_uwwtd,qa


def chempop_geo():
    with tempfile.TemporaryDirectory(prefix="gesui_EIDC_sitegeo_") as tmp:
        d=Path(tmp)
        receipt=csource.retrieve(csource.EIDC_URL,d/"original.zip",csource.SOURCE_DOWNLOAD_LIMIT)
        scores,score_receipt=csource.score_reference(d)
        if score_receipt["SHA256"]!=EXPECTED_BMWP_SCORE_SHA:
            raise RuntimeError("Original published BMWP family scoring dictionary SHA changed")
        with zipfile.ZipFile(d/"original.zip") as z:
            names=z.namelist()
            if len(names)!=326:raise RuntimeError("Original EIDC v2 archive member count changed")
            sites={}
            site_missing=collections.Counter()
            for region in csource.PREFIXES:
                path=csource.find_one_by_region(names,region,r"sitematchedlocation\.csv")
                for row in csource.csv_rows(z,path,["BIO_SITE_ID","BIO_EASTING","BIO_NORTHING"]):
                    site=row["BIO_SITE_ID"].strip()
                    if not site:continue
                    e,n=number(row["BIO_EASTING"]),number(row["BIO_NORTHING"])
                    if not bng_source_valid(e,n):
                        site_missing[region]+=1;continue
                    key=region+"||"+site
                    point=(e,n)
                    if key in sites and float(np.hypot(sites[key][0]-e,sites[key][1]-n))>1e-9:
                        raise RuntimeError("Original ChemPop source site ID maps to conflicting BNG coordinates")
                    sites[key]=point
            raw,weights,orders,prior=cfilter.construct_eligible(z,names,scores)
        if (len(raw)!=EXPECTED_SOURCE_CHEMPOP_SAMPLES or
            prior["source_gate"]["eligible_site_n"]!=EXPECTED_SOURCE_CHEMPOP_SITES):
            raise RuntimeError("Original 2010-18 ChemPop source eligibility drift")
        survey=collections.defaultdict(set)
        for row in raw:
            if row["year"] in TARGET_YEARS:
                key=row["site_region_id"]
                if key not in sites:site_missing[row["region"]+"_eligible_missing"]+=1
                else:survey[row["year"]].add(key)
    qa={
      "original_326_member_EIDC_ZIP_SHA256":receipt["SHA256"],
      "original_BMWP_source_SHA256":score_receipt["SHA256"],
      "original_2010_18_source_eligible_sample_n":len(raw),
      "original_2010_18_source_eligible_site_n":prior["source_gate"]["eligible_site_n"],
      "n_BNG_source_sites_with_coordinates":len(sites),
      "source_BNG_missing_by_region":dict(site_missing),
      "eligible_source_bio_sites_2017":len(survey[2017]),
      "eligible_source_bio_sites_2018":len(survey[2018]),
      "eligible_source_bio_sites_either_2017_or_2018":len(survey[2017]|survey[2018]),
      "native_BIO_ID_source_geo_join_verified":True,
    }
    return sites,survey,qa


def proximity(sources,sites):
    src_keys=sorted(sources)
    site_keys=sorted(sites)
    if not src_keys or not site_keys:
        return {"status":"SOURCE_NO_POINTS","n_original_source_STW":len(src_keys),
                "n_original_bio_sites":len(site_keys)}
    a=np.asarray([sites[k] for k in site_keys],float)
    b=np.asarray([sources[k] for k in src_keys],float)
    distances=np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)
    nearest=distances.min(axis=1)
    out={
        "n_original_source_SA_stations":len(src_keys),
        "n_original_ChemPop_sites":len(site_keys),
        "nearest_STW_straightline_distance_m_p10_p50_p90":[float(x) for x in np.quantile(nearest,[.1,.5,.9])],
        "thresholds":{},
        "straightline_is_not_upstream_flow_or_nutrient_dose":True,
    }
    for metres in SCREENS_METERS:
        edges=distances<=metres
        out["thresholds"][str(metres)+"m"]={
            "bio_site_count_with_any_geographically_near_SA_STW":int(edges.any(axis=1).sum()),
            "SA_STW_count_with_any_geographically_near_BIO_site":int(edges.any(axis=0).sum()),
            "geometric_site_STW_candidate_pairs":int(edges.sum()),
            "does_not_establish_downstream_connectivity":True,
        }
    return out


def main(a):
    out={
      "schema":"GESUI_EA_SA_TP_EIDC_CHEMPOP_BNG_SOURCE_LOCATION_SCREEN_20261120_V1",
      "source_only_no_environmental_or_biological_regressions":True,
      "not_WWTP_measured_total_N_or_P_kg_per_year":True,
      "not_directional_network_route_or_causal_STW_exposure":True,
      "not_source_eligible_BIO_site_equal_to_downstream":True,
      "status":"SOURCE_STOP",
    }
    try:
        ea_year,uwwtd,qa_ea=ea_geo()
        bio,survey,qa_bio=chempop_geo()
        out["EA_SA_tp_source"]=qa_ea
        out["ChemPop_original_biological_site_source"]=qa_bio
        all_source={k:v for year in ea_year.values() for k,v in year.items()}
        out["all_1519_source_biological_sites_screen"]=proximity(all_source,bio)
        by_year={}
        for year in TARGET_YEARS:
            candidate={k:bio[k] for k in survey[year] if k in bio}
            by_year[str(year)]={
                "same_source_observation_year":year,
                "screen_on_all_SA_TP_final_WaterCompany":proximity(ea_year[year],candidate),
                "screen_on_SA_TP_UWWTD_purpose_only":proximity(
                    {k:ea_year[year][k] for k in uwwtd[year] if k in ea_year[year]},candidate),
            }
        out["matched_biology_years"]=by_year
        n_bio_candidates=max(
            by_year[str(y)]["screen_on_all_SA_TP_final_WaterCompany"].get("thresholds",{})
            .get("10000m",{}).get("bio_site_count_with_any_geographically_near_SA_STW",0)
            for y in TARGET_YEARS)
        out["source_gates"]={
            "original_2017_18_watercompany_final_effluent_SA_TP_verified":True,
            "both_years_native_SA_sites_and_BNG_coordinates_stable":True,
            "original_ChemPop_BNG_site_geometry_verified":True,
            "same_year_bio_SA_geographic_screen_has_10_sites_within_10km":
              n_bio_candidates>=MIN_CANDIDATE_SITES_FOR_NETWORK_EXPLORATION,
            "source_to_receiver_directed_hydrological_routing_verified":False,
            "site_source_annual_mass_and_kg_per_year_verified":False,
            "eligible_geographic_proximity_alone_allows_ecological_causal_models":False,
        }
        out["status"]="SOURCE_BNG_PROXIMITY_SCREEN_COMPLETE_NOT_ROUTED"
    except Exception as err:
        out.update({"status":"SOURCE_GEOMETRY_GATE_STOP",
                    "error":type(err).__name__+": "+str(err)[:650]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_EA_SA_TP_CHEMPOP_BNG_SCREEN="+json.dumps(out,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
