#!/usr/bin/env python3
"""Original GBIF source-OPPORTUNITY only at NEW outlet-anchored ±1.5km river points.

128 original official GBIF occurrence-search COUNT queries:
16 real outfalls ×2 up/down sites ×2 temporal windows ×2 broad kingdoms.
Source-based support gate before any numeric EIV, plant community, butterfly
or downstream treatment effect. Coordinates from source-authentic directed
HydroRIVERS geometry, not the prior 48 historical mid-reach GBIF centres.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,time
from pathlib import Path
from urllib.parse import urlencode
import requests

URL="https://api.gbif.org/v1/occurrence/search"
SOURCE_GEOMETRY_SCHEMA="GESUI_20261116_REAL_LOCAL_1500M_HYDROWASTE_HYDRORIVERS_SOURCE_ONLY_V1"
SCHEMA="GESUI_20261116_WWTP_16_LOCAL_UPDOWN_250M_2PERIOD_GBIF_SOURCE_COUNTS_V1"
YEARS=[("2017-2019","2017,2019"),("2020-2022","2020,2022")]
SPECIES_GROUPS=[("plant",6),("nonfocal_Animalia",1)]
ROLES=("UPSTREAM","DOWNSTREAM")
SITE_MIN_PLANT_RECORDS=20
MIN_FACILITIES=10
MIN_COUNTRIES=3
DISTANCE_KM=.25
AGENT="gesui-local250m-outfall-source-counts-science/20261116"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_sites(path):
    gate=json.loads((path/"ORIGINAL16_REAL_WWTP_LOCAL_1500M_RIVER_SITE_GEOMETRY.json").read_text())
    reg=path/"ORIGINAL16_LOCAL_UPSTREAM_DOWNSTREAM_1500M_WWTP_SOURCE_SITES.csv"
    if (gate.get("schema")!=SOURCE_GEOMETRY_SCHEMA or
        gate.get("status")!="PASS_OUTFALL_1500M_GEOGRAPHY_ONLY" or
        gate.get("source_true_chain_outfalls_geometry_pass")!=16 or
        not reg.exists()):
        raise RuntimeError("Previous 16 original outlet-anchored actual river GIS source gate not verified")
    with reg.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    if len(rows)!=16 or len({q["original_outfall_id"] for q in rows})!=16:
        raise RuntimeError("Source actual new local outlet pairs missing")
    if len({q["country"] for q in rows})!=4:
        raise RuntimeError("Source independent-country 16 outlet pilot changed")
    for r in rows:
        if r["source_original_river_true_link_and_orientation_checked"].lower()!="true":
            raise RuntimeError("Actual directed source chain not verified")
        if float(r["original_up_down_river_path_distance_m"])!=3000.:
            raise RuntimeError("Actual upper/lower 1.5km source geometry changed")
        if float(r["candidate_future_GBIF_buffer_radius_m"])!=250.:
            raise RuntimeError("Original frozen 250m buffer changed")
    return rows,{"previous_source_run_original_geometry_schema":gate["schema"],
                 "signed_source_geographic_site_registry_SHA256":sha(reg)}

def query(session,lat,lon,year,taxon_key):
    pars={"taxonKey":taxon_key,"year":year,
          "geoDistance":f"{lat:.7f},{lon:.7f},{DISTANCE_KM:.2f}km",
          "hasCoordinate":"true","hasGeospatialIssue":"false",
          "occurrenceStatus":"PRESENT","limit":0}
    prepared=requests.Request("GET",URL,params=pars).prepare()
    last=None
    for i in range(4):
        try:
            r=session.send(prepared,timeout=(12,75))
            if r.status_code in (429,500,502,503,504):
                raise RuntimeError("Official GBIF temporary HTTP "+str(r.status_code))
            r.raise_for_status()
            obj=r.json()
            if not isinstance(obj,dict) or "count" not in obj:
                raise RuntimeError("Official GBIF query response lacks original count")
            count=obj["count"]
            if not isinstance(count,int) or count<0:
                raise RuntimeError("GBIF count noninteger, invalid source response")
            return {
                "count":count,
                "original_JSON_sha256":hashlib.sha256(r.content).hexdigest(),
                "original_source_url_SHA256":hashlib.sha256(prepared.url.encode()).hexdigest(),
                "GBIF_source_original_query_url":prepared.url}
        except (requests.RequestException,RuntimeError,ValueError) as e:
            last=e
            if i<3:time.sleep(2+2*i)
    raise RuntimeError("Official original GBIF source unavailable for exact study site/year: "+str(last))

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    sites,qa=source_sites(a.geometry)
    session=requests.Session()
    session.headers.update({"Accept":"application/json","User-Agent":AGENT})
    results=[]
    by={}
    for source in sorted(sites,key=lambda q:q["original_outfall_id"]):
        sid=source["original_outfall_id"]
        for role in ROLES:
            lat=float(source[f"new_{role}_1500m_lat"])
            lon=float(source[f"new_{role}_1500m_lon"])
            for p,y in YEARS:
                for cat,key in SPECIES_GROUPS:
                    val=query(session,lat,lon,y,key)
                    item={"original_facility":sid,"country":source["country"],
                        "site_role":role,"GBIF_3year_window":p,
                        "original_GBIF_kingdom_group":cat,"original_GBIF_taxonKey":key,
                        "source_up_down_directed_river_path_total_m":3000,
                        "source_original_1p5km_each_side_of_WWTP":True,
                        "actual_original_site_radius_m":250,
                        "original_GBIF_occurrence_count":val["count"],
                        "GBIF_original_count_JSON_sha256":val["original_JSON_sha256"],
                        "original_API_exact_parameters_hash":val["original_source_url_SHA256"],
                        "original_API_source_query":val["GBIF_source_original_query_url"],
                        "source_0_count_is_not_true_absence":True,
                        "no_trait_or_sewage_response_fitted":True}
                    results.append(item)
                    by[(sid,source["country"],role,p,cat)]=val["count"]
                    time.sleep(.1)
        print("GESUI_LOCAL_GBIF_SOURCE_COUNTS_OUTFALL_COMPLETE="+sid,flush=True)
    if len(results)!=16*2*2*2 or len(by)!=len(results):
        raise RuntimeError("Missing original source count queries 16×2×2×2")
    eligibility={}
    for period,_ in YEARS:
        pairs=[]
        for s in sites:
            sid=s["original_outfall_id"];c=s["country"]
            up=by[(sid,c,"UPSTREAM",period,"plant")]
            down=by[(sid,c,"DOWNSTREAM",period,"plant")]
            if up>=SITE_MIN_PLANT_RECORDS and down>=SITE_MIN_PLANT_RECORDS:
                pairs.append((sid,c))
        eligibility[period]={"n_complete_source_facilities_with_at_least20_PLANT_records_on_both_sides":len(pairs),
            "source_countries":sorted(set(c for _,c in pairs)),
            "source_country_group_count":len({c for _,c in pairs}),
            "source_qualified_original_facility_ids":[x for x,_ in pairs],
            "preoutcome_gate_pass":len(pairs)>=MIN_FACILITIES and len({c for _,c in pairs})>=MIN_COUNTRIES}
    both=set(eligibility[YEARS[0][0]]["source_qualified_original_facility_ids"])&set(
             eligibility[YEARS[1][0]]["source_qualified_original_facility_ids"])
    countries={s["country"] for s in sites if s["original_outfall_id"] in both}
    temporal_pass=len(both)>=MIN_FACILITIES and len(countries)>=MIN_COUNTRIES
    source={"schema":SCHEMA,
       "original_signed_16_directed_outfall_geographic_source":qa,
       "original_observation_query_count":len(results),
       "source_success_all_genuine_GBIF_queries":True,
       "geodistance_radius_km":DISTANCE_KM,
       "predeclared_min_Plantae_occurrences_per_side":SITE_MIN_PLANT_RECORDS,
       "predeclared_min_distinct_outfalls":MIN_FACILITIES,
       "predeclared_min_distinct_country_groups":MIN_COUNTRIES,
       "source_per_period":eligibility,
       "source_two_period_same_facility_qualified_count":len(both),
       "source_two_period_country_count":len(countries),
       "source_two_period_gate_pass":temporal_pass,
       "source_main_2020_2022_gate_pass":eligibility["2020-2022"]["preoutcome_gate_pass"],
       "source_status":"PASS_LOCAL_GBIF_RECORD_OPPORTUNITY_ONLY" if eligibility["2020-2022"]["preoutcome_gate_pass"] else
                       "STOP_LOCAL_GBIF_RECORDING_OPPORTUNITY_250M",
       "no_gbif_species_counts_or_EIV_trait_outcomes_opened":True,
       "no_actual_wwtp_discharge_annual_measured_NP":True,
       "not_individual_plant_or_butterfly_abundance":True,
       "not_same_year_WWTP_water_nutrient_delivery":True,
       "original_source_query_timestamp_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat()}
    (a.output/"LOCAL16_250M_TWO_PERIOD_GBIF_PLANT_AND_NONFOCAL_SOURCE_GATES.json").write_text(
       json.dumps(source,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (a.output/"LOCAL16_250M_TWO_PERIOD_GBIF_OCCURRENCE_SOURCE_COUNTS.csv").open(
        "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]))
        w.writeheader();w.writerows(results)
    print("GESUI_LOCAL16_250M_GBIF_RECORDING_GATE="+json.dumps({
        "status":source["source_status"],
        "primary_2020_2022":eligibility["2020-2022"],
        "early_2017_2019":eligibility["2017-2019"],
        "two_period_same_outfall_n":len(both),
        "source_queries":len(results),
        "no_plant_or_butterfly_effect":True},sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geometry",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
