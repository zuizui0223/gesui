#!/usr/bin/env python3
"""After strictly original 250m STOP, check fixed larger 500m and 1km radius
GBIF Plantae source-opportunity at true ±1.5km physical outlet positions.

Still source-only; records, not biological population. Ruling out overlapping
spatial buffers at 1km prevents source pseudoreplication by concentric circles.
128 new official GBIF count queries and fixed 10 outlet/3 country gate each.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,json,statistics,time
from pathlib import Path
from urllib.parse import urlencode
import requests

SOURCE_GEOMETRY_SCHEMA="GESUI_20261116_REAL_LOCAL_1500M_HYDROWASTE_HYDRORIVERS_SOURCE_ONLY_V1"
SCHEMA="GESUI_20261116_REAL_OUTFALL_500M_1000M_GBIF_2PERIOD_SOURCE_SENSITIVITY_V1"
URL="https://api.gbif.org/v1/occurrence/search"
RADIUS_KM=(.5,1.0)
PERIODS=(("2017-2019","2017,2019"),("2020-2022","2020,2022"))
PLANT_KEY=6
RECORD_THRESHOLD=20
MIN_FAC=10;MIN_COUNTRY=3

def sha(x):return hashlib.sha256(x).hexdigest()
def load_source(a):
    j=json.loads((a/"ORIGINAL16_REAL_WWTP_LOCAL_1500M_RIVER_SITE_GEOMETRY.json").read_text())
    if j.get("schema")!=SOURCE_GEOMETRY_SCHEMA or j.get("source_true_chain_outfalls_geometry_pass")!=16:
        raise RuntimeError("True signed source outlet/river topology original16 is not validated")
    p=a/"ORIGINAL16_LOCAL_UPSTREAM_DOWNSTREAM_1500M_WWTP_SOURCE_SITES.csv"
    with p.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    if len(rows)!=16 or len({q["original_outfall_id"] for q in rows})!=16:
        raise RuntimeError("Original source16 local coordinates missing")
    if any(float(q["new_up_down_Euclidean_distance_m"])<=2000 for q in rows):
        raise RuntimeError("1km up/down radial buffers would overlap, source sensitivity must STOP")
    return rows,sha(p.read_bytes())

def gbif(session,lat,lon,year,radius):
    params={"taxonKey":PLANT_KEY,"year":year,
       "geoDistance":f"{lat:.7f},{lon:.7f},{radius:.2f}km",
       "hasCoordinate":"true","hasGeospatialIssue":"false",
       "occurrenceStatus":"PRESENT","limit":0}
    req=requests.Request("GET",URL,params=params).prepare()
    for attempt in range(4):
        try:
            res=session.send(req,timeout=(12,65))
            if res.status_code in (429,500,502,503,504):
                raise requests.RequestException("Temporary GBIF status "+str(res.status_code))
            res.raise_for_status()
            j=res.json()
            if not isinstance(j,dict) or not isinstance(j.get("count"),int) or j["count"]<0:
                raise ValueError("Official original GBIF count response missing/invalid")
            return int(j["count"]),sha(res.content),sha(req.url.encode())
        except (requests.RequestException,ValueError):
            if attempt==3:raise
            time.sleep(2+2*attempt)

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    source,geo_sha=load_source(a.geometry)
    client=requests.Session()
    client.headers.update({"User-Agent":"gesui-sourceonly-local500m1km-20261116","Accept":"application/json"})
    rows=[]
    counts={}
    for study in sorted(source,key=lambda x:x["original_outfall_id"]):
        f=study["original_outfall_id"]
        for radius in RADIUS_KM:
            for name,year in PERIODS:
                for role in ("UPSTREAM","DOWNSTREAM"):
                    lat=float(study[f"new_{role}_1500m_lat"])
                    lon=float(study[f"new_{role}_1500m_lon"])
                    count,digest,queryhash=gbif(client,lat,lon,year,radius)
                    row={"outfall_id":f,"country":study["country"],"role":role,
                         "GBIF_period":name,"source_radius_km":radius,
                         "original_taxonKey":PLANT_KEY,
                         "source_original_occurrence_count":count,
                         "source_GBIF_original_count_JSON_sha":digest,
                         "source_official_query_parameters_SHA":queryhash,
                         "original_source_up_down_centres_river_distance_m":3000.,
                         "original_source_up_down_straight_distance_m":float(study["new_up_down_Euclidean_distance_m"]),
                         "source_buffers_nonoverlapping":True,
                         "source_no_species_score_or_real_vegetation_effect":True}
                    rows.append(row);counts[(f,study["country"],name,radius,role)]=count
                    time.sleep(.12)
        print("GBIF_LARGER_RADIUS_SOURCE_DONE="+f,flush=True)
    if len(rows)!=128 or len(counts)!=128:raise RuntimeError("Incomplete authentic 16x2x2x2 source count study")
    summary={}
    for radius in RADIUS_KM:
        for period,_ in PERIODS:
            eligible=[(q["original_outfall_id"],q["country"]) for q in source
               if min(counts[(q["original_outfall_id"],q["country"],period,radius,role)]
                      for role in ("UPSTREAM","DOWNSTREAM"))>=RECORD_THRESHOLD]
            allup=[counts[(q["original_outfall_id"],q["country"],period,radius,"UPSTREAM")] for q in source]
            alldown=[counts[(q["original_outfall_id"],q["country"],period,radius,"DOWNSTREAM")] for q in source]
            x={"source_eligible_outfall_count":len(eligible),
               "source_eligible_countries":sorted(set(c for _,c in eligible)),
               "source_eligible_country_count":len({c for _,c in eligible}),
               "source_eligible_facility_IDs":[i for i,_ in eligible],
               "GBIF_UPSTREAM_original_count_median":statistics.median(allup),
               "GBIF_DOWNSTREAM_original_count_median":statistics.median(alldown),
               "GBIF_UPSTREAM_zero_count_sites":sum(z==0 for z in allup),
               "GBIF_DOWNSTREAM_zero_count_sites":sum(z==0 for z in alldown),
               "GBIF_UPSTREAM_highest_count":max(allup),
               "GBIF_DOWNSTREAM_highest_count":max(alldown),
               "source_gate_pass":len(eligible)>=MIN_FAC and len({c for _,c in eligible})>=MIN_COUNTRY}
            summary[f"radius_{radius:.2f}km_period_{period}"]=x
    outcome={"schema":SCHEMA,
         "status":"SOURCE_ONLY_NUMERIC_GBIF_RECORD_OPPORTUNITY_VALIDATED",
         "original_true_16_outfall_geography_SHA256":geo_sha,
         "original_genuine_public_GBIF_source_queries":len(rows),
         "same_source_outfall_ids_and_geographic_area_16":True,
         "fixed_test_radii_km":RADIUS_KM,
         "fixed_year_windows":PERIODS,
         "source_PLANT_original_record_threshold_per_side":RECORD_THRESHOLD,
         "min_distinct_outfall_source_units":MIN_FAC,
         "min_countries_with_supported_outfalls":MIN_COUNTRY,
         "source_250m_original_STOP_is_not_replaced":True,
         "source_same_year_wwtp_NP_actual_dose_unverified":True,
         "no_BMS_butterfly_or_NUTRI_taxon_scores_opened":True,
         "no_causal_wastewater_effect_estimated":True,
         "summary":summary}
    (a.output/"LOCAL16_POSTHOC_500M_1KM_GBIF_PLANT_RECORDING_OPPORTUNITY_QA.json").write_text(
      json.dumps(outcome,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (a.output/"LOCAL16_500M_1KM_PLANT_SOURCE_RECORDING_COUNTS.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print("GESUI_LOCAL16_FIXED_500M_1000M_SOURCE_SENSITIVITY="+json.dumps({
      "status":outcome["status"],
      "main_2020_22_500m":summary["radius_0.50km_period_2020-2022"],
      "main_2020_22_1km":summary["radius_1.00km_period_2020-2022"],
      "early_2017_19_500m":summary["radius_0.50km_period_2017-2019"],
      "early_2017_19_1km":summary["radius_1.00km_period_2017-2019"],
      "no_plant_trait_or_sewage_response_fitted":True
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geometry",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
