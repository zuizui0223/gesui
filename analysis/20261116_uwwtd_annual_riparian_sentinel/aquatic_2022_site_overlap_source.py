#!/usr/bin/env python3
"""Source ONLY: genuine 2022 Austrian/German water and benthic station overlap
at 48 preselected EEA 2022 real-measured N/P source outlet systems.

No aquatic EQR, TP concentration, NDVI vegetation, GBIF plant or butterfly
response values opened. Monitoring-site existence ≠ exposure, direction
or source-specific sewage chemistry.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,json,math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
SRC_SCHEMA="GESUI_20261116_EEA2022_MEASURED_NP_ONE_FW_OUTFALL_DIRECTED_HYDRORIVERS_SOURCE_GEO_V1"
OUT_SCHEMA="GESUI_20261116_EEA2022_48_OUTFALL_REAL_AQUATIC_STATION_SOURCE_JOINT_OVERLAP_V1"
BQE="results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
WISE_BLOB="0f00351e7d907663589ea7a77b0591a9e1efaf36"
EXPECTED_GEOCODED={"DE":214,"AT":87}
EXPECTED_DATE_ELIGIBLE={"DE":242,"AT":100}
RADII_M=(500,1000,2000,5000)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()

def original_source(path):
    gate=json.loads((path/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").read_text())
    target=path/"EEA2022_SATELLITE_SOURCE_BLIND_INDEPENDENT_OUTFALL_PILOT.csv"
    if (gate.get("schema")!=SRC_SCHEMA or
        gate.get("source_status")!="PASS_SOURCE_GEOMETRY_FOR_STAC_METADATA_ONLY" or
        not gate.get("original_source_pilot_pass") or
        gate.get("pilot_source_CSV_SHA256")!=sha(target)):
        raise RuntimeError("Original source 48 satellite/outfall independent site cohort missing")
    with target.open(encoding="utf-8-sig",newline="") as f:rows=list(csv.DictReader(f))
    if len(rows)!=48 or collections.Counter(x["EEA_original_countryCode"] for x in rows)!={"DE":24,"AT":24}:
        raise RuntimeError("EEA original source 48/24/24 cohort changed")
    return rows,gate

def real_bqe(source):
    if blob(source)!=BQE_BLOB:
        raise RuntimeError("Main original 2022 actual EEA benthic site source SHA drift")
    sites=[]
    with source.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.reader(f)
        head=next(reader)
        need=("site_key","countryCode","lat","lon","HYRIV_ID","procedureClassificationSystem")
        if not set(need).issubset(head):
            raise RuntimeError("Original published river BQE site provenance fields missing")
        index={k:head.index(k) for k in need}
        n=0
        for row in reader:
            n+=1
            c=row[index["countryCode"]]
            if c not in ("AT","DE"):continue
            try:
                lat=float(row[index["lat"]]);lon=float(row[index["lon"]])
                h=int(row[index["HYRIV_ID"]])
            except (TypeError,ValueError):continue
            if not (math.isfinite(lat) and math.isfinite(lon)):continue
            sites.append({"id":row[index["site_key"]],"country":c,
               "lat":lat,"lon":lon,"HYRIV_ID":h,
               "method":row[index["procedureClassificationSystem"]]})
    if n!=2475 or collections.Counter(x["country"] for x in sites)!={"AT":68}:
        raise RuntimeError("Original 2022 EEA observed benthic BQE source station country support changed")
    return sites,n

def original_WQ(path):
    if blob(path)!=WISE_BLOB:
        raise RuntimeError("Original national 2022 WQ dated TP monitoring source manifest Gitblob changed")
    j=json.loads(path.read_text(encoding="utf-8-sig"))
    if (j.get("source_original_national_date_eligible_2022_TP_site_total")!=342 or
        j.get("source_original_date_eligible_DE_AT_by_country")!=EXPECTED_DATE_ELIGIBLE):
        raise RuntimeError("Original source TP sample year criteria or 342 date eligibility changed")
    rows=j["nationwide_original_source_2022_station_registry"]
    if len(rows)!=301:
        raise RuntimeError("Original national 342 date-eligible sites have 301 valid-geocoded station records; source registry changed")
    out=[]
    for q in rows:
        if not ("original_WISE6_WGS84_lat" in q and "original_WISE6_WGS84_lon" in q):continue
        try:
            lat=float(q["original_WISE6_WGS84_lat"])
            lon=float(q["original_WISE6_WGS84_lon"])
        except (ValueError,TypeError):continue
        if not (math.isfinite(lat) and math.isfinite(lon)):continue
        if q["countryCode"] not in ("AT","DE"):continue
        out.append({"id":q["original_country_scheme_site_SHA256"],
             "country":q["countryCode"],"lat":lat,"lon":lon,
             "source_dates":int(q["original_2022_TP_independent_sample_dates"])})
    if collections.Counter(q["country"] for q in out)!=EXPECTED_GEOCODED:
        raise RuntimeError("Original 301 official WQ source geocoded TP-date-eligible support changed")
    if len({q["id"] for q in out})!=len(out):
        raise RuntimeError("Original national TP source site hashes reused")
    return out,j

def haversine(lat1,lon1,lat2,lon2):
    p=math.pi/180
    a=(lat2-lat1)*p;b=(lon2-lon1)*p
    x=math.sin(a/2)**2+math.cos(lat1*p)*math.cos(lat2*p)*math.sin(b/2)**2
    return 12742017.6*math.asin(math.sqrt(min(1,max(0,x))))

def match(s,stations):
    lt=float(s["original_source_outfall_lat"])
    ln=float(s["original_source_outfall_lon"])
    cc=s["EEA_original_countryCode"]
    return sorted([(haversine(lt,ln,z["lat"],z["lon"]),z)
       for z in stations if z["country"]==cc],key=lambda x:x[0])

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    pilot,orig=original_source(a.source)
    bqe,bqerows=real_bqe(ROOT/BQE)
    wq,wqsource=original_WQ(a.water)
    outfall_summaries=[]
    for s in sorted(pilot,key=lambda z:z["original_EEA2022_plant_id"]):
        b=match(s,bqe);w=match(s,wq)
        reach0=int(s["source_GIS_original_river_snapped_reach"])
        reach1=int(s["source_true_graph_downstream_reach_id"])
        by_radius={}
        for radius in RADII_M:
            bb=[z for dist,z in b if dist<=radius]
            ww=[z for dist,z in w if dist<=radius]
            hydrology=[z for z in bb if z["HYRIV_ID"] in (reach0,reach1)]
            by_radius[str(radius)]={"EEA_benthic_source_sites_geodesic_radius":len(bb),
                "BQE_sites_on_same_original_source_reach_or_immediate_D1":len(hydrology),
                "2022_TP_source_station_date_valid_nearby":len(ww),
                "has_both_types_nearby_not_same_physical_station":bool(bb and ww),
                "source_monitoring_site_response_data_NOT_OPENED":True}
        outfall_summaries.append({"original_EEA_plant_id":s["original_EEA2022_plant_id"],
            "country":s["EEA_original_countryCode"],
            "original_receiver_HYRIV_ID":reach0,
            "source_immediate_downstream_HYRIV_ID":reach1,
            "nearest_2022_BQE_site_geodesic_m":float(b[0][0]) if b else None,
            "nearest_date_eligible_2022_TP_station_geodesic_m":float(w[0][0]) if w else None,
            "source_all_radius_monitoring_counts":by_radius,
            "source_BQE_EQR_or_TP_numeric_values_NOT_READ":True})
    totals={}
    for radius in RADII_M:
        f=str(radius);totals[f]={
          "outfalls_with_any_EEA2022_BQE_site":sum(q["source_all_radius_monitoring_counts"][f]["EEA_benthic_source_sites_geodesic_radius"]>0 for q in outfall_summaries),
          "outfalls_with_any_EEA2022_BQE_on_outfall_or_D1_reach":sum(q["source_all_radius_monitoring_counts"][f]["BQE_sites_on_same_original_source_reach_or_immediate_D1"]>0 for q in outfall_summaries),
          "outfalls_with_2022_TP_site":sum(q["source_all_radius_monitoring_counts"][f]["2022_TP_source_station_date_valid_nearby"]>0 for q in outfall_summaries),
          "outfalls_with_BQE_and_TP_sites_within_radius":sum(q["source_all_radius_monitoring_counts"][f]["has_both_types_nearby_not_same_physical_station"] for q in outfall_summaries),
          "countries_with_any_2022_BQE":{"AT":sum(q["source_all_radius_monitoring_counts"][f]["EEA_benthic_source_sites_geodesic_radius"]>0 for q in outfall_summaries if q["country"]=="AT"),
                                          "DE":sum(q["source_all_radius_monitoring_counts"][f]["EEA_benthic_source_sites_geodesic_radius"]>0 for q in outfall_summaries if q["country"]=="DE")}}
    # Original 2km joint same-country aquatic monitoring cannot pass if
    # Germany has no published EEA benthic stations in this original panel.
    joint=totals["2000"]["outfalls_with_BQE_and_TP_sites_within_radius"]
    countries=set(q["country"] for q in outfall_summaries
        if q["source_all_radius_monitoring_counts"]["2000"]["has_both_types_nearby_not_same_physical_station"])
    jointpass=(joint>=12 and len(countries)>=2 and
               all(sum(q["country"]==cc and q["source_all_radius_monitoring_counts"]["2000"]["has_both_types_nearby_not_same_physical_station"]
                       for q in outfall_summaries)>=3 for cc in ("AT","DE")))
    if jointpass:
        # Even both nearby categories do not prove actual hydrological
        # downstream-site overlap: therefore geographic support only.
        status="SOURCE_PROXIMITY_COUNTS_PASS_BUT_DOWNSTREAM_NETWORK_EXPOSURE_UNVERIFIED"
    else:
        status="STOP_REAL_TWO_COUNTRY_AQUATIC_BQE_TP_SOURCE_OVERLAP"
    result={"schema":OUT_SCHEMA,"status":status,
       "original_source_48_response_blind_measured_NP_outfalls":len(pilot),
       "original_year_2022_BQE_file_Gitblob":BQE_BLOB,
       "original_year_2022_benthic_BQE_source_total":bqerows,
       "original_BQE_source_2022_DE_AT_country_counts":dict(collections.Counter(x["country"] for x in bqe)),
       "original_WISE6_2022_date_eligible_TP_station_registry_Gitblob":WISE_BLOB,
       "original_2022_national_TP_station_date_eligible":len(wqsource["nationwide_original_source_2022_station_registry"]),
       "original_2022_national_TP_station_with_GPS":len(wq),
       "original_2022_national_TP_GPS_by_country":dict(collections.Counter(x["country"] for x in wq)),
       "predeclared_proximity_distances_m":RADII_M,
       "source_outfall_nearby_station_observations":totals,
       "source_per_outfall_opportunity":outfall_summaries,
       "joint_2km_both_aquatic_data_categories_geodesic_proximity_only":joint,
       "joint_2km_same_2022_wq_and_BQE_source_countries":sorted(countries),
       "independent_AQUATIC_LAND_joint_2country_gate_pass":bool(jointpass),
       "no_same_station_or_original_network_downstream_claim_from_distance":True,
       "EEA_BQE_main_2022_DE0_is_source_scope_not_biological_absence":True,
       "observed_TP_or_benthic_EQR_values_NOT_READ":True,
       "S2_12_source_pixel_two_year_original_GATE_STILL_STOP":True,
       "no_treatment_effect_or_cross_ecosystem_path_fitted":True}
    (a.output/"EEA2022_48_REAL_OUTFALL_TP_BQE_AQUATIC_SITE_SOURCE_OVERLAP.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_EEA2022_48_OUTFALL_AQUATIC_SITE_SOURCE_GATE="+json.dumps({
       "status":status,"original_BQE_DE_AT_source_sites":result["original_BQE_source_2022_DE_AT_country_counts"],
       "original_WQ_DE_AT_sites":result["original_2022_national_TP_GPS_by_country"],
       "source_48_outfall_proximity_counts":totals,
       "joint_2km_outfall_count":joint,
       "joint_2km_country_groups":sorted(countries),
       "no_aquatic_numeric_EQR_or_TP_or_terrestrial_effect":True
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--water",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
