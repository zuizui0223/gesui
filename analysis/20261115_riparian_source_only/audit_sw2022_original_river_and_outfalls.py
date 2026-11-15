#!/usr/bin/env python3
"""Immutable SeBMS 2022 point-site vs ORIGINAL EEA 2022 outfall and HydroRIVERS geometry.

SOURCE / FEASIBILITY ONLY. No species, occurrence extensions, HNS, benthic responses,
or outcome coefficients. Never loosen pre-frozen primary 2022 / 8-cell eligibility
from this output. All derived geometries label uncertain approximations.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,pathlib,urllib.request,zipfile,time
from collections import Counter
import requests
import numpy as np
import pandas as pd
import geopandas as gpd

SEBMS_SITE_GIT_SHA = "a56d143a3474027092a438399e6c2e59066d911f"
PR21_OLD_SITES_GIT_SHA = "bbc1c423dd95eeddda224d8326caba13e61af06f"
FROZEN_PROTOCOL_GIT_SHA = "b9ce415c36d398eee77d8318a20aef4873e0bf44"
RIVERS_URL = "https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_eu_shp.zip"
RIVERS_SHA256 = "500da7d36ceee0aa4c82dd625f5c53d80ba2f7ac14ad022cd3a6601bf55068e5"
REST="https://water.discomap.eea.europa.eu/arcgis/rest/services/WISE_UWWTD"
PLANT=REST+"/TreatmentPlants_2022/MapServer/0"
DCP=REST+"/DischargePoints_2022/MapServer/0"
REPORTER="SE369"
COORD_MAX=100
SNAP_SOURCE_M=500
RIVER_NEAR_M=200
RIVER_INLAND_MIN_M=500
RIVER_INLAND_MAX_M=1500
SITE_FOOTPRINT_M=25
NETWORK_UPDOWN_MIN_KM=1.0
NETWORK_UPDOWN_MAX_KM=5.0
SE_EXTENDED_WGS_BBOX=(7.0,53.5,27.0,71.2)
MAX_ARCHIVE_BYTES=150*1024*1024

def blobsha(b):
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()

def check_git(p,sha):
    data=p.read_bytes()
    got=blobsha(data)
    if got!=sha:raise ValueError(f"Immutable original Git blob SHA changed {p}: {got}")
    return data

def robust_get(url,params):
    last=None
    for i in range(3):
        try:
            r=requests.get(url,params=params,timeout=80,headers={"User-Agent":"gesui-20261115-sourceonly-no-HNS-Swedish-geom"})
            r.raise_for_status()
            obj=r.json()
            if "error" in obj:raise RuntimeError(str(obj["error"])[:350])
            return obj,hashlib.sha256(r.content).hexdigest()
        except Exception as e:
            last=e
            time.sleep(1+i)
    raise RuntimeError(f"EEA original 2022 query failed {url}: {last}")

def source_query(url,where,fields):
    check,sha=robust_get(url+"/query",{"f":"json","where":where,"returnCountOnly":"true"})
    n=int(check["count"])
    if not 0<n<5000:raise ValueError(f"Unexpected national layer source subset cardinality: {n}")
    found=[]
    receipts=[{"type":"count","sha256":sha,"n":n}]
    for offset in range(0,n,1000):
        x,s=robust_get(url+"/query",{"f":"json","where":where,"returnGeometry":"false",
           "outFields":",".join(fields),"orderByFields":"OBJECTID ASC",
           "resultOffset":offset,"resultRecordCount":1000})
        pts=[row["attributes"] for row in x["features"]]
        if not pts:raise ValueError("Empty EEA paged results")
        found+=pts
        receipts.append({"type":"page","offset":offset,"n":len(pts),"sha256":s})
    d=pd.DataFrame.from_records(found)
    if len(d)!=n or d.OBJECTID.duplicated().any() or set(fields)-set(d.columns):
        raise ValueError("EEA source result paging/fields incomplete")
    return d,receipts

def source_outfalls():
    plant_fields=["OBJECTID","repCode","uwwCode","uwwNDischargeMeasured",
                  "uwwPDischargeMeasured","uwwState"]
    dcp_fields=["OBJECTID","repCode","uwwCode","dcpCode","dcpLatitude",
                "dcpLongitude","dcpWaterBodyType","dcpState"]
    plants,r1=source_query(PLANT,
        "repCode = 'SE369' AND uwwNDischargeMeasured > 0 AND uwwPDischargeMeasured > 0",
        plant_fields)
    dcp,r2=source_query(DCP,"repCode = 'SE369'",dcp_fields)
    for table in (plants,dcp):
        table["join_key"]=(table.repCode.astype(str).str.strip().str.upper()
                           +"|"+table.uwwCode.astype(str).str.strip())
    pc=plants.join_key.value_counts()
    dc=dcp.join_key.value_counts()
    plants=plants[plants.join_key.map(pc).eq(1)].copy()
    plants=plants[plants.join_key.map(dc).fillna(0).eq(1)].copy()
    sites=plants.merge(dcp,on="join_key",validate="one_to_one",
                       suffixes=("_plant","_outfall"))
    for field in ["dcpLatitude","dcpLongitude","uwwNDischargeMeasured","uwwPDischargeMeasured"]:
        sites[field]=pd.to_numeric(sites[field],errors="coerce")
    sites=sites[sites.dcpWaterBodyType.astype(str).str.strip().str.upper().eq("FW")]
    sites=sites[sites.dcpLatitude.between(44,72)&sites.dcpLongitude.between(-12,27)]
    sites=sites.dropna(subset=["dcpLatitude","dcpLongitude"])
    # Previous independently audited PR21 original EEA2022 14 single-FW,
    # both-positive registered reported outfalls. Abort on drift, do not redefine.
    if len(sites)!=14:
        raise ValueError("Original 2022 Swedish both-positive N/P single-FW outfalls drifted: "+str(len(sites)))
    out=[]
    for x in sites.itertuples(index=False):
        out.append({"outfall_sha256":hashlib.sha256(str(x.join_key).encode()).hexdigest(),
                    "latitude":float(x.dcpLatitude),"longitude":float(x.dcpLongitude),
                    "reportedN_positive":bool(x.uwwNDischargeMeasured>0),
                    "reportedP_positive":bool(x.uwwPDischargeMeasured>0)})
    return out,{"plants_from_SE369":int(len(plants)),
                "outfalls_all_reported_SE369":int(len(dcp)),
                "qualified_single_freshwater_outfalls":len(out),
                "EEA2022_source_url":PLANT,
                "original_EEA_api_receipts":{"plants":r1,"discharge":r2},
                "year2022_EEA_layer_not_measured_receiver_dose":True,
                "individual_annual_NP_year_not_independently_verified":True}

def hav(lat1,lon1,lat2,lon2):
    p=math.pi/180
    dlat=(lat2-lat1)*p;dlon=(lon2-lon1)*p
    q=math.sin(dlat/2)**2+math.cos(lat1*p)*math.cos(lat2*p)*math.sin(dlon/2)**2
    return 12742*math.asin(min(1,math.sqrt(max(0,q))))

def old_distance_identity(current,old,src):
    old_lookup={(f"{float(x['lat']):.5f}",f"{float(x['lon']):.5f}"):x
                for x in old}
    matched=0;errors=[];matched_site_ids=set()
    for site in current:
        key=(f"{float(site['latitude']):.5f}",f"{float(site['longitude']):.5f}")
        other=old_lookup.get(key)
        if other is None:continue
        matched+=1
        matched_site_ids.add(site["site_location_sha256_16"])
        calc=min(hav(float(site["latitude"]),float(site["longitude"]),
                     x["latitude"],x["longitude"]) for x in src)
        old_d=float(other["nearest_qualified_SWEDEN_2022_FW_WWTP_km"])
        errors.append(abs(calc-old_d))
    if matched<220 or max(errors)>0.00008:
        raise ValueError(f"Original PR21 source positions/distances not recovered; matched={matched},maxdiff={max(errors,default=None)}")
    return {"n_current_sites":len(current),"n_matched_previous_repeat2_geographies":matched,
            "n_new_or_unmatched_geographies":len(current)-matched,
            "old_nearest_source_haversine_max_absolute_difference_km":max(errors),
            "source_analysis_derived_without_any_butterfly_outcome":True}

def download_river(temp):
    target=temp/"HydroRIVERS_v10_eu_shp.zip"
    if not target.exists():
        last=None
        for attempt in range(3):
            try:
                req=urllib.request.Request(RIVERS_URL,headers={
                    "User-Agent":"gesui-sebms-pointsource-river-original-geometry-20261115",
                    "Accept":"application/zip,*/*"})
                with urllib.request.urlopen(req,timeout=90) as r, target.open("wb") as w:
                    total=0
                    while True:
                        b=r.read(1024*1024)
                        if not b:break
                        total+=len(b)
                        if total>MAX_ARCHIVE_BYTES:raise RuntimeError("River archive too large")
                        w.write(b)
                break
            except Exception as e:
                last=e
                target.unlink(missing_ok=True)
                time.sleep(2+attempt)
        if not target.exists():raise RuntimeError("Original HydroRIVERS publisher inaccessible: "+str(last))
    digest=hashlib.sha256(target.read_bytes()).hexdigest()
    if digest!=RIVERS_SHA256:raise ValueError("Authentic HydroRIVERS EU original ZIP SHA changed: "+digest)
    with zipfile.ZipFile(target) as z:
        files=[x for x in z.namelist() if x.lower().endswith(".shp")]
        if len(files)!=1 or "hydrorivers_v10_eu" not in files[0].lower():
            raise ValueError("Original HydroRIVERS EU shapefile member missing or ambiguous")
        base=files[0][:-4]
        members=[n for n in z.namelist() if n.startswith(base) and
                 n.lower().endswith((".shp",".shx",".dbf",".prj"))]
        if len(members)<3:raise ValueError("Original spatial source missing native sidecars")
        for m in members:z.extract(m,temp)
    return temp/files[0],{"original_zip_url":RIVERS_URL,
        "original_zip_SHA256":digest,"original_geodata_file":files[0],
        "original_geom_licensed_not_republished":True}

def near_river_snap(points,river,keep_distance_limit=1700):
    data=pd.DataFrame(points)
    geodata=gpd.GeoDataFrame(data,
        geometry=gpd.points_from_xy(data.longitude,data.latitude),crs="EPSG:4326").to_crs("EPSG:3035")
    g=river[["HYRIV_ID","geometry"]].to_crs("EPSG:3035")
    spatial=gpd.sjoin_nearest(geodata,g,how="left",max_distance=keep_distance_limit,
                             distance_col="distance_m")
    spatial=(spatial.sort_values(["source_index","distance_m","HYRIV_ID"],
                  na_position="last").drop_duplicates("source_index"))
    out={}
    for r in spatial.itertuples(index=False):
        if not hasattr(r,"source_index"):raise RuntimeError("Missing original point ID after spatial join")
        out[str(r.source_index)]={
          "HYRIV_ID":int(r.HYRIV_ID) if pd.notna(r.HYRIV_ID) else None,
          "nearest_river_m":round(float(r.distance_m),3) if pd.notna(r.distance_m) else None}
    if len(out)!=len(points):raise RuntimeError("Nearest river point merge cardinality failed")
    return out

def network_walk(a,b,graph,maxlen=8,max_hops=180):
    if a is None or b is None or a==b:return None
    visited=set();dist=0;current=a
    for _ in range(max_hops):
        if current in visited or current not in graph:return None
        visited.add(current)
        next_,length=graph[current]
        if not math.isfinite(length) or length<0:return None
        dist+=length
        if dist>maxlen:return None
        if next_==b:return round(dist,6)
        if not next_ or next_ not in graph:return None
        current=next_
    return None

def run(args):
    if blobsha(args.freeze.read_bytes())!=FROZEN_PROTOCOL_GIT_SHA:
        raise ValueError("Original frozen v1.0 protocol SHA mismatch; forbidden accidental modification")
    df_raw=check_git(args.sites,SEBMS_SITE_GIT_SHA)
    old_raw=check_git(args.legacy,PR21_OLD_SITES_GIT_SHA)
    sites=list(csv.DictReader(io.StringIO(df_raw.decode("utf-8-sig"))))
    old=list(csv.DictReader(io.StringIO(old_raw.decode("utf-8-sig"))))
    if len(sites)!=264 or sum(s["five_or_more_dates"]=="True" for s in sites)!=127:
        raise ValueError("Original 264 sites / 127 repeated site source drift")
    outfalls,eea=source_outfalls()
    geodesic=old_distance_identity(sites,old,outfalls)
    args.tmpdir.mkdir(parents=True,exist_ok=True)
    shp,provenance=download_river(args.tmpdir)
    # Original authentic geometry. Load Swedish+edge-country reach subset;
    # lack of full graph beyond bbox is flagged, NEVER reclassified as source-zero.
    river=gpd.read_file(shp,bbox=SE_EXTENDED_WGS_BBOX,engine="pyogrio")
    river=river.rename(columns={x:x.upper() for x in river.columns if x!="geometry"})
    fields={"HYRIV_ID","NEXT_DOWN","LENGTH_KM","geometry"}
    if fields-set(river.columns):raise ValueError("Original EU River georeach schema mismatch")
    if river.crs is None:raise ValueError("Original river CRS unknown")
    river=river.to_crs("EPSG:4326")
    river=river[["HYRIV_ID","NEXT_DOWN","LENGTH_KM","geometry"]].copy()
    for f in ["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]:river[f]=pd.to_numeric(river[f],errors="coerce")
    river=river.dropna(subset=["HYRIV_ID","LENGTH_KM"]).sort_values("HYRIV_ID").drop_duplicates("HYRIV_ID")
    if len(river)<10000:raise ValueError("Original Swedish regional river subset implausibly sparse")
    graph={int(r.HYRIV_ID):(int(r.NEXT_DOWN) if pd.notna(r.NEXT_DOWN) else None,
                             float(r.LENGTH_KM))
           for r in river[["HYRIV_ID","NEXT_DOWN","LENGTH_KM"]].itertuples(index=False)}
    sites_in=[{"source_index":s["site_location_sha256_16"],"latitude":float(s["latitude"]),
               "longitude":float(s["longitude"])} for s in sites]
    outfall_in=[{"source_index":q["outfall_sha256"],"latitude":q["latitude"],
                 "longitude":q["longitude"]} for q in outfalls]
    sr=near_river_snap(sites_in,river,keep_distance_limit=1700)
    wr=near_river_snap(outfall_in,river,keep_distance_limit=2000)
    snapped_sources=[x for x in outfalls if wr[x["outfall_sha256"]]["HYRIV_ID"] is not None
                     and wr[x["outfall_sha256"]]["nearest_river_m"]<=SNAP_SOURCE_M]
    if len(snapped_sources)>len(outfalls):raise RuntimeError("Outfall snapping duplicate")
    selected=[];rel=[];category=Counter()
    for s in sites:
        siteid=s["site_location_sha256_16"]
        sm=sr[siteid]
        repeat5=s["five_or_more_dates"]=="True"
        dist=sm["nearest_river_m"]
        if dist is None:
            zone="NO_VALID_RIVER_WITHIN_1700M"
        elif dist+SITE_FOOTPRINT_M<=RIVER_NEAR_M:
            zone="NEAR_0_TO_200M_ENTIRE_25M_RADIUS"
        elif dist-SITE_FOOTPRINT_M>=RIVER_INLAND_MIN_M and dist+SITE_FOOTPRINT_M<=RIVER_INLAND_MAX_M:
            zone="INLAND_500_TO_1500M_ENTIRE_25M_RADIUS"
        else:
            zone="BOUNDARY_OR_INTERMEDIATE_DISTANCE"
        near_source=min((hav(float(s["latitude"]),float(s["longitude"]),
                             o["latitude"],o["longitude"]) for o in outfalls),default=None)
        row={"site_location_sha256_16":siteid,
             "original_2022_may_sep_distinct_visits":int(s["n_distinct_2022_may_sep_visit_dates"]),
             "point_count_5dates":repeat5,
             "original_published_uncertainty_m":float(s["published_uncertainty_m"]),
             "nearest_original_river_m":dist,
             "nearest_original_river_reach_HYRIV_ID":sm["HYRIV_ID"],
             "frozen_distance_band_preliminary":zone,
             "nearest_2022_positive_reported_singleFW_outfall_km":round(near_source,6) if near_source is not None else None,
             "no_river_source_group_assigned_to_HNS":True}
        selected.append(row)
        if repeat5:category[zone]+=1
        if not repeat5 or near_source is None or near_source>7.5 or sm["HYRIV_ID"] is None:
            continue
        for o in snapped_sources:
            d=hav(float(s["latitude"]),float(s["longitude"]),o["latitude"],o["longitude"])
            if d>7.5:continue
            outid=o["outfall_sha256"]
            reach=wr[outid]["HYRIV_ID"]
            if reach==sm["HYRIV_ID"]:
                label="SAME_REACH_ORDER_UNKNOWN"
                downstream=upstream=None
            else:
                downstream=network_walk(reach,sm["HYRIV_ID"],graph)
                upstream=network_walk(sm["HYRIV_ID"],reach,graph)
                label=("SOURCE_TO_SITE_DOWNSTREAM_NETWORK_1TO5KM_APPROX" if downstream is not None and 1<=downstream<=5
                       else "SITE_UPSTREAM_OF_SOURCE_NETWORK_1TO5KM_APPROX" if upstream is not None and 1<=upstream<=5
                       else "ROUTE_OUTSIDE_1TO5_OR_UNRESOLVED")
            rel.append({"site_hash":siteid,"outfall_sha256":outid,
                        "site_2022_5visit_source_only":True,
                        "geom_river_band":zone,
                        "original_network_direction_class":label,
                        "site_geodesic_outfall_km":round(d,6),
                        "source_original_river_snap_m":wr[outid]["nearest_river_m"],
                        "network_full_reach_km_downstream_approx":downstream,
                        "network_full_reach_km_upstream_approx":upstream,
                        "same_reach_flow_order_not_resolved":reach==sm["HYRIV_ID"],
                        "no_biological_outcomes_read":True})
    if len(selected)!=264:raise RuntimeError("Original site count changed")
    countdir=Counter(x["original_network_direction_class"] for x in rel)
    results={
      "schema":"GESUI_SW_SEBMS_2022_POINTS_EEA_SE369_HYDRORIVERS_AUTHENTIC_DIRECTION_SOURCE_ONLY_V1",
      "source_sha":{"frozen_protocol_git_blob":FROZEN_PROTOCOL_GIT_SHA,
        "Swedish_2022_points_git_blob":SEBMS_SITE_GIT_SHA,
        "old_PR21_2022_EEA_proximity_site_table_git_blob":PR21_OLD_SITES_GIT_SHA},
      "source_original_eea2022":eea,
      "source_original_rivers":provenance,
      "older_PR21_nearest_14_2022_outfall_distance_independent_replay":geodesic,
      "study":{"n_2022_point_site_locations":264,"n_point_sites_with_5_dates":127,
        "source_reported_singleFW_positive_NP_outfalls":len(outfalls),
        "source_outfalls_with_native_river_snap_le500m":len(snapped_sources),
        "original_authentic_river_reach_count_expanded_Sweden":len(river),
        "true_native_river_order_not_assumed_polyline_vertex":True,
        "nearest_river_2022_5date_site_distribution":dict(category),
        "source_site_direction_class_counts":dict(countdir),
        "n_2022_5visit_site_source_relations_within_7p5km_and_source_snap":len(rel),
        "n_distinct_2022_5visit_sites_source_relation_within_7p5km":len({x["site_hash"] for x in rel}),
        "n_distinct_registered_outfalls_with_any_relation":len({x["outfall_sha256"] for x in rel}),
        "n_distinct_registered_outfalls_with_downstream1to5km_candidate":len({x["outfall_sha256"] for x in rel if x["original_network_direction_class"]=="SOURCE_TO_SITE_DOWNSTREAM_NETWORK_1TO5KM_APPROX"})},
      "spatial_scope":"Preliminary source only. Approximate NEXT_DOWN full reach length, original 25m event footprints, geographic bbox subset, no actual effluent dilution.",
      "not_yet_source_reference_matched_8cells":True,
      "H7_H6_parent_assignment_NOT_PERFORMED":True,
      "site_near_river_does_NOT_prove_larval_host_or_adult_puddling":True,
      "registered2022_outfall_not_confirmed2022_receiving_water_concentration":True,
      "no_observed_species_counts_or_HNS_read":True,
      "biological_primary_run":False,
      "frozen_24pair_3country_primary_stage":"STOP_NOT_QUALIFIED",
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    with args.siteoutput.open("w",newline="",encoding="utf-8") as f:
        wrt=csv.DictWriter(f,fieldnames=list(selected[0]));wrt.writeheader();wrt.writerows(selected)
    # Summary only per-source pair, public source SHA and location SHA; no raw HydroRIVERS shape redistribution.
    if rel:
        with args.relations.open("w",newline="",encoding="utf-8") as f:
            wrt=csv.DictWriter(f,fieldnames=list(rel[0]));wrt.writeheader();wrt.writerows(rel)
    print("GESUI_SW_HYDRO2022_SOURCE_ONLY="+json.dumps(results,
          sort_keys=True,separators=(",",":"),ensure_ascii=False),flush=True)
    return results

def main():
    p=argparse.ArgumentParser()
    for k in ["sites","legacy","freeze","output","siteoutput","relations","tmpdir"]:
        p.add_argument("--"+k,type=pathlib.Path,required=True)
    run(p.parse_args())

if __name__=="__main__":main()
