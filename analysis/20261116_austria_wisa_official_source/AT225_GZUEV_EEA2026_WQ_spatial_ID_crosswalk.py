#!/usr/bin/env python3
"""Source-only original official Austrian GZÜV national PKINT versus genuine
EEA2026 Water Quality ICM original WISE6 station *spatial identity* crosswalk.

Original 225 source-selected national river sites near 117 positive original
2022 WWTP measured N/P outlets; compare with original 2026 EEA Water Quality
SpatialObjects known original station latitude/longitude. No water, biota
or satellite outcomes opened. Physical station co-location is not proven
receiving-water direction or sampled year.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,sys,zipfile
from pathlib import Path
import numpy as np
import requests
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261116_austria_wisa_official_source"))
from download_ubagzuev_and_match_AT145 import original_public_download,original_station_points

EEA_URL="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_SpatialObjects_DerivedData-csv.zip"
EEA_SIZE=2797036
EEA_SHA="aa505be6f85d06c675e0d0fcb148a9148f80162f8e9590740a887cc95a3e349d"
EEA_MEMBER="Waterbase_v2025_1_S_WISE6_SpatialObject_DerivedData.csv"
EEA_ROWS=90859
AT225_SHA="040d610b6172d165787775406ea50c0f73f34db1db3fcb6c8576898736146dde"
SCHEMA="GESUI_20261116_AT225_GZUEV_TO_EEA_WQ2026_ORIGINAL_STATION_COORDINATE_CROSSWALK_V1"
CUTOFFS=(5.,10.,25.,100.)
RADIUS=6371008.8
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sites(p):
    if sha(p)!=AT225_SHA:raise RuntimeError("Original 225 official PKINT geography source SHA changed")
    with p.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    if len(rows)!=225 or len({r["original_official_GZUEV_PKINT"] for r in rows})!=225:
        raise RuntimeError("Original genuine 225 government station ID rows changed")
    if len({r["source_2022_EEA_UWWTD_original_plant_ID"] for r in rows})!=117:
        raise RuntimeError("Original 117 2022 WWTP source IDs drift")
    return rows
def hav(lat,lon,la,lo):
    x=(np.asarray(la,float)-lat)*math.pi/180
    y=(np.asarray(lo,float)-lon)*math.pi/180
    v=np.sin(x/2)**2+np.cos(lat*math.pi/180)*np.cos(np.asarray(la,float)*math.pi/180)*np.sin(y/2)**2
    return 2*RADIUS*np.arcsin(np.sqrt(np.clip(v,0,1)))
def parse_original_wq(out):
    dest=out/"TEMP_EEA2026_WISE6_SPATIAL_DO_NOT_REHOST.zip"
    r=requests.get(EEA_URL,headers={"User-Agent":"gesui-AT225-EEA2026-station-geometry-source/20261116"},
             stream=True,timeout=(12,100))
    if r.status_code!=200 or not r.url.startswith("https://sdi.eea.europa.eu/"):
        raise RuntimeError("Original official public EEA WQ spatial source inaccessible")
    digest=hashlib.sha256();n=0
    with dest.open("wb") as f:
        for x in r.iter_content(1024*1024):
            if not x:continue
            n+=len(x);digest.update(x);f.write(x)
    if n!=EEA_SIZE or digest.hexdigest()!=EEA_SHA:
        raise RuntimeError("Original EEA WQ spatial archive SHA/size changed")
    info=collections.defaultdict(list)
    invalid=collections.Counter()
    with zipfile.ZipFile(dest) as z:
        if EEA_MEMBER not in z.namelist():raise RuntimeError("Official original spatial CSV member missing")
        reader=csv.DictReader(io.TextIOWrapper(z.open(EEA_MEMBER),encoding="utf-8-sig",newline=""))
        expected=("countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme","lat","lon")
        if any(k not in reader.fieldnames for k in expected):
            raise RuntimeError("Actual original EEA WQ spatial site source fields changed")
        count=0
        for row in reader:
            count+=1
            if row["countryCode"]!="AT":continue
            key=(row["countryCode"],row["monitoringSiteIdentifier"],row["monitoringSiteIdentifierScheme"])
            if not all(key):invalid["incomplete_original_source_key"]+=1;continue
            try:lat=float(row["lat"]);lon=float(row["lon"])
            except (TypeError,ValueError):
                invalid["no_numeric_original_coord"]+=1;continue
            if not (math.isfinite(lat) and math.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180 and
                   45<=lat<=51 and 8<=lon<=19):
                invalid["invalid_original_AT_coordinate"]+=1;continue
            info[key].append((lat,lon))
        if count!=EEA_ROWS:raise RuntimeError("Published official 2026 WQ spatial source row count changed")
    dest.unlink(missing_ok=True)
    positions={}
    ambiguous=[]
    for key,coords in info.items():
        base=coords[0]
        dif=hav(base[0],base[1],np.array([x[0] for x in coords]),np.array([x[1] for x in coords]))
        if np.max(dif)>50:
            ambiguous.append(key);continue
        positions[key]=base
    return positions,{"original_EEA_WQ_spatial_archive_SHA256":EEA_SHA,
        "original_EEA_WQ_spatial_source_rows":EEA_ROWS,
        "source_country_AT_unique_unambiguous_source_station_keys_with_geometry":len(positions),
        "source_original_AT_ambiguous_coords_distinct_keys":len(ambiguous),
        "source_invalid_AT_coordinate_and_original_key_counts":dict(invalid),
        "source_no_2022_any_numeric_concentration_or_benthic_outcome_read":True}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    target=sites(a.crosswalk)
    sess=requests.Session();sess.headers.update({"User-Agent":"gesui-original-GZUEV-to-WISE6-public-station/20261116"})
    gpkg,gov=original_public_download(sess,a.output)
    gdf,govids,qa=original_station_points(gpkg)
    if len(gdf)!=3788 or qa["source_primary_unique_station_identifier_field"]!="PK" or "PKINT" not in gdf:
        raise RuntimeError("Authentic BMLUK source 3788 unique PK national station records invalid")
    byid={str(k):i for i,k in enumerate(gdf["PKINT"])}
    if len(byid)!=3788:raise RuntimeError("Original PKINT mapping duplicates")
    wq,source=parse_original_wq(a.output)
    if not wq:raise RuntimeError("Original actual Austria EEA WQ station coordinate source empty")
    waterkeys=sorted(wq)
    wlat=np.asarray([wq[q][0] for q in waterkeys],float)
    wlon=np.asarray([wq[q][1] for q in waterkeys],float)
    results=[]
    for r in sorted(target,key=lambda q:q["original_official_GZUEV_PKINT"]):
        pk=r["original_official_GZUEV_PKINT"]
        if pk not in byid:raise RuntimeError("Original eligible government GZUEV station not in original GPKG")
        i=byid[pk]
        orig=gdf.geometry.iloc[i];lat=float(orig.y);lon=float(orig.x)
        distance=hav(lat,lon,wlat,wlon)
        order=np.argsort(distance)
        n0=int(order[0]);n1=int(order[1])
        p0=float(distance[n0]);p1=float(distance[n1])
        nearest=waterkeys[n0]
        row={"original_gov_PKINT":pk,
             "original_source_reportedNP_2022_outfall":r["source_2022_EEA_UWWTD_original_plant_ID"],
             "source_GZUEV_official_station_WGS84_lat":lat,"source_GZUEV_official_station_WGS84_lon":lon,
             "nearest_original_EEA_WQ_station_key":nearest[1],
             "nearest_original_EEA_WQ_station_ID_scheme":nearest[2],
             "nearest_same_country_official_WQ_site_geodesic_m":p0,
             "next_nearest_official_WQ_site_geodesic_m":p1,
             "source_10m_unambiguous_nearest":bool(p0<=10 and p1>10),
             "source_10m_coordinate_only_not_verified_station_semantics":True,
             "numeric_2022_TNTP_or_MZB_not_read":True}
        for cutoff in CUTOFFS:
            row[f"source_any_EEA_WQ_site_within_{int(cutoff)}m"]=bool(p0<=cutoff)
        results.append(row)
    # A second source target cannot be allowed to claim identity with the
    # SAME original WQ site merely because nearby point coordinates coincide.
    used=collections.Counter(q["nearest_original_EEA_WQ_station_key"]+"|"+q["nearest_original_EEA_WQ_station_ID_scheme"]
          for q in results if q["source_10m_unambiguous_nearest"])
    for q in results:
        kk=q["nearest_original_EEA_WQ_station_key"]+"|"+q["nearest_original_EEA_WQ_station_ID_scheme"]
        q["source_10m_primary_unique_original_EEA_WQ_station_match"]=bool(
          q["source_10m_unambiguous_nearest"] and used[kk]==1)
    counts={}
    for dist in CUTOFFS:
        count=sum(q[f"source_any_EEA_WQ_site_within_{int(dist)}m"] for q in results)
        fac={q["original_source_reportedNP_2022_outfall"] for q in results
             if q[f"source_any_EEA_WQ_site_within_{int(dist)}m"]}
        counts[str(int(dist))]={"national_original_PKINT_sites_with_any_WQ_source_site":count,
                            "distinct_EEA2022_report_positive_outfalls_with_near_WQ_source_site":len(fac)}
    strict=[q for q in results if q["source_10m_primary_unique_original_EEA_WQ_station_match"]]
    out=a.output/"AT225_GZUEV_ORIGINAL_WISE6_WQ_SPATIAL_ID_AND_COORD_CROSSWALK.csv"
    with out.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    report={"schema":SCHEMA,"original_government_225_source_target_SHA":AT225_SHA,
        "authentic_BMLUK_original_government_FWMST_GPKG_SHA256":gov["original_GZUEV_GeoPackage_source_SHA256"],
        "original_government_n_stations":len(gdf),
        "source_original_WISE6_WQ_spatial":source,
        "official_2022_report_positive_outfall_source_geography_station_n":len(results),
        "fixed_source_coincidence_cutoffs_metres":list(CUTOFFS),
        "source_GZUEV_to_EEA_WQ_different_identifier_spatial_nearest_counts":counts,
        "source_primary_10m_unambiguous_unique_WQ_station_crosswalk_count":len(strict),
        "source_primary_unique_WQ_station_id_and_scheme_pairs":len({(q["nearest_original_EEA_WQ_station_key"],q["nearest_original_EEA_WQ_station_ID_scheme"]) for q in strict}),
        "source_225_derived_PKINT_WQ_station_candidate_crosswalk_SHA256":sha(out),
        "genuine_2022_WQ_TNTP_same_station_samples_UNVERIFIED":True,
        "source_site_coincidence_not_study_year_or_nutrient_load_delivery":True,
        "all_true_same_downstream_reach_or_station_name_UNVERIFIED":True,
        "no_EQR_NDVI_plant_butterfly_ecological_coefficients_fitted":True}
    (a.output/"AT225_GZUEV_OFFICIAL_EEA_WQ_SPATIAL_ID_CROSSWALK_SOURCE_GATE.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    for file in (gpkg,a.output/"ORIGINAL_GZUEV_STATION_SOURCE_DO_NOT_REHOST.zip"):
        file.unlink(missing_ok=True)
    print("GESUI_AT225_GZUEV_EEA_WQ_OFFICIAL_COORDINATE_SOURCE="+json.dumps({
       "source_gov_stations":len(gdf),"source_original_AT_WQ_station_keys":len(wq),
       "n_target_PKINT":len(results),
       "coincidence_counts":counts,
       "unambiguous_unique_10m":len(strict),
       "different_WQ_code_from_national_PKINT_n":sum(
         q["nearest_original_EEA_WQ_station_key"]!=q["original_gov_PKINT"] for q in strict),
       "sampled_2022_NP_not_checked":True},sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--crosswalk",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
