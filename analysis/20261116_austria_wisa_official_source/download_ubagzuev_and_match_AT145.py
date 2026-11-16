#!/usr/bin/env python3
"""Source-only official BMLUK GZÜV freshwater station GeoPackage download.

The download URL is read from the public UBA Nextcloud sharing PAGE'S EXPLICIT
public-page download menu (not guessed/private API). Verify file type and
source SHA, read genuine station ID/CRS. Compare geography to the original
AT145 EEA2022 annual-measured N/P genuine outfalls (from fixed PR116 source
artifact). This does NOT open 2022 observed chemistry, benthic EQR or NDVI.
"""
from __future__ import annotations
import argparse,base64,collections,csv,hashlib,io,json,math,sqlite3,zipfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlparse
import geopandas as gpd
import numpy as np
import requests
from shapely.geometry import Point

SOURCE_SHARE="https://docs.umweltbundesamt.at/s/SNJqq7TEprxJtpJ"
ORIGINAL_DATASET_ID="b054ce75-f0b1-449e-b11c-9603a629c9c5"
ORIGINAL_STATION_FILENAME="FWMST.gpkg"
SOURCE_OUTFALL_SCHEMA="GESUI_20261116_EEA2022_MEASURED_NP_ONE_FW_OUTFALL_DIRECTED_HYDRORIVERS_SOURCE_GEO_V1"
SCHEMA="GESUI_20261116_OFFICIAL_BMLUK_GZUEV_RIVER_GPKG_AT145_OUTFALL_SOURCE_LOCALITY_V1"
RADIUS_M=[500,1000,2000,5000]
MAX_ARCHIVE_BYTES=25_000_000

class HTML(HTMLParser):
    def __init__(self):super().__init__();self.hidden=[]
    def handle_starttag(self,tag,attrs):
        if tag=="input":self.hidden.append(dict(attrs))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def original_public_download(session,out):
    s=session.get(SOURCE_SHARE,timeout=(12,55))
    s.raise_for_status()
    if urlparse(s.url).hostname!="docs.umweltbundesamt.at":
        raise RuntimeError("Source authentic Austrian agency public share unexpectedly redirects")
    page=HTML();page.feed(s.text)
    filename=None;public_menu=None
    for h in page.hidden:
        if h.get("id")=="initial-state-files_sharing-filename":
            filename=json.loads(base64.b64decode(h["value"]).decode())
        if h.get("id")=="initial-state-core-public-page-menu":
            public_menu=json.loads(base64.b64decode(h["value"]).decode())
    if filename!=ORIGINAL_STATION_FILENAME or not isinstance(public_menu,list):
        raise RuntimeError("Original BMLUK source GZUEV station share has changed filename/declared menu")
    entry=next((q for q in public_menu if q.get("id")=="download"),None)
    if not entry or not entry.get("href","").startswith("https://docs.umweltbundesamt.at/"):
        raise RuntimeError("Original public share no public downloadable source file advertised")
    url=entry["href"]
    if SOURCE_SHARE.split("/s/")[1] not in url or "/public.php/dav/files/" not in url or "accept=zip" not in url:
        raise RuntimeError("Source public share declares unexpected non-file URL")
    r=session.get(url,timeout=(15,90),stream=True)
    if r.status_code!=200:raise RuntimeError("Actual public Austrian UBA original gpkg download HTTP "+str(r.status_code))
    total=0
    archive=out/"ORIGINAL_GZUEV_STATION_SOURCE_DO_NOT_REHOST.zip"
    with archive.open("wb") as dest:
        for chunk in r.iter_content(1024*1024):
            if not chunk:continue
            total+=len(chunk)
            if total>MAX_ARCHIVE_BYTES:
                raise RuntimeError("Original GZUEV source exceeds predeclared conservative 25MB file safety cap")
            dest.write(chunk)
    source=out/filename
    # Nextcloud public sharing of one file may honor the public "download"
    # menu but stream the genuine single-file GeoPackage WITHOUT wrapping
    # it in a ZIP, even when its advertised URL includes ?accept=zip.
    # Determine format from authentic file magic; NEVER use a fake ZIP.
    original_download_sha=sha(archive)
    rawmagic=archive.open("rb").read(16)
    if rawmagic.startswith(b"SQLite format 3"):
        archive.replace(source)
        original_delivery="SQLite GeoPackage direct, not ZIP (source URL retained)"
    elif zipfile.is_zipfile(archive):
        with zipfile.ZipFile(archive) as zz:
            candidates=[q for q in zz.infolist() if q.filename.endswith("/"+filename) or q.filename==filename]
            if len(candidates)!=1 or candidates[0].file_size>100_000_000:
                raise RuntimeError("Original published single GeoPackage member missing/oversized")
            with zz.open(candidates[0]) as inp,source.open("wb") as dest:
                while True:
                    chunk=inp.read(1024*1024)
                    if not chunk:break
                    dest.write(chunk)
        original_delivery="ZIP containing official GeoPackage"
    else:
        raise RuntimeError("Official public sharing download has neither genuine SQLite GeoPackage nor valid ZIP, source first16="+str(rawmagic))
    return source,{"original_dataset_id":ORIGINAL_DATASET_ID,
        "original_station_publisher":"BMLUK (Austrian government); served by Austrian Environment Agency",
        "original_public_catalogue_distribution":SOURCE_SHARE,
        "original_public_page_HTML_SHA256":hashlib.sha256(s.content).hexdigest(),
        "original_public_menu_download_href_verbatim":url,
        "original_public_download_SHA256":original_download_sha,
        "original_public_download_bytes":total,
        "original_source_response_delivery":original_delivery,
        "original_GZUEV_GeoPackage_source_SHA256":sha(source),
        "original_gpkg_uncompressed_source_bytes":source.stat().st_size}

def original_station_points(source):
    with sqlite3.connect(f"file:{source}?mode=ro",uri=True) as db:
        check=db.execute("PRAGMA integrity_check").fetchone()[0]
        contents=db.execute("SELECT table_name,data_type,srs_id FROM gpkg_contents").fetchall()
        geometries=db.execute("SELECT table_name,column_name,geometry_type_name,srs_id FROM gpkg_geometry_columns").fetchall()
    if check!="ok" or not contents:
        raise RuntimeError("Official source GZUEV GeoPackage SQLite validation not ok")
    layers=[q[0] for q in geometries if "POINT" in q[2].upper()]
    if len(layers)!=1:
        raise RuntimeError("GZUEV original source does not expose exactly one point station layer: "+str(geometries))
    gdf=gpd.read_file(source,layer=layers[0],engine="pyogrio")
    if len(gdf)<40 or gdf.crs is None or gdf.geometry.isna().any():
        raise RuntimeError("Original official GZUEV source has too few/invalid stations or unknown projection")
    x=gdf.to_crs("EPSG:4326")
    if not x.geometry.geom_type.isin(["Point"]).all():
        raise RuntimeError("Official GZUEV station is not a Point for all source features")
    if not ((x.geometry.x>=9)&(x.geometry.x<=18)&
            (x.geometry.y>=46)&(x.geometry.y<=50)).all():
        raise RuntimeError("Original station location outside Austrian official source extent")
    original_cols=[k for k in x.columns if k!="geometry"]
    # Station identifiers are source-listed, not derived from nearest-river names.
    # The genuine BMLUK station source uses field PK/PKINT rather than
    # an English/id field name. Never turn row order into a fake station ID.
    id_cols=[k for k in ("PK","PKINT") if k in original_cols]
    field_quality={}
    for k in original_cols:
        v=x[k]
        field_quality[k]={"nonmissing":int(v.notna().sum()),
          "distinct_nonmissing":int(v.nunique(dropna=True)),
          "original_source_examples":v.dropna().astype(str).head(4).tolist()}
    station_ids=None;idfield=None
    for k in id_cols:
        v=x[k]
        strings=v.astype(str).str.strip()
        if v.notna().sum()==len(x) and strings.nunique()==len(x) and (strings!="").all():
            station_ids=strings.tolist();idfield=k
            break
    if idfield is None:
        raise RuntimeError("Original GZÜV station PK and PKINT not sufficiently unique; refuse fabricated row-index station IDs")
    # BEGIN_DAT describes the station's registered first operation only; it
    # cannot prove any actual 2022 chemistry or macrozoobenthos sample.
    began=x["BEGIN_DAT"] if "BEGIN_DAT" in x else None
    began2022=None;beganvalid=None
    if began is not None:
        import pandas as pd
        parsed=pd.to_datetime(began,errors="coerce",utc=True)
        beganvalid=int(parsed.notna().sum())
        began2022=int((parsed<pd.Timestamp("2023-01-01T00:00:00Z")).sum())
    return x,station_ids,{
      "original_GeoPackage_issued_point_layers":geometries,
      "original_gpkg_contents":contents,
      "actual_station_layer":layers[0],
      "original_publisher_station_records":len(x),
      "original_native_CRS":str(gdf.crs),
      "original_source_fields":original_cols,
      "original_station_id_candidate_columns":id_cols,
      "source_original_fields_nonmissing_unique_and_examples":field_quality,
      "station_BEGIN_DAT_parseable_count":beganvalid,
      "station_BEGIN_DAT_at_or_before_2022_count_not_measurement_2022":began2022,
      "source_primary_unique_station_identifier_field":idfield,
      "station_coord_bounding_box_lon_lat":list(x.total_bounds),
      "source_duplicate_coordinate_count":int(len(x)-len(x.geometry.drop_duplicates()))}

def originals(folder):
    gate=json.loads((folder/"EEA2022_NP_SOURCEFIRST_LOCAL_DIRECTED_RIVER_GEOMETRY_GATE.json").read_text())
    p=folder/"EEA2022_ALL_SOURCE_QUALIFIED_UWWTD_NP_LOCAL_RIVER_SITES.csv"
    if (gate.get("schema")!=SOURCE_OUTFALL_SCHEMA or
       gate.get("original_all_true_receiving_outfall_geo_rows")!=282 or
       gate.get("original_all_source_valid_geography_CSV_SHA256")!=sha(p)):
        raise RuntimeError("Original true 282 site/year UWWTD/source river framework not original signed PR116 artifact")
    with p.open(encoding="utf-8-sig",newline="") as f:out=list(csv.DictReader(f))
    if len(out)!=282 or collections.Counter(q["EEA_original_countryCode"] for q in out)!={"AT":145,"DE":137}:
        raise RuntimeError("Original true receiving AT145 DE137 source parent source cohort changed")
    return [q for q in out if q["EEA_original_countryCode"]=="AT"],gate

def haversine_distance_vector(lat,lon,slat,slon):
    rad=math.pi/180
    la=(slat-lat)*rad;lo=(slon-lon)*rad
    z=np.sin(la/2)**2+np.cos(lat*rad)*np.cos(slat*rad)*np.sin(lo/2)**2
    return 12742017.6*np.arcsin(np.sqrt(np.minimum(1.,np.maximum(0.,z))))

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"gesui-source-authentic-public-GZUEV-station-academic-20261116/1.0"})
    gpkg,origin=original_public_download(session,args.output)
    stations,ids,qa=original_station_points(gpkg)
    source,orig=originals(args.outfalls)
    xs=np.asarray(stations.geometry.x,float);ys=np.asarray(stations.geometry.y,float)
    per=[]
    for q in sorted(source,key=lambda x:x["original_EEA2022_plant_id"]):
        lat=float(q["original_source_outfall_lat"])
        lon=float(q["original_source_outfall_lon"])
        d=haversine_distance_vector(lat,lon,ys,xs)
        order=np.argsort(d)
        record={"original_EEA_UWWTD_2022_plant":q["original_EEA2022_plant_id"],
                "original_source_report_year":2022,
                "original_annual_positive_reported_N_t":float(q["source_N_t_year_reported_measured_2022"]),
                "original_annual_positive_reported_P_t":float(q["source_P_t_year_reported_measured_2022"]),
                "source_outfall_hydrography_checked":True,
                "source_nearest_official_GZUEV_station_distance_m":float(d[order[0]]),
                "source_nearest_official_GZUEV_station_ID":ids[int(order[0])],
                "source_nearest_station_original_GZUEV_is_not_2022_sample":True}
        for radius in RADIUS_M:
            inrange=np.flatnonzero(d<=radius)
            record[f"source_official_GZUEV_station_location_within_{radius}m_n"]=len(inrange)
            record[f"source_official_GZUEV_station_id_within_{radius}m"]=[ids[i] for i in inrange[:50]]
        per.append(record)
    totals={str(radius):{
        "original_source_2022_NP_report_positive_outfalls_with_any_official_GZUEV_station_within_radius":sum(
             q[f"source_official_GZUEV_station_location_within_{radius}m_n"]>=1 for q in per),
        "unique_source_GZUEV_station_IDs_near_original_outfalls":len({
              z for q in per for z in q[f"source_official_GZUEV_station_id_within_{radius}m"]})}
        for radius in RADIUS_M}
    result={"schema":SCHEMA,
      "original_2022_EEA_nutrient_report_AT_true_outfall_n":len(source),
      "genuine_original_official_GZUEV_geodata_provenance":origin,
      "source_verified_original_GZUEV_station_geodata":qa,
      "original_2022_outfall_geography_sha":orig["original_all_source_valid_geography_CSV_SHA256"],
      "source_official_station_geo_proximity_by_radius":totals,
      "source_AT145_outfall_to_national_GZUEV_site_proximity":per,
      "these_station_locations_include_active_inactive_before2023_not_2022_observations":True,
      "public_2022_same_station_TNTP_and_MZB_not_read_or_verified":True,
      "spatial_proximity_is_not_directed_source_receiving_path":True,
      "no_ecological_effect_NDVI_benthos_or_sewage_association_fitted":True,
      "source_file_not_rehosted_in_public_Git_repository":True}
    (args.output/"BMLUK_GZUEV_RIVER_GPKG_SOURCE_AT145_LOCALITY_AUDIT.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_GZUEV_ORIGINAL_GPKG_AT145_SOURCE="+json.dumps({
        "source_official_GZUEV_file_SHA256":origin["original_GZUEV_GeoPackage_source_SHA256"],
        "publisher_station_count":qa["original_publisher_station_records"],
        "official_station_source_ID_field":qa["source_primary_unique_station_identifier_field"],
        "n_original_AT2022_EEA_source_outfalls":len(source),
        "counts_by_radius":totals,
        "geographic_not_same_year_benthos":True},ensure_ascii=False,sort_keys=True),flush=True)
    # Keep original ZIP/GPKG out of public GitHub Actions artifact to honor
    # catalog source license attribution and limit original third-party bytes.
    for v in (args.output/"ORIGINAL_GZUEV_STATION_SOURCE_DO_NOT_REHOST.zip",gpkg):
        v.unlink(missing_ok=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--outfalls",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
