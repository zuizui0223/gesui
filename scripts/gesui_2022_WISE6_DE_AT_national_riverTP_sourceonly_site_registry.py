#!/usr/bin/env python3
"""Independent NATIONAL DE/AT 2022 WISE6 river-TP monitoring-site metadata.

Re-stream original 102,569,057 EEA Waterbase WISE6 DisaggregatedData rows,
strictly do NOT parse numeric TP concentrations, only 2022 country/site/date
source metadata; restore all national date-eligible stations, not restricted
to originally 42 Lycaenidae/photo source-positive 1degree cells.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from shapely.geometry import Point
import gesui_2026_eea_WISE6_2022_individual_TP_fullstream_sourceonly as prior
import gesui_eea2026_wise6_original_2022_TP_annual_station_coverage as water

EXPECTED_ARCHIVE_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
EXPECTED_SITE_NUM={"DE":242,"AT":100}
EXPECTED_OLD_SITE=187
MIN_DATES=4
MIN_MONTHS=3

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-187-stations",type=Path,required=True)
    p.add_argument("--tmpdir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    old=json.loads(args.source_187_stations.read_text())
    if not old.get("no_original_numeric_river_TP_concentrations_or_LOQ_values_parsed"):
        raise ValueError("Original pre-numeric 187 site source compromised")
    old_sites=set()
    for row in old["original_42_country_1degree_photo_cells_and_WISE6_TP_stations"]:
        for q in row["original_non_numeric_WISE6_2022_station_sampling_support"]:
            old_sites.add(q["monitoring_site_ID_with_scheme_sha256"])
    if len(old_sites)!=EXPECTED_OLD_SITE:raise ValueError("Original source-only 187 site hash baseline changed")
    args.tmpdir.mkdir(parents=True,exist_ok=True)
    archive=args.tmpdir/"EEA2026_original_public_WISE6_DA_2022_nationalTP_source.zip"
    spatialzip=args.tmpdir/"EEA2026_original_public_WISE6_spatial_source.zip"
    source=prior.original_archive(archive)
    if source["archive_sha256"]!=EXPECTED_ARCHIVE_SHA:
        raise ValueError("Original EEA2026 WISE6 2022 individual ZIP actual SHA changed")
    spatial=water.load(water.URL_2026_SPATIAL,spatialzip,water.EXPECTED_SPATIAL_SHA256,6000000)
    positions,spatialQC=water.spatial_source(spatialzip)
    polygons=water.get_geometries()
    sites,originalQC=prior.summarize_original_DE_AT_TP_sample_sites(archive)
    native=Counter();quality=Counter()
    output=[]
    for (country,sid,scheme),data in sites.items():
        nd=len(data["dates"]);nm=len(data["months"])
        if nd<MIN_DATES or nm<MIN_MONTHS:continue
        native[country]+=1
        h=hashlib.sha256((country+"|"+scheme+"|"+sid).encode()).hexdigest()
        pos=positions.get((country,sid,scheme))
        if pos is None:
            quality[country+"_site_missing_or_ambiguous_original_source_coordinates"]+=1
            continue
        lat,lon=pos
        if not polygons[country].covers(Point(lon,lat)):
            quality[country+"_site_not_in_original_national_polygon"]+=1
            continue
        if not(-12<=lon<=26 and 44<=lat<=64):
            quality[country+"_site_outside_old_EU_HydroRIVERS_study_extent"]+=1
            continue
        quality[country+"_valid_geocoded"]+=1
        output.append({
          "countryCode":country,
          "original_country_scheme_site_SHA256":h,
          "original_WISE6_WGS84_lat":float(lat),
          "original_WISE6_WGS84_lon":float(lon),
          "original_2022_TP_independent_sample_dates":nd,
          "original_2022_TP_sample_months":nm,
          "original_2022_TP_record_metadata_rows":data["metadata_rows"],
          "original_source_TP_water_analyte_codes":sorted(data["codes"]),
          "original_source_TP_water_uom":sorted(data["units"]),
          "original_2degree_block":str(math.floor(lat/2))+"_"+str(math.floor(lon/2)),
          "was_in_original_187_site_butterfly_photo_region":h in old_sites,
          "original_2022_water_numeric_TP_concentrations_NOT_READ":True,
        })
    if dict(native)!=EXPECTED_SITE_NUM:
        raise ValueError("National original 2022 TP site count does not match frozen 242 DE/100 AT source metadata")
    if len(set(x["original_country_scheme_site_SHA256"] for x in output))!=len(output):
        raise ValueError("National 2022 site keys not one-to-one")
    overlap=sum(x["was_in_original_187_site_butterfly_photo_region"] for x in output)
    if overlap!=EXPECTED_OLD_SITE:
        raise ValueError("Pre-existing 187 WISE6 site source hashes not all in new national source list")
    result={
       "schema":"gesui_2022_EEA2026_original_WISE6_TP_national_DE_AT_public_station_metadata_source_only_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "original_full_WISE6_disaggregated_2026_zip_sha256":source["archive_sha256"],
       "original_full_source_bytes":source["archive_size_bytes"],
       "original_2026_WISE6_spatial_source_receipt":spatial,
       "original_2022_WISE6_sample_metadata_QC":originalQC,
       "original_2026_WISE6_spatial_metadata_QC":spatialQC,
       "source_original_187_photo_cell_site_metadata_SHA256":sha(args.source_187_stations),
       "source_original_national_date_eligible_2022_TP_site_total":sum(native.values()),
       "source_original_date_eligible_DE_AT_by_country":dict(native),
       "source_original_country_polygon_unambiguous_site_counts":dict(quality),
       "source_original_country_geocoded_valid_stations":len(output),
       "source_original_187_previous_butterfly_photo_cells_station_subset_verified":overlap,
       "source_original_NEW_nonphoto_site_count":len(output)-overlap,
       "nationwide_original_source_2022_station_registry":sorted(output,key=lambda x:(x["countryCode"],x["original_country_scheme_site_SHA256"])),
       "source_only_no_river_TP_numeric_values_parsed_or_selected":True,
       "source_only_no_iNaturalist_butterfly_or_photo_values_read":True,
       "national_environmental_TP_population_DISTINCT_from_42_source_positive_photo_cells":True,
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_EEA2026_ORIGINAL_2022_NATIONAL_DE_AT_RIVER_TP_SOURCE_METADATA_VERIFIED",
        json.dumps({"national_date_eligible":dict(native),"geocoded":dict(quality),
           "total":len(output),"original_photo_region_site_subset":overlap,
           "new_nonphoto_source_sites":len(output)-overlap,
           "no_TP_numeric_results_opened":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
