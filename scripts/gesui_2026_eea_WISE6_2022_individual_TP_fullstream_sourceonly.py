#!/usr/bin/env python3
"""Stream original EEA2026 WISE6 2022 DE/AT river-TP individual SAMPLE METADATA only.

Official original ZIP is ~2.024GB, containing one ~36GB uncompressed CSV.
The official SQL view returned HTTP403, but the publisher ZIP is open.
Download the original ZIP, SHA256 it, STREAM decompress its sole CSV, reject
non-DE/AT rows early. Inspect ONLY original 11 first metadata columns:
(country, station id/scheme, RW, analyte label/code, water matrix/unit,
2022 original sample date/year, and sample id). Numeric resultObservedValue
and all following fields MUST NEVER be parsed or used to select sites.

Join 2022 sites with exact original WISE6 spatial identifiers and
NaturalEarth50m DE/AT polygon; crosswalk the original 42 1degree country
source/photo cells selected without viewing this chemistry. Report complete
source-only support. DO NOT estimate source-WQ-photo association here.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time,urllib.request,zipfile
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path
from shapely.geometry import Point
import gesui_eea2026_wise6_original_2022_TP_annual_station_coverage as old

URL=old.URL_2026_SPATIAL.replace("Waterbase_S_WISE_SpatialObject_DerivedData",
                              "Waterbase_T_WISE6_DisaggregatedData")
SPATIAL_URL=old.URL_2026_SPATIAL
ZIP_MEMBER="waterbase_t_wise6_disaggregateddata.csv"
EXPECTED_ARCHIVE_COMPRESSED_BYTES=2024284753
EXPECTED_MEMBER_UNCOMPRESSED_BYTES=36043700357
MAX_COMPRESSED_BYTES=2300000000
STRICT_YEAR="2022"
COUNTRIES=("DE","AT")
MIN_DAYS=4
MIN_MONTHS=3
MIN_CELL_STATIONS=2
SUPPORT_MIN_TOTAL=20
SUPPORT_MIN_DE=10
SUPPORT_MIN_AT=6
SUPPORT_MIN_2DEG=10
MAX_SOURCE_LINES=440000000

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as r:
        for chunk in iter(lambda:r.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def original_archive(dst):
    headers={
      "User-Agent":"gesui-official-WISE6-2026-original-riverTP-2022-DEAT-metadata-only-stream/1.0",
      "Accept":"application/zip,*/*"}
    request=urllib.request.Request(URL,headers=headers)
    h=hashlib.sha256();size=0
    with urllib.request.urlopen(request,timeout=180) as r,dst.open("wb") as f:
        if r.status!=200:raise ValueError("Original EEA WISE6 2026 individual ZIP inaccessible")
        declared=r.headers.get("Content-Length")
        if declared and int(declared)!=EXPECTED_ARCHIVE_COMPRESSED_BYTES:
            raise ValueError("Original WISE6 2026 source ZIP size changed, do not silently use another release")
        for chunk in iter(lambda:r.read(1024*1024),b""):
            size+=len(chunk)
            if size>MAX_COMPRESSED_BYTES:raise ValueError("Unexpected original zip size")
            h.update(chunk);f.write(chunk)
    if size!=EXPECTED_ARCHIVE_COMPRESSED_BYTES:raise ValueError("Original public WISE6 compressed source length mismatch")
    if not zipfile.is_zipfile(dst):raise ValueError("Original EEA source downloaded is not ZIP")
    return {"original_source_url":URL,"archive_size_bytes":size,"archive_sha256":h.hexdigest()}

def strip(b):
    # EEA original CSV is semicolon delimited; only metadata tokens are decoded.
    return b.strip().strip(b'"').decode("utf-8",errors="replace")

def summarize_original_DE_AT_TP_sample_sites(zippath):
    qc=Counter()
    label_counts=Counter()
    total_tp_codes=Counter()
    total_tp_units=Counter()
    total_tp_matrix=Counter()
    original_site_support={}
    with zipfile.ZipFile(zippath) as z:
        names=z.namelist()
        if names!=[ZIP_MEMBER]:raise ValueError("Original EEA2026 individual ZIP members changed: "+str(names))
        info=z.getinfo(ZIP_MEMBER)
        if info.file_size!=EXPECTED_MEMBER_UNCOMPRESSED_BYTES:
            raise ValueError("Original published 2026 WISE6 uncompressed source changed")
        with z.open(ZIP_MEMBER) as f:
            h=f.readline()
            hdr=h.decode("utf-8-sig",errors="strict").rstrip("\r\n").split(";")
            expected=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                "parameterWaterBodyCategory","observedPropertyDeterminandCode",
                "observedPropertyDeterminandLabel","procedureAnalysedMatrix","resultUom",
                "phenomenonTimeSamplingDate","phenomenonTimeReferenceYear","sampleIdentifier",
                "resultObservedValue"]
            if hdr[:12]!=expected:
                raise ValueError("Published EEA 2026 WISE6 individual CSV meta schema changed")
            # For *every* original record, only country prefix is needed. We
            # never parse resultObservedValue, even when sample metadata matches.
            for line in f:
                qc["original_2026_WISE6_individual_total_source_rows_streamed"]+=1
                if qc["original_2026_WISE6_individual_total_source_rows_streamed"]>MAX_SOURCE_LINES:
                    raise ValueError("Unexpected 2026 source row count / source format drift")
                if not(line.startswith(b"DE;") or line.startswith(b"AT;") or
                       line.startswith(b'"DE";') or line.startswith(b'"AT";')):
                    continue
                # Only split FIRST ELEVEN METADATA FIELDS, leaves everything
                # from resultObservedValue onward as an unused final tail.
                tokens=line.split(b";",11)
                if len(tokens)<12:
                    qc["original_DE_AT_short_records"]+=1
                    continue
                c=strip(tokens[0])
                if c not in COUNTRIES:continue
                qc[f"original_{c}_all_year_individual_records"]+=1
                date=strip(tokens[8]);year=strip(tokens[9])
                if year!=STRICT_YEAR or not date.startswith("2022-"):
                    continue
                qc[f"original_{c}_2022_sample_rows"]+=1
                kind=strip(tokens[3])
                if kind!="RW":
                    qc[f"original_{c}_2022_non_RW_{kind[:12]}"]+=1
                    continue
                qc[f"original_{c}_2022_RW_sample_rows"]+=1
                analyte=strip(tokens[5])
                label_counts[(c,analyte)]+=1
                if analyte.casefold()!="total phosphorus":
                    continue
                qc[f"original_{c}_2022_RW_total_phosphorus_sample_rows"]+=1
                matrix=strip(tokens[6]);uom=strip(tokens[7]);code=strip(tokens[4])
                total_tp_matrix[(c,matrix)]+=1
                total_tp_units[(c,uom)]+=1
                total_tp_codes[(c,code)]+=1
                # exact chemical analyzed matrix must be water, not biota or sediment.
                if matrix not in ("W",):
                    qc[f"original_{c}_2022_RW_total_P_non_water_matrix_rows"]+=1
                    continue
                sid=strip(tokens[1]);scheme=strip(tokens[2])
                if not sid or not scheme:
                    qc[f"original_{c}_2022_RW_total_P_blank_site_scheme_rows"]+=1
                    continue
                sampleid=strip(tokens[10])
                key=(c,sid,scheme)
                if key not in original_site_support:
                    original_site_support[key]={"dates":set(),"months":set(),"metadata_rows":0,
                                                 "unique_sample_ids":set(),"units":set(),"codes":set()}
                v=original_site_support[key]
                v["metadata_rows"]+=1
                if len(date)>=10 and date[4]=="-" and date[7]=="-":
                    v["dates"].add(date[:10]);v["months"].add(date[:7])
                if sampleid:v["unique_sample_ids"].add(sampleid)
                v["units"].add(uom);v["codes"].add(code)
                if qc["original_2026_WISE6_individual_total_source_rows_streamed"]%12000000==0:
                    print("GESUI_EEA2026_ORIGINAL_36GB_METADATA_STREAM_12M_ROWS",
                          json.dumps({"rows_total":qc["original_2026_WISE6_individual_total_source_rows_streamed"],
                            "DE_2022_TP_water":qc["original_DE_2022_RW_total_phosphorus_sample_rows"],
                            "AT_2022_TP_water":qc["original_AT_2022_RW_total_phosphorus_sample_rows"],
                            "water_sites":len(original_site_support)}),flush=True)
    output={
      "all_source_rows_and_country_year_category_counts":dict(qc),
      "country_2022_RW_determinand_label_counts_top":[{"country":c,"analyte":name,"n":n}
                for (c,name),n in label_counts.most_common(100)],
      "country_2022_RW_total_phosphorus_analyte_code_counts":[{"country":c,"code":name,"n":n}
                for (c,name),n in total_tp_codes.most_common()],
      "country_2022_RW_total_phosphorus_original_water_sample_unit_counts":[{"country":c,"unit":name,"n":n}
                for (c,name),n in total_tp_units.most_common()],
      "country_2022_RW_total_phosphorus_original_matrix_counts":[{"country":c,"matrix":name,"n":n}
                for (c,name),n in total_tp_matrix.most_common()],
      "original_2022_RW_TP_site_count_any_water_samples":len(original_site_support),
      "original_2022_RW_TP_sites_ge4_distinct_sample_dates_and_ge3_months":
         sum(len(v["dates"])>=MIN_DAYS and len(v["months"])>=MIN_MONTHS for v in original_site_support.values()),
      "original_source_2026_WISE6_member_uncompressed_bytes":EXPECTED_MEMBER_UNCOMPRESSED_BYTES,
      "NO_NUMERIC_TP_CONCENTRATION_FIELDS_PARSED_OR_READ":True,
    }
    return original_site_support,output

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--2022-original-worldcover",type=Path,required=True,dest="wc")
    a.add_argument("--2023-photo-source",type=Path,required=True,dest="photo")
    a.add_argument("--tmpdir",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    p=a.parse_args()
    world=json.loads(p.wc.read_text());photos=json.loads(p.photo.read_text())
    if world["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED" or not photos["2023_outcome_year_source_gate_pass"]:
        raise ValueError("2021–2023 42 country-cell source is not fixed")
    if photos["original_2022_WorldCover_land_photo_source_sha256"]!=sha(p.wc):
        raise ValueError("Original photo country×cell geometry not locked to WorldCover source")
    available={(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))):
        x for x in photos["source_2023_all_original43_butterfly_photo_count_and_receipts"]
        if x["2023_photo_support_ge20"]}
    original=[x for x in world["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
            if (x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))) in available]
    if len(original)!=42 or Counter(q["countryCode"] for q in original)!=Counter({"DE":27,"AT":15}):
        raise ValueError("Original strict balanced three-year sample unexpectedly changed")
    original_cells={(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))):x for x in original}
    p.tmpdir.mkdir(parents=True,exist_ok=True)
    zipfilelocal=p.tmpdir/"original_EEA2026_WISE6_2022_DA.zip"
    spatialzip=p.tmpdir/"original_EEA2026_WISE6_spatial.zip"
    spatial_receipt=old.load(SPATIAL_URL,spatialzip,old.EXPECTED_SPATIAL_SHA256,6000000)
    positions,spaceqc=old.spatial_source(spatialzip)
    masks=old.get_geometries()
    sourcezip=original_archive(zipfilelocal)
    print("GESUI_EEA2026_WISE6_2022_RIVERTP_SOURCE_FULL_ZIP_SHA_LOCKED",
          json.dumps({"download_sha256":sourcezip["archive_sha256"],"archive_bytes":sourcezip["archive_size_bytes"]}),flush=True)
    stations,sourceqc=summarize_original_DE_AT_TP_sample_sites(zipfilelocal)
    # Exact monitoringSiteIdentifier AND Scheme AND country join, no numeric
    # river values or wastewater-source outcomes used to select.
    sourcecount=Counter();tile_stations=defaultdict(list);all_qc_sites=0
    for (country,sid,scheme),data in stations.items():
        if len(data["dates"])<MIN_DAYS or len(data["months"])<MIN_MONTHS:
            sourcecount[f"{country}_TP_site_less_than4dates_or3months"]+=1;continue
        sourcecount[f"{country}_2022_TP_date_supported_all_national_sites"]+=1
        pos=positions.get((country,sid,scheme))
        if pos is None:
            sourcecount[f"{country}_TP_date_supported_site_unlinked_to_original_WISE_coordinates"]+=1;continue
        lat,lon=pos
        if not masks[country].covers(Point(lon,lat)):
            sourcecount[f"{country}_TP_2022_site_outside_NaturalEarth50m_country"]+=1;continue
        sourcecount[f"{country}_2022_TP_sites_with_actual_original_country_geocoordinates"]+=1
        tile=(country,math.floor(lat),math.floor(lon))
        if tile not in original_cells:continue
        sourcecount[f"{country}_2022_TP_sites_in_original_42_Lycaenidae_photo_cells"]+=1
        all_qc_sites+=1
        tile_stations[tile].append({
            "monitoring_site_ID_with_scheme_sha256":hashlib.sha256((country+"|"+scheme+"|"+sid).encode()).hexdigest(),
            "original_2022_TP_distinct_sample_dates":len(data["dates"]),
            "original_2022_TP_distinct_sample_months":len(data["months"]),
            "original_2022_TP_sample_metadata_rows":data["metadata_rows"],
            "original_2022_TP_unique_sample_ID_count":len(data["unique_sample_ids"]),
            "original_source_TP_analyte_codes":sorted(data["codes"]),
            "original_source_TP_units":sorted(data["units"]),
            "no_original_numeric_TP_concentrations_accessed":True,
        })
    sourcegrid=[]
    for key,cell in sorted(original_cells.items()):
        items=tile_stations[key]
        sourcegrid.append({
          "country":key[0],"1deg_lat_floor":key[1],"1deg_lon_floor":key[2],
          "original_2deg_spatial_validation_block":cell["country_2deg_block"],
          "original_EEA2026_2022_RW_TP_individual_sample_stations_with_ge4dates_3months":len(items),
          "original_non_numeric_WISE6_2022_station_sampling_support":items,
          "original_ge2_TP_monitoring_sites_cell":len(items)>=MIN_CELL_STATIONS,
          "NUMERIC_TP_WATER_CONCENTRATION_NOT_READ":True,
        })
    good=[q for q in sourcegrid if q["original_ge2_TP_monitoring_sites_cell"]]
    by=Counter(q["country"] for q in good)
    blocks=set(q["country"]+"|"+str(q["original_2deg_spatial_validation_block"]) for q in good)
    coveragepass=(len(good)>=SUPPORT_MIN_TOTAL and by["DE"]>=SUPPORT_MIN_DE and
                  by["AT"]>=SUPPORT_MIN_AT and len(blocks)>=SUPPORT_MIN_2DEG)
    result={
      "schema":"gesui_official_EEA2026_WISE6_disaggregated_original_2022_RW_TP_sample_date_coverage_DE_AT_Lycaenidae_sourceonly_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "original_EEA2026_WISE6_full_2GB_individual_zip_SHA256":sourcezip["archive_sha256"],
      "original_EEA2026_WISE6_original_disagg_download_receipt":sourcezip,
      "original_EEA2026_WISE6_spatial_source_receipt":spatial_receipt,
      "original_WISE6_spatial_country_source_QC":spaceqc,
      "original_2022_WorldCover43_SHA256":sha(p.wc),
      "original_2023_frozen_photo_source_SHA256":sha(p.photo),
      "original_2022_TP_sample_country_analyte_and_station_metadata":sourceqc,
      "original_strict_2021_2023_source_photo_cells":len(original_cells),
      "original_2022_TP_country_coordinate_and_cell_join_QC":dict(sourcecount),
      "original_2022_TP_date_supported_stations_in_42_photo_cells":all_qc_sites,
      "original_42_country_1degree_photo_cells_and_WISE6_TP_stations":sourcegrid,
      "TP_source_cells_with_ge2_stations":len(good),
      "TP_source_cells_with_ge2_by_country":dict(by),
      "TP_source_cells_2deg_original_spatial_blocks":len(blocks),
      "source_preoutcome_geography_repeated_site_gate_20cells_DE10_AT6_blocks10":coveragepass,
      "no_original_numeric_river_TP_concentrations_or_LOQ_values_parsed":True,
      "no_butterfly_photo_API_queries_run_in_this_stage":True,
      "no_adult_puddling_or_wastewater_causal_claim":True,
      "source_metadata_requirements":"Original 2022, RW, exact Total phosphorus label, procedureAnalyzedMatrix W, >=4 distinct days >=3 months, exact country/source identifier scheme and original coordinate",
      "note_bulk_2GB_zip_not_committed_or_shared":True,
      "one_country_1degree_photo_reporting_tile_is_not_a_single_river_reach":True,
    }
    p.out.parent.mkdir(parents=True,exist_ok=True)
    p.out.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_EEA2026_OFFICIAL_INDIVIDUAL_RIVER_TP_2022_COUNTRY_SITE_SUPPORT_GATE_COMPLETED",
          json.dumps({"original_raw_source_rows":sourceqc["all_source_rows_and_country_year_category_counts"].get("original_2026_WISE6_individual_total_source_rows_streamed"),
          "2022_DE_AT_TP_source_rows":{k:v for k,v in sourceqc["all_source_rows_and_country_year_category_counts"].items() if "2022_RW_total_phosphorus" in k},
          "national_date_supported_sites":{k:v for k,v in sourcecount.items() if "date_supported_all_national" in k},
          "photo_cell_sites":all_qc_sites,"source_2site_tiles":len(good),"country":dict(by),
          "spatial_blocks":len(blocks),"gate":coveragepass,
          "no_numeric_TP_outcomes_read":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
