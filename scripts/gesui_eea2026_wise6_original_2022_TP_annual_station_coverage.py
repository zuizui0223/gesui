#!/usr/bin/env python3
"""Original 2026 EEA Waterbase WISE6 2022 TP annual monitoring site SUPPORT only.

Download verified published ZIP archives (aggregated 200MB ZIP / 2.7GB
expanded CSV, spatial 3.4MB ZIP). Stream semicolon-delimited original
annual CSV, read ONLY first 12 provenance/measurement-count fields, NO
numeric TP concentration resultMeanValue/resultMedianValue etc. Country
geometry uses SHA-verified original NaturalEarth50m boundaries, and strict
3-year 42 DE/AT Lycaenidae source tiles were fixed previously. No WWTP
effect and no butterfly source association fitted here.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,urllib.request,zipfile,time
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from shapely.geometry import shape,Point

URL_2026_AGG="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_AggregatedData"
URL_2026_SPATIAL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_S_WISE_SpatialObject_DerivedData"
EXPECTED_AGG_SHA256="8492c92702524ce23ac488cba62631f277b1e9bd4e5ea8ece3bea37569307b7c"
EXPECTED_SPATIAL_SHA256="afd64e1d5582492c09453f5a05aae0b1a9d9a38c8814778afdb5c07e43a1c974"
NATURAL_URL="https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_50m_admin_0_countries.geojson"
NATURAL_SHA256="3e458fc036ad0a66411f2c1e6cac49c5d7bfb81cb1123bc513b22511a2b7fdeb"
EXPECTED_3YEAR_SOURCE_CELLS={"DE":27,"AT":15}
MIN_TP_SAMPLES_PER_REPORTED_AGGREGATE=4
MIN_TP_SITES_PER_COUNTRY1DEG_CELL=2
MIN_RECEPTOR_CELLS_TOTAL=20
MIN_RECEPTOR_CELLS_DE=10
MIN_RECEPTOR_CELLS_AT=6
MIN_RECEPTOR_2DEG_BLOCKS=10
EXPECTED_EEA_WISE6_AGG_FILE="waterbase_t_wise6_aggregateddata.csv"
EXPECTED_EEA_WISE6_SPAT_FILE="waterbase_s_wise_spatialobject_deriveddata.csv"
ROWS_HARD_CAP=42_000_000

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for ch in iter(lambda:f.read(1048576),b""):h.update(ch)
    return h.hexdigest()

def load(url,dst,expected_sha,max_size):
    req=urllib.request.Request(url,headers={"User-Agent":"gesui-sourceonly-EEA2026-WISE6-measured-2022-TP-site-countryQA/1.0"})
    h=hashlib.sha256();count=0
    with urllib.request.urlopen(req,timeout=150) as r,dst.open("wb") as f:
        for ch in iter(lambda:r.read(1024*1024),b""):
            count+=len(ch)
            if count>max_size:raise ValueError("Original official zip unexpectedly exceeds size ceiling")
            h.update(ch);f.write(ch)
    if h.hexdigest()!=expected_sha:raise ValueError("Original published EEA 2026 source archive checksum drift: "+url)
    return {"public_original_url":url,"zip_sha256":h.hexdigest(),"size_bytes":count}

def make_key(country,id,scheme):
    return (country.strip('"'),id.strip('"'),scheme.strip('"'))

def get_geometries():
    req=urllib.request.Request(NATURAL_URL,headers={"User-Agent":"gesui-2026-WISE6-NaturalEarth-country-mask/0.1"})
    with urllib.request.urlopen(req,timeout=90) as r:raw=r.read(5000000)
    if hashlib.sha256(raw).hexdigest()!=NATURAL_SHA256:raise ValueError("Original Natural Earth country polygon checksum drift")
    obj=json.loads(raw)
    out={}
    for c,three in (("DE","DEU"),("AT","AUT")):
        features=[f for f in obj["features"] if
              (f.get("properties",{}).get("ISO_A2")==c or
               f.get("properties",{}).get("ADM0_A3")==three) and f.get("geometry")]
        if len(features)!=1:raise ValueError("Original exact country geometry not unique: "+c)
        out[c]=shape(features[0]["geometry"])
    return out

def spatial_source(zip_path):
    positions={}; qc=Counter();example=[]
    with zipfile.ZipFile(zip_path) as z:
        with io.TextIOWrapper(z.open(EXPECTED_EEA_WISE6_SPAT_FILE),encoding="utf-8-sig",
                              newline="",errors="replace") as f:
            rd=csv.DictReader(f,delimiter=";")
            cols=rd.fieldnames or []
            for k in ("countryCode","monitoringSiteIdentifier",
                      "monitoringSiteIdentifierScheme","lat","lon"):
                if k not in cols:raise ValueError("Missing original EEA2026 WISE spatial source field "+k)
            for q in rd:
                qc["spatial_source_records"]+=1
                c=str(q["countryCode"] or "").strip()
                if c not in ("DE","AT"):continue
                qc["spatial_country_"+c+"_rows"]+=1
                try:
                    lat=float(q["lat"]);lon=float(q["lon"])
                except (TypeError,ValueError):qc["spatial_invalid_latlon"]+=1;continue
                if not(-90<=lat<=90 and -180<=lon<=180) or not math.isfinite(lat+lon):
                    qc["spatial_invalid_latlon"]+=1;continue
                key=make_key(c,q["monitoringSiteIdentifier"],q["monitoringSiteIdentifierScheme"])
                if key in positions:
                    prev=positions[key]
                    if prev is None:
                        qc["additional_duplicate_for_already_ambiguous_site_key"]+=1
                    elif abs(prev[0]-lat)>1e-6 or abs(prev[1]-lon)>1e-6:
                        qc["different_positions_same_site_key"]+=1
                        # Drop site only when public spatial rows disagree.
                        positions[key]=None
                    else:
                        qc["duplicate_key_same_location"]+=1
                        # Identical site metadata rows are harmless, retain coordinate.
                    continue
                positions[key]=(lat,lon)
    return positions,dict(qc)

def strip(v):return str(v or "").strip().strip('"')

def station_metadata(zip_path):
    qc=Counter();label_counts=Counter();tp_entries={}; tp_uom=Counter();tp_period=Counter(); tp_codes=Counter()
    categories_2022=Counter();labels_2022_allcategories=Counter();matrices_2022=Counter()
    with zipfile.ZipFile(zip_path) as z:
        with io.TextIOWrapper(z.open(EXPECTED_EEA_WISE6_AGG_FILE),encoding="utf-8-sig",errors="replace",newline="") as f:
            header=f.readline().strip().split(";")
            need=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                  "parameterWaterBodyCategory","observedPropertyDeterminandCode",
                  "observedPropertyDeterminandLabel","procedureAnalysedMatrix","resultUom",
                  "phenomenonTimeReferenceYear","parameterSamplingPeriod","resultNumberOfSamples"]
            if any(v not in header for v in need):raise ValueError("Original 2026 WISE6 annual aggregated source schema drift")
            if len(header)<26:raise ValueError("Wrong original WISE6 year-TP source record")
            ix={k:header.index(k) for k in need}
            max_ix=max(ix.values())
            for line in f:
                qc["official_annual_data_records_checked"]+=1
                if qc["official_annual_data_records_checked"]>ROWS_HARD_CAP:raise ValueError("Bound on official Waterbase records exceeded")
                # Parse only original first 12 fields. Numeric resultMeanValue,
                # resultMedianValue and other concentration columns are not split,
                # accessed, or loaded into any result.
                vals=line.split(";",max_ix+1)
                if len(vals)<=max_ix:qc["invalid_annual_short_rows"]+=1;continue
                c=strip(vals[ix["countryCode"]])
                if c not in ("DE","AT"):continue
                qc["annual_country_"+c+"_all_year_rows"]+=1
                year=strip(vals[ix["phenomenonTimeReferenceYear"]])
                if year!="2022":continue
                qc["annual_country_"+c+"_2022_rows"]+=1
                cat=strip(vals[ix["parameterWaterBodyCategory"]])
                categories_2022[(c,cat)]+=1
                if cat!="RW":
                    # Source metadata only, to distinguish original WISE country
                    # category coding / missingness from true absence of data.
                    labels_2022_allcategories[(c,cat,strip(vals[ix["observedPropertyDeterminandLabel"]]))]+=1
                    continue
                matrices_2022[(c,strip(vals[ix["procedureAnalysedMatrix"]]))]+=1
                qc["annual_country_"+c+"_2022_RW_rows"]+=1
                lab=strip(vals[ix["observedPropertyDeterminandLabel"]])
                label_counts[(c,lab)]+=1
                if lab.casefold()!="total phosphorus":continue
                qc["annual_country_"+c+"_2022_RW_total_phosphorus_rows"]+=1
                tp_codes[(c,strip(vals[ix["observedPropertyDeterminandCode"]]))]+=1
                tp_uom[(c,strip(vals[ix["resultUom"]]))]+=1
                tp_period[(c,strip(vals[ix["parameterSamplingPeriod"]]))]+=1
                matrix=strip(vals[ix["procedureAnalysedMatrix"]])
                try:
                    nn=int(float(strip(vals[ix["resultNumberOfSamples"]])))
                except (ValueError,TypeError):
                    qc["TP_2022_source_missing_or_invalid_sample_count"]+=1
                    nn=0
                if nn<MIN_TP_SAMPLES_PER_REPORTED_AGGREGATE:
                    qc["TP_2022_fewer_than4_samples"]+=1
                    continue
                key=make_key(c,vals[ix["monitoringSiteIdentifier"]],vals[ix["monitoringSiteIdentifierScheme"]])
                if not all(key):continue
                qc["TP_2022_annual_records_with_ge4_samples"]+=1
                prev=tp_entries.get(key)
                # Retain site visit metadata only, not concentration or nutrient
                # values from the original file.
                ent={"countryCode":c,"original_site_key_sha256":hashlib.sha256("|".join(key).encode()).hexdigest(),
                     "original_TP_analyte_code":strip(vals[ix["observedPropertyDeterminandCode"]]),
                     "original_matrix_code":matrix,
                     "original_result_unit":strip(vals[ix["resultUom"]]),
                     "original_TP_2022_aggregate_sampling_period":strip(vals[ix["parameterSamplingPeriod"]]),
                     "original_TP_2022_one_aggregate_n_sample_support":nn}
                if prev is None or nn>prev["original_TP_2022_one_aggregate_n_sample_support"]:
                    tp_entries[key]=ent
            if qc["official_annual_data_records_checked"]<10:raise ValueError("No Waterbase annual records")
    return tp_entries,{
      "source_row_counts":dict(qc),
      "2022_original_WISE6_country_waterbody_category_code_counts":{str(k[0])+"|"+str(k[1]):v for k,v in categories_2022.items()},
      "2022_original_WISE6_country_RW_matrix_counts":{str(k[0])+"|"+str(k[1]):v for k,v in matrices_2022.items()},
      "2022_original_WISE6_top_analyte_labels_OUTSIDE_original_RW_category":[
          {"countryCode":k[0],"original_category":k[1],"label":k[2],"rows":v}
          for k,v in labels_2022_allcategories.most_common(60)],
      "2022_RW_observed_analyte_labels_top35":[{"countryCode":k[0],"label":k[1],"rows":v} for k,v in label_counts.most_common(35)],
      "strict_total_phosphorus_analyte_code_counts":{k[0]+"|"+k[1]:v for k,v in tp_codes.items()},
      "strict_total_phosphorus_original_result_UOM_counts":{k[0]+"|"+k[1]:v for k,v in tp_uom.items()},
      "strict_total_phosphorus_sample_period_counts":{k[0]+"|"+k[1]:v for k,v in tp_period.items()},
      "strict_total_phosphorus_2022_siteIDs_with_ge4_original_annual_samples":len(tp_entries),
      "no_WQ_mean_or_median_TP_numeric_outcome_columns_parsed":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--2022-original-WorldCover-source",type=Path,required=True,dest="wc")
    p.add_argument("--2023-exact-photo-cell-source",type=Path,required=True,dest="photo")
    p.add_argument("--tmpdir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    world=json.loads(a.wc.read_text())
    photo=json.loads(a.photo.read_text())
    if world["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED" or not photo["2023_outcome_year_source_gate_pass"]:
        raise ValueError("Original 2022/23 Lycaenidae source cell panel not locked")
    if photo["original_2022_WorldCover_land_photo_source_sha256"]!=sha(a.wc):
        raise ValueError("Original 42-cell 2023 photo source not linked to original 2022 WorldCover")
    prev={(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))):x for x in photo["source_2023_all_original43_butterfly_photo_count_and_receipts"]}
    original=[x for x in world["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
              if prev[(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"])))]["2023_photo_support_ge20"]]
    if len(original)!=42 or Counter(x["countryCode"] for x in original)!=Counter({"DE":27,"AT":15}):
        raise ValueError("Original strict DE/AT 2021-23 butterfly-photo geography changed")
    original_keys={(x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"]))):x for x in original}
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    aggzip=a.tmpdir/"EEA2026_WISE6_aggregated_original.zip"
    spatialzip=a.tmpdir/"EEA2026_WISE6_spatial_original.zip"
    aggReceipt=load(URL_2026_AGG,aggzip,EXPECTED_AGG_SHA256,250000000)
    spatialReceipt=load(URL_2026_SPATIAL,spatialzip,EXPECTED_SPATIAL_SHA256,6000000)
    positions,sourceSpace=spatial_source(spatialzip)
    sites,quality=station_metadata(aggzip)
    countries=get_geometries()
    geoqc=Counter();tiles=defaultdict(list)
    for key,item in sites.items():
        c=key[0];station=positions.get(key)
        if station is None:
            geoqc["TP_source_station_missing_original_WISE_spatial_metadata_or_position_ambiguous"]+=1
            continue
        lat,lon=station
        if not countries[c].covers(Point(lon,lat)):
            geoqc["official_site_countryCode_not_within_NaturalEarth50m_country_polygon"]+=1
            continue
        geoqc["2022_RiverTP_original_national_site_geo_valid_"+c]+=1
        tile=(c,math.floor(lat),math.floor(lon))
        if tile not in original_keys:continue
        geoqc["2022_RiverTP_original_3year_photo_cell_station_source_linked_"+c]+=1
        # Metadata only; numeric P concentration not opened.
        tiles[tile].append({
           "original_site_SHA256":item["original_site_key_sha256"],
           "original_TP_2022_one_aggregate_n_sample_support":item["original_TP_2022_one_aggregate_n_sample_support"],
           "original_TP_analyte_code":item["original_TP_analyte_code"],
           "original_TP_unit":item["original_result_unit"],
           "original_TP_sampling_period":item["original_TP_2022_aggregate_sampling_period"],
        })
    cells=[]
    for key,x in sorted(original_keys.items()):
        these=tiles[key]
        cells.append({
           "countryCode":key[0],"lat_floor":key[1],"lon_floor":key[2],
           "country_2degree_block":x["country_2deg_block"],
           "official_2026_WISE6_original_TP_annual_year2022_distinct_station_n_ge4samples":len(these),
           "original_2026_waterquality_TP_2022_site_metadata":these,
           "original_station_sample_support_ge2":len(these)>=MIN_TP_SITES_PER_COUNTRY1DEG_CELL,
           "NO_TP_NUMERIC_CONCENTRATION_RESULTS_OPENED":True,
        })
    supported=[x for x in cells if x["original_station_sample_support_ge2"]]
    cby=Counter(x["countryCode"] for x in supported)
    blocks=set(x["countryCode"]+"|"+str(x["country_2degree_block"]) for x in supported)
    complete=(len(supported)>=MIN_RECEPTOR_CELLS_TOTAL and cby["DE"]>=MIN_RECEPTOR_CELLS_DE and
              cby["AT"]>=MIN_RECEPTOR_CELLS_AT and len(blocks)>=MIN_RECEPTOR_2DEG_BLOCKS)
    out={
      "schema":"gesui_EEA2026_original_WISE6_2022_countrymasked_RW_TP_aggregated_station_support_in_exact_42_Lycaenidae_cells_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "source_official_EEA2026_aggregated_zip":aggReceipt,
      "source_official_EEA2026_spatial_zip":spatialReceipt,
      "source_official_EEA2026_WISE_spatial_site_original_QC":sourceSpace,
      "source_WQ2022_metadata_only_original":quality,
      "strict_2021_2023_country_1deg_photo_sample_total":42,
      "original_2022_WorldCover43_source_SHA256":sha(a.wc),
      "original_2023_photo_source_SHA256":sha(a.photo),
      "country_masked_geographic_quality":dict(geoqc),
      "strict_2021_2022_2023_original42_photo_cell_station_metadata":cells,
      "source_original_year_2022_aggregated_total_P_label_exact":"Total phosphorus",
      "source_original_category_RW_only":True,
      "original_aggregated_reported_resultNumberOfSamples_ge4":True,
      "source_2022_TP_cells_with_ge2_annual_aggregated_station_sources":len(supported),
      "source_2022_TP_supported_cells_by_country":dict(cby),
      "source_2022_TP_distinct_country2deg_source_blocks":len(blocks),
      "new_riverTP_vs_lycaenidae_three_layer_geographic_source_gate":complete,
      "no_original_TP_concentration_or_chemical_effect_values_loaded":True,
      "no_numeric_river_WQ_results_used_to_select_country_cells":True,
      "no_photo_or_observer_personal_data_read":True,
      "no_official_EEA2026_dated_disaggregated_data_unzipped_or_needed_here":True,
      "source_zip_200MB_original_not_committed_or_uploaded":True,
      "noncovered_WISE6_disaggregated_data_may_fill_aggregated_missingness":True,
      "WWTP_reported_P_tonnes_not_at_station_P_mgL":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_ORIGINAL_EEA2026_WISE6_2022_TOTALPHOSPHORUS_COUNTRY_SITE_MATCH_SOURCE_ONLY",
      json.dumps({"national_TP_site_QC":dict(geoqc),
          "source_2022_RW_TP_annual_stations":quality["strict_total_phosphorus_2022_siteIDs_with_ge4_original_annual_samples"],
          "supported_original42_cells":len(supported),"bycountry":dict(cby),"2deg_blocks":len(blocks),
          "bridge_support_gate":complete,"numeric_TP_outcomes_not_opened":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
