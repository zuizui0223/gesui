#!/usr/bin/env python3
"""Read numeric EEA2026 2022 river TOTAL phosphorus only at 187 pre-frozen sites.

This is the first reading of numeric 2022 original WISE6 TP concentrations
after the measured station/date support gate PASSED and numeric handling was
committed to the protocol. Output site-year median daily-measurements mg{P}/L,
then country×1° cell median of independent site medians, *not butterfly effects*.
"""
from __future__ import annotations
import argparse,hashlib,json,math,statistics,urllib.request,zipfile
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
ORIGINAL_ZIP_BYTES=2024284753
ORIGINAL_CSV_UNCOMPRESSED_BYTES=36043700357
ORIGINAL_SOURCE_ROWS=102569057
ORIGINAL_TP_ROWS_COUNTRY={"DE":3690,"AT":1281}
ORIGINAL_STATIONS=187
ZIP_MEMBER="waterbase_t_wise6_disaggregateddata.csv"
TP_CODE="CAS_7723-14-0"
TP_UNIT="mg{P}/L"
MIN_VALID_DATES=4
MIN_VALID_MONTHS=3
MAX_BELOWLOQ_FRAC=0.5
MIN_SITE_PER_CELL=2
MIN_CELLS=20
MIN_DE=10
MIN_AT=6
MIN_BLOCKS=10

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for x in iter(lambda:f.read(1048576),b""):h.update(x)
    return h.hexdigest()

def original_zip(p):
    req=urllib.request.Request(URL,headers={
     "User-Agent":"gesui-EEA2026-WISE6-TP2022-numeric-after-source-gate-LOQ/1.0",
     "Accept":"application/zip,*/*"})
    h=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=180) as r,p.open("wb") as out:
        if r.status!=200:raise ValueError("EEA official original water 2026 source not 200")
        for chunk in iter(lambda:r.read(1048576),b""):
            size+=len(chunk)
            if size>ORIGINAL_ZIP_BYTES:raise ValueError("EEA source archive size exceeded frozen original")
            h.update(chunk);out.write(chunk)
    if size!=ORIGINAL_ZIP_BYTES or h.hexdigest()!=SHA:
        raise ValueError("EEA original 2026 numeric source ZIP changed, refuse historical comparison")
    return {"original_url":URL,"original_archive_bytes":size,"original_archive_sha256":h.hexdigest()}

def s(x):return x.strip().strip(b'"').decode("utf8",errors="replace")

def use_original(sourcejson):
    if not sourcejson["source_preoutcome_geography_repeated_site_gate_20cells_DE10_AT6_blocks10"]:
        raise ValueError("Numeric TP must not open before 2022 source monitoring-site 3-layer support gate PASSES")
    expected=sourcejson["original_EEA2026_WISE6_full_2GB_individual_zip_SHA256"]
    if expected!=SHA:raise ValueError("Original water numerical source digest differs from source-only station source")
    q=sourcejson["original_42_country_1degree_photo_cells_and_WISE6_TP_stations"]
    keys={}
    if len(q)!=42:raise ValueError("Source-only original photo cells changed")
    for cell in q:
        c=cell["country"];lat=cell["1deg_lat_floor"];lon=cell["1deg_lon_floor"]
        for item in cell["original_non_numeric_WISE6_2022_station_sampling_support"]:
            idhash=item["monitoring_site_ID_with_scheme_sha256"]
            if idhash in keys:raise ValueError("Same original WISE site assigned to >1 country×1deg region")
            keys[idhash]={"country":c,"lat_floor":lat,"lon_floor":lon,
                     "source_metadata_rows":item["original_2022_TP_sample_metadata_rows"],
                     "source_metadata_valid_dates":item["original_2022_TP_distinct_sample_dates"],
                     "source_metadata_valid_months":item["original_2022_TP_distinct_sample_months"],
                     "source_TP_units":item["original_source_TP_units"],
                     "source_original_codes":item["original_source_TP_analyte_codes"]}
    if len(keys)!=ORIGINAL_STATIONS:raise ValueError("Frozen 2022 sample source station 187 count changed")
    return keys,q

def main():
    arg=argparse.ArgumentParser()
    arg.add_argument("--original-2022-tp-site-metadata",type=Path,required=True)
    arg.add_argument("--tmpdir",type=Path,required=True)
    arg.add_argument("--out",type=Path,required=True)
    a=arg.parse_args()
    src=json.loads(a.original_2022_tp_site_metadata.read_text())
    keys,cells=use_original(src)
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    file=a.tmpdir/"EEA2026_original_2022_TP_Numeric.zip"
    receipt=original_zip(file)
    qc=Counter();data={}
    for k in keys:
        data[k]={"original_sample_rows":0,"below_LOQ_rows":0,"unknown_flag_rows":0,
                 "original_missing_numeric_rows":0,"original_nonfinite_or_negative_rows":0,
                 "original_invalid_bool_flag_rows":0,"nonneg_0_concentration_rows":0,
                 "valid_quantified_by_day":defaultdict(list)}
    with zipfile.ZipFile(file) as z:
        if z.namelist()!=[ZIP_MEMBER]:raise ValueError("Original EEA2026 WISE6 ZIP member change")
        info=z.getinfo(ZIP_MEMBER)
        if info.file_size!=ORIGINAL_CSV_UNCOMPRESSED_BYTES:
            raise ValueError("Original WISE6 36GB numeric data changed")
        with z.open(info) as f:
            hdr=f.readline().decode("utf-8-sig").strip().split(";")
            if hdr[0]!="countryCode" or hdr[4]!="observedPropertyDeterminandCode" or hdr[5]!="observedPropertyDeterminandLabel" or hdr[11]!="resultObservedValue" or hdr[12]!="resultQualityObservedValueBelowLOQ":
                raise ValueError("Original source numeric TP result columns unexpectedly changed")
            for line in f:
                qc["original_2026_WISE6_full_source_rows_checked"]+=1
                if not(line.startswith((b"DE;",b"AT;",b'"DE";',b'"AT";'))):continue
                # At FIRST stage only the original sample station/data metadata
                # parsed. Numeric results parsed only if the original source-
                # frozen country/scheme/site SHA is among exactly 187 sites.
                v=line.split(b";",14)
                if len(v)<14:continue
                c=s(v[0])
                if c not in ("DE","AT"):continue
                if s(v[9])!="2022" or not s(v[8]).startswith("2022-") or s(v[3])!="RW" or s(v[5])!="Total phosphorus" or s(v[6])!="W":continue
                if s(v[4])!=TP_CODE or s(v[7])!=TP_UNIT:
                    raise ValueError("Original 2022 Total-P sample analyte code/unit changed")
                qc["original_2022_totalP_water_matrix_national_matched_sample_rows_"+c]+=1
                stationhash=hashlib.sha256((c+"|"+s(v[2])+"|"+s(v[1])).encode()).hexdigest()
                if stationhash not in keys:continue
                qc["source_only_frozen_187_station_observation_rows"]+=1
                date=s(v[8])[:10]
                item=data[stationhash];item["original_sample_rows"]+=1
                flag=s(v[12]).casefold()
                if flag in ("true","1"):
                    item["below_LOQ_rows"]+=1
                    continue
                if flag not in ("false","0",""):
                    item["original_invalid_bool_flag_rows"]+=1
                    continue
                if flag=="":
                    item["unknown_flag_rows"]+=1
                numstr=s(v[11])
                if not numstr:
                    item["original_missing_numeric_rows"]+=1
                    continue
                try:
                    val=float(numstr)
                except (ValueError,TypeError):
                    item["original_nonfinite_or_negative_rows"]+=1
                    continue
                if not math.isfinite(val) or val<0:
                    item["original_nonfinite_or_negative_rows"]+=1
                    continue
                if val==0:item["nonneg_0_concentration_rows"]+=1
                item["valid_quantified_by_day"][date].append(val)
    if qc["original_2026_WISE6_full_source_rows_checked"]!=ORIGINAL_SOURCE_ROWS:
        raise ValueError("Original 2026 EEA total source row count mismatch")
    for c in ("DE","AT"):
        if qc[f"original_2022_totalP_water_matrix_national_matched_sample_rows_{c}"]!=ORIGINAL_TP_ROWS_COUNTRY[c]:
            raise ValueError("Original country/year TP sample count differs from source-only 2022 QC")
    sites=[];bycell=defaultdict(list);qcr=Counter()
    for idhash,meta in keys.items():
        item=data[idhash]
        if item["original_sample_rows"]!=meta["source_metadata_rows"]:
            raise ValueError("Source-only pre-outcome site-year raw TP count differs after numerical opening")
        n=item["original_sample_rows"]
        belowfrac=item["below_LOQ_rows"]/n if n else 1.
        dates=sorted(item["valid_quantified_by_day"])
        months=set(date[:7] for date in dates)
        qualified=(len(dates)>=MIN_VALID_DATES and len(months)>=MIN_VALID_MONTHS and
                  belowfrac<=MAX_BELOWLOQ_FRAC)
        vals=[statistics.median(item["valid_quantified_by_day"][d]) for d in dates] if dates else []
        val=statistics.median(vals) if qualified else None
        if val is not None and not math.isfinite(val):raise ValueError("Nonfinite independently measured station TP median")
        rec={
          "original_station_siteKey_scheme_country_SHA256":idhash,
          "countryCode":meta["country"],"lat_floor":meta["lat_floor"],"lon_floor":meta["lon_floor"],
          "original_source_metadata_sample_rows":meta["source_metadata_rows"],
          "source_actual_numeric_TP_sample_rows":n,
          "numeric_quantified_days":len(dates),
          "numeric_quantified_months":len(months),
          "numeric_LOQ_below_rows":item["below_LOQ_rows"],
          "numeric_LOQ_censor_fraction":belowfrac,
          "numeric_missing_LOQ_flag_rows":item["unknown_flag_rows"],
          "numeric_missing_value_rows":item["original_missing_numeric_rows"],
          "numeric_invalid_result_rows":item["original_nonfinite_or_negative_rows"],
          "numeric_site_support_ge4days_ge3months_LOQfraction_lehalf":qualified,
          "original_EEA2022_TP_station_median_of_daily_medians_mg_P_per_L":val,
          "not_source_specific_station_2022_P_dose":True,
        }
        sites.append(rec)
        qcr["source_original_site_records"]+=1
        if qualified:
            qcr["numeric_TP_quantified_site_years"]+=1
            qcr["numeric_TP_quantified_site_years_"+meta["country"]]+=1
            bycell[(meta["country"],meta["lat_floor"],meta["lon_floor"])].append(val)
        else:
            qcr["numeric_LOQ_or_daily_sampling_qualification_fail_sites"]+=1
    outputcells=[]
    for orig in cells:
        key=(orig["country"],orig["1deg_lat_floor"],orig["1deg_lon_floor"])
        valid=bycell[key]
        outputcells.append({
          "countryCode":key[0],"lat_floor":key[1],"lon_floor":key[2],
          "original_country_2degree_block":orig["original_2deg_spatial_validation_block"],
          "original_source_only_2022_TP_ge4dates_3months_site_support":orig["original_EEA2026_2022_RW_TP_individual_sample_stations_with_ge4dates_3months"],
          "2022_original_WISE6_quantified_TP_site_years_n":len(valid),
          "2022_WISE6_measured_site_median_TP_cell_median_mg_P_per_L":
              statistics.median(valid) if len(valid)>=MIN_SITE_PER_CELL else None,
          "2022_TP_quantified_ge2site_cell_gate_pass":len(valid)>=MIN_SITE_PER_CELL,
          "site_year_equal_weight_not_original_sample_count_weight":True,
        })
    ok=[c for c in outputcells if c["2022_TP_quantified_ge2site_cell_gate_pass"]]
    country=Counter(c["countryCode"] for c in ok)
    blocks=set(c["countryCode"]+"|"+str(c["original_country_2degree_block"]) for c in ok)
    gate=(len(ok)>=MIN_CELLS and country["DE"]>=MIN_DE and country["AT"]>=MIN_AT and len(blocks)>=MIN_BLOCKS)
    report={
        "schema":"gesui_EEA2026_2022_original_WISE6_River_total_P_numeric_mgP_L_site_year_and_country_cell_outcome_QC_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2022_EEA2026_WISE6_DisaggregatedData_ZIP_verified":receipt,
        "original_PREnumeric_2022_TP_site_metadata_source_SHA256":sha(a.original_2022_tp_site_metadata),
        "original_2022_TP_analyte_original_code":TP_CODE,
        "original_2022_TP_unit_original":TP_UNIT,
        "TP_real_2022_site_qualification":"ge4 original quantified dates ge3 quantified months, below LOQ fraction le0.5, all below LOQ values excluded without substitution; each site-year median of daily medians",
        "raw_original_EEA_WISE6_countryTP_national_sample_rows":dict(qc),
        "real_original_WISE6_numeric_TP_station_sample_qc":dict(qcr),
        "original_source_187_sites_all_numeric_and_censor_QC":sites,
        "strict_original_42_Lycaenidae_country_photo_cells_and_measured_TP_outcomes":outputcells,
        "source_original_42_cells_with_ge2_quantified_TP_sites_n":len(ok),
        "source_original_quantified_TP_cells_by_country":dict(country),
        "source_original_2degree_blocks_with_qualified_numerical_TP":len(blocks),
        "numeric_P_water_TP_to_Lycaenidae_photo_geographic_gate_pass":gate,
        "no_wastewater_attribution_and_no_puddling_claim":True,
        "original_WWTP_P_ton_year_does_not_equal_ambient_WISE6_mg_P_per_L":True,
        "original_numeric_TP_site_record_not_independent_butterfly_population":True,
        "no_original_36GB_CSV_or_observer_immaterial_private_data_committed":True,
        "sample_source_numeric_values_only_after_documented_preoutcome_station_gate":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_REAL_2022_WISE6_RIVER_TP_ORIGINAL_NUMERIC_OUTCOME_FIRST_CHECK",
          json.dumps({"sites_source_original":qcr["source_original_site_records"],
            "sites_quantified":qcr["numeric_TP_quantified_site_years"],
            "2022_TP_qualified_cells":len(ok),"country":dict(country),"original_2deg_blocks":len(blocks),
            "gate":gate,"total_below_LOQ_sample_n":sum(x["numeric_LOQ_below_rows"] for x in sites),
            "full_source_rows_read":qc["original_2026_WISE6_full_source_rows_checked"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
