#!/usr/bin/env python3
"""EEA Waterbase2026 WISE6 original 2022 river TP OBSERVATION STATUS source audit.

This independent WQ metadata inspection reads the original full public ZIP
2026 WISE6 archive SHA, then extracts ONLY three raw status flags (no numeric
resultObservedValue, no photo/observer IDs) for original 187 pre-frozen
WISE6 river TP monitoring stations. It checks site count equivalence with
the source-only original year, no source-effect model. Never interprets
unknown original status codes as accepted or rejected automatically.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,urllib.request,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
SOURCE_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
SOURCE_BYTES=2024284753
SOURCE_ROWS=102569057
SOURCE_MEMBER="waterbase_t_wise6_disaggregateddata.csv"
EXPECTED_MATCHED_SAMPLE_ROWS=2471
FIELDS=("resultObservationStatus","metadata_statusCode","metadata_observationStatus",
        "resultQualityObservedValueBelowLOQ")
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def download(p):
    req=urllib.request.Request(URL,headers={
     "User-Agent":"gesui-original-WISE6-2022-TP-187-station-raw-status-only/1.0",
     "Accept":"application/zip,*/*"})
    h=hashlib.sha256();n=0
    with urllib.request.urlopen(req,timeout=180) as r,p.open("wb") as o:
        if r.status!=200:raise ValueError("Original official EEA2026 source unavailable")
        for b in iter(lambda:r.read(1024*1024),b""):
            n+=len(b)
            if n>SOURCE_BYTES:raise ValueError("Official source exceeds fixed original size")
            h.update(b);o.write(b)
    if n!=SOURCE_BYTES or h.hexdigest()!=SOURCE_SHA:
        raise ValueError("Original independent 2026 WISE6 full ZIP SHA changed")
    return {"download_bytes":n,"original_EEA2026_archive_SHA256":h.hexdigest()}
def sitehash(country,scheme,site):
    return hashlib.sha256((country+"|"+scheme+"|"+site).encode()).hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-2022-TP-sourceonly",type=Path,required=True)
    p.add_argument("--tmpdir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    src=json.loads(a.original_2022_TP_sourceonly.read_text())
    if not src["source_preoutcome_geography_repeated_site_gate_20cells_DE10_AT6_blocks10"]:
        raise ValueError("Original fixed source-only 2022 TP 187 site gate failed")
    sites={}
    for cell in src["original_42_country_1degree_photo_cells_and_WISE6_TP_stations"]:
        for st in cell["original_non_numeric_WISE6_2022_station_sampling_support"]:
            key=st["monitoring_site_ID_with_scheme_sha256"]
            if key in sites:raise ValueError("Original site hash duplicate")
            sites[key]={"country":cell["country"],
             "source_prequality_numeric_TP_sample_rows":st["original_2022_TP_sample_metadata_rows"],
             "original_count":0,"status":{f:Counter() for f in FIELDS}}
    if len(sites)!=187:raise ValueError("Original source 187 station set drift")
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    zipfilelocal=a.tmpdir/"original_EEA2026_WISE6_2022_TP_status_only.zip"
    receipt=download(zipfilelocal)
    original_rows=0;matchrows=0;sampleqc=Counter();source_quality=defaultdict(Counter)
    with zipfile.ZipFile(zipfilelocal) as z:
        if z.namelist()!=[SOURCE_MEMBER]:raise ValueError("Original WISE6 ZIP no longer has expected sole individual CSV")
        with z.open(SOURCE_MEMBER) as r:
            header=r.readline().decode("utf-8-sig").strip().split(";")
            if header[11]!="resultObservedValue" or any(f not in header for f in FIELDS):
                raise ValueError("WISE6 2022 original measured TP source row/quality schema drift")
            index={x:header.index(x) for x in FIELDS}
            for line in r:
                original_rows+=1
                if not line.startswith((b"DE;",b"AT;",b'"DE";',b'"AT";')):continue
                v=line.split(b";",11)
                if len(v)<12:continue
                simple=lambda i:v[i].strip().strip(b'"').decode("utf-8",errors="replace")
                if simple(9)!="2022" or not simple(8).startswith("2022-") or simple(3)!="RW" or simple(5)!="Total phosphorus" or simple(6)!="W":
                    continue
                country=simple(0)
                h=sitehash(country,simple(2),simple(1))
                if h not in sites:continue
                # Full semicolon CSV parse of ONLY the 2471 source-frozen TP
                # samples to safely locate flags past original free-text Remarks.
                rr=list(csv.reader([line.decode("utf-8",errors="replace")],
                     delimiter=";",quotechar='"'))
                if len(rr)!=1 or len(rr[0])!=len(header):
                    sampleqc["source_sample_original_csv_quoting_unreadable"]+=1
                    continue
                row=rr[0]
                matchrows+=1
                st=sites[h]
                st["original_count"]+=1
                for f in FIELDS:
                    value=str(row[index[f]]).strip()
                    st["status"][f][value or "MISSING_OR_BLANK"]+=1
                    source_quality[f][value or "MISSING_OR_BLANK"]+=1
    if original_rows!=SOURCE_ROWS:
        raise ValueError("Original WISE6 all source 102,569,057 rows count changed")
    if matchrows!=EXPECTED_MATCHED_SAMPLE_ROWS:
        raise ValueError(f"Original 2022 numeric TP station 2471 source rows changed {matchrows}")
    diff=[(h,z["source_prequality_numeric_TP_sample_rows"],z["original_count"])
          for h,z in sites.items() if z["source_prequality_numeric_TP_sample_rows"]!=z["original_count"]]
    if diff:raise ValueError("Original pre-numeric 187-site TP sample count mismatch")
    result={
      "schema":"gesui_original_EEA2026_WISE6_2022_riverTP_187_station_observationStatus_quality_codes_QC_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_WISE6_public_source":receipt,
      "original_2022_PREnumeric_TP_187station_source_SHA256":sha(a.original_2022_TP_sourceonly),
      "original_source_WISE6_total_rows_audited":original_rows,
      "original_2022_TP_source187_sites":len(sites),
      "original_2022_TP_source_samples_flagged":matchrows,
      "observation_flag_field_values_overall":{k:dict(v) for k,v in source_quality.items()},
      "original_flag_data_parsing_QC":dict(sampleqc),
      "original_site_quality_flags":[{"site_original_SHA256":h,
                  "countryCode":z["country"],
                  "source_2022_TP_station_sample_rows":z["original_count"],
                  "flag_values_original":{k:dict(v) for k,v in z["status"].items()}
              } for h,z in sites.items()],
      "original_flag_code_semantics_NOT_GUESSED":True,
      "no_numeric_TP_mgP_L_values_parsed_in_this_status_audit":True,
      "no_independent_WWTP_chemical_or_butterfly_effect_fit_in_this_audit":True,
      "no_raw_individual_WISE6_values_or_raw_sample_identifiers_committed":True,
      "next":"Investigate official status code dictionary; if a quality code disqualifies original sample rows, re-evaluate 184 numeric site medians with the same predeclared TP source frame and flag-based sample rules before publication."
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_WISE6_2022_TP_ORIGINAL_STATUS_CODE_AUDIT",
      json.dumps({"source_samples":matchrows,"source_stations":len(sites),
           "original_raw_status_distributions":result["observation_flag_field_values_overall"],
           "no_river_TP_numerical_values_opened":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
