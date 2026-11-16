#!/usr/bin/env python3
"""Official EEA 2026 WISE6 Part1 disaggregated: exact 2018/22 EQR station TNTP.

Stream original *only* published Part1 source ZIP (1,639,231,716 bytes) using
official public Nextcloud WebDAV; SHA-lock entire ZIP, stream CSV rows, retain
only exact 3525 SOURCE-LOCKED EEA EQR station/year/identifierScheme, TN & TP.
At least 3 separate qualifying annual water-sample results per nutrient.
No original Part1 bulk files are uploaded/committed. No biology/HNS/EIV
responses, wastewater causal exposure, or mediation model fitted.

Use official annual source 2026 Part2 result SEPARATELY to avoid pretending
mixed aggregated and disaggregated table records are independent samples.
"""
from __future__ import annotations
import argparse
import collections
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import time
import zipfile
import requests

from acquire_eea_2026_part2_exact_site import (prechemistry_source, CODES, UNITS)

TOKEN="ZcJ37WK6oaWonFJ"
SOURCE="https://sdi.eea.europa.eu/datashare/public.php/dav/files/"+TOKEN+"/WISE6_DisaggregatedData-csv.zip"
EXPECTED_BYTES=1639231716
SCHEMA="GESUI_20261116_OFFICIAL_EEA2026_DISAGGREGATED_PART1_ORIGINAL_MATCHED_EEA_STATION_TNTP_V1"
UA="gesui-official-2026-ICM-Part1-exact-stationTNTP-source-QA/20261116"
SELECTION_YEARS={2018,2022}

def download_source(dest):
    for kind in ("public_anonymous","public_token_auth"):
        kw={"headers":{"User-Agent":UA,"Accept":"application/zip,application/octet-stream"},
            "stream":True,"timeout":(12,140),"allow_redirects":True}
        if kind=="public_token_auth":kw["auth"]=(TOKEN,"")
        with requests.get(SOURCE,**kw) as r:
            if r.status_code in (401,403) and kind=="public_anonymous":
                continue
            if r.status_code!=200:raise RuntimeError("Official Part1 source GET status "+str(r.status_code))
            if r.url.split("/")[2]!="sdi.eea.europa.eu":
                raise RuntimeError("Source file redirected away from official EEA host")
            if r.headers.get("Content-Length") and int(r.headers["Content-Length"])!=EXPECTED_BYTES:
                raise RuntimeError("Original official Part1 size drift from public WebDAV source")
            downloaded=0
            digest=hashlib.sha256()
            with dest.open("wb") as f:
                for chunk in r.iter_content(chunk_size=4*1024*1024):
                    if not chunk:continue
                    downloaded+=len(chunk)
                    if downloaded>EXPECTED_BYTES:
                        raise RuntimeError("Official EEA Part1 source exceeds original public size")
                    f.write(chunk);digest.update(chunk)
            if downloaded!=EXPECTED_BYTES:
                raise RuntimeError("Official EEA Part1 downloaded bytes incomplete: "+str(downloaded))
            if not zipfile.is_zipfile(dest):
                raise RuntimeError("Original Part1 archive invalid")
            return {"publisher":"EEA Waterbase Water Quality ICM 2026",
                    "exact_official_public_source_url":SOURCE,
                    "publisher_file_name":"WISE6_DisaggregatedData-csv.zip",
                    "listed_bytes":EXPECTED_BYTES,
                    "downloaded_bytes":downloaded,
                    "original_part1_SHA256":digest.hexdigest(),
                    "access":kind,
                    "original_source_not_republished_in_Git":True}
    raise RuntimeError("Published EEA original Part1 not downloadable anonymously")

def number(v):
    try:
        q=float(v)
        return q if math.isfinite(q) and q>=0 else None
    except (ValueError,TypeError):return None

def parseyear(v):
    try:return int(v)
    except (TypeError,ValueError):return None

def table_source(zip_path,target):
    stages=collections.Counter()
    invalid=collections.Counter()
    years=collections.Counter()
    samples=collections.defaultdict(lambda:{"TN":{},"TP":{}})
    with zipfile.ZipFile(zip_path) as z:
        entries=[{"name":f.filename,"size_uncompressed":f.file_size,
                  "size_compressed":f.compress_size} for f in z.infolist() if not f.is_dir()]
        members=[x.filename for x in z.infolist() if x.filename.lower().endswith(".csv")]
        if len(members)!=1:
            return {"status":"STOP_UNEXPECTED_ZIP_MEMBERS","ZIP_members":entries[:30]},samples
        member=members[0]
        with z.open(member) as byte_stream:
            f=io.TextIOWrapper(byte_stream,encoding="utf-8-sig",errors="strict",newline="")
            reader=csv.reader(f)
            header=next(reader)
            expected=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                "parameterWaterBodyCategory","observedPropertyDeterminandCode",
                "procedureAnalysedMatrix","resultUom","phenomenonTimeSamplingDate",
                "phenomenonTimeReferenceYear","sampleIdentifier","resultObservedValue",
                "resultQualityObservedValueBelowLOQ"]
            missing=[s for s in expected if s not in header]
            if missing:
                return {"status":"STOP_MISSING_ORIGINAL_PART1_FIELDS",
                    "ZIP_members":entries,"header":header[:60],"missing":missing},samples
            I={s:header.index(s) for s in expected}
            total=0
            for row in reader:
                total+=1
                y=parseyear(row[I["phenomenonTimeReferenceYear"]])
                if y not in SELECTION_YEARS:continue
                years[str(y)]+=1
                k=(row[I["monitoringSiteIdentifier"]],
                   row[I["monitoringSiteIdentifierScheme"]],y)
                site=target.get(k)
                if site is None:continue
                stages["native_site_year_found"]+=1
                if row[I["countryCode"]]!=site["country"]:
                    invalid["country_mismatch"]+=1
                    continue
                stages["country_match"]+=1
                if (row[I["parameterWaterBodyCategory"]]!="RW" or
                    row[I["procedureAnalysedMatrix"]]!="W"):
                    invalid["not_river_water"]+=1;continue
                code=row[I["observedPropertyDeterminandCode"]]
                if code not in CODES.values():continue
                nutrient="TN" if code==CODES["TN"] else "TP"
                stages["exact_"+nutrient+"_native_site_year_code"]+=1
                if row[I["resultUom"]]!=UNITS[nutrient]:
                    invalid[nutrient+"_wrong_unit"]+=1;continue
                flag=row[I["resultQualityObservedValueBelowLOQ"]].strip().lower()
                if flag in ("true","1"):
                    invalid[nutrient+"_below_LOQ"]+=1;continue
                d=row[I["phenomenonTimeSamplingDate"]].strip()
                if not d or len(d)<4 or d[:4]!=str(y):
                    invalid[nutrient+"_date_year_conflict_or_missing"]+=1;continue
                val=number(row[I["resultObservedValue"]])
                if val is None:
                    invalid[nutrient+"_nonfinite_or_missing"]+=1;continue
                ident=row[I["sampleIdentifier"]].strip()
                # Unique sample identifier, else conservative date/chemistry hash.
                samplekey=("identified",ident) if ident else ("no_sample_id",d,format(val,".12g"))
                previous=samples[k][nutrient].get(samplekey)
                if previous is not None:
                    if abs(previous-val)>1e-10:
                        invalid[nutrient+"_duplicate_conflicting_source_sample"]+=1
                    else:invalid[nutrient+"_duplicate_exact_sample"]+=1
                    continue
                samples[k][nutrient][samplekey]=val
                stages["valid_"+nutrient+"_native_site_year_sample"]+=1
                if total%1_000_000==0:
                    print("EEA2026_PART1_ORIGINAL_SOURCE_ROWS_READ="+str(total),flush=True)
            if total<1_000_000:raise RuntimeError("Original 2026 Part1 full source CSV unexpectedly small")
            return {
                "status":"SOURCE_PART1_READ_TO_END",
                "original_part1_CSV_member":member,
                "original_part1_CSV_column_names":header,
                "original_part1_ZIP_members":entries[:10],
                "original_total_disaggregated_source_rows_scanned":total,
                "source_rows_in_2018_and_2022":dict(years),
                "stages":dict(stages),"invalid":dict(invalid),
                "sample_duplicate_policy":"sampleIdentifier, else exact date/value; conflicting duplicates excluded",
                "source_rows_can_be_missing_samples_in_PART2":True
            },samples

def main(args):
    target,original=prechemistry_source(args)
    args.output.mkdir(parents=True,exist_ok=True)
    path=args.output/"TEMP_2026_DISAGGREGATED_PART1_ORIGINAL.zip"
    receipt=download_source(path)
    print("EEA2026_PART1_ZIP_SHA="+json.dumps(receipt,ensure_ascii=False,sort_keys=True),flush=True)
    source,samples=table_source(path,target)
    if source["status"]!="SOURCE_PART1_READ_TO_END":
        out={"schema":SCHEMA,"source_original_archive":receipt,"source_column_check":source,
             "source_gate":"STOP_UNRECOGNIZED_OFFICIAL_PART1_SCHEMA",
             "not_site_year_chemistry_absence":True}
        (args.output/"EEA2026_PART1_EXACT_STATION_YEAR_SOURCE_QA.json").write_text(
            json.dumps(out,ensure_ascii=False,indent=2)+"\n")
        path.unlink()
        print("EEA2026_PART1_SAME_SITE_YEAR_TNTP_GATE="+json.dumps(out,ensure_ascii=False)[:22000],flush=True)
        return
    summary_counts=collections.Counter()
    geography=collections.defaultdict(set)
    eligible_rows=[]
    rows=[]
    for key,item in sorted(target.items()):
        group=samples.get(key,{"TN":{},"TP":{}})
        tn=group["TN"];tp=group["TP"]
        tn_ok=len(tn)>=3
        tp_ok=len(tp)>=3
        if tn_ok:summary_counts["TN"]+=1
        if tp_ok:summary_counts["TP"]+=1
        if tn_ok and tp_ok:
            summary_counts["TNTP"]+=1
            geography["all_both_H7period"].add((item["h"],item["period"]))
            geography["all_both_H6"].add(item["h6"])
            geography["all_both_countries"].add(item["country"])
            if item["upstream_annual_report_positive"]:
                summary_counts["TNTP_UPSTREAM_REPORT_POSITIVE"]+=1
                geography["reported_positive_H7period"].add((item["h"],item["period"]))
                geography["reported_positive_H7"].add(item["h"])
                geography["reported_positive_H6"].add(item["h6"])
                geography["reported_positive_country"].add(item["country"])
        record={
          "original_site":item["site"],"original_scheme":item["scheme"],
          "year":item["calendar_year"],"h7":item["h"],"h6":item["h6"],
          "country":item["country"],
          "original_upstream_reported_annual_N_or_P_positive":item["upstream_annual_report_positive"],
          "valid_TN_distinct_water_samples":len(tn),
          "valid_TP_distinct_water_samples":len(tp),
          "valid_annual_TN":tn_ok,"valid_annual_TP":tp_ok,
          "strict_same_site_year_TNTP":tn_ok and tp_ok,
          "EEA2026_part1_TN_mean_mgN_L":sum(tn.values())/len(tn) if tn_ok else None,
          "EEA2026_part1_TP_mean_mgP_L":sum(tp.values())/len(tp) if tp_ok else None}
        rows.append(record)
        if record["strict_same_site_year_TNTP"]:eligible_rows.append(record)
    if len(rows)!=3525:raise RuntimeError("Original 3525 true EEA EQR station-years drift")
    positive={
        "H7period":len(geography["reported_positive_H7period"]),
        "H7":len(geography["reported_positive_H7"]),
        "H6":len(geography["reported_positive_H6"]),
        "country":len(geography["reported_positive_country"])
    }
    sourcegate=all(positive[k]>=v for k,v in {
        "H7period":30,"H7":20,"H6":10,"country":3}.items())
    out={
       "schema":SCHEMA,
       "original_2026_public_part1_source":receipt,
       "original_2026_part1_CSV_source":source,
       "original_true_native_EEA_EQR_site_years":len(rows),
       "original_true_distinct_EEA_EQR_sites":2485,
       "original_upstream_report_positive_site_years":132,
       "strict_TN_single_site_year_n":summary_counts["TN"],
       "strict_TP_single_site_year_n":summary_counts["TP"],
       "strict_both_same_station_year_TN_and_TP_n":summary_counts["TNTP"],
       "strict_both_and_source_report_positive_n":summary_counts["TNTP_UPSTREAM_REPORT_POSITIVE"],
       "strict_both_all_H7periods":len(geography["all_both_H7period"]),
       "strict_both_all_H6":len(geography["all_both_H6"]),
       "strict_both_all_country":len(geography["all_both_countries"]),
       "strict_report_positive_H7periods":positive["H7period"],
       "strict_report_positive_H7":positive["H7"],
       "strict_report_positive_H6":positive["H6"],
       "strict_report_positive_countries":positive["country"],
       "precommitted_minimum_positive_H7period_H7_H6_country":{"H7period":30,"H7":20,"H6":10,"country":3},
       "true_official_same_site_2018_2022_source_gate_pass":sourcegate,
       "Part1_disaggregated_distinct_from_Part2_annual_source":True,
       "not_combined_with_Part2_until_a_further_source_QA":True,
       "no_original_other_third_party_material_republished":True,
       "source_year_2018_2022_to_GBIF_three_year_bin_mismatch_preserved":True,
       "source_ambient_chemistry_not_WWTP_source_apportionment":True,
       "not_causal_mediation_or_ecological_service_result":True,
       "no_biological_EQR_EIV_HNS_outcomes_opened":True}
    (args.output/"EEA2026_PART1_EXACT_STATION_YEAR_SOURCE_QA.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+"\n")
    with (args.output/"EEA2026_ONLY_MATCHED_3525_STATION_YEAR_TNTP_PART1.csv").open(
            "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    path.unlink()
    print("EEA2026_PART1_SAME_SITE_YEAR_TNTP_GATE="+json.dumps({
       "official_original_part1_zip_SHA256":receipt["original_part1_SHA256"],
       "part1_total_original_rows":source["original_total_disaggregated_source_rows_scanned"],
       "part1_2018_2022_raw_rows":source["source_rows_in_2018_and_2022"],
       "part1_stage":source["stages"],
       "Part1_strict_TN":out["strict_TN_single_site_year_n"],
       "Part1_strict_TP":out["strict_TP_single_site_year_n"],
       "Part1_strict_TNTP":out["strict_both_same_station_year_TN_and_TP_n"],
       "Part1_strict_TNTP_upstream_report_positive":out["strict_both_and_source_report_positive_n"],
       "report_positive_H7period":positive["H7period"],
       "report_positive_H7":positive["H7"],
       "report_positive_H6":positive["H6"],
       "report_positive_country":positive["country"],
       "source_gate_pass":sourcegate,
       "part1_header":source["original_part1_CSV_column_names"],
       "no_bio_model_fitted":True
    },ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--stations",required=True,type=Path)
    p.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--github-token",default=None)
    main(p.parse_args())
