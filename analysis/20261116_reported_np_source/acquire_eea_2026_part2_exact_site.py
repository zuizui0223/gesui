#!/usr/bin/env python3
"""2026 EEA Waterbase ICM Part2: exact 2018/22 original EQR station TN/TP.

Official July 2026 publisher ZIP only, NOT stale Discodata SQL. Download
ONE original 170,422,309-byte aggregated-site CSV ZIP through public
Nextcloud WebDAV (/datashare), compute SHA256, stream CSV only. Strictly
match ORIGINAL EEA EQR station|scheme|year and original 132 reported-positive
station-years; no analysis of EQR/HNS/Plant N traits, no causal SEM.

All original values remain in the GitHub Actions temporary workspace; only
matched-site/year source summary and quality flags are uploaded as artifact.
No third-party bulk file committed or uploaded.
"""
from __future__ import annotations
import argparse
import base64
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

EEA_SHARE="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_AggregatedData-csv.zip"
PUBLIC_SHARE_TOKEN="ZcJ37WK6oaWonFJ"
EXPECTED_ZIP_BYTES=170422309
EEA_SOURCE_NAME="EEA2018_EEA2022_EPA1819_NATIVE_HYBAS7_SOURCE_ONLY_GEOGRAPHY.json"
EEA_STATION_H7_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
H7_PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
BLOBS={
 "EEA18":("0f30bfbc3dac13f6a182c8a910e9899c04778322","2018_measured_upstream_P_t_per_year","2018_measured_upstream_N_t_per_year"),
 "EEA22":("325f199e985d9862ffee5b0edb205c61a7be3833","upstream_2022_EEA_P_measured_t_per_yr","upstream_2022_EEA_N_measured_t_per_yr")
}
CODES={"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"}
UNITS={"TN":"mg{N}/L","TP":"mg{P}/L"}
YEARS={"EEA18":2018,"EEA22":2022}
PERIODS={"EEA18":"2017-2019","EEA22":"2020-2022"}
RUN_SOURCE="GESUI_20261116_EEA_ICM2026_OFFICIAL_CSV_P2_SOURCE_ONLY"
UA="gesui-scientific-EEA-2026-annual-site-TNTP-locked-2018-2022/1"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def gitblob(raw):
    return hashlib.sha1(("blob "+str(len(raw))+"\0").encode()+raw).hexdigest()

def one(path,name):
    found=list(path.rglob(name))
    if len(found)!=1:raise RuntimeError("Original EEA source missing or ambiguous: "+name)
    return found[0]

def get_git_original_blob(blobsha,token):
    url="https://api.github.com/repos/zuizui0223/gesui/git/blobs/"+blobsha
    headers={"User-Agent":UA,"Accept":"application/vnd.github+json"}
    if token:headers["Authorization"]="Bearer "+token
    r=requests.get(url,headers=headers,timeout=(10,75))
    r.raise_for_status()
    ob=r.json()
    if ob.get("encoding")!="base64":raise RuntimeError("Original Git blob content not base64")
    raw=base64.b64decode(ob["content"].replace("\n",""))
    if gitblob(raw)!=blobsha:raise RuntimeError("Original EEA reported station source Git blob SHA changed")
    obj=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    if len(obj) not in (2475,2634):raise RuntimeError("Original EEA site rows not 2475/2634")
    ids={x["site_key"]:x for x in obj}
    if len(ids)!=len(obj):raise RuntimeError("Original EEA reported station IDs duplicated")
    return ids

def prechemistry_source(args):
    if sha(args.panel)!=H7_PANEL_SHA:raise RuntimeError("Original 454 H7-panel SHA drift")
    with args.panel.open(newline="",encoding="utf-8-sig") as f:
        r=list(csv.DictReader(f))
    if (len(r),len({x["h"] for x in r}),len({x["h6"] for x in r}))!=(454,232,90):
        raise RuntimeError("Original 454/232/90 original basin geography changed")
    keymap={(x["h"],x["period"]):(x["h6"],x["country"]) for x in r}
    if len(keymap)!=454:raise RuntimeError("Duplicate original H7×period")
    original=one(args.stations,EEA_SOURCE_NAME)
    if sha(original)!=EEA_STATION_H7_SHA:raise RuntimeError("Original 6330 EEA/EPA geographic source SHA drift")
    raw=json.loads(original.read_text())
    rows=raw.get("candidate_source_site_period_rows")
    if len(rows)!=6330:raise RuntimeError("Original 6330 EEA/EPA source site-rows changed")
    station_reports={}
    for panel,(sha_,_,__) in BLOBS.items():
        station_reports[panel]=get_git_original_blob(sha_,args.github_token)
        if len(station_reports[panel])!=(2634 if panel=="EEA18" else 2475):
            raise RuntimeError("Exact original 2018/22 EEA source station rows changed")
    station_source=[]
    for q in rows:
        pan=q.get("panel")
        if pan not in YEARS:continue
        h=str(q.get("HYBAS_ID_L07"))
        period=q.get("period")
        if period!=PERIODS[pan] or (h,period) not in keymap:continue
        native=q.get("native_id")
        if native!=q.get("physical_id") or not isinstance(native,str):
            raise RuntimeError("Source original site identifiers not native physical")
        bits=native.split("|")
        if len(bits)!=3 or any(not v for v in bits):raise RuntimeError("Original station ID not CC|identifier|scheme")
        if q.get("country")!=bits[0]:
            raise RuntimeError("Original station-country mismatch")
        report=station_reports[pan].get(native)
        if report is None:raise RuntimeError("EEA original reported-source station not found")
        _,p,n=BLOBS[pan]
        pv=float(report[p]);nv=float(report[n])
        if any(not math.isfinite(v) or v<0 for v in (pv,nv)):
            raise RuntimeError("Invalid original upstream yearly reported source mass")
        station_source.append({
            "native_site":native,"scheme":bits[2],"site":bits[1],"country":bits[0],
            "calendar_year":YEARS[pan],"period":period,"h":h,
            "h6":keymap[(h,period)][0],
            "upstream_annual_report_positive":bool(pv>0 or nv>0)})
    if len(station_source)!=3525 or len({x["native_site"] for x in station_source})!=2485:
        raise RuntimeError("Original 3525 source EEA site-year or 2485 station unique count drift")
    if (sum(x["upstream_annual_report_positive"] for x in station_source),
        len({(x["h"],x["period"]) for x in station_source if x["upstream_annual_report_positive"]}))!=(132,56):
        raise RuntimeError("Original 132 report-positive site-year/56 H7period selection drift")
    target={}
    for row in station_source:
        k=(row["site"],row["scheme"],row["calendar_year"])
        if k in target:
            raise RuntimeError("Duplicate monitoringSiteIdentifier+Scheme+year in original station assignment")
        target[k]=row
    if len(target)!=3525:raise RuntimeError("Original exact station target keys duplicated")
    return target,station_source

def download_part2(dest):
    h={"User-Agent":UA,"Accept":"application/zip,application/octet-stream"}
    for mode in ("anonymous","public_share_auth"):
        kw=dict(headers=h,timeout=(15,130),stream=True,allow_redirects=True)
        if mode=="public_share_auth":kw["auth"]=(PUBLIC_SHARE_TOKEN,"")
        with requests.get(EEA_SHARE,**kw) as r:
            if r.status_code in (401,403) and mode=="anonymous":continue
            if r.status_code!=200:
                raise RuntimeError("Official EEA Part2 public source GET "+str(r.status_code)+" mode="+mode)
            if r.url.split("/")[2] not in ("sdi.eea.europa.eu",):
                raise RuntimeError("Publisher download redirected off official host")
            if r.headers.get("Content-Length"):
                length=int(r.headers["Content-Length"])
                if length!=EXPECTED_ZIP_BYTES:
                    raise RuntimeError("Publisher original ZIP differs from independently source-listed size")
            count=0
            dig=hashlib.sha256()
            with dest.open("wb") as out:
                for chunk in r.iter_content(chunk_size=1048576):
                    if not chunk:continue
                    count+=len(chunk)
                    if count>EXPECTED_ZIP_BYTES:
                        raise RuntimeError("Source returned more than original public file listing size")
                    dig.update(chunk);out.write(chunk)
            if count!=EXPECTED_ZIP_BYTES:
                raise RuntimeError("Source ZIP was incomplete: "+str(count))
            if not zipfile.is_zipfile(dest):
                raise RuntimeError("Publisher returned non-ZIP HTML despite 200 and full length")
            return {"publisher_public_share_file":"WISE6_AggregatedData-csv.zip",
                    "public_source_listed_bytes":EXPECTED_ZIP_BYTES,
                    "original_verified_download_bytes":count,
                    "official_2026_original_part2_ZIP_SHA256":dig.hexdigest(),
                    "source_access_mode":mode,
                    "original_publisher_url":EEA_SHARE,
                    "source_checked_byte_length":True}
    raise RuntimeError("Neither anonymous nor public token source access available")

def parse_int(v):
    try:return int(float(v))
    except (TypeError,ValueError):return None

def parse_strict(v,nutrient):
    try:
        measured=float(v["resultMeanValue"])
        number=parse_int(v["resultNumberOfSamples"])
        if not (math.isfinite(measured) and measured>=0 and number is not None and number>=3):
            return None,"VALUE_OR_SAMPLE_COUNT_NOT_ELIGIBLE"
        if str(v["resultUom"])!=UNITS[nutrient]:
            return None,"NONCANONICAL_UNIT"
        if str(v.get("resultQualityMeanBelowLOQ","")).strip().lower() in ("true","1"):
            return None,"BELOW_LOQ"
        return {"mean_mgL":measured,"samples":number,
                "units":str(v["resultUom"])},None
    except (ValueError,TypeError,KeyError):
        return None,"BAD_NUMBER_OR_REQUIRED_COLUMN"

def aggregate_2026(zip_path,target):
    with zipfile.ZipFile(zip_path) as zp:
        info=[{"name":v.filename,"compressed_bytes":v.compress_size,"uncompressed_bytes":v.file_size}
            for v in zp.infolist() if not v.is_dir()]
        names=[f.filename for f in zp.infolist() if f.filename.lower().endswith(".csv")]
        if len(names)!=1:
            # The publisher ZIP must be inventoried before assuming single table.
            return {"status":"STOP_UNKNOWN_ARCHIVE_STRUCTURE","ZIP_member_metadata":info},{}
        filename=names[0]
        with zp.open(filename) as h:
            stream=io.TextIOWrapper(h,encoding="utf-8-sig",errors="strict",newline="")
            head=stream.readline()
            if len(head)>50000:raise RuntimeError("Original CSV header too large")
            delim=";" if head.count(";")>head.count(",") else ","
            cols=next(csv.reader([head],delimiter=delim))
            needed=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                "phenomenonTimeReferenceYear","parameterWaterBodyCategory",
                "observedPropertyDeterminandCode","procedureAnalysedMatrix",
                "resultUom","resultNumberOfSamples","resultQualityMeanBelowLOQ","resultMeanValue"]
            absent=[c for c in needed if c not in cols]
            if absent:
                return {"status":"STOP_REQUIRED_ORIGINAL_COLUMNS_MISSING","archive_members":info,
                        "column_names":cols[:70],"missing_columns":absent},{}
            reader=csv.DictReader(stream,fieldnames=cols,delimiter=delim)
            row_count=0
            bysite=collections.defaultdict(lambda:{"TN":[],"TP":[]})
            stage=collections.Counter()
            status=collections.Counter()
            yrs=collections.Counter()
            matched_countries=collections.Counter()
            for row in reader:
                row_count+=1
                y=parse_int(row.get("phenomenonTimeReferenceYear"))
                if y not in (2018,2022):continue
                yrs[str(y)]+=1
                k=(row.get("monitoringSiteIdentifier"),row.get("monitoringSiteIdentifierScheme"),y)
                old=target.get(k)
                if old is None:continue
                stage["station_year_identified"]+=1
                if str(row.get("countryCode"))!=old["country"]:
                    status["COUNTRY_MISMATCH"]+=1
                    continue
                stage["country_matched"]+=1
                if row.get("parameterWaterBodyCategory")!="RW" or row.get("procedureAnalysedMatrix")!="W":
                    status["NONRIVER_OR_NONWATER"]+=1;continue
                if row.get("observedPropertyDeterminandCode") not in CODES.values():continue
                nutrient="TN" if row["observedPropertyDeterminandCode"]==CODES["TN"] else "TP"
                stage[nutrient+"_exact_code_and_site_year"]+=1
                good,why=parse_strict(row,nutrient)
                if why:
                    status[nutrient+"_"+why]+=1;continue
                # True source observations are stored only for matched EEA
                # sites and not in Git. They will remain an Actions artifact.
                bysite[k][nutrient].append(good)
                matched_countries[old["country"]]+=1
            if row_count<1000:
                raise RuntimeError("Official EEA original 2026 Part2 CSV implausibly small")
            return {
                "status":"SOURCE_CSV_READ_TO_END",
                "original_csv_member":filename,
                "ZIP_member_metadata":info[:20],
                "source_header":cols,
                "original_total_source_csv_rows_scanned":row_count,
                "source_selected_2018_2022_rows_by_year":dict(yrs),
                "original_exact_station_source_gates":dict(stage),
                "source_ineligible_reasons":dict(status),
                "countries_with_qualified_TN_or_TP":dict(matched_countries),
                "original_full_source_sorted_not_assumed":True
            },bysite

def main(args):
    target,original=prechemistry_source(args)
    args.output.mkdir(parents=True,exist_ok=True)
    path=args.output/"PUBLISHER_2026_TEMP_AGGREGATED_ORIGINAL.zip"
    if args.schema_only:
        raise RuntimeError("schema_only mode not supported: original ZIP streaming identity required")
    source=download_part2(path)
    print("EEA2026_SOURCE_ZIP_SHA="+json.dumps(source,ensure_ascii=False,sort_keys=True),flush=True)
    scan,values=aggregate_2026(path,target)
    if scan["status"]!="SOURCE_CSV_READ_TO_END":
        result={"schema":RUN_SOURCE,"download":source,"source_parse":scan,
                "same_site_year_TNTP_status":"STOP_SCHEMA_NOT_CONFIRMED",
                "no_biological_outcome_fitted":True}
        (args.output/"EEA2026_EXACT_STATION_YEAR_TNTP_SOURCE_QA.json").write_text(
            json.dumps(result,indent=2,ensure_ascii=False)+"\n")
        path.unlink()
        print("EEA2026_EXACT_SITE_CHEMISTRY_GATE="+json.dumps(result,ensure_ascii=False)[:22000],flush=True)
        return
    records=[]
    stat=collections.Counter()
    for key,old in sorted(target.items()):
        got=values.get(key,{"TN":[],"TP":[]})
        tn=got["TN"];tp=got["TP"]
        row={"h":old["h"],"h6":old["h6"],
             "country":old["country"],"period":old["period"],
             "year":old["calendar_year"],
             "monitoringSiteIdentifier":old["site"],
             "monitoringSiteIdentifierScheme":old["scheme"],
             "original_site_reported_upstream_annual_N_or_P_positive":old["upstream_annual_report_positive"],
             "qualified_exact_TN_count":len(tn),"qualified_exact_TP_count":len(tp),
             "strict_single_TN":len(tn)==1,
             "strict_single_TP":len(tp)==1,
             "strict_same_station_year_TNTP":len(tn)==1 and len(tp)==1,
             "source_flag_NOT_zero_chemistry":True}
        if len(tn)==1:row["source_original_TN_mgN_L"]=tn[0]["mean_mgL"];row["n_TN_samples"]=tn[0]["samples"]
        else:row["source_original_TN_mgN_L"]=None;row["n_TN_samples"]=None
        if len(tp)==1:row["source_original_TP_mgP_L"]=tp[0]["mean_mgL"];row["n_TP_samples"]=tp[0]["samples"]
        else:row["source_original_TP_mgP_L"]=None;row["n_TP_samples"]=None
        records.append(row)
        if row["strict_single_TN"]:stat["TN"]+=1
        if row["strict_single_TP"]:stat["TP"]+=1
        if row["strict_same_station_year_TNTP"]:
            stat["BOTH"]+=1
            if old["upstream_annual_report_positive"]:stat["BOTH_REPORT_POSITIVE"]+=1
    valid=[r for r in records if r["strict_same_station_year_TNTP"]]
    valid_positive=[r for r in valid if r["original_site_reported_upstream_annual_N_or_P_positive"]]
    h7p=set((r["h"],r["period"]) for r in valid_positive)
    support={"H7_periods":len(h7p),"H7":len({x[0] for x in h7p}),
             "H6":len({r["h6"] for r in valid_positive}),
             "country_groups":len({r["country"] for r in valid_positive})}
    required={"H7_periods":30,"H7":20,"H6":10,"country_groups":3}
    passgate=all(support[x]>=required[x] for x in required)
    summary={"schema":RUN_SOURCE,
             "official_publisher_2026_original":source,
             "actual_annual_CSV_source":scan,
             "n_original_genuine_EEA_site_years":len(records),
             "n_original_distinct_sites":2485,
             "original_reported_positive_site_years":132,
             "original_reported_positive_H7_periods":56,
             "strict_TN_site_years":stat["TN"],
             "strict_TP_site_years":stat["TP"],
             "strict_same_site_year_both":stat["BOTH"],
             "strict_same_site_year_both_and_original_upstream_report_positive":stat["BOTH_REPORT_POSITIVE"],
             "source_supported_upstream_reported_positive_TNTP":support,
             "precommitted_source_eligibility_threshold":required,
             "source_joint_2018_2022_station_year_TNTP_gate_pass":passgate,
             "site_year_TNTP_equals_ambient_river_chemistry_NOT_effluent_delivery":True,
             "source_stations_in_H7_basin_not_exact_GBIF_flora_patch":True,
             "no_larval_host_consumption_or_leaf_nitrogen_observed":True,
             "no_EQR_HNS_NUTRI_W_outcomes_opened":True,
             "no_causal_W_to_Plant_to_Butterfly_model_fitted":True,
             "completeness_includes_AGGREGATED_ONLY_DisaggregatedData_missingness_unknown":True,
             "no_unmatched_report_as_pristine_controls":True}
    (args.output/"EEA2026_EXACT_STATION_YEAR_TNTP_SOURCE_QA.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    with (args.output/"EEA2026_ONLY_MATCHED_3525_STATION_YEAR_TNTP_SOURCE.csv").open(
            "w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]))
        w.writeheader();w.writerows(records)
    # Do not upload the original 170MB archive as an Actions artifact.
    path.unlink()
    print("EEA2026_EXACT_SITE_CHEMISTRY_GATE="+json.dumps({
        "source_ZIP_SHA256":source["official_2026_original_part2_ZIP_SHA256"],
        "original_csv_row_count":scan["original_total_source_csv_rows_scanned"],
        "original_site_years":len(records),"strict_TN":stat["TN"],"strict_TP":stat["TP"],
        "strict_both":stat["BOTH"],"strict_both_report_positive":stat["BOTH_REPORT_POSITIVE"],
        "positive_support":support,"source_gate_pass":passgate,
        "original_CSV_header":scan["source_header"][:65],
        "original_year_rows":scan["source_selected_2018_2022_rows_by_year"],
        "not_fitted":True},sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--stations",required=True,type=Path)
    p.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    p.add_argument("--output",required=True,type=Path)
    p.add_argument("--github-token",default=None)
    p.add_argument("--schema-only",action="store_true")
    main(p.parse_args())
