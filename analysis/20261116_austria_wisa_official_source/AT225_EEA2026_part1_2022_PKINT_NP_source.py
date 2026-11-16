#!/usr/bin/env python3
"""True official EEA WISE6 ICM2026 Part1 96,597,294-source-row water
chemistry source-only 2022 TN/TP sample availability for 225 *exact* original
BMLUK Austrian GZÜV PKINT sites near 117 measured-annual-NP UWWTD outfalls.

Do not open river TP/TN values in saved artifact, plant NUTRI, butterfly HNS,
benthic EQR numbers or Sentinel NDVI. Never assign source EEA WWTP annual
tonnes/year to actual river concentrations. Missing WQ sampling != zero.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
import io
import json
import math
import subprocess
import zipfile
from pathlib import Path

import requests

ORIG_SOURCE="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_DisaggregatedData-csv.zip"
TOKEN="ZcJ37WK6oaWonFJ"
EXPECTED_ZIP_BYTES=1639231716
EXPECTED_ZIP_SHA="ff779521749b40b952cce35a4671240c3b6c47d87f32a8fef9240a4ac503a082"
EXPECTED_SOURCE_ROWS=96597294
CROSSWALK_SHA="040d610b6172d165787775406ea50c0f73f34db1db3fcb6c8576898736146dde"
BQE_REL="results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
SCHEMA="GESUI_20261116_EEA_ICM2026_PART1_GZUEV_225_ORIGINAL_PKINT_2022_TNTP_SOURCE_ONLY_V1"
CODES={"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"}
UNITS={"TN":"mg{N}/L","TP":"mg{P}/L"}
ACTUAL_SCHEME="euMonitoringSiteCode"
YEAR=2022
N_DAYS=3
N_SAMPLES=3


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def gitblob(p):
    b=p.read_bytes()
    return hashlib.sha1(("blob "+str(len(b))+"\0").encode()+b).hexdigest()


def load_crosswalk(p):
    if sha(p)!=CROSSWALK_SHA:
        raise RuntimeError("Original signed BMLUK PKINT×outfall 225-row source crosswalk SHA changed")
    with p.open(newline="",encoding="utf-8-sig") as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=225:
        raise RuntimeError("Original PKINT source crosswalk did not contain exactly 225 source pairs")
    ids={}
    outfalls=collections.defaultdict(set)
    for q in rows:
        pid=q["original_official_GZUEV_PKINT"].strip()
        pk=q["original_official_GZUEV_PK"].strip()
        out=q["source_2022_EEA_UWWTD_original_plant_ID"].strip()
        if not pid.startswith("ATFW") or pid!="AT"+pk:
            raise RuntimeError("Original Austrian government PK↔international PKINT identity does not match")
        if pid in ids:
            raise RuntimeError("One official GZÜV PKINT identified more than once in 225 target records")
        ids[pid]={"outfall":out,"geodesic_distance_m":float(q["public_source_point_distance_to_outfall_m"])}
        if not (0<=ids[pid]["geodesic_distance_m"]<=2000):
            raise RuntimeError("Original 2km source preselected station now out of original bound")
        outfalls[out].add(pid)
    if len(ids)!=225 or len(outfalls)!=117:
        raise RuntimeError("Original 225 physical monitoring IDs or 117 annual measured outlet ID counts changed")
    return ids,outfalls


def load_benthic_original(path):
    if gitblob(path)!=BQE_BLOB:
        raise RuntimeError("Original 2022 EEA benthic monitoring-site source SHA changed")
    with path.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.reader(f)
        header=next(reader)
        idx=header.index("site_key")
        cidx=header.index("countryCode")
        result=[]
        n=0
        for row in reader:
            n+=1
            if row[cidx]!="AT":continue
            parts=row[idx].split("|")
            if len(parts)!=3 or parts[0]!="AT" or parts[2]!=ACTUAL_SCHEME:
                raise RuntimeError("Actual EEA benthic Austrian station ID source format changed")
            result.append(parts[1])
    if n!=2475 or len(result)!=68 or len(set(result))!=68:
        raise RuntimeError("Original EEA2022 2475 source rows 68 Austrian field monitoring identities changed")
    return set(result)


def download_original(out):
    for mode in ("anonymous","public_share_basic"):
        kw={"headers":{"User-Agent":"gesui-AT225-official-EEA-WISE6-2026-PKINT-source-QA/20261116",
                       "Accept":"application/zip,application/octet-stream"},
            "stream":True,"timeout":(12,140)}
        if mode=="public_share_basic":
            kw["auth"]=(TOKEN,"")
        with requests.get(ORIG_SOURCE,**kw) as r:
            if r.status_code in (401,403) and mode=="anonymous":continue
            if r.status_code!=200:
                raise RuntimeError("Original public EEA July2026 Part1 source not downloadable status "+str(r.status_code))
            if not r.url.startswith("https://sdi.eea.europa.eu/"):
                raise RuntimeError("Original EEA WISE6 Part1 official host changed")
            if (r.headers.get("Content-Length") and
                int(r.headers["Content-Length"])!=EXPECTED_ZIP_BYTES):
                raise RuntimeError("Original official Part1 source Content-Length SHA pinned support drift")
            digest=hashlib.sha256()
            n=0
            with out.open("wb") as f:
                for chunk in r.iter_content(4*1024*1024):
                    if not chunk:continue
                    n+=len(chunk)
                    if n>EXPECTED_ZIP_BYTES:
                        raise RuntimeError("Original published EEA full source ZIP exceeds fixed size")
                    f.write(chunk);digest.update(chunk)
            if n!=EXPECTED_ZIP_BYTES or digest.hexdigest()!=EXPECTED_ZIP_SHA:
                raise RuntimeError("Original EEA WISE6 source ZIP original bytes or SHA drift")
            return {"official_source_URL":ORIG_SOURCE,"original_zip_SHA256":digest.hexdigest(),
                    "source_verified_bytes":n,"public_download_method":mode,
                    "source_original_full_zip_not_redistributed":True}
    raise RuntimeError("Original EEA WISE6 public source not accessible")


def good_date(text):
    d=str(text).strip()[:10]
    if len(d)!=10 or d[4]!="-" or d[7]!="-":return None
    try:
        date=datetime.date.fromisoformat(d)
    except ValueError:return None
    return date.isoformat() if date.year==YEAR else None


def parse_csv_stream(source_path,sites):
    cnt=collections.Counter()
    errors=collections.Counter()
    by=collections.defaultdict(lambda:{"TN":{},"TP":{}})
    conflicts=set()
    schemes=collections.Counter()
    with zipfile.ZipFile(source_path) as z:
        names=[q.filename for q in z.infolist() if q.filename.lower().endswith(".csv")]
        if len(names)!=1 or z.getinfo(names[0]).compress_type!=9:
            raise RuntimeError("Original WISE6 source CSV DEFLATE64 ZIP schema changed")
        p=subprocess.Popen(["7z","x","-so",str(source_path),names[0]],
                           stdout=subprocess.PIPE,stderr=subprocess.PIPE,bufsize=1024*1024)
        try:
            with p.stdout:
                reader=csv.reader(io.TextIOWrapper(
                    p.stdout,encoding="utf-8-sig",errors="strict",newline=""))
                header=next(reader)
                req=("countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                     "phenomenonTimeSamplingDate","observedPropertyDeterminandCode",
                     "parameterWaterBodyCategory","procedureAnalysedMatrix","resultUom",
                     "resultObservedValue","resultQualityObservedValueBelowLOQ",
                     "sampleIdentifier")
                if any(q not in header for q in req):
                    raise RuntimeError("Original publisher Part1 source column list changed")
                indices={k:header.index(k) for k in req}
                n=0
                for row in reader:
                    n+=1
                    if n%10000000==0:
                        print("GESUI_EEA_AT225_ORIGINAL_PART1_ROWS_READ="+str(n),flush=True)
                    site=row[indices["monitoringSiteIdentifier"]]
                    if row[indices["countryCode"]]!="AT" or site not in sites:
                        continue
                    date=good_date(row[indices["phenomenonTimeSamplingDate"]])
                    if date is None:continue
                    cnt["target_AT_PKINT_original_2022_any_parameter_sample_rows"]+=1
                    scheme=row[indices["monitoringSiteIdentifierScheme"]]
                    schemes[scheme]+=1
                    if scheme!=ACTUAL_SCHEME:
                        errors["source_station_scheme_not_exact_EEA_euMonitoringSiteCode"]+=1
                        continue
                    cnt["target_exact_PKINT_and_scheme_source_rows"]+=1
                    if row[indices["parameterWaterBodyCategory"]]!="RW" or row[indices["procedureAnalysedMatrix"]]!="W":
                        errors["not_river_water_original_source"]+=1
                        continue
                    code=row[indices["observedPropertyDeterminandCode"]]
                    nutrient=next((k for k,v in CODES.items() if v==code),None)
                    if nutrient is None:continue
                    cnt["original_2022_"+nutrient+"_actual_analyte_rows"]+=1
                    if row[indices["resultUom"]]!=UNITS[nutrient]:
                        errors[nutrient+"_unexpected_chemical_unit"]+=1;continue
                    l=row[indices["resultQualityObservedValueBelowLOQ"]].strip().lower()
                    if l in ("true","1"):
                        errors[nutrient+"_original_below_LOQ_explicit"]+=1;continue
                    try:
                        val=float(row[indices["resultObservedValue"]])
                    except (ValueError,TypeError):
                        errors[nutrient+"_non_numeric"]+=1;continue
                    if not math.isfinite(val) or val<0:
                        errors[nutrient+"_negative_or_nonfinite"]+=1;continue
                    if l=="":
                        errors[nutrient+"_missing_LOQ_flag_not_explicit_true"]+=1
                    samp=row[indices["sampleIdentifier"]].strip()
                    key=(date,samp) if samp else (date,format(val,".12g"))
                    previous=by[site][nutrient].get(key)
                    if previous is not None:
                        if abs(previous[0]-val)>1e-10:
                            conflicts.add((site,nutrient))
                            errors[nutrient+"_contradictory_duplicate_original_sample"]+=1
                        else:errors[nutrient+"_exact_duplicate_sample"]+=1
                        continue
                    by[site][nutrient][key]=(val,l)
                    cnt["original_2022_"+nutrient+"_source_qualified_numeric_sample"]+=1
                if n!=EXPECTED_SOURCE_ROWS:
                    raise RuntimeError("Original public Part1 full-source CSV rows drift "+str(n))
            rc=p.wait(timeout=55)
            if rc!=0:
                raise RuntimeError("True source DEFLATE64 stream failed rc="+str(rc))
        finally:
            if p.poll() is None:p.kill();p.wait(timeout=30)
    return by,conflicts,{"full_original_CSV_source_rows":n,
          "original_source_2022_target_rows_by_scheme":dict(schemes),
          "original_source_selected_2022_stages":dict(cnt),
          "original_source_rejected_2022_issues":dict(errors),
          "original_source_sample_key_date_and_source_sampleIdentifier_checked":True,
          "target_numeric_concentrations_not_exported":True}


def site_count(data):
    out={}
    for nu in ("TN","TP"):
        x=data[nu]
        dates={k[0] for k in x}
        out[nu]={"source_valid_samples":len(x),
                 "source_distinct_sample_days":len(dates),
                 "source_sample_calendar_start":min(dates) if dates else None,
                 "source_sample_calendar_end":max(dates) if dates else None,
                 "source_explicit_LOQ_false_sample_count":sum(v[1] in ("false","0") for v in x.values()),
                 "source_blank_LOQ_flag_sample_count":sum(v[1]=="" for v in x.values())}
    return out


def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    site,outfalls=load_crosswalk(a.crosswalk)
    bqe=load_benthic_original(Path(BQE_REL))
    zp=a.output/"TEMP_EEA2026_PART1_OFFICIAL_DO_NOT_REPUBLISH.zip"
    receipt=download_original(zp)
    print("GESUI_EEA_AT225_OFFICIAL_ARCHIVE_SHA="+json.dumps(receipt,sort_keys=True),flush=True)
    raw,conflicts,scan=parse_csv_stream(zp,site)
    zp.unlink()
    summary=collections.Counter()
    per=[]
    eligible_both=set()
    for pk in sorted(site):
        stats=site_count(raw.get(pk,{"TN":{},"TP":{}}))
        flags={}
        for nu in ("TN","TP"):
            flags[nu]=(stats[nu]["source_valid_samples"]>=N_SAMPLES and
                      stats[nu]["source_distinct_sample_days"]>=N_DAYS and
                      (pk,nu) not in conflicts)
            if flags[nu]:summary["source_true_2022_"+nu+"_station_pass"]+=1
        both=flags["TN"] and flags["TP"]
        if both:
            eligible_both.add(pk)
            summary["source_true_2022_both_TNTP_station_pass"]+=1
        record={
            "official_source_PKINT":pk,"source_EEA2022_original_UWWTD_facility":site[pk]["outfall"],
            "source_original_outfall_to_station_geodesic_m":site[pk]["geodesic_distance_m"],
            "original_EEA2022_benthic_field_site_in_original_68_AT":pk in bqe,
            "original_EEA_EEA2022_national_benthos_coordinate_source_anomaly":pk=="ATFW54110117",
            "original_2022_TN_eligible":flags["TN"],"original_2022_TP_eligible":flags["TP"],
            "original_2022_both_TNTP_eligible":both,
            "TN_valid_distinct_original_samples":stats["TN"]["source_valid_samples"],
            "TP_valid_distinct_original_samples":stats["TP"]["source_valid_samples"],
            "TN_valid_distinct_dates":stats["TN"]["source_distinct_sample_days"],
            "TP_valid_distinct_dates":stats["TP"]["source_distinct_sample_days"],
            "source_TN_first_date":stats["TN"]["source_sample_calendar_start"],
            "source_TN_last_date":stats["TN"]["source_sample_calendar_end"],
            "source_TP_first_date":stats["TP"]["source_sample_calendar_start"],
            "source_TP_last_date":stats["TP"]["source_sample_calendar_end"],
            "original_2022_sample_concentration_values_NOT_EXPORTED":True,
            "does_NOT_prove_downstream_hydrologic_receiving_link":True}
        per.append(record)
    outfall_any={site[pk]["outfall"] for pk in eligible_both}
    bqe_any=set(bqe)&eligible_both
    source_site_n_bqe=len(set(site)&bqe)
    result={
        "schema":SCHEMA,"source_gate":"SOURCE_ONLY_2022_WATER_CHEMISTRY_QUALIFIED_NO_CAUSAL_MODEL",
        "original_official_EEA_source_zip":receipt,
        "original_Part1_csv_read_to_end":scan,
        "original_external_225_PKINT_source_crosswalk_SHA256":CROSSWALK_SHA,
        "preselected_original_official_PKINT_station_n":len(site),
        "preselected_original_2022_measured_NP_outfall_n":len(outfalls),
        "original_actual_2022_totalN_sites_pass_3distinct_days":summary["source_true_2022_TN_station_pass"],
        "original_actual_2022_totalP_sites_pass_3distinct_days":summary["source_true_2022_TP_station_pass"],
        "original_2022_both_same_PKINT_TNTP_3days_each":len(eligible_both),
        "original_2022_annual_reported_NP_outfalls_with_at_least_one_2km_WQ_TNTP_station":len(outfall_any),
        "source_outfalls_with_both_chemistry_original_ids":sorted(outfall_any),
        "original_2022_known_EEA_benthic_sample_site_PKINT_overlap_in225":source_site_n_bqe,
        "original_2022_known_EEA_benthic_sample_site_PKINT_with_both_NP_same_station":len(bqe_any),
        "known_EEA2022_benthic_source_coordinate_discrepancy_ATFW54110117_retained":True,
        "strict_same_site_true_benthos_and_2022_WQ_requires_original_river_reach_source":True,
        "national_WISA_independent_benthic_records_not_checked":True,
        "EEA_source_unknown_schemes_not_merged":True,
        "invalid_or_conflicting_station_nutrient_sample_sets":len(conflicts),
        "water_resource_genuine_2022_measured_at_monitor_not_specific_effluent_input":True,
        "no_plant_NUTRI_butterfly_NDVI_EQR_values_read_or_fitted":True,
        "next_full_eco_source_minimum_20_independent_reaches_3_H6_and_10_reference_stations_not_tested":True}
    p=a.output/"AT225_GZUEV_PKINT_2022_WISE6_PART1_SOURCE_ONLY_TNTP_STATIONS.csv"
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(per[0]))
        w.writeheader();w.writerows(per)
    result["derived_source_225_PKINT_stations_csv_SHA256"]=sha(p)
    (a.output/"AT225_GZUEV_2022_WISE6_PART1_CHEMISTRY_ONLY_GATE.json").write_text(
       json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_AT225_GZUEV_2022_TNTP_WISE6_SOURCE="+json.dumps({
        "status":result["source_gate"],
        "original_source_part1_rows":scan["full_original_CSV_source_rows"],
        "original_target_PKINT":len(site),
        "original_target_outfalls":len(outfalls),
        "TN_sites":result["original_actual_2022_totalN_sites_pass_3distinct_days"],
        "TP_sites":result["original_actual_2022_totalP_sites_pass_3distinct_days"],
        "TNTP_same_site":result["original_2022_both_same_PKINT_TNTP_3days_each"],
        "outfalls_with_qualifying_same_station_water":len(outfall_any),
        "known_original_EEA_benthic_station_ids_in225":source_site_n_bqe,
        "known_original_EEA_benthic_and_TNTP_same_station":len(bqe_any),
        "source_artifact_225_stations_SHA":result["derived_source_225_PKINT_stations_csv_SHA256"],
        "no_wastewater_causal_treatment_or_NDVI_biota_fitted":True
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--crosswalk",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
