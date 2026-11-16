#!/usr/bin/env python3
"""Original EEA2026 annual Part2 2022 original GZÜV 225 PKINT site-only
total N/P availability, separate from Part1 individual sample source.
No original numeric nutrient means, EQR, Sentinel greenness or butterfly
traits exported/analyzed. Exact station IDs, year, country and scheme.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,zipfile
from pathlib import Path
import requests

SCHEMA="GESUI_20261116_EEA2026_PART2_AT225_PKINT_2022_TNTP_SOURCE_ONLY_V1"
URL="https://sdi.eea.europa.eu/datashare/public.php/dav/files/ZcJ37WK6oaWonFJ/WISE6_AggregatedData-csv.zip"
EXPECTED_SIZE=170422309
EXPECTED_SHA="315396800d76e46db415164851fb4f99765b7ec2654828fc42782d902e679a7a"
EXPECTED_ROWS=6450192
CROSSWALK_SHA="040d610b6172d165787775406ea50c0f73f34db1db3fcb6c8576898736146dde"
CODES={"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"}
UNITS={"TN":"mg{N}/L","TP":"mg{P}/L"}


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def sites(path):
    if sha(path)!=CROSSWALK_SHA:raise RuntimeError("Original government signed PKINT target CSV source changed")
    with path.open(newline="",encoding="utf-8-sig") as f:rows=list(csv.DictReader(f))
    if len(rows)!=225:raise RuntimeError("No original 225 unique national sites")
    by={q["original_official_GZUEV_PKINT"]:q for q in rows}
    if len(by)!=225 or len({q["source_2022_EEA_UWWTD_original_plant_ID"] for q in rows})!=117:
        raise RuntimeError("Original crosswalk station/outfall cardinalities changed")
    if any(not k.startswith("ATFW") for k in by):
        raise RuntimeError("GZÜV official site station code changed")
    return by

def download(dest):
    r=requests.get(URL,stream=True,timeout=(12,140),
           headers={"User-Agent":"gesui-AT225-original-EEA-2022-Part2-source/20261116"})
    if r.status_code!=200 or not r.url.startswith("https://sdi.eea.europa.eu/"):
        raise RuntimeError("Original EEA source archive official HTTP/source host changed")
    h=hashlib.sha256();n=0
    with dest.open("wb") as f:
        for c in r.iter_content(1024*1024*4):
            if not c:continue
            n+=len(c)
            if n>EXPECTED_SIZE:raise RuntimeError("Original Part2 source file size exceeded")
            f.write(c);h.update(c)
    if n!=EXPECTED_SIZE or h.hexdigest()!=EXPECTED_SHA:
        raise RuntimeError("True original Part2 bytes or SHA256 changed")
    return {"original_publisher":"EEA WQ ICM2026 public Part2",
            "original_source_SHA256":h.hexdigest(),"source_zip_bytes":n,
            "source_public_URL":URL,"source_no_bulk_data_redistribution":True}

def scan(archive,selection):
    found=collections.defaultdict(lambda:{"TN":[],"TP":[]})
    stats=collections.Counter()
    schemes=collections.Counter()
    with zipfile.ZipFile(archive) as z:
        members=[x.filename for x in z.infolist() if x.filename.lower().endswith(".csv")]
        if len(members)!=1:
            raise RuntimeError("Official Part2 original ZIP no unique annual CSV source")
        f=io.TextIOWrapper(z.open(members[0]),encoding="utf-8-sig",newline="",errors="strict")
        head=f.readline();delimiter=";" if head.count(";")>head.count(",") else ","
        names=next(csv.reader([head],delimiter=delimiter))
        required=("countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
                 "phenomenonTimeReferenceYear","parameterWaterBodyCategory",
                 "procedureAnalysedMatrix","observedPropertyDeterminandCode",
                 "resultNumberOfSamples","resultUom","resultMeanValue","resultQualityMeanBelowLOQ")
        if any(k not in names for k in required):
            raise RuntimeError("Original EEA annual source schema changed")
        reader=csv.DictReader(f,fieldnames=names,delimiter=delimiter)
        n=0
        for row in reader:
            n+=1
            if row["countryCode"]!="AT" or row["phenomenonTimeReferenceYear"]!="2022":
                continue
            key=row["monitoringSiteIdentifier"]
            if key not in selection:continue
            stats["original_AT2022_exact_PKINT_any_parameter_rows"]+=1
            scheme=row["monitoringSiteIdentifierScheme"]
            schemes[scheme]+=1
            if scheme!="euMonitoringSiteCode":
                stats["nonstandard_station_scheme_not_pooled"]+=1
                continue
            if row["parameterWaterBodyCategory"]!="RW" or row["procedureAnalysedMatrix"]!="W":
                stats["not_river_water"]+=1;continue
            code=row["observedPropertyDeterminandCode"]
            nu=next((k for k,v in CODES.items() if v==code),None)
            if nu is None:continue
            stats[nu+"_exact_analyte_rows"]+=1
            if row["resultUom"]!=UNITS[nu]:
                stats[nu+"_wrong_unit"]+=1;continue
            try:
                v=float(row["resultMeanValue"]);count=int(row["resultNumberOfSamples"])
            except (ValueError,TypeError):
                stats[nu+"_invalid_statistic"]+=1;continue
            if not math.isfinite(v) or v<0 or count<3:
                stats[nu+"_nonfinite_or_insufficient_original_samples"]+=1;continue
            if row["resultQualityMeanBelowLOQ"].strip().lower() in ("true","1"):
                stats[nu+"_explicit_belowLOQ"]+=1;continue
            found[key][nu].append({"original_year":2022,
                   "original_source_num_samples":count,
                   "source_quality_eligible_numeric_mean_not_saved":True})
            stats[nu+"_qualified_annual_source_records"]+=1
        if n!=EXPECTED_ROWS:raise RuntimeError("Actual EEA Part2 original rowcount differs")
    return found,{"original_full_source_CSV_rows":n,"original_2022_AT_PKINT_scheme_counts":dict(schemes),
                  "original_target_2022_source_stage_counts":dict(stats)}

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    target=sites(a.crosswalk)
    tmp=a.output/"TEMP_original_EEA_PART2_bulk_DO_NOT_REHOST.zip"
    receipt=download(tmp)
    found,source=scan(tmp,target)
    tmp.unlink()
    out=[]
    for pk in sorted(target):
        q=found.get(pk,{"TN":[],"TP":[]})
        TN=len(q["TN"])==1
        TP=len(q["TP"])==1
        out.append({"official_GZUEV_PKINT":pk,
           "original_2022_EEA_reported_NP_outfall_ID":target[pk]["source_2022_EEA_UWWTD_original_plant_ID"],
           "source_annual_2022_AT_PKINT_TN_exact_valid":TN,
           "source_annual_2022_AT_PKINT_TP_exact_valid":TP,
           "source_annual_2022_both_same_PKINT_NP":TN and TP,
           "source_annual_2022_TN_records":len(q["TN"]),
           "source_annual_2022_TP_records":len(q["TP"]),
           "source_annual_2022_values_not_saved":True,
           "source_station_within_2km_not_proven_receiving_reach":True})
    both={q["official_GZUEV_PKINT"] for q in out if q["source_annual_2022_both_same_PKINT_NP"]}
    totals={
      "TN":sum(q["source_annual_2022_AT_PKINT_TN_exact_valid"] for q in out),
      "TP":sum(q["source_annual_2022_AT_PKINT_TP_exact_valid"] for q in out),
      "both":len(both),
      "outfalls_any_both":len({target[pk]["source_2022_EEA_UWWTD_original_plant_ID"] for pk in both})
    }
    summary={
      "schema":SCHEMA,"status":"ORIGINAL_EEA_PART2_SOURCE_ONLY_READ_TO_END_NO_WWTP_ECOLOGY",
      "official_2026_Part2_source":receipt,
      "original_full_CSV_source_support":source,
      "original_225_GZUEV_site_crosswalk_SHA256":CROSSWALK_SHA,
      "original_2022_225_station_sample_annual_numeric_support":totals,
      "original_WISE6_Part1_same_id_2022_rows_already_verified_zero":True,
      "original_national_H2O_WISA_2022_chemistry_unverified":True,
      "station_2km_proximity_not_hydrologic_delivery":True,
      "no_original_water_TNTP_means_exported":True,
      "no_benthic_NDVI_EIV_butterfly_effect_fitted":True}
    file=a.output/"AT225_GZUEV_PKINT_EEA2026_2022_PART2_ANNUAL_NP_SOURCE.csv"
    with file.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    summary["source_annual_225_csv_SHA256"]=sha(file)
    (a.output/"AT225_GZUEV_PKINT_2022_EEA2026_PART2_ANNUAL_SOURCE_GATE.json").write_text(
        json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_AT225_2022_ORIGINAL_PART2_SOURCE="+json.dumps({
       "source_original_rows":source["original_full_source_CSV_rows"],
       "source_genuine_PKINT_target_count":len(target),
       "qualifying_AT2022_annual_site_TN":totals["TN"],
       "qualifying_AT2022_annual_site_TP":totals["TP"],
       "qualifying_AT2022_annual_site_both":totals["both"],
       "outfalls_with_any_valid_2022_annual_both":totals["outfalls_any_both"],
       "scheme_counts":source["original_2022_AT_PKINT_scheme_counts"],
       "source_observation_count":source["original_target_2022_source_stage_counts"],
       "no_biology_or_pollution_effect":True},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--crosswalk",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
