#!/usr/bin/env python3
"""Recover true original 2018 and 2022 calendar-month measured summer TP for
the source-fixed one-per-24-H6 French water stations; no satellite pixels read.

The 24 station source sample was chosen by SHA256 site IDs before this new
source test; do not exclude stations by TP magnitude or their satellite NDVI.
"""
import argparse,csv,datetime,hashlib,io,json,math,statistics,urllib.request,zipfile
from pathlib import Path
from collections import defaultdict,Counter

ORIGINAL_STAC_SHA="22414dec230bc9bf9c308d823f166303cc843143c6ade922af6cb44b1f2c3f3c"
ORIGINAL_DATE_SHA="2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264"
ORIGINAL_TPZIP_SHA="26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e"
ORIGINAL_TPZIP_MD5="b596dbe9833c49426147a99f10f3fe31"
ORIGINAL_TPZIP_N_BYTES=25003722
ORIGINAL_STATION_CSV_SHA="1824b0d9b84d35a460b4d1d0513659f93a16c59d7f5666a2d8707be1c6da9942"
SUMMER=(6,7,8)
YEARS=(2018,2022)
GATE_MIN_SITE=18
SCHEMA="GESUI_FR263_24H6_ORIGINAL_2018_22_JJA_MONTH_MATCHED_TP_SOURCE_AND_VALUES_V1"
URL="https://ndownloader.figshare.com/files/70896288"

def sha(b):return hashlib.sha256(b).hexdigest()
def verify(path,expected):
    b=path.read_bytes()
    if sha(b)!=expected:raise RuntimeError("Original source artifact SHA mismatch: "+str(path))
    return b

def source_sites(signed_stac,signed_dates):
    stac=json.loads(verify(signed_stac,ORIGINAL_STAC_SHA))
    dates=json.loads(verify(signed_dates,ORIGINAL_DATE_SHA))
    if (stac["schema"]!="GESUI_20261117_FR263_S2_2018_2022_24H6_PILOT_SOURCE_OPPORTUNITY_V1" or
       stac["source_supported_bothyear_pilot_n"]!=24 or len(stac["pilot_sites"])!=24):
        raise RuntimeError("Original signed 24-H6 station source changed")
    ds=dates["original_source_only_exact_same_EEA2018_2022_biological_site_and_annual_measured_TP_support"]
    if len(ds)!=263 or len({z["EEA2022_site_key"] for z in ds})!=263:
        raise RuntimeError("Original 263 EEA+TP source changed")
    keys={z["EEA2022_site_key"]:z for z in ds}
    out=[]
    for x in stac["pilot_sites"]:
        key=x["site_key"]
        if key not in keys:raise RuntimeError("Selected 24 H6 station not in original 263")
        y=keys[key]
        if x["H6"]!=y["HYBAS_ID_L06_parent"] or x["H7"]!=y["HYBAS_ID_L07"] or not(y["source_2018_eligible"] and y["source_2022_eligible"]):
            raise RuntimeError("Original fixed station river/year H6 geography mismatch")
        out.append({"site_key":key,"H6":x["H6"],"H7":x["H7"],"WQ":y["original_WQ_ID"]})
    if len(out)!=24 or len({x["H6"] for x in out})!=24:
        raise RuntimeError("One unique preoriginal station per H6 missing")
    return out

def load_public_zip():
    req=urllib.request.Request(URL,headers={"User-Agent":"gesui-20261117-FR24H6-source-only-original-JJA-Tp-matching","Accept":"application/octet-stream,application/zip,*/*"})
    bb=bytearray()
    with urllib.request.urlopen(req,timeout=120) as resp:
        while chunk:=resp.read(131072):
            bb.extend(chunk)
            if len(bb)>26000000:raise RuntimeError("Original publisher source too large; fail-closed")
    raw=bytes(bb)
    if (len(raw)!=ORIGINAL_TPZIP_N_BYTES or sha(raw)!=ORIGINAL_TPZIP_SHA or
       hashlib.md5(raw).hexdigest()!=ORIGINAL_TPZIP_MD5):
        raise RuntimeError("Original Figshare TP.2 ZIP bytes/MD5/SHA changed")
    return zipfile.ZipFile(io.BytesIO(raw))

def source_and_monthly_true_TP(archive,site):
    path="TP.2_-_Harmonised/"+site["WQ"]+".csv"
    if path not in archive.namelist():
        raise RuntimeError("Original selected WQ source file absent: "+site["WQ"])
    records=csv.DictReader(io.StringIO(archive.read(path).decode("utf-8-sig")))
    dayvals=defaultdict(list)
    annual=defaultdict(list)
    for v in records:
        try:
            day=datetime.date.fromisoformat(str(v["date"]).strip()[:10])
            f=float(v["value"])
        except (KeyError,TypeError,ValueError):continue
        if day.year not in YEARS or not math.isfinite(f) or f<=0:continue
        annual[day.year].append(f)
        if day.month in SUMMER:dayvals[(day.year,day.month)].append(f)
    # Source MONTH availability independent of magnitude signs (already positive flagged)
    common=sorted(m for m in SUMMER if dayvals[(2018,m)] and dayvals[(2022,m)])
    q={**site,"original_source_common_summer_months":common,
        "n_original_2018_2022_common_summer_months":len(common),
        "original_summer_source_eligible":len(common)>=2,
        "original_2018_2022_annual_TP_positive_count":{str(y):len(annual[y]) for y in YEARS}}
    rows=[]
    for m in common:
        a=float(statistics.median(dayvals[(2018,m)]))
        b=float(statistics.median(dayvals[(2022,m)]))
        rows.append({"month":m,"2018_real_TP_mg_L":a,"2022_real_TP_mg_L":b,
                     "delta_ln_TP_2022_minus_2018":math.log(b)-math.log(a),
                     "original_2018_TP_obs_n":len(dayvals[(2018,m)]),
                     "original_2022_TP_obs_n":len(dayvals[(2022,m)])})
    q["original_summer_real_TP_month_pair_rows"]=rows
    if len(common)>=2:
        q["2018_summer_common_month_mean_ln_TP"]=float(statistics.mean(math.log(y["2018_real_TP_mg_L"]) for y in rows))
        q["2022_summer_common_month_mean_ln_TP"]=float(statistics.mean(math.log(y["2022_real_TP_mg_L"]) for y in rows))
        q["delta_summer_common_month_mean_ln_TP"]=float(
            q["2022_summer_common_month_mean_ln_TP"]-q["2018_summer_common_month_mean_ln_TP"])
    q["original_2018_annual_TP_median_mg_L"]=float(statistics.median(annual[2018]))
    q["original_2022_annual_TP_median_mg_L"]=float(statistics.median(annual[2022]))
    return q

def verify_original_numeric(station_csv,rows):
    b=verify(station_csv,ORIGINAL_STATION_CSV_SHA)
    r=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"))))
    d={x["EEA_identical_native_site_key_2018_2022"]:x for x in r}
    if len(d)!=263:raise RuntimeError("Original aquatic station 263 source rows drift")
    for x in rows:
        z=d[x["site_key"]]
        for year in YEARS:
            actual=float(z[f"real_original_{year}_measured_TP_median_mg_L"])
            recomputed=x[f"original_{year}_annual_TP_median_mg_L"]
            if not math.isclose(actual,recomputed,abs_tol=1e-11,rel_tol=1e-11):
                raise RuntimeError(f"Original real TP median source disagrees at {x['site_key']} year {year}")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--stac",required=True,type=Path)
    p.add_argument("--datesource",required=True,type=Path)
    p.add_argument("--original-number",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    sites=source_sites(a.stac,a.datesource)
    archive=load_public_zip()
    rows=[source_and_monthly_true_TP(archive,s) for s in sites]
    verify_original_numeric(a.original_number,rows)
    passing=[x for x in rows if x["original_summer_source_eligible"]]
    counts=Counter(x["n_original_2018_2022_common_summer_months"] for x in rows)
    result={
       "schema":SCHEMA,
       "status":"PASS_FIXED_24_H6_JJA_SAME_MONTH_RIVER_TP_SOURCE_NO_NDVI_FIT" if len(passing)>=GATE_MIN_SITE else "STOP_SOURCE_REAL_SUMMER_TP_SUPPORT_NOT_ECOLOGY_NULL",
       "original_24_H6_source_station_n":len(sites),
       "original_real_JJA_2018_2022_monthmatched_TP_eligible":len(passing),
       "predeclared_min_eligible_source_18_H6":GATE_MIN_SITE,
       "original_source_common_summer_month_counts":dict(counts),
       "original_TP2_Figshare_checksum_SHA256":ORIGINAL_TPZIP_SHA,
       "original_TP2_Figshare_publisher_MD5":ORIGINAL_TPZIP_MD5,
       "all_original_24_annual_TP_recomputed_medians_equal_prior_genuine_source":True,
       "original_source_24_station_records":rows,
       "same_original_S2_NDVI_month_coverage_NOT_OPENED":True,
       "no_source_2021_landover_2018_2022_assumption_necessary_to_recover_TP":True,
       "not_exact_same_day_TP_Sentinel2_or_water_flux":True,
       "no_original_wastewater_specific_or_terrestrial_ecology_model_fit":True
    }
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"FR263_24H6_ORIGINAL_2018_2022_SUMMER_MONTH_TP_SOURCE.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_24H6_SUMMER_ORIGINAL_TP_SOURCE="+json.dumps({
      "n_original":24,"n_eligible":len(passing),"counts_by_JJA_shared_months":dict(counts),
      "status":result["status"]}),flush=True)
if __name__=="__main__":
    main()
