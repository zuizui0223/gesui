#!/usr/bin/env python3
"""Rebuild the immutable FR178 2018/22 same-river aquatic–riparian SOURCE frame.

Only original SHA-locked source geometry and observational month presence are
used; this script never reads TP concentrations, EQR, NDVI or plant responses.
"""
import argparse,csv,hashlib,io,json
from pathlib import Path
from collections import Counter

GEO_SHA="79e91cd2c11776e6a74935bf3eba8fe8d27d924e5196e0b6fab9a0489d45ec7b"
MONTHS_SHA="a02a9f724fb1b1671a9f248ed169620bc25db8b9a5914be566f6728e25b2c4f1"
SOURCE_ID_SHA="76ba45169e58901a2cd82a18244d022773074963a6ac983b527f06302a4aedb8"
FIELDS=["site_key","HYRIV_ID","H7","H6","eaa_lat","eea_lon",
        "original_station_to_original_reach_m","source_reach_length_m",
        "original_along_reach_window_length_m","original_JJA_same_calendar_months",
        "n_original_JJA_same_calendar_months"]
SCHEMA="GESUI_FR178_TRUE_LINE_AND_ORIGINAL_JJA_SOURCE_ONLY_FIXED_V2"

def sha(b):return hashlib.sha256(b).hexdigest()

def run(a):
    original_geom=a.geometry.read_bytes()
    original_month=a.months.read_bytes()
    if sha(original_geom)!=GEO_SHA or sha(original_month)!=MONTHS_SHA:
        raise RuntimeError("Original verified pre-ecology geography/month source SHA changed")
    geo=list(csv.DictReader(io.StringIO(original_geom.decode("utf-8-sig"))))
    m=json.loads(original_month)
    month_rows=m.get("original_263_physically_identical_2018_2022_WQ_and_EEA_source_matched_months")
    if len(geo)!=263 or len(month_rows)!=263:
        raise RuntimeError("Original 263 identical physical river site opportunity rows changed")
    md={x["original_EEA_identical_site_key"]:x for x in month_rows}
    if len(md)!=263 or {x["site_key"] for x in geo}!=set(md):
        raise RuntimeError("Original EEA native IDs do not coincide exactly")
    rows=[]
    audits=Counter()
    for g in geo:
        r=md[g["site_key"]]
        if str(r["original_H6_parent"])!=g["H6"] or str(r["original_HYBAS_ID_L07"])!=g["H7"]:
            raise RuntimeError("Actual original parent basin mismatch")
        actual=sorted(set(r["source_common_original_calendar_months"])&{6,7,8})
        geom=(g["geometry_pass"]=="True" and g["same_original_reach_within_250m"]=="True"
              and g["segment_at_least_900m"]=="True")
        summer=len(actual)>=2
        audits["same_true_reach_and_900m" if geom else "not_source_eligible_geometry"]+=1
        audits["same_years_TP_JJA_2plus" if summer else "source_undersampled_TP_JJA"]+=1
        if not(geom and summer):continue
        row={k:g[k] for k in FIELDS[:9]}
        row["original_JJA_same_calendar_months"]=",".join(map(str,actual))
        row["n_original_JJA_same_calendar_months"]=len(actual)
        rows.append(row)
    rows=sorted(rows,key=lambda r:r["site_key"])
    ids=[x["site_key"] for x in rows]
    ids_sha=sha("\n".join(ids).encode())
    h7={x["H7"] for x in rows}
    h6={x["H6"] for x in rows}
    if len(rows)!=178 or len(h7)!=58 or len(h6)!=23 or ids_sha!=SOURCE_ID_SHA:
        raise RuntimeError("Fixed source-only 178/58H7/23H6 station frame changed")
    if not (len(rows)>=120 and len(h7)>=40 and len(h6)>=15):
        raise RuntimeError("Predefined source geography/time threshold failed")
    a.out.mkdir(parents=True,exist_ok=True)
    path=a.out/"FR178_FIXED_NO_SATELLITE_OUTCOME_STATIONS.csv"
    with path.open("w",newline="",encoding="utf-8") as f:
        wr=csv.DictWriter(f,fieldnames=FIELDS);wr.writeheader();wr.writerows(rows)
    out={"schema":SCHEMA,"status":"PASS_178_PRE_SATELLITE_WATER_RIVER_SOURCE_ONLY",
         "original_263":263,"true_reach_geo_n":192,"original_JJA_TP_dates_2plus_n":241,
         "fixed_join_n":len(rows),"true_H7_n":len(h7),"true_H6_n":len(h6),
         "selected_native_site_IDs_SHA256":ids_sha,
         "source_original_geo_csv_SHA256":GEO_SHA,"source_original_month_only_json_SHA256":MONTHS_SHA,
         "selected_roster_csv_SHA256":sha(path.read_bytes()),
         "source_exclusions":dict(audits),
         "original_river_geometry_from_SHA_pinned_thirdparty_conversion_not_official_ZIP":True,
         "no_new_true_2018_2022_S2_pixel_coverage_claim":True,
         "no_original_TP_level_EQR_or_NDVI_read_to_select_station":True}
    (a.out/"FR178_FIXED_PRE_SATELLITE_SOURCE_QA.json").write_text(
        json.dumps(out,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("FR178_RIVER_JJA_SOURCE="+json.dumps(out,sort_keys=True),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--geometry",type=Path,required=True)
    p.add_argument("--months",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())

if __name__=="__main__":main()
