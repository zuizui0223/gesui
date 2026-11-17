#!/usr/bin/env python3
"""Reproduce independent, source-only FR178 genuine riverline and same JJA water-month
frame before opening ANY new 178-site Sentinel2 vegetation pixels.

Original 263 genuine WQ/EEA sites 2018/22 + original signed hydrological line
distance/length source + authentic positive TP observation-month support.
Only ID, geographic metadata and calendar observation opportunity are used.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter

GEO_SHA="79e91cd2c11776e6a74935bf3eba8fe8d27d924e5196e0b6fab9a0489d45ec7b"
WQMONTH_SHA="3d31c4c0f6f0dd251919b8c1b9eda18729cd4feebe8a713349940a7683c91b67"
SCHEMA="GESUI_FR178_ORIGINAL_HYDRORIVERS_GEOMETRY_AND_JJA_MONTHS_SOURCE_FIRST_V1"
SUMMER={6,7,8}
def sha(raw):return hashlib.sha256(raw).hexdigest()
def records(path,expected):
    raw=path.read_bytes()
    if sha(raw)!=expected:raise RuntimeError("Exact original signed source bytes changed: "+str(path))
    return list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))

def run(args):
    geo=records(args.geo,GEO_SHA)
    obs=records(args.months,WQMONTH_SHA)
    if len(geo)!=263 or len(obs)!=263:raise RuntimeError("Source population no longer original same 263 physical sites")
    keyed={r["EEA_identical_station_key"]:r for r in obs}
    if len(keyed)!=263 or len({r["site_key"] for r in geo})!=263:raise RuntimeError("Duplicated original EEA site key")
    matched=[];fails=[]
    for p in geo:
        key=p["site_key"]
        if key not in keyed:raise RuntimeError("Original same source site key missing")
        q=keyed[key]
        if p["H7"]!=q["HYBAS_ID_L07"] or p["H6"]!=q["HYBAS_ID_L06"]:
            raise RuntimeError("Hydrographical parent H6/H7 source identity disagrees")
        # Intentionally DO NOT use numerical TP/I2M2 EQR outcome columns.
        allmonths={int(v) for v in q["source_original_calendar_months_shared"].split(",") if v}
        common=sorted(SUMMER & allmonths)
        geo_ok=p["geometry_pass"].strip().lower()=="true"
        eligible=geo_ok and len(common)>=2
        if not eligible:
            fails.append({"site_key":key,
                "reason":"bad_original_riverline_geometry" if not geo_ok else "JJA_both_year_TP_months_lt2",
                "bothyear_source_JJA_months":common})
            continue
        matched.append({"site_key":key,"HYRIV_ID":p["HYRIV_ID"],
          "H7":p["H7"],"H6":p["H6"],
          "eea_lat":p["eaa_lat"],"eea_lon":p["eea_lon"],
          "original_station_to_original_reach_m":p["original_station_to_original_reach_m"],
          "source_reach_length_m":p["source_reach_length_m"],
          "original_along_reach_window_length_m":p["original_along_reach_window_length_m"],
          "original_JJA_same_calendar_months":",".join(map(str,common)),
          "n_original_JJA_same_calendar_months":str(len(common))})
    h7={r["H7"] for r in matched};h6={r["H6"] for r in matched}
    if not(len(matched)==178 and len(h7)==58 and len(h6)==23):
        raise RuntimeError("Independent source frame changed unexpectedly: "+str((len(matched),len(h7),len(h6))))
    if len({r["site_key"] for r in matched})!=178:
        raise RuntimeError("Original selected site IDs not unique")
    matched.sort(key=lambda r:r["site_key"])
    args.out.mkdir(parents=True,exist_ok=True)
    outcsv=args.out/"FR178_TRUE_REACH_TWO_YEAR_JJA_MEASURED_TP_SOURCE_ONLY.csv"
    with outcsv.open("w",encoding="utf-8",newline="") as stream:
        w=csv.DictWriter(stream,fieldnames=list(matched[0]))
        w.writeheader();w.writerows(matched)
    report={
      "schema":SCHEMA,
      "status":"PASS_ORIGINAL_178_SOURCE_ONLY_NO_SATELLITE_PIXEL_ECOLOGY",
      "original_full_same_2018_22_WQ_EEA_sites":len(geo),
      "original_genuine_riverline_geometry_pass":sum(r["geometry_pass"].lower()=="true" for r in geo),
      "original_bothyear_JJA_TP_months_at_least2":sum(len(SUMMER&{int(v) for v in r["source_original_calendar_months_shared"].split(",") if v})>=2 for r in obs),
      "original_both_178_sites":len(matched),
      "original_178_unique_H7":len(h7),
      "original_178_unique_H6":len(h6),
      "predefined_source_minima_station_H7_H6":[120,40,15],
      "source_go_for_independent_raster_opportunity":bool(len(matched)>=120 and len(h7)>=40 and len(h6)>=15),
      "original_selection_uses_only_true_riverline_geometry_and_TP_sampling_months":True,
      "old_24H6_joint_primary_gate_n14_below18_remains_STOP":True,
      "this_expanded_posthoc_sample_is_NOT_a_confirmatory_independent_new_cohort":True,
      "FR178_source_csv_SHA256":sha(outcsv.read_bytes()),
      "original_input_geography_SHA256":GEO_SHA,
      "original_input_source_months_SHA256":WQMONTH_SHA,
      "no_2018_2022_NDVI_pixel_values_or_TP_magnitudes_selected":True,
      "not_wwtp_NP_flux_or_ecosystem_service":True
    }
    (args.out/"FR178_GEOMETRY_JJA_TP_SOURCE_GATE.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print("GESUI_FR178_TRUE_RIVER_JJA_SOURCE_FIRST="+json.dumps(report,sort_keys=True),flush=True)
    return report

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--geo",type=Path,required=True)
    p.add_argument("--months",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
