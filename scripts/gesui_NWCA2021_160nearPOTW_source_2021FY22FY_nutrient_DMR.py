#!/usr/bin/env python3
"""Actual EPA FY2021/FY2022 DMR source-only scan for fixed 160 NPDES permits.

Only 2021 calendar-month DMR data. Two FY releases required because FY2021
ends Sep2021 and FY2022 contains Oct-Dec2021. Uses site/permit hashes frozen
BEFORE reading original DMRs; no wetland plant Shannon or TN/TP values read.
Permit-level DMR DOES NOT match the geodesically nearby outfall feature yet.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

BASE=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=BASE/"EPA_NWCA2021_160_SEWAGE_PERMITS_BEFORE_DMR2021_SOURCE.json"
CODE={"00600":"total nitrogen as N","00665":"total phosphorus as P","50050":"effluent flow"}
EFFLUENT={"1","2","EG"}
EARLY_DATE="2021-01-01"
LATE_DATE="2021-12-31"
EXPECTED_FY=(2021,2022)
BASE_URL="https://echo.epa.gov/files/echodownloads/npdes_dmrs_fy"
FIELDS=("EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","PERM_FEATURE_TYPE_CODE",
        "PARAMETER_CODE","MONITORING_LOCATION_CODE","MONITORING_PERIOD_END_DATE",
        "DMR_VALUE_NMBR","DMR_VALUE_STANDARD_UNITS","DMR_UNIT_DESC","STANDARD_UNIT_DESC",
        "DMR_VALUE_QUALIFIER_CODE","NODI_CODE","VALUE_TYPE_CODE",
        "DMR_VALUE_ID","VERSION_NMBR","STATISTICAL_BASE_CODE")
def sha_file(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()
def sha_text(x):return hashlib.sha256(x.encode()).hexdigest()
def norm(x):return str(x or "").strip().upper()
def calendaryear(x):
    s=norm(x)
    if len(s)>=4 and s[:4].isdigit() and s[:4]=="2021":return "2021"
    if s.endswith("2021"):return "2021"
    if s.endswith("/21") or s.endswith("-21"):
        # last two-digit year ambiguous; preserve source and do not guess
        return "UNKNOWN_TWO_DIGIT_YEAR"
    return "OTHER_OR_UNKNOWN"
def scan(root,source_file,fy):
    source=json.loads((root/SOURCE).read_text(encoding="utf-8"))
    if source["schema"]!="gesui_US_NWCA2021_original128near_160municipalNPDESpermitted_2021DMR_PREOUTCOME_source_v01":
        raise ValueError("Original pre-DMR source not frozen")
    if source["NO_FY2021_FY2022_DMR_measurements_accessed_in_this_source_stage"] is not True:
        raise ValueError("Pre-outcome source not blind to DMR values")
    permitset=set(source["exact_source_frozen_near10_registered_municipal_permit_ID_SHA256"])
    if len(permitset)!=160 or len(source["original_site_near10_permits_2021outfallPlausible_from_2026_registry"])!=128:
        raise ValueError("Original near160 and 128 source sites changed")
    expected_member=f"npdes_dmr_fy{fy}.csv"
    if fy not in EXPECTED_FY:raise ValueError("Only fixed FY2021 and FY2022 source permitted")
    if source_file.stat().st_size<50_000_000:
        raise ValueError("EPA original national FY DMR ZIP unexpectedly small")
    seen=Counter()
    selected=[]
    with zipfile.ZipFile(source_file) as z:
        names=[name for name in z.namelist() if name.split("/")[-1].lower()==expected_member]
        if len(names)!=1:
            raise ValueError("Original national EPA DMR CSV ZIP member not found: "+str(z.namelist()[:12]))
        with z.open(names[0]) as binary:
            text=io.TextIOWrapper(binary,encoding="utf-8-sig",errors="replace",newline="")
            reader=csv.reader(text)
            head=next(reader)
            idx={name:i for i,name in enumerate(head)}
            missing=set(FIELDS)-set(idx)
            if missing:raise ValueError("Original DMR national source dictionary mismatch: "+repr(sorted(missing)))
            get=lambda r,name: r[idx[name]].strip() if len(r)>idx[name] else ""
            for r in reader:
                seen["all_national_DMR_table_rows"]+=1
                if len(r)<len(head):
                    seen["malformed_short_row"]+=1
                    continue
                code=get(r,"PARAMETER_CODE")
                if code not in CODE:continue
                seen["fixed_N_or_P_or_flow_parameter_rows"]+=1
                pid=norm(get(r,"EXTERNAL_PERMIT_NMBR"))
                permit_hash=sha_text(pid)
                if permit_hash not in permitset:continue
                seen["source_160_permits_fixed_parameter_rows_all_years"]+=1
                date=get(r,"MONITORING_PERIOD_END_DATE")
                yr=calendaryear(date)
                if yr!="2021":
                    seen["source_160_permits_rows_other_or_unrecognized_year"]+=1
                    if yr=="UNKNOWN_TWO_DIGIT_YEAR":seen["ambiguous_date_format_rows"]+=1
                    continue
                seen["2021_source_160permit_parameter_rows"]+=1
                feature=get(r,"PERM_FEATURE_NMBR")
                featuretype=get(r,"PERM_FEATURE_TYPE_CODE")
                monitoring=get(r,"MONITORING_LOCATION_CODE")
                nodi=norm(get(r,"NODI_CODE"))
                qual=norm(get(r,"DMR_VALUE_QUALIFIER_CODE"))
                vt=norm(get(r,"VALUE_TYPE_CODE"))
                x=get(r,"DMR_VALUE_NMBR")
                try:
                    val=float(x)
                    numeric=(math.isfinite(val) and val>=0)
                except (ValueError,OverflowError):numeric=False
                main_effluent=(norm(featuretype)=="EXO" and norm(monitoring) in EFFLUENT)
                usable=(main_effluent and numeric and not nodi and qual in ("","="))
                if main_effluent:seen["2021_original_160permits_extoutfall_effluent_rows"]+=1
                if usable:seen["2021_numeric_uncensored_extoutfall_final_effluent_rows"]+=1
                record={
                    "source_permit_ID_SHA256":permit_hash,
                    "source_permit_feature_ID_SHA256":sha_text(pid+"|"+norm(feature)),
                    "source_permitted_feature_suffix_hash_SHA256":sha_text(norm(feature)),
                    "original_parameter_code":code,
                    "original_monitoring_period_end_date":date,
                    "original_monitoring_location_code":norm(monitoring),
                    "original_permitted_feature_type_code":norm(featuretype),
                    "original_2021_DMR_VALUE_ID_SHA256":sha_text(norm(get(r,"DMR_VALUE_ID"))+"|"+pid+"|"+norm(feature)),
                    "original_DMR_VERSION_NMBR":get(r,"VERSION_NMBR"),
                    "original_DMR_STATISTICAL_BASE_CODE":get(r,"STATISTICAL_BASE_CODE"),
                    "original_DMR_VALUE_TYPE_CODE":vt,
                    "original_DMR_UNIT_DESC":get(r,"DMR_UNIT_DESC"),
                    "original_STANDARD_UNIT_DESC":get(r,"STANDARD_UNIT_DESC"),
                    "original_value_numeric_nonnegative":numeric,
                    "original_DMR_NODI_CODE":nodi or None,
                    "original_DMR_VALUE_QUALIFIER_CODE":qual or None,
                    "strict_external_outfall_effluent_numeric_uncensored":usable,
                    "DMR_original_reported_numeric_value_NOT_REPUBLISHED":True,
                }
                selected.append(record)
                seen[f"2021_{code}_all_monitoring_rows"]+=1
                if usable:seen[f"2021_{code}_numeric_final_effluent_rows"]+=1
                if nodi:seen["2021_has_NoData_NODI"]+=1
                if qual not in ("","="):seen["2021_has_qualifier_not_equal"]+=1
                if not numeric:seen["2021_DMR_value_missing_non_numeric_invalid"]+=1
    if seen["all_national_DMR_table_rows"]<100000:
        raise ValueError("Original nationwide FY2021/2022 DMR table source unexpectedly truncated")
    group=defaultdict(lambda:{"N":set(),"P":set(),"flow":set(),"any":set(),"features":set()})
    for row in selected:
        if not row["strict_external_outfall_effluent_numeric_uncensored"]:continue
        k=row["source_permit_ID_SHA256"]
        t=CODE[row["original_parameter_code"]]
        ty={"total nitrogen as N":"N","total phosphorus as P":"P","effluent flow":"flow"}[t]
        group[k][ty].add(row["original_monitoring_period_end_date"])
        group[k]["any"].add(row["original_monitoring_period_end_date"])
        group[k]["features"].add(row["source_permit_feature_ID_SHA256"])
    output={
        "schema":"gesui_NWCA2021_160_original_municipal_permits_EPA_national_DMR_fiscal_sourceonly_v01",
        "source_fiscal_year":fy,
        "source_original_FY_DMR_official_url":f"{BASE_URL}{fy}.zip",
        "source_original_FY_DMR_ZIP_SHA256":sha_file(source_file),
        "source_original_160permit_preDMR_manifest_SHA256":sha_file(root/SOURCE),
        "original_CY2021_calendar_end_date_filter_not_equating_fiscal_year_to_calendar":True,
        "original_2021_total_N_code_00600_and_total_P_code_00665_only":True,
        "original_2021_external_outfall_effluent_location_codes":["1","2","EG"],
        "source_rows_national_and_160_permits":dict(seen),
        "source_160_permit_level_any_fy_DMR_numeric_NP_reporting_summary":{
            "permits_any_numeric_final_effluent_P":sum(bool(v["P"]) for v in group.values()),
            "permits_any_numeric_final_effluent_N":sum(bool(v["N"]) for v in group.values()),
            "permits_both_N_and_P_numeric_final_effluent":sum(bool(v["N"]) and bool(v["P"]) for v in group.values()),
            "permits_any_numeric_final_effluent_flow":sum(bool(v["flow"]) for v in group.values()),
        },
        "original_2021_full_source_selected_DMR_rows_no_raw_reported_values":selected,
        "official_EPA_FY_DMR_extract_not_original_2021_wwtp_to_wetland_flow":True,
        "strict_2021_outfall_feature_match_to_geospatial_2026_ECHO_not_yet_verified":True,
        "source_no_plant_TNTP_or_Shannon_observations_read":True,
        "UTC":datetime.now(timezone.utc).isoformat()
    }
    return output

def main():
    p=argparse.ArgumentParser();p.add_argument("--fy",type=int,required=True)
    p.add_argument("--zip",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=scan(root,a.zip,a.fy)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"fy":a.fy,
        "source_rows":result["source_rows_national_and_160_permits"],
        "source_permits":result["source_160_permit_level_any_fy_DMR_numeric_NP_reporting_summary"]},indent=2))
if __name__=="__main__":main()
