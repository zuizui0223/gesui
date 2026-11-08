#!/usr/bin/env python3
"""Source-only EEA UWWTD historical 2018 measured WWTP N/P mass audit.

Uses exactly the public 2017–18 v07/r01 ZIP, full archived reporter+plant IDs,
original measured (not calculated or estimated) N/P discharge fields, original
discharge-point/receiving-water code and reported-year fields. No EQR outcomes,
no choice of historical country based on biology. Archive itself is not mirrored.
"""
from __future__ import annotations
import argparse,hashlib,io,json,re,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd

ROOT="eea_t_uwwtd-waterbase_p_2017-2018_v07_r01/Waterbase_UWWTD_v8_csv/"
TABLES=("UWWTPS.csv","DischargePoints.csv","ReportPeriod.csv","LOV.csv")
META_RUN=42440990345
EXPECTED_MEMBERS=23
BBOX=(-12.,44.,26.,64.)
MEASURED=("uwwNDischargeMeasured","uwwPDischargeMeasured")
FULL_JOIN=("repCode","uwwCode")
POINT_FIELDS=("dcpLatitude","dcpLongitude")
def sha_file(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
    return h.hexdigest()

def check_archive(a,inventory):
    saved=json.loads(inventory.read_text())
    assert saved["source_version"]=="eea_t_uwwtd-waterbase_p_2017-2018_v07_r01"
    meta={x["relative_path"]:(int(x["bytes_uncompressed"]),str(x["source_member_CRC32"]))
           for x in saved["all_public_original_archive_members"]}
    if len(meta)!=EXPECTED_MEMBERS:raise ValueError("Original file count changed")
    output={}
    with zipfile.ZipFile(a) as z:
        current={x.filename:(int(x.file_size),format(x.CRC,"08x")) for x in z.infolist() if not x.is_dir()}
        if current!=meta:raise ValueError("EEA historical ZIP member size/name/CRC changed")
        for name in TABLES:
            fname=ROOT+name
            raw=z.read(fname)
            h=hashlib.sha256(raw).hexdigest()
            output[name]={"name":fname,"source_uncompressed_SHA256":h,
                         "bytes":len(raw),"member_CRC32":meta[fname][1]}
            output[name]["frame"]=pd.read_csv(io.BytesIO(raw),dtype="string",low_memory=False)
    return output

def freq(df,col,maxn=30):
    return {str(k):int(v) for k,v in df[col].fillna("MISSING").value_counts().head(maxn).items()}

def add_key(d):
    for col in FULL_JOIN:
        d[col]=d[col].astype("string").str.strip()
    d["source_key"]=d.repCode.fillna("")+"|"+d.uwwCode.fillna("")
    d.loc[d.repCode.isna()|d.uwwCode.isna()|d.repCode.eq("")|d.uwwCode.eq(""),"source_key"]=pd.NA
    return d

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--inventory",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    files=check_archive(args.archive,args.inventory)
    raw={n:f.pop("frame") for n,f in files.items()}
    uw=add_key(raw["UWWTPS.csv"])
    outlets=add_key(raw["DischargePoints.csv"])
    periods=raw["ReportPeriod.csv"].copy()
    lov=raw["LOV.csv"].copy()
    if set(MEASURED)-set(uw.columns) or set(POINT_FIELDS)-set(outlets.columns):
        raise ValueError("2018 measured load and outlet coordinates unavailable")
    for col in MEASURED:
        uw[col]=pd.to_numeric(uw[col],errors="coerce")
    for col in POINT_FIELDS:
        outlets[col]=pd.to_numeric(outlets[col],errors="coerce")
    uw["both_NP_measured_positive"]=uw[list(MEASURED)].gt(0).all(axis=1)
    uw["N_measured_positive"]=uw[MEASURED[0]].gt(0)
    uw["P_measured_positive"]=uw[MEASURED[1]].gt(0)
    whole=len(uw)
    sourcecodes=uw.source_key.value_counts()
    distinct_ww=uw[uw.source_key.notna() & uw.source_key.map(sourcecodes).eq(1)].copy()
    outs=outlets[outlets.source_key.notna()].copy()
    outlet_counts=outs.source_key.value_counts()
    distinct_ww["n_source_outlets"]=distinct_ww.source_key.map(outlet_counts).fillna(0).astype(int)
    both=distinct_ww[distinct_ww.both_NP_measured_positive].copy()
    if len(both)<10:raise ValueError("Historical reported BOTH measured N/P sources too sparse")
    joined=both[both.n_source_outlets.eq(1)].merge(
        outs.drop_duplicates("source_key"),on="source_key",
        how="inner",suffixes=("_plant","_dcp"),validate="one_to_one")
    joined["outlet_coords_in_EU_bbox"]=(
        joined.dcpLatitude.between(BBOX[1],BBOX[3])&
        joined.dcpLongitude.between(BBOX[0],BBOX[2]))
    europe=joined.loc[joined.outlet_coords_in_EU_bbox].copy()
    if "dcpWaterBodyType" in europe.columns:
        e=europe.dcpWaterBodyType.astype("string").str.strip().str.upper()
        europe["is_direct_freshwater_FW"]=e.eq("FW")
    else:
        europe["is_direct_freshwater_FW"]=False
    source_qc={
      "original_source_version":"eea_t_uwwtd-waterbase_p_2017-2018_v07_r01",
      "original_members_matched_initial_report":True,
      "original_2018_zip_outer_sha256":sha_file(args.archive),
      "source_table_original_SHA256":files,
      "original_plant_rows":whole,
      "original_unique_plant_key_rows":len(distinct_ww),
      "original_report_period_rows":len(periods),
      "original_report_period_2018_reference_flags":{
         col:freq(periods,col,35) for col in ("repSituationAt","repReportedPeriod","repCode","rptMStateKey") if col in periods},
      "plant_national_report_codes":freq(uw,"repCode",50),
      "source_plant_state_codes":freq(uw,"uwwState",35),
      "outlet_waterbody_codes":freq(outlets,"dcpWaterBodyType",40)
               if "dcpWaterBodyType" in outlets else {},
      "outlet_surface_water_codes":freq(outlets,"dcpSurfaceWaters",40)
               if "dcpSurfaceWaters" in outlets else {},
      "source_reference_Codelist_LOV_sample":[dict(x) for _,x in lov.head(5).iterrows()],
      "source_2018_N_measured_positive_plant_rows":int(uw.N_measured_positive.sum()),
      "source_2018_P_measured_positive_plant_rows":int(uw.P_measured_positive.sum()),
      "source_2018_NP_both_measured_positive_plant_rows":int(uw.both_NP_measured_positive.sum()),
      "source_2018_unique_plant_NP_both_measured_positive_rows":len(both),
      "source_2018_positive_both_measured_plus_unique_outlet":len(joined),
      "source_2018_positive_both_measured_unique_outlet_in_EU_BBOX":len(europe),
      "source_2018_positive_both_measured_unique_outlet_direct_FW_in_BBOX":
          int(europe.is_direct_freshwater_FW.fillna(False).sum()),
      "source_2018_positive_both_measured_unique_outlet_waterbody_code":{
          str(k):int(v) for k,v in europe.dcpWaterBodyType.fillna("MISSING").value_counts().items()},
      "source_2018_positive_both_unique_outlet_by_reporter":{
          str(k):int(v) for k,v in europe.repCode.value_counts().items()},
      "source_2018_missing_or_multiple_outlet":{
         "no_outlet":int(both.n_source_outlets.eq(0).sum()),
         "multiple_outlets":int(both.n_source_outlets.gt(1).sum())},
      "N_measured_mass_original_positive_unit_metadata_not_assumed":True,
      "source_loaded_without_Biology_EQR_outcomes":True,
    }
    for col in MEASURED:
        pos=europe.loc[europe.is_direct_freshwater_FW.fillna(False),col].astype(float)
        source_qc["source_"+col+"_positive_mass_quantiles"]=(
          {str(q):float(pos.quantile(q)) for q in (0.1,0.5,0.9)} if len(pos) else {})
    report={
       "schema":"gesui_UWWTD_2018_original_measured_NP_UWWTP_outlet_source_only_QC_v01",
       "date_utc":datetime.now(timezone.utc).isoformat(),
       "first_source_2018_archive_run":META_RUN,
       "source_QC":source_qc,
       "source_mass_values_are_discharge_reported_not_receptor_chemistry":True,
       "original_2018_Biology_EQR_not_opened":True,
       "source_history_and_reporting_reference_year_2018_needs_exact_validation":True,
       "unmeasured_calculated_or_estimated_source_loads_not_substituted":True,
       "limits":[
        "The 2017–18 archive includes published original national reporting, but each reporter's reference year must be verified",
        "The old outlet uses dcpWaterBodyType and dcpSurfaceWaters; direct freshwater must match official code meaning, not nearest river alone",
        "A positive mass field may be reported through noncurrent source or nonoperational plant, must inspect original state",
        "Original 2018-specific units and measurement period must be checked from code dictionary",
        "All-source WWTP annual masses are not river P or N causal doses"
       ]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_2018_UWWTD_MEASURED_NP_SOURCE_GATE",json.dumps({
      k:v for k,v in source_qc.items()
      if k not in ("source_table_original_SHA256","source_reference_Codelist_LOV_sample",
                   "original_report_period_2018_reference_flags","plant_national_report_codes",
                   "source_plant_state_codes","outlet_waterbody_codes","outlet_surface_water_codes")},
      sort_keys=True),flush=True)
    print("GESUI_2018_UWWTD_YEAR_STATE_AND_OUTLET_CODES",json.dumps({
        k:v for k,v in source_qc.items() if k in
        ("original_report_period_2018_reference_flags","source_plant_state_codes",
         "outlet_waterbody_codes","outlet_surface_water_codes")},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
