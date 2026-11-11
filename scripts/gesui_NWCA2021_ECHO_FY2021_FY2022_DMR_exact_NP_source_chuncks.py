#!/usr/bin/env python3
"""Only REPORTING SUPPORT from EPA's ORIGINAL FY2021/FY2022 national DMR CSV.

Source pre-frozen exact 344 municipal permitted feature numbers from registered
2026 ECHO POTW points near 128 EPA NWCA2021 wetlands. Filter official DMR
permit AND feature AND EXO (actual external outfall) to original 2021 *calendar*
monitoring end dates. No NWCA plant Shannon/FQAI/TN/TP outcome columns read.
Counts of reported measurements are not discharge quantity, concentration
at wetlands, or a source-to-receptor flow path.
"""
from __future__ import annotations
import argparse,hashlib,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd

ROOT=Path("results/us_NWCA2021_plant_TNTP_20261110")
CROSS=ROOT/"NWCA2021_568WETLAND_160POTW_EXACT_PERMIT_FEATURE_SOURCE_BEFORE_2021DMR.json"
CODES={"00600":"total_nitrogen_as_N","00665":"total_phosphorus_as_P","50050":"flow_through_plant"}
YEAR=2021
CORE={"EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","PERM_FEATURE_TYPE_CODE",
      "PARAMETER_CODE","MONITORING_PERIOD_END_DATE","DMR_VALUE_NMBR",
      "DMR_VALUE_QUALIFIER_CODE","NODI_CODE","VALUE_TYPE_CODE",
      "STANDARD_UNIT_DESC","DMR_VALUE_ID"}
ADDITIONAL={"DMR_UNIT_DESC","STATISTICAL_BASE_CODE","DMR_VALUE_STANDARD_UNITS","DMR_EVENT_ID","PARAMETER_DESC"}
BATCH=125000

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def norm(s):
    return str(s).strip().upper()

def counted_group():
    return {"original_DMR_2021_event_rows":0,
            "original_DMR_2021_uncensored_numeric_records":0,
            "original_DMR_2021_numeric_positive_records":0,
            "original_monitor_dates_with_valid_uncensored_numeric":set(),
            "original_DMR_reporting_NODI":Counter(),
            "original_DMR_reported_qualifiers":Counter(),
            "original_DMR_unit_type_statistic_combinations":Counter(),
            "original_2021_DMR_value_IDs":set(),
            "original_DMR_source_fiscal_year":Counter()}

def analyze(root, archive,fy):
    origin=json.loads((root/CROSS).read_text())
    if origin["schema"]!="gesui_NWCA2021_568wetland_160municipalPermit_2021DMR_PREOUTCOME_exactPermitFeature_v01":
        raise ValueError("2021 source feature crosswalk not frozen before historical DMR values")
    if origin["2021_real_effluent_TN_TP_FY2021_FY2022_DMR_VALUES_NOT_YET_READ"] is not True:
        raise ValueError("Pre-2021 DMR source selection was not response blind")
    original_pairs=origin["registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology"]
    targets={(norm(x["original_public_NPDES_permit_id"]),norm(x["original_public_NPDES_feature_number"])) for x in original_pairs}
    permits={x[0] for x in targets}
    if len(permits)!=160 or len(targets)!=344:
        raise ValueError("Original 160 permits and 344 feature numbers changed")
    if not archive.is_file():raise ValueError("Original official EPA DMR national archive unavailable")
    source_zip_sha=sha(archive)
    per=defaultdict(counted_group)
    source_info=Counter()
    with zipfile.ZipFile(archive) as z:
        # EPA documentation says singular DMR, but actual official FY2021/22
        # national ZIP members use the plural NPDES_DMRS_FY2021.csv filename.
        csvs=[p for p in z.namelist() if p.split("/")[-1].lower() in (f"npdes_dmr_fy{fy}.csv",f"npdes_dmrs_fy{fy}.csv")]
        if len(csvs)!=1:
            raise ValueError("Official NPDES FY archive missing exactly one national DMR csv; members="+str(z.namelist()[:10]))
        with z.open(csvs[0]) as f:
            columns=list(pd.read_csv(f,nrows=0).columns)
        missing=CORE-set(columns)
        if missing:raise ValueError("Official EPA 2021 DMR source columns missing: "+str(sorted(missing)))
        usecols=CORE|(ADDITIONAL&set(columns))
        with z.open(csvs[0]) as f:
            csv=pd.read_csv(f,low_memory=False,usecols=lambda c:c in usecols,
                            dtype=str,keep_default_na=False,on_bad_lines="error",
                            chunksize=BATCH)
            for chunk_n,block in enumerate(csv,1):
                source_info["whole_original_fiscal_year_DMR_rows"]+=len(block)
                sub=block.loc[block.EXTERNAL_PERMIT_NMBR.str.strip().str.upper().isin(permits)].copy()
                source_info["original_160_candidate_NPDES_permit_DMR_rows_all_features"]+=len(sub)
                if sub.empty:continue
                sub["PERMIT"]=sub.EXTERNAL_PERMIT_NMBR.str.strip().str.upper()
                sub["FEATURE"]=sub.PERM_FEATURE_NMBR.str.strip().str.upper()
                valid=[(p,ft) in targets for p,ft in zip(sub.PERMIT,sub.FEATURE)]
                sub=sub.loc[valid].copy()
                source_info["original_160permit_344exact_feature_DMR_rows_all_years"]+=len(sub)
                if sub.empty:continue
                sub=sub.loc[sub.PERM_FEATURE_TYPE_CODE.str.strip().str.upper().eq("EXO")].copy()
                source_info["only_original_external_outfall_EXO_DMR_rows"]+=len(sub)
                if sub.empty:continue
                dates=pd.to_datetime(sub.MONITORING_PERIOD_END_DATE,errors="coerce",format="mixed")
                sub=sub.loc[dates.ge("2021-01-01")&dates.le("2021-12-31")].copy()
                source_info["EXO_2021_calendar_monitoring_end_DMR_rows"]+=len(sub)
                if sub.empty:continue
                sub["DATE"]=pd.to_datetime(sub.MONITORING_PERIOD_END_DATE,errors="coerce",format="mixed").dt.strftime("%Y-%m-%d")
                sub["CODE"]=sub.PARAMETER_CODE.str.strip().str.zfill(5)
                sub=sub.loc[sub.CODE.isin(CODES)].copy()
                source_info["EXO_2021_calendar_TN_TP_or_flow_rows"]+=len(sub)
                if sub.empty:continue
                # DMR_VALUE_NMBR may be concentration, mass or flow depending
                # on VALUE_TYPE and STANDARD_UNIT_DESC. NO summation.
                sub["NUM"]=pd.to_numeric(sub.DMR_VALUE_NMBR,errors="coerce")
                sub["QUAL"]=sub.DMR_VALUE_QUALIFIER_CODE.str.strip().str.upper()
                sub["NODI"]=sub.NODI_CODE.str.strip().str.upper()
                sub["UNIT"]=sub.STANDARD_UNIT_DESC.str.strip().str.upper()
                sub["VALUE_TYPE"]=sub.VALUE_TYPE_CODE.str.strip().str.upper()
                sub["STATBASE"]=sub.STATISTICAL_BASE_CODE.str.strip().str.upper() if "STATISTICAL_BASE_CODE" in sub else ""
                sub["SOURCE_ID"]=sub.DMR_VALUE_ID.str.strip()
                for r in sub.itertuples(index=False):
                    key=(r.PERMIT,r.FEATURE,r.CODE)
                    a=per[key]
                    a["original_DMR_2021_event_rows"]+=1
                    a["original_DMR_reporting_NODI"][r.NODI if r.NODI else "(blank)"]+=1
                    a["original_DMR_reported_qualifiers"][r.QUAL if r.QUAL else "(blank)"]+=1
                    a["original_DMR_unit_type_statistic_combinations"][(r.UNIT,r.VALUE_TYPE,r.STATBASE)]+=1
                    if r.SOURCE_ID:
                        a["original_2021_DMR_value_IDs"].add(r.SOURCE_ID)
                    a["original_DMR_source_fiscal_year"][str(fy)]+=1
                    # Missing/known not discharged/estimated and censored DMR
                    # cannot masquerade as numerical 2021 effluent.
                    if (pd.notna(r.NUM) and math.isfinite(r.NUM) and r.NUM>=0
                            and r.NODI in ("","0") and r.QUAL in ("","=","E")):
                        a["original_DMR_2021_uncensored_numeric_records"]+=1
                        if r.NUM>0:a["original_DMR_2021_numeric_positive_records"]+=1
                        a["original_monitor_dates_with_valid_uncensored_numeric"].add(r.DATE)
                if chunk_n%80==0:
                    print("DMR_READ_FY",fy,"CHUNK",chunk_n,
                          "nationwide_rows",source_info["whole_original_fiscal_year_DMR_rows"],
                          "exact_feature_2021_target_NPflow",source_info["EXO_2021_calendar_TN_TP_or_flow_rows"],flush=True)
    structured=[]
    for (p,feat,code),v in sorted(per.items()):
        structured.append({
            "original_public_NPDES_permit_id":p,
            "original_public_NPDES_feature_number":feat,
            "DMR_2021_source_pollutant_code":code,
            "DMR_2021_source_pollutant_name":CODES[code],
            "source_2021_FY"+str(fy)+"_EXO_reported_rows":v["original_DMR_2021_event_rows"],
            "source_2021_uncensored_numeric_records":v["original_DMR_2021_uncensored_numeric_records"],
            "source_2021_positive_numeric_records":v["original_DMR_2021_numeric_positive_records"],
            "source_2021_valid_numeric_monitoring_end_dates":sorted(v["original_monitor_dates_with_valid_uncensored_numeric"]),
            "source_2021_original_DMR_valueIDs_unique_n":len(v["original_2021_DMR_value_IDs"]),
            "source_2021_NODI_code_counts":dict(v["original_DMR_reporting_NODI"]),
            "source_2021_qualifier_counts":dict(v["original_DMR_reported_qualifiers"]),
            "source_2021_unit_value_type_statistic_counts":{
               "|".join(key):n for key,n in sorted(v["original_DMR_unit_type_statistic_combinations"].items())},
            "source_original_fiscal_year_rows":dict(v["original_DMR_source_fiscal_year"]),
            "value_units_not_homogeneous_DO_NOT_SUM_as_NP_loading":True
        })
    unique=set((v["original_public_NPDES_permit_id"],v["original_public_NPDES_feature_number"]) for v in structured)
    originalN=set((v["original_public_NPDES_permit_id"],v["original_public_NPDES_feature_number"])
                  for v in structured if v["DMR_2021_source_pollutant_code"]=="00600" and v["source_2021_uncensored_numeric_records"]>0)
    originalP=set((v["original_public_NPDES_permit_id"],v["original_public_NPDES_feature_number"])
                  for v in structured if v["DMR_2021_source_pollutant_code"]=="00665" and v["source_2021_uncensored_numeric_records"]>0)
    return {
        "schema":"gesui_US_NWCA2021_160POTW_344feature_EPA_FY_DMR_NP_reported_source_only_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_FY":fy,
        "original_source_568wetland_160permit_344feature_PRE_DMR_SHA256":sha(root/CROSS),
        "original_official_EPA_ICIS_NPDES_fiscal_DMR_zip_url":f"https://echo.epa.gov/files/echodownloads/npdes_dmrs_fy{fy}.zip",
        "original_official_EPA_ICIS_NPDES_fiscal_DMR_zip_SHA256":source_zip_sha,
        "original_fiscal_DMR_national_row_and_filter_QC":dict(source_info),
        "original_total_160_registered_candidate_permits":len(permits),
        "original_total_344_candidate_permitted_features":len(targets),
        "actual_2021_calendar_year_EXO_344_features_with_any_DMR_target_parameter":len(unique),
        "actual_2021_calendar_year_EXO_344_features_with_measured_TN_00600":len(originalN),
        "actual_2021_calendar_year_EXO_344_features_with_measured_TP_00665":len(originalP),
        "actual_2021_calendar_year_EXO_344_features_with_both_numeric_TN_and_TP":len(originalN&originalP),
        "all_exact_original_registered_permit_FEATURE_DMR_2021_reporting_records":structured,
        "DMR_FY_2021_or_2022_file_MAY_INCLUDE_dates_outside_calendar2021_filtered":True,
        "reported_numerical_DMR_not_verified_outfall_to_wetland_flow_or_plant_nutrient_dose":True,
        "exact_EPA_2021_NODI_missing_not_real_zero":True,
        "original_NWCA2021_plant_H_ALL_FQAI_measured_NP_values_NOT_READ":True,
        "source_FY2021_FY2022_duplicates_not_yet_merged_or_counted_twice":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--dmr-zip",type=Path,required=True)
    p.add_argument("--fy",type=int,choices=(2021,2022),required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    o=analyze(root,a.dmr_zip,a.fy)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(o,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
      "fy":o["source_FY"],
      "filter":o["original_fiscal_DMR_national_row_and_filter_QC"],
      "TN_features":o["actual_2021_calendar_year_EXO_344_features_with_measured_TN_00600"],
      "TP_features":o["actual_2021_calendar_year_EXO_344_features_with_measured_TP_00665"],
      "both_features":o["actual_2021_calendar_year_EXO_344_features_with_both_numeric_TN_and_TP"]
    },indent=2))

if __name__=="__main__":
    main()
