#!/usr/bin/env python3
"""Source-only national UK PRTR 2017/18 WWTP water-release N/P mass census.

IMPORTANT 2026 PRTR original CSV file names do NOT identify ReportingYear
reliably. For example the 2018-named CSV is mostly ReportingYear 2017.
Download year-NAMED 2017, 2018 AND 2019 official government source files,
filter BY original ReportingYear, and explicitly audit overlapping releases
and contradictory facility-pollutant-year values before aggregating anything.

Use actual documented PRTR primary activity "5.(f)" (urban waste-water
treatment plants >=100,000 PE), NACE "37.00" and exact WATER release codes.
TOTAL - PHOSPHORUS and TOTAL - NITROGEN are distinct original substances.
Do not pool NH3 or atmospheric nitrogen with TOTAL - NITROGEN. Official
TotalQuantity E-PRTR is annual mass kg/year; Quantity is waste-transfer
weight and MUST NOT be substituted. M/C/E are distinct basis codes.
Facility coordinates are NOT physical river outfalls or receiving reach.
Missing report, offsite transfer or under-threshold release is not 0 load.

Only publication source hashes, aggregate counts, pollutant/method class
distributions and crossfile conflicts are emitted. No facility identifiers,
raw source rows, coordinates or reported quantities are rehosted.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,math,json,tempfile
from pathlib import Path

import source_only_2017_2018_prtr_schema as source

FILES={
  "2017":source.SOURCES["2017"],
  "2018":source.SOURCES["2018"],
  "2019":"https://assets.publishing.service.gov.uk/media/6a633032243ee000f7127f32/uk_prtr_dataset_2019.csv",
}
VALID_YEARS=("2017","2018")
WWTP_ACTIVITY="5.(f)"
WWTP_NACE="37.00"
TARGET_POLLUTANTS=("TOTAL - PHOSPHORUS","TOTAL - NITROGEN")
MEDIUM="WATER"
REQUIRED={
    "NationalID","ReportingYear","NACEMainEconomicActivityCode",
    "MainEconomicActivityName","AnnexIActivityCode","MediumCode",
    "PollutantCode","MethodBasisCode","MethodTypeCode",
    "TotalQuantity","Quantity","CountryID","LatitudeMeasure","LongitudeMeasure"
}
def numeric(x):
    try:
        y=float(str(x or "").replace(",","").strip())
        return y if math.isfinite(y) else None
    except (ValueError,TypeError):return None

def readcsv(path):
    with path.open("r",encoding="utf-8-sig",errors="replace",newline="") as f:
        rd=csv.DictReader(f)
        fields=set(rd.fieldnames or [])
        if not REQUIRED.issubset(fields):
            raise RuntimeError("Official PRTR CSV columns changed: "+repr(REQUIRED-fields))
        rows=list(rd)
    if len(rows)>300000 or len(rows)<100:
        raise RuntimeError("Unexpected original PRTR original file row count")
    return rows

def native_key(row):
    return (row["ReportingYear"].strip(),
            row["NationalID"].strip(),
            row["PollutantCode"].strip(),
            row["MediumCode"].strip(),
            row["Designation"].strip())

def classify(row):
    act=row["AnnexIActivityCode"].strip()
    nace=row["NACEMainEconomicActivityCode"].strip()
    if act==WWTP_ACTIVITY and nace==WWTP_NACE:return "strict_5f_and_NACE37"
    if act==WWTP_ACTIVITY:return "5f_only_NACE_other"
    if nace==WWTP_NACE:return "NACE37_only_non5f"
    if act.startswith("5.(g)"):return "industrial_wastewater_5g"
    return "other"

def census(a):
    all_rows=[]
    by_file={}
    for labelledyear,url in FILES.items():
        with tempfile.TemporaryDirectory(prefix="gesui_UKPRTR_"+labelledyear+"_") as temp:
            target=Path(temp)/"official.csv"
            receipt=source.download(url,target)
            info={"status":receipt["status"],
                  "url":url,"http":receipt.get("http"),
                  "error":receipt.get("error"),
                  "sha":receipt.get("source_sha256"),
                  "file_bytes":receipt.get("source_bytes")}
            if receipt["status"]!="OFFICIAL_ORIGINAL_CSV_BYTES_OBTAINED":
                by_file[labelledyear]=info
                continue
            rows=readcsv(target)
        info["original_rows"]=len(rows)
        info["actual_ReportingYear_counts"]=dict(sorted(collections.Counter(
            r["ReportingYear"].strip() for r in rows).items()))
        info["2017_2018_water_PN_report_rows"]=collections.Counter()
        for row in rows:
            yr=row["ReportingYear"].strip()
            if yr in VALID_YEARS:
                if row["MediumCode"].strip()==MEDIUM and row["PollutantCode"].strip() in TARGET_POLLUTANTS:
                    info["2017_2018_water_PN_report_rows"][yr+"|"+row["PollutantCode"].strip()]+=1
                all_rows.append((labelledyear,row))
        info["2017_2018_water_PN_report_rows"]=dict(info["2017_2018_water_PN_report_rows"])
        by_file[labelledyear]=info

    result={
       "schema":"GESUI_UK_PRTR_WWTW_2017_18_FACILITY_YEAR_WATER_PN_SOURCE_GATE_V1",
       "official_publisher":"DEFRA UK PRTR csv 2017/2018/2019 published 2026",
       "all_original_files_publication":"https://www.gov.uk/government/publications/download-industrial-emissions-and-waste-transfer-uk-prtr-data",
       "source_only":True,
       "source_lock_internal_reporting_year_not_CSV_filename":True,
       "pollutant_codes":list(TARGET_POLLUTANTS),
       "WWTW_sector_activity":WWTP_ACTIVITY,
       "WWTW_sector_NACE":WWTP_NACE,
       "reported_mass_unit":"kg/year for PRTR TotalQuantity where MediumCode WATER",
       "NOT_annual_STW_outfall_mass":True,
       "NO_source_to_ChemPop_river_routing_or_BMWP_effect":True,
       "NO_zeros_for_missing_PRTR_reporting_threshold":True,
       "source_files":by_file,
       "status":"SOURCE_STOP",
       "by_true_year":{},
    }
    if any(v["status"]!="OFFICIAL_ORIGINAL_CSV_BYTES_OBTAINED" for v in by_file.values()):
        result["reason"]="AT_LEAST_ONE_ORIGINAL_GOVUK_CSV_SOURCE_NOT_RETRIEVED"
        return result

    # Original NationalID represents facility identifier, not outfall ID.
    by_key=collections.defaultdict(list)
    source_report_count=collections.Counter()
    summary={}
    for file_label,row in all_rows:
        yr=row["ReportingYear"].strip()
        if row["MediumCode"].strip()!=MEDIUM or row["PollutantCode"].strip() not in TARGET_POLLUTANTS:
            continue
        cl=classify(row)
        source_report_count[(yr,cl,row["PollutantCode"].strip(),file_label,row["MethodBasisCode"].strip())]+=1
        if not row["NationalID"].strip():
            raise RuntimeError("Original water P/N release missing NationalID")
        by_key[native_key(row)].append((file_label,row))

    keys_status=collections.Counter()
    records_by_yr=collections.defaultdict(list)
    exact_crossfile_duplicates=0
    crossfile_conflict=0
    ambiguous_overlap_examples_hashed=[]
    for key,vals in by_key.items():
        rows=[r for _,r in vals]
        # Some EU releases can be non-water / other designation but this key
        # explicitly locks origin release medium, pollutant and designation.
        # Compare scientific inputs, not source wrapper or facility free text.
        scientific=lambda r:(
            r["AnnexIActivityCode"].strip(),r["NACEMainEconomicActivityCode"].strip(),
            r["MethodBasisCode"].strip(),r["MethodTypeCode"].strip(),
            r["TotalQuantity"].strip(),r["CountryID"].strip()
        )
        fingerprints={scientific(r) for r in rows}
        labels={src for src,_ in vals}
        if len(vals)>1:
            if len(fingerprints)==1:
                exact_crossfile_duplicates+=len(vals)-1
                keys_status["duplicate_orig_full_scientific_record_same_values"]+=1
            else:
                crossfile_conflict+=1
                keys_status["same_year_facility_PN_medium_conflicting_original_release_values"]+=1
                if len(ambiguous_overlap_examples_hashed)<8:
                    ambiguous_overlap_examples_hashed.append(hashlib.sha256(
                        ("|".join(key)).encode()).hexdigest())
                continue
        r=rows[0]
        if len(vals)==1:keys_status["single_raw_record"]+=1
        v=numeric(r["TotalQuantity"])
        missing_mass=not (v is not None and v>=0)
        if missing_mass:keys_status["water_PN_report_without_valid_TotalQuantity"]+=1
        if v is not None and v==0:keys_status["genuine_zero_reported_total_mass"]+=1
        if v is not None and v>0:keys_status["positive_reported_total_mass"]+=1
        records_by_yr[key[0]].append({
            "id":key[1],"target":key[2],"type":classify(r),
            "method":r["MethodBasisCode"].strip(),
            "quantity":v,"inputfile":sorted(labels),
            "country":r["CountryID"].strip(),
            "authority":r.get("CompetentAuthorityPartyName","").strip(),
            "nuts_region":r.get("NutsRegionID","").strip(),
            "river_basin":r.get("RiverBasinDistrictID","").strip(),
            "latitude":numeric(r.get("LatitudeMeasure","")),
            "longitude":numeric(r.get("LongitudeMeasure","")),
            "designation":key[4]
        })
    result["overlap_sourceQA"]={
        "original_water_PN_raw_record_n":sum(len(v) for v in by_key.values()),
        "unique_year_facility_pollutant_medium_designation_native_keys":len(by_key),
        "exact_source_scientific_duplicates_beyond_first":exact_crossfile_duplicates,
        "n_original_release_native_keys_CONFLICTING_across_files":crossfile_conflict,
        "conflicting_native_key_SHA256_first8":ambiguous_overlap_examples_hashed,
        "source_record_status":dict(keys_status),
        "source_reported_rows_by_year_facilitytype_pollutant_file_method":{
            "|".join(keys):value for keys,value in sorted(source_report_count.items())},
    }
    for yr in VALID_YEARS:
        rel=records_by_yr[yr]
        classified={}
        for category in ("strict_5f_and_NACE37","5f_only_NACE_other",
                         "NACE37_only_non5f","industrial_wastewater_5g","other"):
            target={}
            for pollutant in TARGET_POLLUTANTS:
                subset=[x for x in rel if x["type"]==category and x["target"]==pollutant]
                positive=[x for x in subset if x["quantity"] is not None and x["quantity"]>0]
                target[pollutant]={
                    "original_distinct_facilities":len({x["id"] for x in subset}),
                    "original_water_release_rows":len(subset),
                    "nonnegative_TotalQuantity_rows":sum(x["quantity"] is not None and x["quantity"]>=0 for x in subset),
                    "strict_positive_water_PN_report_rows":len(positive),
                    "distinct_facilities_with_positive_water_release":len({x["id"] for x in positive}),
                    "basis_method_M_C_E_original_rows":dict(collections.Counter(x["method"] for x in subset)),
                    "facility_mass_rows_with_M_basis":sum(x["method"]=="M" and x["quantity"] is not None and x["quantity"]>0 for x in subset),
                    "original_CountryID_source":dict(collections.Counter(x["country"] for x in subset)),
                    "original_competent_authority_counts":dict(collections.Counter(x["authority"] for x in subset)),
                    "original_NUTS_region_code_counts":dict(collections.Counter(x["nuts_region"] for x in subset)),
                    "original_river_basin_code_counts":dict(collections.Counter(x["river_basin"] for x in subset)),
                    "n_original_nonmissing_valid_facility_lat_lon":sum(
                        x["latitude"] is not None and x["longitude"] is not None
                        and -12<=x["longitude"]<=5 and 49<=x["latitude"]<=61
                        for x in subset),
                    "facility_lat_lon_is_NOT_permitted_outfall":True,
                    "no_release_sum_until_source_overlap_and_outlet_identity_validated":True,
                }
            p={x["id"] for x in rel if x["type"]==category and x["target"]==TARGET_POLLUTANTS[0] and x["quantity"] is not None and x["quantity"]>0}
            n={x["id"] for x in rel if x["type"]==category and x["target"]==TARGET_POLLUTANTS[1] and x["quantity"] is not None and x["quantity"]>0}
            target["same_native_facility_with_both_positive_water_P_and_N"]=len(p&n)
            target["native_facilities_with_either_positive_P_or_N"]=len(p|n)
            classified[category]=target
        result["by_true_year"][yr]=classified
    result["source_support_gate"]={
        "correct_original_2017_18_reporting_year_recovered":True,
        "all_three_official_file_wrappers_retrieved":True,
        "source_year_report_rows_may_span_wrong_named_file":True,
        "duplicate_original_releases_across_file_wrappers_detected":exact_crossfile_duplicates>0,
        "crossfile_contradictions_absent":crossfile_conflict==0,
        "WWTW_5f_nace37_water_P_and_N_joint_records_in_both_years":
            all(result["by_true_year"][yr]["strict_5f_and_NACE37"]["same_native_facility_with_both_positive_water_P_and_N"]>0 for yr in VALID_YEARS),
        "native_facility_code_not_physical_permitted_outfall":True,
        "PRTR_annual_total_emission_method_M_C_E_not_all_measured":True,
        "not_2017_18_actual_routed_water_PN_load_or_causal_EA_outfall":True,
        "NO_BMWP_BIO_MODEL_AUTHORIZED":True,
    }
    result["status"]="ORIGINAL_2017_2019_WRAPPER_2017_18_TRUE_YEAR_WWTW_WATER_PN_CENSUS_COMPLETED"
    return result

def main(a):
    data=census(a)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(data,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_TRUE_YEAR_UKPRTR_WWTW_PN_SOURCE="+json.dumps({
        "status":data["status"],"reason":data.get("reason"),
        "files":{k:{"status":v["status"],"year_counts":v.get("actual_ReportingYear_counts"),
                    "sha":v.get("sha"),"error":v.get("error"),
                    "year_release_rows":v.get("2017_2018_water_PN_report_rows")}
                 for k,v in data["source_files"].items()},
        "source_conflicts":data.get("overlap_sourceQA"),
        "by_true_year":data.get("by_true_year"),
        "gates":data.get("source_support_gate"),
    },sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
