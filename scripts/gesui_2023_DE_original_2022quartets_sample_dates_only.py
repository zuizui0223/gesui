#!/usr/bin/env python3
"""2023 individual WISE6 river TP *sampling-calendar only* at the ORIGINAL
eight German four-arm WWTP/registered-sourcezero 2022 quartets.

The 2022 8×4 station cohort and source labels were fixed independently of
original 2023 outcomes. Never parse resultObservedValue, below-LOQ numeric
concentration or any 2023 butterfly observation. Verify official original
2026 public ZIP SHA256 and stream first 11 source metadata columns only.
"""
from __future__ import annotations
import argparse,hashlib,json,urllib.request,zipfile
from datetime import datetime,timezone
from collections import Counter
from pathlib import Path

URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
SOURCE_ZIP_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
SOURCE_2022_DATE_BLOB="c09582bee2c09ad038bd9fd16926d837644da075"
ZIP_BYTES=2024284753
CSV_BYTES=36043700357
ROWS=102569057
YEAR="2023"
ANALYTE="CAS_7723-14-0"
UNIT="mg{P}/L"
MIN_SITE_DAYS=4
MIN_SITE_MONTHS=3
MIN_FOURARM_MONTHS=3
MIN_MONTH_QUARTETS=5
MIN_MONTH_BLOCKS=4
MIN_FOURARM_DAYS=3
MIN_DAY_QUARTETS=4
MIN_DAY_BLOCKS=3

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for ch in iter(lambda:f.read(1048576),b""):h.update(ch)
    return h.hexdigest()
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def sitehash(country,scheme,site):
    return hashlib.sha256((country+"|"+scheme+"|"+site).encode()).hexdigest()
def val(t):
    return t.strip().strip(b'"').decode("utf-8",errors="replace")

def source_read(path):
    if blob(path)!=SOURCE_2022_DATE_BLOB:
        raise ValueError("Original pre-2023 2022 source-only 8-quartet date gate altered")
    j=json.loads(path.read_text())
    if (j["source_quartets_ge3_same_month_n"]!=7 or
        j["source_quartets_ge3_same_day_n"]!=0 or
        not j["source_same_calendar_month_fourstation_gate_PASS"] or
        not j["NO_original_TP_numeric_resultObservedValue_or_LOQ_values_read"]):
        raise ValueError("Original eight-quartet 2022 source/date selection not verified")
    q=j["all_original_8_no_replacement_quartets_and_date_coverage"]
    if len(q)!=8:raise ValueError("Unexpected 2022 original number of quartets")
    ids={}
    for x in q:
        arms=x["four_arm_original_WISE_country_scheme_site_hashes"]
        if len(arms)!=4 or set(arms)!={
         "source_upstream","source_downstream",
         "registered_zero_upstream","registered_zero_downstream"}:
            raise ValueError("Unexpected 2022 exact original quartet structure")
        for side,hashid in arms.items():
            if hashid in ids:raise ValueError("Original source quartet reuses a WISE station")
            ids[hashid]={"original_2022_quartet_index":x["quartet_index"],"arm":side}
    if len(ids)!=32:raise ValueError("Not original 32 2022 source-only frozen WISE monitoring sites")
    return q,ids,j

def download(path):
    req=urllib.request.Request(URL,headers={
       "User-Agent":"gesui-2023-original-WISE6-DE-quartet-source-date-only/1.0",
       "Accept":"application/zip,*/*"})
    h=hashlib.sha256();size=0
    with urllib.request.urlopen(req,timeout=200) as r,path.open("wb") as f:
        if r.status!=200:raise ValueError("Official 2026 public WISE6 source not HTTP 200")
        for ch in iter(lambda:r.read(1048576),b""):
            size+=len(ch)
            if size>ZIP_BYTES:raise ValueError("Official original source ZIP length unexpectedly changed")
            h.update(ch);f.write(ch)
    if size!=ZIP_BYTES or h.hexdigest()!=SOURCE_ZIP_SHA:
        raise ValueError("Official 2026 original WISE6 ZIP SHA mismatch")
    return {"official_source_url":URL,"verified_original_archive_SHA256":h.hexdigest(),
             "original_public_compressed_bytes":size}

def scan(path,ids):
    samples={k:{"days":set(),"months":set(),"distinct_sampleids":set(),
                "source_2023_date_only_rows":0} for k in ids}
    count=Counter()
    with zipfile.ZipFile(path) as z:
        if z.namelist()!=["waterbase_t_wise6_disaggregateddata.csv"]:
            raise ValueError("Original full WISE6 ZIP member differs")
        info=z.infolist()[0]
        if info.file_size!=CSV_BYTES:raise ValueError("Original 36GB source member changed")
        with z.open(info) as f:
            head=f.readline().decode("utf-8-sig").strip().split(";")
            expected=("countryCode monitoringSiteIdentifier monitoringSiteIdentifierScheme "
                     "parameterWaterBodyCategory observedPropertyDeterminandCode "
                     "observedPropertyDeterminandLabel procedureAnalysedMatrix resultUom "
                     "phenomenonTimeSamplingDate phenomenonTimeReferenceYear sampleIdentifier resultObservedValue").split()
            if head[:12]!=expected:raise ValueError("Original 29-field WISE6 header changed")
            for line in f:
                count["all_original_EEA2026_sample_rows"]+=1
                if not line.startswith((b"DE;",b'"DE";')):
                    continue
                # ONLY first 11 source/date fields; original resultObservedValue
                # lives in unexamined last token, NOT parsed/decoded.
                tok=line.split(b";",11)
                if len(tok)<12:continue
                if (val(tok[9])!=YEAR or val(tok[3])!="RW" or
                    val(tok[4])!=ANALYTE or val(tok[5])!="Total phosphorus" or
                    val(tok[6])!="W" or val(tok[7])!=UNIT):
                    continue
                day=val(tok[8])[:10]
                if len(day)!=10 or not day.startswith("2023-") or day[4]!="-" or day[7]!="-":
                    count["original_2023_DE_RW_TP_rows_invalid_sampling_day"]+=1
                    continue
                count["original_2023_DE_RW_TP_all_published_sample_rows"]+=1
                h=sitehash("DE",val(tok[2]),val(tok[1]))
                if h not in ids:continue
                q=samples[h]
                q["source_2023_date_only_rows"]+=1
                q["days"].add(day)
                q["months"].add(day[:7])
                sid=val(tok[10])
                if sid:q["distinct_sampleids"].add(sid)
                count["original_frozen32_site_2023_source_sample_rows"]+=1
    if count["all_original_EEA2026_sample_rows"]!=ROWS:
        raise ValueError("Original WISE6 all-source sampled row denominator drift")
    return samples,dict(count)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-2022-eightquartet-source",required=True,type=Path)
    ap.add_argument("--tmpdir",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    quartets,ids,original=source_read(a.original_2022_eightquartet_source)
    a.tmpdir.mkdir(parents=True,exist_ok=True)
    zipp=a.tmpdir/"WISE6_original_public_2023_DE_TP_dates.zip"
    dl=download(zipp)
    sampling,counters=scan(zipp,ids)
    frozen={}
    for h,item in sorted(sampling.items()):
        frozen[h]={
          "original_2022_sitehash_country_scheme_ID":h,
          "original_2022_quartet_index":ids[h]["original_2022_quartet_index"],
          "original_arm_role":ids[h]["arm"],
          "original_2023_total_TP_sample_date_metadata_rows":item["source_2023_date_only_rows"],
          "original_2023_distinct_qualified_source_days":sorted(item["days"]),
          "original_2023_distinct_qualified_source_months":sorted(item["months"]),
          "original_2023_distinct_sample_ID_count":len(item["distinct_sampleids"]),
          "source_2023_site_4days_3months_meta_gate":
              len(item["days"])>=MIN_SITE_DAYS and len(item["months"])>=MIN_SITE_MONTHS,
          "never_accessed_numeric_TP":True,
        }
    quartet_stats=[]
    bymonth=[];byday=[]
    for q in quartets:
        arms=q["four_arm_original_WISE_country_scheme_site_hashes"]
        ds=[sampling[h]["days"] for h in arms.values()]
        ms=[sampling[h]["months"] for h in arms.values()]
        day=sorted(set.intersection(*ds))
        month=sorted(set.intersection(*ms))
        site_ok=all(frozen[h]["source_2023_site_4days_3months_meta_gate"] for h in arms.values())
        okmonth=site_ok and len(month)>=MIN_FOURARM_MONTHS
        okday=site_ok and len(day)>=MIN_FOURARM_DAYS
        row={
         "original_2022_quartet_index":q["quartet_index"],
         "original_2022_same_river_country2degree_source_block":q["source_country2degree_upstream_block"],
         "original_2022_source_interposed_reported_P_tonnes_per_year":
                 q["original_source_interposed_P_tonnes_per_year"],
         "original_WISE_source_only_fixed_4site_hashes":arms,
         "source2023_site_days_n":{role:len(sampling[h]["days"]) for role,h in arms.items()},
         "source2023_site_months_n":{role:len(sampling[h]["months"]) for role,h in arms.items()},
         "original_four_arm_shared_2023_calendar_days":day,
         "original_four_arm_shared_2023_calendar_months":month,
         "2023_four_arms_all_ge4days_and_ge3months":site_ok,
         "source_2023_four_arm_ge3_same_months_gate":okmonth,
         "source_2023_four_arm_ge3_same_days_gate":okday,
         "no_2023_TP_numeric_values_or_WWTP_actual_2023_loads_accessed":True,
        }
        quartet_stats.append(row)
        if okmonth:bymonth.append(row)
        if okday:byday.append(row)
    month_blocks=set(q["original_2022_same_river_country2degree_source_block"] for q in bymonth)
    day_blocks=set(q["original_2022_same_river_country2degree_source_block"] for q in byday)
    month_pass=len(bymonth)>=MIN_MONTH_QUARTETS and len(month_blocks)>=MIN_MONTH_BLOCKS
    day_pass=len(byday)>=MIN_DAY_QUARTETS and len(day_blocks)>=MIN_DAY_BLOCKS
    out={
     "schema":"gesui_2023_original_WISE6_DE_TP_at_2022_frozen_WWTP_source_control_quartets_date_only_gate_v01",
     "utc":datetime.now(timezone.utc).isoformat(),
     "POSTHOC_to_known_2022_same_month_positive_riverTP_gradients":True,
     "source_WISE6_original_verified_2026_ZIP":dl,
     "original_2022_date_only_8quartet_source_Git_blob_SHA1":SOURCE_2022_DATE_BLOB,
     "source_2022_8quartet_source_file_SHA256":sha(a.original_2022_eightquartet_source),
     "original_2022_source_quartets_n":len(quartets),
     "original_2022_source_sites_n":len(ids),
     "year_of_NEW_2023_actual_water_sample_source":YEAR,
     "2023_original_country_all_river_TP_sample_metadata_counts":counters,
     "source_only_original32_site_2023_sampling_date_metadata":frozen,
     "original8_unreplaced_quartets_2023_fourarm_sample_time_coverage":quartet_stats,
     "min_site_days":MIN_SITE_DAYS,"min_site_months":MIN_SITE_MONTHS,
     "strict_4arm_same2023_day_min_days":MIN_FOURARM_DAYS,
     "strict_4arm_same2023_month_min_months":MIN_FOURARM_MONTHS,
     "strict_2023_same_day_quartets_n":len(byday),
     "strict_2023_same_month_quartets_n":len(bymonth),
     "strict_2023_same_day_source_blocks_n":len(day_blocks),
     "strict_2023_same_month_source_blocks_n":len(month_blocks),
     "strict_2023_frozen_4arm_same_DAY_numeric_gate_PASS":day_pass,
     "strict_2023_frozen_4arm_same_MONTH_numeric_gate_PASS":month_pass,
     "never_parsed_2023_resultObservedValue_LOQ_or_WWTP_chemical_load":True,
     "2022_annual_reported_WWTP_source_not_2023_treatment_or_daily_dose":True,
     "no_lycaenidae_abundance_photo_puddling_or_host_quality_accessed":True,
     "if_month_gate_fails_no_new_2023_WISE6_numeric_TP_gradient_permitted":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_2023_DE_WISE6_TP_ORIGINAL32_QUARTETS_SOURCE_DATE_GATE",
      json.dumps({"original_source_2022_quartets":len(quartets),
         "original32_site_2023_sample_rows":counters["original_frozen32_site_2023_source_sample_rows"],
         "quartets_ge3_shared_months":len(bymonth),
         "quartets_ge3_shared_days":len(byday),
         "months_geoblocks":len(month_blocks),
         "day_geoblocks":len(day_blocks),
         "2023_month_numeric_GATE_PASS":month_pass,
         "2023_day_numeric_GATE_PASS":day_pass,
         "no_2023_numeric_TP_opened":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
