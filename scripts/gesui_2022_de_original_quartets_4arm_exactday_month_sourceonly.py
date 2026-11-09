#!/usr/bin/env python3
"""2022 DE WISE6 TP same-DAY / same-calendar-MONTH source eligibility, no TP values.

Read ONLY first 11 original EEA2026 WISE6 public fields, from SHA-verified
2GB ZIP/36GB streaming CSV, for EXACT original 32 station SHA IDs in
precommitted WWTP source-interposed/source-zero quartets.
Numeric resultObservedValue and LOQ columns are never parsed or used.
"""
from __future__ import annotations
import argparse,hashlib,json,urllib.request,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
ARCHIVE_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
SOURCE_QUARTET_BLOB="8bd67be25968311ca16be4026888daa783fbe69f"
SOURCE_NUMERIC32_SHA="28b679184cbebf6376f1c65d7c160f7511f057ff22787e98ff3fa265b617b261"
ZIP_SIZE=2024284753
CSV_SIZE=36043700357
EXPECTED_ROWS=102569057
EXPECTED_SAMPLE32=32
TP="CAS_7723-14-0"
UNIT="mg{P}/L"
MIN_ALLFOUR_DAY=3
MIN_QUARTETS_DAY=4
MIN_BLOCKS_DAY=3
MIN_ALLFOUR_MONTH=3
MIN_QUARTETS_MONTH=5
MIN_BLOCKS_MONTH=4

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def get_original(path):
    r=urllib.request.Request(URL,headers={"User-Agent":"gesui-2022-same-day-river-TP-monitor-meta-gate/1.0"})
    h=hashlib.sha256();n=0
    with urllib.request.urlopen(r,timeout=150) as s,path.open("wb") as f:
        if s.status!=200:raise ValueError("Original official WISE6 ZIP HTTP invalid")
        for chunk in iter(lambda:s.read(1024*1024),b""):
            f.write(chunk);h.update(chunk);n+=len(chunk)
            if n>ZIP_SIZE:raise ValueError("Source ZIP exceeded original source-size lock")
    if n!=ZIP_SIZE or h.hexdigest()!=ARCHIVE_SHA:
        raise ValueError("Source original WISE6 2026 archive SHA mismatch")
    return {"source_url":URL,"sha256":h.hexdigest(),"original_zip_bytes":n}
def val(b):return b.strip().strip(b'"').decode("utf8",errors="replace")
def hkey(c,scheme,site):
    return hashlib.sha256((c+"|"+scheme+"|"+site).encode()).hexdigest()
def read_quartets(q,n):
    if blob(q)!=SOURCE_QUARTET_BLOB:
        raise ValueError("Original pre-TP numeric quartet source-only Git blob changed")
    if sha(n)!=SOURCE_NUMERIC32_SHA:
        raise ValueError("Original frozen 32-site TP numeric QC checksum changed")
    source=json.loads(q.read_text())
    num=json.loads(n.read_text())
    selected=source["selected_original_8_source_vs_clean_zero_monitor_quartets"]
    numeric=num["source_original_32_monitoring_sites_2022_numeric_TP_QC"]
    if len(selected)!=8 or len(numeric)!=32:
        raise ValueError("Original 8 quartets or source 32 TP station count drift")
    quality={r["original_country_scheme_site_SHA256"]:r for r in numeric}
    if len(quality)!=32:raise ValueError("Original 32 numeric TP site hashes not unique")
    quartets={}
    expected=set()
    for i,q in enumerate(selected,1):
        src=q["source"];ctrl=q["registered_source_zero"]
        ids={
          "source_upstream":src["upstream_station_SHA256"],
          "source_downstream":src["downstream_station_SHA256"],
          "registered_zero_upstream":ctrl["upstream_station_SHA256"],
          "registered_zero_downstream":ctrl["downstream_station_SHA256"]}
        for p in ids.values():
            if p in expected:raise ValueError("Original 8 source quartets share a station")
            expected.add(p)
            if p not in quality:raise ValueError("Original source-station ID not included in 32 numeric QC")
        quartets[i]={"quartet_index":i,
           "source_country2degree_upstream_block":src["upstream_country2degree_block"],
           "four_arm_original_WISE_country_scheme_site_hashes":ids,
           "original_source_interposed_P_tonnes_per_year":src["real_exact_source_between"]["strictly_between_original_P_reported_t_year"],
           "original_control_registered_interposed_source_count":ctrl["real_exact_source_between"]["strictly_between_total_original_Fw_outfalls_n"],
           "source_only_match_cost":q["value"],
        }
    if expected!=set(quality):raise ValueError("Quartet source population not 32 frozen original numeric sites")
    return quartets,quality

def scan_metadata(path,ids):
    records={key:{"date_days":set(),"sample_rows":0,"sample_ids":set()} for key in ids}
    qc=Counter()
    with zipfile.ZipFile(path) as z:
        if z.namelist()!=["waterbase_t_wise6_disaggregateddata.csv"]:
            raise ValueError("Unexpected source CSV ZIP member")
        info=z.infolist()[0]
        if info.file_size!=CSV_SIZE:raise ValueError("Original 2026 36GB source member size drift")
        with z.open(info) as f:
            header=f.readline().decode("utf-8-sig").strip().split(";")
            expected=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
             "parameterWaterBodyCategory","observedPropertyDeterminandCode",
             "observedPropertyDeterminandLabel","procedureAnalysedMatrix","resultUom",
             "phenomenonTimeSamplingDate","phenomenonTimeReferenceYear","sampleIdentifier",
             "resultObservedValue"]
            if header[:12]!=expected:raise ValueError("Original WISE6 source field columns drift")
            for line in f:
                qc["original_all_source_rows"]+=1
                if not (line.startswith((b"DE;",b'"DE";'))):continue
                tok=line.split(b";",11)
                if len(tok)<12:continue
                if val(tok[9])!="2022" or val(tok[3])!="RW" or val(tok[4])!=TP or val(tok[5])!="Total phosphorus" or val(tok[6])!="W" or val(tok[7])!=UNIT:
                    continue
                dt=val(tok[8])[:10]
                if len(dt)!=10 or dt[:4]!="2022" or dt[4]!="-" or dt[7]!="-":continue
                qc["DE_2022_RW_TP_source_sample_rows"]+=1
                id=hkey("DE",val(tok[2]),val(tok[1]))
                if id not in ids:continue
                r=records[id]
                r["sample_rows"]+=1
                r["date_days"].add(dt)
                sample=val(tok[10])
                if sample:r["sample_ids"].add(sample)
                qc["original_32_station_meta_only_sample_rows"]+=1
    if qc["original_all_source_rows"]!=EXPECTED_ROWS:
        raise ValueError("Original EEA2026 row count changed")
    return records,dict(qc)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-source-quartets",type=Path,required=True)
    p.add_argument("--original-32site-numeric-QC",type=Path,required=True)
    p.add_argument("--tempdir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    qu,qual=read_quartets(a.original_source_quartets,a.original_32site_numeric_QC)
    a.tempdir.mkdir(parents=True,exist_ok=True)
    file=a.tempdir/"original_EEA2026_WISE6_2022_sampling_dates.zip"
    original_receipt=get_original(file)
    rec,qc=scan_metadata(file,set(qual))
    for id,q in qual.items():
        if rec[id]["sample_rows"]!=q["original_2022_all_sample_rows"]:
            raise ValueError("Original source and numeric QC disagree on raw 2022 sample rows for station "+id)
    site_source={}
    for id,r in rec.items():
        d=sorted(r["date_days"])
        site_source[id]={
          "original_public_WISE6_country_scheme_site_sha256":id,
          "2022_original_nonvalue_sample_rows":r["sample_rows"],
          "2022_original_distinct_sample_dates":d,
          "2022_original_distinct_months":sorted(set(x[:7] for x in d)),
          "2022_original_nonpersonal_distinct_sampleIDs_n":len(r["sample_ids"]),
          "original_WISE6_resultObservedValue_or_LOQ_not_read":True,
        }
    stats=[];four_day=[];four_month=[]
    for i,q in qu.items():
        arms=q["four_arm_original_WISE_country_scheme_site_hashes"]
        day={k:set(rec[x]["date_days"]) for k,x in arms.items()}
        mon={k:set(d[:7] for d in val) for k,val in day.items()}
        d_src=day["source_upstream"]&day["source_downstream"]
        d_ctl=day["registered_zero_upstream"]&day["registered_zero_downstream"]
        d_all=d_src&d_ctl
        m_src=mon["source_upstream"]&mon["source_downstream"]
        m_ctl=mon["registered_zero_upstream"]&mon["registered_zero_downstream"]
        m_all=m_src&m_ctl
        eligible_day=len(d_all)>=MIN_ALLFOUR_DAY
        eligible_month=len(m_all)>=MIN_ALLFOUR_MONTH
        if eligible_day:four_day.append(q)
        if eligible_month:four_month.append(q)
        stats.append({
          **q,
          "four_arm_same_full_sampling_dates":sorted(d_all),
          "four_arm_same_YYYY_MM_months":sorted(m_all),
          "source_pair_same_full_sampling_dates":sorted(d_src),
          "registered_source_zero_pair_same_full_sampling_dates":sorted(d_ctl),
          "source_pair_same_YYYY_MM_months":sorted(m_src),
          "registered_source_zero_pair_same_YYYY_MM_months":sorted(m_ctl),
          "strict_same_day_ge3_full_fourstation_support":eligible_day,
          "strict_same_month_ge3_full_fourstation_support":eligible_month,
          "numerical_TP_concentration_was_not_read_by_this_script":True,
        })
        print("GESUI_2022_DE_SOURCEONLY_TP_FOUR_ARM_SAMEDAY_MONTH_AUDIT",json.dumps({
            "quartet":i,"source_country2deg_block":q["source_country2degree_upstream_block"],
            "src_same_dates_n":len(d_src),"control_same_dates_n":len(d_ctl),
            "all4_same_days_n":len(d_all),"all4_same_months_n":len(m_all),
            "source_interposed_original_P_t_year":q["original_source_interposed_P_tonnes_per_year"]},
            sort_keys=True),flush=True)
    blocksday={q["source_country2degree_upstream_block"] for q in four_day}
    blocksmonth={q["source_country2degree_upstream_block"] for q in four_month}
    daygate=len(four_day)>=MIN_QUARTETS_DAY and len(blocksday)>=MIN_BLOCKS_DAY
    monthgate=len(four_month)>=MIN_QUARTETS_MONTH and len(blocksmonth)>=MIN_BLOCKS_MONTH
    out={
      "schema":"gesui_2022_original_EEA_WISE6_Germany_8_source_zero_TP_quartet_exact_day_and_month_preoutcome_source_gate_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "POSTHOC_after_original_7_quartet_annual_median_positive_was_seen":True,
      "source_original_Git_blob_quartet_SHA1":SOURCE_QUARTET_BLOB,
      "source_original_32site_numeric_QC_SHA256":SOURCE_NUMERIC32_SHA,
      "original_2026_EEA_WISE6_full_unmodified_original_zip_sha256":original_receipt["sha256"],
      "original_EEA_water_source_download":original_receipt,
      "2022_TP_CAS_code":TP,"unit":"mg{P}/L",
      "exact_site_source_32_original_identifiers_hashed_and_calendar":site_source,
      "original_2022_station_source_row_counts":qc,
      "same_day_four_original_sites_min_days":MIN_ALLFOUR_DAY,
      "same_day_min_quartets":MIN_QUARTETS_DAY,
      "same_day_min_source2deg_blocks":MIN_BLOCKS_DAY,
      "same_month_four_original_sites_min_months":MIN_ALLFOUR_MONTH,
      "same_month_min_quartets":MIN_QUARTETS_MONTH,
      "same_month_min_source2deg_blocks":MIN_BLOCKS_MONTH,
      "all_original_8_no_replacement_quartets_and_date_coverage":stats,
      "source_quartets_ge3_same_day_n":len(four_day),
      "source_quartets_ge3_same_day_geoblocks_n":len(blocksday),
      "source_quartets_ge3_same_month_n":len(four_month),
      "source_quartets_ge3_same_month_geoblocks_n":len(blocksmonth),
      "source_exact_fourstation_calendar_date_gate_PASS":daygate,
      "source_same_calendar_month_fourstation_gate_PASS":monthgate,
      "STOP_before_new_numeric_monthly_sourceTP_if_both_gates_fail":not(daygate or monthgate),
      "NO_original_TP_numeric_resultObservedValue_or_LOQ_values_read":True,
      "sample_dates_are_not_actual_WWTP_treatment_or_effluent_doses":True,
      "source_zero_quartets_not_untreated_sites":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_DE_2022_SOURCELOCKED_32STATION_FOUR_ARM_SAMEDAY_CALENDAR_GATE",json.dumps({
       "n_8_quartets":len(qu),"source_ge3_day_quartets":len(four_day),
       "same_day_country2deg_blocks":len(blocksday),"daygate":daygate,
       "source_ge3_month_quartets":len(four_month),
       "same_month_2deg_blocks":len(blocksmonth),"monthgate":monthgate,
       "no_numeric_TP_values_loaded":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
