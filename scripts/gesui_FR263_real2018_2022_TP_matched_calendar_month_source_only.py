#!/usr/bin/env python3
"""Source-only actual 2018 versus 2022 TP measured month overlap in EXACT
263 source-selected Peters WQ/EEA I2M2 sites. Read only dates and whether a
numeric TP observation was positive, NOT TP magnitudes/EQR outcomes.

Original Figshare TP.2 harmonized ZIP (v5) SHA + publisher MD5 verified.
This is posthoc after previous original raw-annual spatial positive/temporal
null, but MONTH SUPPORT is frozen before opening newly balanced TP effects.
"""
import collections,csv,datetime,hashlib,io,json,math,os,pathlib,urllib.request,zipfile
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2018_22_TP_MONTH_SOURCE_ONLY"
ROOT.mkdir(parents=True,exist_ok=True)
PRE="2190ac905326f4f2f4f278eb0debb4a7bfda7c4c71347d4941b54cc19d3de264"
ZIP_SHA="26f2af11f4972987ee6f7c3f497ad902820f145851f71919d4f45b589e1d809e"
ZIP_MD5="b596dbe9833c49426147a99f10f3fe31"
ZIP_SIZE=25003722;F_ID=70896288
prepath=pathlib.Path(os.environ["ORIGINAL_DATEONLY_2018_22_263_SAME_SITE_SOURCE"])
pbytes=prepath.read_bytes()
if hashlib.sha256(pbytes).hexdigest()!=PRE:raise RuntimeError("Frozen original 263 WQ/EEA source records SHA changed")
pre=json.loads(pbytes)
if pre["schema"]!="GESUI_ORIGINAL_2018_2022_MEASURED_PETERS_TP_SOURCE_DATE_263_EEA_SAME_STATION_74H7_v1":
    raise RuntimeError("Original exact 263 2018/22 eligible-source type changed")
ss=pre["original_source_only_exact_same_EEA2018_2022_biological_site_and_annual_measured_TP_support"]
if len(ss)!=263 or len({s["original_WQ_ID"] for s in ss})!=263:
    raise RuntimeError("Previous source stage 263 genuine WQ IDs changed")
if len({s["HYBAS_ID_L07"] for s in ss})!=73 or len({s["HYBAS_ID_L06_parent"] for s in ss})!=24:
    raise RuntimeError("Original source 73 H7 or 24 H6 geography changed")
url=f"https://ndownloader.figshare.com/files/{F_ID}"
request=urllib.request.Request(url,headers={"User-Agent":"gesui-2018-2022-source-only-matched-calendar-months/1","Accept":"application/zip,application/octet-stream,*/*"})
data=bytearray()
with urllib.request.urlopen(request,timeout=110) as f:
    while True:
        chunk=f.read(131072)
        if not chunk:break
        data.extend(chunk)
        if len(data)>27000000:raise RuntimeError("Official TP.2 archive exceeded fixed source maximum")
raw=bytes(data)
if len(raw)!=ZIP_SIZE or hashlib.sha256(raw).hexdigest()!=ZIP_SHA or hashlib.md5(raw).hexdigest()!=ZIP_MD5:
    raise RuntimeError("Original TP.2 published source MD5/SHA/bytes mismatch")
archive=zipfile.ZipFile(io.BytesIO(raw))
names=set(archive.namelist())
out=[];faults=collections.Counter()
for source in ss:
    wq=source["original_WQ_ID"]
    fname="TP.2_-_Harmonised/"+wq+".csv"
    if fname not in names:raise RuntimeError("Source-eligible original site has no official WQ file")
    rr=csv.DictReader(io.StringIO(archive.read(fname).decode("utf-8-sig"),newline=""))
    if not {"date","value"}.issubset(rr.fieldnames or []):
        raise RuntimeError("Original author TP.2 date/value fields missing")
    dct=collections.defaultdict(list)
    invalid=0;dup=collections.Counter()
    for row in rr:
        try:
            day=datetime.date.fromisoformat(str(row["date"]).strip()[:10])
            val=float(row["value"])
        except (ValueError,TypeError):
            invalid+=1;continue
        if not math.isfinite(val) or val<=0:
            invalid+=1;continue
        if day.year not in (2018,2022):continue
        # no raw TP value output or source-selection dependence on magnitude
        dct[day.year].append(day)
        dup[(day.year,day)]+=1
    months={yr:{d.month for d in dct[yr]} for yr in (2018,2022)}
    common=sorted(months[2018]&months[2022])
    if len(dct[2018])<4 or len(dct[2022])<4 or min(len(months[2018]),len(months[2022]))<3:
        raise RuntimeError("Original source-frozen annual measurement eligibility contradicted")
    can=len(common)>=3
    if can:faults["eligible_station_month_overlap_at_least3"]+=1
    faults["source_original_positive_duplicate_date_records"]+=sum(v-1 for v in dup.values() if v>1)
    if not common:faults["no_shared_calendar_month"]+=1
    out.append({
      "WQ_ID":wq,
      "original_EEA_identical_site_key":source["EEA2022_site_key"],
      "original_HYBAS_ID_L07":source["HYBAS_ID_L07"],
      "original_H6_parent":source["HYBAS_ID_L06_parent"],
      "source_2018_original_positive_days_n":len(set(dct[2018])),
      "source_2022_original_positive_days_n":len(set(dct[2022])),
      "source_2018_original_positive_readings_n":len(dct[2018]),
      "source_2022_original_positive_readings_n":len(dct[2022]),
      "source_2018_actual_sample_months":sorted(months[2018]),
      "source_2022_actual_sample_months":sorted(months[2022]),
      "source_common_original_calendar_months":common,
      "source_common_month_count":len(common),
      "source_2018_n_samples_in_common_months":sum(d.month in common for d in dct[2018]),
      "source_2022_n_samples_in_common_months":sum(d.month in common for d in dct[2022]),
      "source_first_day_2018":str(min(dct[2018])),
      "source_last_day_2018":str(max(dct[2018])),
      "source_first_day_2022":str(min(dct[2022])),
      "source_last_day_2022":str(max(dct[2022])),
      "source_eligible_at_least3_common_calendar_months":can,
      "original_TP_numeric_magnitudes_NOT_OPENED_FOR_MONTH_ELIGIBILITY":True,
      "original_EEA_I2M2_numeric_change_NOT_OPENED":True
    })
if len(out)!=263:raise RuntimeError("263 original source WQ IDs not preserved")
eligible=[r for r in out if r["source_eligible_at_least3_common_calendar_months"]]
h7={r["original_HYBAS_ID_L07"] for r in eligible}
h6={r["original_H6_parent"] for r in eligible}
gate=bool(len(eligible)>=100 and len(h7)>=30 and len(h6)>=10)
hist=collections.Counter(str(r["source_common_month_count"]) for r in out)
month2018=collections.Counter(m for r in out for m in r["source_2018_actual_sample_months"])
month2022=collections.Counter(m for r in out for m in r["source_2022_actual_sample_months"])
both=collections.Counter(m for r in out for m in r["source_common_original_calendar_months"])
info={
  "schema":"GESUI_FR263_ORIGINAL_2018_2022_TP_COMMON_CALENDAR_MONTH_ELIGIBILITY_SOURCE_ONLY_v1",
  "source_original_2018_2022_TP_I2M2_station_263_pre_dateonly_SHA256":PRE,
  "official_TP2_ZIP_original_SHA256":ZIP_SHA,
  "official_TP2_ZIP_publisher_verified_MD5":ZIP_MD5,
  "all_original_WQ_EEA_same_2018_2022_sites":len(out),
  "source_2018_2022_years": [2018,2022],
  "source_original_WQ_month_overlap_histogram":dict(sorted(hist.items(),key=lambda x:int(x[0]))),
  "source_original_2018_sampled_month_frequency_stationcounts":dict(month2018),
  "source_original_2022_sampled_month_frequency_stationcounts":dict(month2022),
  "source_original_months_seen_BOTH_years_by_sites":dict(both),
  "source_original_reading_quality_counts":dict(faults),
  "source_monthmatched_primary_eligible_sites_n":len(eligible),
  "source_monthmatched_primary_eligible_unique_H7":len(h7),
  "source_monthmatched_primary_eligible_H6_parents":len(h6),
  "source_fixed_pre_new_outcome_minimums":{"stations":100,"H7":30,"H6_parents":10,"shared_calendar_months":3},
  "source_monthmatched_new_biological_model_GATE_PASS":gate,
  "source_chemical_magnitude_and_EQR_year_change_NOT_used_to_choose_matched_months":True,
  "not_necessarily_same_day_or_flow_condition_across_years":True,
  "source_2018_and2022_2018_2022_EEA_same_method_station_ID_preserved":True
}
dataset={
  "schema":info["schema"],
  "original_263_physically_identical_2018_2022_WQ_and_EEA_source_matched_months":out,
  "original_source_monthmatched_new_model_gate_pass":gate,
  "no_numeric_TP_medians_or_EEA_EQR_values_in_this_source":True
}
by=(json.dumps(dataset,sort_keys=True,separators=(",",":"))+"\n").encode()
(ROOT/"FR263_ORIGINAL_TP_I2M2_SAME_STATION_COMMON_MONTHS_SOURCE_ONLY.json").write_bytes(by)
info["source_2018_2022_month_overlap_fixed_SHA256"]=hashlib.sha256(by).hexdigest()
(ROOT/"FR263_ORIGINAL_TP_I2M2_SAME_STATION_COMMON_MONTHS_QA.json").write_text(json.dumps(info,sort_keys=True,indent=2))
print("FR263_TP_ORIGINAL_COMMON_MONTHS_SOURCE_GATE_START",flush=True)
print(json.dumps(info,sort_keys=True),flush=True)
print("FR263_TP_ORIGINAL_COMMON_MONTHS_SOURCE_GATE_END",flush=True)
# Even if false, preserve QA artifact; STOP new effects in downstream job.
