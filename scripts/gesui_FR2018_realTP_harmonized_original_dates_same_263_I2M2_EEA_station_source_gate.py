#!/usr/bin/env python3
"""SOURCE ONLY original Peters 2018 actually measured river TP datewise records
for same WQ IDs as original 2022 true TP–I2M2 EQR, with 2018 exact
same physical EEA station I2M2 source IDs independently preselected.

Official original Figshare v5 TP.2 harmonized ZIP:
70896288, 25003722 B, publisher supplied MD5 b596dbe...
Does NOT read numeric 2018/2022 TP median or 2018/2022 I2M2 EQR.
Only tests positive measurement validity, distinct dates and sampling months.
"""
import collections,csv,datetime,hashlib,io,json,math,os,pathlib,time,urllib.request,zipfile
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_FR2018_2022_TP_datewise_source"
ROOT.mkdir(parents=True,exist_ok=True)
SOURCE_SHA="e5beb11746d36d9f73b5772371589b68b11eed43b087a828ab78fed2255879dd"
ARTICLE="https://api.figshare.com/v2/articles/30546371"
SOURCE_DOI="10.25400/lincolnuninz.30546371.v5"
ZIP_ID=70896288
SOURCE_ZIP_MD5="b596dbe9833c49426147a99f10f3fe31"
SOURCE_ZIP_SIZE=25003722
src=pathlib.Path(os.environ["ORIGINAL_SOURCE_FR2022_273_TP_EQR_2018_BIOLOGY_QUALITY"])
raw=src.read_bytes()
if hashlib.sha256(raw).hexdigest()!=SOURCE_SHA:
    raise RuntimeError("Original source-only 263 genuine identical 2018/2022 EEA station IDs manifest changed")
report=json.loads(raw)
paired=report["source_2018_2022_EEA_same_station_I2M2_biology_geography_only"]
original2022=report["source_station_pair_registrations_2022"]
if len(original2022)!=273 or len(paired)!=263 or \
    len({z["Peters_2022_WQ_ID"] for z in paired})!=263 or \
    len({z["EEA2022_site_key"] for z in paired})!=263 or \
    len({z["original_HYBAS_ID_L07"] for z in paired})!=73 or \
    len({z["original_H6_parent"] for z in paired})!=24:
    raise RuntimeError("Original source physical 2018 EEA I2M2 station continuity cohort changed")
if any(not z["original_2018_2022_procedure_equal_I2M2"] for z in paired):
    raise RuntimeError("Original 2018 EEA biological methods differ from 2022")
# Obtain official publisher API response and source file checksum, no
# substituted modelled TP or metadata first-last sampling-date approximations.
req=urllib.request.Request(ARTICLE,headers={"Accept":"application/json","User-Agent":"gesui-original-2018-TP-datewise-source-only/1"})
with urllib.request.urlopen(req,timeout=50) as response:body=response.read(350000)
article=json.loads(body)
if article.get("doi")!=SOURCE_DOI:
    raise RuntimeError("Official Figshare original nutrient article DOI version differs")
idx=[x for x in article.get("files",[]) if x.get("name")=="TP.2_-_Harmonised.zip"]
if len(idx)!=1 or idx[0]["id"]!=ZIP_ID or idx[0].get("supplied_md5")!=SOURCE_ZIP_MD5 or idx[0].get("size")!=SOURCE_ZIP_SIZE:
    raise RuntimeError("Original publisher 2018 datewise TP archive original file identity/MD5/size mismatch")
link=idx[0]["download_url"]
req=urllib.request.Request(link,headers={"Accept":"application/zip,application/octet-stream,*/*","User-Agent":"gesui-original-Peters-2018-TP-physically-paired-source/1"})
buf=bytearray()
with urllib.request.urlopen(req,timeout=120) as response:
    while True:
        block=response.read(131072)
        if not block:break
        buf.extend(block)
        if len(buf)>27000000:raise RuntimeError("Original measured TP archive unexpectedly oversized; source STOP")
raw=bytes(buf)
if len(raw)!=SOURCE_ZIP_SIZE or hashlib.md5(raw).hexdigest()!=SOURCE_ZIP_MD5:
    raise RuntimeError("Original measured TP publication bytes failed publisher MD5 or length check")
orig_zip_SHA256=hashlib.sha256(raw).hexdigest()
arc=zipfile.ZipFile(io.BytesIO(raw))
members=set(arc.namelist())
if len(members)<1000:raise RuntimeError("Original TP source ZIP has unexpected station file structure")
# Dates normalized calendar-year UTC YYYY-MM-DD prefixes of source original
# harmonized measured value; source quality only boolean valid positive.
date_stats=collections.Counter()
rows=[]
for item in paired:
    k=item["Peters_2022_WQ_ID"]
    file="TP.2_-_Harmonised/"+k+".csv"
    if file not in members:
        date_stats["source_real_WQ_ID_file_missing"]+=1
        rows.append({**item,"source_2018_2019_wq_file_present":False,
                     "source_2018_eligible":False,
                     "source_2022_eligible":False})
        continue
    b=arc.read(file)
    reader=csv.DictReader(io.StringIO(b.decode("utf-8-sig"),newline=""))
    if not {"date","value"}.issubset(reader.fieldnames or ()):
        raise RuntimeError("Original Peters v5 TP.2 datewise CSV source fields not date/value")
    byyear=collections.defaultdict(list)
    invalid=0
    for q in reader:
        t=str(q.get("date") or "").strip()[:10]
        try:
            d=datetime.date.fromisoformat(t)
            val=float(q["value"])
        except (ValueError,TypeError):
            invalid+=1
            continue
        if not math.isfinite(val) or val<=0:
            invalid+=1;continue
        if d.year in (2018,2022):
            byyear[d.year].append(d)
    elig={}
    for year in (2018,2022):
        arr=byyear[year];months={d.month for d in arr};days=set(arr)
        elig[year]=len(arr)>=4 and len(months)>=3 and len(days)>=4
    if elig[2018]:date_stats["original_WQ_source_2018_supports_4_valid_3months"]+=1
    if elig[2022]:date_stats["original_WQ_source_2022_supports_4_valid_3months"]+=1
    if elig[2018] and elig[2022]:date_stats["original_WQ_source_both_years_eligible"]+=1
    if len(byyear[2018])==0:date_stats["original_WQ_source_no_positive_2018_measurements"]+=1
    rows.append({
      "original_WQ_ID":k,
      "EEA2022_site_key":item["EEA2022_site_key"],
      "original_EEA2018_site_key":item["original_EEA2018_site_key"],
      "HYBAS_ID_L07":item["original_HYBAS_ID_L07"],
      "HYBAS_ID_L06_parent":item["original_H6_parent"],
      "original_TP2_WQ_file_present":True,
      "n_2018_original_positive_measured_TP_dates":len(byyear[2018]),
      "n_2022_original_positive_measured_TP_dates":len(byyear[2022]),
      "n_2018_distinct_positive_measured_TP_months":len({d.month for d in byyear[2018]}),
      "n_2022_distinct_positive_measured_TP_months":len({d.month for d in byyear[2022]}),
      "n_2018_distinct_original_sample_dates":len(set(byyear[2018])),
      "n_2022_distinct_original_sample_dates":len(set(byyear[2022])),
      "n_original_invalid_dates_or_nonpositive_TP":invalid,
      "source_2018_eligible":elig[2018],
      "source_2022_eligible":elig[2022],
      "original_2018_2022_numerical_TP_values_NOT_READ_into_source_decision":True,
      "original_2018_2022_EEA_EQR_values_NOT_READ":True,
      "source_2018_2022_EEA_method_same_I2M2":True
    })
if len(rows)!=263:raise RuntimeError("Original 2018 same physical EEA station list unexpectedly changed")
valid=[r for r in rows if r["source_2018_eligible"] and r["source_2022_eligible"]]
ids={x["HYBAS_ID_L07"] for x in valid}
parents={x["HYBAS_ID_L06_parent"] for x in valid}
gate=(len(valid)>=100 and len(ids)>=30 and len(parents)>=10)
qa={
 "schema":"GESUI_ORIGINAL_2018_2022_MEASURED_PETERS_TP_SOURCE_DATE_263_EEA_SAME_STATION_74H7_v1",
 "original_273_EEA_WQ_2018_biological_source_QA_SHA256":SOURCE_SHA,
 "original_Figshare_article":ARTICLE,
 "original_published_DOI":SOURCE_DOI,
 "publisher_true_original_fileID":ZIP_ID,
 "original_published_measuredTP_zip_size_bytes":SOURCE_ZIP_SIZE,
 "official_published_ZIP_md5_verified":SOURCE_ZIP_MD5,
 "measuredTP_zip_original_bytes_SHA256":orig_zip_SHA256,
 "original_TP_zip_all_station_members":len(members),
 "original_2022_TP_EEA_source_273_stations":len(original2022),
 "original_2018_same_I2M2_benthic_repeated_EEA_sites":len(paired),
 "source_WQ_validity_year_counts":dict(date_stats),
 "n_same_WQ_and_EEA_benthic_sites_with_2018_and_2022_measured_TP_sample_support":len(valid),
 "n_source_2018_2022_eligible_H7":len(ids),
 "n_source_2018_2022_eligible_H6":len(parents),
 "n_original_distinct_WQ_ids_at_two_years":len({x["original_WQ_ID"] for x in valid}),
 "predeclared_longitudinal_measured_source_gate_minima":{
    "same_physical_field_sites":100,"independent_native_H7":30,"native_H6_parents":10},
 "source_measured_2018_and_2022_datewise_TP_I2M2_longitudinal_gate_PASS":gate,
 "source_GATE_STATUS":"ELIGIBLE_TO_OPEN_REAL_2018_TP_NUMERIC_AND_EQR_CHANGE_POSTHOC" if gate else "STOP_ORIGINAL_2018_TP_DATE_SUPPORT_UNDERSUPPORTED",
 "no_actual_numerical_TP_2018_or_2022_median_EEA_2018_2022_EQR_combined_at_source_stage":True,
 "cannot_infer_2018_2022_WWTP_source_P_causation_from_2018_2022_ambient_TP":True,
 "original_source_positive_2022_WWTP_P_index_only_four_sites_in_273":True,
 "no_observed_terrestrial_ecosystem_service_function":True
}
d={"schema":qa["schema"],"original_Figshare_TP2_source_ZIP_SHA256":orig_zip_SHA256,
 "original_source_only_exact_same_EEA2018_2022_biological_site_and_annual_measured_TP_support":rows,
 "source_can_read_original_actual_TP_medians_and_EEA_EQR_only_after_gate":gate}
encoded=(json.dumps(d,sort_keys=True,separators=(",",":"))+"\n").encode()
(ROOT/"FR2018_2022_TRUE_TP_I2M2_SOURCE_ONLY_STATION_DATES_ELIGIBLE.json").write_bytes(encoded)
qa["source_2018_2022_actual_TP_sample_eligibility_SHA256"]=hashlib.sha256(encoded).hexdigest()
(ROOT/"FR2018_2022_SOURCE_DATEWISE_TP_263_I2M2_SITE_SUPPORT_QA.json").write_text(json.dumps(qa,sort_keys=True,indent=2))
print("FR2018_TP_TRUE_MEASUREMENT_SOURCE_GATE_START",flush=True)
print(json.dumps(qa,sort_keys=True),flush=True)
print("FR2018_TP_TRUE_MEASUREMENT_SOURCE_GATE_END",flush=True)
# A source-insufficient result is an analyzable negative data-quality result;
# do not overwrite or source-select alternative 2018 years or stations.
