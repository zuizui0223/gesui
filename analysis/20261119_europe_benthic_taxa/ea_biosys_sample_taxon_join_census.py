#!/usr/bin/env python3
"""Official EA Biosys complete SAMPLE/ANALYSIS/TAXON source-join coverage census.

Streams original 62.4MB archive from Environment Agency. This is a
source-quality census only: no taxon/environment slopes, NO wastewater
exposure, NO taxon trait/imputation, and no raw taxon/source rows rehosted.
"""
from __future__ import annotations

import argparse, collections, csv, datetime, hashlib, io, json, tempfile
import urllib.request, zipfile
from pathlib import Path

URL="https://environment.data.gov.uk/ecology/explorer/downloads/INV_OPEN_DATA.zip"
REFSIZE=62373479
SOURCE_REQUIRED={
    "INV_OPEN_DATA_SITE.csv":{"SITE_ID","CATCHMENT","WFD_WATERBODY_ID"},
    "INV_OPEN_DATA_METRICS.csv":{"SITE_ID","SAMPLE_ID","ANALYSIS_ID","SAMPLE_DATE","BMWP_ASPT"},
    "INV_OPEN_DATA_TAXA.csv":{"SITE_ID","SAMPLE_ID","ANALYSIS_ID","TAXON_LIST_ITEM_KEY","TOTAL_ABUNDANCE","ABUNDANCE_CATEGORY"}
}
MAX_ROWS=25000000

def download_source(filename):
    hash=hashlib.sha256();n=0
    req=urllib.request.Request(URL,headers={"User-Agent":"gesui-EA-Biosys-source-census/20261119"})
    with urllib.request.urlopen(req,timeout=65) as r,filename.open("wb") as f:
        if r.status!=200:raise RuntimeError("Official original source HTTP not 200")
        while True:
            chunk=r.read(1024*1024)
            if not chunk:break
            n+=len(chunk)
            if n>100*1024*1024:raise RuntimeError("Unexpected increased EA original zip size")
            f.write(chunk);hash.update(chunk)
    if n!=REFSIZE:raise RuntimeError("Original EA ZIP size drift, check version before analysis")
    return hash.hexdigest()

def records(z,name):
    with z.open(name) as raw:
        txt=io.TextIOWrapper(raw,encoding="utf-8-sig",errors="replace",newline="")
        reader=csv.DictReader(txt)
        if not SOURCE_REQUIRED[name].issubset(set(reader.fieldnames or [])):
            raise RuntimeError("Missing original join keys: "+name)
        for row in reader:yield row

def year_of(s):
    raw=(s or "").strip()
    if not raw:return None
    for fmt in ("%Y-%m-%d","%d/%m/%Y","%Y/%m/%d","%d-%b-%Y","%d/%m/%y","%Y-%m-%d %H:%M:%S"):
        try:return datetime.datetime.strptime(raw,fmt).year
        except ValueError:pass
    try:return datetime.datetime.fromisoformat(raw.replace("Z","+00:00")).year
    except ValueError:return None

def make(args):
    result={"schema":"GESUI_EA_ORIGINAL_BIOSYS_SITE_SAMPLE_TAXON_SOURCE_CENSUS_V1",
            "status":"SOURCE_STOP","analysis_outcome_unopened":True,
            "official_source_URL":URL,"no_pollution_source_or_trait_score_assigned":True,
            "no_WWTP_exposure_joined":True}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_Biosys_raw_") as temp:
            f=Path(temp)/"original.zip"
            digest=download_source(f)
            result["original_zip_SHA256"]=digest
            with zipfile.ZipFile(f) as z:
                if not set(SOURCE_REQUIRED).issubset(set(z.namelist())):
                    raise RuntimeError("Original official three-relational-table source absent")
                sites={}
                site_duplicate=0
                site_conflicts=0
                for row in records(z,"INV_OPEN_DATA_SITE.csv"):
                    key=row["SITE_ID"].strip()
                    if not key:continue
                    detail={"catchment":row["CATCHMENT"].strip(),
                            "wb":row["WFD_WATERBODY_ID"].strip(),
                            "agency":row.get("AGENCY_AREA","").strip()}
                    if key in sites:
                        site_duplicate+=1
                        if sites[key]!=detail:site_conflicts+=1
                        continue
                    sites[key]=detail
                source_count={
                  "site_records_unique_SITE_ID":len(sites),
                  "site_duplicate_ID_extra_rows":site_duplicate,
                  "site_duplicate_conflicting_metadata":site_conflicts,
                  "site_with_catchment":sum(bool(d["catchment"]) for d in sites.values()),
                  "site_with_WFD_waterbody":sum(bool(d["wb"]) for d in sites.values()),
                  "distinct_catchment":len({d["catchment"] for d in sites.values() if d["catchment"]}),
                }
                metrics={}
                metrics_rows=0
                duplicate_analysis=0
                year_counter=collections.Counter()
                raw_dates=collections.Counter()
                sites_by_year=collections.defaultdict(set)
                methods=collections.Counter()
                sample_types=collections.Counter()
                sample_reasons=collections.Counter()
                site_missing=0
                for row in records(z,"INV_OPEN_DATA_METRICS.csv"):
                    metrics_rows+=1
                    aid=row["ANALYSIS_ID"].strip()
                    sid=row["SITE_ID"].strip()
                    yr=year_of(row["SAMPLE_DATE"])
                    if not aid:continue
                    if aid in metrics:
                        duplicate_analysis+=1
                        continue
                    if sid not in sites:site_missing+=1
                    metrics[aid]=(sid,row["SAMPLE_ID"].strip(),yr)
                    year_counter[str(yr) if yr else "UNKNOWN"]+=1
                    if yr is None:raw_dates[(row["SAMPLE_DATE"] or "")[:24]]+=1
                    if yr is not None:sites_by_year[str(yr)].add(sid)
                    methods[row["SAMPLE_METHOD"].strip()]+=1
                    sample_types[row["SAMPLE_TYPE"].strip()]+=1
                    sample_reasons[row["SAMPLE_REASON"].strip()]+=1
                    if metrics_rows>MAX_ROWS:raise RuntimeError("Unexpected >25m metric rows")
                result["site"] = source_count
                result["metrics"]={
                    "row_count":metrics_rows,"unique_original_ANALYSIS_ID":len(metrics),
                    "duplicate_ANALYSIS_ID_extra_rows":duplicate_analysis,
                    "unique_analysis_missing_SITE":site_missing,
                    "years_metric_counts":dict(sorted(year_counter.items())),
                    "years_unique_site_counts":{k:len(v) for k,v in sorted(sites_by_year.items())},
                    "sample_method_counts_top15":methods.most_common(15),
                    "sample_type_counts_top15":sample_types.most_common(15),
                    "sample_reason_counts_top15":sample_reasons.most_common(15),
                    "undecodable_dates_n":sum(raw_dates.values()),
                    "undecodable_date_examples":list(raw_dates.keys())[:5],
                }
                rows=0;bad_a=0;bad_site=0;bad_sample=0;no_taxon=0;count_number=0;count_not_numeric=0;category_any=0
                by_year=collections.Counter()
                by_year_actual_count=collections.Counter()
                observed_analyses=set()
                original_taxon_ids=set()
                for row in records(z,"INV_OPEN_DATA_TAXA.csv"):
                    rows+=1
                    aid=row["ANALYSIS_ID"].strip()
                    sid=row["SITE_ID"].strip()
                    sample=row["SAMPLE_ID"].strip()
                    taxon=row["TAXON_LIST_ITEM_KEY"].strip()
                    if not taxon:no_taxon+=1
                    else:original_taxon_ids.add(taxon)
                    if aid not in metrics:
                        bad_a+=1
                        continue
                    m_site,m_sample,yr=metrics[aid]
                    if m_site!=sid:bad_site+=1
                    if m_sample!=sample:bad_sample+=1
                    observed_analyses.add(aid)
                    by_year[str(yr) if yr else "UNKNOWN"]+=1
                    c=row["TOTAL_ABUNDANCE"].strip()
                    try:
                        val=float(c)
                        if val>0 and val==val:
                            count_number+=1
                            by_year_actual_count[str(yr) if yr else "UNKNOWN"]+=1
                        else:count_not_numeric+=1
                    except (ValueError,OverflowError):count_not_numeric+=1
                    if row["ABUNDANCE_CATEGORY"].strip():category_any+=1
                    if rows>MAX_ROWS:raise RuntimeError("Unexpected >25m taxon rows")
                result["taxa"]={
                    "row_count":rows,
                    "matched_analysis_id_rows":rows-bad_a,
                    "missing_analysis_id_rows":bad_a,
                    "analysis_to_site_conflicts":bad_site,
                    "analysis_to_sample_conflicts":bad_sample,
                    "unique_analysis_with_taxon_rows":len(observed_analyses),
                    "unique_original_taxon_key_in_taxa":len(original_taxon_ids),
                    "missing_taxon_ID":no_taxon,
                    "positive_numeric_TOTAL_ABUNDANCE_rows":count_number,
                    "nonnumeric_or_zero_TOTAL_ABUNDANCE_rows":count_not_numeric,
                    "nonempty_ABUNDANCE_CATEGORY_rows":category_any,
                    "year_taxa_rows":dict(sorted(by_year.items())),
                    "year_taxa_positive_numeric_rows":dict(sorted(by_year_actual_count.items())),
                }
                result["status"]="SOURCE_RELATIONAL_GRANULARITY_AUDITED"
                result["biological_model_authorized"]=False
                result["source_go_criteria"]={
                    "analysis_to_site_sample_join_100pct":bad_a==0 and bad_site==0 and bad_sample==0,
                    "no_duplicate_analysis_ID":duplicate_analysis==0,
                    "no_conflicting_duplicate_SITE_metadata":site_conflicts==0,
                    "count_CWM_can_be_considered_if_numeric_coverage":count_number/(rows-bad_a) if rows>bad_a else None,
                    "must_independently_map_taxon_names_and_score_coverage":True,
                    "must_source_real_WWTP_or_ambient_NP_with_same_sample":True,
                    "not_a_demographic_occupancy_or_replacement_test":True
                }
    except Exception as e:
        result.update({"status":"SOURCE_STOP","reason":type(e).__name__+": "+str(e)[:420]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_EA_BIOSYS_SOURCE_JOIN_CENSUS="+json.dumps({
        "status":result["status"],"reason":result.get("reason"),
        "site":result.get("site"),
        "metrics":result.get("metrics"),
        "taxa":result.get("taxa"),
        "source_go":result.get("source_go_criteria")
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
