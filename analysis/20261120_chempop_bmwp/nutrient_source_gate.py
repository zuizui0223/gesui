#!/usr/bin/env python3
"""Source-only original ChemPop v2 N/P analyte×macroinvertebrate date coverage.

Lists all N/P candidate analyte summary files by publisher names WITHOUT
looking at correlations or biological outcomes. Reads only published source
water-quality count/LOQ coverage, matched to original scoring-eligible
2010-18 biological site×sample-date keys, prior to any nutrient model.
Never treats mean LOQ-based window-statistic as an exactly simultaneous
water-chemistry field measurement. Publishes no original site/sample rows.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
import re
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261119_europe_benthic_taxa"))
import chempop_v2_source_census as source
import chempop_v2_bmwp_filter as f

EXPECTED_SOURCE_SAMPLES=7836
EXPECTED_SOURCE_SITES=707
EXPECTED_SITEYEARS=3884
EXPECTED_SCORED=92
WQ_PATTERN=re.compile(
    r"^data/cp_(an|mi|ne|nw|so|sw|th)_minv_([a-z0-9_]+)stats_([0-9]{8})-([0-9]{8})[.]csv$",
    re.I
)
NUTRIENT_TOKEN=re.compile(r"phosph|phos|nitra|nitrit|ammon|nitrogen",re.I)
PREOUTCOME_GATES={
   "min_matched_distinct_bio_site_date":500,
   "min_matched_distinct_bio_site":150,
   "min_regions_with_valid_source_sites":4,
   "per_valid_region_min_sites":10,
   "require_numeric_mean_LOQ2_and_source_WQ_COUNT_positive":True,
   "paired_biology_site_dates_are_original_ChemPop_benthic_samples":True,
}

def candidate_names(files):
    out=[]
    for path in files:
        m=WQ_PATTERN.fullmatch(path)
        if not m:continue
        region,analyte,start,end=m.groups()
        if not NUTRIENT_TOKEN.search(analyte):continue
        out.append({"path":path,"region":region.lower(),
                    "analyte":analyte.lower(),"range_start":start,"range_end":end})
    return sorted(out,key=lambda x:(x["analyte"],x["region"],x["path"]))

def paired_source(z, file, original_dates):
    counter=collections.Counter()
    matched=set()
    region=file["region"]
    key_columns={"BIO_SITE_ID","BIO_SAMPLE_DATE","WQ_COUNT",
                 "WQ_MEAN_LOQ1","WQ_MEAN_LOQ2","WQ_MEAN_LOQ3",
                 "WQ_COUNT_L_LOQ","WQ_COUNT_G_LOQ"}
    for row in source.csv_rows(z,file["path"],key_columns):
        counter["water_quality_original_summary_rows"]+=1
        site=str(row["BIO_SITE_ID"]).strip()
        d=f.parsed_date(row["BIO_SAMPLE_DATE"])
        if d is None:
            counter["invalid_original_bio_sample_date"]+=1
            continue
        date=d.date().isoformat()
        count=source.finite_numeric(row["WQ_COUNT"])
        if count is None or count<=0:
            counter["no_valid_source_chemical_measurement_count"]+=1
            continue
        mean2=source.finite_numeric(row["WQ_MEAN_LOQ2"])
        if mean2 is None:
            counter["missing_source_nutrient_LOQ2_mean"]+=1
            continue
        if mean2<0:
            counter["negative_numeric_LOQ2_mean"]+=1
            continue
        counter["valid_original_quality_summary_rows"]+=1
        if (region,site,date) in original_dates:
            matched.add((region,site,date))
            counter["matched_source_WQ_rows"]+=1
            # Track missingness and LOQ rules, not levels/outcome correlations.
            mean1=source.finite_numeric(row["WQ_MEAN_LOQ1"])
            mean3=source.finite_numeric(row["WQ_MEAN_LOQ3"])
            if mean1 is not None and mean3 is not None:
                counter["matched_three_LOQ_mean_versions_present"]+=1
            l=source.finite_numeric(row["WQ_COUNT_L_LOQ"])
            h=source.finite_numeric(row["WQ_COUNT_G_LOQ"])
            if l is not None and h is not None and l+h>0:
                if l>0:counter["matched_some_measurements_below_LOQ"]+=1
    return matched,counter

def aggregate_analytes(file_stats):
    groups=collections.defaultdict(list)
    for q in file_stats:groups[q["analyte"]].append(q)
    out={}
    for name,items in sorted(groups.items()):
        matched=set()
        c=collections.Counter()
        for item in items:
            matched.update(item.pop("_matched"))
            c.update(item["coverage_counters"])
        sites=set((r,s) for r,s,d in matched)
        by_region=collections.defaultdict(set)
        for reg,site,date in matched:by_region[reg].add(site)
        nreg=sum(len(sites)>=PREOUTCOME_GATES["per_valid_region_min_sites"]
                 for sites in by_region.values())
        source_pass=(len(matched)>=PREOUTCOME_GATES["min_matched_distinct_bio_site_date"] and
                     len(sites)>=PREOUTCOME_GATES["min_matched_distinct_bio_site"] and
                     nreg>=PREOUTCOME_GATES["min_regions_with_valid_source_sites"])
        out[name]={
            "publisher_analyte_filename_token":name,
            "original_regions_with_analyte":sorted(set(q["region"] for q in items)),
            "original_source_file_count":len(items),
            "matched_quality_site_date_pairs":len(matched),
            "matched_quality_sites":len(sites),
            "matched_quality_regions":sorted(by_region),
            "regions_with_10_or_more_matched_sites":nreg,
            "region_unique_site_count":{r:len(v) for r,v in sorted(by_region.items())},
            "measured_source_LOQ2_count_stats":dict(c),
            "source_coverage_preoutcome_gate_PASS":bool(source_pass),
            "analyte_unit_LOQ_definition_not_verified":True,
            "WQ_window_stats_not_single_same_day_field_chemistry":True,
        }
    return out

def main(args):
    out={"schema":"GESUI_CHEMPOP_V2_WQ_NUTRIENT_SOURCE_ONLY_V1",
         "source_only":True,"not_causal":True,
         "not_biological_outcome_opened":True,
         "original_DOI":"10.5285/23f194b0-1e21-4f14-be9c-6eabdcc1feb7",
         "source_thresholds":PREOUTCOME_GATES,
         "status":"SOURCE_STOP",
         "candidate_files":[],
         "analytes":{}}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_ChemPop_analyte_source_") as temp:
            d=Path(temp)
            zipinfo=source.retrieve(source.EIDC_URL,d/"official.zip",source.SOURCE_DOWNLOAD_LIMIT)
            values, scoreinfo=source.score_reference(d)
            out["source"]={"original_EIDC_ZIP_SHA256":zipinfo["SHA256"],
                           "original_BMWP_RData_SHA256":scoreinfo["SHA256"]}
            if len(values)!=EXPECTED_SCORED:
                raise RuntimeError("Published British BMWP scoring dictionary changed")
            with zipfile.ZipFile(d/"official.zip") as z:
                files=z.namelist()
                if len(files)!=326:
                    raise RuntimeError("Publisher v2 original 326-file source identity drift")
                raw,weight,orders,source_state=f.construct_eligible(z,files,values)
                if (len(raw)!=EXPECTED_SOURCE_SAMPLES
                    or source_state["source_gate"]["eligible_site_n"]!=EXPECTED_SOURCE_SITES):
                    raise RuntimeError("Original 2010-18 source site/sample gate changed")
                bio_source_date_pairs={
                    (row["region"],row["site"],row["sample_date"]) for row in raw
                }
                if len(bio_source_date_pairs)<EXPECTED_SITEYEARS:
                    raise RuntimeError("Source sample dates unexpectedly fewer than siteyears")
                out["original_source_biology"]={
                    "n_eligible_original_samples":len(raw),
                    "n_distinct_site_date":len(bio_source_date_pairs),
                    "n_distinct_sites":EXPECTED_SOURCE_SITES,
                    "n_original_distinct_site_year_expected":EXPECTED_SITEYEARS}
                candidates=candidate_names(files)
                if not candidates:raise RuntimeError("No original N/P candidate analyte files")
                itemout=[]
                for q in candidates:
                    matches,counters=paired_source(z,q,bio_source_date_pairs)
                    item={**q,"coverage_counters":dict(counters),"_matched":matches}
                    itemout.append(item)
                out["analytes"]=aggregate_analytes(itemout)
                out["candidate_files"]= [{k:v for k,v in x.items() if not k.startswith("_")}
                    for x in itemout]
                # No occurrence/trait responses, concentrations or correlations opened.
                out["status"]="SOURCE_ONLY_CHEMPOP_NP_CANDIDATES_AUDITED"
                out["n_candidate_analytes"]=len(out["analytes"])
                out["n_source_supported_analytes"]=sum(
                    x["source_coverage_preoutcome_gate_PASS"] for x in out["analytes"].values())
    except Exception as ex:
        out.update({"status":"SOURCE_STOP","error":type(ex).__name__+": "+str(ex)[:600]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_CHEMPOP_NP_SOURCE_ONLY="+json.dumps({
       "status":out["status"],"error":out.get("error"),
       "source":out.get("source"),
       "original_source_biology":out.get("original_source_biology"),
       "n_candidate_analytes":out.get("n_candidate_analytes"),
       "n_source_supported_analytes":out.get("n_source_supported_analytes"),
       "analytes":{name:{k:v for k,v in info.items()
                       if k in ("matched_quality_site_date_pairs",
                                "matched_quality_sites",
                                "original_regions_with_analyte",
                                "region_unique_site_count",
                                "source_coverage_preoutcome_gate_PASS")}
                   for name,info in out.get("analytes",{}).items()}
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
