#!/usr/bin/env python3
"""ChemPop-v2 original same-sample taxon/effluent source-quality census.

This is a SOURCE-ONLY, pre-outcome gating script. It does not fit biological
responses, inspect coefficient signs, map wastewater mass to receivers, or
rehost original records. It streams official EIDC original ZIP and published
UK family-score dictionary, verifies stable source identities, counts join,
missingness, year, census coverage, and bounded release documentation.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
import io
import json
import math
import re
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import pyreadr

EIDC_UUID="23f194b0-1e21-4f14-be9c-6eabdcc1feb7"
EIDC_URL="https://data-package.ceh.ac.uk/data/"+EIDC_UUID+".zip"
SCORE_COMMIT="08cea3f06eec93703bb10b53d1015d8fe8c626c4"
SCORE_URL="https://raw.githubusercontent.com/robbriers/biotic/"+SCORE_COMMIT+"/R/sysdata.rda"
PREFIXES=("an","mi","ne","nw","so","sw","th")
TRAITS=["EDF_MEAN","EDF_STD","EDF_Q90","EDF_Q95","BIO_ALTITUDE",
        "BIO_WIDTH","BIO_SILT_CLAY","LC_URBAN_PER","LC_ARABLE_PER",
        "LC_WOODLAND_PER","LC_SEMINATURAL_PER","LC_WOODLAND_AREA",
        "LC_ARABLE_AREA"]
MIN_SAMPLES=500
MIN_SITES=200
MIN_REGIONS=5
MIN_FAMILIES=3
MIN_SCORED_FAMILY_PROP=.70
SOURCE_DOWNLOAD_LIMIT=400*1024*1024

HEADERS={"User-Agent":"gesui-EIDC-ChemPop-v2-authentic-source-census-20261119",
         "Accept":"*/*"}

def retrieve(url, target, maxbytes):
    h=hashlib.sha256()
    n=0
    req=urllib.request.Request(url,headers=HEADERS)
    with urllib.request.urlopen(req,timeout=180) as r,target.open("wb") as fd:
        if r.status!=200:
            raise RuntimeError("Unexpected source status "+str(r.status))
        while True:
            part=r.read(65536)
            if not part:break
            n+=len(part)
            if n>maxbytes:raise RuntimeError("Source exceeds bounded size")
            fd.write(part);h.update(part)
    if n<=0:raise RuntimeError("Empty source ZIP")
    return {"SHA256":h.hexdigest(),"source_archive_bytes":n}

def csv_rows(z,name,required):
    with z.open(name) as body:
        txt=io.TextIOWrapper(body,encoding="utf-8-sig",errors="replace",newline="")
        reader=csv.DictReader(txt)
        cols=set(reader.fieldnames or [])
        if not set(required).issubset(cols):
            raise RuntimeError("Source original schema drift: "+name+" "+repr(sorted(cols)))
        for row in reader: yield row

def score_reference(folder):
    p=folder/"indextable.rda"
    out=retrieve(SCORE_URL,p,300000)
    d=pyreadr.read_r(str(p))
    if "indextable" not in d:raise RuntimeError("No official-published index object")
    table=d["indextable"]
    if not set(["Taxon","BMWP"]).issubset(table.columns):
        raise RuntimeError("Source biotic dictionary absent scores")
    scores={}
    for row in table[["Taxon","BMWP"]].itertuples(index=False):
        tax=str(row.Taxon).strip()
        try:val=float(row.BMWP)
        except (ValueError,TypeError):continue
        if not tax or not math.isfinite(val):continue
        if tax in scores:raise RuntimeError("Duplicate original scoring taxon")
        scores[tax]=val
    if len(scores)!=92:raise RuntimeError("Published BMWP score source count changed")
    return scores,out

def finite_numeric(x):
    try:
        v=float(str(x).strip())
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None

def date_year(s):
    x=str(s or "").strip()
    if not x:return None
    for fmt in ("%Y-%m-%d","%d/%m/%Y","%Y/%m/%d","%Y-%m-%d %H:%M:%S",
                "%d-%b-%Y","%d/%m/%y"):
        try:return datetime.datetime.strptime(x,fmt).year
        except ValueError:pass
    try:return datetime.datetime.fromisoformat(x.replace("Z","+00:00")).year
    except ValueError:return None

def find_one_by_region(files, region, suffix):
    matches=[f for f in files if re.fullmatch(
        "data/cp_"+region+"_minv_"+suffix, f,re.IGNORECASE)]
    if len(matches)!=1:raise RuntimeError("Source region "+region+" "+suffix+" missing/duplicate "+repr(matches))
    return matches[0]

def find_taxa(files,region):
    matches=[f for f in files if re.fullmatch(
        "data/cp_"+region+"_minv_taxonabundance_[0-9]{8}-[0-9]{8}\\.csv",
        f,re.IGNORECASE)]
    if len(matches)!=1:raise RuntimeError("Expected unique original benthic taxon table for "+region)
    return matches[0]

def census(z,files,scores):
    source={"prefixes":list(PREFIXES),
            "n_orig_files":len(files),
            "site_census":{},"taxon_census":{},"regional":{},
            "metadata_context":{},"sample_year_counts":{},
            "year_range":[],
            "sample_scored_ratio":{}}
    record_dates=collections.Counter()
    region_eligible_counts=collections.Counter()
    region_sites_eligible=collections.defaultdict(set)
    dataset_site_count=0
    group_global={}
    taxon_species_groups=set()
    for region in PREFIXES:
        sitefile=find_one_by_region(files,region,"sitevariables\\.csv")
        taxfile=find_taxa(files,region)
        sites={}
        edf=[]
        column_ok=collections.Counter()
        methods=collections.Counter()
        for row in csv_rows(z,sitefile,["BIO_SITE_ID","EDF_MEAN",
                                       "LC_ARABLE_PER","LC_URBAN_PER"]):
            key=str(row["BIO_SITE_ID"]).strip()
            if not key:continue
            if key in sites:raise RuntimeError("Duplicate source siteID within original regional site variables")
            sites[key]={k:finite_numeric(row.get(k)) for k in TRAITS}
            for k,v in sites[key].items():
                if v is not None:column_ok[k]+=1
            if sites[key]["EDF_MEAN"] is not None:edf.append(sites[key]["EDF_MEAN"])
        dataset_site_count+=len(sites)
        family_labels=collections.Counter()
        taxa_rows=0
        n_non_numeric_abundance=0
        n_positive_numeric_abundance=0
        missing_family=0
        missing_order=0
        unknown_site=0
        sample_tuple_mismatch=0
        sample_groups={}
        for row in csv_rows(z,taxfile,["BIO_SITE_ID","BIO_SAMPLE_ID",
                                      "BIO_ANALYSIS_ID","BIO_SAMPLE_DATE",
                                      "BIO_ORD_NAME","BIO_FAM_NAME",
                                      "BIO_ABUNDANCE"]):
            taxa_rows+=1
            site=str(row["BIO_SITE_ID"]).strip()
            aid=str(row["BIO_ANALYSIS_ID"]).strip()
            smp=str(row["BIO_SAMPLE_ID"]).strip()
            ds=str(row["BIO_SAMPLE_DATE"]).strip()
            order=str(row["BIO_ORD_NAME"] or "").strip()
            family=str(row["BIO_FAM_NAME"] or "").strip()
            method=str(row["BIO_ANALYSIS_METHOD"] or "").strip()
            if not order:missing_order+=1
            if not family:missing_family+=1
            if site not in sites:unknown_site+=1
            key=(site,aid)
            native=(smp,ds,method)
            if key not in sample_groups:
                sample_groups[key]={
                    "native":native,"positive_numeric_rows":0,
                    "all_rows":0,"families":set(),"scored":set(),
                    "orders":set(),"source_is_site_matched":site in sites}
            item=sample_groups[key]
            if item["native"]!=native:sample_tuple_mismatch+=1
            item["all_rows"]+=1
            if family:
                item["families"].add(family)
                family_labels[family]+=1
                if family in scores:item["scored"].add(family)
            if order:item["orders"].add(order)
            n=finite_numeric(row["BIO_ABUNDANCE"])
            if n is not None and n>0:
                n_positive_numeric_abundance+=1
                item["positive_numeric_rows"]+=1
            else:n_non_numeric_abundance+=1
            if taxa_rows>9000000:raise RuntimeError("Unexpected regional raw taxa count >9m")
        if sample_tuple_mismatch:raise RuntimeError("Within-analysis site/sample/date/method original identity mismatch")
        eligible=set()
        family_coverage=[]
        family_n=[]
        positive_rate=[]
        n_same_samp_has_site_edf=0
        n_same_samp_year_valid=0
        sample_year_counts=collections.Counter()
        for (site,aid),v in sample_groups.items():
            nfamily=len(v["families"])
            nscore=len(v["scored"])
            coverage=nscore/nfamily if nfamily else 0
            family_coverage.append(coverage)
            family_n.append(nfamily)
            positive_rate.append(v["positive_numeric_rows"]/v["all_rows"])
            yr=date_year(v["native"][1])
            sample_year_counts[str(yr) if yr else "UNKNOWN"]+=1
            record_dates[str(yr) if yr else "UNKNOWN"]+=1
            if yr is not None:n_same_samp_year_valid+=1
            edf_value=(sites.get(site) or {}).get("EDF_MEAN")
            if edf_value is not None:n_same_samp_has_site_edf+=1
            if nfamily>=MIN_FAMILIES and coverage>=MIN_SCORED_FAMILY_PROP and edf_value is not None:
                eligible.add((site,aid))
                region_eligible_counts[region]+=1
                region_sites_eligible[region].add(site)
        source["regional"][region]={
            "n_sites_in_source":len(sites),
            "n_sites_with_finite_EDF_MEAN":len(edf),
            "n_distinct_EDF_MEAN":len(set(round(float(x),6) for x in edf)),
            "EDF_MEAN_min_max":[min(edf),max(edf)] if edf else None,
            "n_raw_taxon_rows":taxa_rows,
            "n_raw_sample_ANALYSIS_ID":len(sample_groups),
            "n_positive_numeric_abundance_rows":n_positive_numeric_abundance,
            "n_nonpositive_or_non_numeric_abundance_rows":n_non_numeric_abundance,
            "n_site_keys_missing_sitevariables":unknown_site,
            "n_empty_family_rows":missing_family,
            "n_empty_order_rows":missing_order,
            "n_unambiguous_sample_date":n_same_samp_year_valid,
            "n_samples_with_EDF_MEAN":n_same_samp_has_site_edf,
            "n_samples_with_3_families_70pct_scored_and_EDF":len(eligible),
            "n_unique_sites_with_eligible_samples":len(region_sites_eligible[region]),
            "n_region_exact_family_names":len(family_labels),
            "n_region_exact_scored_BMWP_families":len(set(family_labels)&set(scores)),
            "family_names_top10":family_labels.most_common(10),
            "source_sitevariable_numeric_coverage":dict(column_ok),
            "year_sample_counts":dict(sorted(sample_year_counts.items())),
            "fraction_samples_3_or_more_families":sum(v>=3 for v in family_n)/len(family_n) if family_n else None,
            "median_family_scored_fraction":sorted(family_coverage)[len(family_coverage)//2] if family_coverage else None,
            "fraction_samples_abundance_allpositive":sum(v==1 for v in positive_rate)/len(positive_rate) if positive_rate else None,
        }
        taxon_species_groups.update(family_labels)
        group_global[region]=[len(sites),len(sample_groups),taxa_rows]
    sample_total=sum(x[1] for x in group_global.values())
    valid=([k for k,v in source["regional"].items() if
            v["n_unique_sites_with_eligible_samples"]>=10 and
            v["n_distinct_EDF_MEAN"]>=5])
    source["aggregated"]={
        "total_sites_in_regional_source":dataset_site_count,
        "total_original_samples":sample_total,
        "total_taxon_rows":sum(x[2] for x in group_global.values()),
        "original_distinct_family_labels":len(taxon_species_groups),
        "n_source_eligible_samples":sum(region_eligible_counts.values()),
        "n_sites_with_source_eligible_samples":sum(len(x) for x in region_sites_eligible.values()),
        "n_regions_with_10_eligible_sites_and_EDF_variation":len(valid),
        "which_regions_pass_primitive_coverage":valid,
        "sample_year_counts":dict(sorted(record_dates.items())),
    }
    return source

def metadata_context(z,files):
    info={"README_EDF_occurrences":[],"README_has_source_metadata":False,
          "rocrate_EDF_occurrences":[]}
    if "readme.html" in files:
        body=z.read("readme.html")[:150000].decode("utf-8","replace")
        info["README_has_source_metadata"]=True
        for mat in list(re.finditer(r"EDF|effluent|sewage|dilution|water quality",body,re.I))[:15]:
            info["README_EDF_occurrences"].append(re.sub(r"<[^>]*>"," ",body[max(0,mat.start()-150):mat.end()+220]).strip()[:350])
    if "ro-crate-metadata.json" in files:
        doc=z.read("ro-crate-metadata.json")[:350000].decode("utf-8","replace")
        for mat in list(re.finditer(r"EDF|effluent|sewage|dilution",doc,re.I))[:10]:
            info["rocrate_EDF_occurrences"].append(doc[max(0,mat.start()-90):mat.end()+160])
    return info

def main(a):
    result={"schema":"GESUI_CHEMPOP_V2_SAMPLE_FAMILY_EDF_SOURCE_CENSUS_V1",
            "source_uuid":EIDC_UUID,"source_only":True,
            "does_not_fit_EEAlink_or_effects":True,
            "no_original_bio_taxon_rows_republished":True,
            "no_confirmed_measured_delivered_WWTP_P_or_N":True,
            "status":"SOURCE_STOP"}
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_chempop_sample_source_") as tmp:
            d=Path(tmp)
            archive=retrieve(EIDC_URL,d/"official.zip",SOURCE_DOWNLOAD_LIMIT)
            scores, score_hash=score_reference(d)
            result["archive_source"]=archive
            result["biotic_reference_source"]=score_hash
            result["biotic_reference_score_n"]=len(scores)
            with zipfile.ZipFile(d/"official.zip") as z:
                files=z.namelist()
                if len(files)!=326 or len(set(files))!=len(files):
                    raise RuntimeError("Original published 326-file v2 archive contents changed")
                for p in PREFIXES:
                    find_one_by_region(files,p,"sitevariables\\.csv")
                    find_taxa(files,p)
                result["metadata_context"]=metadata_context(z,files)
                source=census(z,files,scores)
            result.update(source)
            a0=source["aggregated"]
            result["source_gates"]={
                "min_samples":MIN_SAMPLES,"min_sites":MIN_SITES,"min_regions":MIN_REGIONS,
                "min_scored_families_per_sample":MIN_FAMILIES,
                "min_scored_family_fraction":MIN_SCORED_FAMILY_PROP,
                "enough_source_samples":a0["n_source_eligible_samples"]>=MIN_SAMPLES,
                "enough_source_sites":a0["n_sites_with_source_eligible_samples"]>=MIN_SITES,
                "enough_region_contexts":a0["n_regions_with_10_eligible_sites_and_EDF_variation"]>=MIN_REGIONS,
                "original_edf_treatment_definition_source_locked":False,
                "actual_individual_counts_vs_abundance_class_source_locked":False,
                "sample_analyte_TP_TN_same_original_date_locked":False,
            }
            result["status"]="ORIGINAL_SOURCE_CENSUS_COMPLETED"
    except Exception as e:
        result.update({"status":"SOURCE_STOP","source_error":type(e).__name__+": "+str(e)[:600]})
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_CHEMPOP_V2_SOURCE_CENSUS="+json.dumps({
        "status":result["status"],"reason":result.get("source_error"),
        "aggregated":result.get("aggregated"),
        "regional":{r:{k:v.get(k) for k in [
            "n_sites_in_source","n_sites_with_finite_EDF_MEAN",
            "n_distinct_EDF_MEAN","EDF_MEAN_min_max",
            "n_raw_taxon_rows","n_raw_sample_ANALYSIS_ID",
            "n_source_eligible","n_samples_with_3_families_70pct_scored_and_EDF",
            "n_unique_sites_with_eligible_samples",
            "median_family_scored_fraction","fraction_samples_abundance_allpositive"]}
            for r,v in result.get("regional",{}).items()},
        "source_gates":result.get("source_gates"),
        "metadata_context":result.get("metadata_context")
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
