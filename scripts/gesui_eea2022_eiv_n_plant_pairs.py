#!/usr/bin/env python3
"""Harmonized European vascular plant N-affinity analysis of EEA2022 river pairs.

The 18 EEA P-load river source/control pairs and all 36 GBIF 2022-2024 plant
sites were fixed BEFORE trait values were retrieved. GBIF species taxonomy is
queried by original speciesKey, exact accepted vascular species binomial
matches only to Tichy et al. (2023) European Ellenberg-type NUTRI.
No fuzzy, genus fallback, missing-trait imputation or causal sewage claims.
"""
from __future__ import annotations
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import quote
import argparse,hashlib,io,json,math,re,time
import numpy as np,pandas as pd,requests

EIV_URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
GBIF_PREFIX="https://api.gbif.org/v1/species/"
EXPECTED_KEYS=1006
EXPECTED_SITES=36
EXPECTED_PAIRS=18
MIN_TRAIT_SPECIES_PER_SITE=10
MIN_TRAIT_COVERAGE_FRACTION=.40
MIN_PAIRS_FOR_DESCRIPTIVE_P=10
MIN_BLOCKS_FOR_DESCRIPTIVE_P=4
TRAIT_SOURCE="Tichy et al 2023, European Ellenberg-type NUTRI nutrient-affinity, DOI 10.1111/jvs.13168"
SOURCE_RUN=42377235104
PLANT_RUN=42377454765
WORKERS=6

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def get(url,maxb=1_000_000):
    last=None
    for attempt in range(4):
        try:
            r=requests.get(url,timeout=(12,60),
                 headers={"User-Agent":"gesui-pollution-ecology-trait-coverage/0.1",
                          "Accept":"application/json,text/csv,*/*"})
            r.raise_for_status()
            if len(r.content)>maxb:raise ValueError("Unexpected source size cap exceeded")
            return r.content
        except (requests.RequestException,ValueError) as ex:
            last=ex
            if attempt<3:time.sleep((attempt+1)*2.0)
    raise RuntimeError("Public source unavailable after retry "+url+": "+str(last))

def simple_name(x):
    if not isinstance(x,str):return ""
    return re.sub(r"\s+"," ",x.strip()).casefold()

def load_eiv():
    raw=get(EIV_URL,2_000_000)
    if hashlib.sha256(raw).hexdigest()!=EIV_SHA:
        raise ValueError("European original nutrient indicator source checksum drift")
    df=pd.read_csv(io.BytesIO(raw),encoding="utf-8-sig",low_memory=False)
    if len(df)!=8908 or "NUTRI" not in df or "Unnamed: 0" not in df:
        raise ValueError("Harmonized European plant indicator schema unexpected")
    df["name_norm"]=df["Unnamed: 0"].map(simple_name)
    df["N_numeric"]=pd.to_numeric(df["NUTRI"],errors="coerce")
    valid=df["N_numeric"].notna()
    # Nutrient-affinity is ordinal EIV value; numeric >0 is the meaningful
    # nutrient optimum. EIV 'x' and NA are not valid nutrient scores.
    valid &= df.N_numeric.between(1,12)
    valid &= df.name_norm.map(lambda s:len(s.split())==2 and
        not any(x in s for x in ("aggr"," s.l."," s.s."," cf."," ×"," x ")))
    sel=df.loc[valid,["name_norm","N_numeric"]]
    groups=sel.groupby("name_norm")["N_numeric"].agg(["count","min","max"])
    ambiguous=set(groups.loc[groups["max"]-groups["min"]>1e-10].index)
    reliable=sel.loc[~sel.name_norm.isin(ambiguous)].drop_duplicates("name_norm")
    name2N=dict(zip(reliable.name_norm,reliable.N_numeric.astype(float)))
    report={
        "original_source_url":EIV_URL,"original_sha256":EIV_SHA,
        "original_rows":len(df),"numeric_NUTRI_values_in_source":int(df.N_numeric.notna().sum()),
        "species_rank_source_exact_binomials_with_numeric_N":len(sel),
        "ambiguous_duplicate_names_removed":len(ambiguous),
        "unambiguous_exact_names":len(name2N),
        "EIV_numeric_range":[float(reliable.N_numeric.min()),float(reliable.N_numeric.max())],
        "trait_meaning":"European vascular plant nutrient-affinity indicator; not observed river or effluent N/P"
    }
    return name2N,report

def taxonomy_one(key):
    try:
        b=get(GBIF_PREFIX+str(key),40_000)
        obj=json.loads(b)
        canonical=obj.get("canonicalName")
        accepted=(obj.get("rank")=="SPECIES" and obj.get("taxonomicStatus")=="ACCEPTED")
        vascular=("tracheophyta" in str(obj.get("phylum","")).casefold())
        return int(key),{
            "canonicalName":canonical,
            "scientificName":obj.get("scientificName"),
            "rank":obj.get("rank"),"taxon_status":obj.get("taxonomicStatus"),
            "phylum":obj.get("phylum"),"kingdom":obj.get("kingdom"),
            "exact_accepted_vascular_species":bool(accepted and vascular),
            "GBIF_source_status":"OK"}
    except Exception as ex:
        return int(key),{"GBIF_source_status":"ERROR","error":type(ex).__name__+": "+str(ex)[:160]}

def exact_sign(pos,neg):
    n=pos+neg
    if n==0:return None
    return float(min(1.,2*sum(math.comb(n,i) for i in range(min(pos,neg)+1))/2**n))

def site_summary(site,traits):
    pairs=[]
    for key,records in site["specieskey_counts"].items():
        entry=traits.get(int(key))
        if entry and entry["trait_NUTRI"] is not None:
            pairs.append((int(key),float(entry["trait_NUTRI"]),int(records)))
    n=len(pairs)
    coverage=n/int(site["distinct_species_keys"]) if site["distinct_species_keys"] else None
    if not pairs:
        mean=weighted=None
    else:
        mean=float(np.mean([a[1] for a in pairs]))
        weighted=float(np.average([a[1] for a in pairs],weights=[a[2] for a in pairs]))
    return {
      "station_id":site["station_id"],"pair_id":site["pair_id"],
      "role":site["source_status"],"all_GBIF_distinct_speciesKeys":site["distinct_species_keys"],
      "all_GBIF_speciesKey_resolved_occurrences":site["specieskey_resolved_records"],
      "matched_numeric_EIV_N_speciesKeys":n,
      "trait_coverage_fraction_all_Plantae_speciesKeys":coverage,
      "unweighted_mean_trait_N":mean,
      "GBIF_occurrence_weighted_mean_trait_N_sensitivity":weighted,
      "primary_EIV_site_quality":n>=MIN_TRAIT_SPECIES_PER_SITE
         and coverage is not None and coverage>=MIN_TRAIT_COVERAGE_FRACTION,
      "trait_speciesKey_list":sorted(a[0] for a in pairs)
    }

def summarize_pairs(pairs,site_rows):
    by={x["station_id"]:x for x in site_rows}
    output=[]
    for pair in pairs:
        a=by[pair["source_exposed"]["station_id"]]
        b=by[pair["source_control"]["station_id"]]
        qual=a["primary_EIV_site_quality"] and b["primary_EIV_site_quality"]
        primary=(a["unweighted_mean_trait_N"]-b["unweighted_mean_trait_N"]) if qual else None
        secondary=(a["GBIF_occurrence_weighted_mean_trait_N_sensitivity"]-
                  b["GBIF_occurrence_weighted_mean_trait_N_sensitivity"]) if qual else None
        output.append({"pair_id":pair["pair_id"],"spatial_block":pair["spatial_block"],
          "source_positive_station":a["station_id"],"source_zero_station":b["station_id"],
          "source_exposed_EEA_2022_measured_P_tonnes_per_year":
              pair["source_exposed"]["eea2022_reported_measured_upstream_P_ton_yr"],
          "matched_trait_species_positive_zero":[a["matched_numeric_EIV_N_speciesKeys"],
                                                  b["matched_numeric_EIV_N_speciesKeys"]],
          "matched_trait_coverage_positive_zero":[a["trait_coverage_fraction_all_Plantae_speciesKeys"],
                                                  b["trait_coverage_fraction_all_Plantae_speciesKeys"]],
          "primary_trait_source_positive_minus_zero":primary,
          "record_weighted_sensitivity":secondary,
          "eligible_primary":bool(qual)})
    keep=[x for x in output if x["eligible_primary"]]
    eligible=len(keep)
    blocks=len(set(x["spatial_block"] for x in keep))
    d=[x["primary_trait_source_positive_minus_zero"] for x in keep]
    pos=sum(x>0 for x in d);neg=sum(x<0 for x in d)
    sufficient=(eligible>=MIN_PAIRS_FOR_DESCRIPTIVE_P and blocks>=MIN_BLOCKS_FOR_DESCRIPTIVE_P)
    return output,{
        "eligible_pairs":eligible,"distinct_spatial_blocks":blocks,
        "trait_source_positive_minus_zero_mean":float(np.mean(d)) if d else None,
        "trait_source_positive_minus_zero_median":float(np.median(d)) if d else None,
        "positive_pairs":pos,"negative_pairs":neg,"tie_pairs":eligible-pos-neg,
        "exact_two_sided_sign_p_descriptive":exact_sign(pos,neg) if sufficient else None,
        "passes_source_trait_coverage_gate":sufficient,
        "record_weighted_difference_mean_sensitivity":
            float(np.mean([x["record_weighted_sensitivity"] for x in keep])) if keep else None
    }

def run(args):
    raw=json.loads(args.plants.read_text())
    source=json.loads(args.pairs.read_text())
    rec=json.loads(args.receipt.read_text())
    if sha(args.pairs)!=rec["original_files_SHA256"]["SOURCE_MATCHED_PAIR_SELECTION.json"]:
        raise ValueError("Fixed source plant pair design digest changed")
    if len(source["pairs"])!=EXPECTED_PAIRS or len(raw["plants"])!=EXPECTED_SITES:
        raise ValueError("Original source-only cohort or GBIF site list changed")
    keys=sorted({int(k) for site in raw["plants"] for k in site["specieskey_counts"]})
    if len(keys)!=EXPECTED_KEYS or raw["source_selection_run"]!=SOURCE_RUN:
        raise ValueError("Original 1006 GBIF original species key inventory changed")
    name2N,source_qc=load_eiv()
    taxonomy={}
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        future={pool.submit(taxonomy_one,key):key for key in keys}
        for f in as_completed(future):
            key,answer=f.result()
            taxonomy[key]=answer
    if set(taxonomy)!=set(keys):
        raise ValueError("Some GBIF species taxonomy lookups never returned")
    errors=sum(x["GBIF_source_status"]!="OK" for x in taxonomy.values())
    if errors>max(10,int(.03*len(keys))):
        raise ValueError("Too many incomplete species taxonomy responses")
    tmap={}
    for key in keys:
        v=taxonomy[key]
        candidate=simple_name(v.get("canonicalName"))
        score=name2N.get(candidate) if v.get("exact_accepted_vascular_species") else None
        tmap[key]={"gbif_speciesKey":key,"gbif_canonical_name":v.get("canonicalName"),
          "gbif_phylum":v.get("phylum"),"gbif_rank":v.get("rank"),
          "gbif_taxon_status":v.get("taxon_status"),
          "exact_accepted_vascular_species":v.get("exact_accepted_vascular_species",False),
          "EIV_name_exact_match":score is not None,
          "trait_NUTRI":float(score) if score is not None else None,
          "taxonomy_source_status":v["GBIF_source_status"]}
    site_rows=[site_summary(s,tmap) for s in raw["plants"]]
    paired,contrast=summarize_pairs(source["pairs"],site_rows)
    source_qc["GBIF_taxonomy_keys_checked"]=len(keys)
    source_qc["GBIF_taxon_API_error_count"]=errors
    source_qc["GBIF_accepted_vascular_speciesKeys"]=sum(
        v["exact_accepted_vascular_species"] for v in tmap.values())
    source_qc["numeric_EIV_N_exact_match_speciesKeys"]=sum(
        v["EIV_name_exact_match"] for v in tmap.values())
    source_qc["all_Plant_site_count"]=len(site_rows)
    source_qc["site_trait_indicator_coverage_good_n"]=sum(
        s["primary_EIV_site_quality"] for s in site_rows)
    source_qc["plant_site_coverage_fraction_median"]=float(np.median([
        v["trait_coverage_fraction_all_Plantae_speciesKeys"] for v in site_rows
        if v["trait_coverage_fraction_all_Plantae_speciesKeys"] is not None]))
    report={
       "schema":"gesui_EEA2022_P_source_Europe_EIV_N_all_Plantae_species_v01",
       "date_utc":datetime.now(timezone.utc).isoformat(),
       "source_EIV":TRAIT_SOURCE,"source_EIV_qc":source_qc,
       "fixed_source_pairs_before_GBIF":SOURCE_RUN,
       "original_new_GBIF_run":PLANT_RUN,
       "new_site_n":len(site_rows),"matched_pair_n":len(paired),
       "trait_coverage_gates":{"min_rated_species_per_site":MIN_TRAIT_SPECIES_PER_SITE,
                              "min_all_Plantae_key_coverage_fraction":MIN_TRAIT_COVERAGE_FRACTION,
                              "min_paired_sites_for_exact_sign":MIN_PAIRS_FOR_DESCRIPTIVE_P,
                              "min_independent_spatial_blocks":MIN_BLOCKS_FOR_DESCRIPTIVE_P},
       "primary_outcome":contrast,
       "paired":paired,
       "site_trait_coverage_and_means":site_rows,
       "external_source_values_are_nutrient_affinity_not_site_P":True,
       "source_P_mass_not_measured_plant_P_exposure":True,
       "no_species_name_fuzzy_or_genus_matching":True,
       "limits":[
         "GBIF 2022-24 plant occurrences are not biomass or quadrat census",
         "NUTRI is a European vascular plant ecological nutrient optimum, not actual phosphorus concentration",
         "Exact name join excludes otherwise valid synonym/species complexes, and selection by trait availability can bias results",
         "2022 EEA source mass is only a marker of registered wastewater P source, with source-zero not untreated clean river",
         "Land-use and hydrological conditions are matched incompletely, and underlying river registry is a third-party vector conversion",
         "Neither positive nor negative matched trait difference is a proven effluent-mediated response"
       ]
    }
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"GBIF_SPECIESKEY_EIV_EXACT_NAME_CROSSWALK.json").write_text(
        json.dumps({"source":TRAIT_SOURCE,"source_SHA256":EIV_SHA,
                    "species_taxonomy_and_matched_trait":list(tmap.values())},
              indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    (args.out/"EIV_N_MATCHED_RIVER_PLANT_TRAIT_RESULTS.json").write_text(
        json.dumps(report,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_EUROPE_NUTRIENT_AFFINITY_SOURCE_COVERAGE",
          json.dumps(source_qc,sort_keys=True,ensure_ascii=False),flush=True)
    print("GESUI_EEA2022_PLANT_TRAIT_N_MATCHED_PAIR_TEST",
          json.dumps(contrast,sort_keys=True,ensure_ascii=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--plants",type=Path,required=True)
    ap.add_argument("--pairs",type=Path,required=True)
    ap.add_argument("--receipt",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    run(ap.parse_args())
