#!/usr/bin/env python3
"""Fresh original 2011-24 GBIF all-Plantae data on 247 source-fixed global WPI ports.

No WPI butterfly candidate-host inclusion in all-Plantae query. TaxonKey=6,
year 2011..2024, radius 5km, all basisOfRecord with coordinates, proper facet
pagination. Report raw indexed counts, 50-record rarefied speciesKey richness,
site-level original GBIF facet source receipts/sha; no individual observers.
"""
from __future__ import annotations

import argparse,hashlib,json,math,time,urllib.error,urllib.parse,urllib.request
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

API="https://api.gbif.org/v1/occurrence/search"
SOURCE=Path("results/global_all_Plantae_WWp2015_20261110/ORIGINAL_42COUNTRY_WWp_STRATIFIED_PORTS_PRE_GBIF_ALLFLORA_SOURCE.json")
YEAR="2011,2024"
MIN_RAREFACTION_N=50
RADIUS_KM=5
HEAD={"User-Agent":"gesui-global-AllPlantae-2011-24-42country-WWp-source-pilot/0.1 (research; zuizui0223/gesui)",
      "Accept":"application/json"}
PARAMS_FIXED={"taxonKey":"6","year":YEAR,"hasCoordinate":"true",
              "hasGeospatialIssue":"false","occurrenceStatus":"PRESENT",
              "facet":"speciesKey","facetLimit":"1000","limit":"0"}
MIN_PACE_S=1.15
MAX_FACET_PAGES=50
MAX_ATTEMPTS=4

def sha(data):return hashlib.sha256(data).hexdigest()

def exact_sample_expected_species(counts,n=50):
    """E[S_n] of a WITHOUT-replacement, equal-GBIF-record-count speciesKey subsample."""
    N=sum(counts.values())
    if N<n or n<1:return None
    log_denom=math.lgamma(N+1)-math.lgamma(n+1)-math.lgamma(N-n+1)
    total=0.
    for key,k in counts.items():
        if k<=0 or k>N:raise ValueError("Invalid faceted original speciesKey record count")
        if N-k<n:
            total+=1.
            continue
        log_num=math.lgamma(N-k+1)-math.lgamma(n+1)-math.lgamma(N-k-n+1)
        delta=min(0.,log_num-log_denom)
        total+=-math.expm1(delta)
    if total<-.01 or total>min(n,len(counts))+.01:
        raise ValueError("Original GBIF accepted speciesKey rarefaction outside feasible bounds")
    return min(float(n),max(0.,total))

def request_json(url,last_request):
    for attempt in range(MAX_ATTEMPTS):
        delay=max(0.,MIN_PACE_S-(time.monotonic()-last_request[0]))
        if delay:time.sleep(delay)
        last_request[0]=time.monotonic()
        try:
            request=urllib.request.Request(url,headers=HEAD)
            with urllib.request.urlopen(request,timeout=75) as response:
                raw=response.read()
                if response.status!=200:raise ValueError("Unexpected non-200 GBIF response")
            obj=json.loads(raw)
            if not isinstance(obj,dict) or obj.get("error"):
                raise ValueError("Invalid GBIF original occurrence facet response")
            return obj,sha(raw)
        except (urllib.error.HTTPError,urllib.error.URLError,ValueError,TimeoutError) as e:
            if isinstance(e,urllib.error.HTTPError) and e.code not in (408,425,429,500,502,503,504):
                raise
            if attempt==MAX_ATTEMPTS-1:raise
            time.sleep(3*(attempt+1))
    raise RuntimeError("Unreachable original GBIF request")

def one_port(row,last_request):
    lat=float(row["original_port_lat"]);lon=float(row["original_port_lon"])
    params=dict(PARAMS_FIXED)
    params["geoDistance"]=f"{lat:.7f},{lon:.7f},{RADIUS_KM}km"
    pages=[]
    species={}
    original_record_count=None
    for page in range(MAX_FACET_PAGES):
        params["facetOffset"]=str(page*1000)
        url=API+"?"+urllib.parse.urlencode(params)
        result,blobsha=request_json(url,last_request)
        n=result.get("count")
        if type(n) is not int or n<0:raise ValueError("Original GBIF occurrence query count must be a nonnegative integer")
        if original_record_count is None:original_record_count=n
        if n!=original_record_count:raise ValueError("GBIF original occurrence index count changed across species facet pagination")
        original_facets=[x for x in result.get("facets",[]) if str(x.get("field","")).upper() in ("SPECIES_KEY","SPECIESKEY")]
        if len(original_facets)!=1:raise ValueError("GBIF all-Plantae faceted speciesKey source missing/ambiguous")
        facets=original_facets[0].get("counts",[])
        if not isinstance(facets,list) or len(facets)>1000:raise ValueError("Invalid original GBIF speciesKey facet page")
        for r in facets:
            key=int(r["name"]);v=int(r["count"])
            if v<0 or key<=0:raise ValueError("Non-positive original GBIF speciesKey/count")
            if key in species:raise ValueError("SpeciesKey repeated across GBIF facet pages")
            species[key]=v
        pages.append({"query_url":url,"response_sha256":blobsha,
                      "source_count":n,"num_faceted_species_keys":len(facets),
                      "facetOffset":page*1000})
        if len(facets)<1000:
            break
    else:
        raise ValueError("Original GBIF all-Plantae facet needed more than 50 pages; refuse truncation")
    species_records=sum(species.values())
    if species_records>original_record_count:
        raise ValueError("Original GBIF speciesKey faceted records exceeds all Plantae records")
    rare=exact_sample_expected_species(species,MIN_RAREFACTION_N)
    coverage=species_records/original_record_count if original_record_count>0 else None
    qpass=(original_record_count>=50 and species_records>=50 and
           len(species)>=10 and coverage is not None and coverage>=.7)
    return {
        "original_port_id":int(row["original_port_id"]),
        "source_country":row["original_source_country"],
        "indexed_GBIF_all_Plantae_2011_2024_records":original_record_count,
        "faceted_original_speciesKey_records":species_records,
        "faceted_original_distinct_Plantae_speciesKeys":len(species),
        "facet_record_coverage":coverage,
        "fully_paginated_original_GBIF_species_key_pages_n":len(pages),
        "original_exact_hypergeometric_50_plant_record_rarefied_speciesKey":rare,
        "source_quality_min50_faceted_records_min10_species_coverage70pct":bool(qpass),
        "original_GBIF_query_and_response_SHA_pages":pages,
        "all_Plantae_taxonKey_6_NOT_candidate_butterfly_hosts_only":True,
        "not_native_alien_or_plant_abundance_or_effluent_exposure":True,
    },species

def census(root,cache_path,out):
    source_path=root/SOURCE
    source=json.loads(source_path.read_text(encoding="utf-8"))
    if source["schema"]!="gesui_WPI2019_42country_WWp2015_selected6percountry_AllPlantae_2011_24_PREOUTCOME_FRAME_v01":
        raise ValueError("Source-only global plant-free 42-country cohort is unavailable")
    if source["GBIF_new_all_Plantae_requests_in_this_stage"]!=0:
        raise ValueError("Source selection was not outcome-blind")
    expected=source["all_pilot_port_georeferenced_original_sites"]
    if len(expected)!=247 or source["countries_in_source_pilot"]!=42:
        raise ValueError("Source original 247 ports / 42 countries changed")
    if len({s["original_port_id"] for s in expected})!=247:
        raise ValueError("Original port IDs duplicated")
    prior={}
    if cache_path.exists():
        for line in cache_path.read_text(encoding="utf-8").splitlines():
            if line.strip():
                x=json.loads(line)
                p=int(x["original_port_id"])
                if p in prior:raise ValueError("Cached GBIF port repeated: "+str(p))
                prior[p]=x
    if not set(prior).issubset({x["original_port_id"] for x in expected}):
        raise ValueError("GBIF cache includes non-source port")
    api_utc=datetime.now(timezone.utc).isoformat()
    cache_path.parent.mkdir(parents=True,exist_ok=True)
    last_request=[0.]
    for idx,s in enumerate(expected,1):
        pid=s["original_port_id"]
        if pid in prior:continue
        try:
            val,species=one_port(s,last_request)
            saved={**val,"original_speciesKey_count_pairs":sorted([[k,v] for k,v in species.items()]),
                   "GBIF_status":"COMPLETE","GBIF_error":None}
        except Exception as e:
            saved={"original_port_id":pid,"source_country":s["original_source_country"],
                   "GBIF_status":"QUERY_FAILED_NOT_ZERO",
                   "GBIF_error":str(e)[:250],
                   "source_quality_min50_faceted_records_min10_species_coverage70pct":False,
                   "original_speciesKey_count_pairs":[],
                   "original_GBIF_query_and_response_SHA_pages":[]}
        with cache_path.open("a",encoding="utf-8") as f:
            f.write(json.dumps(saved,sort_keys=True,separators=(",",":"))+"\n")
        prior[pid]=saved
        if idx%20==0 or idx==len(expected):
            print("ALL_PLANTAE_GBIF_SOURCE_QUERY",idx,"/",len(expected),
                  "complete",sum(x["GBIF_status"]=="COMPLETE" for x in prior.values()),
                  "support",sum(x["source_quality_min50_faceted_records_min10_species_coverage70pct"] for x in prior.values()),
                  "API_errors",sum(x["GBIF_status"]!="COMPLETE" for x in prior.values()),flush=True)
    if set(prior)!={x["original_port_id"] for x in expected}:
        raise ValueError("Original global sample GBIF response incomplete")
    summary=[]
    for s in expected:
        row=prior[s["original_port_id"]]
        safe={k:v for k,v in row.items() if k not in ("original_speciesKey_count_pairs",)}
        if row["GBIF_status"]=="COMPLETE":
            counts=dict(row["original_speciesKey_count_pairs"])
            if len(counts)!=row["faceted_original_distinct_Plantae_speciesKeys"]:
                raise ValueError("Facet cache speciesKeys not preserved")
            if sum(counts.values())!=row["faceted_original_speciesKey_records"]:
                raise ValueError("Facet record sums differ")
            observed_rarefy=exact_sample_expected_species(counts,50)
            prev_rarefy=row["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"]
            if observed_rarefy is None:
                if prev_rarefy is not None:raise ValueError("Invalid computed rarefaction on low-count ports")
            elif abs(observed_rarefy-prev_rarefy)>1e-9:
                raise ValueError("Original cache rarefaction response is inconsistent")
        summary.append(safe)
    usable=[x for x in summary if x["source_quality_min50_faceted_records_min10_species_coverage70pct"]]
    bycountry=defaultdict(list)
    for x in usable:bycountry[x["source_country"]].append(x)
    variation={g:len({next(z["source_WWp_log1p"] for z in expected if z["original_port_id"]==v["original_port_id"]) for v in vs})>=2
               for g,vs in bycountry.items() if len(vs)>=2}
    possible=(len(usable)>=120 and len(bycountry)>=25 and sum(variation.values())>=20)
    out={
        "schema":"gesui_2011_24_global_AllPlantae_42country_WWp_stratified247_5km_GBIF_50record_rarefaction_v01",
        "api_first_utc":api_utc,
        "source_only_study_manifest_sha256":hashlib.sha256(source_path.read_bytes()).hexdigest(),
        "source_selected_247_fixed_WPI2019_port_IDs_42_countries":True,
        "original_GBIF_current_index_now_not_time_frozen_GBIF_download":True,
        "GBIF_ALL_PLANTAE_taxonKey":6,
        "GBIF_year_2011_2024_fixed":YEAR,
        "public_port_geoDistance_radius_km":RADIUS_KM,
        "GBIF_filter_2011_2024":PARAMS_FIXED,
        "source_GBIF_original_facet_response_cache_SHA256":hashlib.sha256(cache_path.read_bytes()).hexdigest(),
        "original_species_facet_key_counts_only_in_actions_artifact_not_repository":True,
        "original_site_GBIF_facet_summary":summary,
        "original_GBIF_result_completion_summary":{
            "fixed_ports":len(summary),
            "source_qualifying_all_plantae_ports":len(usable),
            "original_query_success_sites":sum(x["GBIF_status"]=="COMPLETE" for x in summary),
            "original_query_failed_sites_not_zero":sum(x["GBIF_status"]!="COMPLETE" for x in summary),
            "original_country_count_with_qualifying_GBIF_allPlantae_sites":len(bycountry),
            "countries_two_eligible_sites_and_WWp_varied":sum(variation.values())
        },
        "original_source_coverage_gate":{
            "min_all_Plantae_ports_w50_plant_speciesKey_records_and_10_species_coverage70":120,
            "min_countries_covered":25,
            "min_countries_with_2plus_plant_sites_and_distinct_WWp":20,
            "ALL_PLANTAE_WWp_statistical_comparison_ALLOWED":possible
        },
        "GBIF_2011_24_all_plants_are_vegetation_RECORDS_not_floristic_community_sample":True,
        "no_identification_of_local_indigenous_or_introduced_status":True,
        "no_sewage_derived_NP_or_biomass_causal_effect_tested":True
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    value=census(root,a.cache,a.out)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(value,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"quality":value["original_GBIF_result_completion_summary"],
                      "gate":value["original_source_coverage_gate"]},indent=2))
if __name__=="__main__":
    main()
