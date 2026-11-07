#!/usr/bin/env python3
"""NEW GBIF Plantae study at 36 river sites source-blind matched before outcome.

Source pairs are immutable (official EEA reported 2022 measured WWTP P mass,
flow, catchment area, crop and human density). GBIF Plantae speciesKey recorded
2022–2024 only, 1km radius. No response-driven site replacement and no
claims of causal wastewater effects or plant abundance.
"""
from __future__ import annotations
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
import argparse,csv,hashlib,json,math,time
import numpy as np
from scipy.stats import binomtest

GBIF="https://api.gbif.org/v1/occurrence/search"
GBIF_YEAR="2022,2024"
GBIF_TAXON_KEY=6
GBIF_RADIUS_KM=1
DEPTH_PRIMARY=30
DEPTH_SENS=10
SOURCE_PAIRS=18
SOURCE_STUDY_RUN=42377235104
USER_AGENT="gesui-EEA2022-measured-P-new-plantae-matched-rivers-v01"

def gbif(params):
    last=None
    for k in range(4):
        try:
            u=GBIF+"?"+urlencode(params)
            q=Request(u,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
            with urlopen(q,timeout=65) as r: body=r.read()
            ans=json.loads(body)
            if "count" not in ans or "error" in ans:raise ValueError("Missing GBIF occurrence count or error")
            return ans
        except (HTTPError,URLError,TimeoutError,ValueError) as e:
            last=e
            if k<3:time.sleep((k+1)*4)
    raise RuntimeError("GBIF source query unsuccessful: "+str(last))

def site_facets(site):
    lat=float(site["latitude"]);lon=float(site["longitude"])
    params={
     "taxonKey":GBIF_TAXON_KEY,"year":GBIF_YEAR,
     "geoDistance":f"{lat:.7f},{lon:.7f},{GBIF_RADIUS_KM}km",
     "hasCoordinate":"true","hasGeospatialIssue":"false",
     "occurrenceStatus":"PRESENT","facet":"speciesKey","facetLimit":1000,
     "limit":0}
    keys={};total=None;pages=0
    while True:
        params["facetOffset"]=pages*1000
        data=gbif(params)
        n=int(data["count"])
        if total is None:total=n
        if n!=total:raise ValueError("GBIF count changed across facets")
        facet=next((f.get("counts",[]) for f in data.get("facets",[])
                    if str(f.get("field","")).upper() in ("SPECIES_KEY","SPECIESKEY")),[])
        if total>0 and not facet and pages==0:
            raise ValueError("GBIF no speciesKey facet for positive record site")
        for record in facet:
            key=int(record["name"])
            if key in keys:raise ValueError("Duplicate facet speciesKey page")
            keys[key]=int(record["count"])
        pages+=1
        if len(facet)<1000:break
        if pages>12:raise ValueError("Unexpectedly high speciesKey pages")
        time.sleep(.12)
    amount=sum(keys.values())
    if amount>total:raise ValueError("GBIF species counts more than total Plantae records")
    return {"station_id":site["station_id"],"source_status":site["source_status"],
       "pair_id":site["pair_id"],"latitude":lat,"longitude":lon,
       "source_study_run_id":SOURCE_STUDY_RUN,"gbif_years":GBIF_YEAR,
       "site_radius_km":GBIF_RADIUS_KM,"all_Plantae_records":total,
       "specieskey_resolved_records":amount,"distinct_species_keys":len(keys),
       "specieskey_counts":keys,"facet_pages":pages}

def exact_rarefaction(counts,depth):
    a=[int(v) for v in counts if int(v)>0]
    n=sum(a)
    if n<depth:return None
    def logcomb(n,k):
        return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
    denom=logcomb(n,depth)
    return float(sum(1-math.exp(logcomb(n-v,depth)-denom)
                     if n-v>=depth else 1 for v in a))

def evaluate(pairs,obs):
    by={x["station_id"]:x for x in obs}
    pairout=[]
    for pair in pairs:
        a=by[pair["source_exposed"]["station_id"]]
        b=by[pair["source_control"]["station_id"]]
        source_counts=[a["specieskey_resolved_records"],b["specieskey_resolved_records"]]
        ratio=max(source_counts)/min(source_counts) if min(source_counts)>0 else float("inf")
        def contrast(depth,min_species,balanced):
            quality=all(x["specieskey_resolved_records"]>=depth and
                        x["distinct_species_keys"]>=min_species for x in (a,b))
            quality=quality and (not balanced or ratio<=10)
            if not quality:return {"eligible":False,"richness_log_difference_source_minus_control":None}
            S=[exact_rarefaction(x["specieskey_counts"].values(),depth) for x in (a,b)]
            return {"eligible":True,
                    "expected_richness_record_standardized":S,
                    "richness_log_difference_source_minus_control":
                       math.log1p(S[0])-math.log1p(S[1])}
        A=contrast(DEPTH_PRIMARY,6,True)
        B=contrast(DEPTH_SENS,4,False)
        species_A=set(a["specieskey_counts"])
        species_B=set(b["specieskey_counts"])
        jacc=1-len(species_A&species_B)/len(species_A|species_B) if (species_A|species_B) else None
        pairout.append({
         "pair_id":pair["pair_id"],"block":pair["spatial_block"],
         "exposed_id":a["station_id"],"control_id":b["station_id"],
         "record_count_exposed_control":source_counts,
         "distinct_specieskeys_exposed_control":[a["distinct_species_keys"],b["distinct_species_keys"]],
         "record_depth_ratio":ratio if math.isfinite(ratio) else None,
         "source_P_load_tonyr":pair["source_exposed"]["eea2022_reported_measured_upstream_P_ton_yr"],
         "source_matching_distance_km":pair["between_site_distance_km"],
         "source_matching_score":pair["source_matching_score"],
         "jaccard_recorded_species_dissimilarity_descriptive":jacc,
         "primary_depth30":A,"sensitivity_depth10":B,
        })
    def summary(field,required_n=10,required_blocks=4):
        selected=[x for x in pairout if x[field]["eligible"]]
        d=[x[field]["richness_log_difference_source_minus_control"] for x in selected]
        pos=sum(v>0 for v in d);neg=sum(v<0 for v in d);ties=sum(v==0 for v in d)
        n=pos+neg
        valid=(len(selected)>=required_n and
               len({x["block"] for x in selected})>=required_blocks)
        return {
          "support_status":"DESCRIPTIVE_PAIRED_STATISTICS_SUPPORTED" if valid else "INSUFFICIENT_SITE_DEPTH_OR_GEOGRAPHIC_INDEPENDENCE",
          "n_eligible_pairs":len(selected),
          "eligible_spatial_blocks":len({x["block"] for x in selected}),
          "mean_log_richness_source_minus_control":float(np.mean(d)) if d else None,
          "median_log_richness_source_minus_control":float(np.median(d)) if d else None,
          "positive_pairs":pos,"negative_pairs":neg,"equal_pairs":ties,
          "exact_two_sided_sign_test_descriptive_p":
            float(binomtest(pos,n,.5,alternative="two-sided").pvalue)
            if valid and n>0 else None,
          "not_randomized_source_assignment":True
        }
    return {
      "schema":"gesui_EEA2022_measured_P_load_blind_matched_plant_recording_v01",
      "status":"FIRST_NEW_GBIF_ALL_PLANTAE_MATCHED_EEA_REPORTED_P_SOURCE_PILOT",
      "source_only_selection_original_run_id":SOURCE_STUDY_RUN,
      "source_2022_measured_P_WWTP_load_is_not_verified_plant_NP_exposure":True,
      "distinct_source_pairs":len(pairs),
      "distinct_new_GBIF_plant_sites":len(obs),
      "specieskeys_union":len(set(k for x in obs for k in x["specieskey_counts"])),
      "gbif_2022_2024_window":GBIF_YEAR,
      "primary":summary("primary_depth30"),
      "sensitivity":summary("sensitivity_depth10"),
      "paired":pairout,
      "limits":[
        "GBIF Plantae index includes nonvascular plant species and is not plant biomass or coverage",
        "Source-zero sites are not clean wastewater-free rivers or nutrient-free controls",
        "EEA 2022 reported P mass is not direct plant-soil/floodplain P exposure, and source persistence over 2022-2024 is unverified",
        "Nearby stream sampling opportunities and recorder bias may be unmatched after source-only covariates",
        "Pairs may be nested in watersheds and 3-degree blocks are not independent experiments",
        "No causal inference or nutrient-based selection effects claimed"
      ]}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--design",type=Path,required=True)
    ap.add_argument("--receipt",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    src=json.loads(args.design.read_text())
    receipt=json.loads(args.receipt.read_text())
    if hashlib.sha256(args.design.read_bytes()).hexdigest()!=receipt["original_files_SHA256"]["SOURCE_MATCHED_PAIR_SELECTION.json"]:
        raise ValueError("Original matched source-only cohort digest changed")
    qc=src["source_design"];pairs=src["pairs"]
    if len(pairs)!=SOURCE_PAIRS or qc["plant_sites_fixed"]!=36 or not qc["source_only_selection"]:
        raise ValueError("Matched source-only cohort not fixed")
    if qc["observed_Peters_2022_TP_concentrations_used_for_selection"] or qc["previous_GBIF_plant_community_used_for_selection"]:
        raise ValueError("Invalid source response-blind design")
    sites=[]
    for pair in pairs:
        for role,key in (("positive","source_exposed"),("zero_reported","source_control")):
            site={**pair[key],"pair_id":pair["pair_id"],"source_status":role}
            sites.append(site)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"RESPONSE_GATES_FIXED_BEFORE_GBIF.json").write_text(
         json.dumps({"original_source_design_sha256":hashlib.sha256(args.design.read_bytes()).hexdigest(),
            "planned_sites":len(sites),"GBIF_years":GBIF_YEAR,"GBIF_buffer_km":GBIF_RADIUS_KM,
            "primary_min_records":30,"primary_min_species":6,
            "primary_within_pair_record_ratio_max":10,
            "sensitivity_min_records":10,"sensitivity_min_species":4,
            "primary_paired_sign_test_requires_pairs":10,
            "primary_paired_sign_test_requires_spatial_blocks":4,
            "sites_fixed_before_new_plant_records":True},indent=2,sort_keys=True)+"\n")
    print("GESUI_NEW_GBIF_PLANTAE_P_MATCHED_FROZEN",
          json.dumps({"source_pairs":len(pairs),"site_count":len(sites),
                      "gbif_years":GBIF_YEAR,"source_design_digest":hashlib.sha256(args.design.read_bytes()).hexdigest()}),flush=True)
    observations=[]
    for i,site in enumerate(sites,1):
        obs=site_facets(site)
        observations.append(obs)
        print("GESUI_NEW_EEA_P_MATCHED_PLANTAE_SITE",json.dumps({
          "site_index":i,"station_id":obs["station_id"],"source_role":obs["source_status"],
          "recorded_plantae":obs["all_Plantae_records"],
          "resolved_species_keys":obs["distinct_species_keys"]}),flush=True)
        time.sleep(.35)
    output=evaluate(pairs,observations)
    (args.out/"ALL_NEW_PLANT_SPECIES_FACETS.json").write_text(
        json.dumps({"schema":"new_source_only_EEA2022_plant_2022_2024",
                    "source_selection_run":SOURCE_STUDY_RUN,
                    "original_design_sha256":hashlib.sha256(args.design.read_bytes()).hexdigest(),
                    "plants":observations},ensure_ascii=False,separators=(",",":"))+"\n")
    (args.out/"PAIRED_PLANTAE_RESULT.json").write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA2022_NEW_MATCHED_PLANTAE_RESULT",
          json.dumps({k:v for k,v in output.items() if k not in ("paired","limits")},sort_keys=True),flush=True)

if __name__=="__main__":
    main()
