#!/usr/bin/env python3
"""Published 2026 Grand River 8 municipal WWTP source-locked aquatic-to-riparian
microplastic STAGE responses. Exploratory reanalysis of Wardlaw et al 2026,
NOT an independent published discovery, nutrient or biomass experiment.

Primary observations: upstream/downstream at 8 physical treatment facilities
across 3 original numbered visits within ONE connected river network, from
the same published cross-trophic-level data table. True aquatic insect
larvae, emergent adults, riparian Tetragnathidae spiders are measured.

* Lock source and site/taxon/stage/visit support BEFORE reading MP numbers.
* Require 2 genuine analysed specimens on each side within each site/visit.
* Require ≥6 WWTP sites with ≥2 matched visits across larva/adult/spider.
* Biological/facility unit n=8, NOT 906 individuals or 24 independent WWTPs.
* Compare within-taxon concentration units, do not multiply across stages.
* Do not interpret microplastics transfer as N/P nutrition or biomass gain.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,itertools,json,math,statistics
from pathlib import Path

SCHEMA="GESUI_20261117_GRAND8_GENUINE_WWTP_AQUATIC_ADULT_RIPARIAN_SPIDER_MP_WITHIN_FACILITY_V1"
PUBLISHER="https://doi.org/10.5683/SP4/U0MN9S"
PAPER="https://doi.org/10.1016/j.envpol.2026.129116"
ORIGINAL_DATAFILE_ID=1213261
ORIGINAL_SHA256="1b69eb88fe38f601b11ff466d6a9b3177546a9b489a74c61eaf2a64d770dec66"
N_SOURCE=906
SITES=("Ayr","Caledonia","Guelph","Hespeler","Kitchener","Paris","Preston","Waterloo")
TIMES=("1","2","3")
PRIMARY={
  "Hydropsychidae_larvae":"Larval Hydropsychidae",
  "Hydropsychidae_adults":"Adult Hydropsychidae",
  "Chironomidae_larvae":"Larval Chironomidae",
  "Chironomidae_adults":"Adult Chironomidae",
  "Tetragnathidae_spider":"Tetragnathidae",
}
SECONDARY={"Heptageniidae_larvae":"Larval Heptageniidae",
           "Heptageniidae_adults":"Adult Heptageniidae"}
MIN_REPLICATES_PER_SIDE=2
MIN_INDEPENDENT_WWTP_WITH_CHAIN=6
MIN_VISITS_PER_WWTP=2

def sha(data):return hashlib.sha256(data).hexdigest()
def avg(v):return float(statistics.mean(v)) if v else None
def median(v):return float(statistics.median(v)) if v else None
def parse(v):
    s=str(v).strip()
    if not s or s.upper() in ("NA","N/A","NAN"):return None
    try:x=float(s)
    except (ValueError,TypeError):return None
    return x if math.isfinite(x) and x>=0 else None

def source_frame(path):
    data=path.read_bytes()
    if sha(data)!=ORIGINAL_SHA256:
        raise RuntimeError("Original published Borealis cross-stage source SHA256 changed")
    rows=list(csv.DictReader(io.StringIO(data.decode("utf-8-sig")),delimiter="\t"))
    if len(rows)!=N_SOURCE:
        raise RuntimeError("Original true 906 source specimens source record drift")
    for q in ("Site","Location","Time","Life_Stage","Taxa","MP_per_mg","Mass_mg_wet","Rep"):
        if q not in rows[0]:raise RuntimeError("Missing published same-stage key "+q)
    group=collections.defaultdict(list)
    for x in rows:
        site=x["Site"].strip()
        loc=x["Location"].strip().lower()
        time=x["Time"].strip()
        taxon=x["Taxa"].strip()
        if site not in SITES or loc not in ("upstream","downstream") or time not in TIMES:
            raise RuntimeError("Original source experimental site×time/up-down keys inconsistent")
        group[(site,time,loc,taxon)].append(x)
    if set(x["Site"] for x in rows)!=set(SITES):
        raise RuntimeError("Original eight distinct treatment plants not present")
    # IMPORTANT: count physical specimens, before opening MP_per_mg.
    all_groups={**PRIMARY,**SECONDARY}
    ncell={}
    for site,time in itertools.product(SITES,TIMES):
        ncell[(site,time)]={}
        for name,ta in all_groups.items():
            nup=len(group[(site,time,"upstream",ta)])
            ndn=len(group[(site,time,"downstream",ta)])
            ncell[(site,time)][name]={"n_upstream":nup,"n_downstream":ndn,
                 "paired_two_replicates_each":min(nup,ndn)>=MIN_REPLICATES_PER_SIDE}
    qualify={}
    for genus,stages in (("Hydropsychidae",("Hydropsychidae_larvae","Hydropsychidae_adults","Tetragnathidae_spider")),
                         ("Chironomidae",("Chironomidae_larvae","Chironomidae_adults","Tetragnathidae_spider")),
                         ("Heptageniidae",("Heptageniidae_larvae","Heptageniidae_adults","Tetragnathidae_spider"))):
        visits={site:[] for site in SITES}
        for site,time in itertools.product(SITES,TIMES):
            if all(ncell[(site,time)][q]["paired_two_replicates_each"] for q in stages):
                visits[site].append(time)
        qualify[genus]={
          "site_has_at_least_two_source_eligible_original_visits":{
              s:len(v)>=MIN_VISITS_PER_WWTP for s,v in visits.items()},
          "n_sites_two_visit_chain":sum(len(v)>=MIN_VISITS_PER_WWTP for v in visits.values()),
          "n_site_time_all_three_stages_matched":sum(map(len,visits.values())),
          "original_source_matched_visits_by_facility":visits,
          "source_eligibility_pass":sum(len(v)>=MIN_VISITS_PER_WWTP for v in visits.values())>=MIN_INDEPENDENT_WWTP_WITH_CHAIN}
    return rows,group,ncell,qualify,sha(data)

def original_numeric(group,site,time,taxon):
    data={}
    for loc in ("upstream","downstream"):
        r=group[(site,time,loc,taxon)]
        vals=[v for x in r if (v:=parse(x["MP_per_mg"])) is not None]
        data[loc]={"n_original_samples":len(r),
                   "n_original_finite_nonnegative_MP_per_mg":len(vals),
                   "MP_per_mg_mean":avg(vals),
                   "MP_per_mg_median":median(vals),
                   "n_true_zero_source_measurements":sum(x==0 for x in vals)}
    good=min(data[q]["n_original_finite_nonnegative_MP_per_mg"] for q in data)>=MIN_REPLICATES_PER_SIDE
    if good:
        d=data["downstream"]["MP_per_mg_mean"]-data["upstream"]["MP_per_mg_mean"]
    else:d=None
    return {"site":site,"time":time,"taxon":taxon,
            "original_unit":"microplastic_particles_per_mg_wet_sample",
            "original_upstream":data["upstream"],"original_downstream":data["downstream"],
            "matched_original_at_least_two_finite_samples_each_side":good,
            "downstream_minus_upstream_MP_per_mg_mean":d}

def main_analysis(group,ncell,qualify):
    source_pass=qualify["Hydropsychidae"]["source_eligibility_pass"] and qualify["Chironomidae"]["source_eligibility_pass"]
    if not source_pass:return None
    # The entire ORIGINAL valid source cohort selected by source n, not signs.
    allmetrics={}
    for label,taxon in {**PRIMARY,**SECONDARY}.items():
        matches=[]
        for site,time in itertools.product(SITES,TIMES):
            if ncell[(site,time)][label]["paired_two_replicates_each"]:
                x=original_numeric(group,site,time,taxon)
                if x["matched_original_at_least_two_finite_samples_each_side"]:
                    matches.append(x)
        facilities={}
        for site in SITES:
            cells=[x for x in matches if x["site"]==site]
            z=[x["downstream_minus_upstream_MP_per_mg_mean"] for x in cells]
            up=avg([x["original_upstream"]["MP_per_mg_mean"] for x in cells])
            dn=avg([x["original_downstream"]["MP_per_mg_mean"] for x in cells])
            facilities[site]={"n_source_matched_visits":len(z),
                              "original_visits":[x["time"] for x in cells],
                              "source_upstream_per_mg_equal_visit_mean":up,
                              "source_downstream_per_mg_equal_visit_mean":dn,
                              "downstream_to_upstream_ratio_if_both_positive":dn/up if up is not None and up>0 and dn is not None else None,
                              "equal_time_average_downstream_minus_upstream_MP_per_mg":avg(z)}
        a=[r["equal_time_average_downstream_minus_upstream_MP_per_mg"] for r in facilities.values()
           if r["n_source_matched_visits"]>=MIN_VISITS_PER_WWTP]
        positive=sum(x>0 for x in a)
        negative=sum(x<0 for x in a)
        neutral=sum(x==0 for x in a)
        # The eight facilities share one Grand River network. Do not use a
        # random independent-WWTP treatment p-value or call this nutrient flux.
        allmetrics[label]={
          "observed_original_taxon":taxon,"original_matched_site_time_n":len(matches),
          "original_8WWTP_facility_site_summaries":facilities,
          "original_facilities_with_at_least_two_valid_visits":len(a),
          "source_min_six_original_facilities_two_visits_measurement_gate_pass":len(a)>=MIN_INDEPENDENT_WWTP_WITH_CHAIN,
          "source_eligibility_status":"MEASUREMENT_SUPPORT_GO" if len(a)>=MIN_INDEPENDENT_WWTP_WITH_CHAIN else "STOP_INSUFFICIENT_NUMERIC_ASSAY_SUPPORT_NO_CONFIRMATORY_STAGE_EFFECT",
          "original_matched_site_time_aggregate_qa":[{"site":x["site"],"time":x["time"],
            "observed_n_up":x["original_upstream"]["n_original_finite_nonnegative_MP_per_mg"],
            "observed_n_down":x["original_downstream"]["n_original_finite_nonnegative_MP_per_mg"]} for x in matches],
          "original_sites_with_same_year_visit_UPDOWN_source_only":len([r for r in facilities.values() if r["n_source_matched_visits"]>0]),
          "site_equal_exposed_both_sides_mean_upstream_per_mg":avg([r["source_upstream_per_mg_equal_visit_mean"] for r in facilities.values() if r["n_source_matched_visits"]>=MIN_VISITS_PER_WWTP]),
          "site_equal_exposed_both_sides_mean_downstream_per_mg":avg([r["source_downstream_per_mg_equal_visit_mean"] for r in facilities.values() if r["n_source_matched_visits"]>=MIN_VISITS_PER_WWTP]),
          "facility_ratio_summary_only_if_original_upstream_positive":median([r["downstream_to_upstream_ratio_if_both_positive"] for r in facilities.values() if r["n_source_matched_visits"]>=MIN_VISITS_PER_WWTP and r["downstream_to_upstream_ratio_if_both_positive"] is not None]),
          "facility_equal_mean_updown_difference_MP_per_mg":avg(a),
          "facility_equal_median_updown_difference_MP_per_mg":median(a),
          "facility_positive_negative_zero_difference_counts":[positive,negative,neutral],
          "sample_level_replicates_not_used_as_independent_WWTP_treatments":True,
          "no_arbitrary_pseudocount_for_ratio_or_log_difference":True,
          "no_cross_taxon_unit_multiplication":True}
    # Site-level relation of adult insect and riparian spider downstream-upstream
    # contrasts is an opportunity coresponse, NOT demonstrated trophic feeding.
    linked={}
    for genus in ("Hydropsychidae","Chironomidae"):
        g=allmetrics[genus+"_adults"]["original_8WWTP_facility_site_summaries"]
        s=allmetrics["Tetragnathidae_spider"]["original_8WWTP_facility_site_summaries"]
        points=[]
        for site in SITES:
            if min(g[site]["n_source_matched_visits"],s[site]["n_source_matched_visits"])>=2:
                points.append({"site":site,"adult_downup":g[site]["equal_time_average_downstream_minus_upstream_MP_per_mg"],
                               "spider_downup":s[site]["equal_time_average_downstream_minus_upstream_MP_per_mg"]})
        linked[genus]={"n_same_facility_source_comparable_adult_spider":len(points),
            "paired_original_site_facility_adult_and_riparian_spider_contrasts":points,
            "no_prey_consumer_individual_feeding_observed":True}
    return {"all_original_taxonomic_stages":allmetrics,
            "same_WWTP_adult_vs_spider_copatterns":linked}

def run(args):
    rows,group,ncell,qualify,source_sha=source_frame(args.published)
    result={
      "schema":SCHEMA,
      "status":"SOURCE_MATCHED_8_WWTP_AQUATIC_ADULT_SPIDER_MATERIAL_OUTCOMES_EXPLORATORY" if
          qualify["Hydropsychidae"]["source_eligibility_pass"] and qualify["Chironomidae"]["source_eligibility_pass"] else
          "STOP_SITE_TIME_STAGE_SOURCE_INSUFFICIENT_NO_NUMERIC_ANALYSIS",
      "original_Dataverse_published_DOI":PUBLISHER,
      "published_study_DOI":PAPER,
      "published_original_file_id":ORIGINAL_DATAFILE_ID,
      "source_original_SHA256":source_sha,
      "original_source_animal_specimens_n":len(rows),
      "original_independent_WWTP_source_facilities_n":8,
      "n_independent_original_hydrologic_basins":1,
      "original_distinct_sample_times_per_facility":3,
      "original_source_site_year_not_explicit_in_cross_species_table":True,
      "original_effluent_source_year_cross_DOI_matching_NOT_assumed":True,
      "original_site_time_taxon_source_counts":[
        {"site":site,"time":time,"stages":ncell[(site,time)]}
        for site,time in itertools.product(SITES,TIMES)],
      "original_aquatic_larva_adult_spider_source_gate":qualify,
      "explicit_site_stage_SOURCE_FIRST_quality_before_numeric_magnitudes":True,
      "reported_nonnegative_concentration_per_wet_mg_not_river_emergence_biomass":True,
      "no_source_specific_N_or_P_nutrient_effect_or_biomass_claim":True,
      "published_data_reanalysis_not_independent_confirmatory_study":True}
    if result["status"].startswith("SOURCE_MATCHED_"):
        result["numerical_original_per_material_stage_contrasts"]=main_analysis(group,ncell,qualify)
    else:
        result["numerical_original_per_material_stage_contrasts"]=None
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"GrandRiver8_true_WWTP_larval_adult_spider_original_MP_source_gated.json").write_text(
       json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("GRAND_RIVER8_TRUE_SOURCE_MATERIAL="+json.dumps({
      "status":result["status"],"original_six_WWTP_source_EFFLUENT_DOI_NOT_FITTED":True,
      "larva_adult_spider_original_2022_2023_stage_source_sites":{
        k:v["n_sites_two_visit_chain"] for k,v in qualify.items()},
      "stage_mean_down_minus_up_MP_per_mg":{
        k:round(v["facility_equal_mean_updown_difference_MP_per_mg"],9)
        if v["facility_equal_mean_updown_difference_MP_per_mg"] is not None else None
        for k,v in (result["numerical_original_per_material_stage_contrasts"] or {}).get("all_original_taxonomic_stages",{}).items()}
      },ensure_ascii=False),flush=True)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--published",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
if __name__=="__main__":
    main()
