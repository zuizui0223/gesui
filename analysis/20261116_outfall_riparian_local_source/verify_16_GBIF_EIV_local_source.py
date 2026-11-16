#!/usr/bin/env python3
"""Source-only original 16-outfall riparian plant EIV support + flowpath distance gate.

No river-plant effect model, EIV numeric values, butterfly outcomes or outfall
nutrient effects read. Separately prove (i) potential original EIV scoreable
GBIF species support at upstream U2/U1/downstream D1 source-fixed circles,
(ii) U1-to-D1 original *river-path* distance supports an actual local exposure
contrast, and (iii) original source records contain actual same-year measured
WWTP annual nutrient load. This pilot source panel was selected before GBIF
species traits, but source gate thresholds were chosen *posthoc after earlier
plant richness analyses* and are not a prospective preregistered test.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path

GBIF_REGISTRY_SHA="0180bbc9be4ce3dcc5aa0c5eb092110e1b182b387f01095c0cf53eaf956f7d9f"
GBIF_FACET_SHA="ab549e0a8e2ddecd6afa5dbaf12180da64166b8f579c3916c46094da62d5c959"
GBIF_ORIGINAL_PRETRT_SHA="7654dd23a1646573d92972967b6ed0ce44cf4be6d50e04ab4f52778eea6a766f"
EIV_11511_SHA="ab45f36c0dad047442bb958bcef3135b45014bafec8a70811bcafb6109a02679"
TRT_SCHEMA="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
REPORT_SCHEMA="GESUI_20261116_ORIGINAL_16_OUTFALL_48_SITE_DIRECTED_EIV_SOURCE_LOCAL_EXPOSURE_GATE_V1"

# Explicit posthoc feasibility, NOT a pre-registered effect test:
SPECIES_MIN=10
SPECIES_COVERAGE_MIN=.40
SCORED_RECORD_MIN=20
SOURCE_LOCAL_KM=2.0
SOURCE_LOOSE_KM=5.0
REQUIRED_FACILITIES=10
REQUIRED_COUNTRIES=3
RIVER_ROLES=("u2","u1","d1")

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def locate(folder,name):
    files=list(folder.rglob(name))
    if len(files)!=1:raise RuntimeError("Required exactly one byte-authentic source item "+name)
    return files[0]
def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def keys(x):
    return str(x).strip()

def main(args):
    args.output.mkdir(parents=True,exist_ok=True)
    registry=locate(args.gbif,"plantae_site_registry.csv")
    facets=locate(args.gbif,"plantae_site_species_records.csv")
    pre=locate(args.gbif,"preplant_sites_and_estimands.json")
    eiv=locate(args.eiv,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv")
    qa=locate(args.eiv,"EU454_FULL_EIV_SOURCE_GATE.json")
    for p,known in [(registry,GBIF_REGISTRY_SHA),(facets,GBIF_FACET_SHA),
                    (pre,GBIF_ORIGINAL_PRETRT_SHA),(eiv,EIV_11511_SHA)]:
        if sha(p)!=known:raise RuntimeError("Original exact source bytes changed: "+p.name)
    g=json.loads(qa.read_text())
    if (g.get("schema")!=TRT_SCHEMA or not g.get("hard_source_support_pass")
        or not g.get("original_EIV_NUTRI_scored_keys_unchanged")
        or g.get("source_hybrid_score_assigned")!=0):
        raise RuntimeError("Source-original EIV strict eligibility not verified")
    e=rows(eiv);site=rows(registry);facet=rows(facets)
    if len(e)!=11511 or len(site)!=48 or len(facet)!=5394:
        raise RuntimeError("Source-original trait/site species rowcount changed")
    if len({q["speciesKey"] for q in e})!=11511:
        raise RuntimeError("Official source EIV speciesKey duplicate")
    scored=set()
    for q in e:
        if q["trait_status"]=="EXACT_MATCH" and q["EIV_NUTRI_exact"].strip()!="":
            scored.add(q["speciesKey"])
        elif q["EIV_NUTRI_exact"].strip():
            raise RuntimeError("Unexpected numeric score on EIV nonmatched/hybrid species")
    if len(scored)!=4011:
        raise RuntimeError("Exact 4011 accepted vascular numeric EIV source changed")
    pp=json.loads(pre.read_text())
    if pp.get("radius_km")!=1 or pp.get("site_units") is None or len(pp["site_units"])!=48:
        raise RuntimeError("Original source selected site registry differs from original response-blind design")
    if pp.get("year")!=[2011,2024] and pp.get("year")!= "2011,2024":
        # Original schema can store [2011,2024]; no conversion allowed.
        raise RuntimeError("Original 14-year GBIF source opportunity period changed: "+str(pp.get("year")))
    sources={q["site_id"]:q for q in site}
    if len(sources)!=48 or {q["river_role"] for q in site}!=set(RIVER_ROLES):
        raise RuntimeError("Original 48 identified study sites missing roles")
    by=defaultdict(lambda:{"source_scored_keys":set(),"source_all_keys":set(),
                           "total_specieskey_records":0,"scoreable_specieskey_records":0})
    for x in facet:
        name=x["site_id"]
        if name not in sources:raise RuntimeError("Original GBIF source facet site not registered")
        n=int(x["speciesKey_record_count"])
        if n<1:raise RuntimeError("No GBIF count for original source speciesKey")
        k=x["gbif_speciesKey"]
        if k in by[name]["source_all_keys"]:raise RuntimeError("Source original duplicate speciesKey in site facet")
        by[name]["source_all_keys"].add(k)
        by[name]["total_specieskey_records"]+=n
        if k in scored:
            by[name]["source_scored_keys"].add(k)
            by[name]["scoreable_specieskey_records"]+=n
    report=[]
    fac=defaultdict(list)
    for q in site:
        k=q["site_id"];x=by[k]
        scored_species=len(x["source_scored_keys"])
        all_species=len(x["source_all_keys"])
        total=x["total_specieskey_records"];n=x["scoreable_specieskey_records"]
        if all_species!=int(q["distinct_species_keys"]) or total!=int(q["indexed_species_record_count"]):
            raise RuntimeError("Original registered site/facet total GBIF speciesKey counts inconsistent")
        cov=scored_species/all_species if all_species else 0
        good=(scored_species>=SPECIES_MIN and cov>=SPECIES_COVERAGE_MIN
              and n>=SCORED_RECORD_MIN)
        row={"wwtp_id":q["wwtp_id"],"country":q["country"],"river_role":q["river_role"],
             "site_id":k,"raw_original_GBIF_species":all_species,
             "source_exact_EIV_scoreable_species":scored_species,
             "fraction_original_species_with_numeric_EIV":cov,
             "original_speciesKey_record_count":total,
             "records_of_EIV_scored_species":n,
             "source_EIV_quality_pass":good,
             "source_U1_D1_river_path_km":float(q["distance_u1_d1_km"]),
             "source_U2_U1_river_path_km":float(q["distance_u2_u1_km"]),
             "source_original_model_effluent_fraction":float(q["source_fraction"]),
             "original_source_reported_annual_N_mass_present":bool(q["observed_N_mass_kg_year"].strip()),
             "original_model_WASTE_DIS_is_not_annual_NP_report":True,
             "GBIF_2011_to_2024_not_contemporary_2022_or_2023":True}
        report.append(row);fac[row["wwtp_id"]].append(row)
    if len(fac)!=16 or any(set(r["river_role"] for r in q)!=set(RIVER_ROLES) for q in fac.values()):
        raise RuntimeError("Original source 16 freshwater outfall units do not have U2 U1 D1")
    goodfac=[q for q in fac.values() if all(x["source_EIV_quality_pass"] for x in q)]
    short2=[q for q in goodfac if q[0]["source_U1_D1_river_path_km"]<=SOURCE_LOCAL_KM]
    short5=[q for q in goodfac if q[0]["source_U1_D1_river_path_km"]<=SOURCE_LOOSE_KM]
    allfac=[q[0] for q in fac.values()]
    n_input_mass=sum(q["original_source_reported_annual_N_mass_present"] for q in allfac)
    all_distance=[q["source_U1_D1_river_path_km"] for q in allfac]
    if len(all_distance)!=16 or len(goodfac)!=13 or len(short2)!=0 or len(short5)!=2 or n_input_mass!=0:
        raise RuntimeError("Unexpected original source-only 16-outfall trait/distance/annual-load counts")
    if min(all_distance)<=0 or max(all_distance)>100:
        raise RuntimeError("Original river path/source distances invalid")
    result={
       "schema":REPORT_SCHEMA,
       "status":"SOURCE_EIV_TRAIT_GEOGRAPHIC_PASS_BUT_LOCAL_DISCHARGE_EXPOSURE_STOP",
       "posthoc_source_feasibility_not_confirmatory_ecology":True,
       "original_16_HydroWASTE_v1_study_facility_count":16,
       "original_fixed_GBIF_site_count":48,
       "original_GBIF_source_speciesKey_facets":len(facet),
       "original_EIV_crosswalk_speciestaxa":len(e),
       "original_EIV_exact_accepted_vascular_species":len(scored),
       "original_source_trait_site_quality_thresholds":{
         "minimum_scored_species_per_site":SPECIES_MIN,
         "minimum_fraction_of_species_EIV_scored":SPECIES_COVERAGE_MIN,
         "minimum_speciesKey_record_counts_of_scored_species":SCORED_RECORD_MIN},
       "individual_site_EIV_scoreable":sum(x["source_EIV_quality_pass"] for x in report),
       "independent_outfalls_all3_sites_EIV_scoreable":len(goodfac),
       "countries_in_all3_qualified_outfalls":len({q[0]["country"] for q in goodfac}),
       "original_source_U1_D1_river_path_km":{
         "min":min(all_distance),"max":max(all_distance),
         "median":float(sorted(all_distance)[7]+sorted(all_distance)[8])/2},
       "all16_outfalls_U1_D1_within_2km":sum(x<=2. for x in all_distance),
       "all16_outfalls_U1_D1_within_5km":sum(x<=5. for x in all_distance),
       "EIV_scored_triplets_U1_D1_within_2km":len(short2),
       "EIV_scored_triplets_U1_D1_within_5km":len(short5),
       "original_source_site_annual_reported_N_mass_present_outfall_n":n_input_mass,
       "EIV_source_coverage_gate_pass":len(goodfac)>=REQUIRED_FACILITIES and
          len({q[0]["country"] for q in goodfac})>=REQUIRED_COUNTRIES,
       "strict_local_same_river_path_2km_gate_pass":len(short2)>=REQUIRED_FACILITIES,
       "annual_source_wwtp_NP_delivery_at_site_gate_pass":False,
       "local_wastewater_effect_model_fitted":False,
       "GBIF_site_centres_are_reach_centroids_not_physically_outfall_proximal":True,
       "source_GBIF_window_14years_wider_than_water_2018_2022":True,
       "third_party_converted_HydroRIVERS_not_official_topology_byte_verified":True,
       "source_provenance_exact_sha256":{
         "GBIF_48_site_registry":GBIF_REGISTRY_SHA,
         "GBIF_5394_speciesKey_site_facets":GBIF_FACET_SHA,
         "original_outcome_blind_48_site_design":GBIF_ORIGINAL_PRETRT_SHA,
         "official_EIV_11511_species_crosswalk":EIV_11511_SHA}}
    out=args.output/"ORIGINAL_16_OUTFALL_EIV_SOURCE_LOCAL_EXPOSURE_GATE.json"
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    with (args.output/"ORIGINAL_16_OUTFALL_48_SITE_EIV_SOURCE_OPPORTUNITY.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(report[0]))
        w.writeheader();w.writerows(sorted(report,key=lambda q:(q["wwtp_id"],q["river_role"])))
    print("GESUI_ORIGINAL_16_OUTFALL_LOCAL_EIV_GATE="+json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--gbif",required=True,type=Path)
    p.add_argument("--eiv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
