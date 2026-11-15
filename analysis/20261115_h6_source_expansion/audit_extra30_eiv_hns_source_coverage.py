#!/usr/bin/env python3
"""STRICT SOURCE QUALITY ONLY for 30 genuinely new H7 sites x two periods.

Inspect archived original 360 GBIF speciesKey facets, the original preexisting
4011 plant EIV NUTRI exact-score keys, and the old 303 butterfly HNS roster
(186 accepted/scored). No new trait imputation, original W, regression, EQR,
butterfly abundance or causal connection.

This source-only gate is NOT the same original 454 panel. 10 new H6 is merely
an independent geographic opportunity if every requisite source check passes.
"""
from __future__ import annotations
import argparse
import csv
import json
import hashlib
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd

FACET_GROUPS=("Plantae","Hesperiidae","Lycaenidae","Nymphalidae","Papilionidae","Pieridae")
PERIODS=("2017-2019","2020-2022")
EXPECTED_NUMERIC_EIV=4011
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
HNS_SHA="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
WKT_SHA="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
TAXA_SHA="afc524ebfa353e4d1c6ea51d208ba529a85fedb57c4d0f28d83545010bb76189"
PRIOR_SCHEMA="GESUI_ORIGINAL_1219_GBIF_COUNT_SOURCE_OPPORTUNITY_FOR_12_NEW_SINGLETON_H6_V1"
SOURCE_SCHEMA="GESUI_SOURCE_ONLY_30_NEW_H7_2ERA_SIX_ORIGINAL_FACET_ACQUISITION_V1"
MIN_PLANT_N=10
MIN_PLANT_SOURCE_FRAC=.4
MIN_HNS_SOURCE_RECORDS=10
MIN_HNS_SPECIES=2
MIN_HNS_RAW_RECORD_FRAC=.7
MIN_NEW_H6=10

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def one(d,name):
    a=list(d.rglob(name))
    if len(a)!=1:raise RuntimeError("Expected one exact source artifact file "+name+": "+repr([str(x) for x in a]))
    return a[0]
def readcsv(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f))
def writecsv(p,rows):
    if not rows:return
    with p.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)

def main(a):
    meta=json.loads(one(a.prior,"ORIGINAL1219_GBIF_NEW_H7_RECORD_OPPORTUNITY_GATE.json").read_text())
    if meta.get("schema")!=PRIOR_SCHEMA or not meta.get("original_fixed_40_H6_target_source_index_count_gate_pass"):
        raise RuntimeError("Original 40 H7 source count opportunity not authenticated")
    previous=readcsv(one(a.prior,"ORIGINAL1219_EEA_NEW_H7_PLANT_LYCA_ANIMALIA_RECORD_OPPORTUNITIES.csv"))
    if len(previous)!=40:raise RuntimeError("Original source eligible H7 registry drift")
    qualified={r["native_H7"]:r for r in previous if r["both_periods_original_three_taxa_index_tiers_pass"]=="True"}
    if len(qualified)!=30 or len({r["official_native_H6"] for r in qualified.values()})!=10:
        raise RuntimeError("Thirty true source-count supported new H7 in ten H6 changed")
    v2=json.loads(one(a.v2,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (v2.get("schema")!="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or not v2.get("hard_source_support_pass") or
        v2.get("source_EIV",{}).get("original_source_SHA256")!=EIV_SHA or
        v2.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175 or
        v2.get("source_hybrid_score_assigned")!=0):
        raise RuntimeError("Former 4011 numeric EIV source-v2 hash/valid hybrid QA missing")
    cross=pd.read_csv(one(a.v2,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"),dtype={"speciesKey":str})
    if len(cross)!=11511 or cross.speciesKey.duplicated().any():
        raise RuntimeError("Original 11511 source plant roster changed")
    score_eiv=set(cross.loc[cross.EIV_NUTRI_exact.notna(),"speciesKey"])
    if len(score_eiv)!=EXPECTED_NUMERIC_EIV:
        raise RuntimeError("Original verified EIV scored accepted speciesKey count mismatch")
    if sha(a.hns)!=HNS_SHA:raise RuntimeError("Original 303 butterfly HNS score table SHA mismatch")
    taxa=pd.read_csv(a.hns,encoding="utf-8-sig",dtype={"butterflyKey":str})
    if len(taxa)!=303 or taxa.butterflyKey.duplicated().any():
        raise RuntimeError("Old 303-source HNS accepted butterfly speciesKey roster altered")
    alias_bad=taxa.actual_used_alias_conflict.astype(str).str.lower().isin(["true","1","yes"])
    allowed=taxa[(taxa.taxon_status.astype(str)=="ACCEPTED")&
                 (~alias_bad)&taxa.host_nutrient_mean.notna()]
    scored_HNS={str(row.butterflyKey):str(row.family) for row in allowed.itertuples(index=False)}
    if len(scored_HNS)!=186:
        raise RuntimeError("Original accepted host N trait only has 186 scored butterfly species")
    source={}
    shard_audit=[]
    for i in range(3):
        d=a.shards/str(i)
        qa=one(d,f"EXTRA_EEA_H7_REAL_SPECIESKEY_full_SHARD_{i:02d}_QA.json")
        q=json.loads(qa.read_text())
        if (q.get("schema")!=SOURCE_SCHEMA or q.get("source_only_mode")!="full" or
            q.get("source_shard")!=i or q.get("original_source_1219WKT_SHA256")!=WKT_SHA or
            q.get("original_234_fixed_family_taxa_source_SHA256")!=TAXA_SHA):
            raise RuntimeError("New genuine GBIF speciesKey shard provenance incorrect")
        datafile=one(d,f"EXTRA_EEA_H7_REAL_SPECIESKEY_full_SHARD_{i:02d}.json")
        if sha(datafile)!=q.get("original_output_source_JSON_SHA256"):
            raise RuntimeError("New 30 original GBIF speciesKey archive SHA mismatch")
        rows=json.loads(datafile.read_text())
        if len(rows)!=q["source_6_group_2period_facets_requested"] or len(rows)!=q["source_H7_in_shard"]*12:
            raise RuntimeError("New source GBIF shard missed original request opportunities")
        shard_audit.append({"shard":i,"n_H7":q["source_H7_in_shard"],
                            "facets":len(rows),"missing":q["source_facet_missing"],
                            "original_archive_SHA256":sha(datafile)})
        for r in rows:
            key=(str(r["HYBAS_ID_L07"]),str(r["period"]),str(r["group"]))
            if key in source or key[0] not in qualified or key[1] not in PERIODS or key[2] not in FACET_GROUPS:
                raise RuntimeError("New original GBIF speciesKey source contains invalid/duplicate H7-period-family key")
            if r["source_status"]=="ok":
                if not isinstance(r["source_speciesKey_counts"],dict) or not isinstance(r["source_query_occurrence_count"],int):
                    raise RuntimeError("Original source failed or lacks true counts")
                if sum(int(n) for n in r["source_speciesKey_counts"].values())>r["source_query_occurrence_count"]:
                    raise RuntimeError("Original indexed species count exceeds queried source record total")
            elif r["source_status"]=="missing":
                if r["source_speciesKey_counts"] is not None:
                    raise RuntimeError("Missing source was incorrectly coded as a biological zero")
            else:raise RuntimeError("Unknown original GBIF API facet status")
            source[key]=r
    if len(source)!=30*2*6 or sum(x["n_H7"] for x in shard_audit)!=30:
        raise RuntimeError("Full original 30×2×6=360 source facet matrix not complete")
    period_reports=[]
    all_plantae={}
    source_missing=[]
    for h in sorted(qualified,key=int):
        for period in PERIODS:
            facets={family:source[(h,period,family)] for family in FACET_GROUPS}
            if any(r["source_status"]!="ok" for r in facets.values()):
                failed=[g for g,x in facets.items() if x["source_status"]!="ok"]
                source_missing.extend([[h,period,g] for g in failed])
                period_reports.append({
                   "h":h,"h6":qualified[h]["official_native_H6"],"period":period,
                   "all_facets_complete":False,"source_failures":";".join(failed),
                   "known_eiv_plant_scored_n":None,"known_EIV_coverage_of_all_indexed_plants":None,
                   "HNS_scored_species_n":None,"HNS_scored_record_coverage_of_fivefamily_raw":None,
                   "both_plant_and_butterfly_trait_source_support":False})
                continue
            plant=facets["Plantae"]
            pl_counts={str(k):int(v) for k,v in plant["source_speciesKey_counts"].items() if int(v)>0}
            all_plantae[(h,period)]=set(pl_counts)
            pl_scored=set(pl_counts)&score_eiv
            pl_cov=len(pl_scored)/len(pl_counts) if pl_counts else 0.
            raw5=sum(facets[g]["source_query_occurrence_count"] for g in FACET_GROUPS if g!="Plantae")
            fam5=Counter()
            seen_scored=set()
            source5=0
            scored_records=0
            for fam in FACET_GROUPS:
                if fam=="Plantae":continue
                for k,n in facets[fam]["source_speciesKey_counts"].items():
                    if int(n)<0:raise RuntimeError("Negative native butterfly source indexed record count")
                    if k in scored_HNS:
                        if scored_HNS[k]!=fam:raise RuntimeError("Original accepted HNS speciesKey family mismatch")
                        if k in seen_scored:raise RuntimeError("Duplicate original numerical-HNS speciesKey across families")
                        seen_scored.add(k);scored_records+=int(n)
                    source5+=int(n)
                    fam5[fam]+=int(n)
            butterfly_coverage=scored_records/raw5 if raw5 else 0.
            ppass=(plant["source_query_occurrence_count"]>=50 and
                len(pl_scored)>=MIN_PLANT_N and pl_cov>=MIN_PLANT_SOURCE_FRAC)
            bpass=(scored_records>=MIN_HNS_SOURCE_RECORDS and
                len(seen_scored)>=MIN_HNS_SPECIES and
                butterfly_coverage>=MIN_HNS_RAW_RECORD_FRAC)
            period_reports.append({
               "h":h,"h6":qualified[h]["official_native_H6"],"period":period,
               "all_facets_complete":True,"source_failures":"",
               "original_Plant_record_index_count":plant["source_query_occurrence_count"],
               "source_all_Plant_speciesKeys":len(pl_counts),
               "known_eiv_plant_scored_n":len(pl_scored),
               "known_EIV_coverage_of_all_indexed_plants":pl_cov,
               "source_5butterfly_families_all_raw_records":raw5,
               "source_5butterfly_indexed_speciesKey_records":source5,
               "HNS_scored_species_n":len(seen_scored),
               "HNS_scored_records":scored_records,
               "HNS_scored_record_coverage_of_fivefamily_raw":butterfly_coverage,
               "plant_EIV_known_original_4011_subset_source_gate_pass":ppass,
               "butterfly_HNS_old186_score_source_gate_pass":bpass,
               "both_plant_and_butterfly_trait_source_support":ppass and bpass})
    bysite=defaultdict(dict)
    for r in period_reports:bysite[r["h"]][r["period"]]=r
    if len(bysite)!=30 or any(len(periods)!=2 for periods in bysite.values()):
        raise RuntimeError("Thirty new original H7×two periods missing after source trait QA")
    h7_good={h for h,q in bysite.items() if
             all(q[p]["both_plant_and_butterfly_trait_source_support"] for p in PERIODS)}
    h6_good={qualified[h]["official_native_H6"] for h in h7_good}
    unknown_plant_keys=set().union(*all_plantae.values()) if all_plantae else set()
    eiv_uncovered=unknown_plant_keys-set(cross.speciesKey)
    scored_matches=len(unknown_plant_keys&score_eiv)
    full_h6=30+len(h6_good)
    report={
       "schema":"GESUI_SOURCE_ONLY_EXTRA30_H7_BOTH_PERIOD_EIV_HNS_TRAIT_SUPPORT_V1",
       "original_source_eiv_sha":EIV_SHA,"original_HNS_roster_sha":HNS_SHA,
       "original_1219_WKT_SHA":WKT_SHA,
       "source_original_count_candidate_30H7_in_10_H6":True,
       "source_original_speciesKey_facets_expected":360,
       "source_original_speciesKey_facets_received":len(source),
       "original_source_missing_facets":len(source_missing),
       "original_source_missing_facet_ids":source_missing,
       "original_source_shard_provenance":shard_audit,
       "source_two_period_H7_with_facet_source_complete":sum(all(q[p]["all_facets_complete"] for p in PERIODS) for q in bysite.values()),
       "source_two_period_H7_with_plant_nutri_original_4011_lowerbound_and_HNS186_eligibility":len(h7_good),
       "source_distinct_new_H6_with_both_period_trait_support":len(h6_good),
       "source_possible_H6_total_existing30_plus_newtrait":full_h6,
       "original_source_min_new_H6_for_40":MIN_NEW_H6,
       "source_full_strict_40_H6_opportunity_gate_with_old_traits":len(h6_good)>=MIN_NEW_H6,
       "source_all_new_H7_plant_unique_keys":len(unknown_plant_keys),
       "original_previous_11511_plants_accepted_taxonomy_union":len(cross),
       "source_new_plant_speciesKeys_not_in_old_11511_taxonomy":len(eiv_uncovered),
       "source_new_plant_speciesKeys_numeric_NUTRI_scored_in_old_4011":scored_matches,
       "unscored_newplant_keys_are_unknown_not_absent":True,
       "source_new_H7_selected_bothperiod_trait_supported_ids":sorted(h7_good,key=int),
       "source_new_H6_with_old_trait_gates":sorted(h6_good,key=int),
       "source_thresholds":{"EIV_numeric_plant_species_min":MIN_PLANT_N,
            "EIV_allPlantae_speciesKey_fraction_min":MIN_PLANT_SOURCE_FRAC,
            "HNS_scored_record_min":MIN_HNS_SOURCE_RECORDS,
            "HNS_scored_distinct_species_min":MIN_HNS_SPECIES,
            "HNS_scored_fraction_of_fivefamily_raw_min":MIN_HNS_RAW_RECORD_FRAC},
       "source_trait_completeness_not_outcome_effect":True,
       "no_real_NUTRI_HNS_regression_or_EEA_EQR_or_effluent_N_dose_fitted":True,
       "original_study_454_covariates_not_source_matched_for_extra_H7":True
    }
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"EXTRA30_H7_SOURCE_TRAIT_QC.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    writecsv(a.output/"EXTRA30_H7_PERIOD_SOURCE_TRAIT_COVERAGE.csv",period_reports)
    print("GESUI_EXTRA30_H7_SOURCE_TRAIT_QA="+json.dumps(
        {k:v for k,v in report.items() if k!="original_source_shard_provenance"},
        ensure_ascii=False,sort_keys=True),flush=True)
    print("H6_TRAIT_COVERAGE_SOURCE_GATE="+("PASS_SOURCE_TRAITS_BUT_NEEDS_NEW_COVARIATES"
        if len(h6_good)>=MIN_NEW_H6 else "STOP_STRICT_40_H6_UNDER_EXISTING_TRAIT_SCORES"),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--prior",required=True,type=Path)
    p.add_argument("--shards",required=True,type=Path)
    p.add_argument("--v2",required=True,type=Path)
    p.add_argument("--hns",type=Path,default=Path("results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"))
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
