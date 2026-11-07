#!/usr/bin/env python3
"""Audit 48-site GBIF river study using only versioned source JSON and stdlib.

This does not download GBIF or refit a post-hoc ecological model. It checks
source-byte identity, site/taxon units, rarefaction, negative-control contrasts,
Jaccard turnover, and Spearman rho independently of the original SciPy script.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/"results/directed_river_plantae_20261107"
GEOMETRY=ROOT/"results/hydrorivers_eu_source_gate_20261107"

def check_close(a,b,eps=1e-10):
    if a is None or b is None:
        assert a is None and b is None, (a,b)
    else:
        assert math.isclose(float(a),float(b),abs_tol=eps,rel_tol=0),(a,b)

def rank_avg(values):
    indexed=sorted(enumerate(values),key=lambda x:x[1])
    out=[0.0]*len(values)
    i=0
    while i<len(indexed):
        j=i+1
        while j<len(indexed) and indexed[i][1]==indexed[j][1]:
            j+=1
        for idx,_ in indexed[i:j]:
            out[idx]=(i+1+j)/2
        i=j
    return out

def spearman(x,y):
    assert len(x)==len(y)>2
    r1,r2=rank_avg(x),rank_avg(y)
    a=statistics.mean(r1);b=statistics.mean(r2)
    numerator=sum((v-a)*(w-b) for v,w in zip(r1,r2))
    denominator=math.sqrt(sum((v-a)**2 for v in r1)*sum((w-b)**2 for w in r2))
    return numerator/denominator

def rarefy(counts,depth):
    a=[int(n) for n in counts if int(n)>0]
    N=sum(a)
    if N<depth:return None
    def lc(n,k):
        return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
    d=lc(N,depth)
    return sum(1-math.exp(lc(N-v,depth)-d) if N-v>=depth else 1 for v in a)

def jaccard(keys1,keys2):
    a,b=set(keys1),set(keys2)
    if not (a|b):return None
    return 1-len(a&b)/len(a|b)

def audit(result_dir:Path=RESULTS):
    receipt=json.loads((result_dir/"SOURCE_RECEIPT.json").read_text())
    assert receipt["original_success_run"]==42365305227
    assert receipt["physical_source_run"]==42364866919
    for fname,expected in receipt["all_original_files_sha256"].items():
        digest=hashlib.sha256((result_dir/fname).read_bytes()).hexdigest()
        assert digest==expected,(fname,digest,expected)
    saved=json.loads((result_dir/"directional_results.json").read_text())
    raw=json.loads((result_dir/"plantae_raw_species_facets.json").read_text())
    frozen=json.loads((result_dir/"preplant_sites_and_estimands.json").read_text())
    geo=json.loads((GEOMETRY/"candidate_site_design.json").read_text())
    assert frozen["source_geometry_original_run"]==42364866919
    assert frozen["geometry_sha256"]==hashlib.sha256((GEOMETRY/"candidate_site_design.json").read_bytes()).hexdigest()
    assert frozen["qc_sha256"]==hashlib.sha256((GEOMETRY/"source_qc.json").read_bytes()).hexdigest()
    assert geo["metadata"]["sample_selected_without_GBIF_response"] is True
    assert len(geo["selected"])==16 and len(frozen["site_units"])==48
    assert raw["schema"]=="gesui_source_fixed_directed_reach_gbif_plantae_facets_v01"
    assert len(raw["facets"])==48 and saved["site_count"]==48
    assert all(x["observed_N_mass_kg_year"] is None for x in frozen["site_units"])
    sites={x["site_id"]:x for x in raw["facets"]}
    assert len(sites)==48
    per_site_facet={}
    for site_id,v in sites.items():
        species={str(k):int(n) for k,n in v["species_key_record_counts"].items()}
        assert len(species)==v["distinct_species_keys"]
        assert sum(species.values())==v["indexed_species_record_count"]
        assert v["total_plantae_occurrence_records"]>=v["indexed_species_record_count"]
        per_site_facet[site_id]=species
    observed_rows=set()
    with (result_dir/"plantae_site_species_records.csv").open(encoding="utf8",newline="") as f:
        for row in csv.DictReader(f):
            key=(row["site_id"],row["gbif_speciesKey"])
            assert key not in observed_rows
            observed_rows.add(key)
            assert row["gbif_year_range"]=="2011,2024"
            assert int(row["source_run_id"])==42364866919
            assert per_site_facet[key[0]][key[1]]==int(row["speciesKey_record_count"])
    assert len(observed_rows)==receipt["site_species_table_rows"]==5394
    all_keys={key[1] for key in observed_rows}
    assert len(all_keys)==receipt["unique_gbif_speciesKeys"]==1402
    registry={}
    with (result_dir/"plantae_site_registry.csv").open(encoding="utf8",newline="") as f:
        for row in csv.DictReader(f):
            assert row["site_id"] not in registry
            registry[row["site_id"]]=row
            assert int(row["total_plantae_occurrence_records"])==sites[row["site_id"]]["total_plantae_occurrence_records"]
    assert len(registry)==48
    saved_by={x["wwtp_id"]:x for x in saved["wwtp_outcomes"]}
    assert len(saved_by)==16
    primary=[]
    low_depth=[]
    for outfall in geo["selected"]:
        wid=outfall["wwtp_id"]
        item=saved_by[wid]
        s=[sites[str(wid)+"_"+suffix] for suffix in ("u2","u1","d1")]
        rec=[x["indexed_species_record_count"] for x in s]
        ratio=max(rec)/min(rec) if min(rec)>0 else float("inf")
        primary_q=all(x["indexed_species_record_count"]>=30 and x["distinct_species_keys"]>=6 for x in s)
        low_q=all(x["indexed_species_record_count"]>=10 and x["distinct_species_keys"]>=4 for x in s)
        assert primary_q==item["depth30"]["eligible_by_records_and_species"]
        assert low_q==item["depth10_sensitivity"]["eligible_by_records_and_species"]
        assert (primary_q and ratio<=10)==item["strict_primary_eligible"]
        for depth,eligible,label in ((30,primary_q,"depth30"),(10,low_q,"depth10_sensitivity")):
            if not eligible:
                assert item[label]["extra_change_downstream_minus_prior_upstream"] is None
                continue
            S=[rarefy(x["species_key_record_counts"].values(),depth) for x in s]
            early=math.log1p(S[1])-math.log1p(S[0])
            late=math.log1p(S[2])-math.log1p(S[1])
            extra=late-early
            obs=item[label]
            for a,b in zip(S,obs["rarefied_expected_plant_species_u2_u1_d1"]):
                check_close(a,b)
            check_close(early,obs["pre_outfall_upstream_log_richness_change"])
            check_close(late,obs["across_outfall_downstream_log_richness_change"])
            check_close(extra,obs["extra_change_downstream_minus_prior_upstream"])
            j_up=jaccard(s[0]["species_key_record_counts"],s[1]["species_key_record_counts"])
            j_down=jaccard(s[1]["species_key_record_counts"],s[2]["species_key_record_counts"])
            check_close(j_up,obs["jaccard_dissimilarity_u2_u1"])
            check_close(j_down,obs["jaccard_dissimilarity_u1_d1"])
            check_close(j_down-j_up,obs["turnover_increment_across_minus_upstream"])
        if item["strict_primary_eligible"]:primary.append(item)
        if low_q:low_depth.append(item)
    assert len(primary)==saved["primary_depth30_and_balance_n"]==12
    assert len(low_depth)==saved["sensitivity_depth10_n"]==14
    for group,field,depth in ((primary,"primary_source_fraction_association","depth30"),
                              (low_depth,"depth10_association_sensitivity","depth10_sensitivity")):
        x=[v["source_fraction"] for v in group]
        y=[v[depth]["extra_change_downstream_minus_prior_upstream"] for v in group]
        rho=spearman(x,y)
        check_close(rho,saved[field]["spearman_rho"])
    assert saved["status"]=="OBSERVATIONAL_GBIF_ALL_PLANTAE_DIRECTED_REACH_NOT_WASTEWATER_EFFECT"
    return {
        "status":"PASS_INDEPENDENT_STDLIB_RECALCULATION",
        "original_real_GBIF_run":receipt["original_success_run"],
        "verified_source_files":len(receipt["all_original_files_sha256"]),
        "GBIF_site_registry":len(registry),
        "site_by_species_records":len(observed_rows),
        "unique_species_keys":len(all_keys),
        "primary_source_eligible_outfalls":len(primary),
        "sensitivity_depth10_outfalls":len(low_depth),
        "primary_rho":spearman([x["source_fraction"] for x in primary],
                               [x["depth30"]["extra_change_downstream_minus_prior_upstream"] for x in primary]),
        "external_network_calls":0,
    }

if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,sort_keys=True,indent=2))
