#!/usr/bin/env python3
"""Independent built-in-Python audit of 36 new response-blind matched GBIF sites.

Recompute exact individual-count rarefaction, sample qualification and exact
paired sign tests from source-frozen EEA measured P station pairs and original
GBIF all-Plantae speciesKey facet pages. Zero network, zero SciPy dependencies.
"""
from __future__ import annotations
from collections import Counter
from pathlib import Path
import argparse,hashlib,json,math,statistics

SOURCE_RUN=42377235104
GBIF_RUN=42377454765
PRIMARY_DEPTH=30
SECONDARY_DEPTH=10

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def expected(counts,depth):
    vals=[int(x) for x in counts if int(x)>0]
    n=sum(vals)
    if n<depth:return None
    def choose(n,k):
        return math.lgamma(n+1)-math.lgamma(k+1)-math.lgamma(n-k+1)
    denom=choose(n,depth)
    return float(sum(1-math.exp(choose(n-i,depth)-denom) if n-i>=depth else 1
                     for i in vals))

def same(a,b,eps=1e-9):
    if not math.isclose(float(a),float(b),rel_tol=0,abs_tol=eps):
        raise ValueError("Independent source facet numerical discrepancy: "+str((a,b)))

def sign_exact(p,n):
    if n==0:return None
    return min(1.,2*sum(math.comb(n,j) for j in range(min(p,n-p)+1))/2**n)

def audit(source,receipt,rawpath,resultpath,gatepath):
    orig=json.loads(source.read_text())
    rec=json.loads(receipt.read_text())
    raw=json.loads(rawpath.read_text())
    result=json.loads(resultpath.read_text())
    frozen=json.loads(gatepath.read_text())
    src_sha=sha(source)
    if src_sha!=rec["original_files_SHA256"]["SOURCE_MATCHED_PAIR_SELECTION.json"]:
        raise ValueError("Source selected station set changed from historical frozen run")
    if raw["original_design_sha256"]!=src_sha or frozen["original_source_design_sha256"]!=src_sha:
        raise ValueError("Original raw plant recording not based on frozen source cohort")
    if raw["source_selection_run"]!=SOURCE_RUN or result["source_only_selection_original_run_id"]!=SOURCE_RUN:
        raise ValueError("Original no-outcome source run mismatch")
    pairs=orig["pairs"]
    recorded=raw["plants"]
    if not (len(pairs)==18 and len(recorded)==36 and result["distinct_new_GBIF_plant_sites"]==36):
        raise ValueError("Original sample was altered or missing GBIF sites")
    siteids=[x["station_id"] for x in recorded]
    if len(siteids)!=len(set(siteids)):raise ValueError("Duplicated GBIF site")
    lookup={s["station_id"]:s for s in recorded}
    root_pairs=result["paired"]
    out_by={x["pair_id"]:x for x in root_pairs}
    if len(out_by)!=18:raise ValueError("Pair table changed")
    all_species=set()
    total_source_records=0
    no_site_records=0
    results={"primary":[],"sensitivity":[]}
    seen_sites=set()
    for row in recorded:
        if row["gbif_years"]!="2022,2024" or row["site_radius_km"]!=1:
            raise ValueError("Year or radius changed")
        counts=row["specieskey_counts"]
        if len(counts)!=row["distinct_species_keys"]:
            raise ValueError("GBIF species-key count mismatch")
        if sum(int(v) for v in counts.values())!=row["specieskey_resolved_records"]:
            raise ValueError("Original GBIF species-keyed record denominator changed")
        if row["all_Plantae_records"]<row["specieskey_resolved_records"]:
            raise ValueError("Species-key records exceed Plantae records")
        if row["all_Plantae_records"]==0:no_site_records+=1
        total_source_records+=row["all_Plantae_records"]
        all_species.update(counts)
    if len(all_species)!=result["specieskeys_union"] or len(all_species)!=1006:
        raise ValueError("Original 1006 species keys not reproduced")
    for pair in pairs:
        id=pair["pair_id"]
        stored=out_by[id]
        ep=pair["source_exposed"]
        ctrl=pair["source_control"]
        if ep["eea2022_reported_measured_upstream_P_ton_yr"]<=0 or ctrl["eea2022_reported_measured_upstream_P_ton_yr"]!=0:
            raise ValueError("EEA source positive versus zero membership changed")
        if ep["station_id"]==ctrl["station_id"] or ep["block"]!=ctrl["block"]:
            raise ValueError("Source matching group invalid")
        seen_sites|={ep["station_id"],ctrl["station_id"]}
        a,b=lookup[ep["station_id"]],lookup[ctrl["station_id"]]
        for tag,depth,min_sp,balance in (("primary",PRIMARY_DEPTH,6,True),
                                         ("sensitivity",SECONDARY_DEPTH,4,False)):
            original=stored["primary_depth30" if tag=="primary" else "sensitivity_depth10"]
            q=all(x["specieskey_resolved_records"]>=depth and
                  x["distinct_species_keys"]>=min_sp for x in (a,b))
            aa=a["specieskey_resolved_records"];bb=b["specieskey_resolved_records"]
            if balance:q=q and min(aa,bb)>0 and max(aa,bb)/min(aa,bb)<=10
            if bool(original["eligible"])!=q:
                raise ValueError("Original eligibility gate not reproduced for "+id)
            if q:
                sa=expected(a["specieskey_counts"].values(),depth)
                sb=expected(b["specieskey_counts"].values(),depth)
                val=math.log1p(sa)-math.log1p(sb)
                same(val,original["richness_log_difference_source_minus_control"])
                results[tag].append((id,pair["spatial_block"],val))
            else:
                if original["richness_log_difference_source_minus_control"] is not None:
                    raise ValueError("Unqualified source contrast was not missing")
    if seen_sites!=set(siteids):raise ValueError("GBIF sites do not cover all frozen source pairs")
    independent={}
    for tag,field in (("primary","primary"),("sensitivity","sensitivity")):
        vals=[x[2] for x in results[tag]]
        saved=result[field]
        eligible=len(vals);pos=sum(v>0 for v in vals);neg=sum(v<0 for v in vals)
        n=pos+neg
        blocks=len(set(x[1] for x in results[tag]))
        support=(eligible>=10 and blocks>=4)
        if eligible!=saved["n_eligible_pairs"] or blocks!=saved["eligible_spatial_blocks"]:
            raise ValueError("Eligible pair denominator changed")
        if pos!=saved["positive_pairs"] or neg!=saved["negative_pairs"]:
            raise ValueError("Positive/negative pair signs changed")
        if eligible:
            same(statistics.mean(vals),saved["mean_log_richness_source_minus_control"])
            same(statistics.median(vals),saved["median_log_richness_source_minus_control"])
        if support:
            same(sign_exact(pos,n),saved["exact_two_sided_sign_test_descriptive_p"])
        else:
            assert saved["exact_two_sided_sign_test_descriptive_p"] is None
        independent[tag]={"source_valid_pairs":eligible,"spatial_blocks":blocks,
             "positive_pairs":pos,"negative_pairs":neg,
             "mean_log_richness_source_minus_control":statistics.mean(vals) if vals else None,
             "exact_sign_p":sign_exact(pos,n) if support else None,
             "scientific_statistical_gate_pass":support}
    return {"status":"PASS_INDEPENDENT_NO_NETWORK_NO_SCIPY_SOURCE_RECOMPUTATION",
       "frozen_original_source_run":SOURCE_RUN,"original_GBIF_run":GBIF_RUN,
       "source_design_SHA256":src_sha,
       "raw_source_GBIF_SHA256":sha(rawpath),
       "original_results_SHA256":sha(resultpath),
       "original_fixed_outcome_gates_SHA256":sha(gatepath),
       "source_matched_pairs":len(pairs),"new_plant_site_n":len(recorded),
       "gbif_plant_species_key_n":len(all_species),
       "total_all_Plantae_GBIF_occurrence_records":total_source_records,
       "plant_sites_no_records":no_site_records,
       "primary_and_secondary_recomputed":independent,
       "source_2022_EEA_P_load_biologically_causal":False,
       "no_live_source_or_GBIF_queries":True}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--receipt",type=Path,required=True)
    p.add_argument("--raw",type=Path,required=True)
    p.add_argument("--result",type=Path,required=True)
    p.add_argument("--gates",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    obj=audit(a.source,a.receipt,a.raw,a.result,a.gates)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_NEW_EEA_P_PLANT_BLIND_AUDIT",json.dumps(obj,sort_keys=True),flush=True)
