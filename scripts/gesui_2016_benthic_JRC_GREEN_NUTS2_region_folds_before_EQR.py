#!/usr/bin/env python3
"""Freeze independent NUTS2-disjoint validation folds before 2016 EQR values.

Source cohort: 2016 site metadata (2016 EEA river Biology official),
2016 GISCO NUTS2 assignment, official JRC 2014-18 modeled domestic/industry
point N/P. No 2016 EQR values loaded. 5 fixed folds assign entire NUTS2
regions, avoiding same spatial exposure region entering train and test.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

SOURCE=Path("results/eea2016_thirdyear_source_20261111/2016_EEA_BENTHIC_GISCO_NUTS2_GREEN_POINT_NP_SOURCE_GATE_BEFORE_EQR.json")
BQE=("benthic_invertebrates","phytobenthos","fish")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run(root):
    r=json.loads((root/SOURCE).read_text())
    gate=r["precommitted_2016_benthic_NP_modelled_NUTS2_SUPPORT_GATE"]
    if not gate["can_open_2016_EEA_numeric_normalised_EQR_and_test_POSTHOC_correlation"]:
        raise ValueError("Independent source-only 2016 NUTS2 benthic sample gate failed; no future EQR")
    if r["all_original_2016_EEA_BQE_EQR_results_NOT_READ"] is not True:
        raise ValueError("2016 source sample selected after biological response")
    orig=r["original_2016_source_site_and_modelled_NUTS2_PN_values_HASHED_IDs_NO_EQR"]
    if len(orig)<4000:raise ValueError("Original GISCO/JRC site source subset unexpectedly small")
    res={}
    for bio in BQE:
        rows=[x for x in orig if x["original_2016_BQE_group"]==bio]
        key=[(x["original_2016_BQE_site_key_SHA256"],x["original_2016_BQE_group"]) for x in rows]
        if len(set(key))!=len(key):raise ValueError("Duplicated original 2016 BQE monitoring site×element")
        regions=Counter(x["original_GISCO2016_NUTS2"] for x in rows)
        allocate={}
        loads=[0]*5
        bins=[[] for _ in range(5)]
        for region,n in sorted(regions.items(),key=lambda z:(-z[1],z[0])):
            j=min(range(5),key=lambda i:(loads[i],i))
            allocate[region]=j+1
            bins[j].append(region)
            loads[j]+=n
        if sum(loads)!=len(rows):raise ValueError("Cross-validation site count differs source cohort")
        if len(allocate)!=len(regions) or len(allocate)<10:raise ValueError("Original source region fold manifest corrupted")
        # Freeze region-level contribution and EEA assessment methods only.
        byregion=[]
        for region,n in sorted(regions.items()):
            sample=[x for x in rows if x["original_GISCO2016_NUTS2"]==region]
            byregion.append({
                "original_GISCO_NUTS2":region,
                "country":sample[0]["original_2016_EEA_country"],
                "source_original_num_2016_EEA_sites":n,
                "original_fixed_NUTS2_disjoint_fold":allocate[region],
                "source_NUTS2_GREEN_pointN_t":sample[0]["original_JRC2014to18_modelled_PsN_t"],
                "source_NUTS2_GREEN_pointP_t":sample[0]["original_JRC2014to18_modelled_PsP_t"],
                "source_NUTS2_area_km2":sample[0]["original_GISCO2016_NUTS2_km2"]
            })
            if any(x["original_JRC2014to18_modelled_PsN_t"]!=sample[0]["original_JRC2014to18_modelled_PsN_t"] or
                   x["original_JRC2014to18_modelled_PsP_t"]!=sample[0]["original_JRC2014to18_modelled_PsP_t"]
                   for x in sample):
                raise ValueError("NUTS2 regional original modeled N/P source inconsistent within source site")
        res[bio]={
            "n_sites_source_fixed_before_EQR":len(rows),
            "n_original_independent_NUTS2":len(regions),
            "fixed_region_to_fivefold_source_only":allocate,
            "original_5fold_test_site_counts":{str(i+1):n for i,n in enumerate(loads)},
            "original_fivefold_test_regions":{str(i+1):bins[i] for i in range(5)},
            "original_reported_country_count":len({x["original_2016_EEA_country"] for x in rows}),
            "original_NUTS2_region_source_inventory":byregion,
            "source_BQE_EQR_2016_value_not_accessed":True
        }
    return {
        "schema":"gesui_2016_EEA_BQE_benthic_phytobenthos_GISCO2016_NUTS2_disjoint_fivefold_YEAR2016_POSTHOC_noEQR_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_GISCO2016_JRC2014to18_official_source_full_site_manifest_SHA256":sha(root/SOURCE),
        "source_biological_classes":res,
        "fixed_alpha_future_Ridge":10,
        "primary_BQE":"benthic_invertebrates",
        "secondary_BQE":"phytobenthos",
        "fish_min5country_gate_failed":True,
        "numeric_EQR_future_quality_gate":{
            "min_primary_2016_valid_numeric_EQR_monitoring_sites":800,
            "min_primary_2016_with_valid_EQR_countries":8,
            "min_primary_independent_NUTS2_regions":60,
            "min_primary_3degree_geographic_blocks":15,
            "min_secondary_valid_numeric_sites":600,
            "min_secondary_numeric_country_count":8,
            "min_secondary_independent_NUTS2":40
        },
        "fixed_NUTS2_leakage_exclusion":True,
        "original_2016_EEA_numeric_EQR_outcome_VALUES_not_read":True,
        "modelled_wastewater_and_industrial_POINT_NP_not_WWTP_actually_reported2016_mass":True,
        "no_causal_pollution_effect_or_species_diversity_promoted":True
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n")
    print(json.dumps({k:{"source_sites":v["n_sites_source_fixed_before_EQR"],
      "independent_NUTS2":v["n_original_independent_NUTS2"],
      "CV_fold_size":v["original_5fold_test_site_counts"],
      "country_n":v["original_reported_country_count"]}
       for k,v in x["source_biological_classes"].items()},indent=2))
if __name__=="__main__":main()
