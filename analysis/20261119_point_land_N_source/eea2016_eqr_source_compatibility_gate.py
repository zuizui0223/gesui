#!/usr/bin/env python3
"""Pre-NEW-MODEL cohort check: old EQR-available 2016 NUTS2 vs point/soil N.

The original EQR archive was already opened for a previous study. Only
the frozen list of region IDs and cohort counts is used; no mean EQR values
are inspected. This is a posthoc cohort compatibility audit, not blinding.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from datetime import datetime,timezone

RESULT=Path("results/eea2016_thirdyear_source_20261111/2016_EEA_RIVER_BENTHIC_EQR_JRC_GREEN_MODEL_POINTNP_NUTS2_HELDOUT_RESULTS.json")
NEEDED=(55,8,5,8)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for c in iter(lambda:f.read(1<<20),b""):h.update(c)
    return h.hexdigest()

def evaluate(source_csv,gate_json):
    srcgate=json.loads(Path(gate_json).read_text())
    if srcgate["results"]["source_only_GO_for_aquatic_only_next_stage"] is not True:
        raise ValueError("First phase source-only 2016 opportunity STOP")
    if srcgate["results"]["source_complete_NUTS2_shared_n"]!=73:
        raise ValueError("Frozen original 73 source regions changed")
    with open(source_csv,newline="") as f:
        src=list(csv.DictReader(f))
    if len(src)!=73:
        raise ValueError("Original 2016 water-only fixed source geographic roster changed")
    bykey={z["country_NUTS2"]:z for z in src}
    if len(bykey)!=73:raise ValueError("Source country/NUTS2 duplicate")
    old=json.loads(RESULT.read_text())
    if old["results"]["benthic_invertebrates"]["precommitted_numeric_EQR_support_gate_PASS"] is not True:
        raise ValueError("Original source frozen EQR validity check was STOP")
    obs=old["results"]["benthic_invertebrates"]["original_NUTS2_disjoint_holdout_Ridge10"]["source_JRC_green_region_sampled_EEA_EQR_summary"]
    if len(obs)!=73:raise ValueError("Old valid EQR NUTS2 inventory changed")
    if len({z["nuts2"] for z in obs})!=73:
        raise ValueError("Original numeric EQR NUTS2 duplicate")
    overlap=[]
    numeric=set()
    for z in obs:
        key=z["country"]+"|"+z["nuts2"]
        numeric.add(key)
        if key not in bykey:continue
        s=bykey[key]
        if int(z["source_frozen_test_fold"])!=int(s["original_geo_fold"]):
            raise ValueError("Original site fold changed")
        overlap.append({
            "country_NUTS2":key,
            "country":z["country"],
            "nuts2":z["nuts2"],
            "original_geo_fold":int(z["source_frozen_test_fold"]),
            "original_fixed_source_site_n":int(s["frozen_EEA2016_benthic_source_site_n"]),
            "original_numeric_valid_site_n":int(z["n_original_2016_BQE_sites"])
        })
    by_country=defaultdict(list)
    folds=Counter()
    for z in overlap:
        source=bykey[z["country_NUTS2"]]
        by_country[z["country"]].append(source)
        folds[z["original_geo_fold"]]+=1
    varied=[]
    for country,rows in by_country.items():
        if len(rows)<3:continue
        p={round(float(z["original_point_N_model_t_per_km2"]),12) for z in rows}
        d={round(float(z["land_atmos_model_N"]),12) for z in rows}
        if len(p)>=2 and len(d)>=2:varied.append(country)
    checks={
      "numeric_valid_and_source_complete_NUTS2_ge_55":len(overlap)>=NEEDED[0],
      "numeric_valid_and_source_complete_countries_ge_8":len(by_country)>=NEEDED[1],
      "countries_3regions_and_P_D_variation_ge_5":len(varied)>=NEEDED[2],
      "all_source_frozen_folds_ge_8_regions":all(folds[i]>=NEEDED[3] for i in range(1,6))
    }
    result={
      "original_source_only_regions":len(src),
      "original_old_valid_EQR_regions":len(obs),
      "original_numeric_missing_from_soilN_model_source":sorted(numeric-set(bykey)),
      "original_soilN_model_source_missing_from_old_numeric_EQR":sorted(set(bykey)-numeric),
      "overlap_region_n":len(overlap),
      "overlap_country_n":len(by_country),
      "overlap_countries":sorted(by_country),
      "overlap_country_with_interregional_P_D_variation":sorted(varied),
      "overlap_fold_counts":{str(i):folds[i] for i in range(1,6)},
      "overlap_original_site_n_before_numeric":sum(z["original_fixed_source_site_n"] for z in overlap),
      "overlap_original_numeric_valid_site_n":sum(z["original_numeric_valid_site_n"] for z in overlap),
      "requirements":checks,
      "new_aquatic_only_four_model_numeric_GO":all(checks.values()),
      "not_a_new_preoutcome_EEA_numeric_cohort":True,
      "no_numerical_EQR_values_read_for_new_model_gate":True,
      "original_water_land_2018_17to19_2x2_SOURCE_STOP_persist":True,
      "region_ids":[z["country_NUTS2"] for z in sorted(overlap,key=lambda z:z["country_NUTS2"])],
    }
    return result

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--eligible-source",type=Path,required=True)
    p.add_argument("--source-gate",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    data=evaluate(a.eligible_source,a.source_gate)
    result={"schema":"gesui_EEA2016_NUTS2_original_eqr_source_overlap_before_new_model_v01",
            "UTC":datetime.now(timezone.utc).isoformat(),
            "original_2016_numeric_region_summary_SHA256":sha(RESULT),
            "fixed_2016_source_region_CSV_SHA256":sha(a.eligible_source),
            "2016_water_only_source_gate_SHA256":sha(a.source_gate),
            "results":data}
    (a.out/"EEA2016_POINT_LANDN_MODEL_PREMODEL_COMPATIBILITY_GATE.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_2016_PREMODEL_COMPATIBILITY_SOURCE_AUDIT",json.dumps(data,ensure_ascii=False),flush=True)

if __name__=="__main__":main()
