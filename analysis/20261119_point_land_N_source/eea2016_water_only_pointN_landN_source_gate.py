#!/usr/bin/env python3
"""Independent 2016 water-only source overlap opportunity; NO EQR values."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime, timezone

MIN_REGION=55
MIN_COUNTRY=8
MIN_WITHIN_COUNTRY_VARIATION=5
MIN_EACH_FOLD=8
FIELDS=("NUTS2_country_id","water_point_N_t","water_point_P_t",
        "land_atmos_N_mg_m2","land_mineral_N_kg_ha","land_organic_N_kg_ha",
        "land_fixation_N_kg_ha")

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()

def read_sources(csvfile,reportfile):
    report=json.loads(Path(reportfile).read_text())
    if report.get("NO_EEA_EQR_EIV_GBIF_BIOLOGICAL_VALUES_READ") is not True:
        raise ValueError("Original JRC source-only receipt missing")
    if report.get("csv_sha256")!=digest(csvfile):
        raise ValueError("Original 269 NUTS2 CSV SHA drift")
    with open(csvfile,newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f)
        if not set(FIELDS).issubset(reader.fieldnames):
            raise ValueError("Original model fields missing")
        rows=list(reader)
    src={}
    for r in rows:
        key=r["NUTS2_country_id"].upper()
        if key in src:raise ValueError("Duplicate model source region")
        fields={}
        for col in FIELDS[1:]:
            try:
                value=float(r[col])
            except (ValueError,TypeError):
                value=None
            fields[col]=value if value is not None and math.isfinite(value) and value>=0 else None
        if all(v is not None for v in fields.values()):
            src[key]=fields
    countries={k.split("|",1)[0] for k in src}
    if len(src)!=269 or len(countries)!=28:
        raise ValueError("Original 269 NUTS2/28 countries cohort changed")
    return src

def gate(registry, sources):
    j=json.loads(Path(registry).read_text())
    if j.get("original_2016_EEA_numeric_EQR_outcome_VALUES_not_read") is not True:
        raise ValueError("2016 source registry not outcome free")
    x=j["source_biological_classes"]["benthic_invertebrates"]
    ro=x["original_NUTS2_region_source_inventory"]
    folds=x["fixed_region_to_fivefold_source_only"]
    if len(ro)!=79 or x["n_original_independent_NUTS2"]!=79 or x["n_sites_source_fixed_before_EQR"]!=2469:
        raise ValueError("2016 source pre-frozen benthic cohort mismatch")
    if len(folds)!=79:raise ValueError("Original 2016 geographically disjoint fold definitions changed")
    overlap=[]
    total_n=0
    for r in ro:
        country=str(r["country"]).upper()
        nuts=str(r["original_GISCO_NUTS2"]).upper()
        key=country+"|"+nuts
        total_n+=int(r["source_original_num_2016_EEA_sites"])
        if nuts not in folds:raise ValueError("Original 2016 geography fold key absent")
        if key not in sources:continue
        source=sources[key]
        p=float(r["source_NUTS2_GREEN_pointN_t"])
        s=source["water_point_N_t"]
        if abs(p-s)>1e-8*max(1,abs(p),abs(s)):
            raise ValueError("2016 point N from independent JRC source disagrees")
        area=float(r["source_NUTS2_area_km2"])
        if not math.isfinite(area) or area<=0:raise ValueError("2016 native GISCO area invalid")
        overlap.append({
            "country":country,"NUTS2":nuts,"country_NUTS2":key,
            "frozen_EEA2016_benthic_source_site_n":int(r["source_original_num_2016_EEA_sites"]),
            "original_geo_fold":int(folds[nuts]),
            "original_point_N_model_t_per_region":s,
            "original_point_N_model_t_per_km2":s/area,
            "land_atmos_model_N":source["land_atmos_N_mg_m2"],
            "land_mineral_model_N":source["land_mineral_N_kg_ha"],
            "land_organic_model_N":source["land_organic_N_kg_ha"],
        })
    if total_n!=2469:raise ValueError("2016 original source site census changed")
    country=defaultdict(list)
    foldcounts=Counter()
    for z in overlap:
        country[z["country"]].append(z)
        foldcounts[z["original_geo_fold"]]+=1
    varied=[]
    for c,v in sorted(country.items()):
        if len(v)<3:continue
        point={round(z["original_point_N_model_t_per_km2"],12) for z in v}
        deposition={round(z["land_atmos_model_N"],12) for z in v}
        if len(point)>=2 and len(deposition)>=2:varied.append(c)
    checks={
        "shared_NUTS2_ge_55":len(overlap)>=MIN_REGION,
        "shared_countries_ge_8":len(country)>=MIN_COUNTRY,
        "countries_with_three_NUTS2_and_both_P_D_vary_ge_5":len(varied)>=MIN_WITHIN_COUNTRY_VARIATION,
        "all_original_source_folds_ge_8_NUTS2":all(foldcounts[i]>=MIN_EACH_FOLD for i in range(1,6))
    }
    support={
        "2016_original_benthic_NUTS2":79,
        "2016_original_country_n":13,
        "2016_original_benthic_source_site_n":2469,
        "source_complete_NUTS2_shared_n":len(overlap),
        "source_shared_countries":sorted(country),
        "source_shared_country_n":len(country),
        "source_varied_within_country_list":varied,
        "source_varied_within_country_count":len(varied),
        "frozen_NUTS2_fold_region_counts":{str(i):foldcounts[i] for i in range(1,6)},
        "shared_EEA2016_original_benthic_site_records":sum(z["frozen_EEA2016_benthic_source_site_n"] for z in overlap),
        "source_only_gate_items":checks,
        "source_only_GO_for_aquatic_only_next_stage":all(checks.values()),
        "never_executed_biotic_numeric_model":True,
        "does_not_reopen_2018_2017to19_water_plant_2x2_stop":True,
    }
    return overlap,support

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--registry",type=Path,required=True)
    p.add_argument("--source-csv",type=Path,required=True)
    p.add_argument("--source-report",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    src=read_sources(args.source_csv,args.source_report)
    overlap,support=gate(args.registry,src)
    obj={
        "schema":"gesui_2016_aquatic_only_NUTS2_point_landN_source_gate_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2016_EEA_benthic_no_EQR_registry_SHA256":digest(args.registry),
        "original_JRC269_source_CSV_SHA256":digest(args.source_csv),
        "results":support,
        "not_an_ecological_effect_estimate":True,
        "no_EEA_numeric_EQR_read":True,
    }
    (args.out/"EEA2016_WATER_ONLY_POINT_LANDN_SOURCE_GATE.json").write_text(
        json.dumps(obj,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    with (args.out/"EEA2016_WATER_ONLY_POINT_LANDN_ELIGIBLE_NUTS2.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(overlap[0]) if overlap else
                        ["country","NUTS2","country_NUTS2"])
        w.writeheader()
        w.writerows(overlap)
    print("GESUI_2016_AQUATIC_ONLY_SOURCE_GATE",json.dumps(support,ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":main()
