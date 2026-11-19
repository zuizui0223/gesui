#!/usr/bin/env python3
"""Source-only GO/STOP for all original EIV-quality H7 × 5 JRC native layers.

The older NUTS2 water x land biotic 2x2 remains STOP. This script never loads
EQR, NUTRI or HNS. It makes no biological effects or significance claims.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
from datetime import datetime,timezone
from statistics import mean

SHARDS=12
ALL_H7=226
SOURCE_TYPES=("atmos2015","mineralN","organicN","fixationN","pointN")
CORE=("atmos2015","mineralN","organicN","pointN")
BASIN_SINGLE_MIN=90.
MIN_UNIQUE_NATIVE_BASINS=100
MIN_H7=100
MIN_H6=40
MIN_COUNTRIES=8
MIN_VARIATION_COUNTRIES=5
MIN_COMPONENTS_PER_FOLD=4
N_FOLDS=5
ORIGINAL_PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for d in iter(lambda:f.read(1<<20),b""):h.update(d)
    return h.hexdigest()

def load_all(shards):
    records=[]
    manifest=[]
    identical=None
    for idx in range(SHARDS):
        candidates=list(Path(shards).rglob(f"JRC_NATIVE_SOURCE_H7_SHARD_{idx:02d}_OF_12.json"))
        if len(candidates)!=1:
            raise ValueError(f"Original H7 SOURCE shard {idx} missing or duplicated")
        file=candidates[0]
        j=json.loads(file.read_text())
        if j.get("schema")!="gesui_JRC_originalH7_2017to19_native_source_226_full_sharded_v1":
            raise ValueError("Source shard schema changed")
        if j.get("shard")!=idx or j.get("n_shards")!=SHARDS or j.get("original_total_H7_EIV_eligible_first_period")!=ALL_H7:
            raise ValueError("Original source shard fixed allocation or total cohort changed")
        if j.get("no_EQR_NUTRI_HNS_biological_outcome_numeric_read") is not True:
            raise ValueError("Source shard biological values read unexpectedly")
        shas=(j["original_H7_original_panel_SHA256"],j["original_2808_facets_SHA256"],
              j["EIV_451_trait_quality_eligibility_SHA256"])
        if shas[0]!=ORIGINAL_PANEL_SHA or (identical is not None and shas!=identical):
            raise ValueError("Original source cohort SHA not identical across all shards")
        identical=shas
        entries=j.get("source_rows",[])
        if len(entries)!=j["number_model_cells"] or len(entries)!=j["H7_processed"]*len(SOURCE_TYPES):
            raise ValueError("Original source report incomplete inside shard")
        records.extend(entries)
        manifest.append({"shard":idx,"source_file_sha256":sha(file),
                         "H7_processed":j["H7_processed"],
                         "status_counts":j["source_status_counts"]})
    combos={(r["h7"],r["source"]) for r in records}
    H7={x for x,_ in combos}
    if len(H7)!=ALL_H7 or len(records)!=ALL_H7*len(SOURCE_TYPES) or len(combos)!=len(records):
        raise ValueError("Original 226 H7 × five native sources not fully represented")
    return records,manifest,identical

class Partition:
    def __init__(self):self.parent={}
    def find(self,x):
        self.parent.setdefault(x,x)
        if x!=self.parent[x]:self.parent[x]=self.find(self.parent[x])
        return self.parent[x]
    def union(self,a,b):
        a,b=self.find(a),self.find(b)
        if a!=b:self.parent[max(a,b)]=min(a,b)

def source_support(records):
    grouped=defaultdict(dict)
    for r in records:
        if r["source"] not in SOURCE_TYPES:raise ValueError("Unknown JRC native source type")
        if r["source"] in grouped[r["h7"]]:raise ValueError("Duplicate H7 × source")
        grouped[r["h7"]][r["source"]]=r
    eligible={}
    technical_error=[]
    rows=[]
    for h in sorted(grouped):
        a=grouped[h]
        if set(a)!=set(SOURCE_TYPES):raise ValueError("Original five H7 source cells missing")
        if len({str(z["true_h6"]) for z in a.values()})!=1:
            raise ValueError("Original true H6 differs within source cells")
        if len({str(z["original_source_country"]) for z in a.values()})!=1:
            raise ValueError("Original source country differs within source cells")
        if any(z.get("status")!="COMPLETE_OFFICIAL_JRC_SOURCE_ONLY" for z in a.values()):
            technical_error.append(h)
        union_status=all(a[k].get("status")=="COMPLETE_OFFICIAL_JRC_SOURCE_ONLY" and
                         (a[k].get("numeric_geom_coverage_pct") or 0)>=BASIN_SINGLE_MIN
                         for k in CORE)
        strict_status=all(a[k].get("status")=="COMPLETE_OFFICIAL_JRC_SOURCE_ONLY" and
                          (a[k].get("max_numeric_single_native_coverage_pct") or 0)>=BASIN_SINGLE_MIN
                          for k in CORE)
        soil_ids=[a[k].get("max_single_native_id") for k in CORE[:3]]
        soil_same=(len(set(soil_ids))==1 and soil_ids[0] is not None)
        point_id=a["pointN"].get("max_single_native_id")
        # Do not assume that a polygon can be assigned by union coverage alone.
        strict=strict_status and soil_same and point_id is not None
        repaired=any(a[k].get("material_original_jrc_ring_repair",False) for k in CORE)
        # Native Esri rings are primary. If original rings need material repair,
        # an independent geometry validation is required before final GO.
        quality=strict and not repaired
        if quality:eligible[h]=a
        rows.append({
            "original_H7":h,"true_H6":str(a["atmos2015"]["true_h6"]),
            "original_source_country":str(a["atmos2015"]["original_source_country"]),
            "union_90_all4_core":union_status,
            "single_native_polygon90_all4_core":strict_status,
            "same_original_soil_HydroID_across_3_sources":soil_same,
            "single_model_source_strict_geom_eligible":strict,
            "original_native_rings_unrepaired_or_nonmaterial":not repaired,
            "conservative_geometry_quality_pass":quality,
            "soil_native_source_ID":soil_ids[0] if soil_same else None,
            "water_native_point_BasinID":point_id,
            "point_model_N_density_proxy_t_km2":a["pointN"].get("model_pointN_t_per_native_area_km2_proxy"),
            "atmos_model_N_mean":a["atmos2015"].get("soil_model_basin_weighted_mean"),
            "mineral_model_N_mean":a["mineralN"].get("soil_model_basin_weighted_mean"),
            "organic_model_N_mean":a["organicN"].get("soil_model_basin_weighted_mean"),
            "fixation_model_N_mean":a["fixationN"].get("soil_model_basin_weighted_mean"),
            "status_code_all5":";".join(k+":"+str(a[k]["status"]) for k in SOURCE_TYPES),
        })
    # Geographic components are on original H6 + distinct namespaced native
    # point source and distinct native SOIL HydroID, never numerical ID equality.
    uf=Partition()
    for h,a in eligible.items():
        h6="H6:"+str(a["atmos2015"]["true_h6"])
        soil="SOIL_HydroID:"+str(a["atmos2015"]["max_single_native_id"])
        point="WATER_BasinID:"+str(a["pointN"]["max_single_native_id"])
        uf.union(h6,soil);uf.union(h6,point)
    components=defaultdict(list)
    for h,a in eligible.items():
        root=uf.find("H6:"+str(a["atmos2015"]["true_h6"]))
        components[root].append(h)
    ordered=sorted(components.items(),key=lambda z:(-len(z[1]),z[0]))
    foldcounts=[0]*N_FOLDS
    fold_h7=[0]*N_FOLDS
    foldmap={}
    for group,hs in ordered:
        idx=min(range(N_FOLDS),key=lambda i:(fold_h7[i],foldcounts[i],i))
        foldmap[group]=idx+1
        fold_h7[idx]+=len(hs)
        foldcounts[idx]+=1
    for row in rows:
        if not row["conservative_geometry_quality_pass"]:
            row["locked_geo_fold"]=None
        else:
            node="H6:"+row["true_H6"]
            row["locked_geo_fold"]=foldmap[uf.find(node)]
    valid=[z for z in rows if z["conservative_geometry_quality_pass"]]
    unique_soil={z["soil_native_source_ID"] for z in valid}
    unique_point={z["water_native_point_BasinID"] for z in valid}
    uniqueH6={z["true_H6"] for z in valid}
    countries={z["original_source_country"] for z in valid}
    varied={}
    for country in sorted(countries):
        local=[z for z in valid if z["original_source_country"]==country]
        indep_soil={z["soil_native_source_ID"] for z in local}
        indep_point={z["water_native_point_BasinID"] for z in local}
        p={round(float(z["point_model_N_density_proxy_t_km2"]),12) for z in local
           if z["point_model_N_density_proxy_t_km2"] is not None}
        d={round(float(z["atmos_model_N_mean"]),12) for z in local
           if z["atmos_model_N_mean"] is not None}
        varied[country]={"H7":len(local),"unique_soil_native_basins":len(indep_soil),
                         "unique_point_native_basins":len(indep_point),
                         "distinct_pointN_density":len(p),"distinct_atmosN":len(d),
                         "qualifies_for_source_within_country_gradient":len(indep_soil)>=3 and len(indep_point)>=3 and len(p)>1 and len(d)>1}
    varied_n=sum(x["qualifies_for_source_within_country_gradient"] for x in varied.values())
    n_union=sum(z["union_90_all4_core"] for z in rows)
    n_strict=sum(z["single_model_source_strict_geom_eligible"] for z in rows)
    checks={
      "original_all_226_H7_5sources_fully_retrieved_without_technical_stop":len(technical_error)==0,
      "at_least_100_distinct_original_soil_native_basins":len(unique_soil)>=MIN_UNIQUE_NATIVE_BASINS,
      "at_least_100_conservatively_geometry_qualified_H7":len(valid)>=MIN_H7,
      "at_least_40_genuine_H6":len(uniqueH6)>=MIN_H6,
      "at_least_8_original_source_country_labels":len(countries)>=MIN_COUNTRIES,
      "at_least_5_countries_with_3_independent_native_basins_and_source_variation":varied_n>=MIN_VARIATION_COUNTRIES,
      "all_fivefolds_at_least_4_independent_connected_components":all(x>=MIN_COMPONENTS_PER_FOLD for x in foldcounts),
    }
    diagnostics={
       "original_full_H7_cohort":len(rows),
       "original_5_source_cell_status_counts":dict(Counter(x["status"] for a in grouped.values() for x in a.values())),
       "H7_with_all4_model_union_source_coverage_at_least90pct":n_union,
       "H7_with_all4_single_native_polygon_90pct_and_same_soil_id":n_strict,
       "H7_conservative_geometry_qualified":len(valid),
       "H7_missing_any_core_source_or_single_basins":len(rows)-n_strict,
       "H7_material_native_Esri_ring_repair_unverified":sum(not z["original_native_rings_unrepaired_or_nonmaterial"] for z in rows),
       "number_native_soil_HydroID_independent":len(unique_soil),
       "number_native_water_point_BasinID_independent":len(unique_point),
       "true_H6_count":len(uniqueH6),"original_country_count":len(countries),
       "original_country_source_counts":dict(Counter(x["original_source_country"] for x in valid)),
       "countries_with_multiple_original_native_sources_and_variation":varied,
       "countries_source_variation_gate_n":varied_n,
       "independent_H6_x_soil_native_x_point_native_connected_components":len(components),
       "geo_components_by_fixed_fold":foldcounts,
       "eligible_H7_by_fixed_fold":fold_h7,
       "source_only_checks":checks,
       "all_hard_source_checks_pass":all(checks.values()),
       "new_biotic_model_permitted_by_this_source_gate":all(checks.values()),
       "original_2018_waterxland_NUTS2_2x2_SOURCE_STOP_remains":True,
       "no_2018_2022_EEQ_OR_PLANT_EIV_NUTRI_BUTTERFLY_HNS_VALUES_READ":True,
       "interpretation":"Model source coverage and geometry only; not ecological effect, delivered dose or independent experiments.",
    }
    return rows,diagnostics

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--shards",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    records,manifest,cohort=load_all(args.shards)
    rows,gate=source_support(records)
    outcome={
      "schema":"gesui_20261119_original_226_H7_JRCnative_all5_source_only_strict_go_stop_v1",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "original_H7_cohort_SHA256":cohort[0],
      "original_facets_SHA256":cohort[1],
      "original_EIV_quality_eligibility_SHA256":cohort[2],
      "shard_receipts":manifest,
      "gate":gate,
      "strict_original_single_native_basin_90pct_unless_false_not_overridden":True,
      "modelled_raster_derived_basin_mean_is_not_original_pixelwise_reaggregation":True,
      "older_NUTS2_2018water_2017land_2x2_stop_unchanged":True
    }
    report=args.out/"JRC_226_H7_NATIVE_SOURCE_ONLY_GO_STOP.json"
    report.write_text(json.dumps(outcome,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    with (args.out/"JRC_226_H7_SOURCE_COVERAGE_AND_FOLD_METADATA_ONLY.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),extrasaction="raise")
        writer.writeheader();writer.writerows(rows)
    print("GESUI_JRC_FULL_226_H7_SOURCE_ONLY_GATE",json.dumps({
        "original_H7":len(rows),
        "union_90":gate["H7_with_all4_model_union_source_coverage_at_least90pct"],
        "strict_single_90":gate["H7_with_all4_single_native_polygon_90pct_and_same_soil_id"],
        "strict_geom_qualified_H7":gate["H7_conservative_geometry_qualified"],
        "independent_soil_basins":gate["number_native_soil_HydroID_independent"],
        "H6":gate["true_H6_count"],
        "countries":gate["original_country_count"],
        "countries_with_source_variation":gate["countries_source_variation_gate_n"],
        "components":gate["independent_H6_x_soil_native_x_point_native_connected_components"],
        "components_by_fold":gate["geo_components_by_fixed_fold"],
        "hard_source_pass":gate["all_hard_source_checks_pass"],
        "checks":gate["source_only_checks"],
    },ensure_ascii=False,sort_keys=True),flush=True)

if __name__=="__main__":main()
