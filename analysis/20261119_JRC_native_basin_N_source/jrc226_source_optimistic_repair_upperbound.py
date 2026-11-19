#!/usr/bin/env python3
"""SOURCE-only optimistic bound: assume all native Esri ring repairs are valid.

An optimistic assumption can only INCREASE geographic support. If it still
fails the original independent JRC basin gate, do not run ecology models.
No EQR, NUTRI, HNS or other biological outcome is read.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path

EXPECTED_SCHEMA="gesui_20261119_original_226_H7_JRCnative_all5_source_only_strict_go_stop_v1"

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):h.update(x)
    return h.hexdigest()

class Union:
    def __init__(self):self.p={}
    def find(self,x):
        self.p.setdefault(x,x)
        if self.p[x]!=x:self.p[x]=self.find(self.p[x])
        return self.p[x]
    def join(self,x,y):
        x,y=self.find(x),self.find(y)
        if x!=y:self.p[max(x,y)]=min(x,y)

def summarize(rows):
    strict=[r for r in rows if r["single_model_source_strict_geom_eligible"]=="True"]
    full=[r for r in rows if r["conservative_geometry_quality_pass"]=="True"]
    if len(rows)!=226 or len({r["original_H7"] for r in rows})!=226:
        raise ValueError("Original source-only H7 population mismatch")
    if not set(r["original_H7"] for r in full).issubset(r["original_H7"] for r in strict):
        raise ValueError("Optimistic pre-geometry cohort not a superset of strict original gate")
    def blocks(data):
        s={x["soil_native_source_ID"] for x in data}
        p={x["water_native_point_BasinID"] for x in data}
        h6={x["true_H6"] for x in data}
        countries=defaultdict(list)
        for r in data:
            if not r["soil_native_source_ID"] or not r["water_native_point_BasinID"]:
                raise ValueError("Missing original JRC native source basin ID despite single-basin eligibility")
            countries[r["original_source_country"]].append(r)
        group=Union()
        for r in data:
            h="H6:"+r["true_H6"]
            group.join(h,"SOIL:"+r["soil_native_source_ID"])
            group.join(h,"POINT:"+r["water_native_point_BasinID"])
        comp=Counter(group.find("H6:"+r["true_H6"]) for r in data)
        varying={}
        for country,z in sorted(countries.items()):
            soil_n=len({r["soil_native_source_ID"] for r in z})
            point_n=len({r["water_native_point_BasinID"] for r in z})
            point_v=len({round(float(r["point_model_N_density_proxy_t_km2"]),12) for r in z})
            atmos_v=len({round(float(r["atmos_model_N_mean"]),12) for r in z})
            varying[country]={
                "H7":len(z),"independent_soil_native_basins":soil_n,
                "independent_point_native_basins":point_n,
                "distinct_pointN_model_density":point_v,
                "distinct_atmos_model_mean":atmos_v,
                "source_variation_country_eligible":soil_n>=3 and point_n>=3 and point_v>1 and atmos_v>1,
            }
        return {
            "H7":len(data),"unique_true_H6":len(h6),
            "original_source_country_labels":len(countries),
            "country_labels":sorted(countries),
            "unique_soil_native_HydroID":len(s),
            "unique_water_point_native_BasinID":len(p),
            "independent_H6_point_soil_connected_components":len(comp),
            "largest_component_H7_n":max(comp.values()) if comp else 0,
            "countries_with_3_native_basins_and_gradient":sum(v["source_variation_country_eligible"] for v in varying.values()),
            "source_country_detail":varying,
        }
    bstrict=blocks(full)
    optimistic=blocks(strict)
    return {
        "original_226_source_H7":len(rows),
        "JRC_4primary_layer_union90_n":sum(r["union_90_all4_core"]=="True" for r in rows),
        "JRC_4primary_layer_single_native90_n":len(strict),
        "JRC_single_native90_with_unverified_material_ring_repair":len(strict)-len(full),
        "strict_geometry_qualified":bstrict,
        "optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED":optimistic,
        "frozen_min_distinct_native_soil_basins":100,
        "independent_native_source_basins_FAIL_even_if_all_repairs_pass":optimistic["unique_soil_native_HydroID"]<100,
        "original_2x2_NUTS2_source_STOP_remains":True,
        "source_only_not_ecological_outcomes":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--summary",type=Path,required=True)
    p.add_argument("--h7-csv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    orig=json.loads(a.summary.read_text())
    if orig.get("schema")!=EXPECTED_SCHEMA or orig["gate"]["original_full_H7_cohort"]!=226:
        raise ValueError("Original completed 226-H7 source gate not verified")
    if orig["gate"]["no_2018_2022_EEQ_OR_PLANT_EIV_NUTRI_BUTTERFLY_HNS_VALUES_READ"] is not True:
        raise ValueError("Original source-only data not attested")
    with a.h7_csv.open(newline="",encoding="utf8") as f:
        rows=list(csv.DictReader(f))
    summary=summarize(rows)
    if summary["strict_geometry_qualified"]["H7"]!=orig["gate"]["H7_conservative_geometry_qualified"]:
        raise ValueError("Original source-only geometry strict count drift")
    if summary["JRC_4primary_layer_single_native90_n"]!=orig["gate"]["H7_with_all4_single_native_polygon_90pct_and_same_soil_id"]:
        raise ValueError("Original single-source 90 count mismatch")
    if summary["strict_geometry_qualified"]["unique_soil_native_HydroID"]!=orig["gate"]["number_native_soil_HydroID_independent"]:
        raise ValueError("Original source-only native exposure independent count drift")
    result={
        "schema":"gesui_20261119_JRC226_H7_SOURCE_ONLY_REPAIR_ACCEPTED_OPTIMISTIC_BOUND_v1",
        "original_source_gate_json_SHA256":sha(a.summary),
        "original_source_H7_CSV_SHA256":sha(a.h7_csv),
        "assumption":"Temporarily accept EVERY material original Esri ring repair as geometrically valid; most permissive possible source upper bound without altering source-record eligibility or single basin 90% threshold",
        "source_only_diagnostics":summary,
        "no_ecological_outcomes_read":True
    }
    a.out.mkdir(parents=True,exist_ok=True)
    outfile=a.out/"JRC226_GEOMETRY_REPAIR_OPTIMISTIC_SOURCE_BOUND.json"
    outfile.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_JRC_226_SOURCE_OPTIMISTIC_BOUND",json.dumps({
        "strict_H7":summary["strict_geometry_qualified"]["H7"],
        "optimistic_single_native_H7":summary["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["H7"],
        "optimistic_soil_basins":summary["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["unique_soil_native_HydroID"],
        "optimistic_H6":summary["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["unique_true_H6"],
        "optimistic_countries":summary["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["original_source_country_labels"],
        "optimistic_source_gradient_countries":summary["optimistic_IF_ALL_ESRI_RING_REPAIRS_ACCEPTED"]["countries_with_3_native_basins_and_gradient"],
        "independent_native_basin_gate_fails_even_optimistically":summary["independent_native_source_basins_FAIL_even_if_all_repairs_pass"]
    },sort_keys=True),flush=True)

if __name__=="__main__":main()
