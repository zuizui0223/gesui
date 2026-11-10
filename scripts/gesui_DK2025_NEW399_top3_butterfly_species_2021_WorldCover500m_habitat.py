#!/usr/bin/env python3
"""POSTHOC species-specific 2025 Danish eBMS reporting vs actual 2021 500m habitat.

Original 2025 species×site records are reconstructed IN MEMORY by separately
audited original 399-site/five-family parser; source files strictly frozen.
Only aggregated species×WWTP source/region outputs are published. This DOES NOT
identify species host plants, nutrient enrichment or causal sewage dose.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import gesui_DK2025_NEW399_original_taxa_pooled_3vs5_denominator_decomposition as species

D=Path("results/ebms_denmark_2023_2024_WWTP_source_only")
HAB=D/"DK2025_NEW399_286SITE_WORLD_COVER2021_500M_NO_FAMILY_SOURCE.json"
HAB_MODEL=D/"DK2025_NEW399_WORLD_COVER2021_LYCAENIDAE_3VS5_WITHIN6_REGION_HABITAT_QC.json"
SPEC=D/"DK2025_NEW399_ORIGINAL_65SPECIES_3VS5_FAMILY_DECOMPOSITION.json"
TOP={
    "Pieris_brassicae":"1920506",
    "Maniola_jurtina":"5135513",
    "Polyommatus_icarus":"5140214"
}
WC=("built","crop","grass","tree")
NBOOT=3000
SEED=20261110

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def fit(rows,usage,include_land):
    regs=sorted({r["region"] for r in rows})
    if any(len({r["band"] for r in rows if r["region"]==g})!=2 for g in regs):
        return None
    X=[];y=[]
    for r in rows:
        data=[1.0,float(r["band"]=="near5")]
        data.extend(float(r["region"]==g) for g in regs[1:])
        data.append(math.log(r["days"]))
        if include_land:data.extend(r["land"][k] for k in WC)
        X.append(data)
        y.append(math.log((r["species"][usage]+.5)/r["days"]))
    mat=np.asarray(X,dtype=float)
    if len(X)<mat.shape[1]+3 or np.linalg.matrix_rank(mat)!=mat.shape[1]:
        return None
    b=np.linalg.lstsq(mat,np.asarray(y,dtype=float),rcond=None)[0]
    result={"n_original_public_site_rows":len(rows),
            "source_shared_1degree_areas":len(regs),
            "log_rate_source_near_adjusted":float(b[1]),
            "exponentiated_original_site_equal_near_vs_remote_recording_rate_ratio":float(math.exp(b[1])),
            "crop_fraction_plus10pct_recorded_rate_ratio":float(math.exp(b[-3]*.1)) if include_land else None,
            "grass_fraction_plus10pct_recorded_rate_ratio":float(math.exp(b[-2]*.1)) if include_land else None,
            "built_fraction_plus10pct_recorded_rate_ratio":float(math.exp(b[-4]*.1)) if include_land else None,
            "tree_fraction_plus10pct_recorded_rate_ratio":float(math.exp(b[-1]*.1)) if include_land else None,
            "reported_counts_plus_half_continuity_and_observation_source_year_only":True}
    return result

def pct(a,q):
    values=sorted(a);return values[int((len(values)-1)*q)]

def run(root):
    previous=json.loads((root/SPEC).read_text(encoding="utf-8"))
    source,site_species,site_meta,taxonomy=species.calculate(root,return_site_level_for_internal_validation=True)
    if source["original_near_vs_remote_three_family_pooled_observational_OR"]!=previous["original_near_vs_remote_three_family_pooled_observational_OR"]:
        raise ValueError("Original species taxonomy/source cohort changed")
    if source["original_near_vs_remote_five_family_pooled_observational_OR"]!=previous["original_near_vs_remote_five_family_pooled_observational_OR"]:
        raise ValueError("Original family source count changed")
    if len(site_meta)!=399 or len(taxonomy)!=63:
        raise ValueError("Original source-locked 399 sites and 63 canonical species changed")
    land=json.loads((root/HAB).read_text())
    model=json.loads((root/HAB_MODEL).read_text())
    if not land["no_2025_butterfly_families_or_numbers_in_this_source"]:
        raise ValueError("Habitat source selected using butterfly outcome")
    if not land["fixed_before_habitat_adjustment_quality_gate"]["habitat_adjusted_2025_butterfly_analysis_may_proceed"]:
        raise ValueError("Original WorldCover source quality failed")
    if model["source_original_SHA256"]["worldcover2021"]!=sha(root/HAB):
        raise ValueError("Original six-region family habitat model source hash differs")
    ww={x["source_original_public_site_SHA256"]:x for x in land["original_site_geo_worldcover_source_counts"]}
    if len(ww)!=286 or any(x["source_QC"]!="COMPLETE" for x in ww.values()):
        raise ValueError("WorldCover pixel source is not complete")
    rows=[]
    for sid,w in ww.items():
        s=site_meta[sid]
        if not s["pre_2025_family_any_occurrence_linked_site_source_gate_pass"]:
            raise ValueError("WorldCover extracted an originally unsupported site")
        if w["source_band"]!=s["source_band"] or w["original_1degree_region"]!=s["original_1degree_region"]:
            raise ValueError("Original source status and area changed")
        if int(s["original_linked_occurrence_dates"])!=w["original_2025_any_published_Occurrence_linked_dates"]:
            raise ValueError("Source reporting days changed")
        rows.append({
            "sid":sid,
            "region":s["original_1degree_region"],
            "band":s["source_band"],
            "days":int(s["original_linked_occurrence_dates"]),
            "land":w["fractions"],
            "species":site_species[sid]
        })
    all_status=Counter(r["band"] for r in rows)
    if all_status!=Counter({"near5":70,"remote10":92,"intermediate":124}):
        raise ValueError("Original 2021 true landcover and 2025 source group cohort changed")
    shared=sorted(set(r["region"] for r in rows if r["band"]=="near5")&
                  set(r["region"] for r in rows if r["band"]=="remote10"))
    sub=[r for r in rows if r["region"] in shared and r["band"] in ("near5","remote10")]
    if len(shared)!=6 or len(sub)!=139 or Counter(r["band"] for r in sub)!=Counter({"near5":61,"remote10":78}):
        raise ValueError("Original six comparable country×1° area cohort changed")
    if shared!=model["within_1degree_comparable_geographic_regions"]:
        raise ValueError("Original 2021 habitat model and source geography mismatch")
    top_set=set(TOP.values())
    if not top_set.issubset(set(taxonomy)):raise ValueError("Selected GBIF source type absent from frozen taxonomy")
    results={}
    for species_name,usage in TOP.items():
        info=taxonomy[usage]
        v={
            "original_GBIF_canonical_usageKey":usage,
            "species_scientificName":info["canonical_species"],
            "original_family":info["family"],
            "source_near_and_remote_2025_original_all286_sites":{},
            "source_region_reporting_corrected_models":{},
            "source_6region_Landcover_adjusted_leave_one_original_region_out":{}
        }
        for band in ("near5","remote10","intermediate"):
            group=[r for r in rows if r["band"]==band]
            x=sum(r["species"][usage] for r in group)
            d=sum(r["days"] for r in group)
            v["source_near_and_remote_2025_original_all286_sites"][band]={
                "original_individualCount":x,"public_monitor_site_rows_with_species":sum(r["species"][usage]>0 for r in group),
                "original_source_linked_survey_dates":d,
                "original_species_recorded_per_published_date":x/d}
            prev_by_taxon=next(x for x in previous["all_accepted_five_family_species_2025_individual_count_by_source_band"]
                              if x["gbif_accepted_usageKey"]==usage)
            if x!=prev_by_taxon["per_group"][band]["original_2025_eBMS_individualCount"]:
                raise ValueError("Original independent site individualCount sums do not match previous species result")
        for key,include_land in (("six_region_survey_days",False),("six_region_survey_days_and_actual_2021_landcover",True)):
            fit0=fit(sub,usage,include_land)
            if fit0 is None:raise ValueError("Region+landcover source model not estimable")
            v["source_region_reporting_corrected_models"][key]=fit0
        for region in shared:
            reduced=[r for r in sub if r["region"]!=region]
            estimate=fit(reduced,usage,True)
            v["source_6region_Landcover_adjusted_leave_one_original_region_out"][region]=(
                estimate["exponentiated_original_site_equal_near_vs_remote_recording_rate_ratio"] if estimate else None)
        results[species_name]=v

    # Whole ORIGINAL 1° areas are resampled; plot no pseudo-independent visits.
    byreg={z:[r for r in sub if r["region"]==z] for z in shared}
    state=SEED
    boot={name:{group:[] for group in ("six_region_survey_days","six_region_survey_days_and_actual_2021_landcover")} for name in TOP}
    for _ in range(NBOOT):
        sampled=[]
        for __ in shared:
            state=(state*1664525+1013904223)&0xFFFFFFFF
            r=shared[int((state/2**32)*len(shared))]
            sampled.extend(byreg[r])
        for name,usage in TOP.items():
            for spec,landflag in (("six_region_survey_days",False),("six_region_survey_days_and_actual_2021_landcover",True)):
                x=fit(sampled,usage,landflag)
                if x:boot[name][spec].append(x["log_rate_source_near_adjusted"])
    for name in TOP:
        for key,record in results[name]["source_region_reporting_corrected_models"].items():
            vals=boot[name][key]
            if len(vals)<NBOOT*.9:
                raise ValueError("Too many original spatial area resamples rank-deficient")
            record["n_valid_original_1degree_block_resamples"]=len(vals)
            record["descriptive_95exp_ratio_from_six_geo_regions"]=[
                math.exp(pct(vals,.025)),math.exp(pct(vals,.975))]
            record["positive_near_logratio_fraction"]=sum(v>0 for v in vals)/len(vals)
            record["only_descriptive_bootstrap_not_WWTP_causal_pvalue"]=True

    return {
        "schema":"gesui_DK2025_new399_top3_species_actual_2021_WorldCover_500m_source_geo_site_rate_POSTHOC_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_SHA256":{
            "original_species_decomposition_63canonical":sha(root/SPEC),
            "source_2021_WorldCover286_pixel_fracs":sha(root/HAB),
            "original_family_2025_SixRegion_model":sha(root/HAB_MODEL)
        },
        "posthoc_selected_3_top_species_after_family_and_species_reporting_composition_patterns_seen":True,
        "original_2025_all399_accepted_species_and_all_family_source_site_counts_reverified":True,
        "original_source_2021_WorldCover286_site_pixel_quality_and_names_reverified":True,
        "original_2025_286_linked_site_counts_by_source":dict(all_status),
        "original_six_shared_1degree_region_ids":shared,
        "original_six_region_original_2025_site_n":len(sub),
        "original_six_region_near_remote_site_counts":dict(Counter(r["band"] for r in sub)),
        "source_species_study_results":results,
        "all_3_species_3000_whole_1degree_region_bootstraps_descriptive":True,
        "habitat_fractions_original_2021_larval_foodplants_NOT_measured":True,
        "source_EEA_2022_NP_annual_report_not_2025_local_stream_exposure":True,
        "none_of_these_species_coefficient_estimates_are_randomized_WWTP_treatments":True,
        "no_individual_observation_rows_or_private_site_species_coordinates_exported":True,
    }

def main():
    q=argparse.ArgumentParser()
    q.add_argument("--out",required=True,type=Path)
    a=q.parse_args()
    root=Path(__file__).resolve().parents[1]
    output=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "original_six_region_sites":output["original_six_region_original_2025_site_n"],
        "species":{name:{
            "near":v["source_near_and_remote_2025_original_all286_sites"]["near5"],
            "remote":v["source_near_and_remote_2025_original_all286_sites"]["remote10"],
            "models":v["source_region_reporting_corrected_models"]
        } for name,v in output["source_species_study_results"].items()}
    },indent=2,ensure_ascii=False))

if __name__=="__main__":
    main()
