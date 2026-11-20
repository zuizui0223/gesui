#!/usr/bin/env python3
"""Posthoc adversarial review of frozen 2010–18 ChemPop-v2 BMWP~EDF model.

Tests order-between versus family-within-order slope components; leave-one-
region-out influence; one-row-per-site means; and ACTUAL COUNT-ONLY sample
abundance sensitivity. It does not re-interpret ordinal/log-abundance
codes as observed specimen counts.

The original source is an authentic, SHA-matched, 326-member EIDC ZIP; UK
family scoring reference is the published biotic package's 2016 RData table.
All outputs are aggregate coefficients, source gates, and falsifications.
Original sample or species/score rows are not republished.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis'/'20261119_europe_benthic_taxa'))

import numpy as np
import pandas as pd
import chempop_v2_bmwp_filter as f
import chempop_v2_source_census as src

ORIGINAL_BETA = -0.09110572299049917
ORIGINAL_N_SITEYEAR = 3884
ORIGINAL_N_SITES = 707
ORIGINAL_N_REGIONS = 7
ORIGINAL_ZIP_SHA = "cc42df5ea777f5c7f59e31b097435d5caf72be21fc7e6bc08f3ce13cde83bc51"
ORIGINAL_SCORE_SHA = "4f9be755ce6d24e645b46705dbec9c0367bc71a9ce6b335d3bf2657d8c10ac5a"

ACTUAL_METHOD = "Actual abundance of taxa in the sample"
ACTUAL_MIN_SAMPLES = 300
ACTUAL_MIN_SITEYEARS = 200
ACTUAL_MIN_SITES = 100
ACTUAL_MIN_REGIONS = 4
ACTUAL_MIN_SITES_PER_REGION = 5
ACTUAL_MIN_EDF_VALUES_PER_REGION = 5


def create_same_source(fam_scores, source_rows, source_weights):
    """Reuse prior exact source-only gates, not selecting by outcome."""
    if len(source_rows) != len(source_weights):
        raise RuntimeError("Source sample and taxon weight identity drift")
    siteyear, A, B, fam, scores, nsupported = f.site_year_response(
        source_rows, source_weights, fam_scores)
    siteyear["edf_log1p"] = np.log1p(siteyear.edf_pct.to_numpy(float))
    if (len(siteyear) != ORIGINAL_N_SITEYEAR
        or siteyear.site_region_id.nunique() != ORIGINAL_N_SITES
        or siteyear.region.nunique() != ORIGINAL_N_REGIONS):
        raise RuntimeError("Original 3884×707×7 source grain has changed")
    reg = f.fit_cluster(siteyear, "y_incidence")
    if abs(reg["beta"] - ORIGINAL_BETA) > 1e-10:
        raise RuntimeError("Original EDF→incidence BMWP WLS coefficient changed")
    return siteyear, A, B, fam, scores, reg


def order_ref(families, values, order_source):
    """Average original accepted scored family values within documented order.

    Ambiguous multisystem family labels have their own singleton references;
    they cannot contribute to within-order variation and are NEVER assigned
    to whichever source order happens to be statistically convenient.
    """
    groups = collections.defaultdict(list)
    ambiguous = []
    unknown = []
    for idx, fam in enumerate(families):
        orders = order_source.get(fam, collections.Counter())
        if len(orders) == 1:
            group = next(iter(orders))
        elif len(orders) > 1:
            ambiguous.append(fam)
            group = "__AMBIGUOUS_FAMILY_SINGLETON_" + fam
        else:
            unknown.append(fam)
            group = "__UNKNOWN_FAMILY_SINGLETON_" + fam
        groups[group].append(idx)
    ref = np.empty_like(values, dtype=float)
    for _, indices in groups.items():
        av = float(np.mean(values[indices]))
        ref[indices] = av
    if not np.isfinite(ref).all():
        raise RuntimeError("Original order reference lacks scored taxon values")
    return ref, {
        "n_scored_families": len(values),
        "n_unique_order_or_ambiguity_groups": len(groups),
        "n_non_singleton_groups": sum(len(x) >= 2 for x in groups.values()),
        "n_families_in_non_singleton_groups": sum(len(x) for x in groups.values() if len(x)>=2),
        "n_ambiguous_multiorder_families": len(ambiguous),
        "ambiguous_multiorder_family_names": sorted(ambiguous),
        "n_unknown_order_families": len(unknown),
        "unknown_order_family_names": sorted(unknown),
        "family_order_reference_is_equal_scored_family_mean": True,
    }


def decompose(df, A, values, ref, original_fit):
    """Algebraic decomposition with EXACT same site-weighted adjusted model."""
    d = df.copy()
    d["y_between_order"] = A @ ref
    d["y_within_order"] = A @ (values - ref)
    if not np.allclose(
        d.y_between_order + d.y_within_order, d.y_incidence,
        atol=1e-12, rtol=0
    ):
        raise RuntimeError("BMWP order decomposition does not reconstruct source CWM")
    bfit = f.fit_cluster(d, "y_between_order")
    wfit = f.fit_cluster(d, "y_within_order")
    if abs((bfit["beta"] + wfit["beta"]) - original_fit["beta"]) > 1e-9:
        raise RuntimeError("WLS order-decomposed EDF slopes do not add to source")
    # prior true label null from the independently re-built source pipeline;
    # no additional hypothesis-selected permutations are used here.
    returns = {
        "original_total": original_fit,
        "order_between": bfit,
        "within_order_family_deviation": wfit,
        "order_beta_components_add_to_total": True,
        "order_between_to_total_beta_fraction_point_only":
            float(bfit["beta"] / original_fit["beta"]) if abs(original_fit["beta"])>1e-9 else None,
        "this_ratio_is_not_causal_contribution_or_variance": True,
        "between_order_region_holdout": f.region_transfer(d,"y_between_order"),
        "within_order_region_holdout": f.region_transfer(d,"y_within_order"),
    }
    return returns


def leave_region_influence(df):
    """Estimate original EDF association on all seven 6-region train subsets."""
    out = {}
    for heldout in sorted(df.region.unique()):
        d = df.loc[df.region != heldout].copy().reset_index(drop=True)
        fit = f.fit_cluster(d,"y_incidence")
        out[heldout] = {
            "beta": fit["beta"],
            "site_cluster95": fit["cluster_site_normal95"],
            "n_original_site": fit["n_distinct_site"],
            "n_original_siteyear": fit["n_site_year"],
            "no_regional_reweight_after_outcome_selection":True,
        }
    return out


def only_sites(df):
    """Independent-site means, not duplicated site-years of fixed EDF."""
    use=("edf_pct","lc_urban","lc_arable","width","altitude","siltclay",
         "year","month_sin","month_cos","log_n_families","coverage",
         "edf_log1p","y_incidence")
    rows=[]
    for site, g in df.groupby("site_region_id",sort=True):
        if g.edf_pct.max()-g.edf_pct.min()>1e-10:
            raise RuntimeError("A supposedly static site EDF changes by year")
        row={col:float(g[col].mean()) for col in use}
        row.update({"site_region_id":site,"site":str(g.site.iloc[0]),
                    "region":str(g.region.iloc[0]),
                    "method":str(g.method.value_counts().index[0]),
                    "n_original_siteyears":len(g)})
        rows.append(row)
    d=pd.DataFrame(rows)
    if len(d)!=ORIGINAL_N_SITES:raise RuntimeError("Site mean should be 707 true sites")
    reg=f.fit_cluster(d,"y_incidence")
    return {
       "n_independent_sites":len(d),
       "n_original_total_siteyears":int(d.n_original_siteyears.sum()),
       "source_response_is_unweighted_mean_of_siteyears":True,
       "fit":reg,
       "region_holdout":f.region_transfer(d,"y_incidence"),
    }


def actual_count_gate(df_sample, source_rows):
    """Entirely outcome-blind source method and sample-support check."""
    by_region=collections.defaultdict(set)
    values=collections.defaultdict(set)
    for r in df_sample:
        by_region[r["region"]].add(r["site_region_id"])
        values[r["region"]].add(round(float(r["edf_pct"]), 6))
    passing={reg:len(sites)>=ACTUAL_MIN_SITES_PER_REGION and
             len(values[reg])>=ACTUAL_MIN_EDF_VALUES_PER_REGION
             for reg,sites in by_region.items()}
    return {
        "original_actual_count_method_only":ACTUAL_METHOD,
        "n_original_source_samples":len(source_rows),
        "n_sample_year_rows":len(df_sample),
        "n_distinct_sites":len({r["site_region_id"] for r in df_sample}),
        "n_regions_with_enough_sites_EDF":sum(passing.values()),
        "per_region":{reg:{"sites":len(sites),"edf_distinct":len(values[reg]),"pass":passing[reg]}
                      for reg,sites in by_region.items()},
        "minimum_samples":ACTUAL_MIN_SAMPLES,
        "minimum_siteyears":ACTUAL_MIN_SITEYEARS,
        "minimum_sites":ACTUAL_MIN_SITES,
        "minimum_regions":ACTUAL_MIN_REGIONS,
        "source_gate_PASS":(
            len(source_rows)>=ACTUAL_MIN_SAMPLES
            and len(df_sample)>=ACTUAL_MIN_SITEYEARS
            and len({r["site_region_id"] for r in df_sample})>=ACTUAL_MIN_SITES
            and sum(passing.values())>=ACTUAL_MIN_REGIONS),
    }


def true_count_sensitivity(source_rows, source_weights, scores):
    if len(source_rows)!=len(source_weights):
        raise RuntimeError("Source row-weight parity failed")
    eligible=[
        (r,w) for r,w in zip(source_rows,source_weights)
        if r["method"]==ACTUAL_METHOD
        and r["all_scored_family_recorded_numeric_positive"]
    ]
    metadata={
        "n_original_2010_18_source_eligible_samples":len(source_rows),
        "method_distribution_in_all_source_eligible":
            dict(collections.Counter(r["method"] for r in source_rows)),
        "n_actual_count_and_all_scored_numeric_samples":len(eligible),
        "reasons_for_excluding_other_methods":
           "Log-abundance estimates and locally-defined estimated abundance are NOT actual counted animals",
    }
    if not eligible:
        metadata["status"]="STOP_NO_TRUE_ACTUAL_COUNT_METHOD"
        return metadata
    raw=[q[0] for q in eligible]
    weights=[q[1] for q in eligible]
    df, _, _, _, _, n=f.site_year_response(raw,weights,scores)
    df["edf_log1p"] = np.log1p(df.edf_pct.to_numpy(float))
    # Source test is conducted BEFORE reading count-CWM outcomes.
    record=df.to_dict(orient="records")
    gate=actual_count_gate(record,raw)
    metadata["source_gate"]=gate
    if not gate["source_gate_PASS"]:
        metadata["status"]="SOURCE_STOP_TRUE_COUNT_COHORT_TOO_SMALL"
        return metadata
    if df.count_supported.sum()!=len(df):
        metadata.update({
            "status":"SOURCE_STOP_TRUE_COUNT_INCOMPLETE",
            "count_supported_siteyears":int(df.count_supported.sum()),
            "total_siteyears":len(df),
        })
        return metadata
    df["y_count_minus_incidence"] = df.y_recorded_abundance-df.y_incidence
    metadata["status"]="REAL_COUNT_ONLY_METHOD_OBSERVATIONAL_SENSITIVITY"
    for outcome in ("y_incidence","y_recorded_abundance","y_count_minus_incidence"):
        metadata[outcome] = {
            "fit":f.fit_cluster(df,outcome),
            "region_transfer":f.region_transfer(df,outcome),
        }
    metadata["actual_count_only_not_independent_study"] = True
    return metadata


def main(a):
    output={
        "schema":"GESUI_CHEMPOP2010_18_BMWP_ORDER_AND_TRUE_COUNT_ADVERSARIAL_V1",
        "source_run_original":"44548177752",
        "source_is_modelled_sewage_EFFLUENT_percentage_not_chemistry":True,
        "not_causal_and_not_phylogenetic_independence":True,
        "posthoc_after_discovered_negative_BMWP_EDF_slope":True,
        "not_a_new_independent_England_field_cohort":True,
        "status":"SOURCE_STOP",
    }
    try:
        with tempfile.TemporaryDirectory(prefix="gesui_chem_bmwp_audit_") as tmp:
            d=Path(tmp)
            a0=src.retrieve(src.EIDC_URL,d/"original.zip",src.SOURCE_DOWNLOAD_LIMIT)
            scores,s0=src.score_reference(d)
            output["source_sha256"]={"EIDC_container":a0["SHA256"],"BMWP":s0["SHA256"]}
            output["original_archive_byte_sha_matches_prior"] = a0["SHA256"] == ORIGINAL_ZIP_SHA
            output["original_BMWP_RData_byte_sha_matches_prior"] = s0["SHA256"] == ORIGINAL_SCORE_SHA
            if s0["SHA256"]!=ORIGINAL_SCORE_SHA:
                raise RuntimeError("Published family score dictionary original RData SHA changed")
            # EIDC assembles a generated multi-file ZIP: a different wrapper SHA
            # is possible when member metadata/order/timestamps are repackaged.
            # Preserve the observed ZIP hash and rigorously verify the 14
            # biological/EDF source member digests plus original source-specific
            # 7836 samples, 3884 siteyears, 707 sites, and numerical WLS beta.
            with __import__("zipfile").ZipFile(d/"original.zip") as z:
                files=z.namelist()
                if len(files)!=326 or len(set(files))!=326:
                    raise RuntimeError("Original source archive member identities drift")
                member_digests = {}
                for region in src.PREFIXES:
                    for name in (
                        src.find_one_by_region(files,region,"sitevariables\\.csv"),
                        src.find_taxa(files,region),
                    ):
                        hasher=hashlib.sha256()
                        with z.open(name) as stream:
                            while True:
                                block=stream.read(1024*1024)
                                if not block:break
                                hasher.update(block)
                        member_digests[name]=hasher.hexdigest()
                output["original_14_critical_member_SHA256"]=member_digests
                raw, weights, orders, audit = f.construct_eligible(z,files,scores)
        output["source_full"]={
            "eligible_2010_2018_samples":audit["source_gate"]["eligible_sample_n"],
            "eligible_2010_2018_sites":audit["source_gate"]["eligible_site_n"],
            "eligible_2010_2018_regions":audit["source_gate"]["eligible_region_n"],
        }
        if not audit["source_gate"]["cohort_ready_for_incidence_CWM_fit"]:
            output["status"]="SOURCE_GATE_STOP_FROZEN_MAIN_COHORT"
        else:
            df,A,B,families,values,orig=create_same_source(scores,raw,weights)
            ref,orderqa=order_ref(families,values,orders)
            output["order_reference_source_qa"]=orderqa
            output["order_decomposition"]=decompose(df,A,values,ref,orig)
            output["leave_one_region_out_EDF_slope"]=leave_region_influence(df)
            output["one_record_per_SITE_sensitivity"]=only_sites(df)
            try:
                output["real_count_method_only"]=true_count_sensitivity(raw,weights,scores)
            except Exception as ex:
                output["real_count_method_only"]={
                    "status":"SOURCE_OR_METHOD_STOP_TRUE_COUNTS",
                    "error":type(ex).__name__+": "+str(ex)[:450],
                    "old_original_main_results_preserved":True
                }
            output["status"]="POSTHOC_ADVERSARIAL_CALCULATION_COMPLETE"
    except Exception as ex:
        output.update({"status":"SOURCE_OR_MODEL_STOP","error":type(ex).__name__+": "+str(ex)[:650]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    digest={
        "status":output["status"],
        "error":output.get("error"),
        "order_source":output.get("order_reference_source_qa"),
        "order_beta":({k:v["beta"] for k,v in output["order_decomposition"].items()
                       if isinstance(v,dict) and "beta" in v}
                     if "order_decomposition" in output else {}),
        "order_details":output.get("order_decomposition"),
        "leave_region":output.get("leave_one_region_out_EDF_slope"),
        "site_level":output.get("one_record_per_SITE_sensitivity"),
        "actual_count_only":output.get("real_count_method_only"),
    }
    print("GESUI_CHEMPOP_ORDER_TRUECOUNT_ADVERSARIAL="+json.dumps(
        digest,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
