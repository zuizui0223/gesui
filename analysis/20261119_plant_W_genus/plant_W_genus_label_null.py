#!/usr/bin/env python3
"""Source-verified EU451 W -> plant NUTRI: between/within accepted genus tests.

Posthoc follow-up to original 2026-11-15 observational W -> plant CWM fit.
Does not test butterfly HNS, nutrient delivery, field biomass, or WWTP causality.
Retains exactly the original 451 H7-periods, 90 genuine H6, 4011 strictly
accepted scored vascular plant species and original GBIF occurrence facets.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis" / "20261115_eu454_plant_nutri"))
import h6_model_if_gate as m

FACETS_SHA = "251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
EIV_SHA = "ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
MODEL_W_BETA_INC_EXPECTED = 0.2152
MODEL_W_BETA_RECORD_EXPECTED = 0.2954
N_H7PERIOD, N_H7, N_H6, N_SCORED, N_GENUS = 451, 232, 90, 4011, 981
BASE = m.P_GEOG

def one(folder: Path, name: str) -> Path:
    files = list(folder.rglob(name))
    if len(files) != 1:
        raise RuntimeError("Source file missing/duplicate: %s: %r" % (name, files))
    return files[0]

def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load_immutable_source(args):
    gate = json.loads(one(args.eiv, "EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (gate.get("schema") != "GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
        or gate.get("source_facets_SHA256") != FACETS_SHA
        or gate.get("source_H7_panel_SHA256") != m.ORIGINAL_PANEL_SHA
        or gate.get("source_EIV", {}).get("original_source_SHA256") != EIV_SHA
        or not gate.get("hard_source_support_pass")
        or gate.get("official_confirmed_valid_unscorable_nonbinomial_hybrids") != 175
        or gate.get("source_hybrid_score_assigned") != 0
        or gate.get("eligible_H7_period") != N_H7PERIOD
        or gate.get("eligible_true_H6") != N_H6):
        raise RuntimeError("Official original source-v2 score/hybrid quality gate changed")
    panel_path = ROOT / "results" / "20261115_sync" / "hns_macro" / "H7_454_HNS_panel.csv"
    if sha(panel_path) != m.ORIGINAL_PANEL_SHA:
        raise RuntimeError("Original 454 H7-period analysis panel SHA failed")
    original = pd.read_csv(panel_path, dtype={"h": str, "h6": str})
    source_eiv = one(args.eiv, "EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv")
    cross = pd.read_csv(source_eiv, dtype={"speciesKey": str})
    if len(cross) != 11511 or cross.speciesKey.duplicated().any():
        raise RuntimeError("Original GBIF 11511 species crosswalk drift")
    scored = cross[cross.EIV_NUTRI_exact.notna()].sort_values("speciesKey", kind="stable")
    if len(scored) != N_SCORED or not (scored.trait_status == "EXACT_MATCH").all():
        raise RuntimeError("Original strictly scored 4011 accepted vascular species lost")
    vals = scored.EIV_NUTRI_exact.to_numpy(float)
    if not np.isfinite(vals).all() or (vals < 1).any() or (vals > 9).any():
        raise RuntimeError("Original NUTRI trait score range changed")
    group = []
    for name in scored.gbif_canonical_name:
        parts = str(name).split()
        if len(parts) != 2 or any(not re.fullmatch(r"[A-Za-z][A-Za-z-]+", x) for x in parts):
            raise RuntimeError("True accepted exact-binomial genus missing")
        group.append(parts[0].casefold())
    genus = np.asarray(group, dtype=str)
    if len(set(genus)) != N_GENUS:
        raise RuntimeError("Original 981 scored genera changed")
    species_to_index = {k: i for i, k in enumerate(scored.speciesKey)}
    eligibility = pd.read_csv(one(args.eiv, "EU454_NUTRI_SOURCE_ELIGIBILITY_PER_H7.csv"),
                              dtype={"h": str, "h6": str})
    means = pd.read_csv(one(args.eiv, "EU454_H7PERIOD_NUTRI_TRAIT_MEANS_SOURCEONLY.csv"),
                        dtype={"h": str})
    if len(eligibility) != 454 or eligibility.duplicated(["h","period"]).any():
        raise RuntimeError("Source 454 H7-period eligibility changed")
    valid = eligibility[eligibility.full_source_trait_eligible.astype(str) == "True"].copy()
    if len(valid) != N_H7PERIOD:
        raise RuntimeError("Source valid 451 H7-period eligibility changed")
    if (valid.NUTRI_coverage_all_Plantae_speciesKeys.lt(.4).any()
        or valid.exact_NUTRI_scored_vascular_accepted_species.lt(10).any()):
        raise RuntimeError("Original source plant score eligibility threshold changed")
    df = original.merge(
        valid[["h","period","NUTRI_coverage_all_Plantae_speciesKeys",
               "exact_NUTRI_scored_vascular_accepted_species"]],
        on=["h","period"], how="inner", validate="one_to_one",
    ).merge(
        means[["h","period","Plant_N_incidence","Plant_N_record","n_scored_species"]],
        on=["h","period"], how="inner", validate="one_to_one",
    )
    if (len(df),df.h.nunique(),df.h6.nunique()) != (N_H7PERIOD,N_H7,N_H6):
        raise RuntimeError("Selected original 451 H7-period, 232 H7, 90 H6 drift")
    if df.country_holdout.nunique() != gate["eligible_majority_country_groups"]:
        raise RuntimeError("Original eligible country grouping changed")
    df["log1p_eiv_scored_taxa"] = np.log1p(
        df.exact_NUTRI_scored_vascular_accepted_species.to_numpy(float))
    df["country_method"] = (
        df.country.fillna("__UNKNOWN__").astype(str) + "::" +
        df.method.fillna("__UNKNOWN__").astype(str)
    )
    if df[["Plant_N_incidence","Plant_N_record","W"]].isna().any().any():
        raise RuntimeError("Original W/plant response missing from eligible source")
    facet_path = one(args.facets, "GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json")
    if sha(facet_path) != FACETS_SHA:
        raise RuntimeError("Genuine 2808 source GBIF facets SHA mismatch")
    facets = json.loads(facet_path.read_text())
    if len(facets) != 2808 or any(x.get("source_status") != "ok" for x in facets):
        raise RuntimeError("Original GBIF facets incomplete")
    plants = {
        (str(x["HYBAS_ID_L07"]),str(x["period"])): x
        for x in facets if x.get("group") == "Plantae"
    }
    if len(plants) != 468:
        raise RuntimeError("Original plant source facets count changed")
    inc = np.zeros((len(df),N_SCORED),dtype=np.float64)
    rec = np.zeros_like(inc)
    missing = []
    for i,row in enumerate(df.itertuples(index=False)):
        key=(str(row.h),str(row.period))
        if key not in plants:
            raise RuntimeError("Original H7-period Plantae facet missing")
        src=plants[key]["source_speciesKey_counts"]
        scored_items=[(species_to_index[str(k)],int(n)) for k,n in src.items()
                      if str(k) in species_to_index and int(n)>0]
        if len(scored_items) != int(row.n_scored_species):
            raise RuntimeError("Source site scored species count mismatch")
        indices=np.asarray([a for a,b in scored_items],int)
        counts=np.asarray([b for a,b in scored_items],float)
        if len(indices)<10 or not (counts>0).all():
            raise RuntimeError("Source taxon quality/positive counts failed")
        inc[i,indices]=1.0/len(indices)
        rec[i,indices]=counts/counts.sum()
    if not np.allclose(inc.sum(axis=1),1.,rtol=0,atol=1e-12):
        raise RuntimeError("Equal species normalization is not 1")
    if not np.allclose(rec.sum(axis=1),1.,rtol=0,atol=1e-12):
        raise RuntimeError("GBIF record weight normalization is not 1")
    maxerr = {
        "species_equal_NUTRI": float(np.max(np.abs(inc@vals-df.Plant_N_incidence))),
        "record_weight_NUTRI":float(np.max(np.abs(rec@vals-df.Plant_N_record))),
    }
    if max(maxerr.values()) >= 1e-10:
        raise RuntimeError("Original source score mean not reproduced " + repr(maxerr))
    return df, vals, genus, inc, rec, maxerr, gate


def decompose(vals,genus):
    groups=defaultdict(list)
    for i,g in enumerate(genus):
        groups[g].append(i)
    gbar=np.zeros(len(vals))
    multi=[]
    for g,idx in groups.items():
        gbar[idx]=float(np.mean(vals[idx]))
        if len(idx)>1:multi.append(np.asarray(idx,dtype=int))
    if not np.allclose(vals, gbar + (vals-gbar),atol=1e-14):
        raise RuntimeError("Genus trait partition not exact")
    return gbar,multi,groups


def regression_fwl_weights(df,weights_by_taxon):
    X,names=m.design_full(df,BASE+["W"])
    focal=names.index("W")
    x=X[:,focal]
    C=np.delete(X,focal,axis=1)
    w=m.H6_weights(df)
    sw=np.sqrt(w)
    if np.linalg.matrix_rank(C) != C.shape[1]:
        raise RuntimeError("Original source geographical WLS reduced model rank-deficient")
    co=np.linalg.lstsq(sw[:,None]*C,sw*x,rcond=None)[0]
    xres=x-C@co
    den=float(np.dot(w*xres,xres))
    if den<1e-10:
        raise RuntimeError("W entirely aliased with geographic background")
    c=(w*xres)@weights_by_taxon/den
    return c


def random_label_test(coeff,vals,genusmean,multigenus,B,seed):
    rng=np.random.default_rng(seed)
    residual=vals-genusmean
    observed_total=float(coeff@vals)
    observed_within=float(coeff@residual)
    global_null=np.zeros(B,float)
    within_null=np.zeros(B,float)
    for i in range(B):
        global_null[i]=coeff@rng.permutation(vals)
        replacement=residual.copy()
        for k in multigenus:
            replacement[k]=rng.permutation(residual[k])
        within_null[i]=coeff@replacement
    def summary(null,observed):
        center=float(null.mean())
        exceed=int(np.sum(np.abs(null-center)>=np.abs(observed-center)-1e-12))
        return {
            "observed_beta":float(observed),
            "null_center":center,
            "two_sided_centered_tail_p":float((exceed+1)/(B+1)),
            "null_percentile95":list(map(float,np.quantile(null,[.025,.975]))),
            "n_permutations":B,
        }
    return {
        "all_4011_scored_plant_species_labels_global":summary(global_null,observed_total),
        "plant_labels_permuted_only_within_accepted_genus":summary(within_null,observed_within),
        "within_genus_null_is_for_within_component_not_total_exposure":True,
    }


def fit(df,A,values,gbar,label,B,seed,multi):
    total=A@values
    between=A@gbar
    within=A@(values-gbar)
    if np.max(np.abs(total-between-within))>1e-10:
        raise RuntimeError("NUTRI genus decomposition identity failed")
    source_name="Plant_N_incidence" if label=="species_equal" else "Plant_N_record"
    if np.max(np.abs(total-df[source_name].to_numpy(float)))>1e-10:
        raise RuntimeError("Raw source CWM means not reproduced")
    full=df.copy()
    for name,arr in (("total",total),("between_genus",between),("within_genus",within)):
        full["computed_"+name]=arr
    regressions={}
    predictions={}
    for name in ("total","between_genus","within_genus"):
        key="computed_"+name
        regressions[name]=m.source_cluster_regression(full,key,BASE+["W"],"W")
        predictions[name]=m.cv_h6(full,key,BASE,"W")
    beta={k:r["adjusted_beta_per_1SD_focal"] for k,r in regressions.items()}
    if abs(beta["total"]-beta["between_genus"]-beta["within_genus"])>1e-9:
        raise RuntimeError("W regression genus coefficients do not add")
    expected=MODEL_W_BETA_INC_EXPECTED if label=="species_equal" else MODEL_W_BETA_RECORD_EXPECTED
    if abs(beta["total"]-expected)>.001:
        raise RuntimeError("Original published +W source regression not reproduced")
    coeff=regression_fwl_weights(full,A)
    if abs(coeff@values-beta["total"])>1e-9:
        raise RuntimeError("FWL source original total coefficient != WLS")
    if abs(coeff@(values-gbar)-beta["within_genus"])>1e-9:
        raise RuntimeError("FWL source within coefficient != WLS")
    if abs(coeff@gbar-beta["between_genus"])>1e-9:
        raise RuntimeError("FWL source between coefficient != WLS")
    orig_cv=m.cv_h6(df,source_name,BASE,"W")
    if abs(orig_cv["H6_equal_added_MSE_gain_pct"] -
           predictions["total"]["H6_equal_added_MSE_gain_pct"])>1e-9:
        raise RuntimeError("Original H6 fold predictions not reproduced")
    return {
        "community_weight_mode":label,
        "reconstructed_original_CWM_variance":float(np.var(total)),
        "between_genus_CWM_variance":float(np.var(between)),
        "within_genus_CWM_variance":float(np.var(within)),
        "between_within_CWM_correlation":float(np.corrcoef(between,within)[0,1]),
        "W_beta_decomposition_source_scale":beta,
        "fraction_of_W_total_slope_from_between_genus_point_only":
            (float(beta["between_genus"]/beta["total"]) if abs(beta["total"])>1e-9 else None),
        "fraction_is_not_contribution_to_causality_or_variance":True,
        "W_adjusted_associations":regressions,
        "true_H6_out_W_added_prediction":predictions,
        "original_CV_exactly_reproduced":True,
        "genus_taxon_labels_permutation":random_label_test(
            coeff,values,gbar,multi,B,seed),
    }


def main(args):
    df,values,genus,Ainc,Arec,maxerr,gate=load_immutable_source(args)
    gbar,multi,groups=decompose(values,genus)
    output={
        "schema":"GESUI_EU451_W_PLANT_EIV_NUTRI_GENUS_LABEL_NULL_V1",
        "date":"2026-11-19",
        "exploratory_posthoc":True,
        "original_EU_W_not_measured_local_wastewater_or_nutrient_supply":True,
        "not_Phylogenetically_corrected_permutation":True,
        "GBIF_occurrence_frequency_not_true_individual_abundance":True,
        "genus_equal_reference_across_all_scored_species":True,
        "original_source_sha256":{
            "H7_panel":m.ORIGINAL_PANEL_SHA,
            "GBIF_2808_facets":FACETS_SHA,
            "EIV_original":EIV_SHA},
        "source_quality_gate":{k:gate[k] for k in [
            "eligible_H7_period","eligible_true_H6","eligible_majority_country_groups",
            "official_confirmed_valid_unscorable_nonbinomial_hybrids"]},
        "n_scored_species":len(values),
        "n_scored_genera":len(groups),
        "n_species_in_multispecies_genera":sum(len(x) for x in multi),
        "n_multispecies_genera":len(multi),
        "source_CWM_replay_max_abs_error":maxerr,
        "analysis":{
            "species_equal":fit(df,Ainc,values,gbar,"species_equal",args.permutations,20261119,multi),
            "record_weight":fit(df,Arec,values,gbar,"record_weight",args.permutations,20261219,multi),
        }
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    digest={}
    for mode,v in output["analysis"].items():
        digest[mode]={
            "W_beta_decomposition":v["W_beta_decomposition_source_scale"],
            "H6_ci":{k:r["H6_cluster_descriptive95"] for k,r in v["W_adjusted_associations"].items()},
            "H6_out_gain_pct":{k:r["H6_equal_added_MSE_gain_pct"] for k,r in v["true_H6_out_W_added_prediction"].items()},
            "H6_out_gain95":{k:r["fixed_out_of_fold_error_H6_resampling95_pct"] for k,r in v["true_H6_out_W_added_prediction"].items()},
            "genus_label_null":v["genus_taxon_labels_permutation"],
        }
    print("GESUI_EU451_PLANT_W_GENUS_TEST="+json.dumps({
        "n_H7period":len(df),"n_H6":df.h6.nunique(),
        "n_EIV_species":len(values),"n_EIV_genera":len(groups),
        "source_CWM_replay_max_abs_error":maxerr,"results":digest,
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--eiv",required=True,type=Path)
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    p.add_argument("--permutations",type=int,default=4999)
    args=p.parse_args()
    main(args)
