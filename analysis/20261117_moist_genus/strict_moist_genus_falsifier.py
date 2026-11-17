#!/usr/bin/env python3
"""Source-sha checked 2026-11-17 MOIST between/within-genus adversarial.

This is a POSTHOC observational analysis of genus composition, not source
water chemistry, actual butterfly feeding, or aquatic-terrestrial causality.
Original sites, GBIF facets, paired plant EIV and strict 162 source HOST scores
are independently replayed before new models are computed.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
HOST_DIR = ROOT / "analysis/20261117_host_moisture"
MODEL_DIR = ROOT / "analysis/20261115_eu454_plant_nutri"
sys.path.insert(0, str(HOST_DIR))
sys.path.insert(0, str(MODEL_DIR))

from strict_FerrerParis_pair_H7_water_land_model import (  # noqa: E402
    reconstructed_original_gbif_community,
    load_original_2014_host_score,
    BASE as ORIGINAL_BASE,
)
from rebuild_H7_hostN_F_and_water_land_provisional import (  # noqa: E402
    genuine_original_data, FACET_SHA, FACET_NAME,
    SCORES as ORIGINAL_ROSTER,
)
from h6_model_if_gate import (  # noqa: E402
    cv_h6, source_cluster_regression, design_full, H6_weights,
)

SCHEMA = "GESUI_STRICT_162_HOSTF_EU447_PLANTF_GENUS_ADVERSARIAL_V1"
F_SHA = "ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
F_URL = "https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
DIMENSIONS = ["LIGHT", "TEMPE", "MOIST", "REACT", "NUTRI", "SALIN"]
RANGES = {"LIGHT": (1,9), "TEMPE": (1,12), "MOIST": (1,12),
          "REACT": (1,9), "NUTRI": (1,9), "SALIN": (0,9)}
N_SOURCE_EIV = 8908
N_ALL_CROSS = 11511
N_PAIR_SPECIES = 3722
N_H7_PERIOD = 447
N_H6 = 89
N_STRICT_BFLY = 162
N_WITHIN_SHUFFLES = 4999
N_GLOBAL_SHUFFLES = 999
SEED = 20261117
BFLY_Y = "Bfly_same_host_F_incidence"
BETWEEN = "Plant_F_between_genus"
WITHIN = "Plant_F_within_genus"
PLANT_F = "Plant_F_incidence"
PLANT_N = "Plant_N_incidence"

def checksum_bytes(b: bytes):
    return hashlib.sha256(b).hexdigest()

def find_one(folder: Path, name: str):
    hits = list(folder.rglob(name))
    if len(hits) != 1:
        raise RuntimeError("Expected exactly one original source file: " + name + " got " + str(len(hits)))
    return hits[0]

def binomial(name):
    ss = name.split()
    return len(ss) == 2 and all(re.fullmatch(r"[a-z][a-z-]+", p) for p in ss)

def normalize(n):
    return re.sub(r"\s+", " ", str(n).strip()).casefold() if pd.notna(n) else ""

def source_pairs(cross_root: Path):
    name = "EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv"
    cross = pd.read_csv(find_one(cross_root, name), dtype={"speciesKey":str})
    gate = json.loads(find_one(cross_root, "EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (len(cross) != N_ALL_CROSS or cross.speciesKey.duplicated().any() or
        gate.get("schema") != "GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA" or
        not gate.get("hard_source_support_pass") or
        gate.get("source_hybrid_score_assigned") != 0 or
        gate.get("official_confirmed_valid_unscorable_nonbinomial_hybrids") != 175 or
        gate.get("source_EIV",{}).get("original_source_SHA256") != F_SHA):
        raise RuntimeError("Original GBIF 11511-taxonomy / EIV source-v2 gate STOP")
    response = requests.get(F_URL,timeout=(15,120),
        headers={"User-Agent":"gesui-20261117-plant-F-genus-original-source-crosscheck"})
    response.raise_for_status()
    if checksum_bytes(response.content) != F_SHA:
        raise RuntimeError("Official original EIV 8908 CSV SHA mismatch")
    raw = pd.read_csv(io.BytesIO(response.content),encoding="utf-8-sig")
    if len(raw) != N_SOURCE_EIV or list(raw.columns) != ["Unnamed: 0"] + DIMENSIONS:
        raise RuntimeError("Published original EIV six traits schema mismatch")
    raw["name"] = raw["Unnamed: 0"].map(normalize)
    raw = raw.loc[raw.name.map(binomial)].copy()
    for c in DIMENSIONS:
        raw[c] = pd.to_numeric(raw[c],errors="coerce")
        low,high = RANGES[c]
        raw.loc[~raw[c].between(low,high),c] = np.nan
    ambiguous = set()
    for c in DIMENSIONS:
        z = raw.groupby("name")[c].agg(["min","max"])
        ambiguous |= set(z.index[(z["max"]-z["min"]) > 1e-10])
    raw = raw.loc[~raw.name.isin(ambiguous)].drop_duplicates("name")
    original_traits = raw.set_index("name")[["NUTRI","MOIST"]].to_dict("index")
    paired = []
    numeric_n = 0
    for row in cross.itertuples(index=False):
        if pd.isna(row.EIV_NUTRI_exact):
            continue
        numeric_n += 1
        name = normalize(row.gbif_canonical_name)
        traits = original_traits.get(name)
        if not binomial(name) or traits is None or pd.isna(traits["NUTRI"]):
            raise RuntimeError("Original numeric EIV N species no longer exact-joined: " + row.speciesKey)
        if abs(float(row.EIV_NUTRI_exact)-float(traits["NUTRI"])) > 1e-10:
            raise RuntimeError("Original GBIF accepted N EIV source score drift: "+row.speciesKey)
        if pd.notna(traits["MOIST"]):
            paired.append((str(row.speciesKey), name.split()[0],
                           float(traits["NUTRI"]), float(traits["MOIST"])))
    if numeric_n != 4011 or len(paired) != N_PAIR_SPECIES:
        raise RuntimeError(f"Expected 4011 accepted numeric N, 3722 same-species N+F; got {numeric_n}, {len(paired)}")
    paired.sort(key=lambda z:z[0])
    skey = [r[0] for r in paired]
    name = [r[1] for r in paired]
    N = np.array([r[2] for r in paired],dtype=float)
    F = np.array([r[3] for r in paired],dtype=float)
    unique_genus = sorted(set(name))
    groups = [np.array([j for j,g in enumerate(name) if g==x],dtype=int)
              for x in unique_genus]
    gmean = np.empty(len(paired))
    for group in groups:
        gmean[group] = float(F[group].mean())
    return skey,N,F,gmean,groups,{
        "original_11511_taxa": len(cross),
        "numeric_N_source_species": numeric_n,
        "same_N_F_numeric_original_species":len(paired),
        "accepted_genera":len(groups),
        "multi_species_genera":sum(len(g)>=2 for g in groups),
        "numeric_F_variable_genera":sum(len(g)>=2 and np.ptp(F[g])>0 for g in groups),
        "within_genus_exchangeable_species":sum(len(g) for g in groups if len(g)>=2),
        "EIV_source_SHA256":F_SHA,
        "names_and_scores_do_not_get_published":True
    }

def reconstruct_strict(a):
    merged,_,originalQA = genuine_original_data(a)
    old = pd.read_csv(ORIGINAL_ROSTER, dtype={"butterflyKey":str},
                      encoding="utf-8-sig").reset_index(drop=True)
    scores,sourcegate = load_original_2014_host_score(a.samehosts)
    true_host = {
        str(z.butterflyKey):{
          "N":float(z.source_same_exact_host_NUTRI_mean),
          "F":float(z.source_same_exact_host_MOIST_mean),
          "HC":float(z.source_2014_host_taxa_paired_trait_coverage)}
        for z in scores.itertuples(index=False)
        if str(z.eligible_strict_same_HOST_N_F_and_accepted_butterfly).lower()=="true"
    }
    if len(true_host) != N_STRICT_BFLY:
        raise RuntimeError("Original same-individual-host plant strict butterfly source STOP")
    all447 = merged.copy()
    all447,strictQA = reconstructed_original_gbif_community(a.facets,all447,old,true_host)
    quality = (
        (all447.source_matched_host_butterfly_scored_species>=2) &
        (all447.source_matched_host_butterfly_scored_records>=10) &
        np.isfinite(all447[[BFLY_Y,"Bfly_same_host_N_incidence",
          "Bfly_same_host_F_record",PLANT_N,PLANT_F,
          "EQR_original_site_equal_basin_period_mean"]]).all(axis=1)
    )
    use = all447.loc[quality].copy().reset_index(drop=True)
    if len(use) != N_H7_PERIOD or use.h6.nunique() != N_H6:
        raise RuntimeError("Original 447 H7period and 89 H6 strict source support drift")
    return use,originalQA,strictQA,sourcegate

def reconstruct_plant_source(use,folder,skey,N,F):
    p = find_one(folder, FACET_NAME)
    if checksum_bytes(p.read_bytes()) != FACET_SHA:
        raise RuntimeError("Original source GBIF 2808 facets SHA changed")
    records = json.loads(p.read_text())
    if len(records)!=2808:
        raise RuntimeError("Original GBIF facets not 2808")
    flora = {(str(r["HYBAS_ID_L07"]),r["period"]):r
             for r in records if r["group"]=="Plantae"}
    if len(flora)!=468:
        raise RuntimeError("Original 468 H7×period plant facets missing")
    lookup = {k:j for j,k in enumerate(skey)}
    A = np.zeros((len(use),len(skey)),dtype=float)
    paired_n = []
    for i,row in enumerate(use.itertuples(index=False)):
        data = flora.get((str(row.h),str(row.period)))
        if data is None or data["source_status"]!="ok":
            raise RuntimeError("Original plant facet missing/failed status")
        raw = {str(k):int(v) for k,v in data["source_speciesKey_counts"].items()
               if int(v)>0}
        inds = [lookup[k] for k in raw if k in lookup]
        paired_n.append(len(inds))
        if len(inds)<10 or len(inds)/max(len(raw),1)<.4:
            raise RuntimeError("Plant original 447 pre-scored pair eligibility violated")
        A[i,inds] = 1.0/len(inds)
    if np.max(abs(A.sum(axis=1)-1)) > 1e-10:
        raise RuntimeError("Plant source incidence matrix not normalized")
    f_error = float(np.max(abs(A@F-use[PLANT_F].to_numpy(float))))
    n_error = float(np.max(abs(A@N-use[PLANT_N].to_numpy(float))))
    nsource = np.expm1(use["plant_pair_log_scored_species"].to_numpy(float))
    count_error = float(np.max(abs(nsource-np.array(paired_n))))
    if max(f_error,n_error,count_error)>1e-8:
        raise RuntimeError("Source original same plant N/F means did not match: "+str((f_error,n_error,count_error)))
    return A,{
        "source_original_2808_facet_SHA256":FACET_SHA,
        "site_periods_plant_score_reproduced":len(use),
        "plant_F_max_absolute_replay_error":f_error,
        "plant_N_max_absolute_replay_error":n_error,
        "original_site_pair_species_count_replay_error":count_error
    }

def FWL_setup(use,A,F,base):
    X,labels = design_full(use,base)
    # WLS uses H6-balanced weighting; unchanged original method/country/period.
    root = np.sqrt(H6_weights(use))
    Q,_ = np.linalg.qr(root[:,None]*X,mode="reduced")
    yr = root*use[BFLY_Y].to_numpy(float)
    yres = yr-Q@(Q.T@yr)
    W_A = root[:,None]*A
    obs = W_A@F
    o = obs-Q@(Q.T@obs)
    yr2 = float(yres@yres)
    if yr2<1e-10 or float(o@o)<1e-10:
        raise RuntimeError("Original weighted source full partial variation STOP")
    r = float((yres@o)/np.sqrt(yr2*(o@o)))
    return W_A,Q,yres,yr2,r,{
        "observed_plant_F_conditioned_weighted_partial_r":r,
        "source_H6_equal_weights":True,
        "original_fixed_numeric_and_country_method_period_adjustments":base,
        "original_control_design_rank":X.shape[1],
    }

def nulls(WA,Q,yres,y2,F,groups,ndraw,rng):
    out = np.empty(ndraw,float)
    batch = 28
    for start in range(0,ndraw,batch):
        count = min(batch,ndraw-start)
        traits = np.empty((len(F),count),dtype=float)
        for j in range(count):
            f = F.copy()
            if groups is None:
                f = rng.permutation(F)
            else:
                for g in groups:
                    if len(g)>1:
                        f[g] = rng.permutation(F[g])
            traits[:,j]=f
        z=WA@traits
        z-=Q@(Q.T@z)
        norms=np.sum(z*z,axis=0)
        if not np.all(norms>1e-12):
            raise RuntimeError("Null genus permutation variance degeneracy")
        out[start:start+count]=(yres@z)/np.sqrt(y2*norms)
    return out

def summarize_null(vals,obs):
    n=len(vals)
    med=float(np.median(vals))
    mean=float(np.mean(vals))
    return {
        "n_trait_label_shuffles":n,
        "null_mean":mean,
        "null_q025_q50_q975":np.quantile(vals,[.025,.5,.975]).tolist(),
        "observed_partial_r":obs,
        "upper_tail_p_MC":float((1+np.sum(vals>=obs-1e-12))/(1+n)),
        "centered_two_sided_MC_descriptive_tail":float(
            (1+np.sum(abs(vals-med)>=abs(obs-med)-1e-12))/(1+n)),
        "permutation_not_independent_field_replication":True
    }

def self_test():
    f=np.array([2.,5.,8.,10.,12.],dtype=float)
    g=[np.array([0,1,2]),np.array([3,4])]
    m=np.zeros(len(f))
    for ix in g:m[ix]=f[ix].mean()
    A=np.array([[.5,.5,0,0,0],[0,.25,.25,.25,.25]])
    assert np.allclose(A@f, A@m+A@(f-m))
    assert abs((A@f)[0]-3.5)<1e-12
    return {"status":"PASS_SYNTHETIC_ALGEBRA_IDENTITY", "n":len(f)}

def run(a):
    use,oldQA,strictQA,hostgate=reconstruct_strict(a)
    key,N,F,Fgenus,genus,sourceQA=source_pairs(a.crosswalk)
    A,plantQA=reconstruct_plant_source(use,a.facets,key,N,F)
    between=A@Fgenus
    within=A@(F-Fgenus)
    err=float(np.max(abs(between+within-use[PLANT_F].to_numpy(float))))
    if err>1e-10:
        raise RuntimeError("Exact original 447 community plant MOIST genus decomposition fails")
    use[BETWEEN]=between
    use[WITHIN]=within

    baseline=list(ORIGINAL_BASE)+[PLANT_N]
    if len(set(baseline))!=len(baseline):
        raise RuntimeError("Original PR119 baseline duplicated")
    # Fully independently reproduce original focal source model, after
    # reconstructing exact same original strict host and EIV plant scores.
    orig_ols=source_cluster_regression(use,BFLY_Y,baseline+[PLANT_F],PLANT_F)
    orig_cv=cv_h6(use,BFLY_Y,baseline,PLANT_F)
    original_beta=orig_ols["adjusted_beta_per_1SD_focal"]
    original_gain=orig_cv["H6_equal_added_MSE_gain_pct"]
    if abs(original_beta-0.13357)>0.0001 or abs(original_gain-22.0050)>0.005:
        raise RuntimeError(f"PR119 baseline source match FAIL beta={original_beta}, heldout={original_gain}")

    between_cv=cv_h6(use,BFLY_Y,baseline,BETWEEN)
    within_cv=cv_h6(use,BFLY_Y,baseline+[BETWEEN],WITHIN)
    between_cluster=source_cluster_regression(
        use,BFLY_Y,baseline+[BETWEEN,WITHIN],BETWEEN)
    within_cluster=source_cluster_regression(
        use,BFLY_Y,baseline+[BETWEEN,WITHIN],WITHIN)
    root=np.sqrt(H6_weights(use))
    X, _ = design_full(use,baseline)
    Q,_=np.linalg.qr(root[:,None]*X,mode="reduced")
    yres=root*use[BFLY_Y].to_numpy(float)-Q@(Q.T@(root*use[BFLY_Y].to_numpy(float)))
    bres=root*between-Q@(Q.T@(root*between))
    wres=root*within-Q@(Q.T@(root*within))
    reg=np.column_stack([bres,wres])
    bw_coef=np.linalg.lstsq(reg,yres,rcond=None)[0]
    SSE_base=float(yres@yres)
    qb=float((yres@bres)/(bres@bres))
    SSE_between=float(np.dot(yres-qb*bres,yres-qb*bres))
    SSE_both=float(np.dot(yres-reg@bw_coef,yres-reg@bw_coef))
    if SSE_base < SSE_between-1e-9 or SSE_between<SSE_both-1e-9:
        raise RuntimeError("Nested genus within/between source SSE invalid")

    WA,Q,yres,y2,observed,regQA=FWL_setup(use,A,F,baseline)
    rng=np.random.default_rng(SEED)
    within_null=nulls(WA,Q,yres,y2,F,genus,N_WITHIN_SHUFFLES,rng)
    global_null=nulls(WA,Q,yres,y2,F,None,N_GLOBAL_SHUFFLES,rng)
    # Local within-H6 site contrast strictly source gate only.
    nH7_byH6=use.groupby("h6")["h"].nunique()
    multih6=int((nH7_byH6>=2).sum())
    gate40=multih6>=40

    return {
      "schema":SCHEMA,
      "status":"SOURCE_VERIFIED_POSTHOC_GENUS_ADVERSARIAL_FITTED",
      "source_pr119_original_2026_11_17_verification":{
        "strict_original_butterfly_species":N_STRICT_BFLY,
        "strict_original_H7period":int(len(use)),
        "strict_original_H6":int(use.h6.nunique()),
        "strict_same_larval_HOSTs_original_gate":hostgate.get("schema"),
        "source_447_trueH7_oldGBIF_QA":oldQA,
        "source_162_strict_QA":strictQA,
        "source_EIV_pair_species_QA":sourceQA,
        "source_plant_H7period_N_F_replay":plantQA,
        "replay_original_PR119_plant_F_beta_per_SD":original_beta,
        "replay_original_PR119_plant_F_heldout_H6_MSEgain_pct":original_gain,
      },
      "source_original_F_CWM_decomposition_max_abs_error":err,
      "source_plantF_between_and_within_component_SD":{
        "total":float(np.std(A@F,ddof=0)),
        "between":float(np.std(between,ddof=0)),
        "within":float(np.std(within,ddof=0)),
      },
      "source_genus_component_after_full_original_covariates":{
        "between_genus_added_to_original_baseline_H6_holdout":between_cv,
        "within_genus_ADDED_AFTER_between_H6_holdout":within_cv,
        "joint_two_component_cluster_beta_between":between_cluster,
        "joint_two_component_cluster_beta_within":within_cluster,
        "in_sample_H6_equal_weighted_SSEgain_between_pct":
          100*(SSE_base-SSE_between)/SSE_base,
        "in_sample_H6_equal_weighted_SSEgain_within_AFTER_between_pct":
          100*(SSE_between-SSE_both)/SSE_between,
        "no_causal_sewage_or_independent_field_inference":True,
      },
      "original_butterfly_vs_plant_F_partial_correlation":regQA,
      "source_matched_plant_F_within_GENUS_label_null":
        summarize_null(within_null,observed),
      "source_matched_plant_F_global_label_null":
        summarize_null(global_null,observed),
      "posthoc_original_multi_H7_within_H6_source_gate":{
        "n_parent_H6_with_at_least_two_distinct_H7":multih6,
        "prelocked_minimum_trueH6":40,
        "status":"PASS_SOURCE_ONLY_SUPPLEMENT_NOT_FITTED" if gate40
                 else "STOP_WITHIN_H6_LOCAL_SPATIAL_REPLICATION",
        "no_within_H6_FE_fitted":True
      },
      "positive_spatial_pattern_not_actual_larval_use":True,
      "EIV_butterfly_host_and_plant_share_reference":True,
      "original_gbif_occurrences_not_butterfly_counts":True,
      "no_water_to_land_observed_wastewater_effect":True,
      "original_AT67_2022_S2_riparian_smoke_STOP_unchanged":True,
      "posthoc_not_independently_preregistered":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--self-test",action="store_true")
    p.add_argument("--facets",type=Path)
    p.add_argument("--crosswalk",type=Path)
    p.add_argument("--plants",type=Path)
    p.add_argument("--samehosts",type=Path)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    if a.self_test:
        print("GESUI_MOIST_GENUS_SMOKE="+json.dumps(self_test(),sort_keys=True),flush=True)
        return
    if not all((a.facets,a.crosswalk,a.plants,a.samehosts,a.output)):
        p.error("--facets --crosswalk --plants --samehosts --output all required")
    a.output.mkdir(parents=True,exist_ok=True)
    result=run(a)
    out=a.output/"STRICT_162_HOST_MOIST_EU447_BETWEEN_WITHIN_GENUS_ADVERSARIAL.json"
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    focus={
       "status":result["status"],
       "source_strict_H7":result["source_pr119_original_2026_11_17_verification"]["strict_original_H7period"],
       "genus":result["source_pr119_original_2026_11_17_verification"]["source_EIV_pair_species_QA"],
       "original_F_beta":result["source_pr119_original_2026_11_17_verification"]["replay_original_PR119_plant_F_beta_per_SD"],
       "original_F_H6_MSE_gain":result["source_pr119_original_2026_11_17_verification"]["replay_original_PR119_plant_F_heldout_H6_MSEgain_pct"],
       "between_H6_MSE_gain":result["source_genus_component_after_full_original_covariates"]["between_genus_added_to_original_baseline_H6_holdout"]["H6_equal_added_MSE_gain_pct"],
       "within_after_between_H6_MSE_gain":result["source_genus_component_after_full_original_covariates"]["within_genus_ADDED_AFTER_between_H6_holdout"]["H6_equal_added_MSE_gain_pct"],
       "within_after_between_H6_fixederror95":result["source_genus_component_after_full_original_covariates"]["within_genus_ADDED_AFTER_between_H6_holdout"]["fixed_out_of_fold_error_H6_resampling95_pct"],
       "within_genus_shuffle":result["source_matched_plant_F_within_GENUS_label_null"],
       "global_shuffle":result["source_matched_plant_F_global_label_null"],
       "local_H6_gate":result["posthoc_original_multi_H7_within_H6_source_gate"]
    }
    print("GESUI_MOIST_GENUS_STRICT162_ORIGINAL_H7_RESULT="+json.dumps(
       focus,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    main()
