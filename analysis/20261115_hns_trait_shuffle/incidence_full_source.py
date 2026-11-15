#!/usr/bin/env python3
"""Independent source-restored incidence vs recorded-dominance HNS label permutations.

Extends the PR95 full-454 within-family slope test with *equal observed species*
incidence using original 2,808 GBIF facets. It does NOT estimate adult/larval
abundance, actual N flux, source treatment, or causal W -> plant -> butterfly.

The 2026-11-15 record-weighted xorshift32 primary remains separately frozen;
this is an explicitly post-result secondary sensitivity using PCG64.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
ORIGINAL_PANEL_SHA="1296d1ded7f1ee185854cad6d1ece1dbb1f94b4e2f48856c8c487a3f3d002fca"
ORIGINAL_TAXON_SHA="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
ORIGINAL_FACET_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
FAMS=["Hesperiidae","Lycaenidae","Nymphalidae","Papilionidae","Pieridae"]
COMMON=["crop","urban","riverflow","naturalout","area","nEEA","nAnimal",
        "nPlant","nButterfly","primary_record_coverage_raw","primary_host_completeness"]

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def quant(v,p):return float(np.quantile(v,p))
def describe(draws,observed):
    ext=int(np.count_nonzero(np.abs(draws)>=abs(observed)-1e-12))
    return {"observed_beta":float(observed),"extreme_permutations":ext,
            "permutation_p_two_sided":float((1+ext)/10000),
            "null_2p5_to_97p5":[quant(draws,.025),quant(draws,.975)],
            "null_median":quant(draws,.5),"null_sd":float(np.std(draws))}

def main(args):
    panel=ROOT/"results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"
    roster=ROOT/"results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"
    if digest(panel)!=ORIGINAL_PANEL_SHA or digest(roster)!=ORIGINAL_TAXON_SHA:
        raise RuntimeError("Genuine original panel or original taxon trait checksum changed")
    df=pd.read_csv(panel,dtype={"h":str,"h6":str})
    if (len(df),df.h.nunique(),df.h6.nunique(),df.country.nunique())!=(454,232,90,15):
        raise RuntimeError("Source geography changed")
    if df.duplicated(["h","period"]).any():raise RuntimeError("Repeated H7-period")
    t=pd.read_csv(roster,dtype={"butterflyKey":str},encoding="utf-8-sig")
    if len(t)!=303 or t.butterflyKey.duplicated().any():raise RuntimeError("Source roster invalid")
    families=t.family.astype(str).to_numpy()
    scores=t.host_nutrient_mean.to_numpy(float)
    bad_alias=t.actual_used_alias_conflict.astype(str).str.lower().isin(["true","1","yes"]).to_numpy()
    scores[(t.taxon_status.astype(str).to_numpy()!="ACCEPTED")|bad_alias]=np.nan
    known=np.isfinite(scores)
    if known.sum()!=186:raise RuntimeError("Missing-score roster differs from the frozen original")
    source=list(args.facets.rglob("GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"))
    if len(source)!=1 or digest(source[0])!=ORIGINAL_FACET_SHA:
        raise RuntimeError("Genuine original 2808 source-facet SHA mismatch")
    facets=json.loads(source[0].read_text())
    if len(facets)!=2808 or any(x.get("source_status")!="ok" for x in facets):
        raise RuntimeError("Genuine source-facet status/denominator changed")
    lookup={(str(x["HYBAS_ID_L07"]),str(x["period"]),x["group"]):x for x in facets}
    if len(lookup)!=2808:raise RuntimeError("Duplicate original facets")
    ix={k:i for i,k in enumerate(t.butterflyKey)}
    B=np.zeros((len(df),len(t)),dtype=np.int64)
    raw_total=0
    for j,row in enumerate(df.itertuples(index=False)):
        for fam in FAMS:
            facet=lookup[(row.h,row.period,fam)]
            raw_total+=int(facet["source_query_occurrence_count"])
            for k,n in facet["source_speciesKey_counts"].items():
                if k not in ix or families[ix[k]]!=fam or B[j,ix[k]]!=0:
                    raise RuntimeError("Original speciesKey family/count reconstruction mismatch")
                B[j,ix[k]]=int(n)
    if int(B.sum())!=4318269 or raw_total!=4400147:
        raise RuntimeError("Original GBIF source counted totals changed")
    g=df.h6.astype(str)
    w=len(df)/(g.nunique()*g.map(g.value_counts()).to_numpy(float))
    df["country_method"]=(df.country.fillna("__UNKNOWN__").astype(str)+"::"+
                          df.method.fillna("__UNKNOWN__").astype(str))
    X=[np.ones(len(df))];names=["intercept"]
    for k in COMMON+["W"]:
        v=df[k].to_numpy(float)
        if not np.isfinite(v).all():raise RuntimeError("Missing fixed design "+k)
        sd=float(np.std(v))
        X.append((v-v.mean())/(sd if sd>1e-12 else 1.))
        names.append(k)
    for k in ["country_method","period"]:
        s=df[k].astype(str)
        for v in sorted(s.unique().tolist())[1:]:
            X.append((s==v).to_numpy(float))
            names.append(k+"="+v)
    X=np.column_stack(X)
    if np.linalg.matrix_rank(X)!=X.shape[1]:
        raise RuntimeError("Original WLS design not full-rank")
    k=names.index("W")
    W=X[:,k]
    X0=np.delete(X,k,axis=1)
    sq=np.sqrt(w)
    resid=W-X0@np.linalg.lstsq(X0*sq[:,None],W*sq,rcond=None)[0]
    denominator=float(np.sum(w*resid**2))
    if not denominator>1e-10:raise RuntimeError("No residual W exposure")
    scored_idx=np.flatnonzero(known)
    rec_den=B[:,known].sum(axis=1)
    inc_den=(B[:,known]>0).sum(axis=1)
    if (rec_den<1).any() or (inc_den<1).any():raise RuntimeError("Zero scored response denominator")
    prefRec=B[:,known]/rec_den[:,None]
    prefInc=(B[:,known]>0)/inc_den[:,None]
    beta_factors_rec=(w*resid)@prefRec/denominator
    beta_factors_inc=(w*resid)@prefInc/denominator
    f= families[known]
    score =scores[known]
    ref=np.array([float(np.mean(score[f==family])) for family in f])
    centered=score-ref
    obs_rec=float(beta_factors_rec@centered)
    obs_inc=float(beta_factors_inc@centered)
    obs_rec_total=float(beta_factors_rec@score)
    obs_inc_total=float(beta_factors_inc@score)
    if abs(obs_rec-.12708367495390768)>1e-9:
        raise RuntimeError("Fails independent source PR95 WLS within-family coefficient")
    if abs(obs_inc_total-.06998)>.0002:
        raise RuntimeError("Fails original PR94 incidence equal source benchmark")
    # Independent record-factor audit from frozen PR95 source CSV.
    source_path=ROOT/"results/20261115_hns_family_taxon/SPECIES_W_SLOPE_ATTRIBUTION.csv"
    with source_path.open(newline="",encoding="utf-8") as handle:
        saved={str(x["speciesKey"]):float(x["W_species_record_composition_slope_factor"])
               for x in csv.DictReader(handle)}
    if len(saved)!=186 or max(abs(saved[str(t.butterflyKey.iloc[i])]-beta_factors_rec[j])
                            for j,i in enumerate(scored_idx))>1e-10:
        raise RuntimeError("Independent projected record-weight factors do not match archived source")
    rng=np.random.default_rng(20261115)
    positions=[np.flatnonzero(f==fam) for fam in FAMS]
    rec_dist=np.empty(9999);inc_dist=np.empty(9999)
    for b in range(9999):
        draw=np.empty_like(centered)
        for ixg in positions:
            draw[ixg]=rng.permutation(centered[ixg])
        rec_dist[b]=beta_factors_rec@draw
        inc_dist[b]=beta_factors_inc@draw
    delta_dist=rec_dist-inc_dist
    output={
       "schema":"GESUI_EU_H7_WITHIN_FAMILY_HNS_INCIDENCE_VS_RECORD_LABEL_NULL_V1",
       "original_panel_SHA256":ORIGINAL_PANEL_SHA,"original_taxon_SHA256":ORIGINAL_TAXON_SHA,
       "original_2808_facets_SHA256":ORIGINAL_FACET_SHA,
       "n_original_H7_periods":454,"n_H6":90,"n_scored_species":186,
       "n_raw_5family_GBIF_records":raw_total,"n_speciesKey_records":int(B.sum()),
       "full_record_HNS_beta_same_fixed_design":obs_rec_total,
       "full_incidence_HNS_beta_same_fixed_design":obs_inc_total,
       "main_record_within_family_beta_recovered":obs_rec,
       "incidence_within_family_beta_new":obs_inc,
       "nine_thousand_ninehundred_ninety_nine_shuffles":9999,
       "rng":"numpy PCG64 default_rng(20261115)",
       "trait_reference_within_family":"unweighted scored species mean fixed within each family",
       "record_within_family_posthoc_parallel":describe(rec_dist,obs_rec),
       "incidence_within_family_posthoc_secondary":describe(inc_dist,obs_inc),
       "record_minus_incidence_within_family_secondary":describe(delta_dist,obs_rec-obs_inc),
       "same_permutations_used_for_both_weightings":True,
       "frozen_original_model_and_454_cohort":True,
       "caveats":[
         "Post-result conditional trait-label exchangeability not independently pre-registered",
         "Equal observed-species incidence is not detection-corrected ecological occupancy",
         "HNS is modelled plant habitat preference, not actual larval plant use or field leaf N",
         "No wastewater exposure delivery, stream-to-land food-web path or causal mediation estimated",
         "Original score uncertainties, recording effort and phylogenetic constraints not randomized"
       ]
    }
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"HNS_WITHINFAMILY_INCIDENCE_VS_RECORD_PERMUTATION.json").write_text(
        json.dumps(output,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_HNS_INCIDENCE_VS_RECORD_PERM="+json.dumps({
       "record_observed":obs_rec,"record_perm":output["record_within_family_posthoc_parallel"],
       "incidence_observed":obs_inc,"incidence_perm":output["incidence_within_family_posthoc_secondary"],
       "delta_observed":obs_rec-obs_inc,"delta_perm":output["record_minus_incidence_within_family_secondary"]
       },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--facets",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    main(ap.parse_args())

