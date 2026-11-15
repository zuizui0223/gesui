#!/usr/bin/env python3
"""Restore original GBIF 454 H7 x period x five family source and attribute HNS.

STRICTLY noncausal, posthoc. All calculated variables are derived from occurrence
records and a static larval-host nutrient habitat preference table, NOT counts
of individual butterflies, local plants, nutrient delivery, or feeding events.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "analysis/20261115_multilayer/scripts"))
sys.path.insert(0, str(ROOT / "analysis/20261115_hns/scripts"))
sys.path.insert(0, str(ROOT / "analysis/20261115_hns_family_decomposition"))
from analyse import verify_source, design, COMMON
from score_core import aggregate_traits
from multilayer import cv_add_value, _matrix

FAMS = ["Hesperiidae", "Lycaenidae", "Nymphalidae", "Papilionidae", "Pieridae"]
SCORES_SHA256 = "00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
SOURCE_ZIP_SHA256 = "4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c"
SOURCE_JSON_SHA256 = "251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
ORIGINAL_MEMBER = "GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json"
MIN_SCORED_RECORDS = 10
MIN_SCORED_SPECIES = 2
MIN_SCORE_COVERAGE = 0.70


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def find_source(folder: Path):
    found = []
    for p in sorted(folder.rglob("*.zip")):
        try:
            with zipfile.ZipFile(p) as z:
                m = next((x for x in z.namelist() if x.endswith(ORIGINAL_MEMBER)), None)
                if m:
                    found.append((p, m, json.loads(z.read(m)), sha(p)))
        except zipfile.BadZipFile:
            continue
    for p in sorted(folder.rglob("*.json")):
        if p.name == ORIGINAL_MEMBER:
            found.append((p, None, json.loads(p.read_text()), sha(p)))
    if len(found) != 1:
        raise RuntimeError("Expected one source JSON in original archived artifact, found " +
                           str([(str(x[0]), x[1]) for x in found]))
    p, member, data, digest = found[0]
    if digest != (SOURCE_ZIP_SHA256 if member else SOURCE_JSON_SHA256):
        raise RuntimeError("Original source-frozen field-facet bytes SHA256 mismatch: " + digest)
    if len(data) != 2808 or any(x.get("source_status") != "ok" for x in data):
        raise RuntimeError("Original 234×2×6 facet counts/statuses invalid")
    return p, member, data, digest


def load_reconstruction(panel: Path, score_path: Path, input_dir: Path):
    df = verify_source(panel)
    if sha(score_path) != SCORES_SHA256:
        raise ValueError("Original 303 host score file SHA mismatch")
    taxa = pd.read_csv(score_path, dtype={"butterflyKey": str}, encoding="utf-8-sig")
    if len(taxa) != 303 or taxa.butterflyKey.duplicated().any():
        raise ValueError("Original 303 butterfly taxonomic roster mismatch")
    p, member, raw, original_sha = find_source(input_dir)
    lookup = {(str(x["HYBAS_ID_L07"]), x["period"], x["group"]): x for x in raw}
    if len(lookup) != 2808:
        raise ValueError("Duplicated family/period/source H7 facet key")
    keys = taxa.butterflyKey.astype(str).to_numpy()
    ix = {key: i for i, key in enumerate(keys)}
    families = taxa.family.astype(str).to_numpy()
    alias_bad = taxa.actual_used_alias_conflict.astype(str).str.lower().isin(
        ["true", "1", "yes"]).to_numpy()
    score = taxa.host_nutrient_mean.to_numpy(dtype=float)
    score[(taxa.taxon_status.astype(str).to_numpy() != "ACCEPTED") | alias_bad] = np.nan
    known = np.isfinite(score)
    if known.sum() != 186:
        raise ValueError("Unexpected number of original numeric accepted HNS species")
    B = np.zeros((len(df), len(taxa)), dtype=np.int64)
    raw_family = np.zeros((len(df), len(FAMS)), dtype=np.int64)
    for j, rec in enumerate(df.itertuples(index=False)):
        seen = set()
        for k, family in enumerate(FAMS):
            facet = lookup[(rec.h, rec.period, family)]
            raw_family[j, k] = int(facet["source_query_occurrence_count"])
            for key, count in facet["source_speciesKey_counts"].items():
                if key not in ix or key in seen:
                    raise ValueError(f"Unmatched or duplicated original butterfly speciesKey: {key}")
                i = ix[key]
                if families[i] != family:
                    raise ValueError("SpeciesKey/family mismatch")
                B[j, i] = int(count)
                seen.add(key)
    total_raw = raw_family.sum(axis=1)
    if (int(B.sum()), int(total_raw.sum())) != (4318269, 4400147):
        raise ValueError("Original 454 cohort counts mismatch")
    q = aggregate_traits(B, score, families,
                         taxa.trait_coverage.to_numpy(dtype=float), total_raw)
    for col, z in [("primary_record_mean", q["record_mean"]),
                   ("primary_incidence_mean", q["incidence_mean"]),
                   ("primary_between_family", q["between_family"]),
                   ("primary_within_family", q["within_family"])]:
        np.testing.assert_allclose(z, df[col].to_numpy(dtype=float), atol=1e-10, rtol=0)
    for i, f in enumerate(q["families"]):
        np.testing.assert_allclose(q["record_family_shares"][:,i],
                                   df["primary_rf_"+f], atol=1e-10, rtol=0)
    return df, taxa, B, score, known, families, raw_family, q, {
        "original_archive_name": p.name,
        "original_archive_member": member,
        "source_archive_sha256": original_sha,
        "score_table_sha256": sha(score_path),
        "original_454_source_sum_speciesKey": int(B.sum()),
        "original_454_all_family_record_count": int(total_raw.sum())
    }


def cluster_fit(df, y):
    """Same exact source design as PR94 without family proportions."""
    import statsmodels.api as sm
    X, labels = design(df, include_family=False)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError("Source full-rank design failed")
    h6 = df.h6.astype(str)
    w = len(df) / (h6.nunique() * h6.map(h6.value_counts()).to_numpy(float))
    model = sm.WLS(np.asarray(y, float), X, weights=w).fit(
        cov_type="cluster", cov_kwds={"groups": h6.to_numpy(), "use_correction":True},
        use_t=False)
    j = labels.index("W")
    b, se = float(model.params[j]), float(model.bse[j])
    return {"beta_W_per_1sd_this_cohort": b,
            "cluster95": [b-1.96*se, b+1.96*se],
            "cluster_se": se,
            "posthoc_nominal_p": float(2*norm.sf(abs(b/se))) if se>0 else None,
            "H6": int(h6.nunique()), "n_H7_period": len(df),
            "W_source_SD": float(df.W.std(ddof=0)),
            "adjusted_for": "PR94: original W, urban/crop, flow/area, effort, host score completeness, country×method, period",
            "not_causal": True}


def country_prediction(df, ycol, numeric=COMMON, categorical=("country_method","period")):
    """Leave-one-source-country out; purge any overlapping H6 from training."""
    groups = df.country_holdout.astype(str).to_numpy()
    h6 = df.h6.astype(str).to_numpy()
    y = df[ycol].to_numpy(float)
    country_sizes = {g:int(sum(groups==g)) for g in sorted(set(groups))}
    predictions = {flag: np.full(len(df),np.nan) for flag in (False,True)}
    for group in country_sizes:
        test = groups==group
        train = (groups!=group) & ~np.isin(h6,h6[test])
        if train.sum()<30 or len(set(h6[train]))<5:
            continue
        for withW in (False,True):
            num = list(numeric)+(["W"] if withW else [])
            a,b = _matrix(df.iloc[np.where(train)[0]],
                          df.iloc[np.where(test)[0]],num,list(categorical))
            ct = pd.Series(h6[train]).value_counts()
            w = np.array([1/ct[g] for g in h6[train]],float)
            w /= w.mean()
            gram = a.T@(w[:,None]*a)
            penalty = np.eye(a.shape[1])*10.0
            penalty[0,0] = 0.0
            coef = np.linalg.solve(gram+penalty,a.T@(w*y[train]))
            predictions[withW][test] = b@coef
    country = []
    for g in country_sizes:
        ids = groups==g
        if country_sizes[g]<10 or len(set(h6[ids]))<3:
            continue
        if not (np.isfinite(predictions[False][ids]).all() and
                np.isfinite(predictions[True][ids]).all()):
            continue
        country.append([g,int(ids.sum()),int(len(set(h6[ids]))),
                        float(np.mean((y[ids]-predictions[False][ids])**2)),
                        float(np.mean((y[ids]-predictions[True][ids])**2))])
    if not country:
        return {"status":"no_eligible_country_evaluation"}
    a=np.array([x[3] for x in country]);b=np.array([x[4] for x in country])
    rng=np.random.default_rng(20261115)
    draws=rng.integers(0,len(country),size=(2000,len(country)))
    bs=100*(a[draws].mean(1)-b[draws].mean(1))/a[draws].mean(1)
    return {"country_equal_H6_purged_MSE_gain_pct":float(100*(a.mean()-b.mean())/a.mean()),
            "fixed_error_country_bootstrap95":np.percentile(bs,[2.5,97.5]).tolist(),
            "n_country_evaluation_groups":len(country),
            "country_groups":[x[0] for x in country],
            "no_country_or_H6_leakage":True,
            "descriptive_fixed_errors_not_refit":True,
            "not_independent_replicate_studies":True}


def build_family(df, B, score, known, families, raw_fam, q):
    total_scored = (B*known).sum(1)
    means=dict(zip(q["families"].tolist(),q["family_reference_means"].tolist()))
    if (total_scored==0).any(): raise ValueError("Original zero HNS-scored denominator")
    for k, fam in enumerate(FAMS):
        mask = (families==fam)&known
        cnt = B[:,mask]
        nrec = cnt.sum(1)
        spp = (cnt>0).sum(1)
        r_cov = np.divide(nrec, raw_fam[:,k],
                          out=np.zeros(len(df)), where=raw_fam[:,k]>0)
        score_sum = cnt@score[mask]
        familymean = np.divide(score_sum,nrec, out=np.full(len(df),np.nan),where=nrec>0)
        speciesmean=np.divide((cnt>0)@score[mask],spp,
                              out=np.full(len(df),np.nan),where=spp>0)
        component = (cnt@(score[mask]-means[fam]))/total_scored
        df["component_"+fam]=component
        df["record_mean_"+fam]=familymean
        df["incidence_mean_"+fam]=speciesmean
        df["count_records_"+fam]=nrec
        df["count_species_"+fam]=spp
        df["source_score_coverage_"+fam]=r_cov
        df["support_"+fam]=((nrec>=MIN_SCORED_RECORDS)&
                           (spp>=MIN_SCORED_SPECIES)&
                           (r_cov>=MIN_SCORE_COVERAGE))
    np.testing.assert_allclose(df[["component_"+f for f in FAMS]].sum(1),
                               df.primary_within_family, atol=1e-10,rtol=0)
    return df,means


def main(args):
    df,taxa,B,score,known,families,raw_family,q,receipt = load_reconstruction(
        args.panel,args.scores,args.inputdir)
    df,means=build_family(df,B,score,known,families,raw_family,q)
    out={"schema":"GESUI_H7_454_FAMILY_SPECIES_HNS_ATTRIBUTION_V1",
         "source_receipts":receipt,
         "n_periods":454,"n_H7":232,"n_true_H6":90,
         "n_species_keys_all":303,"n_valid_fixed_score_keys":186,
         "n_original_untreated_basin_controls":0,
         "family_main_score_reference":means,
         "family_support_gate":{"min_scored_records":MIN_SCORED_RECORDS,
                                "min_scored_species":MIN_SCORED_SPECIES,
                                "min_scored_record_fraction":MIN_SCORE_COVERAGE},
         "family_components":{}, "family_within_means":{},
         "posthoc_not_causal":True,
         "family_facet_records_are_GBIF_occurrence_not_abundance":True,
         "host_N_is_Ellenberg_habitat_affinity_not_leaf_N":True}
    for fam in FAMS:
        col="component_"+fam
        fit=cluster_fit(df,df[col].to_numpy(float))
        cv=cv_add_value(df,col,COMMON,["W"])
        out["family_components"][fam]={"W_association":fit,"H6_holdout_W_added_value":cv,
          "n_orig_scored_species":int(sum((families==fam)&known)),
          "family_mean_HNS_reference":float(means[fam]),
          "source_records_all_periods":int(raw_family[:,FAMS.index(fam)].sum()),
          "source_rows_with_scored_HNS":int(sum(df["count_records_"+fam]>0))}
        supported=df[df["support_"+fam]].reset_index(drop=True)
        summary={"n_supported_H7period":len(supported),
                 "H7":int(supported.h.nunique()),
                 "H6":int(supported.h6.nunique()),
                 "n_rows_with_any_scored_record":int(sum(df["count_records_"+fam]>0)),
                 "median_score_coverage_supported":float(
                     supported["source_score_coverage_"+fam].median()) if len(supported) else None}
        if len(supported)>=100 and supported.h6.nunique()>=25:
            for name in ["record_mean_"+fam,"incidence_mean_"+fam]:
                fitmean=cluster_fit(supported,supported[name].to_numpy(float))
                cvmean=cv_add_value(supported,name,COMMON,["W"])
                ctry=country_prediction(supported,name)
                summary[name]={"W_association":fitmean,
                               "H6_holdout_W_added_value":cvmean,
                               "country_h6_purged_W_added_value":ctry,
                               "score_coverage_selection_posthoc":True}
        else:
            summary["insufficient_support_for_multiresponse"]=True
        out["family_within_means"][fam]=summary
    # Additive decomposition of original PR94 within-family W-slope using same X/w.
    full=cluster_fit(df,df.primary_within_family.to_numpy(float))
    source_beta=float(full["beta_W_per_1sd_this_cohort"])
    total=sum(out["family_components"][f]["W_association"]["beta_W_per_1sd_this_cohort"] for f in FAMS)
    if abs(total-source_beta)>1e-9 or abs(source_beta-0.12708367495390768)>1e-6:
        raise ValueError("Family coefficients do not reproduce original full-source within-HNS beta")
    out["within_total_original_W_association"]=full
    out["additivity_full_source_W_beta_maxabs_error"]=abs(total-source_beta)
    # Species-specific arithmetic contributions to observed W slope, not independently tested trends.
    X,names=design(df,False)
    k=names.index("W")
    h6=df.h6.astype(str)
    w=len(df)/(h6.nunique()*h6.map(h6.value_counts()).to_numpy(float))
    X0=np.delete(X,k,axis=1)
    W=X[:,k]
    rw=W-X0@np.linalg.lstsq(X0*np.sqrt(w[:,None]),
                             W*np.sqrt(w),rcond=None)[0]
    denom=np.sum(w*rw**2)
    scores_centered=np.array([(score[i]-means[families[i]]) if known[i] else 0.
                              for i in range(len(score))])
    p=B/(B[:,known].sum(1)[:,None])
    record_weight=(w*rw)@p/denom
    contributions=record_weight*scores_centered
    if abs(contributions.sum()-source_beta)>1e-9:
        raise ValueError("Species contribution sum not equal PR94 within-family beta")
    species=[]
    for i in range(len(score)):
        if not known[i]:continue
        species.append({"speciesKey":str(taxa.butterflyKey.iloc[i]),
                        "accepted_name":str(taxa.accepted_name.iloc[i]),
                        "family":str(families[i]),
                        "host_N_fixed_score":float(score[i]),
                        "source_454_record_count":int(B[:,i].sum()),
                        "within_family_reference_HNS":float(means[families[i]]),
                        "W_species_record_composition_slope_factor":float(record_weight[i]),
                        "W_beta_component_algebraic":float(contributions[i]),
                        "not_single_species_causal_slope":True})
    species.sort(key=lambda z:-abs(z["W_beta_component_algebraic"]))
    for f in FAMS:
        algebra=sum(x["W_beta_component_algebraic"] for x in species if x["family"]==f)
        assert abs(algebra-out["family_components"][f]["W_association"]["beta_W_per_1sd_this_cohort"])<1e-9
    out["top25_species_absolute_W_slope_contributors"]=species[:25]
    out["species_attribution_reconstruction_sum"]=float(contributions.sum())
    out["all_species_attribution_is_posthoc_and_not_significance"]=True
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"FAMILY_HNS_ATTRIBUTION_RESULTS.json").write_text(
        json.dumps(out,ensure_ascii=False,allow_nan=False,indent=2)+"\n")
    pd.DataFrame(species).to_csv(args.output/"SPECIES_W_SLOPE_ATTRIBUTION.csv",index=False)
    tab=[]
    for f in FAMS:
        x=out["family_components"][f]
        c=x["W_association"]
        a=x["H6_holdout_W_added_value"]
        z=out["family_within_means"][f]
        tab.append({
            "family":f,"source_fixed_numeric_species":x["n_orig_scored_species"],
            "source_records_raw":x["source_records_all_periods"],
            "component_beta":c["beta_W_per_1sd_this_cohort"],
            "component_beta_low":c["cluster95"][0],
            "component_beta_high":c["cluster95"][1],
            "component_h6_cv_MSE_gain_percent":a["H6_equal_MSE_gain_percent"],
            "component_h6_cv95_low":a["H6_heldout_fixed_errors_bootstrap95_percent"][0],
            "component_h6_cv95_high":a["H6_heldout_fixed_errors_bootstrap95_percent"][1],
            "supported_rows_for_within_family_mean":z["n_supported_H7period"],
            "supported_H6":z["H6"],
            "within_family_rec_beta":z.get("record_mean_"+f,{}).get("W_association",{}).get(
                "beta_W_per_1sd_this_cohort"),
            "within_family_rec_h6_gain":z.get("record_mean_"+f,{}).get(
                "H6_holdout_W_added_value",{}).get("H6_equal_MSE_gain_percent"),
            "within_family_rec_country_gain":z.get("record_mean_"+f,{}).get(
                "country_h6_purged_W_added_value",{}).get(
                    "country_equal_H6_purged_MSE_gain_pct")
        })
    pd.DataFrame(tab).to_csv(args.output/"FAMILY_HNS_ATTRIBUTION_SUMMARY.csv",index=False)
    print("GESUI_H7_FAMILY_TAXON_RESULT="+json.dumps(
        {"source_receipts":receipt,
         "within_total_beta":source_beta,
         "family_components":{f:out["family_components"][f]["W_association"]["beta_W_per_1sd_this_cohort"] for f in FAMS},
         "family_record_means":{f:{"support":out["family_within_means"][f]["n_supported_H7period"],
                                   "beta":out["family_within_means"][f].get("record_mean_"+f,{}).get("W_association",{}).get("beta_W_per_1sd_this_cohort"),
                                   "gain":out["family_within_means"][f].get("record_mean_"+f,{}).get("H6_holdout_W_added_value",{}).get("H6_equal_MSE_gain_percent")}
                                for f in FAMS},
         "top5_species":species[:5]},allow_nan=False),flush=True)
    print("FULL_H7_FAMILY_TAXON_JSON="+json.dumps(out,ensure_ascii=False,allow_nan=False,separators=(',',':')),flush=True)
    print("FULL_H7_SPECIES_ATTRIBUTION_CSV="+(args.output/"SPECIES_W_SLOPE_ATTRIBUTION.csv").read_text().replace("\n","\\n"),flush=True)
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--panel",type=Path,default=Path("results/20261115_sync/hns_macro/H7_454_HNS_panel.csv"))
    ap.add_argument("--scores",type=Path,default=Path("results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"))
    ap.add_argument("--inputdir",type=Path,default=Path("analysis/20261115_hns_family_taxon/download"))
    ap.add_argument("--output",type=Path,default=Path("results/20261115_hns_family_taxon"))
    main(ap.parse_args())
