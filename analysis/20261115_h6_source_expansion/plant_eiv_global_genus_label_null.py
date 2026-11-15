#!/usr/bin/env python3
"""Posthoc original-source EIV label shuffles for H6-FE plant–butterfly HNS.

All 183 H7/40 genuine H6 and 13 source controls retained. Exact original
Plantae speciesKey facet incidence for each of two GBIF periods is used to
recreate source plant_N_incidence before any shuffling. Two precommitted
source conditional nulls: 4999 globally across all 4011 scored vascular
species and 4999 within accepted-binomial genus. No taxonomic imputations.
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

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261115_h6_source_expansion"))
from test_183_40_within_H6 import (
    clean_source,features,weights,CONTROLS,PREDICTOR,Y_PRIMARY,
    PANEL_SOURCE_SHA,SOURCE_SCHEMA
)

OLD_FACETS_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
NEW_FACETS_SHA="e076f6c8e2d93282be19aea3c75ae8431302318f068b9d75e278e2a21c20aafb"
SCORE_H7_SOURCE_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
SOURCE_QA_SCHEMA="GESUI_EU454_FULLSOURCE_PLANT_NUTRI_TRAIT_GATE_V2_HYBRID_QA"
N_SCORED=4011
N_DRAW=4999
SEED=20261115
PERIODS=["2017-2019","2020-2022"]
ORIGINAL_BETA=0.11807615189599471

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def source_one(folder,name):
    hits=list(folder.rglob(name))
    if len(hits)!=1:
        raise RuntimeError("Exactly one original genuine source required: "+name+" "+repr([str(x) for x in hits]))
    return hits[0]

def load_facets(folder,name,expected,n):
    file=source_one(folder,name)
    if sha(file)!=expected:
        raise RuntimeError("Original source GBIF facets SHA mismatch "+name)
    records=json.loads(file.read_text())
    if len(records)!=n or any(r.get("source_status")!="ok" for r in records):
        raise RuntimeError("Original GBIF source facet missing/invalid status")
    return {(str(r["HYBAS_ID_L07"]),str(r["period"]),str(r["group"])):r
            for r in records}

def source_frame(a):
    p=source_one(a.source183,"GENUINE_TERRESTRIAL_ONLY_183_H7_40_H6_PAIRED_H7_MEANS.csv")
    qa=json.loads(source_one(a.source183,"GENUINE_TERRESTRIAL_183_H7_40_H6_SOURCE_RECONSTRUCTION_QA.json").read_text())
    if qa.get("schema")!=SOURCE_SCHEMA:
        raise RuntimeError("Wrong old/new 183 H7 source build provenance")
    frame=clean_source(p,qa).sort_values(["h6","h"],kind="stable").reset_index(drop=True)
    if frame.h.nunique()!=183 or frame.h6.nunique()!=40:
        raise RuntimeError("Genuine original geographic H7/H6 source drift")
    cross=source_one(a.eiv,"EU454_PLANT_11511_SOURCE_EIV_EXACT_CROSSWALK.csv")
    gate=json.loads(source_one(a.eiv,"EU454_FULL_EIV_SOURCE_GATE.json").read_text())
    if (gate.get("schema")!=SOURCE_QA_SCHEMA or not gate.get("hard_source_support_pass")
        or gate.get("source_hybrid_score_assigned")!=0
        or gate.get("official_confirmed_valid_unscorable_nonbinomial_hybrids")!=175
        or gate.get("source_EIV",{}).get("original_source_SHA256")!=SCORE_H7_SOURCE_SHA):
        raise RuntimeError("Official full 4011 EIV source-v2 checksum/typed hybrid gate failed")
    taxa=pd.read_csv(cross,dtype={"speciesKey":str})
    if len(taxa)!=11511 or taxa.speciesKey.duplicated().any():
        raise RuntimeError("Original source union 11511 GBIF speciesKey IDs changed")
    rated=taxa.loc[taxa.EIV_NUTRI_exact.notna()].sort_values("speciesKey",kind="stable")
    if len(rated)!=N_SCORED:
        raise RuntimeError("Source 4011 numeric accepted vascular plant keys changed")
    ix={k:i for i,k in enumerate(rated.speciesKey)}
    vals=rated.EIV_NUTRI_exact.to_numpy(float)
    if not np.isfinite(vals).all() or np.any(vals<1) or np.any(vals>9):
        raise RuntimeError("Original published EIV scores invalid or missing")
    gnames=[]
    for name in rated.gbif_canonical_name:
        parts=str(name).split()
        if len(parts)!=2 or not all(re.fullmatch(r"[A-Za-z][A-Za-z-]+",p) for p in parts):
            raise RuntimeError("A source-score species has no true accepted binomial")
        gnames.append(parts[0].casefold())
    genera=defaultdict(list)
    for i,g in enumerate(gnames):genera[g].append(i)
    exchangeable=[np.array(z,int) for z in genera.values() if len(z)>=2]
    old=load_facets(a.oldfacets,"GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json",
                    OLD_FACETS_SHA,2808)
    new=load_facets(a.newfacets,"MERGED_ORIGINAL_30H7_TWO_PERIOD_SIX_FAMILY_360_SOURCE_FACETS.json",
                    NEW_FACETS_SHA,360)
    if set(old)&set(new):raise RuntimeError("Original old/new GBIF sources duplicate H7-period")
    records={**old,**new}
    A=np.zeros((len(frame),len(rated)),float)
    n_source_periods=0
    for i,row in enumerate(frame.itertuples(index=False)):
        for t in PERIODS:
            k=(str(row.h),t,"Plantae")
            if k not in records:raise RuntimeError("Original 183 H7×two period Plantae absent "+str(k))
            f=records[k]
            present=[ix[str(key)] for key,val in f["source_speciesKey_counts"].items()
                     if int(val)>0 and str(key) in ix]
            if len(present)<10:
                raise RuntimeError("Scored original source site-period unexpectedly <10 plant species "+str(k))
            A[i,present]+=0.5/len(present)
            n_source_periods+=1
    if n_source_periods!=366 or not np.allclose(A.sum(axis=1),1.,atol=1e-12,rtol=0):
        raise RuntimeError("Each source H7 must average two complete plant score periods")
    original=A@vals
    error=float(np.max(np.abs(original-frame[PREDICTOR].to_numpy(float))))
    if error>1e-10:raise RuntimeError("Genuine original 183 source CWM NUTRI reconstitution FAILED: "+str(error))
    return frame,A,vals,exchangeable,{
      "n_eiv_numeric_species":len(vals),
      "n_source_plant_facets_used":n_source_periods,
      "n_exact_accepted_genera":len(genera),
      "n_scored_species_exchangeable_within_genera":int(sum(len(v) for v in exchangeable)),
      "n_genera_multiple_scored_species":len(exchangeable),
      "plantCWM_replay_max_abs_error":error,
      "original_183_source_SHA256":PANEL_SOURCE_SHA,
      "original_old_facets_SHA256":OLD_FACETS_SHA,
      "original_new_facets_SHA256":NEW_FACETS_SHA,
      "original_published_EIV_SHA256":SCORE_H7_SOURCE_SHA
    }

def residualize_setup(frame,A,trait):
    X,labels,center=features(frame)
    j=labels.index(PREDICTOR)
    root=np.sqrt(weights(frame))
    X0=np.delete(X,j,axis=1)
    if np.linalg.matrix_rank(X0)!=X0.shape[1]:
        raise RuntimeError("Real H6-FE and 13-controls reduced model not full rank")
    Q,_=np.linalg.qr(root[:,None]*X0,mode="reduced")
    yr=root*frame[Y_PRIMARY].to_numpy(float)
    yres=yr-Q@(Q.T@yr)
    wa=A*root[:,None]
    source_x=wa@trait
    xres=source_x-Q@(Q.T@source_x)
    y2=float(yres@yres);x2=float(xres@xres)
    if y2<1e-12 or x2<1e-12:
        raise RuntimeError("Observed source residuals have zero variation")
    obs=float((yres@xres)/np.sqrt(x2*y2))
    source_sd=float(frame[PREDICTOR].std(ddof=0))
    fwl_beta=source_sd*float((xres@yres)/(xres@xres))
    if abs(fwl_beta-ORIGINAL_BETA)>1e-8:
        raise RuntimeError("Independent QR source beta not identical original H6 FE coefficient "+str(fwl_beta))
    return Q,wa,yres,y2,obs,fwl_beta,{"original_source_plant_N_FWL_beta_1sd":fwl_beta,
       "observed_H6_FE_13controls_weighted_partial_r":obs,
       "observed_H6_FE_13controls_partial_R2":obs**2,
       "source_original_183_plant_CWM_full_cohort_sd":source_sd,
       "unconditional_source_13_controls_and_40H6_FE_kept_fixed":True}

def permutation_group(A,Q,yres,y2,traits,groups,ndraw,rng):
    vals=np.empty(ndraw,float)
    batch=48
    for offset in range(0,ndraw,batch):
        stop=min(offset+batch,ndraw)
        B=np.empty((len(traits),stop-offset),float)
        for col in range(stop-offset):
            s=traits.copy()
            if groups is None:s=rng.permutation(traits)
            else:
                for g in groups:
                    s[g]=rng.permutation(traits[g])
            B[:,col]=s
        Z=A@B
        Z-=Q@(Q.T@Z)
        d=np.sum(Z*Z,axis=0)
        if not np.all(d>1e-12):
            raise RuntimeError("Trait label shuffle created unidentified zero site gradient")
        vals[offset:stop]=(yres@Z)/np.sqrt(y2*d)
    return vals

def summary_null(a,observed):
    extreme=int(np.count_nonzero(np.abs(a)>=abs(observed)-1e-12))
    return {"n_original_label_shuffles":len(a),
      "two_sided_permutation_p":float((extreme+1)/(len(a)+1)),
      "null_exceeding_observed_abs_partial_r":extreme,
      "null_q025_median_q975":np.quantile(a,[.025,.5,.975]).tolist(),
      "null_mean":float(a.mean()),"null_sd":float(a.std()),
      "null_not_causal_sewage_or_field_intervention":True}

def main(a):
    frame,A,vals,genera,qa=source_frame(a)
    Q,WA,yres,y2,observed,beta,obsqa=residualize_setup(frame,A,vals)
    rng=np.random.default_rng(SEED)
    full=permutation_group(WA,Q,yres,y2,vals,None,N_DRAW,rng)
    genus=permutation_group(WA,Q,yres,y2,vals,genera,N_DRAW,rng)
    out={"schema":"GESUI_183H7_40H6_PLANT_EIV_ORIGINAL_GLOBAL_AND_GENUS_LABEL_NULL_V1",
         "source_original_receipts":qa,
         "observed_source_FWL":obsqa,
         "shuffles":N_DRAW,"rng":"numpy PCG64 default_rng(20261115)",
         "all_scored_species_labels_unrestricted":summary_null(full,observed),
         "labels_shuffled_only_within_exact_GBIF_accepted_genus":summary_null(genus,observed),
         "new_taxonomic_scores_or_imputations":0,
         "original_eiv_missingness_and_2period_incidence_preserved":True,
         "original_butterfly_HNS_and_13_controls_and_H6_FE_untouched":True,
         "original_40H6_geographic_cohort_and_all_H7_unchanged":True,
         "posthoc_not_pre_registered_or_phylogenetically_independent":True,
         "not_actual_wastewater_N_or_plants_feeding_causal_mediation":True}
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/"183H7_40H6_PLANT_EIV_GLOBAL_GENUS_LABEL_NULL_RESULTS.json").write_text(
        json.dumps(out,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    report={"n_H7":len(frame),"n_H6":frame.h6.nunique(),
         "n_eiv_scored_species":qa["n_eiv_numeric_species"],
         "n_genus_exchangeable_species":qa["n_scored_species_exchangeable_within_genera"],
         "CWM_original_replay_max_abs":qa["plantCWM_replay_max_abs_error"],
         "source_beta_fwl":beta,"partial_r":observed,
         "all_species_trait_shuffle":out["all_scored_species_labels_unrestricted"],
         "within_genus_trait_shuffle":out["labels_shuffled_only_within_exact_GBIF_accepted_genus"]}
    print("GESUI_PLANT_EIV_LABEL_NULL_183_H7_40_H6="+json.dumps(
        report,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source183",required=True,type=Path)
    p.add_argument("--oldfacets",required=True,type=Path)
    p.add_argument("--newfacets",required=True,type=Path)
    p.add_argument("--eiv",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    main(p.parse_args())
