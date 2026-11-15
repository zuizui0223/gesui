#!/usr/bin/env python3
"""External state-block EPA site×taxon trait-model robustness; source only.

Same untouched three EPA NRSA eras and strict original trait-quality samples
as gesui #99. No new taxonomic/field values are synthesized or imputed.
State blocks and target HUC8 are both absent from each training fold.
Focal interaction: original EPA clinger habit within EPT order × XEMBED.
All three penalties 0.005, 0.025, 0.100 fixed BEFORE looking at new outcomes.
Separate full-scored-taxon and shared-training-taxon conditional CE diagnostics.

This is posthoc observational source geography, not true unseen-continental
transfer, causal trait-mediated selection, or sewage-derived nutrient effects.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
import joint_multinomial as model

CYCLES = model.ERAS
PENALTIES=(0.005,0.025,0.100)
N_STATE_FOLDS=5
BOOT=2000
SEED=20261115

def fold_states(d):
    counts=d.state.astype(str).value_counts()
    if (counts<=0).any() or len(counts)<20:
        raise ValueError("Original state source support unexpectedly small")
    ordered=sorted(((str(k),int(v)) for k,v in counts.items()),
                   key=lambda kv:(-kv[1],kv[0]))
    bins=[{"states":[],"n":0} for _ in range(N_STATE_FOLDS)]
    for name,n in ordered:
        index=min(range(N_STATE_FOLDS),key=lambda k:(bins[k]["n"],k))
        bins[index]["states"].append(name)
        bins[index]["n"]+=n
    code={name:i for i,b in enumerate(bins) for name in b["states"]}
    fold=d.state.astype(str).map(code).to_numpy(dtype=int)
    if len(set(fold))!=N_STATE_FOLDS:
        raise ValueError("Deterministic state block partition failed")
    return fold,bins

def heldout_losses(F,X,Q,logprev,B):
    """Nats/standardized sample, with probabilities over all known EPA taxa."""
    logits=logprev[None,:]+(X@B)@Q.T
    all_logp=logits-logsumexp(logits,axis=1)[:,None]
    return -(F*all_logp).sum(axis=1),logits

def conditional_common(F,logits,seen):
    """Secondary: evaluate only original taxa recorded in training.

    Excludes unseen taxa from a CONDITIONAL log-loss and re-normalizes both
    the observed composition and the predicted probabilities. It does not
    claim accuracy for any novel taxa; report excluded proportions.
    """
    raw=F[:,seen]
    mass=raw.sum(axis=1)
    valid=mass>1e-12
    result=np.full(len(F),np.nan)
    if valid.any():
        obs=raw[valid]/mass[valid,None]
        lp=logits[valid][:,seen]
        lp=lp-logsumexp(lp,axis=1)[:,None]
        result[valid]=-(obs*lp).sum(axis=1)
    return result,valid,mass

def huc_equal_summary(d,lossbase,lossplus,mask=None):
    H=d.huc8.astype(str).to_numpy()
    valid=np.ones(len(d),dtype=bool) if mask is None else np.asarray(mask,bool)
    if valid.sum()<100:
        return {"status":"INSUFFICIENT_SUPPORTED_SHARED_TAXA",
                "source_eligible_n":int(valid.sum())}
    groups,idx=np.unique(H[valid],return_inverse=True)
    if len(groups)<70:
        return {"status":"INSUFFICIENT_HUC8",
                "source_eligible_n":int(valid.sum()),"HUC8":len(groups)}
    num=np.bincount(idx)
    a=np.bincount(idx,weights=lossbase[valid])/num
    b=np.bincount(idx,weights=lossplus[valid])/num
    delta=a-b
    observed=float(delta.mean())
    rng=np.random.default_rng(SEED)
    samples=rng.integers(0,len(groups),size=(BOOT,len(groups)))
    ds=delta[samples].mean(axis=1)
    meanbase=a[samples].mean(axis=1)
    return {
      "status":"OBSERVATIONAL_CV_COMPLETE",
      "source_eligible_n":int(valid.sum()),"HUC8":len(groups),
      "cross_entropy_base_nats":float(a.mean()),
      "cross_entropy_trait_nats":float(b.mean()),
      "delta_nats_per_scored_site":observed,
      "gain_percent":float(100*observed/a.mean()),
      "fixed_error_HUC8_boot95_nats":list(map(float,np.quantile(ds,[.025,.975]))),
      "fixed_error_HUC8_boot95_percent":list(map(float,
                          np.quantile(100*ds/meanbase,[.025,.975]))),
      "no_refit_bootstrap":True
    }

def evaluate_era(d, F, Q, typ, penalty):
    state_folds, state_blocks=fold_states(d)
    group=d.huc8.astype(str).to_numpy()
    CEs={"M0":np.full(len(d),np.nan),"M1":np.full(len(d),np.nan)}
    shared={"M0":np.full(len(d),np.nan),"M1":np.full(len(d),np.nan)}
    common_valid=np.zeros(len(d),bool)
    mass_to_known=np.full(len(d),np.nan)
    fold_receipts=[]
    for fold in range(N_STATE_FOLDS):
        test=state_folds==fold
        target_huc=set(group[test])
        training=(state_folds!=fold)&(~np.isin(group,list(target_huc)))
        if training.sum()<160 or test.sum()<30:
            raise ValueError(f"State fold insufficient original source support {fold}, {training.sum()},{test.sum()}")
        if set(d.loc[training,"state"]) & set(d.loc[test,"state"]):
            raise ValueError("State group train/test leakage")
        if set(group[training]) & set(group[test]):
            raise ValueError("Source hydrological HUC8 train/test leakage")
        X,T=model.covariates(d.loc[training],d.loc[test])
        logprev,w=model.prevalence(F[training],group[training])
        B0,s0=model.optimize(F[training],X,Q,logprev,w,False,lam=penalty)
        B1,s1=model.optimize(F[training],X,Q,logprev,w,True,lam=penalty)
        l0,z0=heldout_losses(F[test],T,Q,logprev,B0)
        l1,z1=heldout_losses(F[test],T,Q,logprev,B1)
        CEs["M0"][test]=l0;CEs["M1"][test]=l1
        seen=F[training].sum(axis=0)>0
        novel=F[test][:,~seen].sum(axis=1)
        cc0,val0,mass=conditional_common(F[test],z0,seen)
        cc1,val1,mass1=conditional_common(F[test],z1,seen)
        if not np.array_equal(val0,val1) or not np.allclose(mass,mass1):
            raise ValueError("Baseline / interaction source common-taxon evaluations differ")
        shared["M0"][test]=cc0;shared["M1"][test]=cc1
        common_valid[test]=val0
        mass_to_known[test]=mass
        fold_receipts.append({
          "fold":fold,"heldout_states":sorted(d.loc[test,"state"].astype(str).unique().tolist()),
          "test_sites":int(test.sum()),"train_sites":int(training.sum()),
          "test_HUC8":len(target_huc),"train_HUC8":len(set(group[training])),
          "source_disjoint_states":True,"source_disjoint_HUC8":True,
          "additional_source_HUC8_purged_from_other_state_training":
                    int(((state_folds!=fold)&np.isin(group,list(target_huc))).sum()),
          "number_known_scored_EPT_taxa_in_training":int(seen.sum()),
          "taxon_roster_all_known_source_only":int(len(seen)),
          "heldout_scored_taxon_mass_unseen_in_train":float(novel.mean()),
          "heldout_sites_with_no_shared_taxon":int((~val0).sum()),
          "trait_cross_interaction_theta":float(B1[0,0]),
          "base_fit":s0,"focal_fit":s1
        })
    if not all(np.isfinite(x).all() for x in CEs.values()):
        raise ValueError("Full genuine EPA state heldout loss missing")
    result={
        "state_block_design":state_blocks,
        "original_source_n_sites":len(d),
        "state_holdout_fold_receipts":fold_receipts,
        "same_original_sample_both_models":True,
        "training_state_and_HUC8_leakage":False,
        "theta_all_five_folds_negative":all(f["trait_cross_interaction_theta"]<0 for f in fold_receipts),
        "full_taxon_roster_from_original_publisher_source":True,
        "test_taxa_ever_unseen_training_are_not_imputed_as_absent":True,
        "fraction_heldout_response_mass_from_taxa_unseen_in_train":float(np.mean(1-mass_to_known)),
        "fraction_heldout_sites_with_no_shared_taxa":float(np.mean(~common_valid)),
        "all_scored_taxa_unconditional":huc_equal_summary(d,CEs["M0"],CEs["M1"]),
        "only_taxa_observed_in_train_conditional":huc_equal_summary(d,shared["M0"],shared["M1"],common_valid),
        "negative_control_label":"Both metrics are same taxon roster and specimen source; not independent observations",
        "penalty":penalty,"outcome":typ
    }
    return result

def run(out):
    sites,counts=model.source()
    results={"schema":"GESUI_EPA3ERA_STATE_HUC8_PURGED_PENALTY_SENSITIVITY_V1",
             "source":"original EPA native 3 survey cycles/2209 eligible site IDs",
             "parent_PR99_huc8_holdout_results_are_unchanged":True,
             "analysis_is_posthoc_exploratory":True,
             "source_field_data_untouched":True,
             "all_source_nulls_and_negative_results_preserved":True,
             "true_field_WWTP_chemical_attribution":False,
             "model_class":"fixed-penalty conditional multinomial fourth-corner, not a JSDM",
             "source_penalty_grid":PENALTIES,
             "state_group_holdout":N_STATE_FOLDS,
             "nested_HUC8_leakage_purge":True,
             "results":{}}
    for era in CYCLES:
        d,Y,Fcount,Finc,Q,meta=model.matrices_for_era(sites,counts,era)
        results["results"][era]={"source_quality":meta,"by_penalty":{}}
        for penalty in PENALTIES:
            k=f"{penalty:.3f}"
            value={
               "count":evaluate_era(d,Fcount,Q,"individual_count_composition",penalty),
               "incidence":evaluate_era(d,Finc,Q,"equal_scored_taxon_incidence",penalty)
            }
            results["results"][era]["by_penalty"][k]=value
            print("STATE_BLOCK_TAXON",era,k,json.dumps({
               name:{
                  "CE_full":v["all_scored_taxa_unconditional"]["gain_percent"],
                  "CE_full95":v["all_scored_taxa_unconditional"]["fixed_error_HUC8_boot95_percent"],
                  "CE_common":v["only_taxa_observed_in_train_conditional"].get("gain_percent"),
                  "CE_common95":v["only_taxa_observed_in_train_conditional"].get("fixed_error_HUC8_boot95_percent"),
                  "unseen_taxa_fraction":v["fraction_heldout_response_mass_from_taxa_unseen_in_train"],
                  "theta_all_folds_negative":v["theta_all_five_folds_negative"],
                  "thetas":[r["trait_cross_interaction_theta"] for r in v["state_holdout_fold_receipts"]]
               } for name,v in value.items()
            },sort_keys=True,allow_nan=False),flush=True)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(results,indent=2,sort_keys=True,allow_nan=False))
    print("STATE_HUC8_SOURCE_AUDIT_DONE",str(out),flush=True)
    return results

def self_test():
    import pandas as pd
    rng=np.random.default_rng(71)
    n=200
    states=np.array([f"S{i}" for i in range(26) for _ in range(8)])[:n]
    d=pd.DataFrame({"state":states,"huc8":[str(i//3) for i in range(n)]})
    f,b=fold_states(d)
    if len(set(f))!=5:raise AssertionError("State-block synthetic failed")
    for z in range(5):
        tr=(f!=z)&(~np.isin(d.huc8,d.huc8[f==z]))
        te=f==z
        assert not (set(d.huc8[tr])&set(d.huc8[te]))
        assert not (set(d.state[tr])&set(d.state[te]))
    F=np.stack([rng.dirichlet(np.ones(11)) for _ in range(n)])
    logits=rng.normal(size=(n,11));seen=np.r_[np.ones(8,bool),np.zeros(3,bool)]
    loss,valid,m=conditional_common(F,logits,seen)
    assert valid.all() and np.isfinite(loss).all()
    print("STATE_BLOCK_HUC8_ZERO_LEAK_TEST_PASS",flush=True)

if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--output",type=Path,default=Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_ept_state_transfer"/"STATE_AND_PENALTY_RESULTS.json")
    a.add_argument("--self-test",action="store_true")
    args=a.parse_args()
    if args.self_test:self_test()
    else:run(args.output)
