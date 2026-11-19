#!/usr/bin/env python3
"""Posthoc country-heldout transfer of W->plant NUTRI genus components.

Source quality gates are based on country and true H6 counts only. Each test
country removes every overlapping H6 from its training countries, and the
existing P_GEOG / original country_method one-hot encoding is train-only.
No source selected based on sign. Unknown target country_method categories
naturally use zeros, not a pooled/masked target category. Point/CI are
descriptive fixed-error country-block resamples, not a new independent study.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import h6_model_if_gate as m
import plant_W_genus_label_null as g


MIN_COUNTRY_TEST_ROWS = 10
MIN_COUNTRY_H6 = 3
MIN_TRAIN_ROWS = 100
MIN_TRAIN_H6 = 35
MIN_TRAIN_COUNTRIES = 4
MIN_ELIGIBLE_COUNTRIES = 5
MIN_MULTI_H7_H6_FOR_LOCAL_WITHIN_H6 = 40
MIN_COUNTRY_TEST_TRAIN_NONOVERLAP = True
N_COUNTRY_BOOT = 2000


def original_decomposition(df, values, genus, inc, rec):
    gbar,_,groups = g.decompose(values,genus)
    result=df.copy()
    for mode, mat in (("incidence",inc),("records",rec)):
        result[mode+"_total"]=mat@values
        result[mode+"_between"]=mat@gbar
        result[mode+"_within"]=mat@(values-gbar)
        if np.max(np.abs(result[mode+"_total"]-
                    result[mode+"_between"]-result[mode+"_within"]))>1e-10:
            raise RuntimeError("Same immutable source genotype/trait partition lost")
    return result


def aggregate_H6_losses(test, y, pred0, pred1):
    h=test.h6.astype(str).to_numpy()
    errs0=(y-pred0)**2
    errs1=(y-pred1)**2
    a=[];b=[]
    for key in sorted(set(h)):
        idx=(h==key)
        a.append(float(np.mean(errs0[idx])))
        b.append(float(np.mean(errs1[idx])))
    a=np.asarray(a);b=np.asarray(b)
    if a.mean()<1e-12:
        raise RuntimeError("Trivial zero-test-MSE")
    return {"H6_equal_baseline_MSE":float(a.mean()),
            "H6_equal_added_W_MSE":float(b.mean()),
            "H6_equal_W_added_gain_pct":float(100*(a.mean()-b.mean())/a.mean())}


def country_out(df,eligible,country,outcome):
    rows_country=df.country_holdout.astype(str)==country
    test=df.loc[rows_country].copy().reset_index(drop=True)
    h6_test=set(test.h6.astype(str))
    train=df.loc[~rows_country &
                 (~df.h6.astype(str).isin(h6_test))].copy().reset_index(drop=True)
    source_n=(len(test),test.h6.nunique(),len(train),train.h6.nunique(),
              train.country_holdout.nunique())
    source={
        "test_n_H7period":source_n[0],
        "test_true_H6":source_n[1],
        "train_n_H7period":source_n[2],
        "train_true_H6":source_n[3],
        "train_other_countries":source_n[4],
        "test_H6_all_removed_from_training":bool(not(set(train.h6)&set(test.h6))),
        "no_target_outcome_used_for_source_quality_selection":True
    }
    if (source_n[0]<MIN_COUNTRY_TEST_ROWS or source_n[1]<MIN_COUNTRY_H6 or
        source_n[2]<MIN_TRAIN_ROWS or source_n[3]<MIN_TRAIN_H6 or
        source_n[4]<MIN_TRAIN_COUNTRIES):
        return {"status":"SOURCE_SUPPORT_STOP", "source":source}
    y=test[outcome].to_numpy(float)
    baseline=m.predict_ridge(train,test,outcome,m.P_GEOG)
    withW=m.predict_ridge(train,test,outcome,m.P_GEOG+["W"])
    if not (np.isfinite(baseline).all() and np.isfinite(withW).all()):
        raise RuntimeError("Out-of-country predictions not finite")
    return {"status":"HELDOUT_COUNTRY_AND_H6_DISJOINT",
            "source":source,
            **aggregate_H6_losses(test,y,baseline,withW)}


def make(args):
    df,values,genus,inc,rec,err,gate=g.load_immutable_source(args)
    full=original_decomposition(df,values,genus,inc,rec)
    groups=full.groupby("h6").h.nunique()
    multi=groups[groups>=2]
    within_source={
        "n_original_H6":int(len(groups)),
        "H6_with_two_or_more_distinct_H7":int(len(multi)),
        "distinct_H7_in_multih7_H6":int(full[full.h6.isin(multi.index)].h.nunique()),
        "source_only_gate_at_least_40_H6":bool(len(multi)>=MIN_MULTI_H7_H6_FOR_LOCAL_WITHIN_H6),
        "model_execution_if_gate_fails":"STOP"
    }
    countries=sorted(full.country_holdout.astype(str).unique().tolist())
    if len(countries)!=11:raise RuntimeError("Original 11 country holdout groups changed")
    result={
        "schema":"GESUI_EU451_W_PLANT_N_GENUS_COUNTRY_HELDOUT_V1",
        "exploratory_after_within_genus_positive":True,
        "source_only_before_opening_country_outcome":True,
        "not_actual_wastewater_N_or_causal":True,
        "quality_thresholds":{
            "min_country_test_rows":MIN_COUNTRY_TEST_ROWS,
            "min_country_test_H6":MIN_COUNTRY_H6,
            "min_train_rows":MIN_TRAIN_ROWS,
            "min_train_H6":MIN_TRAIN_H6,
            "min_train_country_groups":MIN_TRAIN_COUNTRIES,
            "min_eligible_country_groups_for_summary":MIN_ELIGIBLE_COUNTRIES,
            "min_H6_for_local_within_parent_W_test":MIN_MULTI_H7_H6_FOR_LOCAL_WITHIN_H6,
        },
        "original_n_H7period":len(df),
        "original_n_H6":df.h6.nunique(),
        "original_country_groups":countries,
        "within_H6_source_only":within_source,
        "by_country":{},
        "summary":{},
    }
    # Authorize entire country cohort BEFORE reading its endpoint predictions.
    source_gate={}
    for country in countries:
        test=full[full.country_holdout.astype(str)==country]
        train=full[(full.country_holdout.astype(str)!=country) &
                   (~full.h6.astype(str).isin(set(test.h6.astype(str))))]
        eligible_country=(len(test)>=MIN_COUNTRY_TEST_ROWS and
                          test.h6.nunique()>=MIN_COUNTRY_H6 and
                          len(train)>=MIN_TRAIN_ROWS and
                          train.h6.nunique()>=MIN_TRAIN_H6 and
                          train.country_holdout.nunique()>=MIN_TRAIN_COUNTRIES)
        source_gate[country]={
           "eligible":bool(eligible_country),
           "test_n":len(test),"test_H6":test.h6.nunique(),
           "train_n":len(train),"train_H6":train.h6.nunique()}
    authorized=[c for c in countries if source_gate[c]["eligible"]]
    result["source_gate_by_country"]=source_gate
    result["source_eligible_countries"]=authorized
    result["country_level_primary_criterion_gate_pass"]=len(authorized)>=MIN_ELIGIBLE_COUNTRIES

    for label in ("incidence_total","incidence_between","incidence_within",
                  "records_total","records_between","records_within"):
        result["by_country"][label]={}
        for country in countries:
            if not source_gate[country]["eligible"]:
                result["by_country"][label][country]={"status":"SOURCE_SUPPORT_STOP",
                                                      "source":source_gate[country]}
            else:
                result["by_country"][label][country]=country_out(full,True,country,label)
        seen=[result["by_country"][label][c] for c in authorized]
        if len(seen)<MIN_ELIGIBLE_COUNTRIES:
            result["summary"][label]={
                "status":"STOP_UNDER_5_COUNTRIES",
                "n_source_eligible_countries":len(seen)}
            continue
        a=np.asarray([x["H6_equal_baseline_MSE"] for x in seen],float)
        b=np.asarray([x["H6_equal_added_W_MSE"] for x in seen],float)
        rng=np.random.default_rng(20261119)
        ix=rng.integers(0,len(a),size=(N_COUNTRY_BOOT,len(a)))
        gains=100*(a[ix].mean(axis=1)-b[ix].mean(axis=1))/a[ix].mean(axis=1)
        result["summary"][label]={
            "status":"COUNTRY_AND_H6_DISJOINT_PREDICTION",
            "n_source_eligible_countries":len(seen),
            "country_equal_MSE0":float(a.mean()),
            "country_equal_MSE1":float(b.mean()),
            "country_equal_W_gain_pct":float(100*(a.mean()-b.mean())/a.mean()),
            "fixed_loss_country_block_bootstrap95":list(map(float,np.quantile(gains,[.025,.975]))),
            "countries_with_positive_added_W_gain":int(np.sum(a>b)),
            "point_estimate_is_not_trained_after_heldout_country":True,
            "interval_does_not_cover_model_refitting_or_country_sampling":True
        }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_PLANT_W_COUNTRY_TRANSFER="+json.dumps({
        "source_gate":source_gate,
        "within_H6_source_gate":within_source,
        "authorized_countries":authorized,
        "summary":result["summary"]
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--eiv",required=True,type=Path)
    p.add_argument("--facets",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
