#!/usr/bin/env python3
"""Posthoc robustness of original ChemPop EDF after matched 3 river nutrients.

LOQ1/2/3 source substitutions and one-site-one-response sensitivity, all
from identical EIDC v2 sample/date sources. No original biological rows
redistributed; refuse alternate-LOQ analyses if source coverage differs by
more than a prespecified 5%, and compare strict common cohorts when it does.
Original main benchmark effect must be reproduced before new diagnostics.
"""
from __future__ import annotations

import argparse
import collections
import json
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import wq_adjusted_edf as m

MAIN_EDF_CONDITIONAL_BETA=0.0070911445162768275
MAIN_EDF_AFTER_NUTRIENTS_REGION_MSE_GAIN=-0.5185647610452906
EXPECTED_N_MATCHED=5000
EXPECTED_N_SITE=591
EXPECTED_N_SITEYEAR=2624
MIN_ALTERNATE_COHORT_FRACTION=.95


def source_and_primary():
    """Single exact original source construction used by all three LOQ modes."""
    with tempfile.TemporaryDirectory(prefix="gesui_chempop_loq_") as temp:
        d=Path(temp)
        arch=m.source.retrieve(m.source.EIDC_URL,d/"original.zip",m.source.SOURCE_DOWNLOAD_LIMIT)
        score,score_info=m.source.score_reference(d)
        if score_info["SHA256"]!=m.EXPECTED_SCORE_SHA or len(score)!=92:
            raise RuntimeError("Published UK BMWP official source fingerprint drift")
        with zipfile.ZipFile(d/"original.zip") as z:
            files=z.namelist()
            if len(files)!=326:raise RuntimeError("Original publisher v2 ZIP members drift")
            raw,weights,orders,sourceinfo=m.original.construct_eligible(z,files,score)
            keys={(r["region"],r["site"],r["sample_date"]) for r in raw}
            wq,src=m.chemical_summaries(z,files,keys)
            pair,gate=m.source_join_raw(raw,weights,wq)
    if (len(pair)!=EXPECTED_N_MATCHED
        or gate["matched_distinct_sites"]!=EXPECTED_N_SITE
        or gate["matched_distinct_site_year"]!=EXPECTED_N_SITEYEAR
        or not gate["source_gate_PASS"]):
        raise RuntimeError("Source benchmark 5000/591/2624 NP original sample cohort changed")
    return pair,score,{"original_EIDC_archive_SHA256":arch["SHA256"],
                      "source_BMWP_reference_SHA256":score_info["SHA256"],
                      "original_three_chemistry_source_files":src}


def benchmark(pair,score):
    d,qa=m.site_year_matched(pair,score)
    full=m.cluster_EDF_coeff(d,True)["coefficient"]["EDF_log1p_std"]["beta"]
    if abs(full-MAIN_EDF_CONDITIONAL_BETA)>1e-9:
        raise RuntimeError("Original paired full 3-chemical EDF BMWP adjusted coefficient changed")
    external=m.region_transfer_four_models(d)
    if external["status"]!="REGION_DISJOINT_FOUR_NESTED_MODELS":
        raise RuntimeError("Original 7 region holdout source support drift")
    gain=external["contrasts"]["EDF_after_nutrients"]["MSE_additional_gain_percent"]
    if abs(gain-MAIN_EDF_AFTER_NUTRIENTS_REGION_MSE_GAIN)>1e-7:
        raise RuntimeError("Original paired 3-chemical EDF added prediction changed")
    return d,qa,{
        "beta_EDF_after_3chem":full,
        "EDF_added_after_chem_region_loss_MSE_gain_pct":gain,
        "all_original_7_regions":len(external["valid_regions"]),
        "source_reconstruction_verified":True,
    }


def source_only_alt_subset(pair,variant):
    allowed=[]
    status=collections.Counter()
    for row,w,chem in pair:
        replacements={}
        for analyte in m.NUTRIENTS:
            newval=chem[analyte].get(variant)
            if newval is None or not np.isfinite(newval) or newval<0:
                status["missing_or_invalid_"+analyte]+=1
                break
            repl=dict(chem[analyte])
            repl["primary"]=float(newval)
            replacements[analyte]=repl
        else:
            allowed.append((row,w,replacements))
    ratio=len(allowed)/len(pair) if pair else 0
    decision={"alternative":variant,
              "n_alt_fully_scored_original_samples":len(allowed),
              "prior_original_matched_samples":len(pair),
              "fraction_same_source_kept":ratio,
              "missing_reasons":dict(status),
              "prespecified_min_cohort_fraction":MIN_ALTERNATE_COHORT_FRACTION,
              "strict_common_pairable":ratio>=MIN_ALTERNATE_COHORT_FRACTION}
    return allowed,decision


def identical_primary_on_subset(pair, allowed):
    """Compare LOQ2 and LOQ alternative ON IDENTICAL original samples."""
    keys=set((x[0]["site_region_id"],x[0]["sample_date"],x[0]["year"],x[0]["method"])
             for x in allowed)
    primary_subset=[q for q in pair
                    if (q[0]["site_region_id"],q[0]["sample_date"],q[0]["year"],q[0]["method"]) in keys]
    if len(primary_subset)!=len(allowed):
        raise RuntimeError("Original sample ID method+date not unique in paired LOQ cohort")
    return primary_subset


def summarize_models(df):
    V=m.VIF_for_original_controls(m.full_design(df,True))
    if max(V.values())>m.SOURCE_GATES["max_full_design_VIF_EDF_or_one_nutrient"]:
        return {"status":"SOURCE_STOP_COLLINEAR","VIF":V}
    c=m.cluster_EDF_coeff(df,True)
    o=m.cluster_EDF_coeff(df,False)
    transfer=m.region_transfer_four_models(df)
    return {"status":"SOURCE_COHORT_VALID_MODEL",
            "n_site_year":len(df),
            "n_sites":df.site_region_id.nunique(),
            "n_regions":df.region.nunique(),
            "VIF":V,
            "EDF_with_nutrients":c["coefficient"]["EDF_log1p_std"],
            "EDF_without_nutrients":o["coefficient"]["EDF_log1p_std"],
            "four_model_region_transfer":transfer.get("contrasts"),
            "region_transfer_status":transfer["status"]}


def one_site_only(df):
    """Give each source physical river site exactly one equal-weight record."""
    numeric=(tuple(m.original.COVARIATES)+("edf_pct","edf_log1p")+
             tuple("chem_"+k for k in m.NUTRIENTS)+("y_incidence",))
    rows=[]
    for sid,g in df.groupby("site_region_id",sort=True):
        if abs(g.edf_pct.max()-g.edf_pct.min())>1e-10:
            raise RuntimeError("Source EDF unexpectedly changes within site")
        row={col:float(g[col].mean()) for col in numeric}
        row.update({"site_region_id":sid,
            "site":str(g.site.iloc[0]),"region":str(g.region.iloc[0]),
            "method":str(g.method.value_counts().index[0]),
            "n_site_year":len(g)})
        rows.append(row)
    frame=pd.DataFrame(rows)
    if frame.site_region_id.nunique()!=EXPECTED_N_SITE:
        raise RuntimeError("Same original 591 source sites changed")
    return frame


def make(args):
    out={"schema":"GESUI_20261120_CHEMPOP_NP_EDF_LOQ_AND_SITELEVEL_V1",
         "status":"SOURCE_STOP",
         "posthoc_after_original_EDF_after3NP_prediction_was_negative":True,
         "no_original_WQ_or_taxon_raw_rows_republished":True,
         "not_causal":True}
    try:
        pair,score,sourceinfo=source_and_primary()
        out["source"]=sourceinfo
        primary,qa,orig=benchmark(pair,score)
        out["main_source_QA"]=qa
        out["original_main_exact_numeric_replay"]=orig
        site=one_site_only(primary)
        out["single_original_site_per_row"]=summarize_models(site)
        out["single_original_site_per_row"]["n_original_site_years_collapsed"]=len(primary)
        out["alternative_original_LOQ"]= {}
        for variant in ("lower_LOQ_version","upper_LOQ_version"):
            subset,gate=source_only_alt_subset(pair,variant)
            rec={"source_gate":gate}
            if not gate["strict_common_pairable"]:
                rec["status"]="SOURCE_STOP_ALTERNATE_LOQ_TOO_MUCH_MISSINGNESS"
            else:
                originals=identical_primary_on_subset(pair,subset)
                origdf,_=m.site_year_matched(originals,score)
                altdf,_=m.site_year_matched(subset,score)
                if (len(origdf)!=len(altdf)
                    or list(origdf.site_region_id)!=list(altdf.site_region_id)
                    or list(origdf.year)!=list(altdf.year)):
                    raise RuntimeError("Original LOQ2 and alternative source siteyear cohorts changed")
                rec["status"]="SOURCE_COHORT_VALID_POSTHOC_LOQ"
                rec["LOQ2_on_same_subset"]=summarize_models(origdf)
                rec["alternative_LOQ_on_same_subset"]=summarize_models(altdf)
            out["alternative_original_LOQ"][variant]=rec
        out["status"]="POSTHOC_LOQ_AND_SITE_ONLY_COMPLETED"
    except Exception as ex:
        out.update({"status":"SOURCE_OR_MODEL_STOP","error":type(ex).__name__+": "+str(ex)[:650]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    excerpt=lambda r:(None if not isinstance(r,dict) else {
        "EDF_after_NP":r.get("EDF_with_nutrients"),
        "region_transfer":r.get("four_model_region_transfer"),
        "status":r.get("status"),"n_sites":r.get("n_sites")})
    print("GESUI_CHEMPOP_NP_LOQ_SITELEVEL_STRESS="+json.dumps({
        "status":out["status"],"error":out.get("error"),
        "main":out.get("original_main_exact_numeric_replay"),
        "site_only":excerpt(out.get("single_original_site_per_row")),
        "LOQ":{k:{"source":v["source_gate"],
                 "primary":excerpt(v.get("LOQ2_on_same_subset")),
                 "alternative":excerpt(v.get("alternative_LOQ_on_same_subset")),
                 "status":v["status"]} for k,v in out.get("alternative_original_LOQ",{}).items()}
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
