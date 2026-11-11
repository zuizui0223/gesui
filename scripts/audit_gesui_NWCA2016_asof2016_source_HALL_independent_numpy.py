#!/usr/bin/env python3
"""Separate NumPy verification of 2016 historical-permit NWCA 212 wetland results.

Same original 2016 source-frozen 50 NEAR/162 FAR source and CSV sha, but
compute coefficients/country-state CR1 directly using numeric SVD, not
original statsmodels formula. Source 2026 POTW permit date is not 2016 DMR.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from pathlib import Path
from gesui_NWCA2016_historicalyear_NPDES_near10_FAR20_actual_HALL import load
from audit_gesui_NWCA2011_2016_POTW_HALL_original_country_cluster_numpy import one

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
MAIN=D/"NWCA2016_ASOF2016_ECHO_REGISTERED_POTW_HALL_50NEAR162FAR_POSTHOC_QA.json"
SOURCE=D/"NWCA2016_2026_NPDES_PERMIT_EFF_ASOF2016_SOURCE_ONLY_NEAR10_FAR20.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify(root,dir):
    orig=json.loads((root/MAIN).read_text())
    source=json.loads((root/SOURCE).read_text())
    if orig["original_year2016_source_only_permit2026_near_far_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Source year2016 status not same as reported biological model")
    if not source["original_2016_historical_source_vegmetric_plant_H_ALL_NOT_READ"]:
        raise ValueError("Source population contaminated with H_ALL records")
    df,meta=load(root,dir)
    if Counter(df.NEAR)!=Counter({1:50,0:162}) or df.STATE_ECO.nunique()!=23:
        raise ValueError("Original fixed 50 near,162 far before H_ALL missing")
    expected=orig["original_2016_state_ecoregion_NWCA_plant_diversity_by_yearvalid_REGISTERED_POTW_source"]["HALL"]
    out={}
    for name,weighted,refname in (("site_equal",False,"site_equal_OLS"),
                                   ("survey_weighted",True,"original_NWCA_design_weight_WLS")):
        numeric=one(df,weighted)
        ref=expected[refname]
        if abs(numeric["beta_near_minus_far"]-ref["near_minus_far_beta"])>1e-8:
            raise ValueError("2016 historical-permit NumPy independently computed point β disagrees")
        for a,b in zip(numeric["state_cluster95"],ref["state_clustered_descriptive95"]):
            if abs(a-b)>2e-5:
                raise ValueError("2016 historical-permit NumPy state-cluster95 disagrees: "+
                                 repr(numeric["state_cluster95"])+" "+repr(ref["state_clustered_descriptive95"]))
        out[name]=numeric
    dropped={}
    for state in sorted(df.PSTLCODE.unique()):
        d=df[df.PSTLCODE!=state]
        dropped[str(state)]={"weighted":one(d,True)["beta_near_minus_far"],
                             "equal":one(d,False)["beta_near_minus_far"],
                             "n_source_sites_remaining":len(d)}
    return {
      "schema":"gesui_NWCA2016_ECHO_asof2016_POTW_actual_212site_numpy_OLS_WLS_cluster95_independent_QA_v01",
      "2016_original_HALL_50near162far_source_first_SHA256":sha(root/SOURCE),
      "2016_original_HALL_reporting_model_SHA256":sha(root/MAIN),
      "original_212_original_full_HALL_source_siteIDs_and_permit_dates_unchanged":True,
      "independent_2016_cluster_beta_and_95CI_vs_original_statsmodels":out,
      "n_original_US_states_dropped":len(dropped),
      "entire_state_omission_results":dropped,
      "all_state_deletion_count_weighted_negative":sum(r["weighted"]<0 for r in dropped.values()),
      "all_state_deletion_count_site_equal_negative":sum(r["equal"]<0 for r in dropped.values()),
      "weighted_delete_state_beta_minmax":[min(x["weighted"] for x in dropped.values()),max(x["weighted"] for x in dropped.values())],
      "site_equal_delete_state_beta_minmax":[min(x["equal"] for x in dropped.values()),max(x["equal"] for x in dropped.values())],
      "2016_current_ECHO_registered_permit_date_not_actual_2016_wastewater_release":True,
      "no_causal_sewage_benefit_or_harm_inferred":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    result=verify(root,args.source_dir)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"independent":result["independent_2016_cluster_beta_and_95CI_vs_original_statsmodels"],
        "states":result["n_original_US_states_dropped"],
        "weighted_negative_LOO":result["all_state_deletion_count_weighted_negative"],
        "weighted_range":result["weighted_delete_state_beta_minmax"],
        "equal_range":result["site_equal_delete_state_beta_minmax"]},indent=2))
if __name__=="__main__":main()
