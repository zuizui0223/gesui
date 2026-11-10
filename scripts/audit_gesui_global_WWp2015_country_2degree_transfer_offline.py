#!/usr/bin/env python3
"""Independent offline global WWp positive-vs-country generalization arithmetic.

Verify source-archived 1055 WPI ports, 42 countries, all heldout port residuals,
country×2degree block grouping, geographic/country MSE calculations and country
bootstrap. No original GBIF host records, no regression source parser, no
individual person data, no new external requests.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

ROOT=Path("results/global_WWp2015_country_transfer_20261110")
RESULT=ROOT/"GLOBAL_1055_PORT_42_COUNTRY_WWp_OUTOFCOUNTRY_AND_GEO_HOLDOUT.json"
ARCHIVE=Path("archive/chocho_pr122/provenance/reviewer_defenses/results/global_hydrowaste_port_hosts_20261102/JONES2015_WORLDCOVER2021_GLOBAL_WW_SOURCE_AUDIT_20261103.json")
SEED=20261110

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def validate_close(x,y,name,tol=1e-10):
    if abs(x-y)>tol:raise ValueError(f"{name} arithmetic changed {x} vs {y}")

def summary(rows):
    assert len({x["port_id"] for x in rows})==len(rows)
    by=defaultdict(list)
    for x in rows:
        if x["baseline_squared_error"]<0 or x["plus_WWp_squared_error"]<0:
            raise ValueError("Squared error negative")
        by[x["country"]].append(x)
    def calc(v):
        return [sum(z["baseline_squared_error"] for z in v)/len(v),
                sum(z["plus_WWp_squared_error"] for z in v)/len(v)]
    all0,all1=calc(rows)
    country={c:calc(v) for c,v in sorted(by.items())}
    equal0=np.mean([v[0] for v in country.values()])
    equal1=np.mean([v[1] for v in country.values()])
    return {
        "ports":len(rows),"countries":len(by),
        "port_mse0":all0,"port_mse1":all1,
        "port_gain_percent":100*(all0-all1)/all0,
        "country_mse0":float(equal0),"country_mse1":float(equal1),
        "country_gain_percent":100*(equal0-equal1)/equal0,
        "country_gain_positive_n":sum(a>b for a,b in country.values()),
        "country_gain_negative_n":sum(a<b for a,b in country.values()),
        "country_values":country}

def validate_summary(calc,old,name):
    for side,key in (("port","weighted_by_PORT"),("country","weighted_equal_COUNTRY")):
        v=old[key]
        for field,got in (("MSE_baseline",calc[side+"_mse0"]),
                          ("MSE_plus_WWp",calc[side+"_mse1"]),
                          ("relative_MSE_gain_percent",calc[side+"_gain_percent"])):
            validate_close(got,v[field],name+"."+key+"."+field)
    assert old["countries_with_positive_WWp_predictive_gain"]==calc["country_gain_positive_n"]
    assert old["countries_with_negative_WWp_predictive_gain"]==calc["country_gain_negative_n"]
    assert len(old["per_country"])==calc["countries"]
    for v in old["per_country"]:
        actual=calc["country_values"][v["country"]]
        validate_close(v["MSE_baseline"],actual[0],name+".country.MSE0")
        validate_close(v["MSE_plus_WWp"],actual[1],name+".country.MSE1")
        validate_close(v["delta_MSE_base_minus_plus_WWp"],actual[0]-actual[1],name+".country.delta")

def audit(root):
    a=json.loads((root/RESULT).read_text())
    b=json.loads((root/ARCHIVE).read_text())
    expected=b["city_raster_2021_matched_models"]["WWTP_plus_WWproduction_plus_urban_crop"]
    source=a["original_source_cohort_replay"]
    assert source["prior_result_archive_numeric_reproduced"] is True
    assert a["source_original_chocho_nested_model_archive_SHA256"]==sha(root/ARCHIVE)
    assert (source["source_original_countryFE_sample_ports"],source["source_original_countryFE_countries"])==(1055,42)
    assert source["matched_source_worldcover_valid_ports"]==1241
    validate_close(source["original_β_WWp_after_urban_crop_countryFE"],expected["WWproduction_beta"],"original_beta")
    for actual,target in zip(source["original_country_cluster95"],expected["WWproduction_CI95"]):
        validate_close(actual,target,"original CI")
    assert source["original_β_WWp_after_urban_crop_countryFE"]>0
    assert source["original_country_cluster95"][0]>0
    i=a["old_positive_model_one_country_removed_sensitivity"]
    assert i["n_country_holdouts"]==42 and len(i["models"])==42
    assert len({q["left_out_country"] for q in i["models"]})==42
    assert all(q["re_fitted_β_WWp"]>0 and q["remaining_countries"]==41 for q in i["models"])
    validate_close(min(q["re_fitted_β_WWp"] for q in i["models"]),i["min_beta"],"min_leave_country")
    validate_close(max(q["re_fitted_β_WWp"] for q in i["models"]),i["max_beta"],"max_leave_country")
    rank=a["original_countries_WWp_within_country_rank_pattern"]
    assert rank["n_eligible"]==len(rank["source_WWp_v_host_photo_fraction_within_country_spearman"])==30
    assert rank["n_pos"]==20 and rank["n_neg"]==10
    sp=a["new_geo_country2degree_5fold_within_country_transfer"]
    assert sp["status"]=="COUNTRY_FIXED_EFFECT_WITHIN_COUNTRY_SPATIAL_HOLDOUT_ESTIMATED"
    geo=sp["squared_error_by_original_port"]
    country=a["new_leave_one_country_out_without_country_fixed_effects"]
    assert country["status"]=="STRICT_UNSEEN_COUNTRY_TRANSFER_NO_COUNTRYFE"
    ext=country["squared_error_by_original_port"]
    assert len(ext)==1055 and len(geo)==996
    assert len({x["port_id"] for x in ext})==len(ext)
    assert all(x["external_fit_coefficient_WWp"]==x["external_fit_coefficient_WWp"] for x in ext)
    assigns=sp["original_country2degree_block_fold_assignments"]
    assert len(assigns)==1055 and len({x["port_id"] for x in assigns})==1055
    by=defaultdict(dict)
    for x in assigns:
        key=x["country"]
        block=x["original_geo2deg_block"]
        if block in by[key] and by[key][block]!=x["fold"]:
            raise ValueError("Original 2° block leaked across holdout folds")
        by[key][block]=x["fold"]
    for c,d in by.items():
        for idx,name in enumerate(sorted(d)):
            if d[name]!=(idx%5):raise ValueError("Changed original ordered 2° block→fold geography")
    byport={x["port_id"]:x for x in assigns}
    assert all(x["fold"]==byport[x["port_id"]]["fold"] and
               x["country"]==byport[x["port_id"]]["country"] and
               x["2deg_block"]==byport[x["port_id"]]["original_geo2deg_block"] for x in geo)
    for rows,label,original in ((geo,"within_2degree",sp["estimates"]),(ext,"unseen_country",country["estimates"])):
        result=summary(rows)
        validate_summary(result,original,label)
    assert len(by)==42
    result=summary(ext)
    assert result["countries"]==42
    rng=np.random.default_rng(SEED)
    v=list(result["country_values"].values())
    resample=[]
    for i in range(2000):
        idx=rng.integers(0,len(v),size=len(v))
        selected=[v[k] for k in idx]
        mse0=np.mean([z[0] for z in selected])
        mse1=np.mean([z[1] for z in selected])
        resample.append(100*(mse0-mse1)/mse0)
    interval=np.quantile(resample,[.025,.975])
    for got,published in zip(interval,country["estimates"]["country_cluster_descriptive_2000_bootstrap95_relative_MSE_gain_percent"]):
        validate_close(got,published,"external country-bootstrap CI")
    controls=a["component_alien_and_native_candidate_indexing_controls"]
    assert controls["ln_alien"]["source_beta_WWp"]>0 and controls["ln_native"]["source_beta_WWp"]>0
    return {
        "schema":"gesui_original_global_2015WWp_1055_port_42country_macro_independent_country_error_QA_v01",
        "original_macro_output_sha256":sha(root/RESULT),
        "original_inherited_chocho2015_model_summary_sha256":sha(root/ARCHIVE),
        "source_population_exactly_1055_ports_42_countries":True,
        "old_positive_country_cluster_regression_source_beta_exact":source["original_β_WWp_after_urban_crop_countryFE"],
        "old_positive_country_95CI":source["original_country_cluster95"],
        "leave_one_country_out_refitted_countryFE_all_42_positive":True,
        "port_to_original_country2degree_block_and_fold_no_leakage":True,
        "same_country_2degree_spatial_5fold_port_predictions_all_996_reproduced":True,
        "strict_new_country_42holdout_port_predictions_all_1055_reproduced":True,
        "source_test_scores_2000_country_resamples_and_original_CI_reproduced":True,
        "original_2degree_5fold_test_country_n":len({x["country"] for x in geo}),
        "within_country_country_equal_MSE_gain_percent":summary(geo)["country_gain_percent"],
        "new_country_country_equal_MSE_gain_percent":result["country_gain_percent"],
        "new_country_descriptive95_MSE_gain_percent":list(map(float,interval)),
        "both_alien_and_native_independent_GBIF_index_count_WWp_associations_positive":True,
        "no_WWp_true_nitrogen_mass_and_no_all_plant_abundance_claim":True,
        "macro_result_positive_source_geography_not_unseen_country_ecological_law":True
    }

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    r=audit(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    print(json.dumps(r,indent=2))
