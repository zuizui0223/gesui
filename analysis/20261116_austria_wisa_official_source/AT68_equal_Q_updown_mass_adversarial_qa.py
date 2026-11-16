#!/usr/bin/env python3
"""Post-source adversarial Austrian 2022 wastewater-report P vs river TP/BQE:
identical Q dilution upstream and downstream, leave-top3-P/Q stress, and
augmented-lstsq independent Ridge numeric replay. Still noncausal.
"""
from __future__ import annotations
import argparse,collections,json,math,sys
from pathlib import Path
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import (
    cohort,cv,full_ols,REGIONS,ALPHA,P_CONVERSION
)

SCHEMA="GESUI_AT2022_68_SAME_RIVER_2022_UP_DOWN_2022_ANNUAL_P_EQUAL_FLOW_PLACEBO_V1"
ORIGINAL_SOURCE_RESULTS={
    "logTP_M1_raw_up_gain_pct":5.2916573389747965,
    "logTP_M4_flow_up_gain_pct":16.30205276938712,
    "EQR_M3_observed_stream_TP_gain_pct":33.896408233758514,
    "EQR_M6_flow_up_gain_pct":19.26905930019493
}

def normalized_data():
    d=cohort()
    for r in d:
        q=math.expm1(r["logQ"])
        if q<=0:raise RuntimeError("Original source mean river Q missing")
        r["logPdown_flow_proxy"]=math.log1p(P_CONVERSION*r["Pdown_t_yr"]/q)
    return d

def original_models(d):
    baseline=["logQ","logArea","lat","lon"]
    water={"M0":baseline,
       "M1_raw_up":baseline+["logPup"],
       "M2_up_flow":baseline+["logPup_flow_proxy"],
       "M3_down_flow":baseline+["logPdown_flow_proxy"],
       "M4_up_and_down_flow":baseline+["logPup_flow_proxy","logPdown_flow_proxy"]}
    eco={"M0":baseline,
       "M1_actual_stream_TP":baseline+["logTP"],
       "M2_up_flow":baseline+["logPup_flow_proxy"],
       "M3_down_flow":baseline+["logPdown_flow_proxy"],
       "M4_actual_stream_TP_plus_up_flow":baseline+["logTP","logPup_flow_proxy"],
       "M5_actual_stream_TP_plus_down_flow":baseline+["logTP","logPdown_flow_proxy"],
       "M6_TP_plus_both_flow_P":baseline+["logTP","logPup_flow_proxy","logPdown_flow_proxy"]}
    return baseline,water,eco

def source_predict(train,test,features,outcome):
    # INDEPENDENT from the first model's normal equations:
    # augment observed weighted design by identity sqrt-alpha penalty.
    X=np.array([[x[k] for k in features] for x in train],float)
    Xh=np.array([[x[k] for k in features] for x in test],float)
    mu=X.mean(axis=0);sd=X.std(axis=0,ddof=0)
    sd[sd<1e-12]=1.
    X=(X-mu)/sd;Xh=(Xh-mu)/sd
    y=np.array([x[outcome] for x in train],float)
    mat=np.column_stack([np.ones(len(train)),X])
    check=np.column_stack([np.ones(len(test)),Xh])
    penalty=np.diag([0.]+[math.sqrt(ALPHA)]*len(features))
    joined=np.vstack([mat,penalty])
    target=np.concatenate([y,np.zeros(penalty.shape[0])])
    b=np.linalg.lstsq(joined,target,rcond=None)[0]
    return check@b

def model_cv(d,y,models):
    blocks=sorted(set(x["block"] for x in d))
    if blocks!=sorted(REGIONS):raise RuntimeError("Lost original 5 geographic blocks")
    preds={}
    for lab,fields in models.items():
        squared={}
        for block in blocks:
            tr=[x for x in d if x["block"]!=block]
            te=[x for x in d if x["block"]==block]
            if len(te)<4:raise RuntimeError("Too few original source stations in a region")
            yh=np.array([x[y] for x in te])
            yhat=source_predict(tr,te,fields,y)
            squared[block]=float(np.mean((yh-yhat)**2))
        preds[lab]=squared
    metric={k:float(np.mean(list(v.values()))) for k,v in preds.items()}
    return {"outcome":y,"region_mean_MSE":metric,"MSE_by_region":preds}

def comparison(cvdata,baseline,addition):
    a=cvdata["region_mean_MSE"][baseline]
    b=cvdata["region_mean_MSE"][addition]
    gains=[(cvdata["MSE_by_region"][baseline][k],cvdata["MSE_by_region"][addition][k])
           for k in sorted(REGIONS)]
    first=np.array([p[0] for p in gains])
    second=np.array([p[1] for p in gains])
    if np.any(first<=0):raise RuntimeError("Zero prediction source baseline error")
    # Exhaustive integer resampling of five fixed heldout region error sums:
    # 5^5=3125 pseudo samples, fixed original model forecasts (no refit).
    ix=np.array(np.meshgrid(*([range(5)]*5),indexing="ij")).reshape(5,-1).T
    fixed_sample=100*(first[ix].mean(axis=1)-second[ix].mean(axis=1))/first[ix].mean(axis=1)
    return {
       "source_5block_equal_MSE_added_gain_pct":float(100*(a-b)/a),
       "source_5block_fixed_error_exhaustive_resampling_2p5_97p5_pct":
          [float(x) for x in np.quantile(fixed_sample,[.025,.975])],
       "source_5block_individually_added_gain_pct":{
          k:float(100*(cvdata["MSE_by_region"][baseline][k]-
                       cvdata["MSE_by_region"][addition][k])/
                   cvdata["MSE_by_region"][baseline][k]) for k in sorted(REGIONS)},
       "original_source_region_count":len(REGIONS),
       "source_5region_bootstrap_not_independent_coefficient_CI":True}

def independent_check(d):
    baseline=["logQ","logArea","lat","lon"]
    water={
       "M0":baseline,
       "M1_raw_up":baseline+["logPup"],
       "M2_up_flow":baseline+["logPup_flow_proxy"]}
    eco={
       "M0":baseline,
       "M1_TP":baseline+["logTP"],
       "M2_up_flow":baseline+["logPup_flow_proxy"]}
    wa=model_cv(d,"logTP",water)
    er=model_cv(d,"eqr",eco)
    refs=[
       ("logTP_M1_raw_up_gain_pct",comparison(wa,"M0","M1_raw_up")),
       ("logTP_M4_flow_up_gain_pct",comparison(wa,"M0","M2_up_flow")),
       ("EQR_M3_observed_stream_TP_gain_pct",comparison(er,"M0","M1_TP")),
       ("EQR_M6_flow_up_gain_pct",comparison(er,"M0","M2_up_flow"))]
    res=[]
    for key,val in refs:
        ref=ORIGINAL_SOURCE_RESULTS[key]
        actual=val["source_5block_equal_MSE_added_gain_pct"]
        err=abs(actual-ref)
        if err>1e-8:
            raise RuntimeError("Independent source Rigde solver FAILED original result reproduction "+
                    str((key,err,actual,ref)))
        res.append({"metric":key,"original_normal_equation_gain":ref,
                    "independent_augmented_lstsq_gain":actual,
                    "absolute_percentage_point_error":err,
                    "independent_numeric_identity":True})
    return res

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    d=normalized_data()
    if len(d)!=68 or sum(q["Pup_t_yr"]>0 for q in d)!=53:
        raise RuntimeError("Original signed Austria68 source cohort changed")
    qa=independent_check(d)
    baseline,wat_model,eco_model=original_models(d)
    water=model_cv(d,"logTP",wat_model)
    eqr=model_cv(d,"eqr",eco_model)
    water_comp={
      "UPSTREAM_normalized_source_P_added_vs_B":comparison(water,"M0","M2_up_flow"),
      "DOWNSTREAM_normalized_source_P_added_vs_B":comparison(water,"M0","M3_down_flow"),
      "UPSTREAM_normalized_source_P_added_after_DOWNSTREAM":
         comparison(water,"M3_down_flow","M4_up_and_down_flow"),
      "DOWNSTREAM_normalized_source_P_added_after_UPSTREAM":
         comparison(water,"M2_up_flow","M4_up_and_down_flow")}
    eq_comp={
      "UPSTREAM_normalized_source_P_added_vs_B":comparison(eqr,"M0","M2_up_flow"),
      "DOWNSTREAM_normalized_source_P_added_vs_B":comparison(eqr,"M0","M3_down_flow"),
      "actual_2022_stream_TP_added_vs_B":comparison(eqr,"M0","M1_actual_stream_TP"),
      "UPSTREAM_normalized_source_P_after_actual_stream_TP":
         comparison(eqr,"M1_actual_stream_TP","M4_actual_stream_TP_plus_up_flow"),
      "DOWNSTREAM_normalized_source_P_after_actual_stream_TP":
         comparison(eqr,"M1_actual_stream_TP","M5_actual_stream_TP_plus_down_flow"),
      "UPSTREAM_normalized_source_P_AFTER_TP_and_DOWNSTREAM":
         comparison(eqr,"M5_actual_stream_TP_plus_down_flow","M6_TP_plus_both_flow_P"),
      "DOWNSTREAM_normalized_source_P_AFTER_TP_and_UPSTREAM":
         comparison(eqr,"M4_actual_stream_TP_plus_up_flow","M6_TP_plus_both_flow_P")}
    sorted_rows=sorted(d,key=lambda x:(x["logPup_flow_proxy"],x["station"]),reverse=True)
    excluded=[r["station"] for r in sorted_rows[:3]]
    trim=[r for r in d if r["station"] not in excluded]
    if len(trim)!=65 or set(x["block"] for x in trim)!=set(REGIONS):
        raise RuntimeError("Source fixed top3 WWTP P/Q no longer preserves 65 Austrian stations and five blocks")
    trim_water=model_cv(trim,"logTP",wat_model)
    trim_eqr=model_cv(trim,"eqr",eco_model)
    source={
      "schema":SCHEMA,
      "status":"PASS_INDEPENDENT_QA_FIVE_REGION_FLOW_MATCHED_NEGATIVE_CONTROL_DIAGNOSTIC",
      "original_source_68_benthic_actual_2022_river_TP_sites":len(d),
      "original_genuine_5_regions":REGIONS,
      "original_2022_upstream_eligible_reported_P_positive":sum(q["Pup_t_yr"]>0 for q in d),
      "original_2022_downstream_eligible_reported_P_positive":sum(q["Pdown_t_yr"]>0 for q in d),
      "source_mass_t_per_year_to_max_unattenuated_concentration_mgP_L_factor":P_CONVERSION,
      "same_flow_dilution_upstream_and_downstream_and_covariate_controls":True,
      "source_numeric_independent_solver_checks":qa,
      "observed_68_water_TP_original_mass_controls":water,
      "observed_68_benthic_EQR_original_mass_controls":eqr,
      "reported_2022_water_TP_source_dilution_placebo_contrasts":water_comp,
      "reported_2022_benthic_EQR_source_dilution_placebo_contrasts":eq_comp,
      "posthoc_exclusion_of_3_largest_original_upstream_P_per_flow_site_ids":excluded,
      "source_trimmed_65_same_5blocks_water_TP":{
        "UP_normalized":comparison(trim_water,"M0","M2_up_flow"),
        "DOWN_normalized":comparison(trim_water,"M0","M3_down_flow")},
      "source_trimmed_65_same_5blocks_EQR":{
        "UP_normalized":comparison(trim_eqr,"M0","M2_up_flow"),
        "DOWN_normalized":comparison(trim_eqr,"M0","M3_down_flow"),
        "UP_after_measured_TP":comparison(trim_eqr,"M1_actual_stream_TP","M4_actual_stream_TP_plus_up_flow")},
      "geographic_5block_fixed_loss_resampling_is_descriptive_not_confirmatory":True,
      "posthoc_numeric_test_not_independent_field_replication":True,
      "source_annual_TP_from_anyWWTP_not_unique_source_delivery":True,
      "no_2022_rootzone_N_or_riparian_insect_actual_consumer_response":True,
      "no_wwtp_causal_transfer_ecosystem_service_or_butterfly_effect":True}
    (a.output/"AUSTRIA68_2022_UPDOWN_P_FLOW_IDENTICAL_PLACEBO_AND_SOURCE_LEVERAGE.json").write_text(
       json.dumps(source,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    succinct={
      "status":source["status"],
      "original_n":len(d),
      "source_independent_numeric_QA":len(qa),
      "water_TP":{
         k:{"gain":v["source_5block_equal_MSE_added_gain_pct"],
            "interval":v["source_5block_fixed_error_exhaustive_resampling_2p5_97p5_pct"],
            "gains_per_region":v["source_5block_individually_added_gain_pct"]}
         for k,v in water_comp.items()},
      "benthic_EQR":{
         k:{"gain":v["source_5block_equal_MSE_added_gain_pct"],
            "interval":v["source_5block_fixed_error_exhaustive_resampling_2p5_97p5_pct"],
            "regions_positive":sum(x>0 for x in v["source_5block_individually_added_gain_pct"].values())}
         for k,v in eq_comp.items()},
      "trim_n":len(trim),
      "trim_water_UPgain":source["source_trimmed_65_same_5blocks_water_TP"]["UP_normalized"]["source_5block_equal_MSE_added_gain_pct"],
      "trim_water_DOWNgain":source["source_trimmed_65_same_5blocks_water_TP"]["DOWN_normalized"]["source_5block_equal_MSE_added_gain_pct"],
      "trim_EQR_UPgain":source["source_trimmed_65_same_5blocks_EQR"]["UP_normalized"]["source_5block_equal_MSE_added_gain_pct"],
      "trim_EQR_DOWNgain":source["source_trimmed_65_same_5blocks_EQR"]["DOWN_normalized"]["source_5block_equal_MSE_added_gain_pct"],
      "no_sewage_causal_transfer":True}
    print("GESUI_AT68_TRUE_UPDOWN_MATCHED_Q_DILUTION_NEGATIVE_CONTROL="+json.dumps(
         succinct,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
