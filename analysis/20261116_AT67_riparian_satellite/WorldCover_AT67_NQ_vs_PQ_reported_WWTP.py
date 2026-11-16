#!/usr/bin/env python3
"""Posthoc source nutrient-identity falsification:
same Austria67 actual river TP + benthic EQR, observed 2022,
with true 2021 WorldCover landuse:
upstream reported WWTP N/Q vs P/Q and downstream Q-matched placebos.

No measured dissolved/total N at river station, no delivered wastewater
fraction, no causal mediation or terrestrial biomass observations.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,sys,collections
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr,pearsonr

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
sys.path.insert(0,str(HERE.parents[0]/"20261116_austria_wisa_official_source"))
from WorldCover_AT67_WWTP_report_presence_vs_Pdose import original_and_covariates
from AT68_2022_real_WWTP_P_to_ambient_TP_to_benthic_EQR import P_CONVERSION
from AT68_equal_Q_updown_mass_adversarial_qa import model_cv,comparison

SCHEMA="GESUI_20261116_AT67_OFFICIAL_UWWTD2022_REPORTED_N_VS_P_Q_LANDCOVER_CONTROL_V1"
BQE_SOURCE=ROOT/"results/eea2022_biology_bqe_20261108/BENTHIC_INVERTEBRATES_2022_SOURCE_EQR.csv"
BQE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
ORIGINAL67={"15_3":9,"15_4":14,"15_5":21,"16_4":13,"16_5":10}
BASE=["logQ","logArea","lat","lon",
      "wc1km_40_fraction","wc1km_50_fraction",
      "wc1km_10_fraction","wc1km_30_fraction"]
ORIGINAL_EQR_2022_PQ_GAIN=8.388318614326595
ORIGINAL_WATER_2022_PQ_GAIN=4.58218306901158

def crosswalk_source_N_fields(rows):
    data=BQE_SOURCE.read_bytes()
    if hashlib.sha1(("blob "+str(len(data))+"\0").encode()+data).hexdigest()!=BQE_BLOB:
        raise RuntimeError("Source-original EEA2022 genuine water BQE P/N registry blob changed")
    with BQE_SOURCE.open(newline="",encoding="utf-8-sig") as f:
        allrows=list(csv.DictReader(f))
    if len(allrows)!=2475:
        raise RuntimeError("Original published 2022 BQE source row count changed")
    index={q["site_key"]:q for q in allrows if q["countryCode"]=="AT"}
    if len(index)!=68:
        raise RuntimeError("Original Austria 2022 true EEA BQE N/P source site IDs changed")
    for x in rows:
        original=index.get("AT|AT"+x["station"]+"|euMonitoringSiteCode")
        if original is None:raise RuntimeError("Missing official source station "+x["station"])
        qq=math.expm1(x["logQ"])
        if qq<=0:raise RuntimeError("Original mean river Q nonpositive")
        upN=float(original["upstream_2022_EEA_N_measured_t_per_yr"])
        dnN=float(original["downstream_2022_EEA_N_t_per_yr_placebo"])
        if not all(map(math.isfinite,[upN,dnN])) or min(upN,dnN)<0:
            raise RuntimeError("Original 2022 N source annual mass nonfinite/negative")
        x["source_report_N_up_t_yr"]=upN
        x["source_report_N_down_t_yr"]=dnN
        x["source_report_N_UP_Q"]=math.log1p(P_CONVERSION*upN/qq)
        x["source_report_N_DOWN_Q"]=math.log1p(P_CONVERSION*dnN/qq)
        x["source_original_N_year_mass_log"]=math.log1p(upN)
        x["source_original_P_year_mass_log"]=x["logPup"]
        x["source_report_P_DOWN_Q"]=math.log1p(P_CONVERSION*x["Pdown_t_yr"]/qq)
        x["source_original_eligible_N_presence"]=float(upN>0)
    if len(rows)!=67 or collections.Counter(x["block"] for x in rows)!=ORIGINAL67:
        raise RuntimeError("Original genuine 67 source sites changed")
    return rows

def specs(y):
    bg=BASE+([] if y=="logTP" else ["logTP"])
    return {
       "B":bg,
       "P_UP":bg+["logPup_flow_proxy"],
       "N_UP":bg+["source_report_N_UP_Q"],
       "P_N_UP":bg+["logPup_flow_proxy","source_report_N_UP_Q"],
       "P_DOWN":bg+["source_report_P_DOWN_Q"],
       "N_DOWN":bg+["source_report_N_DOWN_Q"],
       "P_UP_RAW":bg+["source_original_P_year_mass_log"],
       "N_UP_RAW":bg+["source_original_N_year_mass_log"],
       "N_PRESENCE":bg+["source_original_eligible_N_presence"],
       "P_UP_AFTER_N_PRESENCE":bg+["source_original_eligible_N_presence","logPup_flow_proxy"]
    }

def eval_source(x,y):
    m=model_cv(x,y,specs(y))
    pairs=[
       ("up_reported_P_over_Q","B","P_UP"),
       ("up_reported_N_over_Q","B","N_UP"),
       ("up_P_Q_added_after_N_Q","N_UP","P_N_UP"),
       ("up_N_Q_added_after_P_Q","P_UP","P_N_UP"),
       ("down_P_Q_placebo","B","P_DOWN"),
       ("down_N_Q_placebo","B","N_DOWN"),
       ("raw_reported_P_year_mass","B","P_UP_RAW"),
       ("raw_reported_N_year_mass","B","N_UP_RAW"),
       ("source_reporting_presence_only","B","N_PRESENCE"),
       ("up_reported_P_Q_after_source_presence","N_PRESENCE","P_UP_AFTER_N_PRESENCE")]
    comp={name:comparison(m,base,extended) for name,base,extended in pairs}
    return m,comp

def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    d,qa=original_and_covariates()
    if d is None or len(d)!=67:
        raise RuntimeError("Original WorldCover true physical source 67/67 failed, no posthoc resampling")
    d=crosswalk_source_N_fields(d)
    upP=np.array([x["logPup_flow_proxy"] for x in d])
    upN=np.array([x["source_report_N_UP_Q"] for x in d])
    upmP=np.array([x["Pup_t_yr"] for x in d])
    upmN=np.array([x["source_report_N_up_t_yr"] for x in d])
    meanP=np.array([x["logTP"] for x in d])
    env={"original_same_eligible_source_N_P_report_presence_overlap_n":
          sum((x["Pup_t_yr"]>0)==(x["source_report_N_up_t_yr"]>0) for x in d),
         "upstream_reported_P_and_N_annual_mass_spearman":float(spearmanr(upmP,upmN).statistic),
         "upstream_reported_P_and_N_potential_Q_diluted_signal_spearman":float(spearmanr(upP,upN).statistic),
         "upstream_reported_P_and_N_potential_Q_diluted_signal_pearson":float(pearsonr(upP,upN).statistic),
         "original_upstream_positive_P_n":int((upmP>0).sum()),
         "original_upstream_positive_N_n":int((upmN>0).sum()),
         "river2022_real_lnTP_vs_sourceN_Q_spearman":float(spearmanr(upN,meanP).statistic)}
    w,wcmp=eval_source(d,"logTP")
    e,ecmp=eval_source(d,"eqr")
    a0=wcmp["up_reported_P_over_Q"]["source_5block_equal_MSE_added_gain_pct"]
    b0=ecmp["up_reported_P_over_Q"]["source_5block_equal_MSE_added_gain_pct"]
    if abs(a0-ORIGINAL_WATER_2022_PQ_GAIN)>1e-8 or abs(b0-ORIGINAL_EQR_2022_PQ_GAIN)>1e-8:
        raise RuntimeError("Independent original Austria67 2022 P/Q source reference not reproduced")
    report={"schema":SCHEMA,"status":"PASS_EXPLORATORY_SOURCE_P_VS_N_COLLINEARITY_ORIGINAL_67_FIELD_SITES",
       "original_source_region5_genuine_observation_n":67,
       "source_true_WWTP2022_eligible_PN_annual_report_N_positive":env["original_upstream_positive_N_n"],
       "source_original_2021_WorldCover_background_valid_stations":qa["source_valid_n"],
       "source_original_EEA2022_benthic_BQE_PN_file_gitblob":BQE_BLOB,
       "original_N_P_reported_same_eligible_plant_pool":True,
       "source_P_N_observational_correlation":env,
       "source_real_2022_river_TP_five_geographic_block_models":w,
       "source_real_2022_river_TP_N_v_P_contrasts":wcmp,
       "source_real_2022_benthic_EQR_five_geographic_block_models":e,
       "source_real_2022_benthic_EQR_N_v_P_contrasts_adjusted_for_real_TP":ecmp,
       "source_original_P_Q_gain_from_landcover_audit_reproduced":b0,
       "same_source_N_P_high_collinearity_prevents_nutrient_specific_causality":True,
       "source_N_t_year_not_actual_river_TN":True,
       "source_Q_based_nutrient_signal_not_measured_effluent_PN_delivered":True,
       "source_no_local_riparian_NDVI_effect_per_prior_smoke_5_of_10_STOP":True,
       "no_aquatic_biomass_riparian_plant_butterfly_or_EcosystemService_effect_fitted":True}
    (a.output/"AUSTRIA67_SOURCE_REPORTED_UPSTREAM_N_VS_P_Q_TO_REAL_WATER_TP_EQR.json").write_text(
       json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_AT67_ORIGINAL_SOURCE_REPORTED_WWTP_N_VS_P_Q="+json.dumps({
      "status":report["status"],
      "n":67,
      "N_P_source_correlation":env,
      "water_logTP_source_gains":{
         k:round(x["source_5block_equal_MSE_added_gain_pct"],4)
         for k,x in wcmp.items()},
      "EQR_conditioned_real_TP_source_gains":{
         k:round(x["source_5block_equal_MSE_added_gain_pct"],4)
         for k,x in ecmp.items()},
      "no_causal_N_or_P_effect":True},ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
