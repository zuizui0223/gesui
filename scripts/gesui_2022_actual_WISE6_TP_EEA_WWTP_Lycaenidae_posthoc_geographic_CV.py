#!/usr/bin/env python3
"""EEA2026 original observed 2022 in-river TP vs EEA WWTP reported source P and
iNaturalist public Lycaenidae photo composition (POSTHOC observational geography).

All original 2022 TP numeric site medians were aggregated and Git SHA-frozen
after source-only repeated-site gate before entering this model. All original
2021–2023 butterfly photo composition and WorldCover2021 land fractions are
source-locked. This is not a mechanistic pathway or mediation model: EEA WWTP
outfall location has not been hydrologically connected to a specific
water-quality monitoring site within each 1degree photo reporting cell.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

ORIGINAL_WISE_ZIP_SHA="2a2ef2a43a16ec068b7974a80ae5fca49e44697660a46ce3b91fb1813abe62d2"
MIN_CELLS=20
MIN_COUNTRY={"DE":10,"AT":6}
MIN_BLOCKS=10
SOURCE_INPUTS={
 "WWTP_P_reported_t_per_year":"registered_P_t_year",
 "WWTP_N_reported_t_per_year":"registered_N_t_year",
 "WWTP_qualified_source_count":"registered_2022_EEA_both_measured_FW_outfalls",
}
EARLY_LYCAENID=["2021_Lycaenidae_photo_n","2021_three_family_photo_total"]
RNG_SEED=20261109
BOOTSTRAP_N=999

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def key(x):return (x["countryCode"],int(float(x["lat_floor"])),int(float(x["lon_floor"])))
def group(x):return x["countryCode"]+"|"+str(x["country_2deg_block"])
def target(x,mode):
    if mode=="actual_TP_2022":
        # Original WISE6 source unit mg{P}/L and original median of daily water samples.
        return math.log1p(float(x["actual_2022_riverTP_cell_median_mgP_per_L"]))
    if mode in ("photo_2022","photo_2023"):
        yy=2022 if mode=="photo_2022" else 2023
        n=float(x[f"{yy}_three_family_photo_total"])
        l=float(x[f"{yy}_Lycaenidae_photo_n"])
        return math.log((l+0.5)/(n-l+0.5))
    raise ValueError("Unknown outcome "+mode)

def arr(rows,fn):return np.asarray([fn(x) for x in rows],float)

def design(train,test,outcome,add_feature):
    if add_feature is None:extra=[]
    else:extra=([add_feature] if isinstance(add_feature,str) else list(add_feature))
    base=[
        ("country_is_AT",lambda r:float(r["countryCode"]=="AT")),
        ("latitude",lambda r:float(r["lat_floor"])+.5),
        ("longitude",lambda r:float(r["lon_floor"])+.5),
        ("2021_WorldCover_built_fraction",lambda r:float(r["WC2021_512_built_land_fraction"])),
        ("2021_WorldCover_crop_fraction",lambda r:float(r["WC2021_512_crop_land_fraction"])),
        ("2021_WorldCover_country_terrestrial_km2",lambda r:math.log1p(float(r["WC2021_512_country_terrestrial_area_km2"]))),
    ]
    if outcome!="actual_TP_2022":
        base.extend([
           ("2021_Plantae_photo_effort",lambda r:math.log1p(float(r["2021_Plantae_photo_n"]))),
           ("2021_Lycaenidae_photo_logit",lambda r:math.log((float(r["2021_Lycaenidae_photo_n"])+.5)/
                            (float(r["2021_three_family_photo_total"])-float(r["2021_Lycaenidae_photo_n"])+.5))),
           ("2021_target_butterfly_photo_effort",lambda r:math.log1p(float(r["2021_three_family_photo_total"]))),
        ])
    X=[np.ones(len(train))];T=[np.ones(len(test))]
    for name,fn in base:
        aa=arr(train,fn);bb=arr(test,fn)
        if name=="country_is_AT":
            X.append(aa);T.append(bb);continue
        m=float(aa.mean());sd=float(aa.std())
        if sd<1e-13:raise ValueError("Degenerate training baseline "+name)
        X.append((aa-m)/sd);T.append((bb-m)/sd)
    baseline=np.column_stack(X);testing=np.column_stack(T)
    for feature in extra:
        if feature in SOURCE_INPUTS:
            source=SOURCE_INPUTS[feature]
            fn=lambda r:math.log1p(float(r[source]))
        elif feature=="observed_2022_river_TP":
            fn=lambda r:math.log1p(float(r["actual_2022_riverTP_cell_median_mgP_per_L"]))
        else:raise ValueError("Unknown source feature "+feature)
        aa=arr(train,fn);bb=arr(test,fn)
        ztrain=np.zeros(len(train));ztest=np.zeros(len(test))
        for country in ("DE","AT"):
            mask=np.asarray([x["countryCode"]==country for x in train])
            mt=np.asarray([x["countryCode"]==country for x in test])
            if mask.sum()<3:raise ValueError("Insufficient country-fold training records")
            mu=float(aa[mask].mean());std=float(aa[mask].std())
            if std<1e-13:raise ValueError("Source within-country training variance zero "+feature+" "+country)
            ztrain[mask]=(aa[mask]-mu)/std
            ztest[mt]=(bb[mt]-mu)/std
        X.append(ztrain);T.append(ztest)
    extended=np.column_stack(X);extend_test=np.column_stack(T)
    return baseline,testing,extended,extend_test

def ols(X,T,y):
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Rank-deficient original 2degree block CV training design")
    return T@np.linalg.lstsq(X,y,rcond=None)[0]

def heldout(rows,outcome,source,with_preds=False):
    regions=sorted(set(group(r) for r in rows))
    pred=[]
    for region in regions:
        train=[x for x in rows if group(x)!=region]
        test=[x for x in rows if group(x)==region]
        X0,T0,X1,T1=design(train,test,outcome,source)
        yt=arr(train,lambda r:target(r,outcome))
        yr=arr(test,lambda r:target(r,outcome))
        p0=ols(X0,T0,yt);p1=ols(X1,T1,yt)
        for x,y,a,b in zip(test,yr,p0,p1):
            pred.append({
                "country":x["countryCode"],"original_2degree_block":region,
                "original_1degree_source_latfloor":x["lat_floor"],
                "original_1degree_source_lonfloor":x["lon_floor"],
                "original_2022_TP_valid_site_years":x["actual_2022_riverTP_qualified_stations_n"],
                "2022_measured_river_TP_mgP_per_L":x["actual_2022_riverTP_cell_median_mgP_per_L"],
                "observed_target_value":float(y),
                "baseline_pred":float(a),
                "source_pred":float(b),
                "baseline_SE":float((y-a)**2),
                "added_source_SE":float((y-b)**2),
                "baseline_abs_err":float(abs(y-a)),
                "source_abs_err":float(abs(y-b)),
            })
    baseline=float(np.mean([x["baseline_SE"] for x in pred]))
    sourceMSE=float(np.mean([x["added_source_SE"] for x in pred]))
    ans={
       "target":outcome,"added_source_predictor":source,
       "n_independent_country_1degree_photo_cells":len(rows),
       "n_original_country_2degree_holdout_regions":len(regions),
       "baseline_MSE":baseline,
       "source_MSE":sourceMSE,
       "MSE_improvement_percent":100*(baseline-sourceMSE)/baseline,
       "baseline_MAE":float(np.mean([x["baseline_abs_err"] for x in pred])),
       "source_MAE":float(np.mean([x["source_abs_err"] for x in pred])),
       "n_hyperparameters_tuned_on_test_data":0,
       "country_results":{},
       "not_direct_2022_nutrient_chemical_mediation":True,
    }
    for c in ("DE","AT"):
        pp=[x for x in pred if x["country"]==c]
        if not pp:continue
        b=float(np.mean([x["baseline_SE"] for x in pp]))
        a=float(np.mean([x["added_source_SE"] for x in pp]))
        ans["country_results"][c]={"n":len(pp),"baseline_MSE":b,
                                   "source_MSE":a,"MSE_improvement_percent":100*(b-a)/b}
    if with_preds:
        ans["original_country2degree_block_heldout_cell_predictions"]=pred
        rng=np.random.default_rng(RNG_SEED)
        by=defaultdict(list)
        for x in pred:by[x["original_2degree_block"]].append(x)
        groups=list(by)
        sims=[]
        for _ in range(BOOTSTRAP_N):
            sample=rng.choice(groups,size=len(groups),replace=True)
            obs=[x for key in sample for x in by[key]]
            b=float(np.mean([x["baseline_SE"] for x in obs]))
            v=float(np.mean([x["added_source_SE"] for x in obs]))
            if b>0:sims.append(100*(b-v)/b)
        ans["source_gain_block_error_bootstrap_95pct_interval"]={
          "n_draws":len(sims),"seed":RNG_SEED,
          "descriptive_interval":[float(v) for v in np.quantile(sims,[.025,.975])],
          "not_randomized_causal_inference":True}
    return ans

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-WISE6-numeric-TP",type=Path,required=True)
    p.add_argument("--original-WorldCover-2022-photo",type=Path,required=True)
    p.add_argument("--original-2021-photo",type=Path,required=True)
    p.add_argument("--original-2023-photo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    w=json.loads(a.original_WISE6_numeric_TP.read_text())
    wc=json.loads(a.original_WorldCover_2022_photo.read_text())
    h=json.loads(a.original_2021_photo.read_text())
    p23=json.loads(a.original_2023_photo.read_text())
    if not w["numeric_P_water_TP_to_Lycaenidae_photo_geographic_gate_pass"]:
        raise ValueError("Nonnumeric source passed but original quantified numeric TP source gate failed: do not fit")
    if w["original_2022_EEA2026_WISE6_DisaggregatedData_ZIP_verified"]["original_archive_sha256"]!=ORIGINAL_WISE_ZIP_SHA:
        raise ValueError("New real water chemistry original file SHA changed")
    if w["original_2022_TP_analyte_original_code"]!="CAS_7723-14-0" or w["original_2022_TP_unit_original"]!="mg{P}/L":
        raise ValueError("TP analyte/unit not unambiguous")
    if wc["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED":
        raise ValueError("Original frozen 2022 landcover+photo source not final")
    if w["original_2022_WorldCover43_SHA256"]!=sha(a.original_WorldCover_2022_photo) or p23["original_2022_WorldCover_land_photo_source_sha256"]!=sha(a.original_WorldCover_2022_photo):
        raise ValueError("2021/2022/2023 photo records not source-locked to exact country/1degree grid")
    if h["original_2021_WorldCover2022_photo_source_full_SHA256"]!=sha(a.original_WorldCover_2022_photo):
        raise ValueError("2021 photo history not matched to original source grid")
    source={key(x):x for x in wc["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]}
    prior={key(x):x for x in h["source_2021_same43_cell_historical_butterfly_photo_records"]}
    future={key(x):x for x in p23["source_2023_all_original43_butterfly_photo_count_and_receipts"]}
    rows=[]
    for q in w["strict_original_42_Lycaenidae_country_photo_cells_and_measured_TP_outcomes"]:
        if not q["2022_TP_quantified_ge2site_cell_gate_pass"]:continue
        r=source.get(key(q))
        if r is None:raise ValueError("Frozen new 2022 river TP source cell not original butterfly geography")
        k=key(q)
        z=prior[k];f=future[k]
        if not z["2021_photo_support_min15"] or not f["2023_photo_support_ge20"]:
            raise ValueError("Original TP cell not in valid 2021/2023 photo history strict support")
        row=dict(r)
        row["2021_Lycaenidae_photo_n"]=z["2021_Lycaenidae_photo_n"]
        row["2021_three_family_photo_total"]=z["2021_three_family_photo_total"]
        row["2022_three_family_photo_total"]=r["2022_total_3family_photos"]
        row["2023_Lycaenidae_photo_n"]=f["2023_Lycaenidae_photo_n"]
        row["2023_three_family_photo_total"]=f["2023_three_family_photo_total"]
        row["actual_2022_riverTP_qualified_stations_n"]=q["2022_WISE6_quantified_TP_site_years_n"]
        row["actual_2022_riverTP_cell_median_mgP_per_L"]=q["2022_WISE6_measured_site_median_TP_cell_median_mg_P_per_L"]
        if not isinstance(row["actual_2022_riverTP_cell_median_mgP_per_L"],(float,int)):
            raise ValueError("Numerically invalid preselected TP region")
        rows.append(row)
    country=Counter(r["countryCode"] for r in rows)
    blocks=set(group(r) for r in rows)
    if len(rows)<MIN_CELLS or country["DE"]<MIN_COUNTRY["DE"] or country["AT"]<MIN_COUNTRY["AT"] or len(blocks)<MIN_BLOCKS:
        raise ValueError("Original numerical water sample quantification did not meet regional support gate")
    env={}
    for name in SOURCE_INPUTS:
        env[name]=heldout(rows,"actual_TP_2022",name,with_preds=True)
        print("GESUI_2022_REAL_WISE6_RIVER_TP_VS_2022_EEA_WWTP_SOURCE_SPATIAL_CV",
              json.dumps({"source":name,"cells":len(rows),
                "MSE_base":env[name]["baseline_MSE"],"MSE_source":env[name]["source_MSE"],
                "gain_percent":env[name]["MSE_improvement_percent"],
                "countries":env[name]["country_results"]},sort_keys=True),flush=True)
    butterfly={}
    for targetyear in ("photo_2022","photo_2023"):
        butterfly[targetyear]={
          "river_TP_above_baseline":heldout(rows,targetyear,"observed_2022_river_TP",with_preds=True),
          "reported_WWTP_source_P_above_baseline":heldout(rows,targetyear,"WWTP_P_reported_t_per_year",with_preds=True),
          "river_TP_above_reported_WWTP_P_and_baseline":heldout(rows,targetyear,
                 ("WWTP_P_reported_t_per_year","observed_2022_river_TP"),with_preds=True),
        }
        for name,result in butterfly[targetyear].items():
            print("GESUI_2022_REAL_TP_TO_LYCAENIDAE_GEOGRAPHIC_PHOTOGRAPHIC_REPORTING_CV",
                  json.dumps({"target":targetyear,"model":name,"cells":len(rows),
                  "baseline_MSE":result["baseline_MSE"],"source_MSE":result["source_MSE"],
                  "gain_percent":result["MSE_improvement_percent"],
                  "countries":result["country_results"]},sort_keys=True),flush=True)
    correlation={}
    for c in ("DE","AT"):
        subset=[r for r in rows if r["countryCode"]==c]
        sourceP=[math.log1p(float(r["registered_P_t_year"])) for r in subset]
        riverP=[math.log1p(float(r["actual_2022_riverTP_cell_median_mgP_per_L"])) for r in subset]
        photoL=[target(r,"photo_2022") for r in subset]
        envR,_=spearmanr(sourceP,riverP)
        rpr,_=spearmanr(riverP,photoL)
        sourcePhoto,_=spearmanr(sourceP,photoL)
        correlation[c]={
            "n_country_1degree_cells":len(subset),
            "source_reported_P_mass_to_real_2022_river_TP_Spearman":float(envR),
            "riverTP_2022_to_Lycaenidae_public_photo_logit_Spearman":float(rpr),
            "source_P_to_Lycaenidae_public_photo_logit_Spearman":float(sourcePhoto),
            "no_mediation_or_causality_identified":True}
    measured=[float(r["actual_2022_riverTP_cell_median_mgP_per_L"]) for r in rows]
    report={
      "schema":"gesui_2022_original_WISE6_actual_riverTP_and_reported_EEA_WWTP_P_vs_2021_2023_Lycaenidae_photo_geographic_v01",
      "UTC":datetime.now(timezone.utc).isoformat(),
      "POSTHOC_comparison_after_positive_source_photo_geography":True,
      "original_EEA2026_numeric_TP_site_and_country_cell_SHA256":sha(a.original_WISE6_numeric_TP),
      "original_ESAWC2021_2022_family_photo_WorldCover_SHA256":sha(a.original_WorldCover_2022_photo),
      "original_2021_family_photo_history_source_SHA256":sha(a.original_2021_photo),
      "original_2023_photo_source_SHA256":sha(a.original_2023_photo),
      "original_country_1degree_source_cells_for_real_TP_numeric_QC":len(rows),
      "country_1degree_cells_by_country":dict(country),
      "original_country_2degree_block_count":len(blocks),
      "original_country_1deg_TP_region_median_mg_P_per_L_range":[min(measured),max(measured)],
      "original_country_1deg_TP_region_median_mg_P_per_L_overall_median":float(np.median(measured)),
      "real_2022_TP_vs_official_EEA2022_reported_WWTP_source_country2degree_spatial_heldout":env,
      "real_2022_TP_and_2022_reported_WWTP_P_vs_Lycaenidae_public_photographic_report_country2degree_spatial_heldout":butterfly,
      "country_nonparametric_cell_source_TP_Lycaenidae_photo_recording_correlations":correlation,
      "source_annual_reported_P_N_tonnes_vs_actual_riverTP_mg_P_per_L_are_different_physical_quantities":True,
      "TP_cell_median_not_direct_chemical_dose_at_any_lycaenid_mud_puddling_location":True,
      "source_EEA_WWTP_reporting_outfalls_not_verified_upstream_of_any_of_the_WISE6_monitoring_sites":True,
      "no_formal_statistical_mediation_or_sewage_ecological_causal_claim":True,
      "iNaturalist_photo_records_not_actual_butterfly_individual_counts":True,
      "source_selected_regions_prior_to_first_TP_numeric_outcome":True,
      "no_new_WISE6_numeric_TP_or_photo_API_download_in_this_model":True,
      "limitations":[
       "Source outfall mass is a cell-level public geographic total, NOT hydrologic upstream load delivered to an individual river monitoring station",
       "Both sampled river stations and photographed butterflies are nonrandom observations in different ecological media",
       "The 2021 pre-year butterfly photographic composition may encode ecosystem and human observer persistence",
       "WorldCover built cropland land area and country longitude latitude cannot fully control urban observation effort or source selection",
       "Only 1 degree cells, not fish/butterfly sites or WWTP-registered chemically exposed river reaches, are modeled",
       "EEA2026 WISE6 2022 individual observations with belowLOQ were excluded with no substitution; missing site chemistry not zero",
       "P and N WWTP reported annual masses collinear; no independent biochemical P mechanism tested"
      ]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA2026_REAL_2022_RIVER_TP_EEAWWTP_AND_LYCAENIDAE_PHOTO_SOURCE_LINK_GATED_DONE",
          json.dumps({"cells":len(rows),"country":dict(country),
           "TP_sourceP_MSE_gain_percent":env["WWTP_P_reported_t_per_year"]["MSE_improvement_percent"],
           "Lycaenidae_2022_photoTP_MSE_gain_percent":butterfly["photo_2022"]["river_TP_above_baseline"]["MSE_improvement_percent"],
           "Lycaenidae_2022_photoTP_incremental_over_Pmass_MSE_gain_percent":
               butterfly["photo_2022"]["river_TP_above_reported_WWTP_P_and_baseline"]["MSE_improvement_percent"],
           "no_sewage_mechanism":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
