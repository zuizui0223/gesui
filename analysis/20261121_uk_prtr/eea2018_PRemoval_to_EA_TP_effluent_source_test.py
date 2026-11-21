#!/usr/bin/env python3
"""Exploratory 2018 UK UWWTD P-removal flag vs observed EA final effluent TP.

Source-natively reported 2018 WWTP PRemoval boolean, original 2018 UK
freshwater outfall point, current EA WQE Water Company final-treated SA
sample point unique reciprocal within 250m, and actual 2018 field/lab
TP0348 mg/L, aggregated one 2018 annual median per unique native SA.
Boolean flags are documented as water treatment technology in the Eionet
UWWTD Article15 data dictionary (uwwPRemoval / uwwNRemoval).

This is a GEOGRAPHICALLY matched, NONRANDOMIZED source descriptive
comparison. UWWTD facility identity has NOT been independently joined by
licence/native ID to EA SA samplepoint; do NOT call this proof of phosphorus
removal effectiveness, delivered river P loads or any BWMP causal response.

Store only aggregate group counts, observed medians/ratios, site-bootstrap
descriptive CI and fixed covariate adjusted conditional coefficient.
Official original plant/site/TP raw source rows, native IDs and coordinates
are NEVER published to git or in job logs.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,sys,tempfile,urllib.request,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261121_uk_prtr"))
import eea2018_physical_outlet_vs_EA_SA_source_spatial as eea_geo
import eea_UWWTD2018_UK_outlet_source_gate as eea
import eea2018_outlet_EA_TP_2018_ChemPop_river_gate as geo
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as ea
import ea_WQE_national_SA_actual2017_18_TP_chemistry_bio_screen as wq
import ea_wqe_2017_2018_SA_tp_siteyear_census as original_wqe

REPORT_YEAR=2018
N_EXPECTED_UWWTP=1799
N_EXPECTED_MUTUAL_SA=930
MIN_GROUP_N=30
MIN_REGIONS=3
MIN_PER_REGION=10
MIN_CAPACITY_VALID_FRACTION=.70
N_BOOT=2000
RNG_SEED=20261121
MAX_RADIUS_M=250

def plant_records():
    req=urllib.request.Request(eea.URL,headers={
        "User-Agent":"gesui-source-verified-EEA2018-UK-phosphorus-removal-vs-final-TP"})
    with tempfile.TemporaryDirectory(prefix="gesui_uk2018_PRemoval_source_") as folder:
        p=Path(folder)/"original.zip"
        n=0
        with urllib.request.urlopen(req,timeout=100) as r,p.open("wb") as file:
            if r.status!=200:raise RuntimeError("Original EEA UK UWWTD archive not accessible")
            sha=hashlib.sha256()
            while True:
                b=r.read(131072)
                if not b:break
                n+=len(b)
                if n>eea.MAX_BYTES:raise RuntimeError("Original official UWWTD archive exceeds cap")
                file.write(b);sha.update(b)
        with zipfile.ZipFile(p) as z:
            original=z.read(eea.BASE+"UWWTPS.csv")
            if hashlib.sha256(original).hexdigest()!=eea.EXPECT["UWWTPS.csv"]:
                raise RuntimeError("Original UK UWWTPS scientific source fingerprint drift")
            data=list(csv.DictReader(io.StringIO(original.decode("utf-8-sig","replace"))))
    selected=[x for x in data if x["repCode"].strip()=="UK5"]
    if len(selected)!=N_EXPECTED_UWWTP:raise RuntimeError("Original 2018 UK plants no longer 1799")
    plants={}
    origin=collections.Counter()
    for r in selected:
        key=r["uwwCode"].strip()
        if not key or key in plants:raise RuntimeError("Original 2018 UK plant IDs blank or duplicate")
        pflag=str(r.get("uwwPRemoval","")).strip()
        nflag=str(r.get("uwwNRemoval","")).strip()
        if pflag not in ("0","1") or nflag not in ("0","1"):
            origin["original_2018_UK_unknown_P_or_N_treatment_flag"]+=1
            continue
        # eea.finite in the OUTLET source audit is correctly restricted
        # to coordinate magnitudes <10000: never reuse it for capacity (PE).
        try:
            cap=float(str(r.get("uwwCapacity","")).strip())
            if not math.isfinite(cap) or cap<=0:cap=None
        except (ValueError,TypeError):
            cap=None
        if cap is None:origin["original_2018_UK_invalid_original_capacity"]+=1
        region=str(r.get("uwwNUTS","")).strip()
        original_region=region[:3] if len(region)>=3 and region.startswith("UK") else None
        if not original_region:origin["original_2018_UK_unknown_region"]+=1
        plants[key]={"PRemoval":int(pflag),"NRemoval":int(nflag),
                     "capacity":float(cap) if cap is not None and cap>0 else None,
                     "original_region":original_region,
                     "source_report_year":REPORT_YEAR}
        origin["2018_P_removal_"+pflag]+=1
        origin["2018_N_removal_"+nflag]+=1
    return plants,{
        "original_2018_UK_UWWTP_source_n":len(selected),
        "source_2018_valid_treatment_rows":len(plants),
        "original_WWTW_capacity_and_region_source_quality":dict(origin),
        "EEA_original_UWWTPS_SHA256":eea.EXPECT["UWWTPS.csv"],
        "official_Eionet_type":"uwwPRemoval and uwwNRemoval boolean source treatment fields",
        "not_measured_NP_load_or_randomized_treatment":True
    }

def original_TP_dates(names):
    """Original WQE page and source code verification before outcome aggregation."""
    records=collections.defaultdict(lambda:collections.defaultdict(list))
    stats=collections.Counter()
    batches=0
    pagecount=0
    for lo in range(0,len(names),wq.BATCH_N):
        codes=names[lo:lo+wq.BATCH_N]
        batches+=1
        for index in range(wq.MAX_BATCH_PAGES):
            rows,rec=wq.page(REPORT_YEAR,codes,index*wq.PAGE_LIMIT)
            pagecount+=1
            for r in rows:
                stats["returned_source_rows"]+=1
                if r["determinand.notation"].strip()!="0348":
                    raise RuntimeError("WQE returned wrong original TP determinand")
                if not r["phenomenonTime"].strip().startswith("2018"):
                    raise RuntimeError("WQE actual source row has wrong original reported sampling year")
                if r["samplingPoint.samplingPointType"].strip()!=original_wqe.SA:
                    stats["non_WaterCompany_final_type"]+=1;continue
                if r["sampleMaterialType"].strip()!="FINAL SEWAGE EFFLUENT":
                    stats["non_final_treated_sample_material"]+=1;continue
                if r["unit"].strip()!="MILLIGRAM PER LITRE":
                    stats["non_mg_per_L"]+=1;continue
                try:tp=float(r["result"].strip())
                except (ValueError,TypeError):stats["censored_or_nonnumeric"]+=1;continue
                if not math.isfinite(tp) or tp<0:
                    stats["negative_or_nonfinite_original_TP"]+=1;continue
                sa_code=r["samplingPoint.notation"].strip()
                day=r["phenomenonTime"].strip()[:10]
                records[sa_code][day].append(tp)
                stats["valid_numeric_TP_records"]+=1
            if len(rows)<wq.PAGE_LIMIT:break
            if index==wq.MAX_BATCH_PAGES-1:
                raise RuntimeError("Source EA WQE 250-per-page capped TP records truncated")
    siteyear={}
    for code,days in records.items():
        daily=[float(np.median(v)) for v in days.values()]
        if not daily:continue
        siteyear[code]={"annual_median_TP_mgL":float(np.median(daily)),
                        "n_source_dates":len(days),
                        "n_raw_numeric_TP":sum(len(v) for v in days.values()),
                        "n_source_months":len({d[:7] for d in days}),
                        "this_is_final_effluent_TP_mgL_not_kg_yr":True}
    stats.update({"distinct_numeric_TP_SA_stations":len(siteyear),
                  "official_WQE_250row_page_n":pagecount,
                  "official_source_batches":batches,
                  "TP_per_SA_year_is_median_of_source_day_medians":True})
    return siteyear,dict(stats)

def original_sampling_selection_QA(matches,plants,tps):
    """Audit P-removal measurement-selection BEFORE interpreting TP outcomes."""
    prior=collections.Counter()
    actual=collections.Counter()
    for code,item in matches.items():
        plant=plants.get(item["outlet_id"][1])
        if plant is None:
            prior["original_plant_ID_missing"]+=1
            continue
        k=str(plant["PRemoval"])
        prior[k]+=1
        if code in tps:
            actual[k]+=1
    return {
      "official_2018_EEA_plant_PRemoval_yes_of_all_1799":418,
      "original_all_2018_UWWTPS_plant_rows":1799,
      "geographical_reciprocal250m_source_before_WQE_TP_by_PRemoval":dict(prior),
      "actually_2018_TP_measured_source_by_PRemoval":dict(actual),
      "sampling_fraction_with_actual_TP_by_original_PRemoval":{
        k:float(actual[k]/prior[k]) if prior[k] else None for k in ("0","1")},
      "PRemoval_flag_is_nonrandomly_related_to_regulatory_TP_monitoring_opportunity_POSSIBLE":True,
      "still_not_direct_native_2018_plant_ID_crosswalk":True,
    }


def make_frame(matches,plants,tps):
    rows=[]
    missing=collections.Counter()
    for sa,info in matches.items():
        native_plant=info["outlet_id"][1]
        plant=plants.get(native_plant)
        if plant is None:missing["source_plant_row_not_MAPPED_by_official_UWWTP_ID"]+=1;continue
        obs=tps.get(sa)
        if obs is None:missing["EA_SA_no_actual_2018_TP_numeric_observations"]+=1;continue
        row={
            "P_removal":plant["PRemoval"],
            "N_removal":plant["NRemoval"],
            "capacity":plant["capacity"],
            "region":plant["original_region"],
            "annual_TP":obs["annual_median_TP_mgL"],
            "n_days":obs["n_source_dates"],
            "n_months":obs["n_source_months"],
            "nearest_outlet_SA_distance":info["original2018_FW_outlet_to_SA_distance_m"]
        }
        rows.append(row)
    df=pd.DataFrame(rows)
    if df.empty:raise RuntimeError("True 2018 source-matched UWWTD treatment vs WQE TP frame empty")
    qa={
      "source_250m_mutual_original_outlet_SA_pairs":len(matches),
      "original_2018_eligible_real_final_TP_SA_sample_sites":len(tps),
      "2018_actual_TP_source_matched_plants":len(df),
      "source_missing_reason_counts":dict(missing),
      "PRemoval_source_flag_0_1_counts":{
         str(k):int(v) for k,v in df.P_removal.value_counts().items()},
      "NRemoval_source_flag_0_1_counts":{
         str(k):int(v) for k,v in df.N_removal.value_counts().items()},
      "source_true_2018_TP_sampling_months_median":
         float(df.n_months.median()),
      "source_months_at_least2":int((df.n_months>=2).sum()),
      "region_3letter_original_NUTS_counts":{
         str(k):int(v) for k,v in df.region.fillna("MISSING").value_counts().items()},
      "capacity_source_valid_by_PRemoval":{
        str(p):{
         "n_source":len(d),
         "n_original_positive_capacity":int(d.capacity.notna().sum()),
         "source_capacity_complete_fraction":float(d.capacity.notna().mean())
        } for p,d in df.groupby("P_removal")},
      "one_2018_sampoint_source_median_TP_per_original_EA_native_SA_ID":True,
      "historic_outlet_to_2018_SA_monitoring_point_mutual_unique_geography_not_ID":True,
    }
    qualified={
      "PRemoval_0_at_least30_2018_true_SA_TP":(df.P_removal==0).sum()>=MIN_GROUP_N,
      "PRemoval_1_at_least30_2018_true_SA_TP":(df.P_removal==1).sum()>=MIN_GROUP_N,
      "at_least3_original_NUTS3_regions_with10_sources":
       sum(df.region.value_counts()>=MIN_PER_REGION)>=MIN_REGIONS,
      "capacity_ge70pct_in_each_PRemoval_source_group":all(
        q.capacity.notna().mean()>=MIN_CAPACITY_VALID_FRACTION for _,q in df.groupby("P_removal")),
      "both_source_flags_are_actual_0_and_1":set(df.P_removal)=={0,1}
    }
    return df,qa,qualified

def modelling(df,qual):
    if not all(qual.values()):
        return {"status":"SOURCE_SUPPORT_STOP_BEFORE_MODELLING",
                "source_flags":qual,
                "no_PRemoval_TP_result_opened":True}
    d=df.copy()
    sub=d[d.capacity.notna() & d.region.notna()].copy()
    if min((sub.P_removal==p).sum() for p in (0,1))<MIN_GROUP_N:
        return {"status":"SOURCE_MODEL_STOP_CAPACITY_COMPLETECASE_TOO_SMALL",
                "n_after_source_covariate_listwise_deletion":len(sub)}
    source_group={
        p:sub[sub.P_removal==p].annual_TP.to_numpy(float)
        for p in (0,1)
    }
    med0=float(np.median(source_group[0]))
    med1=float(np.median(source_group[1]))
    # Keep zeros as real numeric nondetect; ratio is undefined when 0.
    diff=med1-med0
    ratio=med1/med0 if med0>0 else None
    rng=np.random.default_rng(RNG_SEED)
    bs0=source_group[0];bs1=source_group[1]
    bootdiff=[];bootratio=[]
    for _ in range(N_BOOT):
        x=np.median(rng.choice(bs1,len(bs1),replace=True))
        y=np.median(rng.choice(bs0,len(bs0),replace=True))
        bootdiff.append(float(x-y))
        bootratio.append(float(x/y) if y>0 else np.nan)
    bootratio=np.asarray(bootratio,float)
    X=pd.DataFrame({
        "const":np.ones(len(sub)),
        "PRemoval":sub.P_removal.to_numpy(float),
        "NRemoval":sub.N_removal.to_numpy(float),
        "log_capacity":np.log1p(sub.capacity.to_numpy(float))
    },index=sub.index)
    # Raw uwwNUTS area geography, not nominal region labels inferred
    # based on response quality. Omit the first fixed category.
    for key in sorted(sub.region.dropna().unique())[1:]:
        X["region_"+key]=(sub.region==key).astype(float)
    if np.linalg.matrix_rank(X.to_numpy(float))<X.shape[1]:
        model={"status":"MODEL_STOP_EXACT_DESIGN_COLLINEARITY"}
    else:
        fit=sm.OLS(np.log1p(sub.annual_TP.to_numpy(float)),X).fit(cov_type="HC3")
        model={"status":"NONCAUSAL_CONDITIONAL_OLS",
               "conditional_log1p_TP_PRemoval_beta":float(fit.params["PRemoval"]),
               "heteroskedastic_HC3_normal95":[float(x) for x in fit.conf_int().loc["PRemoval"]],
               "PRemoval_nominal_p":float(fit.pvalues["PRemoval"]),
               "n_geographically_matched_unique_2018_SA":len(sub),
               "n_source_regions":sub.region.nunique(),
               "original_source_controls":["PRemoval","NRemoval","log_capacity","original_uwwNUTS3"],
               "UWWTD_WWTW_treatment_is_not_randomized":True,
               "permit_facility_ID_still_not_joined":True}
    return {
      "status":"EXPLORATORY_SOURCE_TP_DIFFERENCE_CALCULATED",
      "n_source_PRemoval_0":len(source_group[0]),
      "n_source_PRemoval_1":len(source_group[1]),
      "PRemoval_0_median_yearly_final_TP_mgL":med0,
      "PRemoval_1_median_yearly_final_TP_mgL":med1,
      "PRemoval_1_minus_0_median_TP_mgL":diff,
      "PRemoval_1_divide_0_median_TP_ratio":ratio,
      "source_equal_SA_site_boot2000_median_diff95":[float(x) for x in np.quantile(bootdiff,[.025,.975])],
      "source_equal_SA_site_boot2000_median_ratio95":
        [float(x) for x in np.quantile(bootratio[np.isfinite(bootratio)],[.025,.975])] if np.isfinite(bootratio).any() else None,
      "original_source_treatment_specific_TP_quantile25_median_75":{
         str(p):[float(v) for v in np.quantile(source_group[p],[.25,.5,.75])]
         for p in (0,1)},
      "adjusted_model":model,
      "source_outcome_is_observed_final_effluent_total_P_concentration_not_river_P_mass":True
    }

def json_source_scalar(value):
    """Serialize original NumPy scalar support flags without changing their scientific meaning."""
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,np.ndarray):return value.tolist()
    raise TypeError("Unsupported original source audit JSON value: "+str(type(value)))


def main(args):
    out={"schema":"GESUI_20261121_UWWTD2018_UK_PRemoval_FLAG_VS_EA_FINAL_EFFLUENT_TP_2018_SOURCE_GEOMATCH_V1",
         "status":"SOURCE_STOP","posthoc_and_not_randomized":True,
         "outlet_to_SA_source_geo250m_not_verified_native_ID_permit_crosswalk":True,
         "source_only_before_BIO_effects":True,
         "not_WWTP_annual_P_mass_kg":True}
    try:
        historic,out_counts,geoqa=eea_geo.exact_source()
        plants,plantqa=plant_records()
        if len(plants)!=N_EXPECTED_UWWTP:
            raise RuntimeError("Original 1799 UK PRemoval Boolean source cohort shifted")
        current,ea_qa=ea.collect()
        matching,qa=geo.select_current_SA_2018_UWWTD_historic_outlets(
            historic,out_counts,current)
        if len(matching)!=N_EXPECTED_MUTUAL_SA:
            raise RuntimeError("Original mutual 2018 FW outlet-current SA 250m 930 source cohort drift")
        TP_annual,tpqa=original_TP_dates(sorted(matching))
        treatment_sampling_qa=original_sampling_selection_QA(matching,plants,TP_annual)
        cohort,sourceqa,sourcegates=make_frame(matching,plants,TP_annual)
        out.update({
            "status":"SOURCE_PASS_TREATMENT_FLAG_TP_COHORT_AUDITED",
            "source_2018_UWWTP_treatment_QA":plantqa,
            "source_geo_QA":qa,
            "official_EA_WQE_2018_TP_source_QA":tpqa,
            "original_treatment_group_TP_monitoring_selection_QA":treatment_sampling_qa,
            "source_matched_2018_cohort_QA":sourceqa,
            "preoutcome_source_gates":sourcegates
        })
        out["exploratory_TP_flag_response"]=modelling(cohort,sourcegates)
        if out["exploratory_TP_flag_response"]["status"]=="EXPLORATORY_SOURCE_TP_DIFFERENCE_CALCULATED":
            out["status"]="EXPLORATORY_UWWTD2018_PRemoval_EA_TP_CONCENTRATION_COMPARISON_NOT_CAUSAL"
    except Exception as exc:
        out.update({"status":"SOURCE_OR_MODEL_STOP","source_error":type(exc).__name__+": "+str(exc)[:650]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False,default=json_source_scalar)+"\n")
    print("GESUI_UWWTD2018_PREMOVAL_EA_FINAL_TP="+json.dumps({
       "status":out["status"],"error":out.get("source_error"),
       "original_PRemoval":out.get("source_2018_UWWTP_treatment_QA"),
       "TP_monitoring_selection_by_PRemoval":out.get("original_treatment_group_TP_monitoring_selection_QA"),
       "matched_source_cohort":out.get("source_matched_2018_cohort_QA"),
       "source_gates":out.get("preoutcome_source_gates"),
       "source_2018_totalP_comparison":out.get("exploratory_TP_flag_response"),
    },ensure_ascii=False,sort_keys=True,allow_nan=False,default=json_source_scalar),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
