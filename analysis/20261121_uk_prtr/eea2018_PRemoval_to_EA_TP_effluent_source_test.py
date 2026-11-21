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
import argparse,collections,csv,datetime,hashlib,io,json,math,sys,tempfile,urllib.request,zipfile
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
EXPECTED_SA_COHORT_SHA256='d55a7773bd85c59bf137dda4dac835cd0a1f77d8912a1164a422ac48ee7cee07'
MIN_GROUP_N=30
MIN_REGIONS=3
MIN_PER_REGION=10
MIN_CAPACITY_VALID_FRACTION=.70
N_BOOT=2000
RNG_SEED=20261121
MAX_RADIUS_M=250
UK_NUTS1_PREFIXES=frozenset(('UKC','UKD','UKE','UKF','UKG','UKH',
                           'UKI','UKJ','UKK','UKL','UKM','UKN'))

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
        # The original uwwNUTS field is finer-scale; a three-character
        # prefix is NUTS1, never NUTS3. Validate the actual UK prefix.
        original_region=region[:3] if region[:3] in UK_NUTS1_PREFIXES else None
        if not original_region:origin["original_2018_UK_unknown_region"]+=1
        load_text=str(r.get("uwwLoadEnteringUWWTP","")).strip()
        try:
            load=float(load_text)
            if not math.isfinite(load) or load<0:load=None
        except (TypeError,ValueError):load=None
        state_fields={k:str(v).strip() for k,v in r.items()
                      if 'state' in k.lower() or 'status' in k.lower()}
        plants[key]={"PRemoval":int(pflag),"NRemoval":int(nflag),
                     "capacity":float(cap) if cap is not None and cap>0 else None,
                     "original_region":original_region,
                     "original_full_uwwNUTS":region,
                     "source_state_fields":state_fields,
                     "load_entering_PE":load,
                     "source_report_year":REPORT_YEAR}
        origin["2018_P_removal_"+pflag]+=1
        origin["2018_N_removal_"+nflag]+=1
    return plants,{
        "original_2018_UK_UWWTP_source_n":len(selected),
        "source_2018_valid_treatment_rows":len(plants),
        "original_WWTW_capacity_and_region_source_quality":dict(origin),
        "EEA_original_UWWTPS_SHA256":eea.EXPECT["UWWTPS.csv"],
        "EEA_retrieved_archive_SHA256":sha.hexdigest(),
        "EEA_retrieved_archive_bytes":n,
        "region_control_level":"NUTS1 prefix of original uwwNUTS",
        "full_original_NUTS_code_length_counts":dict(collections.Counter(
            str(len(p['original_full_uwwNUTS'])) for p in plants.values())),
        "original_source_state_field_counts":{
            k:dict(collections.Counter(p['source_state_fields'].get(k,'') or 'MISSING'
                for p in plants.values()))
            for k in sorted({k for p in plants.values() for k in p['source_state_fields']})},
        "incoming_load_PE_source_quality_by_PRemoval":{
            str(g):{'n':sum(p['PRemoval']==g for p in plants.values()),
                'nonnegative_finite_reported':sum(p['PRemoval']==g and p['load_entering_PE'] is not None for p in plants.values())}
            for g in (0,1)},
        "incoming_load_PE_is_not_measured_P_or_N_mass":True,
        "official_Eionet_type":"uwwPRemoval and uwwNRemoval boolean source treatment fields",
        "not_measured_NP_load_or_randomized_treatment":True
    }

def original_TP_dates(names):
    """Original WQE page and source code verification before outcome aggregation."""
    if len(names)!=len(set(names)):
        raise RuntimeError('Duplicate requested SA candidate')
    records=collections.defaultdict(lambda:collections.defaultdict(list))
    stats=collections.Counter()
    private_site_stats={name:collections.Counter() for name in names}
    receipts=[]
    exact_rows_seen=set()
    batches=0
    pagecount=0
    for lo in range(0,len(names),wq.BATCH_N):
        codes=names[lo:lo+wq.BATCH_N]
        batches+=1
        for index in range(wq.MAX_BATCH_PAGES):
            rows,rec=wq.page(REPORT_YEAR,codes,index*wq.PAGE_LIMIT)
            receipts.append(dict(rec,batch_index=batches-1))
            pagecount+=1
            for r in rows:
                stats["returned_source_rows"]+=1
                sa_code=r['samplingPoint.notation'].strip()
                if sa_code not in codes:
                    raise RuntimeError('WQE returned a station outside the fixed candidate batch')
                fate=private_site_stats[sa_code]
                fate['returned_TP_rows']+=1
                if r["determinand.notation"].strip()!="0348":
                    raise RuntimeError("WQE returned wrong original TP determinand")
                if not r["phenomenonTime"].strip().startswith("2018"):
                    raise RuntimeError("WQE actual source row has wrong original reported sampling year")
                if r["samplingPoint.samplingPointType"].strip()!=original_wqe.SA:
                    stats["non_WaterCompany_final_type"]+=1;fate['wrong_type']+=1;continue
                if r["sampleMaterialType"].strip()!="FINAL SEWAGE EFFLUENT":
                    stats["non_final_treated_sample_material"]+=1;fate['wrong_material']+=1;continue
                if r["unit"].strip()!="MILLIGRAM PER LITRE":
                    stats["non_mg_per_L"]+=1;fate['wrong_unit']+=1;continue
                exact_key=tuple(sorted(r.items()))
                if exact_key in exact_rows_seen:
                    stats['exact_export_duplicate_rows_excluded']+=1
                    fate['exact_export_duplicate_rows_excluded']+=1
                    continue
                exact_rows_seen.add(exact_key)
                text_result=r['result'].strip()
                if text_result.startswith(('<','>','≤','≥')):
                    stats['bounded_censored_results']+=1;fate['bounded_censored_results']+=1;continue
                try:tp=float(text_result)
                except (ValueError,TypeError):
                    stats['malformed_or_nonnumeric_results']+=1
                    fate['malformed_or_nonnumeric_results']+=1;continue
                if not math.isfinite(tp) or tp<0:
                    stats["negative_or_nonfinite_original_TP"]+=1
                    fate['negative_or_nonfinite_original_TP']+=1;continue
                day=r["phenomenonTime"].strip()[:10]
                datetime.date.fromisoformat(day)
                records[sa_code][day].append(tp)
                stats["valid_numeric_TP_records"]+=1
                fate['eligible_numeric_TP_records']+=1
                if tp==0:stats['numeric_zero_values']+=1
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
        stats['site_days_with_multiple_distinct_export_rows']+=sum(len(v)>1 for v in days.values())
    stats.update({"distinct_numeric_TP_SA_stations":len(siteyear),
                  "official_WQE_250row_page_n":pagecount,
                  "official_source_batches":batches,
                  "TP_per_SA_year_is_median_of_source_day_medians":True})
    stats['censored_or_nonnumeric']=stats['bounded_censored_results']+stats['malformed_or_nonnumeric_results']
    source_fates={}
    for code,counts in private_site_stats.items():
        label=('eligible_numeric_TP' if code in siteyear else
               'no_eligible_numeric_TP' if counts['returned_TP_rows'] else
               'no_returned_TP_record')
        source_fates[code]=dict(counts,fate=label)
    qa=dict(stats)
    qa.update(source_page_receipts=receipts,
              candidate_site_fate_counts=dict(collections.Counter(v['fate'] for v in source_fates.values())),
              page_size_requested=wq.PAGE_LIMIT,
              all_batches_reached_a_terminal_short_page=True,
              numeric_zero_is_recorded_value_not_inferred_nondetect=True,
              duplicate_audit_scope='export rows passing SA type, final-effluent material and mg/L unit checks',
              exact_export_duplicates_do_not_establish_independent_lab_replicates=True)
    return siteyear,qa,source_fates

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
            "outlet_easting":float(info['outlet_coord'][0]),
            "outlet_northing":float(info['outlet_coord'][1]),
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
      "at_least3_original_NUTS1_prefix_regions_with10_sources":
       sum(df.region.value_counts()>=MIN_PER_REGION)>=MIN_REGIONS,
      "capacity_ge70pct_in_each_PRemoval_source_group":all(
        q.capacity.notna().mean()>=MIN_CAPACITY_VALID_FRACTION for _,q in df.groupby("P_removal")),
      "both_source_flags_are_actual_0_and_1":set(df.P_removal)=={0,1}
    }
    return df,qa,qualified

def modelling(df,qual):
    if min((df.P_removal==p).sum() for p in (0,1))<MIN_GROUP_N:
        return {"status":"SOURCE_SUPPORT_STOP_BEFORE_MODELLING",
                "source_flags":qual,
                "no_PRemoval_TP_result_opened":True}
    d=df.copy()
    complete=(d.capacity.notna() & (d.capacity>0) & d.region.notna()
              & np.isfinite(d.outlet_easting) & np.isfinite(d.outlet_northing))
    sub=d[complete].copy()
    model_supported=all(qual.values()) and min((sub.P_removal==p).sum() for p in (0,1))>=MIN_GROUP_N
    source_group={
        p:d[d.P_removal==p].annual_TP.to_numpy(float)
        for p in (0,1)
    }
    med0=float(np.median(source_group[0]))
    med1=float(np.median(source_group[1]))
    # Keep source-recorded numeric zero; it does not itself prove nondetection.
    # The median ratio is undefined if the denominator is zero.
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
        "log_capacity":np.log1p(sub.capacity.to_numpy(float)),
        "historical_outlet_easting_100km":
            (sub.outlet_easting.to_numpy(float)-sub.outlet_easting.mean())/100000,
        "historical_outlet_northing_100km":
            (sub.outlet_northing.to_numpy(float)-sub.outlet_northing.mean())/100000,
    },index=sub.index)
    # Three-character prefix of original uwwNUTS = NUTS1, not NUTS3.
    for key in sorted(sub.region.dropna().unique())[1:]:
        X["region_"+key]=(sub.region==key).astype(float)
    if not model_supported:
        model={"status":"SOURCE_MODEL_STOP_ADJUSTMENT_SUPPORT",
               "source_flags":qual,
               "n_after_source_covariate_listwise_deletion":len(sub)}
    elif np.linalg.matrix_rank(X.to_numpy(float))<X.shape[1]:
        model={"status":"MODEL_STOP_EXACT_DESIGN_COLLINEARITY"}
    else:
        fit=sm.OLS(np.log1p(sub.annual_TP.to_numpy(float)),X).fit(cov_type="HC3")
        model={"status":"NONCAUSAL_CONDITIONAL_OLS",
               "conditional_log1p_TP_PRemoval_beta":float(fit.params["PRemoval"]),
               "heteroskedastic_HC3_normal95":[float(x) for x in fit.conf_int().loc["PRemoval"]],
               "PRemoval_nominal_p":float(fit.pvalues["PRemoval"]),
               "n_geographically_matched_unique_2018_SA":len(sub),
               "n_source_regions":sub.region.nunique(),
               "original_source_controls":["PRemoval","NRemoval","log_capacity",
                   "original_uwwNUTS1_prefix","historical_outlet_easting_100km",
                   "historical_outlet_northing_100km"],
               "restores_original_protocol_coordinates_omitted_in_historical_run":True,
               "UWWTD_WWTW_treatment_is_not_randomized":True,
               "permit_facility_ID_still_not_joined":True}
    # Posthoc diagnostic, explicitly separate from the pre-outcome frozen
    # full adjustment. This sequence shows whether composition by geography
    # and plant size explains the unadjusted treatment-flag contrast; it is
    # NOT a mediation sequence and is not a causal effect decomposition.
    nested={}
    y=np.log1p(sub.annual_TP.to_numpy(float))
    diagnostic_specs=(
        ("flag_only",()),
        ("flag_and_capacity",("capacity",)),
        ("flag_and_region",("region",)),
        ("flag_capacity_region",("capacity","region")),
        ("flag_capacity_region_NRemoval",("capacity","region","NRemoval")),
        ("flag_capacity_region_NRemoval_BNG",("capacity","region","NRemoval","BNG"))
    ) if model_supported else ()
    for label,parts in diagnostic_specs:
        Z=pd.DataFrame({"const":np.ones(len(sub)),
                        "PRemoval":sub.P_removal.to_numpy(float)},index=sub.index)
        if "capacity" in parts:
            Z["log_capacity"]=np.log1p(sub.capacity.to_numpy(float))
        if "NRemoval" in parts:
            Z["NRemoval"]=sub.N_removal.to_numpy(float)
        if "BNG" in parts:
            Z['historical_outlet_easting_100km']=(sub.outlet_easting-sub.outlet_easting.mean())/100000
            Z['historical_outlet_northing_100km']=(sub.outlet_northing-sub.outlet_northing.mean())/100000
        if "region" in parts:
            for rk in sorted(sub.region.dropna().unique())[1:]:
                Z["region_"+rk]=(sub.region==rk).astype(float)
        if np.linalg.matrix_rank(Z.to_numpy(float))<Z.shape[1]:
            nested[label]={"status":"MODEL_STOP_COLLINEAR"}
        else:
            mfit=sm.OLS(y,Z).fit(cov_type="HC3")
            nested[label]={
                "status":"POSTHOC_NONCAUSAL_ADJUSTMENT_DIAGNOSTIC",
                "PRemoval_beta_log1p_TP":float(mfit.params["PRemoval"]),
                "PRemoval_site_HC3_normal95":[float(v) for v in mfit.conf_int().loc["PRemoval"]],
                "nominal_p":float(mfit.pvalues["PRemoval"]),
                "n":int(len(sub)) }
    return {
      "status":"EXPLORATORY_SOURCE_TP_DIFFERENCE_CALCULATED",
      "posthoc_sequential_source_adjustment_diagnostics":nested,
      "primary_descriptive_n_all_eligible_TP":len(d),
      "adjusted_model_completecase_n":len(sub),
      "excluded_from_adjusted_model_missing_covariate_n":len(d)-len(sub),
      "median_contrast_and_log1p_mean_coefficients_are_different_estimands":True,
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
    out={"schema":"GESUI_20261121_UWWTD2018_UK_PRemoval_FLAG_VS_EA_FINAL_EFFLUENT_TP_2018_SOURCE_GEOMATCH_V2",
         "status":"SOURCE_STOP","posthoc_and_not_randomized":True,
         "outlet_to_SA_source_geo250m_not_verified_native_ID_permit_crosswalk":True,
         "source_only_before_BIO_effects":True,
         "not_WWTP_annual_P_mass_kg":True}
    safe_run_metadata=dict(out)
    stage='historical_outlet_source'
    try:
        historic,out_counts,geoqa=eea_geo.exact_source()
        stage='historical_plant_source'
        plants,plantqa=plant_records()
        if len(plants)!=N_EXPECTED_UWWTP:
            raise RuntimeError("Original 1799 UK PRemoval Boolean source cohort shifted")
        stage='current_EA_registry'
        current,ea_qa=ea.collect()
        stage='candidate_source_cohort'
        matching,qa=geo.select_current_SA_2018_UWWTD_historic_outlets(
            historic,out_counts,current)
        if len(matching)!=N_EXPECTED_MUTUAL_SA:
            raise RuntimeError("Original mutual 2018 FW outlet-current SA 250m 930 source cohort drift")
        cohort_hash=hashlib.sha256('\n'.join(sorted(matching)).encode()).hexdigest()
        if cohort_hash!=EXPECTED_SA_COHORT_SHA256:
            raise RuntimeError('Original 930 SA station-set fingerprint changed')
        qa['baseline_run_44850406458_SA_station_set_SHA256_verified']=True
        qa['complete_pair_mapping_has_no_retrospective_baseline_fingerprint']=True
        stage='original_2018_TP_records'
        TP_annual,tpqa,site_fates=original_TP_dates(sorted(matching))
        stage='source_support_aggregation'
        treatment_sampling_qa=original_sampling_selection_QA(matching,plants,TP_annual)
        from premoval_source_support import build_source_support
        support=build_source_support(matching,plants,TP_annual,site_fates)
        cohort,sourceqa,sourcegates=make_frame(matching,plants,TP_annual)
        out.update({
            "status":"SOURCE_PASS_TREATMENT_FLAG_TP_COHORT_AUDITED",
            "source_2018_UWWTP_treatment_QA":plantqa,
            "source_2018_historical_outlet_QA":geoqa,
            "source_current_EA_registry_QA":ea_qa,
            "source_geo_QA":qa,
            "source_selection_and_covariate_support":support,
            "official_EA_WQE_2018_TP_source_QA":tpqa,
            "original_treatment_group_TP_monitoring_selection_QA":treatment_sampling_qa,
            "source_matched_2018_cohort_QA":sourceqa,
            "preoutcome_source_gates":sourcegates
        })
        stage='source_TP_model'
        out["exploratory_TP_flag_response"]=modelling(cohort,sourcegates)
        if out["exploratory_TP_flag_response"]["status"]=="EXPLORATORY_SOURCE_TP_DIFFERENCE_CALCULATED":
            out["status"]="EXPLORATORY_UWWTD2018_PRemoval_EA_TP_CONCENTRATION_COMPARISON_NOT_CAUSAL"
        stage='serialize_aggregate'
        payload=json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False,default=json_source_scalar)+'\n'
    except Exception as exc:
        # Third-party exception messages can contain native IDs or records.
        # Publish only a controlled stage and exception class; clear potentially
        # non-serializable partial model output before writing the STOP receipt.
        out=dict(safe_run_metadata,status='SOURCE_OR_MODEL_STOP',
                 source_error=type(exc).__name__,source_error_stage=stage,
                 source_error_detail_redacted=True)
        payload=json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+'\n'
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(payload)
    print("GESUI_UWWTD2018_PREMOVAL_EA_FINAL_TP="+json.dumps({
       "status":out["status"],"error":out.get("source_error"),
       "original_PRemoval":out.get("source_2018_UWWTP_treatment_QA"),
       "TP_monitoring_selection_by_PRemoval":out.get("original_treatment_group_TP_monitoring_selection_QA"),
       "matched_source_cohort":out.get("source_matched_2018_cohort_QA"),
       "source_gates":out.get("preoutcome_source_gates"),
       "source_2018_totalP_comparison":out.get("exploratory_TP_flag_response"),
    },ensure_ascii=False,sort_keys=True,allow_nan=False,default=json_source_scalar),flush=True)
    # Expected source-support STOP remains an explicit data result. Unexpected
    # retrieval, schema, coding or model failures must fail the Actions step.
    return 1 if 'source_error' in out else 0

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    raise SystemExit(main(p.parse_args()))

