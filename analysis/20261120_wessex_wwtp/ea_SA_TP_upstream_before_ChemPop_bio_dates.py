#!/usr/bin/env python3
"""Source-only true sampling DATE order for official treated STW TP and ChemPop.

Crosswalk exact 2017/18 original ChemPop EIDC native site×BIO_SAMPLE_DATE
(source-eligibility only, NO actual BMWP trait outcome) against original
EA WQE water-company final treated effluent SA total P observations 0348
at upstream 50km HydroRIVERS NEXT_DOWN source reach, retaining current
national 5722 SA registry and sampling point != actual physical outfall.

Compare pre-BIO 30/90/365 day TP observation records and after-BIO90
comparison. Do NOT infer measured wastewater P mass/year, hydraulic
delivery, causal vegetation/biotic response, source N/P origin or ecology.
"""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import sys
import tempfile
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as nat
import ea_SA_TP_to_ChemPop_location_sourcegate as local
import ea_SA_TP_ChemPop_HydroRIVERS_directed_route_gate as route
import ea_SA_TP_ChemPop_networkfirst_allSA_history as net
import ea_WQE_national_SA_actual2017_18_TP_chemistry_bio_screen as hist

sys.path.insert(0,str(ROOT/"analysis"/"20261119_europe_benthic_taxa"))
import chempop_v2_source_census as eidc
import chempop_v2_bmwp_filter as srcbio

YEARS=(2017,2018)
WINDOWS_DAYS=(30,90,365)
FUTURE_WINDOW_DAYS=90
EXPECTED_INITIAL_BIO_UNION=148
EXPECTED_INITIAL_BIO_SITEYEARS=237
EXPECTED_ORIGINAL_BIO_SAMPLES=7836
EXPECTED_ORIGINAL_BIO_SITES=707
EXPECTED_SCORE_SHA="4f9be755ce6d24e645b46705dbec9c0367bc71a9ce6b335d3bf2657d8c10ac5a"


def origin_bio_sample_dates():
    """Fresh official-source reconstruction; no scored-community BMWP outcomes."""
    with tempfile.TemporaryDirectory(prefix="gesui_BIO_source_date_") as tmp:
        d=Path(tmp)
        archive=eidc.retrieve(eidc.EIDC_URL,d/"original.zip",eidc.SOURCE_DOWNLOAD_LIMIT)
        score,dictionary=eidc.score_reference(d)
        if dictionary["SHA256"]!=EXPECTED_SCORE_SHA:
            raise RuntimeError("Original BMWP source dictionary changed")
        with zipfile.ZipFile(d/"original.zip") as z:
            if len(z.namelist())!=326:raise RuntimeError("Original ChemPop ZIP dataset membership drift")
            rows,weights,orders,audit=srcbio.construct_eligible(z,z.namelist(),score)
    if (len(rows)!=EXPECTED_ORIGINAL_BIO_SAMPLES or
        audit["source_gate"]["eligible_site_n"]!=EXPECTED_ORIGINAL_BIO_SITES):
        raise RuntimeError("Original 7836/707 score-eligible BIO source cohort changed")
    d=collections.defaultdict(set)
    for row in rows:
        yr=row["year"]
        if yr not in YEARS:continue
        key=(yr,row["site_region_id"])
        original=row.get("sample_date")
        if not original:
            raise RuntimeError("EIDC original BIO source actual sample date missing")
        day=dt.date.fromisoformat(str(original)[:10])
        if day.year!=yr:raise RuntimeError("Original biological sample year/date conflict")
        d[key].add(day)
    return d,{
        "original_EIDC_326file_source_zip_SHA256":archive["SHA256"],
        "original_BMWP_source_reference_SHA256":dictionary["SHA256"],
        "n_original_BIO_2017_18_unique_source_site_dates":sum(len(x) for x in d.values()),
        "n_2017_BIO_siteyears":sum(key[0]==2017 for key in d),
        "n_2018_BIO_siteyears":sum(key[0]==2018 for key in d),
        "all_native_dates_truly_from_original_EIDC_BIO_sample":True,
        "no_BIO_response_BMWP_or_abundance_CWM_calculated":True,
    }


def upstream_source_ids(rel,measured):
    bybio=collections.defaultdict(set)
    for year in YEARS:
        for sa,v in rel[year].items():
            if sa not in measured[year]:continue
            for biokey in v["upstream"]:
                bybio[(year,biokey)].add(sa)
    return bybio


def historical_sa_tp_dates(measured):
    all_days=collections.defaultdict(set)
    for year in YEARS:
        for code,v in measured[year].items():
            for sampleday in v["valid_days"]:
                d=dt.date.fromisoformat(str(sampleday)[:10])
                if d.year!=year:raise RuntimeError("Genuine STW TP source sample year mismatch")
                all_days[code].add(d)
    return all_days


def evaluate_windows(bio_dates,src_by_bio,sa_tp_days):
    """Count ACTUAL upstream TP assays preceding all source-qualified BIO dates."""
    flags={
       "pre30":30,"pre90":90,"pre365":365,"post90":90,"same_day":0,
    }
    year_result={}
    aggregate={k:{"union_bio":set(),"site_year_count":0,
                  "source_bio_site_dates_supported":0,
                  "all_bio_dates_supported_siteyears":0,
                  "n_original_bio_site_dates":0}
               for k in flags}
    total_bio_source_dates=0
    for yr in YEARS:
        eligible_keys=sorted(k for k in src_by_bio if k[0]==yr)
        checked={}
        for name in flags:
            checked[name]={"siteyear_any":0,"siteyear_all":0,
                           "n_supported_bio_site_dates":0,
                           "n_total_bio_site_dates":0,
                           "distinct_bio_sites":set()}
        for year,bid in eligible_keys:
            days=sorted(bio_dates.get((year,bid),set()))
            if not days:raise RuntimeError("Network BIO site lacked ORIGINAL sample-date support")
            source_ids=src_by_bio[(year,bid)]
            if not source_ids:raise RuntimeError("Empty genuine upstream STW source set")
            measured_days=sorted(set().union(*(sa_tp_days.get(s,set()) for s in source_ids)))
            if not measured_days:raise RuntimeError("An original TP-measured STW source has zero official sample dates")
            total_bio_source_dates+=len(days)
            for flag,daysize in flags.items():
                supported=0
                for bio_day in days:
                    if flag.startswith("pre"):
                        # Source timestamps and BIO survey records do not
                        # establish the WITHIN-DAY order. Same-calendar-day
                        # TP is NOT evidence of a preceding source sample.
                        passed=any(0<(bio_day-date).days<=daysize
                                   for date in measured_days)
                    elif flag=="same_day":
                        passed=bio_day in measured_days
                    else:
                        passed=any(0<(date-bio_day).days<=daysize
                                   for date in measured_days)
                    supported+=int(passed)
                checked[flag]["n_total_bio_site_dates"]+=len(days)
                checked[flag]["n_supported_bio_site_dates"]+=supported
                if supported:
                    checked[flag]["siteyear_any"]+=1
                    checked[flag]["distinct_bio_sites"].add(bid)
                if supported==len(days):
                    checked[flag]["siteyear_all"]+=1
        year_result[str(yr)]={}
        for flag,v in checked.items():
            year_result[str(yr)][flag]={
                "original_same_year_upstream_candidate_bio_siteyears":len(eligible_keys),
                "site_year_with_any_BIO_date_preceded_or_followed_by_an_original_TP_sample":
                    v["siteyear_any"],
                "site_year_where_ALL_BIO_dates_have_an_original_TP_sample_in_window":
                    v["siteyear_all"],
                "distinct_BIO_sites_in_window":len(v["distinct_bio_sites"]),
                "original_BIO_site_date_count":v["n_total_bio_site_dates"],
                "BIO_site_dates_with_any_STW_TP_sample_in_window":v["n_supported_bio_site_dates"],
            }
            agg=aggregate[flag]
            agg["union_bio"].update(v["distinct_bio_sites"])
            agg["site_year_count"]+=v["siteyear_any"]
            agg["all_bio_dates_supported_siteyears"]+=v["siteyear_all"]
            agg["source_bio_site_dates_supported"]+=v["n_supported_bio_site_dates"]
            agg["n_original_bio_site_dates"]+=v["n_total_bio_site_dates"]
    overall={flag:{
        "unique_BIO_sites_with_any_date_in_window":len(v["union_bio"]),
        "site_years_with_any_true_BIO_date_in_window":v["site_year_count"],
        "site_years_where_ALL_true_BIO_dates_in_window":v["all_bio_dates_supported_siteyears"],
        "supported_BIO_site_date_pairs":v["source_bio_site_dates_supported"],
        "total_BIO_site_date_pairs":v["n_original_bio_site_dates"],
        "this_is_TP_monitoring_temporal_opportunity_not_wastewater_delivery":True,
    } for flag,v in aggregate.items()}
    if any(v["site_years_with_any_true_BIO_date_in_window"]>
           EXPECTED_INITIAL_BIO_SITEYEARS for v in overall.values()):
        raise RuntimeError("Temporal original source evidence exceeds total network-upstream BIO siteyears")
    return year_result,overall


def main(a):
    out={
        "schema":"GESUI_20261120_GENUINE_TP_STW_SA_MEASURED_BEFORE_BIO_SAMPLED_SOURCE_GATE_V1",
        "source_only_no_actual_BMWP_or_biomass_analysed":True,
        "current_2026_SA_not_historical_STW_outfall_or_annual_kgP":True,
        "years_with_official_TP_sources":list(YEARS),
        "pre_BIO_TP_sample_windows_in_days":list(WINDOWS_DAYS),
        "post_BIO_future_comparison_days":FUTURE_WINDOW_DAYS,
        "2017_left_truncated_pre_365_day_TP_history_missing_2016":True,
        "not_actual_treatment_plant_effluent_outlet_route":True,
        "status":"SOURCE_STOP"
    }
    try:
        registry,regqa=nat.collect()
        bio,yr,bioqa=local.chempop_geo()
        with tempfile.TemporaryDirectory(prefix="gesui_temporal_TP_hydro_") as temp:
            rivers,graph,hydroqa=route.get_rivers(Path(temp))
            b,bqa=route.snap({k:bio[k] for k in yr[2017]|yr[2018]},
                             rivers,"original genuine eligible ChemPop river sites")
            sa,saqa=route.snap(registry,rivers,"current WQE WaterCompany SA sampling sites")
            eligible,rel,prior=net.relations(sa,b,yr,graph)
            measured,history=hist.candidate_history(sorted(eligible))
            original_year,original=net.validate_measured(measured,rel,graph,b,yr)
        if (original["n_DISTINCT_true_BIO_site_union_network_upstream_measured_source"]!=EXPECTED_INITIAL_BIO_UNION
            or original["n_SITE_YEAR_bio_network_upstream_measured_source"]!=EXPECTED_INITIAL_BIO_SITEYEARS):
            raise RuntimeError("Rebuilt official 148/237 historical TP SA upstream source gates changed")
        true_bio_dates,biodateqa=origin_bio_sample_dates()
        for year in YEARS:
            compare=set(k for yr0,k in true_bio_dates if yr0==year)
            if compare != set(yr[year]):
                raise RuntimeError("Native original BIO survey year/sample dates disagree with source registration")
        bybio=upstream_source_ids(rel,measured)
        if len(bybio)!=EXPECTED_INITIAL_BIO_SITEYEARS:
            raise RuntimeError("True source route 237 site-year keys not recovered from original EA source")
        tps=historical_sa_tp_dates(measured)
        peryear,total=evaluate_windows(true_bio_dates,bybio,tps)
        out.update({
            "source_recapture":{
                "original_candidate_current_SA_source_n":len(eligible),
                "EA_total_P_history":{k:v for k,v in history.items() if k not in ("source_page_SHA256",)},
                "original_148_distinct_BIO_237_SITEYEAR_upstream_source_geometry_reproduced":True,
                "native_BIO_sampling_DATE_qa":biodateqa,
                "origin_WQE_SA_TP_dates_from_original_real_2017_18_sources":True
            },
            "true_temporal_order_before_BIO_by_year":peryear,
            "true_temporal_order_overall":total,
            "source_gates":{
                "same_year_source_reach_count_reproduced":True,
                "original_BIO_dates_from_actual_sample_date_not_year_proxy":True,
                "pre90_any_BIO_source_observation_siteyears_ge_100":
                    total["pre90"]["site_years_with_any_true_BIO_date_in_window"]>=100,
                "pre90_distinct_BIO_sites_ge_50":
                    total["pre90"]["unique_BIO_sites_with_any_date_in_window"]>=50,
                "pre90_all_BIO_dates_source_observation_siteyears_ge_100":
                    total["pre90"]["site_years_where_ALL_true_BIO_dates_in_window"]>=100,
                "same_day_measurement_not_assumed_to_precede_BIO":True,
                "pre365_year2017_left_censored_due_to_unretrieved2016_TP":True,
                "true_STW_outfall_location_and_actual_river_flow_delivered":False,
                "real_WWTP_P_mass_kg_year":False,
                "allow_BMWP_causal_model":False
            },
            "status":"SOURCE_NETWORK_AND_TP_SAMPLE_DATES_COMPLETED_NOT_CAUSAL"
        })
    except Exception as e:
        out.update({"status":"SOURCE_OR_TEMPORAL_GATE_STOP",
                    "error":type(e).__name__+": "+str(e)[:750]})
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_TP_STW_MEASUREMENT_BEFORE_CHEMPOP_BIO="+json.dumps({
        "status":out["status"],"error":out.get("error"),
        "source":out.get("source_recapture"),
        "year":out.get("true_temporal_order_before_BIO_by_year"),
        "overall":out.get("true_temporal_order_overall"),
        "gates":out.get("source_gates"),
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
