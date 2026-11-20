#!/usr/bin/env python3
"""Historical 2017/2018 measured treated-effluent TP in nearby national SA sites.

Source-only + strict year validation. Source geographical candidates are
the complete 2026 WQE SA registered sampling-points within two km of any
ChemPop original source-quality-eligible BIO site sampled in 2017 or 2018.
The official observation API must provide original TOTAL PHOSPHORUS as P
0348, sample point type Water Company FINAL/TREATED effluent (SA), sample
material FINAL SEWAGE EFFLUENT, numeric TP, unit mg/L, and source year
equal to BIO sample year. On this background, quantify spatial opportunity
at 0.5, 1 and 2 km from each genuine measured TP SA point in the SAME YEAR.

NO biological outcome value or slope; NO direct topology or permit-outfall
direction, NO source-specific effluent P mass, and no unrecorded discharge
is imputed zero. No original site names, IDs, raw samples or concentrations
are published in results, just aggregate source counts and SHA receipts.

Official EA 2026 API supports comma-separated pointNotation, NOT repeated
pointNotation key: confirmed using two independently observed original
Bath-area effluent monitoring points (run #44678384278).
"""
from __future__ import annotations

import argparse,collections,csv,hashlib,io,json,math,sys,time,urllib.error,urllib.parse,urllib.request
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis"/"20261120_wessex_wwtp"))
import ea_national_SA_registry_to_ChemPop_source_ceiling as national
import ea_SA_TP_to_ChemPop_location_sourcegate as local
import ea_wqe_2017_2018_SA_tp_siteyear_census as actual

WQE_OBS="https://environment.data.gov.uk/water-quality/data/observation"
UA="gesui-original-EA-WQE-2017-18-treated-SA-TP-near-ChemPop-20261120"
HEADERS={
    "User-Agent":UA,"API-Version":"1","Accept":"text/csv",
    "CSV-Header":"present","Content-Type":"application/json"}
YEARS=(2017,2018)
SEARCH_RADIUS_M=2000
SCREENS_M=(500,1000,2000)
BATCH_N=40
PAGE_LIMIT=250  # EA official WQE API max 250 observations per page (2025-11 documented)
MAX_BATCH_PAGES=12
MAX_BATCH_REPLY_BYTES=750000
EXPECTED_SOURCE_N=5722
EXPECTED_BIO_BY_YEAR={2017:389,2018:366}
MIN_SA_OBS_DATES=1
MULTI_SITE_CONTRACT_RUN="44678384278"


def source_candidates(source,survey,biosites):
    sites=sorted(survey[2017] | survey[2018])
    point_names=sorted(source)
    if len(sites)!=463:raise RuntimeError("Original 2017/18 scored ChemPop site union changed")
    b=np.asarray([biosites[k] for k in sites],float)
    s=np.asarray([source[k] for k in point_names],float)
    distance=np.linalg.norm(s[:,None,:]-b[None,:,:],axis=2)
    inside=distance.min(axis=1)<=SEARCH_RADIUS_M
    names=[point_names[i] for i in np.flatnonzero(inside)]
    return names,{
        "registered_SA_current_source_n":len(point_names),
        "original_bio_site_n_either_2017_or_2018":len(sites),
        "candidate_current_SA_sites_within_2km_source_screen":len(names),
        "criterion_is_year_agnostic_geography_only":True,
        "source_candidate_sampling_code_SHA256":hashlib.sha256(
            "\n".join(names).encode()).hexdigest(),
    }


def page(year,notations,skip):
    request_params={
        "dateFrom":f"{year}-01-01","dateTo":f"{year}-12-31",
        "determinand":"0348","pointNotation":",".join(notations),
        "limit":PAGE_LIMIT,"skip":skip
    }
    request=urllib.request.Request(
        WQE_OBS+"?"+urllib.parse.urlencode(request_params),
        method="POST",data=b"",headers=HEADERS)
    receipt={
        "year":year,"n_requested_native_SA_points":len(notations),
        "skip":skip,"status":None,"n_observations":None,
        "sha256":None,"content_type":None}
    try:
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request,timeout=70) as response:
                    receipt["status"]=int(response.status)
                    receipt["content_type"]=response.headers.get("Content-Type","")
                    raw=response.read(MAX_BATCH_REPLY_BYTES+1)
                break
            except urllib.error.HTTPError as exc:
                if exc.code in (502,503,504) and attempt<2:
                    time.sleep(1.5*(attempt+1))
                    continue
                raise
        if len(raw)>MAX_BATCH_REPLY_BYTES:raise RuntimeError("Official 250-observation WQE page over capped 750 KB")
        receipt["source_bytes"]=len(raw)
        receipt["sha256"]=hashlib.sha256(raw).hexdigest()
        if receipt["status"]==204 or not raw:
            receipt["n_observations"]=0
            return [],receipt
        if receipt["status"]!=200:
            raise RuntimeError("Unexpected original WQE source status")
        reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")))
        required={"samplingPoint.notation","samplingPoint.samplingPointType",
          "sampleMaterialType","determinand.notation","phenomenonTime",
          "result","unit","samplingPurpose"}
        if not required.issubset(set(reader.fieldnames or [])):
            raise RuntimeError("Official EA original native observation schema changed")
        rows=list(reader)
        receipt["n_observations"]=len(rows)
        if len(rows)>PAGE_LIMIT:
            raise RuntimeError("Official source response exceeds official WQE 250 observations per page")
        requested=set(notations)
        if any(x["samplingPoint.notation"].strip() not in requested for x in rows):
            raise RuntimeError("Official pointNotation comma filter is not strict")
        return rows,receipt
    except urllib.error.HTTPError as err:
        raise RuntimeError(
           f"WQE original source HTTP {err.code} on fixed source batch, year={year}, skip={skip}: {str(err)[:200]}") from err


def candidate_history(names):
    out={yr:collections.defaultdict(lambda:{
        "valid_days":set(),"valid_months":set(),"observations":0,
        "UWWTD_observations":0}) for yr in YEARS}
    receipts=[]
    status=collections.Counter()
    all_batches=0
    for i in range(0,len(names),BATCH_N):
        batch=names[i:i+BATCH_N]
        if len(batch)>100:raise RuntimeError("WQE source max 100 official pointNotations violated")
        all_batches+=1
        for yr in YEARS:
            for page_num in range(MAX_BATCH_PAGES):
                rows,rec=page(yr,batch,page_num*PAGE_LIMIT)
                receipts.append(rec)
                status["all_source_observations_returned"]+=len(rows)
                for row in rows:
                    code=row["samplingPoint.notation"].strip()
                    if row["determinand.notation"].strip()!="0348":
                        raise RuntimeError("Official source returned different determinand")
                    date=row["phenomenonTime"].strip()[:10]
                    if not date.startswith(str(yr)):
                        raise RuntimeError("Original WQE source year mismatch")
                    if row["samplingPoint.samplingPointType"].strip()!=actual.SA:
                        status["source_record_not_WaterCompany_SA"]+=1;continue
                    if row["sampleMaterialType"].strip()!="FINAL SEWAGE EFFLUENT":
                        status["source_record_not_final_sewage_effluent"]+=1;continue
                    if row["unit"].strip()!="MILLIGRAM PER LITRE":
                        status["source_record_not_mgL"]+=1;continue
                    try:tp=float(row["result"].strip())
                    except (ValueError,TypeError):
                        status["source_TP_censored_or_non_numeric"]+=1;continue
                    if not math.isfinite(tp) or tp<0:
                        status["source_TP_nonfinite_or_negative"]+=1;continue
                    value=out[yr][code]
                    value["observations"]+=1
                    value["valid_days"].add(date)
                    value["valid_months"].add(date[:7])
                    if "UWWTD" in row["samplingPurpose"].upper():
                        value["UWWTD_observations"]+=1
                    status["valid_WaterCompany_final_TP_mgL"]+=1
                if len(rows)<PAGE_LIMIT:break
                if page_num==MAX_BATCH_PAGES-1:
                    raise RuntimeError("Source per 40 STW/one year pagination truncated")
    summary={
        "total_source_registered_SA_candidates_screened":len(names),
        "source_batch_size_native_pointNotations":BATCH_N,
        "official_EA_WQE_observation_max_250_per_page_enforced":True,
        "source_batches":all_batches,
        "source_WQE_requested_years":list(YEARS),
        "source_page_count":len(receipts),
        "source_page_SHA256":[r["sha256"] for r in receipts],
        "source_record_fates":dict(status),
        "historically_measured_STW_sampling_points_by_year":{
            str(yr):len(out[yr]) for yr in YEARS
        },
        "historically_measured_STW_sampling_points_both_years":
             len(set(out[2017])&set(out[2018])),
        "historically_measured_STW_sampling_points_by_year_with_2_months":{
            str(yr):sum(len(v["valid_months"])>=2 for v in out[yr].values())
            for yr in YEARS},
        "historically_measured_STW_sampling_points_by_year_with_2_dates":{
            str(yr):sum(len(v["valid_days"])>=2 for v in out[yr].values())
            for yr in YEARS},
        "years_are_same_original_as_biological_study_and_NOT_2026_registry_activity":True,
        "source_WQE_is_concentration_mgL_not_annual_P_kg":True,
    }
    return out,summary


def compare_year(national,obs,year,geom,survey):
    year_set=sorted(survey[year])
    stw=sorted(obs[year])
    if len(year_set)!=EXPECTED_BIO_BY_YEAR[year]:
        raise RuntimeError("Original year-specific biological source status changed")
    if not stw:
        return {"n_bio_sites":len(year_set),"n_STW_historic_measured_TP":0,
                "status":"NO_MEASURED_STW_SOURCE"}
    a=np.asarray([geom[k] for k in year_set],float)
    b=np.asarray([national[k] for k in stw],float)
    dist=np.linalg.norm(a[:,None,:]-b[None,:,:],axis=2)
    out={
        "source_observation_year":year,
        "n_biology_original_sourceeligible_sites":len(year_set),
        "n_official_STW_current_sampling_points_WITH_ACTUAL_TP_measurement_same_year":len(stw),
        "n_total_candidate_source_SA_site_pairs_near_bio":int((dist<=SEARCH_RADIUS_M).sum()),
        "nearest_measured_TP_SA_point_distance_median_m":float(np.median(dist.min(axis=1))),
        "threshold_source_opportunity":{
            f"{th}m":{
              "bio_sites_with_at_least_one_STW_SA_measured_TP_same_year_within":int((dist<=th).any(axis=1).sum()),
              "STW_SA_measured_TP_same_year_with_a_bio_site_nearby":int((dist<=th).any(axis=0).sum()),
              "undirected_STW_bio_candidate_pairs":int((dist<=th).sum()),
              "still_not_downstream_or_facility_effluent_dose":True,
            } for th in SCREENS_M
        },
        "source_concentration_chemistry_not_biological_exposure":True,
    }
    for min_months in (2,):
        valid=[s for s,v in obs[year].items() if len(v["valid_months"])>=min_months]
        if valid:
            sub=np.asarray([national[k] for k in sorted(valid)],float)
            dm=np.linalg.norm(a[:,None,:]-sub[None,:,:],axis=2)
            out[f"same_year_STW_with_at_least_{min_months}_sample_months"]={
                "n_STW_measured_TP":len(valid),
                "bio_near_500m":int((dm<=500).any(axis=1).sum()),
                "bio_near_1000m":int((dm<=1000).any(axis=1).sum()),
                "bio_near_2000m":int((dm<=2000).any(axis=1).sum()),
                "not_downstream":True
            }
    return out


def make(args):
    report={
        "schema":"GESUI_20261120_NATIONAL_WQE_2017_18_ACTUAL_SA_TP_CHEMPOP_YEAR_MATCHED_SOURCE_V1",
        "source_only":True,"date":"2026-11-20",
        "status":"SOURCE_STOP",
        "current_registration_does_not_mean_measured_2017_18_TP":True,
        "actual_measurements_only_official_0348_total_phosphorus_as_P_mgL":True,
        "year_specific_WaterCompany_final_treated_not_storm_or_river":True,
        "geography_still_undirected_BNG_NOT_river_transport":True,
        "no_annual_STW_P_mass_kg_year":True,
        "no_biological_outcome_or_effect_test_opened":True,
        "source_candidate_radius_m":SEARCH_RADIUS_M,
    }
    try:
        registry,registryqa=national.collect()
        bio,years,bioqa=local.chempop_geo()
        candidates,candidateqa=source_candidates(registry,years,bio)
        if not candidates or len(candidates)>3000:
            raise RuntimeError("National geometry source candidate set absent or unexpectedly large")
        measured,obsqa=candidate_history(candidates)
        report.update({
            "national_current_SA_registry_QA":registryqa,
            "original_ChemPop_source_year_QA":bioqa,
            "pre_original_outcome_geographic_candidate_QA":candidateqa,
            "original_EA_WQE_actual_2017_18_TP_source_QA":obsqa,
            "verified_same_year_TP_source_to_BIO_geography":{
                str(yr):compare_year(registry,measured,yr,bio,years)
                for yr in YEARS
            },
            "source_gates":{
                "original_current_registry_5722_SA_and_BNG_verified":True,
                "original_2017_2018_STW_TP_measured_sources_separated_from_unobserved_current_registry":True,
                "official_same_year_source_by_native_WQE_pointNotation_verified":True,
                "directed_river_network_outfall_to_ChemPop_bio_site_verified":False,
                "measured_STW_P_kg_year_or_annual_flow_verified":False,
                "cannot_promote_same_year_geometric_proximity_to_actual_BMWP_sewage_effect":True,
            },
        })
        report["status"]="SOURCE_YEAR_MATCHED_SA_TP_BIO_PROXIMITY_VERIFIED_NOT_ROUTED"
    except Exception as ex:
        report.update({"status":"SOURCE_REGISTRY_OR_SAMPLE_YEAR_GATE_STOP",
                       "source_error":type(ex).__name__+": "+str(ex)[:750]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+"\n")
    summary={
        "status":report["status"],"source_error":report.get("source_error"),
        "candidate_source_QA":report.get("pre_original_outcome_geographic_candidate_QA"),
        "measured_STW_source_QA":{
            k:v for k,v in report.get("original_EA_WQE_actual_2017_18_TP_source_QA",{}).items()
            if k not in ("source_page_SHA256",)
        },
        "year_matched":report.get("verified_same_year_TP_source_to_BIO_geography"),
        "source_gates":report.get("source_gates")
    }
    print("GESUI_NATIONAL_SA_ACTUAL_2017_18_TP_BIO_GEO="+json.dumps(
        summary,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
