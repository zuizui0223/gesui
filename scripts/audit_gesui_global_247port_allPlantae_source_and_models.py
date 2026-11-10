#!/usr/bin/env python3
"""Independent no-GBIF-requery data quality audit of source-locked all Plantae.

Verify all 247 original port IDs, country/source locations, every public GBIF
facet query URL, 50-record source minimum, records vs indexed speciesKey,
source sample support gates, and any measured WWp macro model identity.
Do NOT merge candidate butterfly host plants into an all-Plantae response.
"""
from __future__ import annotations
import argparse,hashlib,json,math,urllib.parse
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

DIR=Path("results/global_all_Plantae_WWp2015_20261110")
SOURCE=DIR/"ORIGINAL_42COUNTRY_WWp_STRATIFIED_PORTS_PRE_GBIF_ALLFLORA_SOURCE.json"
PLANT=DIR/"ORIGINAL_247PORT_42COUNTRY_GBIF_ALL_PLANTAE_2011_24_RAREFACTION_SOURCE.json"
MODEL=DIR/"GLOBAL_42COUNTRY_ALL_PLANTAE_50RECORD_RAREFIED_WWp2015_COUNTRYFE_RESULT.json"

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def independently_recompute_OLS_cluster_CI(sources,gbif,model):
    """Different QR/linalg path from original statsmodels country formula."""
    source={x["original_port_id"]:x for x in sources}
    available={x["original_port_id"]:x for x in gbif if x["source_quality_min50_faceted_records_min10_species_coverage70pct"]}
    groups=defaultdict(list)
    for pid in available:groups[source[pid]["original_source_country"]].append(pid)
    countries=sorted(k for k,ids in groups.items()
                     if len(ids)>=2 and len({source[v]["source_WWp_log1p"] for v in ids})>=2)
    rows=sorted([pid for pid in available if source[pid]["original_source_country"] in countries],
                key=lambda pid:(source[pid]["original_source_country"],pid))
    expected=model["source_per_country_WWp_variation_model_sample_n"]
    if len(rows)!=expected["port_sites"] or len(countries)!=expected["countries"]:
        raise ValueError("Independent rarefied source site case list differs")
    if rows and (len(rows),len(countries))!=(161,34):
        raise ValueError("Unexpected macro first pilot size after source QA")
    base=[];targets=defaultdict(list)
    for pid in rows:
        p=source[pid];z=available[pid]
        country=p["original_source_country"]
        lat=float(p["original_port_lat"]);lon=float(p["original_port_lon"])
        Xrow=[1.,float(p["source_WWp_log1p"]),
              math.log1p(float(p["source_WWTP_count_25km"])),
              float(p["source_ESA2021_built_fraction"]),
              float(p["source_ESA2021_crop_fraction"]),
              abs(lat),math.sin(math.radians(lon)),math.cos(math.radians(lon))]
        Xrow.extend(float(country==g) for g in countries[1:])
        base.append(Xrow)
        targets["log1p_50rarefied_allPlantae_main"].append(
            math.log1p(float(z["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"])))
        targets["log1p_raw_Plantae_speciesKey_main"].append(
            math.log1p(float(z["faceted_original_distinct_Plantae_speciesKeys"])))
        targets["log1p_Plantae_indexed_records_main"].append(
            math.log1p(float(z["indexed_GBIF_all_Plantae_2011_2024_records"])))
        targets["log1p_50rarefied_allPlantae_plus_reporting_records_control"].append(
            math.log1p(float(z["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"])))
    groups_by_obs=[source[pid]["original_source_country"] for pid in rows]
    out={}
    for name,y in targets.items():
        X=np.asarray(base,dtype=float)
        if name.endswith("plus_reporting_records_control"):
            X=np.column_stack((X,np.asarray(targets["log1p_Plantae_indexed_records_main"])))
        yy=np.asarray(y)
        coeff=np.linalg.lstsq(X,yy,rcond=None)[0]
        if np.linalg.matrix_rank(X)!=X.shape[1]:raise ValueError("Independent country-FE source design rank deficient")
        residual=yy-X@coeff
        bread=np.linalg.inv(X.T@X)
        meat=np.zeros((X.shape[1],X.shape[1]))
        G=len(countries);N=len(rows);K=X.shape[1]
        for c in countries:
            idx=np.asarray([v==c for v in groups_by_obs])
            influence=X[idx].T@residual[idx]
            meat+=np.outer(influence,influence)
        correction=(G/(G-1))*((N-1)/(N-K))
        cv=bread@meat@bread*correction
        se=math.sqrt(float(cv[1,1]))
        beta=float(coeff[1])
        ci=[beta-1.959963984540054*se,beta+1.959963984540054*se]
        ref=model["macro_country_models"][name]
        if abs(beta-ref["beta_log_WWp"])>1e-8:
            raise ValueError("Independent raw numpy OLS coefficient does not match original country FE: "+name)
        for actual,target in zip(ci,ref["country_cluster95"]):
            if abs(actual-target)>1e-5:
                raise ValueError("Independent country-sandwich 95% CI differs: "+name)
        out[name]={"independent_beta":beta,"independent_cluster95":ci,
                   "independent_source_sites":N,"independent_source_countries":G,
                   "nominal_cluster_SE":se,"independent_country_FE_regressors":K}
    return out

def audit(root):
    source=json.loads((root/SOURCE).read_text())
    gbif=json.loads((root/PLANT).read_text())
    model=json.loads((root/MODEL).read_text())
    if gbif["source_only_study_manifest_sha256"]!=sha(root/SOURCE):
        raise ValueError("Original source manifest changed after GBIF data query")
    if model["original_source_port_WWp_selection_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Model not aligned to original source cohort")
    if model["independent_GBIF_allPlantae_2011_24_247port_facet_summary_SHA256"]!=sha(root/PLANT):
        raise ValueError("Model not aligned to new all-Plantae recording opportunity")
    src=source["all_pilot_port_georeferenced_original_sites"]
    observations=gbif["original_site_GBIF_facet_summary"]
    if len(src)!=247 or len(observations)!=247:raise ValueError("Original globally fixed 247 port count mismatch")
    a={x["original_port_id"]:x for x in src}
    b={x["original_port_id"]:x for x in observations}
    if len(a)!=247 or len(b)!=247 or set(a)!=set(b):
        raise ValueError("Changing original country/port cohort")
    good=Counter()
    bad=Counter()
    countries=defaultdict(set)
    positive=Counter()
    paged=0
    for siteid,z in b.items():
        s=a[siteid]
        if z["source_country"]!=s["original_source_country"]:
            raise ValueError("Changed original source country")
        if z["GBIF_status"]=="COMPLETE":
            n=int(z["indexed_GBIF_all_Plantae_2011_2024_records"])
            k=int(z["faceted_original_speciesKey_records"])
            richness=int(z["faceted_original_distinct_Plantae_speciesKeys"])
            coverage=z["facet_record_coverage"]
            if not(0<=k<=n) or richness>k or richness<0:
                raise ValueError("Indexed all-Plantae source count does not respect records/species")
            if (coverage is None and n>0) or (coverage is not None and abs(coverage-(k/n if n else 0))>1e-9):
                raise ValueError("Original speciesKey facet coverage not reconciled")
            pages=z["original_GBIF_query_and_response_SHA_pages"]
            if len(pages)!=z["fully_paginated_original_GBIF_species_key_pages_n"] or len(pages)>50:
                raise ValueError("GBIF faceted speciesKey pagination corrupted")
            if any(len(p["response_sha256"])!=64 for p in pages):
                raise ValueError("Original faceted GBIF per-page SHA missing")
            if richness>len(pages)*1000:
                raise ValueError("Original speciesKey facet page count truncated")
            for idx,rec in enumerate(pages):
                u=urllib.parse.urlparse(rec["query_url"])
                pars=urllib.parse.parse_qs(u.query)
                checks={
                    "taxonKey":"6",
                    "year":"2011,2024",
                    "hasCoordinate":"true",
                    "hasGeospatialIssue":"false",
                    "occurrenceStatus":"PRESENT",
                    "facet":"speciesKey",
                    "facetLimit":"1000",
                    "limit":"0",
                    "facetOffset":str(1000*idx)
                }
                if u.scheme!="https" or u.netloc!="api.gbif.org" or u.path!="/v1/occurrence/search":
                    raise ValueError("Original GBIF facet source URL not official GBIF")
                if any(pars.get(key,[None])!=[v] for key,v in checks.items()):
                    raise ValueError("Original all-Plantae source GBIF query differs from precommitted parameters")
                coord=pars.get("geoDistance",[""])[0]
                if coord!=f"{s['original_port_lat']:.7f},{s['original_port_lon']:.7f},5km":
                    raise ValueError("Original all-Plantae request coordinates do not match frozen port")
                if rec["source_count"]!=n or rec["facetOffset"]!=1000*idx:
                    raise ValueError("GBIF query total count changed across faceted pages")
            isgood=(k>=50 and richness>=10 and coverage is not None and coverage>=.7 and n>=50)
            if bool(isgood)!=z["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
                raise ValueError("Predeclared 50record/10species/70percent observation gate changed")
            if isgood:
                S=z["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"]
                if S is None or not(1<=S<=min(50,richness)):
                    raise ValueError("Exact rarefied Plantae species richness outside allowed range")
                good[s["original_source_country"]]+=1
                countries[s["original_source_country"]].add(s["source_WWp_log1p"])
                positive["source_eligible_ports"]+=1
            else:bad["low_original_2011_24_GBIF_Plantae_species_records_or_coverage"]+=1
            paged+=len(pages)
        else:
            bad["failed_new_GBIF_api_requests_NOT_zeros"]+=1
            if z["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
                raise ValueError("Failed new GBIF source called an eligible all-plants port")
    total=positive["source_eligible_ports"]
    possible=(total>=120 and len(good)>=25 and
              sum(a>=2 and len(countries[k])>=2 for k,a in good.items())>=20)
    if possible!=gbif["original_source_coverage_gate"]["ALL_PLANTAE_WWp_statistical_comparison_ALLOWED"]:
        raise ValueError("Original global sample/GBIF source gate decisions are not reproducible")
    if total!=gbif["original_GBIF_result_completion_summary"]["source_qualifying_all_plantae_ports"]:
        raise ValueError("Original source accepted study-site count changed")
    if model["model_possible"]!=possible:
        raise ValueError("Reported macro model disregards original source-only 247port quality gate")
    independent_coefficients=independently_recompute_OLS_cluster_CI(src,observations,model) if possible else None
    if possible:
        m=model["macro_country_models"]
        if m is None or len(m)!=4:raise ValueError("Four original agreed all-Plantae source models missing")
        old=model["source_per_country_WWp_variation_model_sample_n"]
        if old["countries"]<20 or old["port_sites"]>total:
            raise ValueError("Country source-supported model cohort inconsistent")
        for outcome,estimate in m.items():
            if estimate["original_new_ALL_Plantae_ports_n"]!=old["port_sites"] or estimate["original_new_ALL_Plantae_countries_n"]!=old["countries"]:
                raise ValueError("Individual source-model cases mismatched")
            if len(estimate["country_cluster95"])!=2:
                raise ValueError("Source model lacks clustered country confidence interval")
    elif model["macro_country_models"] is not None:
        raise ValueError("Model results reported after original source-effort gate failure")
    return {
        "schema":"gesui_42country_247port_allPlantae_WWp2015_independent_GBIF_queries_and_model_QA_v01",
        "source_sha256":{"original42country_port_freeze":sha(root/SOURCE),
                         "original247port_2011_24_GBIF_allPlantae":sha(root/PLANT),
                         "countryFE_50rarefied_allPlantae_model":sha(root/MODEL)},
        "exact247_original_port_and_all42_source_countries_preserved":True,
        "GBIF_all_Plantae_faceted_speciesKey_query_options_verified_every_page":True,
        "GBIF_all_Plantae_speciesKey_50record_rarefaction_exact_values_valid_range":True,
        "GBIF_original_facet_query_response_SHA_pages_n":paged,
        "GBIF_source_supported_pilot_sites":total,
        "original_GBIF_source_supported_countries":len(good),
        "GBIF_query_or_opportunity_missingness":dict(bad),
        "original_all_Plantae_source_gate_evaluated_independently":possible,
        "all_macro_models_source_case_count_reconciled":True,
        "four_new_ALL_Plantae_model_source_coefficients_and_country_cluster_CI_independently_reproduced":independent_coefficients,
        "not_proof_of_2015_direct_wastewater_discharge_or_NP_uptake":True,
        "new_GBIF_Plantae_2011_24_dataset_NOT_candidate_host_only":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    v=audit(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(v,sort_keys=True,indent=2)+"\n")
    print(json.dumps(v,indent=2))
if __name__=="__main__":
    main()
