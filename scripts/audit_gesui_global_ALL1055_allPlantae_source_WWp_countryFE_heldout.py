#!/usr/bin/env python3
"""Independent source and OLS country-sandwich audit of ALL1055 all-Plantae.

Rebuild all original 1055 site/source ID joins and rarefaction/GBIF query
metadata, original 247 unchanged and new808 counted. Recompute full countryFE
OLS beta and 95% cluster-CI via raw numerical normal-equation sandwich, not
original statsmodels/patsy implementation. All predictive squared errors
are also arithmetically rechecked (no new GBIF requests).
"""
from __future__ import annotations
import argparse,hashlib,json,math,urllib.parse
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

D=Path("results/global_all_Plantae_WWp2015_20261110")
S=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
G=D/"ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json"
M=D/"GLOBAL_ALL1055_WWp2015_ALL_PLANTAE_50RAREFIED_42COUNTRY_MACRO_AND_COUNTRY_HOLDOUT.json"
O=D/"ORIGINAL_247PORT_42COUNTRY_GBIF_ALL_PLANTAE_2011_24_RAREFACTION_SOURCE.json"
REF=1.959963984540054

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def close(x,y,where,tol=1e-7):
    if not math.isfinite(x) or not math.isfinite(y) or abs(x-y)>tol:
        raise ValueError("Mismatch at "+where+": "+str(x)+" vs "+str(y))

def sample_case(source,records):
    dd=[]
    for pid,p in source.items():
        z=records[pid]
        if not z["source_quality_min50_faceted_records_min10_species_coverage70pct"]:continue
        if z["GBIF_status"]!="COMPLETE":raise ValueError("Invalid quality-supported source record")
        dd.append((p,z))
    group=defaultdict(list)
    for p,z in dd:group[p["original_source_country"]].append((p,z))
    valid={k for k,vs in group.items() if len(vs)>=2 and len({v[0]["source_WWp_log1p"] for v in vs})>=2}
    arr=sorted([(p,z) for p,z in dd if p["original_source_country"] in valid],
               key=lambda x:(x[0]["original_source_country"],x[0]["original_port_id"]))
    return arr,sorted(valid)

def linear_source_fit(arr,cols,extra_effort=False):
    countries=sorted({p["original_source_country"] for p,z in arr})
    X=[];Y=[]
    for p,g in arr:
        lat=float(p["original_port_lat"]);lon=float(p["original_port_lon"])
        row=[1.,float(p["source_WWp_log1p"]),math.log1p(float(p["source_WWTP_count_25km"])),
             float(p["source_ESA2021_built_fraction"]),float(p["source_ESA2021_crop_fraction"]),
             abs(lat),math.sin(math.radians(lon)),math.cos(math.radians(lon))]
        row.extend(float(p["original_source_country"]==k) for k in countries[1:])
        if extra_effort:row.append(math.log1p(g["indexed_GBIF_all_Plantae_2011_2024_records"]))
        X.append(row)
        if cols=="rarefied":Y.append(math.log1p(g["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"]))
        elif cols=="richness":Y.append(math.log1p(g["faceted_original_distinct_Plantae_speciesKeys"]))
        elif cols=="records":Y.append(math.log1p(g["indexed_GBIF_all_Plantae_2011_2024_records"]))
        else:raise ValueError(cols)
    x=np.asarray(X,float)
    y=np.asarray(Y,float)
    coef=np.linalg.lstsq(x,y,rcond=None)[0]
    if np.linalg.matrix_rank(x)!=x.shape[1]:
        raise ValueError("CountryFE full all1055 original source design rank-deficient")
    resid=y-x@coef
    bread=np.linalg.inv(x.T@x)
    cluster=np.zeros((x.shape[1],x.shape[1]))
    for c in countries:
        mask=np.array([p["original_source_country"]==c for p,z in arr])
        t=x[mask].T@resid[mask]
        cluster+=np.outer(t,t)
    N=x.shape[0];K=x.shape[1];G=len(countries)
    factor=(G/(G-1))*((N-1)/(N-K))
    covariance=bread@cluster@bread*factor
    se=math.sqrt(covariance[1,1])
    beta=float(coef[1])
    return {"n_ports":N,"n_countries":G,"beta":beta,
            "country_cluster95":[beta-REF*se,beta+REF*se],
            "source_regressors":K}

def measure_holdouts(raw,expected,description):
    if len({r["port_id"] for r in raw})!=len(raw):raise ValueError("Duplicate source holdout WPI port record "+description)
    groups=defaultdict(list)
    for r in raw:
        groups[r["country"]].append(r)
        if r["baseline_error_sq"]<0 or r["plus2015WWp_error_sq"]<0:
            raise ValueError("Negative squared prediction error")
    def mse(rows):
        return [sum(x["baseline_error_sq"] for x in rows)/len(rows),
                sum(x["plus2015WWp_error_sq"] for x in rows)/len(rows)]
    a,b=mse(raw)
    ca=[];cb=[]
    for g,rows in groups.items():
        q,r=mse(rows)
        ca.append(q);cb.append(r)
    ga=sum(ca)/len(ca);gb=sum(cb)/len(cb)
    for key,v in (
        ("port_pooled_MSE_baseline",a),("port_pooled_MSE_plusWWp",b),
        ("country_equal_MSE_baseline",ga),("country_equal_MSE_plusWWp",gb),
        ("port_pooled_MSE_relative_gain_percent",100*(a-b)/a),
        ("country_equal_MSE_relative_gain_percent",100*(ga-gb)/ga)):
        close(v,expected[key],description+"."+key,tol=1e-9)
    if expected["n_ports"]!=len(raw) or expected["n_countries"]!=len(groups):
        raise ValueError("Source holdout country sample count changed")
    return {"port_count":len(raw),"countries":len(groups),
            "port_equal_relative_mse_gain_percent":100*(a-b)/a,
            "country_equal_relative_mse_gain_percent":100*(ga-gb)/ga}

def run(root):
    files={k:root/p for k,p in (("frozen1055",S),("sourceGBIF1055",G),("macro_model",M),("old247",O))}
    input={k:json.loads(f.read_text()) for k,f in files.items()}
    src=input["frozen1055"];gbif=input["sourceGBIF1055"];m=input["macro_model"]
    assert src["GBIF_new808_species_records_NOT_OBSERVED_in_stage0"] is True
    assert gbif["original_all1055_WWp_WPI_country_source_frame_SHA256"]==sha(root/S)
    assert gbif["old_original_247_WPI_GBIF_all_Plantae_source_result_SHA256"]==sha(root/O)
    assert m["original_1055_WPI_WWp2015_full_source_SHA256"]==sha(root/S)
    assert m["original_247old_plus_808new_allPlantae_GBIF_1055_result_SHA256"]==sha(root/G)
    p={x["original_port_id"]:x for x in src["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]}
    q={x["original_port_id"]:x for x in gbif["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    old={x["original_port_id"]:x for x in input["old247"]["original_site_GBIF_facet_summary"]}
    if not(len(p)==len(q)==1055 and len(old)==247 and set(q)==set(p)):
        raise ValueError("All original 1055 WPI fixed source site IDs not observed")
    for pid in old:
        if q[pid]!=old[pid]:
            raise ValueError("Old 247 source flora exactly conserved between pilot and full follow-up")
    if len(gbif["four_new_2026_11_10_GBIF_allPlantae_source_shards_SHA256"])!=4:
        raise ValueError("Source GBIF new808 not exactly four independent 202site shard artifacts")
    stats=Counter()
    countries=defaultdict(set)
    for pid,x in q.items():
        if x["source_country"]!=p[pid]["original_source_country"]:
            raise ValueError("Original source country not identical")
        if x["GBIF_status"]!="COMPLETE":
            if x["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
                raise ValueError("Failed GBIF request treated as eligible allPlantae site")
            stats["api_fail"]+=1;continue
        stats["api_success"]+=1
        n=x["indexed_GBIF_all_Plantae_2011_2024_records"]
        k=x["faceted_original_speciesKey_records"]
        s=x["faceted_original_distinct_Plantae_speciesKeys"]
        if not (n>=k>=s>=0):raise ValueError("GBIF source facet records/species keys inconsistent")
        pages=x["original_GBIF_query_and_response_SHA_pages"]
        if len(pages)!=x["fully_paginated_original_GBIF_species_key_pages_n"]:
            raise ValueError("Original GBIF page count incomplete")
        for i,r in enumerate(pages):
            url=urllib.parse.urlparse(r["query_url"])
            pars=urllib.parse.parse_qs(url.query)
            if url.netloc!="api.gbif.org" or url.path!="/v1/occurrence/search":
                raise ValueError("Site response not original official GBIF API URL")
            for k0,v0 in (("taxonKey","6"),("year","2011,2024"),("facet","speciesKey"),
                          ("facetLimit","1000"),("limit","0"),("facetOffset",str(i*1000))):
                if pars.get(k0)!=[v0]:raise ValueError("Original source query mutated "+k0)
            if pars.get("geoDistance")!=[f"{p[pid]['original_port_lat']:.7f},{p[pid]['original_port_lon']:.7f},5km"]:
                raise ValueError("Original GBIF site location changed from frozen WPI")
            if len(r["response_sha256"])!=64 or r["source_count"]!=n:
                raise ValueError("Original GBIF raw response SHA or page count missing")
        if x["source_quality_min50_faceted_records_min10_species_coverage70pct"]:
            if n<50 or k<50 or s<10 or x["facet_record_coverage"]<.7:
                raise ValueError("Original precommitted allPlantae rarefaction sample gate changed")
            val=x["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"]
            if not 1<=val<=min(50,s):raise ValueError("Invalid all1055 source 50-record rarefaction")
            stats["source_valid"]+=1
            countries[p[pid]["original_source_country"]].add(p[pid]["source_WWp_log1p"])
        else:stats["source_invalid_plant_sampling_support"]+=1
    variation=sum(len(v)>=2 for v in countries.values())
    expected=(stats["source_valid"]>=550 and len(countries)>=35 and variation>=30)
    if expected!=gbif["precommitted_full_1055_AllPlantae_response_support_GATE"]["ALL1055_AllPlantae_WWp2015_macro_regression_ALLOWED"]:
        raise ValueError("Full 1055 original GBIF plant source quality gate mismatch")
    if expected!=m["ALL_Plantae_original_macro_model_eligible"]:
        raise ValueError("Full macro model improperly fit after failed source gate")
    coeffs={}
    cv={}
    if expected:
        sample,cs=sample_case(p,q)
        baseline=m["original_all_Plantae_full1055_model_actual_sample"]
        if len(sample)!=baseline["n_original_port_records"] or len(cs)!=baseline["n_original_countries"]:
            raise ValueError("Country-with-two-ports macro cohort changed")
        mname={
            "rarefied50_Plantae":("rarefied",False),
            "raw_Plantae_speciesKey":("richness",False),
            "all_Plantae_GBIF_record_count":("records",False),
            "rarefied50_plus_original_plant_record_effort":("rarefied",True)
        }
        for key,(response,effort) in mname.items():
            value=linear_source_fit(sample,response,effort)
            oldfit=m["all1055_2015WWp_and_2011_24_all_Plantae_countryFE_models"][key]
            close(value["beta"],oldfit["beta_on_log1p_modeled_WWp2015"],"FULL_WWp_beta_"+key)
            for a,b in zip(value["country_cluster95"],oldfit["country_cluster95"]):
                close(a,b,"FULL_WWp_clusterCI_"+key,1e-6)
            coeffs[key]=value
        regional=m["same_country_2degree_spatial_block_5fold_predictions"]
        if regional.get("status")=="ESTIMATED_SPATIALLY_OUTOFBLOCK_WITHIN_ALREADY_REPRESENTED_COUNTRIES":
            errors=regional["source_test_port_errors_original"]
            cv["within_country_2deg"]=measure_holdouts(errors,regional["original_same_country_2degree_5fold_source_independent_mse"],"within_country_2deg")
            folds=defaultdict(set)
            for row in errors:
                pair=(row["country"],row["source_2degree_block"])
                folds[pair].add(row["fold"])
            if any(len(x)!=1 for x in folds.values()):
                raise ValueError("Source 2° country geography leaked across heldout folds")
        unseen=m["STRICT_unseen_country_without_country_fixed_effects"]
        if unseen.get("status")=="STRICT_UNSEEN_COUNTRY_ALL_PLANTAE_NO_COUNTRY_FE":
            cv["unseen_country"]=measure_holdouts(unseen["source_test_port_errors_original"],
                         unseen["unseen_country_modelled_WWp_pilot_forecast_gain"],"unseen_country")
            if cv["unseen_country"]["port_count"]!=len(sample):
                raise ValueError("Unseen country model failed to score every original source-qualified port")
    return {
        "schema":"gesui_2026_11_global1055_WWp2015_allPlantae_countryFE_2degree_countryHoldout_independent_arithmetic_QA_v01",
        "source_inputSHA256":{name:sha(f) for name,f in files.items()},
        "old247_sources_also_in_full1055_exactly_unchanged":True,
        "original_four_new808_GBIF_source_shards_202sites_each":True,
        "original_all1055_query_urls_years_taxon_radii_and_sha_checked":True,
        "original_full_GBIF_all_Plantae_record_and_speciesKey_QC":dict(stats),
        "quality_supported_countries":len(countries),
        "independent_source_support_gate":expected,
        "independent_four_full_global_WWp_countryFE_beta_countrycluster95":coeffs,
        "independent_spatial_country_holdout_MSE":cv,
        "new808_outcome_posthoc_after_pilot_known":True,
        "modeled2015_WWp_is_not_real_wastewater_nutrient_source_to_terrestrial_plants":True,
        "allplantae_GBIF_species_record_not_local_floristic_abundance_or_non_native_status":True,
    }

def main():
    a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True);args=a.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:x[k] for k in ("original_full_GBIF_all_Plantae_record_and_speciesKey_QC",
                 "quality_supported_countries","independent_source_support_gate",
                 "independent_four_full_global_WWp_countryFE_beta_countrycluster95",
                 "independent_spatial_country_holdout_MSE")},indent=2))
if __name__=="__main__":
    main()
