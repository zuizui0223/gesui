#!/usr/bin/env python3
"""Independent NumPy (no statsmodels) source + country-cluster model check.

Uses only already-archived GBIF indexed recording totals at 1055 original ports.
Reproduces all cross-kingdom WWp source coefficients and 95% country-cluster
CIs, independent of patsy/statsmodels in primary model. Does not requery GBIF.
"""
from __future__ import annotations
import argparse,hashlib,json,math,urllib.parse
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

D=Path("results/global_all_Plantae_WWp2015_20261110")
S=D/"ORIGINAL_ALL1055_42COUNTRY_WWp2015_BEFORE_NEW808_GBIF_SOURCE.json"
P=D/"ORIGINAL_ALL1055_42COUNTRY_GBIF_PLANTAE_2011_24_50RECORD_RAREFACTION_SOURCE.json"
C=D/"GLOBAL1055_WWp2015_GBIF_ANIMALIA_FUNGI_COUNT_SOURCE_BEFORE_EFFECT.json"
M=D/"GLOBAL1055_WWp_PLANTAE_ANIMALIA_FUNGI_RECORDING_COUNTRYFE_AND_PLANT_RICHNESS_QC.json"
REF=1.959963984540054

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def numeric_design(sample,classes):
    countries=sorted({x["country"] for x in sample})
    X=[]
    for x in sample:
        o=x["port"];a=x["Animalia"];f=x["Fungi"]
        lat=float(o["original_port_lat"]);lon=float(o["original_port_lon"])
        X.append([
          1.,float(o["source_WWp_log1p"]),
          math.log1p(float(o["source_WWTP_count_25km"])),
          float(o["source_ESA2021_built_fraction"]),
          float(o["source_ESA2021_crop_fraction"]),
          abs(lat),math.sin(math.radians(lon)),math.cos(math.radians(lon)),
          *[float(x["country"]==c) for c in countries[1:]],
          *([math.log1p(a),math.log1p(f)] if classes else [])
        ])
    return np.asarray(X,dtype=float),countries

def targets(x):
    p=x["plant"];a=x["Animalia"];f=x["Fungi"]
    counts=p["indexed_GBIF_all_Plantae_2011_2024_records"]
    rare=p["original_exact_hypergeometric_50_plant_record_rarefied_speciesKey"]
    richness=p["faceted_original_distinct_Plantae_speciesKeys"]
    return {
        "Plantae_record_count":math.log1p(counts),
        "Animalia_record_count":math.log1p(a),
        "Fungi_record_count":math.log1p(f),
        "Combined_Animalia_Fungi_record_count":math.log1p(a+f),
        "Plantae_vs_Animalia_log_record_ratio":math.log1p(counts)-math.log1p(a),
        "Plantae_vs_Fungi_log_record_ratio":math.log1p(counts)-math.log1p(f),
        "Plantae_raw_speciesKey":math.log1p(richness),
        "Plantae_rarefied_50record":math.log1p(rare),
        "Plantae_raw_speciesKey_plus_Animalia_Fungi_records":math.log1p(richness),
        "Plantae_rarefied_50record_plus_Animalia_Fungi_records":math.log1p(rare),
        "Plantae_GBIF_record_count_plus_Animalia_Fungi_records":math.log1p(counts)
    }

def fit(sample,name,expected):
    extra="plus_Animalia_Fungi_records" in name
    X,countries=numeric_design(sample,extra)
    Y=np.asarray([targets(x)[name] for x in sample],dtype=float)
    coef=np.linalg.lstsq(X,Y,rcond=None)[0]
    if np.linalg.matrix_rank(X)!=X.shape[1]:raise ValueError("Original source model has rank deficient WPI×country design")
    residual=Y-X@coef
    B=np.linalg.inv(X.T@X)
    meat=np.zeros((X.shape[1],X.shape[1]))
    for country in countries:
        idx=np.asarray([s["country"]==country for s in sample])
        v=X[idx].T@residual[idx]
        meat+=np.outer(v,v)
    n=len(sample);k=X.shape[1];g=len(countries)
    factor=(g/(g-1))*((n-1)/(n-k))
    se=math.sqrt(float((B@meat@B*factor)[1,1]))
    beta=float(coef[1])
    ci=[beta-REF*se,beta+REF*se]
    if abs(beta-expected["beta_2015_modelled_WWp"])>1e-8:
        raise ValueError("Independent WWp regression coefficient mismatch: "+name)
    for a,b in zip(ci,expected["country_cluster95"]):
        if abs(a-b)>1e-6:
            raise ValueError("Independent country-cluster 95CI failed: "+name)
    if expected["source_port_count"]!=n or expected["source_country_count"]!=g:
        raise ValueError("Independent regression fitted different original source sites")
    return {"n_original_sites":n,"n_original_countries":g,
            "independent_beta_WWp":beta,"independent_country_cluster95":ci,
            "independent_GBIF_source_regressors":k,
            "same_country_cluster_sandwich_but_independent_linear_algebra":True}

def run(root):
    files={"WPI_WWp":root/S,"original_plant":root/P,"new_anim_fungi":root/C,"primary_model":root/M}
    data={k:json.loads(f.read_text()) for k,f in files.items()}
    macro=data["primary_model"]
    original=data["WPI_WWp"]
    plant=data["original_plant"]
    counts=data["new_anim_fungi"]
    assert macro["source_files_SHA256"]["WPI_1055_WWp_source"]==sha(root/S)
    assert macro["source_files_SHA256"]["original_GBIF_1055_Plantae"]==sha(root/P)
    assert macro["source_files_SHA256"]["new_Animalia_Fungi_count_source"]==sha(root/C)
    assert counts["precommitted_SOURCE_QC_before_effect"]["cross_kingdom_recording_control_statistical_models_PERMITTED"]
    src={x["original_port_id"]:x for x in original["source_all_1055_public_country_port_site_metadata_no_new_GBIF_response"]}
    flora={x["original_port_id"]:x for x in plant["original_1055_port_all_Plantae_source_GBIF_facet_summaries"]}
    taxa={x["original_port_id"]:x for x in counts["all_original_source_port_animal_and_fungi_indexed_counts"]}
    if not(len(src)==len(flora)==len(taxa)==1055) or set(src)!=set(flora) or set(src)!=set(taxa):
        raise ValueError("Original WPI port IDs not identical for kingdom comparisons")
    complete=Counter()
    sample=[]
    for pid in sorted(src):
        port=src[pid];q=taxa[pid];p=flora[pid]
        if port["original_source_country"]!=q["original_source_country"]:
            raise ValueError("GBIF source port changed countries")
        for label,gbifkey in (("Animalia",1),("Fungi",5)):
            r=q[label]
            if r["GBIF_status"]=="COMPLETE":
                count=r["original_indexed_record_count"]
                if type(count) is not int or count<0:raise ValueError("Invalid GBIF count")
                url=urllib.parse.urlparse(r["request_url"])
                params=urllib.parse.parse_qs(url.query)
                if url.netloc!="api.gbif.org" or url.path!="/v1/occurrence/search":
                    raise ValueError("Nonofficial GBIF count response")
                expected={
                    "taxonKey":str(gbifkey),"year":"2011,2024",
                    "geoDistance":f"{port['original_port_lat']:.7f},{port['original_port_lon']:.7f},5km",
                    "hasCoordinate":"true","hasGeospatialIssue":"false",
                    "occurrenceStatus":"PRESENT","limit":"0"
                }
                if any(params.get(k)!=[v] for k,v in expected.items()):
                    raise ValueError("Comparator query did not match original Plantae source footprint")
                if len(r["response_sha256"])!=64:raise ValueError("Original recorded count response SHA missing")
                complete[label+"_API_complete"]+=1
            elif r["GBIF_status"]=="API_ERROR_NOT_ZERO":
                if r["original_indexed_record_count"] is not None:raise ValueError("Failed GBIF source treated as zero")
            else:raise ValueError("GBIF request result status unknown")
        if not p["source_quality_min50_faceted_records_min10_species_coverage70pct"]:continue
        complete["previous_Plant752"]+=1
        if any(q[t]["GBIF_status"]!="COMPLETE" for t in ("Animalia","Fungi")):continue
        complete["both_controls_complete_in_original_plant"]+=1
        sample.append({"pid":pid,"port":port,"plant":p,"country":port["original_source_country"],
                       "Animalia":q["Animalia"]["original_indexed_record_count"],
                       "Fungi":q["Fungi"]["original_indexed_record_count"]})
    assert complete["previous_Plant752"]==752
    allowed={c for c in {s["country"] for s in sample} if len([x for x in sample if x["country"]==c])>=2 and len({x["port"]["source_WWp_log1p"] for x in sample if x["country"]==c})>=2}
    sample=[s for s in sample if s["country"] in allowed]
    assert len(sample)==macro["source_case_support"]["n_sites"]
    assert len(allowed)==macro["source_case_support"]["n_original_countries"]
    results={name:fit(sample,name,estimate) for name,estimate in macro["same_source_2015_WWp_countryFE_crosskingdom_recording_models"].items()}
    if len(results)!=11:raise ValueError("Expected original 11 source-locked response/control models")
    return {
        "schema":"gesui_global1055_original_WWp_Plantae_Animalia_Fungi_count_independent_statistical_QA_v01",
        "source_sha256":{k:sha(v) for k,v in files.items()},
        "GBIF_kingdom_count_source_query_URL_original_site_year_radius_verified":True,
        "source_site_population_1055_verified_without_posthoc_selection":True,
        "source_original_gbif_plant_752_support_verified":True,
        "source_and_API_completeness":dict(complete),
        "number_of_independent_crosskingdom_models_verified":len(results),
        "independent_numpy_countryFE_OLS_and_cluster_sandwich_replication":results,
        "no_species_richness_of_animalia_fungi_measured":True,
        "no_actual_treatment_sewage_plant_vertebrate_or_fungi_causal_effect_claimed":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "source":out["source_and_API_completeness"],
        "verified_models":out["number_of_independent_crosskingdom_models_verified"],
        "independent_models":{k:{"beta":v["independent_beta_WWp"],
            "ci":v["independent_country_cluster95"]} for k,v in out["independent_numpy_countryFE_OLS_and_cluster_sandwich_replication"].items()}
    },indent=2))

if __name__=="__main__":
    main()
