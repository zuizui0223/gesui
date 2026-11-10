#!/usr/bin/env python3
"""POSTHOC Danish EBMS NEW399 habitat confounding check (2025 × 2021 WorldCover).

Original ~286 public eBMS sites and 2021 raw WorldCover500m fractions were
frozen BEFORE this code opens licensed 2025 butterfly FAMILY count aggregates.
Absolute recorded rate and 3 vs 5 family composition depend on exact
denominator, WWTP source status, environment and site inclusion, NOT inferred
causal environmental sewage effects. Six 1degree comparison regions only.
"""
from __future__ import annotations

import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

ROOT="results/ebms_denmark_2023_2024_WWTP_source_only/"
HAB=ROOT+"DK2025_NEW399_286SITE_WORLD_COVER2021_500M_NO_FAMILY_SOURCE.json"
FAMILY=ROOT+"DENMARK_2025_NEW399_SITE_LYCAENIDAE_FAMILY_POST_SOURCE_GATE_COUNTS.json"
FROZEN=ROOT+"DENMARK_2025_NEW399_SITE_EVENT_ANY_OCC_WWTP20_SOURCE_GATE_BEFORE_FAMILY.json"
YEAR=2025
CONTROL_FRACTIONS=("built","crop","grass","tree")
BASE_COUNTS=("Lycaenidae","Pieridae","Papilionidae","Nymphalidae","Hesperiidae")
ANALYSES=("logLycaenidaePerDate","logitThree","logitFive",
          "logPieridaePerDate","logNymphalidaePerDate","logHesperiidaePerDate")
B=3500
SEED=20261110

def sha(b):return hashlib.sha256(b).hexdigest()

def sign(x):return (1 if x>0 else -1 if x<0 else 0)

def read(root):
    files={"worldcover2021":root/HAB,"original_family_source":root/FAMILY,
           "original_sourceonly_before_outcome":root/FROZEN}
    raw={k:json.loads(v.read_text(encoding="utf-8")) for k,v in files.items()}
    h=raw["worldcover2021"]
    p=raw["original_family_source"]
    source=raw["original_sourceonly_before_outcome"]
    if not h["fixed_before_habitat_adjustment_quality_gate"]["habitat_adjusted_2025_butterfly_analysis_may_proceed"]:
        raise ValueError("2021 habitat source quality gate did not pass")
    if not h["no_2025_butterfly_families_or_numbers_in_this_source"]:
        raise ValueError("Original pre-adjustment habitat sample selected on butterfly outcome")
    if not p["source_family_model_sample_gate"]["original_2025_family_contrast_reporting_gate_PASS"]:
        raise ValueError("Original biological family count reporting support failed")
    if h["original_source_only_new399_sha256"]!=sha(files["original_sourceonly_before_outcome"].read_bytes()):
        raise ValueError("Original source ID cohort changed since WorldCover pixel sampling")
    if p["source_blinded_to_2025_family_frozen_PREOUTCOME_JSON_SHA256"]!=sha(files["original_sourceonly_before_outcome"].read_bytes()):
        raise ValueError("Original 2025 family recording cohort has different source manifest")
    if h["2021_worldcover_valid_source_sites"]!=286 or h["original_286_source_eligible_monitoring_site_count"]!=286:
        raise ValueError("Source site count incomplete")
    world={x["source_original_public_site_SHA256"]:x for x in h["original_site_geo_worldcover_source_counts"]}
    if len(world)!=286 or any(x["source_QC"]!="COMPLETE" for x in world.values()):
        raise ValueError("Original WorldCover 2021 site-level raw 10m quality gate changed")
    prior=[x for x in p["2025_original_399_site_accepted_family_counts"]
           if x["source_gate_2plus_any_Occurrence_dates"]]
    if len(prior)!=286 or len({x["original_public_site_SHA256"] for x in prior})!=286:
        raise ValueError("Original 2025 all-family source sample not 286 distinct sites")
    matched=[]
    for x in prior:
        w=world[x["original_public_site_SHA256"]]
        if w["source_band"]!=x["source_band"] or w["original_1degree_region"]!=x["original_geo_1deg_cell"]:
            raise ValueError("Source distance band or original 1degree location changed")
        if w["original_2025_any_published_Occurrence_linked_dates"]!=x["linked_visit_days"]:
            raise ValueError("Original recording-denominator changed in habitat merge")
        if any(not(0<=w["fractions"][field]<=1) for field in CONTROL_FRACTIONS):
            raise ValueError("Original pixel fraction is invalid")
        fam=x["family_individualCount"]
        if any(int(fam[k])!=fam[k] or fam[k]<0 for k in BASE_COUNTS):
            raise ValueError("Noninteger original accepted butterfly family counts")
        matched.append({
            "site_sha":x["original_public_site_SHA256"],"group":x["source_band"],
            "region":x["original_geo_1deg_cell"],"days":x["linked_visit_days"],
            "fractions":w["fractions"],"counts":fam,
            "WWTP_outfall_group_SHA256":w["source_near_WWTP_id_SHA256"]})
    return files,matched,raw


def response(row,name):
    c=row["counts"];L=c["Lycaenidae"];d=row["days"]
    if name in ("logLycaenidaePerDate","logPieridaePerDate","logNymphalidaePerDate","logHesperiidaePerDate"):
        f={"logLycaenidaePerDate":"Lycaenidae","logPieridaePerDate":"Pieridae",
           "logNymphalidaePerDate":"Nymphalidae","logHesperiidaePerDate":"Hesperiidae"}[name]
        return math.log((c[f]+.5)/d)
    denom3=L+c["Pieridae"]+c["Papilionidae"]
    denom5=denom3+c["Nymphalidae"]+c["Hesperiidae"]
    if name=="logitThree":
        return math.log((L+.5)/(c["Pieridae"]+c["Papilionidae"]+.5)) if denom3>0 else None
    if name=="logitFive":
        return math.log((L+.5)/(c["Pieridae"]+c["Papilionidae"]+c["Nymphalidae"]+c["Hesperiidae"]+.5)) if denom5>0 else None
    raise ValueError("Unrecognized original family outcome")


def fit(rows,name,include_land):
    valid=[r for r in rows if response(r,name) is not None]
    regions=sorted({r["region"] for r in valid})
    if len(valid)<=len(regions)+len(CONTROL_FRACTIONS)+5:return None
    if any(sum(r["group"]=="near5" for r in valid if r["region"]==k)==0 or
           sum(r["group"]=="remote10" for r in valid if r["region"]==k)==0 for k in regions):
        # A region with no supported original family record in one source class
        # cannot identify within-region butterfly contrast.
        return None
    design=[];y=[]
    for r in valid:
        v=[1.,float(r["group"]=="near5")]
        v += [float(r["region"]==region) for region in regions[1:]]
        v.append(math.log(r["days"]))
        if include_land:v += [r["fractions"][key] for key in CONTROL_FRACTIONS]
        design.append(v);y.append(response(r,name))
    x=np.asarray(design,dtype=float);val=np.asarray(y,dtype=float)
    if np.linalg.matrix_rank(x)!=x.shape[1]:return None
    beta,resid,rank,sv=np.linalg.lstsq(x,val,rcond=None)
    return {"sample_site_n":len(valid),"source_regions_with_both_classes":len(regions),
            "regressors":int(x.shape[1]),"coefficient_log_near_minus_remote":float(beta[1]),
            "exponentiated_continuity_adjusted_ratio":float(math.exp(beta[1])),
            "site_equal_not_individual_count_weighted":True}


def percentile(seq,q):
    v=sorted(seq);return v[int(math.floor((len(v)-1)*q))]


def summarize_habitat(rows):
    def stats(a,key):
        vals=sorted(r["fractions"][key] for r in a)
        return {"n":len(vals),"mean":sum(vals)/len(vals),"median":vals[(len(vals)-1)//2],
                "q25":percentile(vals,.25),"q75":percentile(vals,.75)}
    result={}
    for g in ("near5","remote10","intermediate"):
        x=[r for r in rows if r["group"]==g]
        result[g]={key:stats(x,key) for key in (*CONTROL_FRACTIONS,"grass_plus_tree_fraction_of_land","water_fraction_of_all_circle_pixels")}
    return result


def run(root):
    files,allrows,_=read(root)
    nearremote=[r for r in allrows if r["group"] in ("near5","remote10")]
    grouped=defaultdict(list)
    for row in nearremote:grouped[row["region"]].append(row)
    shared=sorted(k for k,a in grouped.items() if {r["group"] for r in a}=={"near5","remote10"})
    if len(shared)!=6:raise ValueError("Unchanged original comparison cohort must have exactly six overlapping areas")
    chosen=[r for r in nearremote if r["region"] in shared]
    if len(chosen)!=139 or Counter(r["group"] for r in chosen)!=Counter({"near5":61,"remote10":78}):
        raise ValueError("Original six shared-area 2025 sample unexpectedly changed")
    data={
        "schema":"gesui_DK2025_NEW399_2021_WorldCover500m_habitat_6region_siteEqual_logit_rates_POSTHOC_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "source_original_SHA256":{k:sha(v.read_bytes()) for k,v in files.items()},
        "site_source_286_original_worldcover2021_samples_all_complete":True,
        "2025_original_286_eligible_documented_survey_sites_by_band":dict(Counter(r["group"] for r in allrows)),
        "2021_source_pixel_habitat_distributions_on_all_original_286_sites":summarize_habitat(allrows),
        "within_1degree_comparable_geographic_regions":shared,
        "within_shared_6_original_regions_site_counts":dict(Counter(r["group"] for r in chosen)),
        "original_shared_6_regions_total_site_n":len(chosen),
        "estimate_family_source_geometry_is_static_2022_not_WWTP_nutrient_dose":True,
        "analyses_explicitly_posthoc_after_3vs5_family_sign_inversion_known":True,
        "model":"unweighted site-level OLS, response log((accepted Lycaenidae count +0.5)/(Pieridae+Papilionidae or plus Nymphalidae+Hesperiidae+0.5)) only if positive original family denominator, or log((family count+0.5)/original ANY-OCC linked survey date); fixed original country×1° region dummy, original log linked days and original 2021 500m fractions built, crop, grass, tree. Exponentiated slope on near5 vs remote10 is a descriptive transformed recording ratio, not true abundance, dose or population effect.",
        "models":{},
        "bootstrap_cluster_source_regions":6,"bootstrap_2021_habitat_draws_target":B,
        "bootstrap_seed":SEED,
        "no_independent_source_instrument_or_randomized_WWTP_effect":True,
        "year2021_worldcover_not_species_level_larval_host_measurement":True,
    }
    for outcome in ANALYSES:
        data["models"][outcome]={
            "region_plus_logSurveyDays":fit(chosen,outcome,False),
            "region_plus_logSurveyDays_plus_WorldCover2021":fit(chosen,outcome,True),
        }
        if any(v is None for v in data["models"][outcome].values()):
            raise ValueError("Site-level effect model rank or source-control support fails: "+outcome)
    # Entire 1degree blocks, not site visits or local WWTPs, are the replication unit.
    state=SEED
    boot={kind:{model:[] for model in ("region_plus_logSurveyDays","region_plus_logSurveyDays_plus_WorldCover2021")} for kind in ANALYSES}
    for _ in range(B):
        sample=[]
        for __ in shared:
            state=(1664525*state+1013904223)&0xffffffff
            region=shared[int((state/(2**32))*len(shared))]
            sample.extend(grouped[region])
        for kind in ANALYSES:
            for model,land in (("region_plus_logSurveyDays",False),
                               ("region_plus_logSurveyDays_plus_WorldCover2021",True)):
                result=fit(sample,kind,land)
                if result is not None:boot[kind][model].append(result["coefficient_log_near_minus_remote"])
    for kind in ANALYSES:
        for model in data["models"][kind]:
            vals=boot[kind][model]
            if len(vals)<B*.90:raise ValueError("More than 10 percent geographic resamples rank-deficient")
            data["models"][kind][model]["descriptive_1degree_region_cluster_bootstrap95_exp_ratio"]=[
                float(math.exp(percentile(vals,.025))),
                float(math.exp(percentile(vals,.975)))]
            data["models"][kind][model]["bootstrap_n_valid"]=len(vals)
            data["models"][kind][model]["bootstrap_fraction_log_slope_positive"]=sum(x>0 for x in vals)/len(vals)
            data["models"][kind][model]["confidence_interval_is_not_confirmatory_randomization"]=True
    data["between_original_1degree_region_habitat_class_balance"]={}
    for g in shared:
        a=grouped[g]
        out={}
        for key in CONTROL_FRACTIONS:
            near=[r["fractions"][key] for r in a if r["group"]=="near5"]
            remote=[r["fractions"][key] for r in a if r["group"]=="remote10"]
            out[key]={"near_n":len(near),"remote_n":len(remote),
                      "near_mean":sum(near)/len(near),"remote_mean":sum(remote)/len(remote),
                      "near_minus_remote_fraction":sum(near)/len(near)-sum(remote)/len(remote)}
        data["between_original_1degree_region_habitat_class_balance"][g]=out
    return data


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=run(root)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({
        "near_remote_worldcover_mean":{g:{k:v["mean"] for k,v in x.items()} for g,x in out["2021_source_pixel_habitat_distributions_on_all_original_286_sites"].items()},
        "six_1deg_regions":out["within_1degree_comparable_geographic_regions"],
        "n":out["within_shared_6_original_regions_site_counts"],
        "models":out["models"],
    },indent=2))

if __name__=="__main__":
    main()
