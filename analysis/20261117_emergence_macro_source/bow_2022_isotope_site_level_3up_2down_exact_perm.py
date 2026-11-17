#!/usr/bin/env python3
"""Bow River 2022 genuine aquatic→riparian isotopes: SOURCE-GATED SITE-LEVEL audit.

One Bow River, FIVE sites, THREE upstream unexposed vs TWO downstream exposed
according to published author metadata, NOT 5 independent WWTPs.
The 184 specimens are biological SUBSAMPLES; inference n=5 sites and one
river, with 10 exact permutations of two exposed-site labels among five.
Do not call isotope enrichment a measured effluent nutrient flux or spider
diet mixing proportion. Do not treat 2023 lab-analysis Date as field visit.
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math,statistics,subprocess
from collections import Counter,defaultdict
from pathlib import Path

SOURCE="https://github.com/Karen-Kidd-Lab-Ecotoxicology/Bow-River-Microbiome-2022"
COMMIT="32b36abfcf31117ac1bd9d7f3fdbc17cdcedd795"
ISO_SHA="e9ca82f50bc6f42a43213df8a0daa6d118b01cbaf16953e9e084d603e3097076"
META_SHA="f0a77f4efb51e0f0f91b87fbabb8ffcc6ef43685bbceaeec29da5879023ba5e8"
ISO_FILE="stable_isotope_analysis/StableIsotope2022.csv"
META_FILE="microbiome_analysis/Raw data/metadata_noKK2078_samples_ordered.csv"
GROUPS=("Aquatic Larvae","Aquatic Adult","Spider","Basal Food Source","Terrestrial Prey")
N=5
EXPECTED_ISO_GROUPS={"Aquatic Larvae":53,"Aquatic Adult":55,
 "Spider":30,"Basal Food Source":30,"Terrestrial Prey":16}
UPSTREAM=("Cochrane","Cushing Bridge","Sunalta")
DOWNSTREAM=("Graves Bridge","Policeman Flats")
SCHEMA="GESUI_20261117_BOW_2022_EXACT_AUTHOR_SITE_SOURCE_ISOTOPE_PAIRED_BIOTIC_V1"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
def csv_source(root,path,expected):
    p=root/path
    if sha(p)!=expected:raise RuntimeError("Original author raw CSV changed: "+path)
    with p.open(encoding="utf-8-sig",newline="") as f:
        return list(csv.DictReader(f))
def parse(x):
    try:
        n=float(str(x).strip())
        return n if math.isfinite(n) else None
    except (ValueError,TypeError):return None
def site_from_original_metadata(rows):
    if len(rows)!=380:raise RuntimeError("Original microbiome metadata 380 source count changed")
    d=defaultdict(lambda:defaultdict(set))
    month=defaultdict(set)
    counts=Counter()
    for r in rows:
        if r.get("Stream_Type")!="Bow River":continue
        if r.get("Sample_Type")!="Sample":continue
        site=str(r.get("Site","")).strip()
        if not site:continue
        d[site]["exposure"].add(r.get("Exposure"))
        d[site]["stream_type"].add(r.get("Stream_Type"))
        month[site].add(r.get("Month"))
        counts[site]+=1
    if set(d)!=set(UPSTREAM)|set(DOWNSTREAM):
        raise RuntimeError("Field observation five original Bow sites changed")
    for site in d:
        status=d[site]["exposure"]
        expected={"Unexposed"} if site in UPSTREAM else {"Exposed"}
        if status!=expected or d[site]["stream_type"]!={"Bow River"}:
            raise RuntimeError("Published upstream/downstream treatment status drift at "+site)
        if month[site]!={"July"}:
            raise RuntimeError("Author 2022 Bow collection season differs at "+site)
    return {k:{"source_original_exposure":next(iter(v["exposure"])),
                "original_2022_field_collection_month":"July",
                "original_microbiome_sample_rows":counts[k]}
            for k,v in sorted(d.items())}
def original_isotopes(rows):
    if len(rows)!=184:raise RuntimeError("Original stable isotope 184 specimens changed")
    counts=Counter(r["sample_group"] for r in rows)
    if dict(counts)!=EXPECTED_ISO_GROUPS:raise RuntimeError("Original five source stage sample sizes changed")
    if Counter(r["Site"] for r in rows)!={"Graves Bridge":45,"Policeman Flats":41,
          "Sunalta":40,"Cochrane":29,"Cushing Bridge":29}:
        raise RuntimeError("Original five site sample sizes changed")
    return rows
def med(vals):
    vals=[x for x in vals if x is not None]
    return float(statistics.median(vals)) if vals else None
def stat(rows,site,group,taxon=None):
    r=[x for x in rows if x["Site"]==site and x["sample_group"]==group and
       (taxon is None or x.get("Taxa")==taxon)]
    return {"specimen_n":len(r),
      "valid_d15N_n":sum(parse(x.get("d15N")) is not None for x in r),
      "valid_d13C_n":sum(parse(x.get("d13C")) is not None for x in r),
      "median_d15N":med(parse(x.get("d15N")) for x in r),
      "median_d13C":med(parse(x.get("d13C")) for x in r)}
def difference(rows,exp,n_exp=2):
    keys=[r["site"] for r in rows]
    if len(set(keys))!=len(keys):raise RuntimeError("Site was treated as repeat")
    y={r["site"]:r["value"] for r in rows}
    source={r["site"]:r["exposure"] for r in rows}
    good=sorted(y)
    if any(y[k] is None for k in good):
        return {"status":"NO_COMPLETE_SITE_PAIRING_NEEDED_FOR_EXACT_TEST",
                "sites_with_measured_y":sum(y[k] is not None for k in good)}
    observed=statistics.mean(y[k] for k in good if source[k]=="Exposed")-\
             statistics.mean(y[k] for k in good if source[k]=="Unexposed")
    if sum(source[k]=="Exposed" for k in good)!=n_exp:
        raise RuntimeError("Unexpected original published exposed-site count")
    perm=[]
    for ex in itertools.combinations(good,n_exp):
        treated=set(ex);untreated=set(good)-treated
        val=statistics.mean(y[k] for k in treated)-statistics.mean(y[k] for k in untreated)
        perm.append(val)
    p_two=sum(abs(x)>=abs(observed)-1e-12 for x in perm)/len(perm)
    p_pos=sum(x>=observed-1e-12 for x in perm)/len(perm)
    return {"status":"EXACT_SOURCE_SITE_LABEL_PERMUTATION_DESCRIPTIVE_ONLY",
      "n_independent_source_rivers":1,"n_observation_sites":len(good),
      "n_exposed_downstream_sites":n_exp,
      "n_unexposed_upstream_sites":len(good)-n_exp,
      "downstream_mean_site_median_minus_upstream_mean_site_median":float(observed),
      "mean_upstream_site_medians":float(statistics.mean(y[k] for k in good if source[k]=="Unexposed")),
      "mean_downstream_site_medians":float(statistics.mean(y[k] for k in good if source[k]=="Exposed")),
      "exact_permutation_assignments":len(perm),
      "exact_twosided_label_permutation_p":float(p_two),
      "exact_positive_tail_label_permutation_p":float(p_pos),
      "min_attainable_one_sided_p":1/len(perm),
      "permutations_do_not_remove_downstream_spatial_confounding":True}
def run(a):
    root=a.original.resolve()
    remote=subprocess.check_output(["git","-C",str(root),"remote","get-url","origin"],text=True).strip()
    actual=subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"],text=True).strip()
    if actual!=COMMIT or SOURCE.lower() not in remote.lower():
        raise RuntimeError("Original published repo commit/remote not source locked")
    orig=original_isotopes(csv_source(root,ISO_FILE,ISO_SHA))
    meta=site_from_original_metadata(csv_source(root,META_FILE,META_SHA))
    sites=sorted(meta)
    allsite=[]
    for site in sites:
        raw={"site":site,"source_exposure":meta[site]["source_original_exposure"],
             "source_original_microbiome_n":meta[site]["original_microbiome_sample_rows"]}
        for group in GROUPS:
            raw[group]=stat(orig,site,group)
        for taxon in ("Biofilm","Terrestrial Leaves","Araneidae","Tetragnathidae"):
            group="Basal Food Source" if taxon in ("Biofilm","Terrestrial Leaves") else "Spider"
            raw[f"Taxon_{taxon}"]=stat(orig,site,group,taxon)
        # Isotopic offsets are baseline-adjusted source indicators, not nutrient
        # transfer fractions.  Basal biofilm is n=3 original rows at each site.
        bio=raw["Taxon_Biofilm"]
        for group in ("Aquatic Larvae","Aquatic Adult","Spider"):
            for k in ("d15N","d13C"):
                target=raw[group]["median_"+k];baseline=bio["median_"+k]
                raw[f"{group}_minus_biofilm_{k}"]=(
                    float(target-baseline) if target is not None and baseline is not None else None)
        allsite.append(raw)
    for site in allsite:
        if site["Spider"]["specimen_n"]!=6 or site["Taxon_Biofilm"]["specimen_n"]!=3:
            raise RuntimeError("No valid six-spider and three-biofilm per original site")
        if any(site[q]["valid_d15N_n"]<3 or site[q]["valid_d13C_n"]<3
               for q in ("Aquatic Larvae","Aquatic Adult","Spider","Taxon_Biofilm")):
            raise RuntimeError("Insufficient original carbon/nitrogen isotope support at one site")
    comparisons=[]
    patterns=[
      *[(g,k) for g in ("Aquatic Larvae","Aquatic Adult","Spider","Taxon_Biofilm",
                         "Taxon_Terrestrial Leaves","Taxon_Araneidae",
                         "Taxon_Tetragnathidae","Terrestrial Prey")
                     for k in ("median_d15N","median_d13C")],
      *[(None,k) for k in ("Aquatic Larvae_minus_biofilm_d15N",
                           "Aquatic Adult_minus_biofilm_d15N",
                           "Spider_minus_biofilm_d15N",
                           "Aquatic Larvae_minus_biofilm_d13C",
                           "Aquatic Adult_minus_biofilm_d13C",
                           "Spider_minus_biofilm_d13C")]]
    for group,k in patterns:
        values=[{"site":x["site"],"exposure":x["source_exposure"],
                 "value":x[group][k] if group else x[k]} for x in allsite]
        matched=[x for x in values if x["value"] is not None]
        item={"source_biological_group":group or "within_site_delta_vs_biofilm",
              "source_isotope_variable":k,
              "n_source_site_medians_supported":len(matched)}
        if len(matched)==5:
            item["contrast"]=difference(values)
        else:item["contrast"]={"status":"INCOMPLETE_SOURCE_SITE_SUPPORT_DO_NOT_FIT_FIVE_SITE_EXACT_LABEL_TEST"}
        comparisons.append(item)
    output={"schema":SCHEMA,
      "status":"PASS_SITE_LEVEL_ISOTOPE_SOURCE_NUMERIC_DESCRIPTIVE_NOT_WWTP_CAUSAL",
      "original_public_author_repo":SOURCE,
      "original_public_repo_commit":COMMIT,
      "original_isotope_sha256":ISO_SHA,"original_microbiome_metadata_sha256":META_SHA,
      "original_isotope_specimen_n":184,
      "original_microbiome_metadata_rows":380,
      "original_observation_river_n":1,
      "original_river_field_sites_n":5,
      "original_unexposed_sites_n":3,
      "original_exposed_sites_n":2,
      "source_group_counts":dict(EXPECTED_ISO_GROUPS),
      "original_bow_field_collection_month":"July",
      "original_isotope_lab_analysis_date_field_has_2023_values_NOT_2023_field_visits":True,
      "no_independent_five_WWTP_replicates":True,
      "site_source_exposure":meta,
      "site_level_original_iso_medians":allsite,
      "site_label_exact_permutation_comparisons":comparisons,
      "source_exposure_labels_not_randomly_assigned_treatments":True,
      "observed_isotope_enrichment_not_chemical_source_identification_without_effluent_baseline":True,
      "spider_and_emergent_adult_medians_not_true_emergence_biomass_or_consumed_flux":True,
      "all_permutation_p_descriptive_and_not_valid_independent_causal_treatment_p":True,
      "no_new_cross_site_multiriver_meta_analysis_or_causal_transfer_identified":True}
    a.out.mkdir(parents=True,exist_ok=True)
    file=a.out/"BOW_2022_GENUINE_3UP_2DOWN_SITE_MEDIAN_ISOTOPES.json"
    file.write_text(json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_BOW_REAL_ISOTOPE_ANALYSIS="+json.dumps({
        "status":output["status"],"sites":5,"iso_specimens":184,
        "source_spider_n":30,
        "group_d15N_site_exposure_differences":{
            x["source_biological_group"] or x["source_isotope_variable"]:
                round(x["contrast"].get("downstream_mean_site_median_minus_upstream_mean_site_median"),5)
            for x in comparisons
            if (x["source_isotope_variable"]=="median_d15N" and
                x["source_biological_group"] in ("Aquatic Larvae","Aquatic Adult","Spider")) or
                x["source_isotope_variable"]=="Spider_minus_biofilm_d15N"
        }},ensure_ascii=False,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--original",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
