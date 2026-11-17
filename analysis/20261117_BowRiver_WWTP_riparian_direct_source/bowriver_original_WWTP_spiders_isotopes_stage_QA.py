#!/usr/bin/env python3
"""Source-pinned original Bow River 2022 wastewater-exposure isotope stage audit.

REANALYSIS of Diesbourg et al. (2025), not a new independent study. Three
reported WWTPs along one river are NOT five independent WWTP treatments.
Exposure labels belong to original published 2022 metadata, with two
exposed and three unexposed river sampling sites.

Do not fit a pseudoreplicated individual-level 184-sample treatment p-value.
Report source-authentic SITE means and within-site isotope baselines only.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, itertools, json, statistics
from collections import Counter
from pathlib import Path

UPSTREAM_REPO="Karen-Kidd-Lab-Ecotoxicology/Bow-River-Microbiome-2022"
UPSTREAM_COMMIT="32b36abfcf31117ac1bd9d7f3fdbc17cdcedd795"
ORIGINAL_GITBLOB_SHA1={
 "raw":"849ad337ecdb689e2a8265a3e3db5ef8671a506a",
 "subset":"404e3136eb0e5cd0196b9d7b205778e71280171f",
 "metadata":"638a729d8eecfca225f671e805e26d242662f983"
}
SITES=("COCH","SUN","CUSHBR","GRVBR","PMF")
EXPOSURES={"COCH":"Unexposed","SUN":"Unexposed","CUSHBR":"Unexposed",
           "GRVBR":"Exposed","PMF":"Exposed"}
TAXA=("Tetragnathidae","Araneidae","Terrestrial Leaves","Biofilm","Chironomidae","Hydropsychidae")
SCHEMA="GESUI_20261117_BOWRIVER_ORIGINAL_WWTP_SITE_SPIDER_INSECT_ISOTOPES_SOURCE_REANALYSIS_V1"

def git_blob(data:bytes)->str:
    return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
def loadcsv(path,expected):
    b=path.read_bytes()
    if git_blob(b)!=expected:raise RuntimeError("Git blob of 2025 publication source changed: "+str(path))
    return list(csv.DictReader(io.StringIO(b.decode("utf-8-sig")))),hashlib.sha256(b).hexdigest()
def num(row,col):
    x=str(row.get(col,"")).strip()
    if not x or x.upper() in ("N/A","NA"):return None
    try:return float(x)
    except ValueError:return None
def mean(rows,col):
    x=[num(r,col) for r in rows]
    values=[v for v in x if v is not None]
    return float(statistics.mean(values)) if values else None
def avg(vals):
    v=[x for x in vals if x is not None]
    return float(statistics.mean(v)) if v else None
def one_type(rows,field,label):
    q=[r for r in rows if r[field]==label]
    return {"n":len(q),"d13C_mean":mean(q,"d13C"),"d15N_mean":mean(q,"d15N")}
def minus(a,b):return None if a is None or b is None else a-b

def run(args):
    raw,h_raw=loadcsv(args.raw,ORIGINAL_GITBLOB_SHA1["raw"])
    sub,h_sub=loadcsv(args.subset,ORIGINAL_GITBLOB_SHA1["subset"])
    md,h_md=loadcsv(args.metadata,ORIGINAL_GITBLOB_SHA1["metadata"])
    if (len(raw),len(sub),len(md))!=(184,155,380):
        raise RuntimeError("Exact real original 184/155/380 source rows changed")
    if len({x["SampleName"] for x in raw})!=184 or len({x["SampleName"] for x in sub})!=155:
        raise RuntimeError("Original 2022 isotopic specimen IDs not unique")
    if not {x["SampleName"] for x in sub}.issubset({x["SampleName"] for x in raw}):
        raise RuntimeError("Original authors' excluded-subset not nested in raw data")
    if set(x["Site_Code"] for x in raw)!=set(SITES):
        raise RuntimeError("Original isotope 5-site source cohort changed")
    source_exp={}
    for x in md:
        if x["Site_code"] in SITES and x["Exposure"]:
            source_exp.setdefault(x["Site_code"],set()).add(x["Exposure"])
    if {k:list(sorted(v)) for k,v in source_exp.items()}!={k:[v] for k,v in EXPOSURES.items()}:
        raise RuntimeError("Original published WWTP exposure source status differs")
    per_site=[]
    original_site_counts={"COCH":29,"SUN":40,"CUSHBR":29,"GRVBR":45,"PMF":41}
    for site in SITES:
        r=[x for x in raw if x["Site_Code"]==site]
        if len(r)!=original_site_counts[site]:
            raise RuntimeError("Original raw isotope count drift at "+site)
        m=[x for x in md if x["Site_code"]==site and x["Exposure"]]
        taxa={t:one_type(r,"Taxa",t) for t in TAXA}
        stages={t:one_type(r,"Stage",t) for t in ("Larvae","Adult","Spider","Shrub","Biofilm")}
        if taxa["Tetragnathidae"]["n"]!=3 or taxa["Araneidae"]["n"]!=3 or taxa["Terrestrial Leaves"]["n"]!=3:
            raise RuntimeError("Original matched spider/terrestrial baseline per site lacks three specimens")
        water={}
        for p in ("Total_Phosphorus","Ammonia","Conductivity","Dissolved_oxygen","Total_Organic_Carbon"):
            v=sorted({num(x,p) for x in m if num(x,p) is not None})
            water[p]={"n_distinct_reported_original_values":len(v),"only_value":v[0] if len(v)==1 else None}
        tt=taxa["Tetragnathidae"]
        leaves=taxa["Terrestrial Leaves"]
        per_site.append({
            "site":site,"original_exposure_label":EXPOSURES[site],
            "n_raw_isotope_specimens":len(r),
            "n_author_155_subset_specimens":sum(x["Site_Code"]==site for x in sub),
            "n_original_2022_metadata_rows":len(m),
            "original_source_unique_water_chemistry":water,
            "original_d13C_d15N_by_taxon":taxa,
            "original_d13C_d15N_by_stage":stages,
            "tetragnathid_d13C_minus_terrestrial_leaves_d13C_permil":
                minus(tt["d13C_mean"],leaves["d13C_mean"]),
            "tetragnathid_d15N_minus_terrestrial_leaves_d15N_permil":
                minus(tt["d15N_mean"],leaves["d15N_mean"]),
            "tetragnathid_d13C_minus_local_aquatic_adult_mean_d13C_permil":
                minus(tt["d13C_mean"],stages["Adult"]["d13C_mean"])})
    exposed=[x for x in per_site if x["original_exposure_label"]=="Exposed"]
    unexposed=[x for x in per_site if x["original_exposure_label"]=="Unexposed"]
    if (len(exposed),len(unexposed))!=(2,3):
        raise RuntimeError("Original independently classified exposure sites not 2/3")
    comparisons={}
    properties={
      "tetragnathid_d13C_permil":lambda x:x["original_d13C_d15N_by_taxon"]["Tetragnathidae"]["d13C_mean"],
      "araneid_d13C_permil":lambda x:x["original_d13C_d15N_by_taxon"]["Araneidae"]["d13C_mean"],
      "tetragnathid_d15N_permil":lambda x:x["original_d13C_d15N_by_taxon"]["Tetragnathidae"]["d15N_mean"],
      "terrestrial_leaves_d15N_permil":lambda x:x["original_d13C_d15N_by_taxon"]["Terrestrial Leaves"]["d15N_mean"],
      "tetragnathid_d13C_minus_local_terrestrial_leaves_d13C_permil":
         lambda x:x["tetragnathid_d13C_minus_terrestrial_leaves_d13C_permil"],
      "tetragnathid_d15N_minus_local_terrestrial_leaves_d15N_permil":
         lambda x:x["tetragnathid_d15N_minus_terrestrial_leaves_d15N_permil"],
      "tetragnathid_d13C_minus_local_aquatic_adult_d13C_permil":
         lambda x:x["tetragnathid_d13C_minus_local_aquatic_adult_mean_d13C_permil"]}
    for field,f in properties.items():
        E=avg([f(x) for x in exposed])
        U=avg([f(x) for x in unexposed])
        if E is None or U is None:raise RuntimeError("Critical original isotope endpoint was missing")
        comparisons[field]={
          "exposed_two_SITE_mean":E,
          "unexposed_three_SITE_mean":U,
          "exposed_minus_unexposed_site_equal_difference":E-U,
          "leave_one_of_five_SITE_OUT_difference_range":[
             min((avg([f(q) for q in exposed if q["site"]!=d["site"]])-
                  avg([f(q) for q in unexposed if q["site"]!=d["site"]]))
                 for d in per_site if
                   [q for q in exposed if q["site"]!=d["site"]] and
                   [q for q in unexposed if q["site"]!=d["site"]]),
             max((avg([f(q) for q in exposed if q["site"]!=d["site"]])-
                  avg([f(q) for q in unexposed if q["site"]!=d["site"]]))
                 for d in per_site if
                   [q for q in exposed if q["site"]!=d["site"]] and
                   [q for q in unexposed if q["site"]!=d["site"]])],
          "no_individual_specimen_level_p_value":True,
          "not_randomized_or_causally_identified":True}
    report={
      "schema":SCHEMA,
      "status":"PASS_ORIGINAL_SITE_EQUAL_WWTP_ISOTOPE_REANALYSIS_NO_NEW_CAUSAL_INFERENCE",
      "publication":"Diesbourg et al. Environmental Pollution 2025 DOI 10.1016/j.envpol.2025.125948",
      "data_citation":"https://github.com/"+UPSTREAM_REPO,
      "original_pinned_commit":UPSTREAM_COMMIT,
      "original_git_blob_sha1":ORIGINAL_GITBLOB_SHA1,
      "original_raw_file_sha256":{"raw":h_raw,"subset":h_sub,"metadata":h_md},
      "n_original_raw_isotopic_specimens":184,
      "n_original_author_analysis_subset":155,
      "n_original_excluded_from_author_subset":29,
      "n_original_biological_and_water_metadata":380,
      "n_sites_with_original_isotopes":5,
      "n_original_exposure_labeled_sites":2,
      "n_original_unexposed_labeled_sites":3,
      "site_month_water_TP_unique_values_NOT_five_independent_timepoints":True,
      "GRVBR_exposed_site_has_no_original_TP_chemistry_value":True,
      "only_one_exposed_site_with_original_water_TP":True,
      "n_independent_randomized_WWTP_treatments":0,
      "original_2022_sampling_vs_2023_lab_analysis_dates_differ":True,
      "source_qa_by_site":per_site,
      "site_equal_exposed_vs_unexposed_DESCRIPTIVE_ONLY":comparisons,
      "original_aquatic_prey_and_riparian_spider_source_both_measured":True,
      "source_no_2022_measured_annual_wwtp_P_mass_joined":True,
      "source_does_NOT_identify_local_animal_diet_fraction_or_nutrient_flux_from_simple_isotope_differences":True,
      "source_dataset_reanalysis_not_a_new_independent_replicate_of_published_study":True,
      "source_no_freshwater_macro_wwtp_to_riparian_predator_H6_contrast":True}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"BowRiver2022_source_pinned_isotope_pathway_SITE_EQUAL_QA.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    rows=[]
    for x in per_site:
        for tax,values in x["original_d13C_d15N_by_taxon"].items():
            rows.append({
             "source_site":x["site"],"original_site_exposure":x["original_exposure_label"],
             "taxon":tax,"n_source_original_specimens":values["n"],
             "mean_d13C_permil":values["d13C_mean"],"mean_d15N_permil":values["d15N_mean"],
             "site_is_independent_unit_not_individual_specimen":True,
             "site_original_source_blob_SHA1":ORIGINAL_GITBLOB_SHA1["raw"]})
    with (args.out/"BowRiver2022_site_taxon_published_isotope_means_DESCRIPTIVE.csv").open(
         "w",encoding="utf-8",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print("GESUI_BOWRIVER_2022_ORIGINAL_WWTP_ISOTOPE_PATHWAY="+json.dumps({
      "status":report["status"],"n_original":184,"n_author_subset":155,
      "site_n":5,"exposed":2,"controls":3,
      "tetrag_d13C_exposed_minus_controls_permil":
        comparisons["tetragnathid_d13C_permil"]["exposed_minus_unexposed_site_equal_difference"],
      "tetrag_d13C_adjusted_relative_to_leaves_exposed_minus_controls_permil":
        comparisons["tetragnathid_d13C_minus_local_terrestrial_leaves_d13C_permil"]["exposed_minus_unexposed_site_equal_difference"],
      "tetrag_d15N_exposed_minus_controls_permil":
        comparisons["tetragnathid_d15N_permil"]["exposed_minus_unexposed_site_equal_difference"],
      "tetrag_d15N_adjusted_relative_to_leaves_exposed_minus_controls_permil":
        comparisons["tetragnathid_d15N_minus_local_terrestrial_leaves_d15N_permil"]["exposed_minus_unexposed_site_equal_difference"]
    }),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--raw",required=True,type=Path)
    p.add_argument("--subset",required=True,type=Path)
    p.add_argument("--metadata",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    run(p.parse_args())
