#!/usr/bin/env python3
"""USGS Kotalik 2022 one industrial outfall: source-gated two PFAS chemicals
across effluent, water, caddisfly larvae, emergent adults and riparian spiders.

Uses original 2022 USGS Table1, one August 1–2 field event, not independent
outfalls/years or direct organism feeding. Chemicals PFOS and 6:2 FTS are
selected from the *previously published* 2025 paper's contrasting mechanisms,
not posthoc best-sign candidates in this source.
Protect against nondetect misclassification, mixed ng/L vs ng/g, and missing
individual specimen IDs. No mixing model, p values, or causal transfer rates.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,time,urllib.request,urllib.error
from pathlib import Path
URL="https://www.sciencebase.gov/catalog/item/68af0ddcd4be0253bb2c7c6b?format=json"
SOURCE_SHA="05f09778feb45fad755c64216ac37b95246a80dd3c15ebcbbbf31155b78f8970"
FILE="SoyoilFoodWeb_Table1_2022Results.txt"
SCHEMA="GESUI_20261118_INDUSTRIAL_2022_AUG_ONE_EVENT_PFOS_6_2FTS_SOURCE_GATED_STAGE_QA_V1"
CHEM=("PFOS","6:2 FTS")
GROUPS={
 "effluent":("Outfall","Effluent","NA"),
 "river_upstream":("Upstream","Surface Water","NA"),
 "river_downstream":("Downstream","Surface Water","NA"),
 "trichoptera_larvae_up":("Upstream","Trichoptera","Larvae"),
 "trichoptera_larvae_down":("Downstream","Trichoptera","Larvae"),
 "trichoptera_adult_up":("Upstream","Trichoptera","Adult"),
 "trichoptera_adult_down":("Downstream","Trichoptera","Adult"),
 "spider_up":("Upstream","Tetragnathidae","NA"),
 "spider_down":("Downstream","Tetragnathidae","NA")
}
BIO_PAIRS=(("trichoptera_larvae_up","trichoptera_larvae_down"),
           ("trichoptera_adult_up","trichoptera_adult_down"),
           ("spider_up","spider_down"))
DATE_SCOPE={"20220801","20220802"}
MIN_ANALYTES_PER_EACH_BIO_GROUP=2  # before reading values
def sha(b):return hashlib.sha256(b).hexdigest()
def fetch(url):
    if not url.startswith("https://www.sciencebase.gov/"):raise RuntimeError("Untrusted source")
    last=None
    for k in range(3):
        try:
            request=urllib.request.Request(url,headers={"User-Agent":"gesui-20261118-source-locked-industrial-PFAS-tissue-only/1","Accept":"text/plain,application/json,*/*"})
            with urllib.request.urlopen(request,timeout=75) as r:
                b=r.read(1_000_001)
                if len(b)>1_000_000:raise RuntimeError("Too many original bytes")
                return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            last=e
            if k==2:raise RuntimeError("USGS original access failure: "+str(last)[:250])
            time.sleep(2+2*k)
def read():
    item=json.loads(fetch(URL))
    if item.get("id")!="68af0ddcd4be0253bb2c7c6b":raise RuntimeError("Original USGS item wrong")
    meta={f["name"]:f for f in item["files"]}
    raw=fetch(meta[FILE]["downloadUri"])
    if sha(raw)!=SOURCE_SHA:raise RuntimeError("Original fixed full file SHA changed")
    rows=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter="\t"))
    if len(rows)!=2839:raise RuntimeError("Source full 2022 table expected 2839")
    return rows

def parse(x):
    try:
        v=float(x)
        return v if math.isfinite(v) and v>=0 else None
    except (ValueError,TypeError):return None

def run(args):
    rows=read()
    original={}
    support={}
    for key,(site,typ,stage) in GROUPS.items():
        original[key]={}
        for chem in CHEM:
            q=[r for r in rows if r["Site"].strip()==site and r["Sample_type"].strip()==typ
               and r["Lifestage"].strip()==stage and r["Parameter"].strip()==chem
               and str(r["Sample_datetime"])[:8] in DATE_SCOPE]
            units=sorted(set(r["Unit"].strip() for r in q))
            dates=sorted(set(str(r["Sample_datetime"])[:8] for r in q))
            remark_counts=dict(collections.Counter(r["Remark"].strip() for r in q))
            records={"source_original_assay_rows":len(q),
               "original_location":site,"original_type":typ,"original_lifestage":stage,
               "original_year":2022,"original_field_dates":dates,
               "original_unit_as_published":units,
               "original_Remark_censor_codes":remark_counts,
               "unique_organism_specimen_ids_absent_from_USGS_2022_header":True}
            original[key][chem]=(q,records)
    flags=[]
    for chem in CHEM:
        for group,record in original.items():
            size=len(record[chem][0])
            need=1 if group in ("effluent","river_upstream","river_downstream") else MIN_ANALYTES_PER_EACH_BIO_GROUP
            ok=size>=need
            flags.append({"chem":chem,"group":group,"n":size,"required":need,"supported":ok})
        for a,b in BIO_PAIRS:
            sa=original[a][chem][1]["original_unit_as_published"]
            sb=original[b][chem][1]["original_unit_as_published"]
            if len(sa)!=1 or sa!=sb:
                flags.append({"chem":chem,"group":a+"=>"+b,"supported":False,
                     "failure":"original_upstream_downstream_same_taxon_units_missing_or_different",
                     "unit_upstream":sa,"unit_downstream":sb})
    source_pass=all(x["supported"] for x in flags)
    result={"schema":SCHEMA,"status":"GO_POSTHOC_SINGLE_EVENT_STAGE_SOURCE_FOR_DESCRIPTIVE" if source_pass else
                "STOP_SOURCE_STAGE_OR_UNIT_INELIGIBLE_NO_NUMERIC_PFAS_PATH",
       "original_2022_full_source_SHA256":SOURCE_SHA,
       "original_2022_PFAS_table_2839_analyte_assay_rows_not_2839_independent_animals":True,
       "one_actual_industrial_outfall_and_one_reach":True,
       "physical_Bow_or_Grand_independent_WWTP_replication_NOT_present":True,
       "only_original_2022_aug_1_2_field_event":True,
       "original_chemical_candidates_from_published_contrasting_transfer_mechanisms":list(CHEM),
       "source_primary_taxon_lifestage_and_units_gate_before_magnitude":flags,
       "source_group_by_chemical_eligibility":{
         g:{chem:q[1] for chem,q in entries.items()} for g,entries in original.items()},
       "no_report_level_censored_observations_imputed_to_zero":True,
       "same_aquatic_larvae_adult_individual_paired_metamerphosis_NOT_TRACKED":True,
       "sample_organism_individual_ID_not_present":True,
       "no_claim_of_direct_consumption_or_trophic_transfer_factor":True}
    if source_pass:
        # descriptive chemistry only, all same-source physical matrices in
        # their own original tissue/unit scale.
        chemobj={}
        for group,chem_map in original.items():
            chemobj[group]={}
            for chem,(q,prior) in chem_map.items():
                detected=[v for r in q if str(r["Remark"]).strip()!="BDL" and
                         (v:=parse(r["Result_value"])) is not None]
                n_above=len(detected)
                chemobj[group][chem]={
                    "n_original_chemical_assay_rows":len(q),
                    "n_source_detected_non_censored_numeric":n_above,
                    "n_source_BDL":sum(r["Remark"].strip()=="BDL" for r in q),
                    "original_units":prior["original_unit_as_published"],
                    "observed_detected_only_arithmetic_mean":sum(detected)/len(detected) if detected else None,
                    "observed_detected_only_median":sorted(detected)[len(detected)//2] if detected else None,
                    "not_all_source_samples_detected":n_above!=len(q)}
        result["observed_only_not_complete_censored_stage_chemistry"]=chemobj
        contrasts={}
        for chem in CHEM:
          contrasts[chem]={}
          for up,down in (("river_upstream","river_downstream"),*BIO_PAIRS):
            x=chemobj[up][chem];y=chemobj[down][chem]
            z={"original_unit_same_matrix":x["original_units"],
               "n_detected_upstream":x["n_source_detected_non_censored_numeric"],
               "n_detected_downstream":y["n_source_detected_non_censored_numeric"],
               "detected_only_mean_difference":None,"censor_free_mean_ratio":None,
               "cannot_compute_true_unbiased_concentration_ratio_with_source_BDL":True}
            if (x["n_source_detected_non_censored_numeric"]>=2 and
                y["n_source_detected_non_censored_numeric"]>=2):
                u=x["observed_detected_only_arithmetic_mean"]
                d=y["observed_detected_only_arithmetic_mean"]
                z["detected_only_mean_difference"]=d-u
                if u and u>0:z["detected_only_mean_ratio_descriptive_selection_bias"]=d/u
            contrasts[chem][up+"_to_"+down]=z
        result["within_matrix_detected_only_descriptive_contrasts"]=contrasts
    else:
        result["observed_only_not_complete_censored_stage_chemistry"]=None
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"USGS_2022_INDUSTRIAL_AUG1_2_PFOS_6_2FTS_BIOLOGICAL_STAGE_SOURCE_QA.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("USGS_INDUSTRIAL_ORIGINAL_2022_TRACER="+json.dumps({
       "status":result["status"],
       "source_groups":{g:{c:q[1]["source_original_assay_rows"] for c,q in rs.items()}
          for g,rs in original.items()},
       "source_2022_cross_biota_numeric_provisional":result.get("within_matrix_detected_only_descriptive_contrasts")
    },ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
