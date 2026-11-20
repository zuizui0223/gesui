#!/usr/bin/env python3
"""Source-only authentic Borealis Grand River 8 municipal WWTPs, two genuine
published datasets (effluent/biofilm/benthic vs larvae/adults/spiders/isotope).

Before reading numeric microplastic outcome magnitudes:
 - Re-fetch exactly the original publisher datafile IDs + SHA256.
 - Validate true WWTP IDs, up/down, time/year, life stage, measured units,
   taxon support, duplicated isotope sample IDs.
 - Never treat particles themselves as animal replicates or two papers as
   independent treatment source experiments.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,math,time,urllib.error,urllib.request
from collections import Counter,defaultdict
from pathlib import Path

BASE="https://borealisdata.ca/api/access/datafile/"
SCHEMA="GESUI_20261117_BOREALIS_GRAND8_TWO_DOI_SOURCE_LOCKED_WWTP_SITE_YEAR_STAGE_GATE_V1"
FILES={
 "effluent":{"id":1030117,"sha":"f79dd728f93e8c69efe5390e8730ede88024f9ed3b95c7c22d3215386746ae0b",
     "title":"Effluent_Summary.tab","doi":"10.5683/SP3/TLDY6W"},
 "biofilm":{"id":1030368,"sha":"d8cb0b90321cec6053118290b7e63058be2e44efb8f7c631083a415ba7afdf42",
     "title":"Biofilm_Summary.tab","doi":"10.5683/SP3/TLDY6W"},
 "sediment":{"id":1030371,"sha":"070a87d5d48dca54e8ceaabd104515d8f4ea0d48dc1e61e27e89f5a938291590",
     "title":"Sediment_Summary.tab","doi":"10.5683/SP3/TLDY6W"},
 "aquatic_2022":{"id":1057471,"sha":"aa9a27fed3df116ceeb4bdca954bdf75e5e6594f615094de6af3e1fcf6b85d19",
     "title":"AquaticInsect_Summary.tab","doi":"10.5683/SP3/TLDY6W"},
 "cross_species":{"id":1213261,"sha":"1b69eb88fe38f601b11ff466d6a9b3177546a9b489a74c61eaf2a64d770dec66",
     "title":"22.23.MP.Invert.Summary.tab","doi":"10.5683/SP4/U0MN9S"},
 "isotopes":{"id":1219236,"sha":"2cc350defbec3210bccaba49bab80be4286b6a916626c315c2f3e6f40181645e",
     "title":"Stable_Isotope_data.tab","doi":"10.5683/SP4/U0MN9S"}
}
EXPECTED={"effluent":77,"biofilm":240,"sediment":240,"aquatic_2022":426,"cross_species":906,"isotopes":918}
FIELDS={"effluent":["Site","time","MPs.per.L"],"biofilm":["Site","Location","Time","MP_per_g"],
 "sediment":["Site","Location","Time","MP_per_kg"],
 "aquatic_2022":["Site","Location","Time","taxa","MP_per_mg"],
 "cross_species":["Site","Location","Time","Life_Stage","Taxa","MP_per_mg"],
 "isotopes":["Site","Location","Year","Timepoint","Taxa","Sample_ID","d13C","d15N"]}
MAX_BYTES=2_500_000
def fetch(fileid):
    req=urllib.request.Request(BASE+str(fileid),headers={"User-Agent":"gesui-original-Borealis-eight-WWTP-source-stage-audit/20261117",
         "Accept":"text/tab-separated-values,text/plain,*/*"})
    with urllib.request.urlopen(req,timeout=90) as res:
        raw=res.read(MAX_BYTES+1)
        if len(raw)>MAX_BYTES:raise RuntimeError("Borealis original source larger than 2.5MB")
        return raw
def sha(b):return hashlib.sha256(b).hexdigest()
def scan(name,meta,raw):
    if sha(raw)!=meta["sha"]:
        raise RuntimeError("Original published Borealis DataFile SHA differs "+name)
    rr=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter="\t"))
    if len(rr)!=EXPECTED[name]:
        raise RuntimeError("Original source record count changed "+name)
    header=list(rr[0])
    if any(k not in header for k in FIELDS[name]):
        raise RuntimeError("Critical original published Site/Year/stage/metric column absent "+name)
    # Original publication site IDs and upstream/downstream taxonomy are
    # inherently pre-numeric. Do not inspect numeric response values for
    # source selection; quality denominator is only site/stage/time support.
    keys=[k for k in ["Site","Location","Life_Stage","Family","Taxa","taxa","Year",
        "Time","time","Timepoint","Place","Rep","rep","Sample_ID"] if k in header]
    levels={k:Counter(str(r.get(k) or "").strip() for r in rr) for k in keys}
    cells=Counter()
    for r in rr:
        key=tuple(str(r.get(x) or "").strip() for x in keys if x not in ("Rep","rep","Sample_ID","Place"))
        cells[key]+=1
    metric=FIELDS[name][-1]
    nonempty=sum(str(r.get(metric) or "").strip() not in ("","NA","NaN","N/A") for r in rr)
    return {"published_DOI":meta["doi"],"published_source_file_id":meta["id"],
      "original_filename":meta["title"],"source_sha256":meta["sha"],
      "original_real_record_n":len(rr),"original_columns":header,
      "original_source_fields_examined_pre_numeric":keys,
      "original_sample_site_ids":sorted(x for x in levels["Site"] if x),
      "original_site_count":len([x for x in levels["Site"] if x]),
      "original_category_levels":{
        k:dict(sorted(v.items())) for k,v in levels.items()
         if k!="Sample_ID" and len(v)<=40},
      "source_stage_x_year_x_WWTP_site_x_updown_counts":[
        {"keys":dict(zip([x for x in keys if x not in ("Rep","rep","Sample_ID","Place")],tuple_key)),
         "original_sample_records":int(n)}
        for tuple_key,n in sorted(cells.items())],
      "source_numeric_measurement_column":metric,
      "source_nonblank_numeric_field_record_n":nonempty,
      "source_microplastic_magnitudes_NOT_READ_OR_FITTED":True,
      "original_specimen_or_particle_duplicates_not_assumed_independent":True,
      "original_duplicate_lab_sample_ids":len(rr)-len([x for x in levels.get("Sample_ID",{}) if x])
          if "Sample_ID" in levels else None}
def run(a):
    result={"schema":SCHEMA,"status":"PASS_ORIGINAL_8_WWTP_TWO_DOI_SOURCE_METADATA_WHEN_FULL",
      "original_2026_published_two_DOIs":["10.5683/SP3/TLDY6W","10.5683/SP4/U0MN9S"],
      "source_file_count_in_strict_original_comparison":len(FILES),"source_files":{}}
    for name,meta in FILES.items():
        b=fetch(meta["id"])
        out=scan(name,meta,b)
        result["source_files"][name]=out
        print("BOREALIS_GRAND8_ORIGINAL_SOURCE="+json.dumps({
           "endpoint":name,"n":out["original_real_record_n"],
           "sites":out["original_site_count"],
           "original_levels":{k:v for k,v in out["original_category_levels"].items()
                              if k in ("Location","Life_Stage","Year","Time","time","Timepoint","Taxa","taxa")}},
           ensure_ascii=False,sort_keys=True)[:6500],flush=True)
    sites=[set(x["original_sample_site_ids"]) for x in result["source_files"].values()]
    if len(set.union(*sites))!=8 or any(len(s)!=8 for s in sites):
        raise RuntimeError("Original two-DOI shared WWTP site ID support not complete across six true matrices")
    if len(set.intersection(*sites))!=8:
        raise RuntimeError("Original physical source site keys are not exactly shared across data types")
    result["source_exact_8WWTP_site_names_common_to_six_original_files"]=sorted(set.intersection(*sites))
    result["original_unique_source_plant_count"]=8
    result["cross_DOI_site_name_join_pass"]=True
    result["same_original_year_and_biological_stage_cells_NOT_YET_QUALIFIED"]=True
    result["source_effect_or_particle_flux_NOT_YET_ESTIMATED"]=True
    result["no_treatment_replication_across_different_global_river_basins"]=True
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"GRAND8_REAL_EFFLUENT_BENTHIC_EMERGENCE_SPIDER_SITE_YEAR_SOURCE_QA.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("BOREALIS_GRAND8_SOURCE_FULL="+json.dumps(
        {"sites":8,"file_sources":len(FILES),"status":result["status"],
         "shared_site_ids":result["source_exact_8WWTP_site_names_common_to_six_original_files"]},
        ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    run(p.parse_args())
