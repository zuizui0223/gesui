#!/usr/bin/env python3
"""Source-only public Bow River 2022 municipal wastewater aquatic-riparian foodweb
repository audit: stable isotopes/insects/spiders real group counts and columns.

No claims of independent WWTP replicates without original spatial treatment
metadata. Git clone is done separately; do not redistribute third-party raw.
"""
import argparse,csv,hashlib,json,os,re,subprocess
from pathlib import Path
from collections import Counter

OWNER="Karen-Kidd-Lab-Ecotoxicology"
PROJECT="Bow-River-Microbiome-2022"
SOURCE="https://github.com/"+OWNER+"/"+PROJECT
SCHEMA="GESUI_20261117_ORIGINAL_BOW_RIVER_WWTP_AQUATIC_RIPARIAN_DATA_SOURCE_FIRST_V1"
SUFFIXES={".csv",".tsv",".txt",".r",".R",".md",".xlsx",".rds",".RDS",".json"}
MAX_NFILE=250
MAX_TEXT_BYTES=6_000_000

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        while b:=f.read(1<<20):h.update(b)
    return h.hexdigest()
def rel(root,path):return str(path.relative_to(root))
def audit_csv(root,path):
    raw=path.read_bytes()
    text=raw.decode("utf-8-sig","replace")
    delim="\t" if path.suffix.lower()==".tsv" else ","
    recs=csv.DictReader(text.splitlines(),delimiter=delim)
    cols=list(recs.fieldnames or [])
    if not cols:return {"rows":0,"column_names":[],"error":"no_header"}
    cats={k:Counter() for k in cols if re.search(r"(site|location|river|reach|wwtp|outfall|down|upstream|treat|exposure|stream.*type|taxon|taxa|family|order|stage|aqu.*terr|sample.*type|matrix|tissue|year|month|date|city|group|habitat|consumer)",k,re.I)}
    n=0
    cross_site=Counter()
    cross_station_group=Counter()
    nonempty=Counter()
    for row in recs:
        if not any(v for v in row.values()):continue
        n+=1
        for k in cols:
            v=row.get(k,"")
            if v and str(v).strip():nonempty[k]+=1
        site=str(row.get("Site") or row.get("Site_code") or "").strip()
        exposure=str(row.get("Exposure") or "").strip()
        stream=str(row.get("Stream_Type") or "").strip()
        if site:cross_site[(site,exposure,stream)]+=1
        group=str(row.get("sample_group") or row.get("Stage") or "").strip()
        if site and group:cross_station_group[(site,group)]+=1
        for k in cats:
            value=str(row.get(k) or "").strip()
            if value and len(cats[k])<150:cats[k][value]+=1
    return {
      "original_real_data_records":n,
      "column_names":cols[:120],
      "column_count":len(cols),
      "candidate_true_biological_site_and_treatment_columns":list(cats),
      "category_distinct_count":{k:len(v) for k,v in cats.items()},
      "category_observed_n":{k:dict(v.most_common(30)) for k,v in cats.items() if len(v)<=35},
      "original_site_by_published_exposure_and_stream_type":[{"site":k[0],"source_exposure":k[1],
        "source_stream_type":k[2],"n_samples":v}
        for k,v in sorted(cross_site.items())],
      "original_site_by_biological_sample_group":[{"site":k[0],"group":k[1],"n_samples":v}
        for k,v in sorted(cross_station_group.items())],
      "nonempty_rows_per_candidate_column":{k:nonempty[k] for k in cats},
      "numeric_fields_not_opened_or_tested":True}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=a.original.resolve()
    if not (root/".git").exists():raise RuntimeError("Original public Git repository clone missing")
    actual=subprocess.check_output(["git","-C",str(root),"rev-parse","HEAD"],text=True).strip()
    remote=subprocess.check_output(["git","-C",str(root),"remote","get-url","origin"],text=True).strip()
    if SOURCE.lower() not in remote.lower():
        raise RuntimeError("Original repository remote not published author source")
    files=sorted(p for p in root.rglob("*") if p.is_file() and ".git" not in p.parts)
    if len(files)>MAX_NFILE:raise RuntimeError("Original source repository file count out of bound")
    out=[]
    for f in files:
        size=f.stat().st_size
        if f.suffix not in SUFFIXES and f.suffix.lower() not in {x.lower() for x in SUFFIXES}:continue
        r={"path":rel(root,f),"bytes":size,
           "sha256":sha(f),"kind":f.suffix.lower()}
        if f.suffix.lower() in (".csv",".tsv") and size<=MAX_TEXT_BYTES:
            r["source_data_shape"]=audit_csv(root,f)
        if f.suffix.lower()==".r" and size<MAX_TEXT_BYTES:
            t=f.read_text(encoding="utf-8",errors="replace")
            r["R_script_references_to_original_source_columns"]={
              "mentions_isotope":bool(re.search(r"(isotope|delta|d15N|d13C|15N|13C)",t,re.I)),
              "mentions_spiders":bool(re.search(r"(spider|Tetragnath)",t,re.I)),
              "mentions_WWTP_or_effluent":bool(re.search(r"(effluent|WWTP|wastewater)",t,re.I))}
        out.append(r)
    result={"schema":SCHEMA,
      "source_public_repo":SOURCE,
      "actual_source_commit_SHA":actual,
      "source_files_n":len(files),
      "source_file_schemas":out,
      "source_tables_with_original_real_data":sum("source_data_shape" in x for x in out),
      "source_tables_with_isotope_fields":[x["path"] for x in out if "source_data_shape" in x and
        any(re.search(r"(delta|d13c|d15n|isotope|13C|15N)",k,re.I)
            for k in x["source_data_shape"]["column_names"])],
      "source_updown_independent_WWTP_site_replicate_count_NOT_YET_ESTABLISHED":True,
      "no_aquatic_riparian_effect_or_treatment_model_fitted":True}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"BOW_RIVER_2022_ORIGINAL_PUBLIC_SOURCE_FIELD_CENSUS.json").write_text(
      json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_ORIGINAL_BOW_RIVER_SOURCE="+json.dumps({
      "sha":actual,"source_files":len(files),"n_real_tables":result["source_tables_with_original_real_data"],
      "isotope_table_candidates":result["source_tables_with_isotope_fields"]},ensure_ascii=False),flush=True)

if __name__=="__main__":main()
