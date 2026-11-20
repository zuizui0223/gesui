#!/usr/bin/env python3
"""Original USGS Farmington public PFAS assay to stable-isotope specimen ID audit.

Tests whether the *same biological specimen* has a supported field year,
matrix and taxon in the separate 3192-assay PFAS and 76 stable isotope rows.
Five focal streams in ONE Farmington watershed, NOT municipal effluent sites.
No chemical concentrations or isotope values read; no biological effects.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,re,time,urllib.request,urllib.error
from pathlib import Path

ORIGINAL="https://www.sciencebase.gov/catalog/item/69b80121b66b010d85582e8b?format=json"
SHA_META="a1ffa182ccc51c5f0018cda1cca99bda1a2407346f3d03fd9b716ae20f832524"
FILES={
 "Farmington_PFAS_FWstudy.csv":"f2d6cc636782506ec645b75d06c839e0a74fe01ee2fd526146536425f4ecc0d3",
 "Farmington_PFAS_Stable_Isotopes.csv":"963111a4d4ae2e3cbb510765c0b5acd3c20aa66cde2d779210c73ff5ca893872"
}
CANONICAL=("Burr Pond Brook","Hop Brook","Pequabuck River","Ratlum Brook","Russell Brook")
SCHEMA="GESUI_20261118_USGS_FARMINGTON_5STREAM_ASSAY_ISOTOPE_TRUE_SPECIMEN_YEAR_JOIN_GATE_V1"
def sha(x):return hashlib.sha256(x).hexdigest()
def norm(x):return re.sub(r"[^a-z0-9]","",str(x or "").lower())
def fetch(url):
    if not url.startswith("https://www.sciencebase.gov/"):raise RuntimeError("not USGS original")
    for trial in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={
                "User-Agent":"gesui-farmington-original-specimen-source-id-join/2026-11-18"}),timeout=75) as r:
                b=r.read(1_500_001)
                if len(b)>1_500_000:raise RuntimeError("original source unexpectedly oversized")
                return b
        except (urllib.error.URLError,urllib.error.HTTPError,TimeoutError):
            if trial==2:raise
            time.sleep(2+trial)
def source():
    b=fetch(ORIGINAL)
    if sha(b)!=SHA_META:raise RuntimeError("metadata changed")
    meta=json.loads(b)
    allfiles={x["name"]:x for x in meta["files"]}
    rows={}
    for name,expected in FILES.items():
        raw=fetch(allfiles[name]["downloadUri"])
        if sha(raw)!=expected:raise RuntimeError("original source sha mismatched: "+name)
        rr=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
        rows[name]=rr
    if len(rows["Farmington_PFAS_FWstudy.csv"])!=3192 or len(rows["Farmington_PFAS_Stable_Isotopes.csv"])!=76:
        raise RuntimeError("source counts changed")
    return rows
def typchem(s):
    s=str(s)
    if s=="Tetragnathidae":return "Spider"
    if s.startswith("Larval_"):return "Larval"
    if s.startswith("Emergent_"):return "Adult"
    if s=="Biofilm":return "Biofilm"
    return "Other"
def typiso(r):
    s=r["sample_type"]
    if s=="Spider":return "Spider"
    if s.startswith("Larval_"):return "Larval"
    if s.startswith("Emergent_"):return "Adult"
    if s=="Biofilm":return "Biofilm"
    return "Other"

def run(args):
    data=source()
    chem=data["Farmington_PFAS_FWstudy.csv"]
    isot=data["Farmington_PFAS_Stable_Isotopes.csv"]
    # Collapse each 12-analyte block to one genuine specimen; never infer that
    # one analyte row is one independent insect/spider replicate.
    specimens={}
    for r in chem:
        if r["site"]=="Field blank":continue
        site=r["site"]
        idkey=(site,r["field_id"])
        if idkey not in specimens:specimens[idkey]=r
        elif any(specimens[idkey][q]!=r[q] for q in ("sample_type","site")):
            raise RuntimeError("Original chemistry sample ID conflicts across assay rows")
    if len(specimens)!=258 or set(k[0] for k in specimens)!=set(CANONICAL):
        raise RuntimeError("Original chemical true specimen n258 or site names changed")
    # Distinct year/sample type in isotope source.
    iso_by_field=collections.defaultdict(list)
    for row in isot:
        if row["site"] not in CANONICAL:raise RuntimeError("unexpected isotope original field site")
        for candidate in ("u_conn_sample_id","usgs_sample_id"):
            key=norm(row.get(candidate))
            if key:
                iso_by_field[(row["site"],key)].append((candidate,row))
    pairs=[]
    ambiguity=0
    for (site,field_id),r in specimens.items():
        idvariants=[("field_id",r.get("field_id")),("cese_id",r.get("cese_id"))]
        found={}
        for nm,val in idvariants:
            nv=norm(val)
            if not nv:continue
            for what,iso in iso_by_field.get((site,nv),[]):
                found[id(iso)]=(nm,what,iso)
        if len(found)>1:ambiguity+=1
        if len(found)==1:
            from_ch,from_iso,iso=list(found.values())[0]
            pairs.append({"site":site,"chemical_stage":typchem(r["sample_type"]),
              "isotope_stage":typiso(iso),"isotope_real_year":iso["year"],
              "matching_fields":[from_ch,from_iso]})
    coverage=collections.Counter((p["site"],p["chemical_stage"],p["isotope_real_year"]) for p in pairs)
    stage_rows=[]
    for site in CANONICAL:
        for year in ("2022","2023"):
            c={}
            for cat in ("Larval","Adult","Spider","Biofilm"):
                c[cat]=int(coverage[(site,cat,year)])
            stage_rows.append({"site":site,"year":year,
             "exact_individual_specimen_id_matches":c,
             "three_biotic_stages_at_least_one_exact_original_match":all(c[k]>0 for k in ("Larval","Adult","Spider")),
             "three_biotic_stages_at_least_two_exact_original_matches":all(c[k]>=2 for k in ("Larval","Adult","Spider"))})
    levels=collections.Counter((r["site"],r["year"],typiso(r)) for r in isot)
    isotope_presence=[]
    for site in CANONICAL:
        for year in ("2022","2023"):
            x={k:int(levels[site,year,k]) for k in ("Larval","Adult","Spider","Biofilm")}
            isotope_presence.append({"site":site,"year":year,"original_isotope_source_group_counts":x,
              "original_sameyear_spider_and_larvae_and_adult_at_least_one":
                all(x[k]>0 for k in ("Larval","Adult","Spider")),
              "original_sameyear_spider_and_larvae_and_adult_at_least_two":
                all(x[k]>=2 for k in ("Larval","Adult","Spider"))})
    n1=sum(x["original_sameyear_spider_and_larvae_and_adult_at_least_one"] for x in isotope_presence)
    n2=sum(x["original_sameyear_spider_and_larvae_and_adult_at_least_two"] for x in isotope_presence)
    q1=sum(x["three_biotic_stages_at_least_one_exact_original_match"] for x in stage_rows)
    q2=sum(x["three_biotic_stages_at_least_two_exact_original_matches"] for x in stage_rows)
    report={
      "schema":SCHEMA,"status":"ORIGINAL_SOURCE_SPECIMEN_ID_JOIN_AUDIT_COMPLETED_NO_ECOLOGY",
      "n_published_original_PFAS_analyte_assay_rows_including_blanks":3192,
      "n_original_chemical_true_biotic_and_abiotic_field_specimens_excluding_blanks":258,
      "n_original_stable_isotope_data_rows":76,
      "n_exact_original_isotope_chemical_individual_ID_matches":len(pairs),
      "n_ambiguous_many_to_many_original_specimen_matches":ambiguity,
      "n_siteyear_isotope_any_3_biotic_stage_sameyear":n1,
      "n_siteyear_isotope_at_least2_all_3_stages_sameyear":n2,
      "n_siteyear_exact_specimen_linked_any_3_biotic_stage_sameyear":q1,
      "n_siteyear_exact_specimen_linked_at_least2_all_3_stages_sameyear":q2,
      "original_actual_isotope_siteyear_stage_source_counts":isotope_presence,
      "original_exact_specimenID_chemical_isotope_join_counts":stage_rows,
      "exact_join_stages_and_source_fields":dict(collections.Counter(
         ("/".join(p["matching_fields"]) for p in pairs)).most_common()),
      "no_real_annual_PFAS_response_value_or_isotope_magnitude_opened":True,
      "source_year_of_PFAS_insects_not_imputed_from_isotope_without_ID_match":True,
      "actual_source_independent_water_catchments":1,
      "source_USGS_not_outfall_exposed_unexposed_case_control":True,
      "original_input_sha256":FILES}
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"FARMINGTON_USGS_EXACT_BIOTIC_PFAS_ISOTOPE_SAMPLE_ID_YEAR_SOURCE_QA.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_USGS_FARMINGTON_ORIGINAL_EXACT_SPECIMEN_YEAR="+json.dumps({
        "status":report["status"],"exact_matched_specimens":len(pairs),
        "sameyear_all_3_biota_isotope_min1":n1,
        "exact_samples_sameyear_min1":q1,"exact_samples_sameyear_min2":q2},ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    run(p.parse_args())
