#!/usr/bin/env python3
"""Public USGS Farmington 5-stream PFAS source opportunity with exact originals.

True source DOI 10.5066/P139C2SN, ScienceBase 69b80121b66b010d85582e8b.
This is one river watershed with five sampled tributary streams, NOT five
independent wastewater treatment catchments or an effluent-up/down treatment.
Inspect sources, site/year/stage/assay grain and explicit one-site alias ONLY;
do not load concentrations, assign site exposure, or claim transfer quantity.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,collections,re,time,urllib.parse,urllib.request,urllib.error
from pathlib import Path

SCIENCEBASE="https://www.sciencebase.gov/catalog/item/69b80121b66b010d85582e8b?format=json"
SHA_META="a1ffa182ccc51c5f0018cda1cca99bda1a2407346f3d03fd9b716ae20f832524"
FILES={
 "Farmington_InsectMasses.csv":"76e7e0ebd4e2e77be4b7bc8bbcb93de5460714fe140a663682495db709b08283",
 "Farmington_PFAS_FWstudy.csv":"f2d6cc636782506ec645b75d06c839e0a74fe01ee2fd526146536425f4ecc0d3",
 "Farmington_PFAS_Stable_Isotopes.csv":"963111a4d4ae2e3cbb510765c0b5acd3c20aa66cde2d779210c73ff5ca893872",
 "SiteLocations.csv":"8ab6d062d199d3bbb1aa45f6b1e4920b81469fdb0d187298131a56c2cd07dc1e"
}
NAME_ALIAS={"Pequabuck Brook":"Pequabuck River"}
CANONICAL=("Burr Pond Brook","Hop Brook","Pequabuck River","Ratlum Brook","Russell Brook")
SCHEMA="GESUI_20261118_USGS_FARMINGTON_5STREAM_REAL_CHEM_INSECT_MASS_ISOTOPE_SITE_YEAR_STAGE_SOURCE_GATE_V1"
MAX_BYTES=1_500_000
def sha(x):return hashlib.sha256(x).hexdigest()
def get(url,maxbytes=MAX_BYTES):
    if not url.startswith("https://www.sciencebase.gov/"):raise RuntimeError("Untrusted non-USGS original data URL")
    for trial in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(url,headers={
                "User-Agent":"gesui-source-locked-usgs-Farmington-five-stream/2026-11-18",
                "Accept":"application/json,text/csv,*/*"}),timeout=70) as res:
                b=res.read(maxbytes+1)
                if len(b)>maxbytes:raise RuntimeError("USGS published file exceeded source safety cap")
                return b
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
            if trial==2 or isinstance(e,urllib.error.HTTPError) and e.code not in (429,500,502,503,504):
                raise
            time.sleep(2+trial*2)
    raise RuntimeError("Unreachable")

def csv_rows(data,file):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))
def standard(name):return NAME_ALIAS.get(name,name)
def nonblank(x):return x is not None and str(x).strip() not in ("","NA","N/A","NaN","nan")
def calendar(v):
    s=str(v or "").strip()
    y=re.search(r"(?:19|20)\d\d",s)
    return y.group(0) if y else "unresolved_date"

def aggregate(chem,isotope,masses,locations):
    if [len(chem),len(isotope),len(masses),len(locations)]!=[3192,76,128,22]:
        raise RuntimeError("Original four quantitative-source census sizes changed")
    avail={v["Site"] for v in locations}
    if not set(CANONICAL).issubset(avail):
        raise RuntimeError("Original five natural stream names not in independent official location file")
    if {v["site"] for v in masses}!={*CANONICAL,"Pequabuck Brook"}-{ "Pequabuck River"}:
        raise RuntimeError("The original single literature site-name alias changed")
    for file,rs in (("PFAS",chem),("ISOTOPE",isotope),("MASSES",masses)):
        for r in rs:
            t=standard(r.get("site",""))
            # Published instrument quality-control blanks are not stream
            # sites or extra river replicates. Preserve and count separately.
            if file=="PFAS" and t=="Field blank":continue
            if t not in CANONICAL:raise RuntimeError(f"Original {file} includes unexpected sampled stream {t}")
    for r in chem:
        if r.get("units") not in ("ng/g ww","ng/mL"):
            raise RuntimeError("Original chemical tissue versus water units changed")
    metric_fields=set(chem[0])
    sample_class_fields=["site","sample_type","collected","units","field_id","cese_id","pfas_type"]
    if not set(sample_class_fields).issubset(metric_fields):raise RuntimeError("Original FWstudy fields missing")
    actual_types=sorted({v["sample_type"] for v in chem})
    unique_sample_groups=collections.defaultdict(set)
    assay_rows=collections.Counter()
    assay_compounds=collections.defaultdict(set)
    blankrows=sum(str(r["site"]).strip()=="Field blank" for r in chem)
    for r in chem:
        if str(r["site"]).strip()=="Field blank":continue
        site=standard(r["site"])
        year=calendar(r["collected"])
        # Field-ID is specimen/sampling unit. Compounds analysed on the same
        # specimen are NOT counted as independent biological samples.
        key=(site,year,r["sample_type"],r["units"])
        specimen=str(r["field_id"]).strip()
        if nonblank(specimen):unique_sample_groups[key].add(specimen)
        assay_rows[key]+=1
        assay_compounds[key].add(r["pfas_type"])
    chem_source=[{"site":x[0],"year":x[1],"sample_type":x[2],"units":x[3],
          "source_independent_field_id_n":len(unique_sample_groups[x]),
          "original_compound_assay_rows":n,
          "source_distinct_PFAS_analyte_count":len(assay_compounds[x])}
         for x,n in sorted(assay_rows.items())]
    mass_groups=collections.Counter((standard(x["site"]),x["year"],x["life_stage"],x.get("order","")) for x in masses)
    isotope_groups=collections.Counter((standard(x["site"]),x["year"],x["sample_type"],x["lifestage"],x.get("order",""))
                                       for x in isotope)
    return {
       "site_alias_from_published_biotic_masses_to_all_other_USGS_original_tables":NAME_ALIAS,
       "canonical_site_count":len(CANONICAL),
       "independent_water_catchment_count":1,
       "source_chemical_3192_rows_true_field_ids_not_false_replicates":True,
       "original_PFAS_sample_type_levels":actual_types,
       "original_USGS_field_blank_assay_rows_excluded_from_five_biological_rivers":blankrows,
       "original_all_PFAS_rows_including_QA_field_blanks":len(chem),
       "original_nonblank_river_PFAS_assays":sum(assay_rows.values()),
       "original_PFAS_site_year_stage_field_sample_support":chem_source,
       "original_insect_mass_site_year_stage_order_count":[{"site":k[0],"year":k[1],"life_stage":k[2],"order":k[3],"n":v}
         for k,v in sorted(mass_groups.items())],
       "original_isotope_site_year_sample_type_stage_order_count":[
         {"site":k[0],"year":k[1],"sample_type":k[2],"life_stage":k[3],"order":k[4],"n":v}
         for k,v in sorted(isotope_groups.items())],
       "source_wwtp_exposure_upsream_downstream_NOT_present":True,
       "source_nonconcentration_numerical_effects_opened":False
    }

def run(args):
    args.out.mkdir(parents=True,exist_ok=True)
    report={"schema":SCHEMA,"source_DOI":"10.5066/P139C2SN",
      "source_2026_published_article":"10.1021/acs.est.5c16942",
      "expected_source_metadata_SHA256":SHA_META,
      "source_file_expected_SHA256":FILES,
      "public_data_CC0_USGS_release":True}
    try:
        md=get(SCIENCEBASE,maxbytes=7_000_000)
        if sha(md)!=SHA_META:raise RuntimeError("USGS original source metadata changed from first verified file collection")
        meta=json.loads(md)
        files={f["name"]:f for f in meta["files"]}
        if not all(p in files for p in FILES):raise RuntimeError("Original USGS four data filenames missing")
        payloads={}
        for name,expected in FILES.items():
            url=files[name]["downloadUri"]
            b=get(url)
            if sha(b)!=expected:raise RuntimeError("Original source SHA256 mutated for "+name)
            payloads[name]=csv_rows(b,name)
            print("FARMINGTON_TRUE_SOURCE_FILE_VERIFIED="+json.dumps(
              {"name":name,"bytes":len(b),"source_original_n_rows":len(payloads[name])},ensure_ascii=False),flush=True)
        report.update(aggregate(payloads["Farmington_PFAS_FWstudy.csv"],
          payloads["Farmington_PFAS_Stable_Isotopes.csv"],
          payloads["Farmington_InsectMasses.csv"],payloads["SiteLocations.csv"]))
        report["status"]="USGS_ORIGINAL_SOURCE_FOUR_TABLE_TRUE_SITE_YEAR_STAGE_GATE_AUDITED_NO_PFAS_EFFECT"
    except Exception as ex:
        report["status"]="STOP_PUBLISHED_SOURCE_ID_OR_SCHEMA_NOT_ECOLOGY_NULL"
        report["source_error_type"]=type(ex).__name__
        report["source_error_summary"]=str(ex)[:450]
    (args.out/"USGS_FARMINGTON_FIVE_STREAM_PFAS_ORIGINAL_CHEM_BIOTA_SITE_YEAR_STAGE_GATE.json").write_text(
       json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_FARMINGTON_FIVE_STREAM_STAGE="+json.dumps({
      "status":report["status"],"n_source_site":report.get("canonical_site_count"),
      "source_types":report.get("original_PFAS_sample_type_levels"),
      "error":report.get("source_error_summary")},ensure_ascii=False),flush=True)
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--out",required=True,type=Path)
    run(a.parse_args())
