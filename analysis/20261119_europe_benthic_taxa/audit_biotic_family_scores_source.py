#!/usr/bin/env python3
"""Source-only auditing the published UK family BMWP/WHPT score table.

Original EA Biosys Taxon Info has nomenclature but NO per-family BMWP scores.
This audits a public CRAN biotic package R source dictionary by immutable
Git commit, not a site-by-year biological fit. It does not assume a sample's
BMWP_ASPT is a valid per-taxon trait, nor silently assign a genus a score.
"""
from __future__ import annotations
import argparse, csv, hashlib, io, json, tempfile, urllib.request, zipfile
from pathlib import Path

import pyreadr

REF="08cea3f06eec93703bb10b53d1015d8fe8c626c4"
SOURCE="https://raw.githubusercontent.com/robbriers/biotic/"+REF+"/R/sysdata.rda"
EA_TAXA="https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip"

def fetch(url,maxsize):
    req=urllib.request.Request(url,headers={"User-Agent":"gesui-20261119-family-tolerance-source-only"})
    with urllib.request.urlopen(req,timeout=40) as r:
        raw=r.read(maxsize+1)
        if len(raw)>maxsize:raise RuntimeError("Source exceeds bounded exact family dictionary archive gate")
        return raw

def make(args):
    result={"schema":"GESUI_EA_UK_BENTHIC_FAMILY_TRAIT_REFERENCE_SOURCE_ONLY_V1",
            "status":"SOURCE_STOP",
            "official_EA_TaxonInfo_does_not_have_BMWP_per_family":True,
            "biotic_package_reference":"CRAN biotic 0.1.x 2016",
            "biotic_GPL3_source_commit":REF,
            "no_EA_biotic_observation_or_environment_effect_fitted":True,
            "no_automatic_species_to_family_inference":True}
    try:
        data=fetch(SOURCE,300000)
        familyzip=fetch(EA_TAXA,6000000)
        result["original_dictionary_rda_sha256"]=hashlib.sha256(data).hexdigest()
        result["original_EA_taxon_zip_sha256"]=hashlib.sha256(familyzip).hexdigest()
        with tempfile.TemporaryDirectory(prefix="gesui_trait_") as tmp:
            fp=Path(tmp)/"sysdata.rda";fp.write_bytes(data)
            x=pyreadr.read_r(str(fp))
        result["rdata_object_keys"]=list(x.keys())
        if "indextable" not in x:
            raise RuntimeError("Original biotic package lacked indextable")
        table=x["indextable"]
        if not {"Taxon","BMWP","WHPT"}.issubset(table.columns):
            raise RuntimeError("Original score table lacks named columns")
        if table["Taxon"].duplicated().any():
            raise RuntimeError("Scoring family dictionary has duplicate Taxon labels")
        with zipfile.ZipFile(io.BytesIO(familyzip)) as z:
            if sorted(z.namelist())!=["OPEN_DATA_TAXON_INFO.csv"]:
                raise RuntimeError("Official EA taxon ZIP file member changed")
            with z.open("OPEN_DATA_TAXON_INFO.csv") as f:
                txt=io.TextIOWrapper(f,encoding="utf-8-sig",errors="replace",newline="")
                rows=list(csv.DictReader(txt))
        if not all("TAXON_RANK" in r and "PREFERRED_TAXON_NAME" in r for r in rows):
            raise RuntimeError("EA taxonomy source lacks rank/name columns")
        family={r["PREFERRED_TAXON_NAME"].strip() for r in rows
                if r["TAXON_RANK"].strip().upper()=="FAMILY" and r["PREFERRED_TAXON_NAME"].strip()}
        score=table.loc[table["BMWP"].notna()]
        scored=set(str(k).strip() for k in score.Taxon)
        exact=sorted(scored & family)
        result.update({
            "status":"EXTERNAL_PUBLISHED_TRAIT_REFERENCE_AND_EA_NOMENCLATURE_VERIFIED",
            "n_EA_taxon_rows":len(rows),
            "n_EA_family_taxa_distinct":len(family),
            "n_published_scored_BMWP_taxa":len(scored),
            "n_EA_exact_family_scored":len(exact),
            "n_EA_family_scored_fraction":round(len(exact)/len(family),5) if family else None,
            "published_scoring_columns":list(table.columns),
            "n_published_scoring_table_rows":len(table),
            "not_original_EA_2026_official_version_of_scoring":True,
            "not_a_true_scored_sample_coverage":True,
            "taxon_resolution_mapping_and_indices_need_independent_validation":True,
            "no_unscored_trait_imputation":True
        })
    except Exception as e:
        result.update({"status":"SOURCE_STOP","error":type(e).__name__+": "+str(e)[:400]})
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("GESUI_EA_TRAIT_REFERENCE_SOURCE_ONLY="+json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    make(p.parse_args())
