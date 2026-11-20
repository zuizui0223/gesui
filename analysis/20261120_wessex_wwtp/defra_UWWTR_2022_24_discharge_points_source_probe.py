#!/usr/bin/env python3
"""DEFRA UK UWWTR Article15 official ODS source-only download/schema gate.

Original water companies' Department for Environment, Food & Rural Affairs
CC-BY UWWTR Water Treatment in England datasets, publicly listed by
data.gov.uk. Official 2022/2024 Article15 ODS contains treatment plants,
discharge points and receiving waters in its dataset description.

CURRENT public source 2022 or 2024 MUST NOT be backfilled as a verified
2017–18 physical outlet without time-native plant/source identifiers and
independent stability checks. This source-only job checks URLs/downloads,
ODS ZIP member list and source identity, NOT permit joins/WWTP exposure,
nutrient concentrations, mass nor ecological outcomes.
"""
from __future__ import annotations
import argparse,hashlib,json,tempfile,urllib.error,urllib.request,zipfile
from pathlib import Path

DATASETS={
    "DEFRA_UWWTR_ARTICLE15_2022":{
        "url":"https://s3.eu-west-1.amazonaws.com/data.defra.gov.uk/Urban%2BWaste%2BWater%2BTreatment%2BDirective/Article15_13Dec2022.ods",
        "version_year":2022},
    "DEFRA_UWWTR_ARTICLE15_2024":{
        "url":"https://s3.eu-west-1.amazonaws.com/data.defra.gov.uk/Urban%2BWaste%2BWater%2BTreatment%2BDirective/UWWTR_Art15_24thOct2024.ods",
        "version_year":2024},
}
OFFICIAL_PAGE="https://www.data.gov.uk/dataset/d7e2c57b-110a-462b-97a0-9833e7d26cc2/wastewater-treatment-in-england"
LIMIT=50*1024*1024
HEADERS={
   "User-Agent":"gesui-original-DEFRA-UWWTR-Article15-outlet-source-gate/20261120",
   "Accept":"application/vnd.oasis.opendocument.spreadsheet,*/*",
}
def source(d,info,out):
    obj={"official_source_url":info["url"],"source_reference_year":info["version_year"],
         "status":"SOURCE_STOP","HEAD_status":None,"GET_status":None,
         "error":None}
    req=urllib.request.Request(info["url"],headers=HEADERS,method="HEAD")
    try:
        with urllib.request.urlopen(req,timeout=35) as r:
            obj["HEAD_status"]=r.status
            obj["HEAD_content_length"]=r.headers.get("Content-Length")
            obj["HEAD_content_type"]=r.headers.get("Content-Type")
    except Exception as e:
        obj["HEAD_error"]=type(e).__name__+": "+str(e)[:250]
    try:
        req=urllib.request.Request(info["url"],headers=HEADERS,method="GET")
        h=hashlib.sha256();n=0
        with urllib.request.urlopen(req,timeout=75) as r,(out/"source.ods").open("wb") as f:
            obj["GET_status"]=r.status
            if r.status!=200:raise RuntimeError("Original DEFRA source HTTP!=200")
            while True:
                chunk=r.read(65536)
                if not chunk:break
                n+=len(chunk)
                if n>LIMIT:raise RuntimeError("Original ODS larger than 50MB source-only cap")
                h.update(chunk);f.write(chunk)
        obj["source_bytes"]=n
        obj["source_SHA256"]=h.hexdigest()
        if not zipfile.is_zipfile(out/"source.ods"):
            raise RuntimeError("Official dataset URL did not return a ZIP/ODS package")
        with zipfile.ZipFile(out/"source.ods") as z:
            names=z.namelist()
            obj["zip_member_count"]=len(names)
            obj["original_member_names"]=sorted(names)
            obj["original_content_xml_uncompressed_bytes"]=(
                z.getinfo("content.xml").file_size if "content.xml" in names else None)
            typ=z.read("mimetype")[:150].decode("ascii","replace") if "mimetype" in names else None
            obj["original_mimetype"]=typ
            if "application/vnd.oasis.opendocument.spreadsheet"!=typ:
                raise RuntimeError("ODS official data mime declaration incorrect")
        obj["status"]="OFFICIAL_DEFRA_CC_BY_ODS_SOURCE_PACKAGE_OBTAINED_NOT_JOINED"
    except urllib.error.HTTPError as e:
        obj.update({"GET_status":e.code,"error":"HTTPError "+str(e)[:250]})
    except Exception as e:
        obj["error"]=type(e).__name__+": "+str(e)[:450]
    return obj

def main(a):
    result={
      "schema":"GESUI_DEFRA_ENGLISH_2022_2024_UWWTR_ART15_OUTLET_SOURCE_ONLY_V1",
      "official_govuk_source_metadata":OFFICIAL_PAGE,
      "official_dataset_licence":"Creative Commons Attribution",
      "original_publication_2022_2024_not_2017_18":True,
      "does_NOT_establish_2017_18_STW_outfall_location":True,
      "not_fitted_to_BMWP_or_EA_TP_measurements":True,
      "not_source_loads_NP_kgyr":True,
      "source_results":{}
    }
    with tempfile.TemporaryDirectory(prefix="gesui_DEFRA_ODS_source_") as folder:
        p=Path(folder)
        for name,entry in DATASETS.items():
            result["source_results"][name]=source(name,entry,p)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True,allow_nan=False)+"\n")
    print("GESUI_DEFRA_UWWTR_DISCHARGE_ODS_SOURCE_ONLY="+json.dumps({
        "status_by_year":{name:{"status":v["status"],"HEAD_status":v.get("HEAD_status"),
                           "GET_status":v.get("GET_status"),
                           "GET_error":v.get("error"),
                           "size":v.get("source_bytes"),
                           "SHA256":v.get("source_SHA256"),
                           "members":v.get("zip_member_count"),
                           "content_xml_bytes":v.get("original_content_xml_uncompressed_bytes")}
                          for name,v in result["source_results"].items()},
        "no_backfill_of_2017_18_outfalls":True,
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
