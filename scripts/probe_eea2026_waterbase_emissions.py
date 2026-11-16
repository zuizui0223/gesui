#!/usr/bin/env python3
"""Official 2026 EEA Waterbase Emissions metadata: agricultural/urban N and P loads.

No biological EQR scores or downstream effect outcomes used. Source ZIP is
only inventoried and scientific table headers are fingerprinted; full
original EEA data are not redistributed.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,zipfile
from datetime import datetime,timezone
from pathlib import Path
BASE="https://www.eea.europa.eu/en/datahub/datahubitem-view/5bd950b1-52b9-42f3-b0d7-bf51ff3056b7"
DOWNLOAD="https://sdi.eea.europa.eu/datashare/s/KkqFMjiif9eS7bn/download"
MAX_ENTRIES=250
MAX_PREVIEW=32768
def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1048576),b""):
            h.update(buf)
    return h.hexdigest()
def read_schema(z,info):
    name=info.filename
    output={"file":name,"size_uncompressed":info.file_size,
            "size_compressed":info.compress_size,"CRC32":format(info.CRC,"08x")}
    if info.filename.lower().endswith((".csv",".txt",".sql")):
        with z.open(info) as f:buf=f.read(MAX_PREVIEW)
        output["sample_32kb_SHA256"]=hashlib.sha256(buf).hexdigest()
        lines=buf.decode("utf-8-sig","replace").splitlines()
        if lines:
            first=lines[0]
            sep=max(("\t",";",","),key=lambda d:first.count(d))
            try:output["first_record_as_column_names"]=next(csv.reader([first],delimiter=sep))[:150]
            except Exception as e:output["header_error"]=str(e)[:90]
        output["source_delimiter_inferred"]=sep if lines else None
    return output
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    if not zipfile.is_zipfile(a.archive):raise ValueError("Official EEA Emissions download not a valid ZIP")
    zsha=digest(a.archive)
    with zipfile.ZipFile(a.archive) as z:
        entries=[x for x in z.infolist() if not x.is_dir()]
        if len(entries)>MAX_ENTRIES:raise ValueError("Source archive contains unexpectedly many files")
        records=[read_schema(z,x) for x in entries]
    likely=[x for x in records if any(k in x["file"].casefold()
       for k in ("emission","riverineinput","definition","spatial","README".casefold()))]
    obj={
      "schema":"gesui_eea_Waterbase_Emissions_2026_source_inventory_v01",
      "EEA_public_dataset":BASE,"EEA_exact_download":DOWNLOAD,
      "EEA_version":"eea_t_waterbase-emissions_p_1977-2025_v01_r00",
      "retrieved_utc":datetime.now(timezone.utc).isoformat(),
      "archive_container_sha256":zsha,
      "original_archive_bytes":a.archive.stat().st_size,
      "members":records,
      "candidate_emissions_and_spatial_members":likely,
      "scientific_N_or_P_year2022_values_read":False,
      "actual_2022_Emissions_agriculture_load_available_not_yet_verified":True,
      "original_EEA_full_archive_redistributed":False,
      "biology_EQR_outcomes_or_model_residuals_accessed":False,
      "limitations":[
       "Waterbase Emissions source data can be annual country/catchment aggregates instead of individual receiving reaches",
       "The source may have a subset of countries or no 2022 observations for some diffuse agriculture emissions",
       "Diffuse source codes and measurement estimation techniques require original codebook verification",
       "If agricultural emissions are assigned via ecological condition, conditioning on them can cause circularity",
       "Do not label emissions as river-sampled measured TN/TP concentrations or fixed WWTP effluent contributions"
      ]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2026_EEA_EMISSIONS_SOURCE_INVENTORY",json.dumps({
       "original_zip_sha256":zsha,"archive_bytes":obj["original_archive_bytes"],
       "total_members":len(records),"candidate_files":[x["file"] for x in likely][:35],
       "no_ecological_model_run":True},sort_keys=True),flush=True)
    for info in likely[:24]:
        print("GESUI_2026_EEA_EMISSIONS_SOURCE_SCHEMA",json.dumps(
          {"name":info["file"],"uncompressed_bytes":info["size_uncompressed"],
           "columns":info.get("first_record_as_column_names",[])},sort_keys=True),flush=True)
if __name__=="__main__":
    main()
