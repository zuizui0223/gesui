#!/usr/bin/env python3
"""Inspect OFFICIAL 2026 EEA WISE6 36GB individual-chemistry CSV header by
2MiB HTTP-range prefix from its 2GB published ZIP, without downloading it
or inspecting any numeric measurement / station sample records.
"""
from __future__ import annotations
import argparse,hashlib,json,struct,urllib.request,zlib
from datetime import datetime,timezone
from pathlib import Path
URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
EXPECTED_NAME="waterbase_t_wise6_disaggregateddata.csv"
MAX_PREFIX=2097152

def main():
    a=argparse.ArgumentParser();a.add_argument("--out",required=True,type=Path);p=a.parse_args()
    req=urllib.request.Request(URL,headers={
       "User-Agent":"gesui-2026-original-WISE6-2GB-ZIP-HEADER-only-range/1.0",
       "Range":f"bytes=0-{MAX_PREFIX-1}"})
    with urllib.request.urlopen(req,timeout=95) as r:
        if r.status!=206:raise ValueError("Original source did not honor bounded exact byte range")
        content_range=r.headers.get("Content-Range")
        data=r.read(MAX_PREFIX+1)
    if len(data)!=MAX_PREFIX or data[:4]!=b"PK\x03\x04":
        raise ValueError("Original public ZIP local-file header or source byte-range changed")
    compression=struct.unpack_from("<H",data,8)[0]
    n=struct.unpack_from("<H",data,26)[0]
    extra_n=struct.unpack_from("<H",data,28)[0]
    name=data[30:30+n].decode("utf8")
    if name!=EXPECTED_NAME or compression not in (8,):
        raise ValueError("Original ZIP file or compression not expected: "+name)
    start=30+n+extra_n
    h=zlib.decompressobj(-15)
    # Bound uncompressed inspection to first 64 KB of original CSV, not data.
    preview=h.decompress(data[start:],max_length=65536)
    header=preview.splitlines()[0].decode("utf-8-sig",errors="replace")
    sep=";" if header.count(";")>header.count(",") else ","
    names=header.split(sep)
    wanted=["countryCode","monitoringSiteIdentifier","monitoringSiteIdentifierScheme",
      "parameterWaterBodyCategory","observedPropertyDeterminandLabel","resultUom",
      "phenomenonTimeSamplingDate","resultObservedValue",
      "resultQualityObservedValueBelowLOQ","procedureLOQValue"]
    if len(names)<15:raise ValueError("Original Waterbase individual WISE6 source header truncated")
    output={
      "schema":"gesui_EEA2026_original_WISE6_36GB_disaggregated_sampling_2022TP_source_CSV_HEADER_ONLY_v01",
      "retrieved_UTC":datetime.now(timezone.utc).isoformat(),
      "official_published_WISE6_DisaggregatedData_URL":URL,
      "source_HTTP_206_content_range":content_range,
      "source_first_2MiB_compressed_range_SHA256":hashlib.sha256(data).hexdigest(),
      "original_ZIP_compression_method":compression,
      "original_ZIP_member_name":name,
      "original_WISE6_DisaggregatedData_CSV_field_delimiter":sep,
      "original_WISE6_DisaggregatedData_CSV_field_names":names,
      "required_2022_TP_filtered_retrieval_fields_present":{x:x in names for x in wanted},
      "original_numeric_concentrations_observation_rows_NOT_PARSED":True,
      "source_is_EEA2026_public_bulk_download_not_an_approximate_TP_dataset":True,
      "original_whole_Waterbase_DisaggregatedData_ZIP_NOT_DOWNLOADED":True,
      "original_sample_station_COUNTS_not_YET_available":True,
      "no_butterfly_observations_or_personal_data_read":True,
    }
    p.out.parent.mkdir(parents=True,exist_ok=True)
    p.out.write_text(json.dumps(output,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA_WISE6_ORIGINAL_36GB_INDIVIDUAL_TP_SOURCE_CSV_HEADER_ONLY",
      json.dumps({"column_count":len(names),"all_core_fields_present":all(output["required_2022_TP_filtered_retrieval_fields_present"].values()),
        "names":names,"full_2GB_original_data_not_downloaded":True},sort_keys=True),flush=True)
if __name__=="__main__":main()
