#!/usr/bin/env python3
"""EEA 2026 WISE6 2022 individual water-quality SOURCE ZIP capability audit.

Original public WISE6 DisaggregatedData bulk ZIP is >2GB compressed. This
script makes HEAD and *bounded* HTTP byte-range probes and decodes ZIP central
directory size metadata without downloading all original observations.
No numeric water chemistry, station IDs, personal data or butterfly data used.
Do not treat access/size barrier as biological absence.
"""
from __future__ import annotations
import argparse,hashlib,json,re,struct,urllib.request,urllib.error
from datetime import datetime,timezone
from pathlib import Path
URL="https://discodata.eea.europa.eu/download/WISE_SOE/latest/Waterbase_T_WISE6_DisaggregatedData"
MAX_RESPONSE_READ=140000
TAIL_BYTES=131072

def ask(headers=None,method="GET"):
    hdr={"User-Agent":"gesui-EEA2026-original-WISE6-individual-TP-bulk-ZIP-safe-range/1.0",
         "Accept":"application/zip,*/*"}
    if headers:hdr.update(headers)
    resp={"request_range":hdr.get("Range"),"method":method,"original_url":URL}
    try:
        req=urllib.request.Request(URL,headers=hdr,method=method)
        with urllib.request.urlopen(req,timeout=50) as r:
            body=r.read(MAX_RESPONSE_READ) if method=="GET" else b""
            resp.update({"http":r.status,
             "content_length_header":r.headers.get("Content-Length"),
             "content_range_header":r.headers.get("Content-Range"),
             "accept_ranges_header":r.headers.get("Accept-Ranges"),
             "content_type":r.headers.get("Content-Type"),
             "body_bytes_read":len(body),
             "bounded_body_sha256":hashlib.sha256(body).hexdigest() if body else None})
        return resp,body
    except urllib.error.HTTPError as e:
        resp.update({"http":e.code,"error":str(e)});return resp,b""
    except Exception as e:
        resp["error"]=type(e).__name__+": "+str(e);return resp,b""

def zipcentral(blob):
    out={"range_response_contains_ZIP_EOCD":False,"zip_member_metadata":[]}
    i=blob.rfind(b"PK\x05\x06")
    if i<0 or i+22>len(blob):return out
    out["range_response_contains_ZIP_EOCD"]=True
    (_sig,disk,disk2,on_disk,entries,cdsize,cdoffset,comment)=struct.unpack_from("<IHHHHIIH",blob,i)
    out.update({"ZIP_center_total_original_entries":entries,
        "ZIP_CD_size_bytes":cdsize,"ZIP_CD_uncompressed_offset":cdoffset,
        "ZIP64_or_unusually_large_classic_offsets_possible":bool(cdsize==0xffffffff or cdoffset==0xffffffff)})
    for m in re.finditer(re.escape(b"PK\x01\x02"),blob):
        pos=m.start()
        if pos+46>len(blob):continue
        cmp=struct.unpack_from("<I",blob,pos+20)[0]
        unc=struct.unpack_from("<I",blob,pos+24)[0]
        n=struct.unpack_from("<H",blob,pos+28)[0]
        e=struct.unpack_from("<H",blob,pos+30)[0]
        com=struct.unpack_from("<H",blob,pos+32)[0]
        if pos+46+n+e+com>len(blob):continue
        fn=blob[pos+46:pos+46+n].decode("utf8",errors="replace")
        extra=blob[pos+46+n:pos+46+n+e]
        z64={}
        k=0
        while k+4<=len(extra):
            h,sz=struct.unpack_from("<HH",extra,k)
            val=extra[k+4:k+4+sz]
            if h==1:
                ptr=0
                if unc==0xffffffff and ptr+8<=len(val):
                    z64["uncompressed"]=struct.unpack_from("<Q",val,ptr)[0];ptr+=8
                if cmp==0xffffffff and ptr+8<=len(val):
                    z64["compressed"]=struct.unpack_from("<Q",val,ptr)[0];ptr+=8
            k+=4+sz
        out["zip_member_metadata"].append({"name":fn,
           "compressed_bytes":z64.get("compressed",cmp),
           "uncompressed_bytes":z64.get("uncompressed",unc),
           "from_original_range_response":True})
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    head,_=ask(method="HEAD")
    start,part0=ask({"Range":"bytes=0-1023"})
    tail,part_end=ask({"Range":f"bytes=-{TAIL_BYTES}"})
    midata=zipcentral(part_end) if tail.get("http")==206 else {
       "range_response_contains_ZIP_EOCD":False,
       "reason":"Server did not return 206 partial-content for tail probe"}
    obj={
        "schema":"gesui_official_EEA2026_WISE6_original_2022_individual_TP_source_zip_range_and_size_gate_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "EEA2026_original_published_full_individual_measurement_ZIP_URL":URL,
        "HEAD":head,"bounded_range_start":start,"bounded_range_tail":tail,
        "original_ZIP_directory_when_range_supported":midata,
        "no_full_2GB_2026_individual_chemistry_ZIP_downloaded":True,
        "no_2022_TP_concentration_or_site_data_rows_read":True,
        "no_2022_butterfly_photo_or_abundance_record_accessed":True,
        "source_only_failure_not_absence_of_DE_AT_2022_River_TP":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_EEA2026_ORIGINAL_INDIVIDUAL_WATERQUALITY_BOUNDED_ZIP_PROBE",
       json.dumps({"HEAD_bytes":head.get("content_length_header"),
          "start_http":start.get("http"),"start_content_range":start.get("content_range_header"),
          "tail_http":tail.get("http"),"tail_content_range":tail.get("content_range_header"),
          "ZIP_metadata":midata},sort_keys=True),flush=True)

if __name__=="__main__":main()
