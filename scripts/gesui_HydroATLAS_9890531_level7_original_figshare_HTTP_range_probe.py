#!/usr/bin/env python3
"""Original HydroATLAS (HydroSHEDS Figshare) shapefile ZIP Source technical
availability. NEVER download entire 3.98 GB ZIP, only HEAD and 1k+64k Range.

If not fully HTTP206 source, SOURCE STOP and do not assume land cover controls.
"""
import hashlib,json,os,pathlib,time,urllib.request,urllib.error
BASE="https://api.figshare.com/v2/articles/9890531"
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HydroATLAS_original_cropland_source"
ROOT.mkdir(parents=True,exist_ok=True)
def request(url,method="GET",ranger=None):
    heads={"User-Agent":"gesui-original-HydroATLAS-L7-covariate-source-20261112/1",
          "Accept":"application/json,application/zip,application/octet-stream,*/*"}
    if ranger:heads["Range"]=ranger
    req=urllib.request.Request(url,headers=heads,method=method)
    with urllib.request.urlopen(req,timeout=65) as res:
        sizebound=66000 if ranger=="bytes=-65536" else 4096
        raw=res.read(sizebound) if method!="HEAD" else b""
        hdr={k.lower():v for k,v in res.headers.items()}
        return {"http_status":res.status,"resolved_url":res.url,
           "content_length_header":hdr.get("content-length"),
           "content_type":hdr.get("content-type"),
           "content_range":hdr.get("content-range"),
           "accept_ranges":hdr.get("accept-ranges"),
           "bytes_sampled":len(raw),
           "sampled_SHA256":hashlib.sha256(raw).hexdigest()}
result={"schema":"HYDROATLAS_2019_FIGSHARE_OFFICIAL_L7_PARTIAL_ZIP_BYTES_SOURCE_QA_v1",
 "original_publisher_article_API":BASE,"original_article_DOI":"10.6084/m9.figshare.9890531",
 "original_file_metadata":None,"partial_access_results":{},
 "NO_FULL_4GB_ARCHIVE_REQUESTED":True,
 "no_raw_landcover_urban_or_agriculture_values_fitted":True}
try:
    req=urllib.request.Request(BASE,headers={"User-Agent":"gesui-source-only-figshare-archive-inventory/1","Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=50) as r:
        raw=r.read(400000)
    info=json.loads(raw)
    result["figshare_article_source_SHA256"]=hashlib.sha256(raw).hexdigest()
    files=[f for f in info.get("files",[]) if f.get("name")=="BasinATLAS_Data_v10_shp.zip"]
    if len(files)!=1:raise RuntimeError("Cannot select unique official HydroATLAS BasinATLAS Data v10 shapefile")
    f=files[0]
    source_url=f.get("download_url")
    if not source_url or int(f.get("size") or 0)<2000000000:
        raise RuntimeError("Published original BasinATLAS source URL/size mismatch")
    result["original_file_metadata"]={k:f.get(k) for k in ["id","name","size","download_url","supplied_md5","computed_md5","is_link_only"]}
    for label,kw in [("HEAD",("HEAD",None)),("byte_0_1023",("GET","bytes=0-1023")),("suffix_65536",("GET","bytes=-65536"))]:
        try:
            result["partial_access_results"][label]=request(source_url,method=kw[0],ranger=kw[1])
        except Exception as e:
            result["partial_access_results"][label]={"error":type(e).__name__,"detail":str(e)[:250]}
    z=result["partial_access_results"]
    result["server_byte_ranges_reliable"]=bool(
       z.get("byte_0_1023",{}).get("http_status")==206
       and z.get("suffix_65536",{}).get("http_status")==206
       and z.get("byte_0_1023",{}).get("bytes_sampled")==1024
       and z.get("suffix_65536",{}).get("bytes_sampled")==65536)
    result["status"]="PASS_TO_REMOTE_DBf_ONLY" if result["server_byte_ranges_reliable"] else "STOP_4GB_SOURCE_NOT_RANDOM_RANGE_READABLE"
except Exception as e:
    result["status"]="STOP_original_publisher_metadata_inaccessible"
    result["source_error"]=type(e).__name__+":"+str(e)[:250]
path=ROOT/"HYDROATLAS_OFFICIAL_FIGSHARE_H7_DBF_ONLY_RANGE_PREFLIGHT.json"
path.write_text(json.dumps(result,indent=2,sort_keys=True))
print("HYDROATLAS_SOURCE_RANGE_RESULT_START",flush=True)
print(json.dumps(result,sort_keys=True),flush=True)
print("HYDROATLAS_SOURCE_RANGE_RESULT_END",flush=True)
