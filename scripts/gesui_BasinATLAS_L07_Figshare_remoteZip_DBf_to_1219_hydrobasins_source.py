#!/usr/bin/env python3
"""Remote ZIP range-read only original HydroATLAS BasinATLAS Level07 DBF.
Publisher Figshare file 20087237, 4.276GB; never fetch full archive.

CRC-verify native original DBF through python zipfile and record its SHA256.
Crosswalk only the originally source-frozen 1219 HydroBASINS Level07 IDs,
without reading ecosystem biological effects.
"""
import collections,hashlib,io,json,os,pathlib,re,time,urllib.error,urllib.request,zipfile
import shapefile

ART_SHA="87ed43e0dd02e6fba1228566878b17e1ad35d9a903f1ef4d26eb9dc2a9a45e38"
ZIP_SIZE=4276492333
ZIP_ID=20087237
ZIP_MD5="ce5013294bdcb6f884ff6cdf1a1acae0"
URL=f"https://ndownloader.figshare.com/files/{ZIP_ID}"
ORIGINAL_H7_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_HydroATLAS_level07_raw_DBf"
ROOT.mkdir(parents=True,exist_ok=True)
p=pathlib.Path(os.environ["HYBAS7_FROZEN_1219_COHORT"])
if hashlib.sha256(p.read_bytes()).hexdigest()!=ORIGINAL_H7_SHA:
    raise RuntimeError("Prior 1219 monitored River Level07 source IDs SHA has changed")
meta=json.loads(p.read_text())
target_ids={str(r["HYBAS_ID_L07"]) for r in meta["candidate_source_site_period_rows"] if r["HYBAS_ID_L07"] is not None}
if len(target_ids)!=1219:raise RuntimeError("Original 1219 study basin population changed")
source_prev=json.loads(pathlib.Path(os.environ["HYDROATLAS_OFFICIAL_RANGE_QA"]).read_text())
if source_prev["status"]!="PASS_TO_REMOTE_DBf_ONLY":
    raise RuntimeError("Original exact publisher archive Range source gate did not pass")
file_meta=source_prev["original_file_metadata"]
if file_meta["id"]!=ZIP_ID or file_meta["size"]!=ZIP_SIZE or file_meta["computed_md5"]!=ZIP_MD5:
    raise RuntimeError("Original massive Figshare source ID size or MD5 changed")
class HttpRangeReader(io.RawIOBase):
    """Lightweight seekable ZIP random access; each block uses a fresh official
    Figshare public redirect, avoiding a cached expiring 10-second signed link.
    """
    def __init__(self,url,total,block=4*1024*1024):
        self.url=url;self.total=total;self.block=block
        self.position=0;self.cache_start=-1;self.cache=b"";self.calls=0;self.bytes=0
        self.max_chunk=block
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        if whence==0:nextpos=offset
        elif whence==1:nextpos=self.position+offset
        elif whence==2:nextpos=self.total+offset
        else:raise ValueError("Unsupported ZIP seek reference")
        if nextpos<0 or nextpos>self.total:raise ValueError("Original ZIP seek outside source file")
        self.position=nextpos
        return self.position
    def _fetch(self,at):
        begin=(at//self.block)*self.block
        end=min(self.total-1,begin+self.block-1)
        h={"User-Agent":"gesui-HydroATLAS-remote-ZIP-only-L07-attributes/1",
           "Accept":"application/zip,application/octet-stream,*/*",
           "Range":f"bytes={begin}-{end}"}
        err=None
        for attempt in range(4):
            try:
                req=urllib.request.Request(self.url,headers=h)
                with urllib.request.urlopen(req,timeout=70) as response:
                    if response.status!=206:raise RuntimeError("Original source refused random Access Range")
                    got=response.headers.get("Content-Range","")
                    if got!=f"bytes {begin}-{end}/{self.total}":
                        raise RuntimeError("Unexpected returned Range for source")
                    raw=response.read(end-begin+2)
                if len(raw)!=(end-begin+1):raise RuntimeError("Original Range response incomplete")
                self.cache=raw;self.cache_start=begin
                self.calls+=1;self.bytes+=len(raw)
                if self.calls%10==0:
                    print("HydroATLAS original ranges fetched",self.calls,"bytes",self.bytes,flush=True)
                return
            except Exception as e:
                err=e
                time.sleep(1+2*attempt)
        raise RuntimeError("Official HydroATLAS original random ZIP source read failed "+str(err))
    def read(self,n=-1):
        if n is None or n<0:
            # zipfile._EndRecData() uses no-arg read() only AFTER it seeks
            # to the ZIP footer (22 bytes from EOF). Safely permit just the
            # remaining footer bytes, never any 4GB unbounded retrieval.
            n=self.total-self.position
            if n>self.max_chunk:
                raise RuntimeError("No unbounded read() of 4GB official HydroATLAS archive allowed")
        if n==0:return b""
        n=min(n,self.total-self.position)
        chunks=[]
        while n>0:
            if not (self.cache_start<=self.position<self.cache_start+len(self.cache)):
                self._fetch(self.position)
            pos=self.position-self.cache_start
            take=min(len(self.cache)-pos,n)
            chunks.append(self.cache[pos:pos+take])
            self.position+=take;n-=take
        return b"".join(chunks)
r=HttpRangeReader(URL,ZIP_SIZE)
with zipfile.ZipFile(r,mode="r",allowZip64=True) as z:
    members=z.infolist()
    allnam=[{
       "name":m.filename,"original_zip_member_size":m.file_size,
       "zip_member_compressed_size":m.compress_size,"ZIP_CRC32":f"{m.CRC:08x}"
    } for m in members if not m.is_dir() and
           (m.filename.lower().endswith(".dbf") or "readme" in m.filename.lower())]
    selected=[m for m in members if not m.is_dir()
       and m.filename.lower().endswith(".dbf")
       and re.search(r"(?:^|/|\\|_)lev0?7(?:_|\.|/|\\)",m.filename.lower())]
    if len(selected)!=1:
        raise RuntimeError("Official all-level BasinATLAS Level07 DBF ambiguous; exact archived candidates "+
            json.dumps(allnam[:25]))
    dbf=selected[0]
    if dbf.file_size>500_000_000 or dbf.compress_size>350_000_000:
        raise RuntimeError("Original level7 DBF source unexpectedly huge; stop instead of retrieving full 4GB")
    print("Selected only original Level07 DBF:",dbf.filename,"source uncompressed",dbf.file_size,
          "source compressed",dbf.compress_size,flush=True)
    digest=hashlib.sha256()
    parts=[];size=0
    with z.open(dbf) as content:
        while True:
            chunk=content.read(1048576)
            if not chunk:break
            parts.append(chunk);digest.update(chunk);size+=len(chunk)
    if size!=dbf.file_size:raise RuntimeError("Original Level07 DBF downloaded length failed source check")
    raw=b"".join(parts)
    # ZipExtFile verifies the archived CRC32 at EOF, plus original member SHA
    reader=shapefile.Reader(dbf=io.BytesIO(raw),encoding="latin-1")
    fields=[x[0] for x in reader.fields[1:]]
    field_lower={x.lower():x for x in fields}
    if "hybas_id" not in field_lower:
        raise RuntimeError("Original Level07 DBF lacks real HYBAS_ID, do not join by points or alias")
    desired=["crp_pc_sse","crp_pc_use","urb_pc_sse","urb_pc_use",
            "glc_pc_s22","glc_pc_u22","dis_m3_pyr","hft_ix_s09",
            "pop_ct_s15","pop_ct_u15"]
    found={x:field_lower.get(x) for x in desired}
    if not (found["crp_pc_sse"] and found["urb_pc_sse"] and found["glc_pc_s22"]):
        raise RuntimeError("Source original BasinATLAS missing fixed cropland/artificial cover fields; "+
                           repr([x for x in fields if "crp_" in x.lower() or "glc_pc_s22" in x.lower()]))
    foundrows={}
    n_all=0
    for original in reader.iterRecords():
        n_all+=1
        d=dict(zip(fields,original))
        hid=str(d[field_lower["hybas_id"]])
        if hid in target_ids:
            if hid in foundrows:raise RuntimeError("Duplicate native BasinATLAS H7 ID source")
            foundrows[hid]={x:d[full] if full is not None else None for x,full in found.items()}
    if len(foundrows)!=1219:
        raise RuntimeError("Original BasinATLAS target H7 basin source coverage below 100%: "+
            str(len(foundrows))+"/1219; report source rather than artificially fill")
summary={
 "schema":"GESUI_ORIGINAL_FIGSHARE_BASINATLAS_L07_1219_MONITORED_BASIN_LANDUSE_SOURCE_v1",
 "source_figure_article":"10.6084/m9.figshare.9890531",
 "source_global_original_ZIP_fileID":ZIP_ID,
 "source_archive_4GB_publisher_MD5_not_independently_full_download_verified":ZIP_MD5,
 "source_original_zip_size_bytes":ZIP_SIZE,
 "original_hydrology_BasinATLAS_Level07_DBf_name":dbf.filename,
 "native_original_dbf_source_SHA256":digest.hexdigest(),
 "native_original_dbf_zip_CRC32":f"{dbf.CRC:08x}",
 "native_original_level07_DBf_row_count":n_all,
 "source_original_target_1219_L7_basins_covered":len(foundrows),
 "original_expected_source_H7_geographic_cohort_SHA256":ORIGINAL_H7_SHA,
 "source_fields":found,
 "source_zip_member_size_bytes":dbf.file_size,
 "source_zip_member_compressed_size_bytes":dbf.compress_size,
 "HTTP_range_requests_made":r.calls,"HTTP_bytes_transferred":r.bytes,
 "did_not_download_global_full_4GB_zip":True,
 "original_raw_landuse_units_and_missing_sentinels_not_modified":True,
 "raw_numeric_values_inspected_for_regression":False,
 "no_ecological_effect_calculated":True,
 "source_necessary_field_coverage":{
  x:{"nonmissing":sum(v[x] is not None for v in foundrows.values()),
      "observed_type_names":dict(collections.Counter(type(v[x]).__name__ for v in foundrows.values())),
      "raw_min":min((v[x] for v in foundrows.values() if isinstance(v[x],(int,float))),default=None),
      "raw_max":max((v[x] for v in foundrows.values() if isinstance(v[x],(int,float))),default=None)}
    for x in found},
 "original_zip_dbf_candidate_members":allnam,
}
jsonpath=ROOT/"HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json"
data={"schema":summary["schema"],"source_original_DBf_SHA256":digest.hexdigest(),
 "original_1219_H7_covariate_by_real_ID":foundrows,
 "source_values_static_not_2018_2022_anthropogenic_change":True}
payload=(json.dumps(data,sort_keys=True,separators=(",",":"))+"\n").encode()
jsonpath.write_bytes(payload)
summary["native_1219_original_covariate_manifest_SHA256"]=hashlib.sha256(payload).hexdigest()
(ROOT/"HYDROATLAS_V1_L07_1219_SOURCE_GEO_AND_CROPLAND_URBAN_QA.json").write_text(json.dumps(summary,indent=2,sort_keys=True))
print("BASINATLAS_L07_ORIGINAL_SOURCE_RESULT_START",flush=True)
short={k:v for k,v in summary.items() if k!="original_zip_dbf_candidate_members"}
print(json.dumps(short,sort_keys=True),flush=True)
print("BASINATLAS_L07_ORIGINAL_SOURCE_RESULT_END",flush=True)
