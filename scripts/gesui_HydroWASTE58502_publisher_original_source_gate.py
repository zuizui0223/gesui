#!/usr/bin/env python3
"""HydroWASTE v1 exact original official source-gate for 58502 WWTP points.

Official HydroSHEDS project links Figshare original published binary. No
wastewater-ecosystem effect estimates; this only verifies source coordinates,
modelled flow units, HydroRIVERS IDs, missingness and checksum.
"""
import csv,hashlib,io,json,math,os,pathlib,statistics,time,urllib.request,urllib.error,zipfile
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_hydrowaste_original_source"
ROOT.mkdir(parents=True,exist_ok=True)
URLS=[
 "https://figshare.com/ndownloader/files/31910714",
 "https://ndownloader.figshare.com/files/31910714"
]
def download(url):
    last=None
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-HydroSHEDS-official-HydroWASTE-source-20261112/1",
                     "Accept":"application/zip,application/octet-stream,*/*"})
            with urllib.request.urlopen(req,timeout=70) as r:
                b=r.read(30000000)
                actual=r.url
            if not zipfile.is_zipfile(io.BytesIO(b)):
                raise RuntimeError("Original publisher link did not deliver ZIP")
            return b,actual
        except Exception as e:
            last=e;time.sleep(1+attempt*2)
    raise RuntimeError("Official HydroWASTE public source unavailable "+url+": "+str(last))
error={}
original=None
for url in URLS:
    try:
        data,actual=download(url)
        original={"source_url":url,"resolved_url":actual,"data":data}
        break
    except Exception as e:error[url]=str(e)[:300]
if original is None:
    result={"status":"STOP_OFFICIAL_HYDROWASTE_SOURCE_UNAVAILABLE","attempted_official_publisher_urls":error,
            "modelled_WWTP_count_not_assumed_as_observed":True}
    (ROOT/"HYDROWASTE_GLOBAL_ORIGINAL_SOURCE_GATE.json").write_text(json.dumps(result,indent=2))
    print("HYDROWASTE_SOURCE_GATE_RESULT_START")
    print(json.dumps(result))
    print("HYDROWASTE_SOURCE_GATE_RESULT_END")
    raise RuntimeError("Unable to verify genuine 58502 official treatment plants")
data=original["data"]
if len(data)>30000000:raise RuntimeError("Unexpected original HydroWASTE source too large")
original_sha=hashlib.sha256(data).hexdigest()
z=zipfile.ZipFile(io.BytesIO(data));z.testzip()
members=[t for t in z.infolist() if not t.is_dir()]
source_csv=[t for t in members if t.filename.lower().endswith(".csv") and ("hydrowaste" in t.filename.lower())]
if len(source_csv)!=1:raise RuntimeError("Official HydroWASTE archive CSV ambiguous "+str([x.filename for x in members]))
f=source_csv[0];raw=z.read(f)
if len(raw)>50000000:raise RuntimeError("Original global raw WWTP source unexpectedly big")
reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig"),newline=""))
original_columns=reader.fieldnames or []
cols={x.lower():x for x in original_columns}
# Preserve original fields no rename of scientific models.
req=["lon_out","lat_out","waste_dis","hyriv_id"]
if not all(x in cols for x in req):
    raise RuntimeError("Expected original HydroWASTE source fields missing "+str(original_columns))
def v(row,k):
    return row.get(cols[k],"")
def ffloat(row,k):
    try:
        x=float(v(row,k))
        return x if math.isfinite(x) else None
    except (ValueError,TypeError,KeyError):return None
n=0;unique=set();dup=0;locgood=0;outfallgood=0;flowgood=0;nonpositiveflow=0
levels={};countries={};qual_waste={}
rows=[]
for r in reader:
    n+=1
    x=ffloat(r,"lon_out");y=ffloat(r,"lat_out")
    good= x is not None and y is not None and -180<=x<=180 and -90<=y<=90
    locgood+=good
    hv=v(r,"hyriv_id")
    outfallgood+=bool(hv.strip() and hv.strip().lower() not in ("na","null","nan"))
    discharge=ffloat(r,"waste_dis")
    flowgood+=discharge is not None
    nonpositiveflow+=discharge is not None and discharge<=0
    for key,dct in [("level",levels),("country",countries),("qual_waste",qual_waste)]:
        if key in cols:
            s=v(r,key)
            dct[s]=dct.get(s,0)+1
    if "wwtp_id" in cols:
        key=v(r,"wwtp_id")
        if key in unique:dup+=1
        unique.add(key)
    # Do not commit 58502 public raw data. Artifact contains original
    # download only for source check; no biological outcome in source.
    rows.append(r)
if n!=58502:raise RuntimeError("Original global plant count expected 58502, got "+str(n))
passgate=locgood/n>=.99 and outfallgood/n>=.95 and flowgood/n>=.95 and dup==0
original_member={
 "name":f.filename,"size_bytes":len(raw),
 "source_raw_SHA256":hashlib.sha256(raw).hexdigest(),
 "original_columns":original_columns,
 "source_N_rows":n,
 "latitude_longitude_valid":locgood,
 "HYRIV_ID_present":outfallgood,
 "modeled_WWTP_discharge_numeric":flowgood,
 "modeled_WWTP_discharge_nonpositive":nonpositiveflow,
 "source_WWTP_id_duplicate":dup,
 "treatment_level_frequency":levels,
 "country_frequency_n_distinct":len(countries),
 "quality_waste_flow_frequency":qual_waste
}
report={
 "schema":"GESUI_ORIGINAL_GLOBAL_HYDROWASTE58K_SOURCE_ONLY_v1",
 "download_source":"HydroSHEDS original CC BY 4.0 license, linked Figshare publisher file",
 "source_link":original["source_url"],
 "actual_redirect_source":original["resolved_url"],
 "source_original_ZIP_SHA256":original_sha,
 "original_ZIP_size_bytes":len(data),
 "original_source_members":[{"name":x.filename,"size":x.file_size} for x in members],
 "original_CSV_original_source":original_member,
 "source_gate_pass":passgate,
 "not_measured_annual_2018_2022_sewage_NP":True,
 "not_annual_time_varying_flow":True,
 "no_biological_EQR_EPT_plant_or_butterfly_effects":True,
 "no_receiver_specific_flowpath_yet":True
}
(ROOT/"HYDROWASTE_GLOBAL_ORIGINAL_SOURCE_GATE.json").write_text(json.dumps(report,indent=2,sort_keys=True))
(ROOT/"HYDROWASTE_v10_original_publisher_archive.zip").write_bytes(data)
print("HYDROWASTE_SOURCE_GATE_RESULT_START",flush=True)
small={k:v for k,v in report.items() if k!="original_CSV_original_source"}
small["original_CSV_original_source"]={k:v for k,v in original_member.items() if k!="original_columns"}
small["original_CSV_original_source"]["original_columns"]=original_columns
print(json.dumps(small,sort_keys=True),flush=True)
print("HYDROWASTE_SOURCE_GATE_RESULT_END",flush=True)
if not passgate:raise RuntimeError("Fixed original HydroWASTE source quality gate failed; no match or slopes")
