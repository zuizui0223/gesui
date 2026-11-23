#!/usr/bin/env python3
"""Authorized ORIGINAL CDSE S3 preflight for 2018+2022 10m TCD and confidence.

No passwords printed/committed. Read ONLY TIFF header range from four genuinely
published original assets for one source-frozen France EPSG3035 tile. Missing
GitHub Actions secrets -> explicit SOURCE STOP. No EEA biological outcome read.
"""
import argparse,csv,hashlib,io,json,os,pathlib,re,urllib.request
import boto3
from botocore.config import Config
from botocore.exceptions import BotoCoreError,ClientError

S3_ENDPOINT="https://eodata.dataspace.copernicus.eu"
BUCKET="eodata"
TILE="E36N28"
DATASETS={
 "TCD":("clms_vlcc_tree-cover-density_europe_10m_yearly_v1","5a6f0dc9f67ec7e3671f5984e9cf0425988e23dbb59ce1f6a83a18dc8661d6c9"),
 "TCDCL":("clms_vlcc_tree-cover-density-confidence-layer_europe_10m_yearly_v1","347053a15fe8e3cf5df7fa261daf7ae0cf9dcab4cf1fe5f81237a28f735b4222")
}
CSV_ROOT="https://s3.waw3-1.cloudferro.com/swift/v1/CatalogueCSV/landcover_landuse/tree_cover_density/"
def original_catalog(dataset,expected_sha):
 url=CSV_ROOT+dataset+"/"+dataset+"_cog.csv"
 with urllib.request.urlopen(urllib.request.Request(url,headers={
    "User-Agent":"gesui-official-CDSE-private-key-verified-original-catalogue-QA/20261124"}),timeout=60) as r:
  raw=r.read(8000000)
 if hashlib.sha256(raw).hexdigest()!=expected_sha:
  raise RuntimeError("Original public publisher S3 catalogue SHA changed")
 index={}
 for row in csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter=";"):
  for yr in ("2018","2022"):
   if f"_S{yr}_" in row["name"] and f"_{TILE}_03035_" in row["name"]:
    if yr in index:raise RuntimeError("Ambiguous original annual tile")
    index[yr]=row
 if len(index)!=2:raise RuntimeError("Original publisher source annual tile missing")
 return index
def source_key(row):
 orig=row["s3_path"]
 if not orig.startswith("s3://eodata/CLMS/landcover_landuse/tree_cover_density/"):
  raise RuntimeError("Original publisher source S3 key outside CLMS tree cover")
 tile=row["name"]
 if not re.fullmatch(r"CLMS_HRLVLCC_TCD(?:CL)?_S20(?:18|22)_R10m_E36N28_03035_V\d\d_R\d\d",tile):
  raise RuntimeError("Source publisher tile name not pinned expected 2018 or2022")
 key=orig[len("s3://eodata/"):].strip("/")+"/"+tile+".tif"
 return key
def run(args):
 reports={"schema":"GESUI_FR381_CLMS_TCD2018_2022_OFFICIAL_S3_AUTHENTICATED_SOURCE_PREFLIGHT_v1",
  "original_source_tile":TILE,"original_s3_endpoint":S3_ENDPOINT,"original_bucket":BUCKET,
  "fixed_years":[2018,2022],
  "original_canopy_raster_VALUES_read":False,
  "EEA_2018_2022_EQR_biological_outcomes_read":False,
  "temporary_credentials_written_to_disk":False,
  "secrets_never_printed":True}
 access=os.environ.get("CDSE_S3_ACCESS_KEY","")
 secret=os.environ.get("CDSE_S3_SECRET_KEY","")
 if not access or not secret:
  reports.update({"source_status":"STOP_NO_CDSE_S3_SECRET_CONFIGURED",
   "suggested_action":"Add CDSE_S3_ACCESS_KEY and CDSE_S3_SECRET_KEY in GitHub Actions repository Secrets; re-run workflow_dispatch. Never post credentials to chat or git."})
  return reports
 # Original source discovery always from official public CSV, not guessed filenames.
 try:
  catalogues={key:original_catalog(*value) for key,value in DATASETS.items()}
  cli=boto3.client("s3",endpoint_url=S3_ENDPOINT,
    aws_access_key_id=access,aws_secret_access_key=secret,
    region_name="default",config=Config(signature_version="s3v4",
     connect_timeout=15,read_timeout=35,retries={"max_attempts":2},
     s3={"addressing_style":"path"}))
  result={}
  for name,years in catalogues.items():
   result[name]={}
   for year,row in years.items():
    key=source_key(row)
    try:
     r=cli.get_object(Bucket=BUCKET,Key=key,Range="bytes=0-31")
     payload=r["Body"].read(32)
     first4=payload[:4].hex()
     native_tiff=first4 in ("49492a00","4d4d002a","49492b00","4d4d002b")
     result[name][year]={
      "original_publication_source_name":row["name"],
      "original_publisher_file_bytes":int(row["content_length"]),
      "original_publisher_content_hash":row["checksum_value"],
      "source_sha256_of_first32_returned_bytes":hashlib.sha256(payload).hexdigest(),
      "source_object_etag":r.get("ETag"),
      "public_native_TIFF_header_confirmed":native_tiff,
      "first_four_bytes_hex":first4,
      "s3_explicit_range_only":True,
      "actual_canopy_pixel_sample_values_READ":False}
    except (BotoCoreError,ClientError,TimeoutError) as e:
     status=(e.response.get("Error",{}).get("Code") if isinstance(e,ClientError) else type(e).__name__)
     result[name][year]={"source_access":"STOP", "error_type_or_code":str(status)[:90],
       "original_source_filename":row["name"],"actual_raster_values_READ":False}
  reports["original_CLMS_TCD_and_TCDCL_four_auth_source_receipts"]=result
  allpass=all(bool(result[p][yr].get("public_native_TIFF_header_confirmed"))
              for p in DATASETS for yr in ("2018","2022"))
  reports["source_status"]="PASS_CDSE_AUTHENTICATED_SOURCE_TIFF_HEADERS_ONLY_NOT_381SITE_PIXELS" if allpass else "STOP_CDSE_ORIGINAL_S3_SOURCE_ACCESS"
 except Exception as e:
  reports["source_status"]="STOP_CDSE_ORIGINAL_SOURCE_QA"
  reports["error_type"]=type(e).__name__
  reports["error_detail_nonsecret"]=str(e)[:150]
 return reports
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--out",type=pathlib.Path,required=True)
 a=p.parse_args()
 out=run(a)
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print("FR381_TCD2018_22_CDSE_S3_AUTHORIZED_PREFLIGHT",json.dumps(out,sort_keys=True),flush=True)
