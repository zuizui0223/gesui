#!/usr/bin/env python3
"""Original open Impact Observatory / Microsoft 10m Annual LULC V2 2018+2022.

Source-only independent public STAC item/TIFF actual source. No aquatic EQR,
no land change fractions, and NO raw image bytes committed. Never log a SAS
signature/token. Reject rendered RGB images rather than true native class COG.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,pathlib,urllib.parse,urllib.request,urllib.error,time

COLL="io-lulc-annual-v02"
STAC="https://planetarycomputer.microsoft.com/api/stac/v1"
SAS="https://planetarycomputer.microsoft.com/api/sas/v1/token/"+COLL
AWS_STAC="https://api.impactobservatory.com/stac-aws"
S3_BUCKET="io-10m-annual-lulc"
EXAMPLE_BBOX=(2.0,48.4,3.0,49.0)
YEAR=(2018,2022)
USER_AGENT="gesui 2018/2022 original IO annual LULC source-only QA/20261124"
def get(url,cap=1500000,range=None):
 hdr={"User-Agent":USER_AGENT,"Accept":"application/json,application/octet-stream,*/*"}
 if range:hdr["Range"]=range
 req=urllib.request.Request(url,headers=hdr)
 try:
  with urllib.request.urlopen(req,timeout=27) as r:
   body=r.read(cap+1)
   u=urllib.parse.urlparse(r.url)
   clean=urllib.parse.urlunparse((u.scheme,u.netloc,u.path,"","",""))
   return {"status":"HTTP_SUCCESS","code":r.status,"clean_url":clean,
      "type":r.headers.get("Content-Type"),
      "content_range":r.headers.get("Content-Range"),
      "body_size":len(body),"body_capped":len(body)>cap,
      "body_sha256":hashlib.sha256(body[:cap]).hexdigest(),"_bytes":body[:cap]}
 except urllib.error.HTTPError as e:
  return {"status":"HTTP_ERROR","code":e.code,"error_type":"HTTPError"}
 except Exception as e:
  return {"status":"NETWORK_ERROR","error_type":type(e).__name__,"detail":str(e)[:130]}
def public(r):
 return {k:v for k,v in r.items() if k not in ("_bytes","_json")}
def jsonget(url):
 r=get(url,cap=3000000)
 if r["status"]=="HTTP_SUCCESS":
  try:r["_json"]=json.loads(r["_bytes"])
  except ValueError:r["parse_error"]="not valid JSON"
 return r
def query(q,base=STAC):
 return base+"/search?"+urllib.parse.urlencode(q)
def source_item_features(year):
 p={"collections":COLL,"bbox":",".join(map(str,EXAMPLE_BBOX)),
    "datetime":f"{year}-01-01T00:00:00Z/{year}-12-31T23:59:59Z",
    "limit":"20"}
 r=jsonget(query(p))
 return r
def try_tiff(href,token):
 if not href:return {"status":"NO_ASSET"}
 candidates=[("original_STAC_asset",href)]
 if href.startswith("https://") and token:
  sep="&" if "?" in href else "?"
  candidates.append(("official_public_SAS",href+sep+token))
 if href.startswith("s3://"+S3_BUCKET+"/"):
  candidates.append(("original_open_AWS_S3","https://"+S3_BUCKET+".s3.us-west-2.amazonaws.com/"+href[len("s3://"+S3_BUCKET+"/"):]))
 output=[]
 for mode,url in candidates:
  r=get(url,cap=64,range="bytes=0-63")
  magic=r.get("_bytes",b"")[:4].hex()
  valid=magic in ("49492a00","4d4d002a","49492b00","4d4d002b")
  output.append({"mode":mode,"response":public(r),"first4_original_tiff_magic":magic,
                 "genuine_native_TIFF_header":valid})
  if valid:break
 return output
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--out",type=pathlib.Path,required=True)
 a=p.parse_args()
 coll=jsonget(STAC+"/collections/"+COLL)
 old=jsonget(AWS_STAC+"/collections")
 collection=coll.get("_json") or {}
 if coll["status"]!="HTTP_SUCCESS" or collection.get("id")!=COLL:
  raise RuntimeError("Official open original Microsoft Planetary Computer STAC collection not authenticated")
 if collection.get("license") not in ("CC-BY-4.0","CC-BY-4.0-ESRI","CC-BY-4.0-ImpactObservatory"):
  # Field may be "CC-BY-4.0"; fail the source contract rather than broadening
  raise RuntimeError("Original publisher expected CC BY 4.0 license not confirmed: "+str(collection.get("license")))
 src_codes=[]
 for asset in (collection.get("item_assets") or {}).values():
  if isinstance(asset,dict) and asset.get("file:values"):
   src_codes+=asset["file:values"]
 # SAS token is a short-lived **official public** read credential, not a secret
 # from user. NEVER serialize token or signed query string in artifacts.
 tokenr=jsonget(SAS)
 token=(tokenr.get("_json") or {}).get("token")
 if not isinstance(token,str):token=""
 features={}
 for year in YEAR:
  q=source_item_features(year)
  items=(q.get("_json") or {}).get("features",[])
  # This STAC service can return earlier item IDs for a requested annual
  # datetime window because metadata time extents differ from item YEAR.
  # Record that discrepancy; NEVER treat a 2017 tile as a 2018 pixel.
  mismatched=[i.get("id") for i in items if not str(i.get("id","")).endswith("-"+str(year))]
  matched=[i for i in items if str(i.get("id","")).endswith("-"+str(year))]
  if not matched and items:
   # Reuse the actual returned UTM/MGRS supercell prefix for the exact-year
   # SOURCE item lookup; require real HTTP200 with exact original ID.
   prefix=str(items[0].get("id","")).rsplit("-",1)[0]
   if not __import__("re").fullmatch(r"[0-9]{1,2}[C-X]",prefix):
    raise RuntimeError("Original STAC source item MGRS ID invalid "+prefix[:30])
   expected=prefix+"-"+str(year)
   direct=jsonget(STAC+"/collections/"+COLL+"/items/"+expected)
   item=direct.get("_json") or {}
   if direct["status"]=="HTTP_SUCCESS" and item.get("id")==expected:
    matched=[item]
  data=[]
  for item in matched[:3]:
   itemdate=item.get("properties",{}).get("datetime")
   ident=item.get("id","")
   attrs=item.get("assets",{})
   target=attrs.get("data",{})
   href=target.get("href")
   if not ident.endswith("-"+str(year)):
    raise RuntimeError("Original exact year item identity mismatch")
   content={
    "year":year,"item_id":ident,"collection":item.get("collection"),
    "source_item_datetime":itemdate,"source_grid_epsg":item.get("properties",{}).get("proj:epsg"),
    "original_bbox":item.get("bbox"),
    "native_asset_type":target.get("type"),"original_asset_roles":target.get("roles"),
    "native_asset_public_path_no_query":urllib.parse.urlparse(href or "").path,
    "native_asset_public_hostname":urllib.parse.urlparse(href or "").netloc,
    "native_asset_public_scheme":urllib.parse.urlparse(href or "").scheme,
    "actual_tiff_http_checks":try_tiff(href,token) if len(data)<2 else []}
   data.append(content)
  features[str(year)]={"source_http":public(q),"original_item_count_in_example_bbox":len(items),
                        "year_mismatched_item_count_in_first20":len(mismatched),
                        "first_wrong_year_original_item_ids":mismatched[:4],
                        "original_verified_matching_year_source_items":data,
                        "original_native_2018_2022_year_asset_go":any(
       any(v.get("genuine_native_TIFF_header") for v in x["actual_tiff_http_checks"]) for x in data)}
 genuine=all(features[str(year)]["original_native_2018_2022_year_asset_go"] for year in YEAR)
 results={
  "schema":"GESUI_FRANCE_IO_LULC_V2_ANNUAL2018_2022_ORIGINAL_STAC_REAL_TIFF_SOURCE_v1",
  "original_collection":COLL,"source_publisher":"Impact Observatory, Esri and Microsoft",
  "actual_source_collection_public_HTTP":public(coll),
  "original_collection_id":collection.get("id"),"source_license":collection.get("license"),
  "source_title":collection.get("title"),"source_temporal_extent":collection.get("extent",{}).get("temporal"),
  "source_original_class_codes":src_codes[:30],
  "original_AWS_STAC_discovery":public(old),
  "genuine_public_SAS_endpoint_responded":tokenr.get("status")=="HTTP_SUCCESS" and bool(token),
  "genuine_public_SAS_endpoint_code":tokenr.get("code"),
  "original_source_example_bbox":EXAMPLE_BBOX,
  "source_2018_2022_native_asset_source":features,
  "original_2018_and_2022_native_COG_TIFF_headers_publicly_readable":genuine,
  "original_EEA_benthos_EQR_or_annual_vegetation_values_READ":False,
  "source_go_for_381_France_before_pixel_support":genuine,
  "source_status":"SOURCE_PUBLIC_NATIVE_COG2018_2022_GO_NEED_PIXEL_SITE_QUALITY" if genuine else "SOURCE_STOP_NO_PUBLIC_NATIVE_BOTH_YEAR_COG_TIFF"
 }
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
 print("IO_LULC_V2_FRANCE_2018_2022_REAL_COG_SOURCE",json.dumps({
   "collection":results["source_title"],"license":results["source_license"],
   "SAS":results["genuine_public_SAS_endpoint_responded"],
   "source_years":features,
   "source_status":results["source_status"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
