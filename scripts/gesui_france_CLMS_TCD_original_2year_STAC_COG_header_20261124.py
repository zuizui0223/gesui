#!/usr/bin/env python3
"""Official Copernicus annual TCD2018 and 2022 COG sample tile pixel access.

Ground source item names for French EPSG3035 tile E36N28 from complete original
provider TCD CSV, verify exact STAC metadata and anonymous HTTP raster byte
range, not numerical canopy pixel or EEA biological response. If authentication
required, explicitly SOURCE STOP rather than claiming change measured.
"""
from __future__ import annotations
import argparse,csv,datetime,hashlib,io,json,pathlib,urllib.request,urllib.error,urllib.parse,re
IDENT="clms_vlcc_tree-cover-density_europe_10m_yearly_v1"
CSV="https://s3.waw3-1.cloudferro.com/swift/v1/CatalogueCSV/landcover_landuse/tree_cover_density/"+IDENT+"/"+IDENT+"_cog.csv"
S3_ENDPOINT="https://eodata.dataspace.copernicus.eu/"
STAC="https://stac.dataspace.copernicus.eu/v1"
TILE="E36N28"
EXPECTED_CSV_SHA="5a6f0dc9f67ec7e3671f5984e9cf0425988e23dbb59ce1f6a83a18dc8661d6c9"
def fetch(url,limit=180000,headers=None,method="GET"):
 hdr={"User-Agent":"gesui official annual CLMS TCD original pixel access pilot/20261124",
      "Accept":"application/json,application/octet-stream,*/*"}
 hdr.update(headers or {})
 try:
  req=urllib.request.Request(url,headers=hdr,method=method)
  with urllib.request.urlopen(req,timeout=40) as f:
   body=f.read(limit+1) if method=="GET" else b""
   return {"status":"HTTP_SUCCESS","code":f.status,"final_url":f.url,
     "content_type":f.headers.get("Content-Type"),
     "content_length":f.headers.get("Content-Length"),
     "content_range":f.headers.get("Content-Range"),
     "etag":f.headers.get("ETag"),
     "sample_bytes":len(body),
     "bytes_over_limit":len(body)>limit,
     "first8_hex":body[:8].hex(),
     "sha256_fetched":hashlib.sha256(body[:limit]).hexdigest(),
     "_body":body[:limit]}
 except urllib.error.HTTPError as e:
  return {"status":"HTTP_ERROR","code":e.code,"error":str(e)[:160],
      "www_authenticate":e.headers.get("WWW-Authenticate"),
      "content_type":e.headers.get("Content-Type")}
 except Exception as e:
  return {"status":"NETWORK_ERROR","error":type(e).__name__+":"+str(e)[:200]}
def slim(x):return {k:v for k,v in x.items() if k!="_body"}
def data(url):
 r=fetch(url)
 try:
  if r["status"]=="HTTP_SUCCESS":r["_json"]=json.loads(r["_body"])
 except Exception as e:r["json_error"]=str(e)[:120]
 return r
def uri_to_public_s3(href):
 if href.startswith("s3://"):
  return S3_ENDPOINT+href[len("s3://"):]
 return href if href.startswith("https://") else None
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--out",type=pathlib.Path,required=True)
 a=p.parse_args()
 r=fetch(CSV,limit=3100000)
 if r["status"]!="HTTP_SUCCESS" or hashlib.sha256(r["_body"]).hexdigest()!=EXPECTED_CSV_SHA:
  raise RuntimeError("CLMS full original source catalogue SHA changed")
 reader=csv.DictReader(io.StringIO(r["_body"].decode("utf-8-sig")),delimiter=";")
 rows={}
 for t in reader:
  if (TILE+"_03035_" in t["name"] and
      ("_S2018_" in t["name"] or "_S2022_" in t["name"])):
   yr="2018" if "_S2018_" in t["name"] else "2022"
   if yr in rows:raise RuntimeError("Duplicate original annual TCD French example tile")
   rows[yr]=t
 if set(rows)!={"2018","2022"}:raise RuntimeError("Real original matched French tile pair absent")
 collections=[IDENT+"_cog",IDENT]
 probes={}
 try:
  # Official STAC first page is paginated; we probe exact ID, not assume missing
  # simply because TCD absent from the first page.
  for coll in collections:
   url=STAC+"/collections/"+urllib.parse.quote(coll,safe="")
   v=data(url)
   probes[coll]={"collection":slim(v),"declared_id":(v.get("_json") or {}).get("id"),
      "item_assets":list((v.get("_json") or {}).get("item_assets",{}).keys())[:30]}
 except Exception as e:
  probes["probe_error"]=type(e).__name__+":"+str(e)[:150]
 valid_coll=next((k for k,v in probes.items() if isinstance(v,dict) and v.get("declared_id")==k),None)
 items={}
 assets=[]
 for yr,t in rows.items():
  name=t["name"]
  sourceitem={"original_name":name,"official_product_UUID":t["id"],
   "official_published_COG_zip_bytes":int(t["content_length"]),
   "original_publisher_MD5":t["checksum_value"],
   "original_s3_path":t["s3_path"],
   "original_bbox":t["bbox"]}
  if valid_coll:
   u=STAC+"/collections/"+urllib.parse.quote(valid_coll,safe="")+"/items/"+urllib.parse.quote(name,safe="")
   result=data(u)
   j=result.get("_json") or {}
   sourceitem["item_STAC_response"]=slim(result)
   sourceitem["item_STAC_ID"]=j.get("id")
   sourceitem["item_STAC_collection"]=j.get("collection")
   sourceitem["item_STAC_geometry_type"]=(j.get("geometry") or {}).get("type")
   for key,v in (j.get("assets") or {}).items():
    if not isinstance(v,dict):continue
    href=v.get("href","")
    candidate={"asset_key":key,"stac_href":href[:750],
     "stac_media_type":v.get("type"),
     "roles":v.get("roles"),"auth":v.get("auth:refs"),
     "storage":v.get("storage:refs"),
     "file_type":str(v.get("title",""))[:150]}
    sourceitem.setdefault("source_assets",[]).append(candidate)
    if ("cog" in str(v.get("type","")).lower() or
        "tif" in str(v.get("type","")).lower() or
        "tree cover density" in key.lower() or key.lower() in ("data","tcd")):
     public_uri=uri_to_public_s3(href)
     if public_uri and len(assets)<4:
      assets.append({"year":yr,"asset":candidate,"original_public_uri_to_test":public_uri})
  # Original provider OData service: verify only product file access headers;
  # no guessed TIFF path derived from the ZIP S3 prefix.
  od="https://catalogue.dataspace.copernicus.eu/odata/v1/Products("+t["id"]+")"
  p_od=data(od+"?$expand=Assets")
  body=p_od.get("_json") or {}
  sourceitem["odata_exact_source_product"]=slim(p_od)
  sourceitem["odata_exact_source_key"]=body.get("Id")
  sourceitem["odata_original_assets_names"]=[
   {"Name":e.get("Name"),"S3Path":str(e.get("S3Path",""))[:250],"ContentType":e.get("ContentType")}
    for e in (body.get("Assets") or []) if isinstance(e,dict)][:12]
  items[yr]=sourceitem
 checks=[]
 for a0 in assets[:4]:
  out={"year":a0["year"],"source_asset":a0["asset"],
       "uri":a0["original_public_uri_to_test"]}
  v=fetch(out["uri"],limit=1024,headers={"Range":"bytes=0-1023"})
  out["raw_byte_range"]=slim(v)
  out["first_bytes_are_TIFF"]=v["status"]=="HTTP_SUCCESS" and v.get("first8_hex","")[:8] in ("49492a00","4d4d002a","49492b00","4d4d002b")
  checks.append(out)
 passed_by_year={yr:any(i["year"]==yr and i["first_bytes_are_TIFF"] for i in checks) for yr in ("2018","2022")}
 output={
  "schema":"gesui_official_TCD2018_2022_France_real_COG_anonymous_pixel_access_pilot_v1",
  "official_CLMS_TCD10m_DOI":"10.2909/e677441e-fb94-431c-b4f9-304f10e4dfd8",
  "official_original_catalogue_SHA":EXPECTED_CSV_SHA,
  "original_sample_France_EPSG3035_100km_tile":TILE,
  "source_original_STAC_exact_collection_probes":probes,
  "source_original_items_2018_2022":items,
  "source_official_COG_data_byte_range_checks":checks,
  "2018_and_2022_genuine_anonymous_raster_TIFF_access":passed_by_year,
  "both_years_real_TIFF_source_pass":all(passed_by_year.values()),
  "source_status":"PASS_BOTH_YEAR_REAL_TIFF_HEADERS_NOT_YET_SITE_PIXEL_QA" if all(passed_by_year.values()) else "STOP_NO_PUBLIC_SOURCE_TIFF_FOR_BOTH_YEARS",
  "biological_EEA2018_2022_outcomes_READ":False,
  "TCD_original_pixel_values_READ":False,
  "next":"Require legitimate individual 381 France source sites, 10m pixel in both annuli, original TCD quality flag; no source change inference before those pass."}
 p=a.out;p.parent.mkdir(parents=True,exist_ok=True)
 p.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("CLMS_TCD2018_2022_FR_REAL_PIXEL_ACCESS",json.dumps({
  "source_collection":valid_coll,
  "original_STAC_collection_receipts":{k:v.get("collection",{}).get("code") for k,v in probes.items() if isinstance(v,dict)},
  "two_year_tile_names":{k:v["original_name"] for k,v in items.items()},
  "assets":{k:v.get("source_assets") for k,v in items.items()},
  "head_checks":checks,
  "source_status":output["source_status"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
