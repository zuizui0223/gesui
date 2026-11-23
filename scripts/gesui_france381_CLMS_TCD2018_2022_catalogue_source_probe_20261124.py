#!/usr/bin/env python3
"""Original CLMS TCD annual 2018 & 2022 catalogue/pixel access SOURCE ONLY.

Source discovery and public access receipts: no EEA biological observations,
no EQR, no vegetation numeric rasters read. Truthful GO/STOP; avoid guessing
TCD asset URLs or creating a surrogate vegetation change series.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,pathlib,time,urllib.request,urllib.error,csv,io,re
from urllib.parse import urlencode
IDENTIFIER="clms_vlcc_tree-cover-density_europe_10m_yearly_v1"
URLS={
 "official_clms_product":"https://land.copernicus.eu/en/products/high-resolution-layer-forests-and-tree-cover/tree-cover-density-2018-present-raster-10-m-europe-yearly",
 "official_cdse_stac_collections":"https://stac.dataspace.copernicus.eu/v1/collections",
 "official_cdse_csv_directory":"https://csv.dataspace.copernicus.eu/CLMS/landcover_landuse/tree_cover_density/"+IDENTIFIER+"/",
 "official_cdse_csv_cog_inventory":"https://s3.waw3-1.cloudferro.com/swift/v1/CatalogueCSV/landcover_landuse/tree_cover_density/"+IDENTIFIER+"/"+IDENTIFIER+"_cog.csv",
 "official_cdse_odata_one_record":"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?"+urlencode({
  "$count":"true","$top":"1",
  "$filter":"Collection/Name eq 'CLMS' and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'datasetIdentifier' and att/OData.CSC.StringAttribute/Value eq '"+IDENTIFIER+"')"
 })}
def get(url,limit,range_header=None,request_method="GET"):
 headers={"User-Agent":"gesui-ecology-research-original-CLMS-TCD-source-QA/20261124",
          "Accept":"application/json,text/csv,text/plain,*/*"}
 if range_header:headers["Range"]=range_header
 request=urllib.request.Request(url,headers=headers,method=request_method)
 try:
  with urllib.request.urlopen(request,timeout=35) as r:
   b=r.read(limit+1) if request_method=="GET" else b""
   return {"status":"HTTP_SUCCESS","code":r.status,"effective_url":r.url,
    "content_type":r.headers.get("Content-Type"),
    "content_length":r.headers.get("Content-Length"),
    "content_range":r.headers.get("Content-Range"),
    "accept_ranges":r.headers.get("Accept-Ranges"),
    "etag":r.headers.get("ETag"),
    "last_modified":r.headers.get("Last-Modified"),
    "retrieved_bytes":len(b),"capped":len(b)>limit,
    "sha256_of_retrieved_bytes":hashlib.sha256(b[:limit]).hexdigest(),
    "_body":b[:limit]}
 except urllib.error.HTTPError as e:
  return {"status":"HTTP_ERROR","code":e.code,"reason":str(e)[:220]}
 except Exception as e:
  return {"status":"NETWORK_ERROR","reason":type(e).__name__+": "+str(e)[:250]}
def safe(j):
 return {k:v for k,v in j.items() if k!="_body"}
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--out",type=pathlib.Path,required=True)
 a=ap.parse_args()
 receipts={}
 for name in ("official_cdse_csv_cog_inventory","official_cdse_stac_collections","official_cdse_odata_one_record"):
  url=URLS[name]
  if name=="official_cdse_csv_cog_inventory":
   receipts[name+"_HEAD"]=get(url,0,request_method="HEAD")
   receipts[name+"_RANGE"]=get(url,192000,range_header="bytes=0-191999")
  elif name=="official_cdse_stac_collections":
   receipts[name]=get(url,3000000)
  else:
   receipts[name]=get(url,64000)
 observations={}
 c=receipts["official_cdse_csv_cog_inventory_RANGE"]
 if c["status"]=="HTTP_SUCCESS":
  text=c["_body"].decode("utf-8-sig",errors="replace")
  lines=text.splitlines()
  columns=lines[0].split(",") if lines else []
  candidate=[str(x) for x in columns if any(y in x.lower() for y in (
    "year","name","date","href","url","s3","bbox","filename","path","time","geometr"))]
  counts={"2018":0,"2022":0}
  matches=[]
  for line in lines[1:]:
   for yr in counts:
    if re.search(r"(?<!\d)"+yr+r"(?!\d)",line):
     counts[yr]+=1
     if len(matches)<5:matches.append({"year":yr,"raw_source_catalog_line_max_180chars":line[:180]})
  observations["csv_header_n_cols"]=len(columns)
  observations["csv_headers"]=columns[:65]
  observations["csv_source_keys_of_interest"]=candidate[:65]
  observations["preview_lines"]=len(lines)
  observations["preview_2018_2022_string_hits"]=counts
  observations["preview_selected_year_rows"]=matches
  observations["range_http_partial_206"]=c.get("code")==206
  observations["source_complete_inventory_not_yet_verified"]=True
 st=receipts["official_cdse_stac_collections"]
 if st["status"]=="HTTP_SUCCESS":
  try:
   data=json.loads(st["_body"])
   collections=data.get("collections",[])
   observations["stac_collections_found"]=len(collections)
   observations["stac_total_collections_source_truncated"]=bool(st["capped"])
   matches=[{"id":x.get("id"),"title":x.get("title"),
             "time":x.get("extent",{}).get("temporal",{}).get("interval")}
    for x in collections if any(k in (str(x.get("id"))+" "+str(x.get("title"))).lower()
                             for k in ("tree-cover-density","tcd","tree cover density"))]
   observations["stac_tcd_collection_matches"]=matches[:50]
  except Exception as e:observations["stac_json_parse_error"]=type(e).__name__+":"+str(e)[:100]
 od=receipts["official_cdse_odata_one_record"]
 if od["status"]=="HTTP_SUCCESS":
  try:
   data=json.loads(od["_body"])
   observations["odata_provider_count"]=data.get("@odata.count")
   observations["odata_one_item_name"]=(data.get("value") or [{}])[0].get("Name")
   observations["odata_provider_query_tcd_exact_filter"]=True
  except Exception as e:observations["odata_json_parse_error"]=type(e).__name__
 output={
  "schema":"gesui_Copernicus_CLMS_TCD2018_TCD2022_France381_catalogue_source_probe_v1",
  "timestamp_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
  "original_publisher_TCD_DOI":"10.2909/e677441e-fb94-431c-b4f9-304f10e4dfd8",
  "fixed_years":[2018,2022],"fixed_original_pixel_size_m":10,
  "publisher_dataset_id":IDENTIFIER,
  "source_only_does_NOT_read_2018_or_2022_EQR_or_TCD_values":True,
  "official_endpoints":URLS,
  "source_receipts":{k:safe(v) for k,v in receipts.items()},
  "source_metadata":observations,
  "original_two_year_TCD_tile_pixel_data_proven":False,
  "original_tcd_annual_quality_raster_data_proven":False,
  "source_status":"SOURCE_CATALOGUE_DISCOVERY_ONLY_NOT_ORIGINAL_PIXEL_GO",
  "reason":"Only original CLMS catalogue HEAD/range/metadata. Native 2018 and 2022 TCD GeoTIFF + quality band with actual pixels and same 381 France original site source ring support required before ECOLOGY."}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("CLMS_TCD2018_2022_ORIGINAL_SOURCE_CATALOGUE",json.dumps({
  "status":output["source_status"],
  "source_receipts":output["source_receipts"],
  "metadata":observations},sort_keys=True),flush=True)
if __name__=="__main__":main()
