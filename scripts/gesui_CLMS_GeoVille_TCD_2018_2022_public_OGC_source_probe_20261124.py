#!/usr/bin/env python3
"""Source-only CLMS GeoVille public OGC WCS/WMS numeric original TCD opportunity.

No aquatic EQR and no original S3/CDSE authentication bypass; uses openly
advertised WMS/WCS endpoints. Do not confuse display maps with actual TCD values.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,pathlib,urllib.parse,urllib.request,urllib.error,xml.etree.ElementTree as ET,re,math,concurrent.futures
BASE="https://geoserver.vlcc.geoville.com/geoserver/ows"
LAYER={2018:"HRL_TCF:TCD_S2018",2022:"HRL_TCF:TCD_S2022"}
def source(url,max_bytes=3000000):
 request=urllib.request.Request(url,headers={
   "Accept":"application/xml,text/xml,application/json,text/plain,*/*",
   "User-Agent":"gesui-official-CLMS-public-OGC-source-only/20261124"})
 try:
  with urllib.request.urlopen(request,timeout=12) as r:
   b=r.read(max_bytes+1)
   return {"status":"HTTP_SUCCESS","code":r.status,"content_type":r.headers.get("Content-Type"),
    "response_bytes_capped":len(b)>max_bytes,"response_body_bytes":min(len(b),max_bytes),
    "response_SHA256":hashlib.sha256(b[:max_bytes]).hexdigest(),
    "effective_url":r.url,"text":b[:max_bytes]}
 except urllib.error.HTTPError as e:
  body=e.read(1200)
  return {"status":"HTTP_ERROR","code":e.code,"reason":str(e)[:200],
    "response_body_sample":body.decode("utf-8",errors="replace")[:350]}
 except Exception as e:
  return {"status":"NETWORK_ERROR","reason":type(e).__name__+":"+str(e)[:300]}
def slim(r):
 return {k:v for k,v in r.items() if k!="text"}
def req(params):
 return BASE+"?"+urllib.parse.urlencode(params)
def parse_xml(r):
 if r["status"]!="HTTP_SUCCESS":return None
 try:return ET.fromstring(r["text"])
 except Exception as e:return {"parse_error":type(e).__name__+":"+str(e)[:140]}
def descendants_names(root):
 return [e.text.strip() for e in root.iter()
         if e.tag.rsplit("}",1)[-1] in ("Name","CoverageId","Identifier") and e.text]
def diagnostic(root,kind):
 if root is None:return {"ok":False}
 if isinstance(root,dict):return root
 ids=descendants_names(root)
 hits=[x for x in ids if any(k in x.upper() for k in ("TCD_S2018","TCD_S2022","TCDCL","CONFIDENCE"))]
 formats=[e.text.strip() for e in root.iter()
          if e.tag.rsplit("}",1)[-1] in ("Format","SupportedFormat") and e.text]
 return {"root_tag":root.tag.rsplit("}",1)[-1],
  "original_supported_source_names_found":hits[:35],
  "source_named_TCD_2018":any("TCD_S2018" in q for q in hits),
  "source_named_TCD_2022":any("TCD_S2022" in q for q in hits),
  "source_format_capability_preview":list(dict.fromkeys(formats))[:30],
  "all_names_count":len(ids),"full_response_not_truncated":True}
def getfeatureinfo(year):
 # Study-independent original public source point in French mainland, not a
 # site selected on ecological outcomes (longitude=0.9 latitude=48.2).
 lon,lat=0.9,48.2
 x=6378137*math.radians(lon);y=6378137*math.log(math.tan(math.pi/4+math.radians(lat)/2))
 url=req({"service":"WMS","version":"1.1.1","request":"GetFeatureInfo",
  "layers":LAYER[year],"query_layers":LAYER[year],
  "srs":"EPSG:3857","bbox":f"{x-256},{y-256},{x+256},{y+256}",
  "width":101,"height":101,"x":50,"y":50,
  "styles":"","format":"image/png","info_format":"text/plain","feature_count":10})
 v=source(url,max_bytes=8000)
 result=slim(v)
 if v["status"]=="HTTP_SUCCESS":
  result["body_first_350chars"]=v["text"].decode("utf-8",errors="replace")[:350]
  result["contains_observed_percent_or_raw_band_numeric"]=bool(
    re.search(r"(?i)(?:TCD|VALUE|GRAY_INDEX|BAND_1|RASTER)\s*[:=]\s*[\d.]+",
     result["body_first_350chars"]))
 return result
def main():
 parser=argparse.ArgumentParser()
 parser.add_argument("--out",required=True,type=pathlib.Path)
 a=parser.parse_args()
 receipts={}
 parsed={}
 for service,version in (("WMS","1.3.0"),("WCS","2.0.1"),("WCS","1.0.0")):
  key=service+"_"+version.replace(".","_")
  v=source(req({"service":service,"version":version,"request":"GetCapabilities"}),
           max_bytes=5000000)
  receipts[key]=slim(v)
  root=parse_xml(v)
  parsed[key]=diagnostic(root,key)
  if v["status"]=="HTTP_SUCCESS" and "text" in v and v["text"]:
   receipts[key]["first_response_130chars"]=v["text"][:130].decode("utf-8",errors="replace")
 probe={year:getfeatureinfo(year) for year in LAYER}
 tcd_years=any(q["source_named_TCD_2018"] and q["source_named_TCD_2022"]
       for k,q in parsed.items() if k.startswith("WMS") and "source_named_TCD_2018" in q)
 wcs_years=any(q["source_named_TCD_2018"] and q["source_named_TCD_2022"]
       for k,q in parsed.items() if k.startswith("WCS") and "source_named_TCD_2018" in q)
 ga={str(y):v for y,v in probe.items()}
 result={
  "schema":"GESUI_FRANCE_PUBLIC_GEOVILLE_WMS_WCS_TCD2018_2022_SOURCE_ONLY_v1",
  "requested_at_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
  "public_original_publisher_source_layer_preview":"https://geoserver.vlcc.geoville.com/geoserver/web/wicket/bookmarkable/org.geoserver.web.demo.MapPreviewPage?0=&filter=false",
  "public_OGC_endpoint":BASE,
  "candidate_original_layers":{str(k):v for k,v in LAYER.items()},
  "original_OGC_capability_transport":receipts,
  "genuine_getcapabilities_parsed":parsed,
  "fixed_france_original_nonbiological_source_GetFeatureInfo":ga,
  "original_public_layer_both_years_WMS_support":tcd_years,
  "original_public_layer_both_years_WCS_numerical_coverage_ADVERTISED":wcs_years,
  "real_original_unstyled_numeric_TCD_raster_access_verified":False,
  "both_years_native_TCD_confidence_layer_numeric_access_verified":False,
  "actual_2018_and_2022_FRANCE381_TCD_pixel_change_CALCULATED":False,
  "source_status":"DISCOVERY_ONLY_NO_GEO_TIFF_QUALITY_PIXEL_GO",
  "source_wms_colormap_not_a_percentage_substitute":True,
  "no_EEA_river_biological_outcomes_opened":True}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
 print("CLMS_PUBLIC_GEOVILLE_2018_2022_OGC_SOURCE",json.dumps({
   "transport":receipts,"capabilities":parsed,
   "getfeatureinfo":ga,
   "WMS_2018_2022_both":tcd_years,"WCS_2018_2022_both":wcs_years,
   "data_source_status":result["source_status"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
