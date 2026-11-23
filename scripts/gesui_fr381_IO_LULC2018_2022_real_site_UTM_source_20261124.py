#!/usr/bin/env python3
"""Source-only geospatial pairing of original France381 EEA/RZ site universe
with real IO LULC 9-class V2 2018+2022 yearly STAC supercells.

2018 and 2022 response EQR *numeric columns are NEVER read*. Original real
EEA France381 site IDs remain private on runner; only aggregated geo receipts
and unchanged source folds are saved in public GitHub artifact.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,sys,urllib.parse,urllib.request,re,math
import numpy as np
sys.path.insert(0,"scripts")
import gesui_FR381_TCD2018_2022_quality_tile_matched_source_20261124 as prior

COLL="io-lulc-annual-v02"
STAC="https://planetarycomputer.microsoft.com/api/stac/v1"
BBOX=(-5.5,41.0,10.1,51.6)
PAD_LAT=0.0047
PAD_LON=0.009
YEARS=(2018,2022)
def httpjson(url):
 req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui LULC 2018 2022 original France381 spatial source audit/20261124"})
 with urllib.request.urlopen(req,timeout=35) as r:
  data=json.loads(r.read(5000000))
  if r.status!=200:raise RuntimeError("Official public STAC HTTP bad")
  return data
def source_items():
 url=STAC+"/search?"+urllib.parse.urlencode({
  "collections":COLL,"bbox":",".join(map(str,BBOX)),
  "datetime":"2018-01-01T00:00:00Z/2018-12-31T23:59:59Z","limit":"150"})
 data=httpjson(url)
 features=data.get("features",[])
 # Avoid a source paginated result silently missing real site classes.
 if not features:raise RuntimeError("Original STAC France land-use source empty")
 if any(x.get("rel")=="next" for x in data.get("links",[])):
  raise RuntimeError("Original source STAC result paginated at limit150: cannot assure all France tiles")
 original_ids=sorted(x.get("id","") for x in features)
 prefixes=sorted({s.rsplit("-",1)[0] for s in original_ids
                  if re.fullmatch(r"\d{1,2}[C-X]-20\d\d",s)})
 if len(prefixes)<3:raise RuntimeError("Original French STAC grid MGRS supercells too few")
 by_year={}
 for prefix in prefixes:
  by_year[prefix]={}
  for y in YEARS:
   originalid=prefix+"-"+str(y)
   try:
    item=httpjson(STAC+"/collections/"+COLL+"/items/"+originalid)
   except Exception as exc:
    by_year[prefix][str(y)]={"status":"SOURCE_ITEM_NOT_ACCESSIBLE","error_type":type(exc).__name__}
    continue
   if item.get("id")!=originalid or item.get("collection")!=COLL:
    raise RuntimeError("Original source per-year item ID mismatch")
   dat=item.get("assets",{}).get("data",{})
   if "tiff" not in str(dat.get("type","")).lower():
    raise RuntimeError("Original year geo item not TIFF")
   b=item.get("bbox")
   if not b or len(b)!=4:
    raise RuntimeError("Source original STAC item WGS bbox missing")
   epsg=(item.get("properties") or {}).get("proj:epsg")
   if epsg is None or not str(epsg).startswith("326"):
    raise RuntimeError("Original IO annual source original UTM grid not EPSG326")
   link=dat.get("href","")
   if not link.startswith("https://") or "io-lulc" not in link:
    raise RuntimeError("Original scientific COG href missing")
   by_year[prefix][str(y)]={"status":"ORIGINAL_SOURCE_STAC_EXACT_ITEM",
     "source_item":originalid,"source_publisher_bbox":b,
     "source_utm_epsg":epsg,
     "original_COG_path":urllib.parse.urlparse(link).path,
     "original_cog_hostname":urllib.parse.urlparse(link).netloc,
     "source_item_geometry_type":(item.get("geometry") or {}).get("type"),
     "source_native_10m_class_tiff":True}
 return by_year,{"source_query_URL":url,"source_query_total_features":len(features),
  "source_original_MGRS_supercell_prefixes":prefixes,
  "source_first_records":original_ids[:10],"source_no_pagination":True}
def covers(b,lon,lat,padding=False):
 pad_lon=PAD_LON if padding else 0.
 pad_lat=PAD_LAT if padding else 0.
 return b[0]+pad_lon<=lon<=b[2]-pad_lon and b[1]+pad_lat<=lat<=b[3]-pad_lat
def map_sites(fr,by_year):
 count=collections.Counter();matched=[]
 by_fold=collections.Counter()
 failed=collections.Counter()
 fixed_by_fold={}
 for _,row in fr.iterrows():
  lon=float(row.lon_2018);lat=float(row.lat_2018)
  fold=int(row.pre_2018_source_only_spatial_fold)
  # Both annual products must exist in EXACT same native UTM zone as
  # the river site. STAC 01U/60U antimeridian item bbox sometimes reports
  # [-180,180] across *all* longitude! Selecting by bbox alone misassigns
  # French sites to Pacific MGRS supercells. Fix strictly by longitude-based
  # UTM zone (pre-outcome geographic source only), then original 500m bbox.
  original_UTM_zone=int(math.floor((lon+180)/6))+1
  if original_UTM_zone not in (30,31,32):
   raise RuntimeError("Actual French mainland site falls outside expected UTM zones 30,31,32")
  original_site_epsg=32600+original_UTM_zone
  pair=[]
  for stem,ys in by_year.items():
   A=ys.get("2018",{});B=ys.get("2022",{})
   if A.get("status")!="ORIGINAL_SOURCE_STAC_EXACT_ITEM" or B.get("status")!="ORIGINAL_SOURCE_STAC_EXACT_ITEM":
    continue
   if int(A["source_utm_epsg"])!=original_site_epsg:
    continue
   if A["source_utm_epsg"]!=B["source_utm_epsg"]:
    continue
   if A["source_publisher_bbox"]!=B["source_publisher_bbox"]:
    raise RuntimeError("Original yearly source MGRS grid bbox differ: "+stem)
   if covers(A["source_publisher_bbox"],lon,lat,True):
    pair.append(stem)
  if not pair:
   failed[fold]+=1
   continue
  if len(pair)>1:
   # overlapping STAC boxes require a single deterministic source assignment
   # based on original source STAC bbox, not ecological outcomes.
   pair.sort()
  selected=pair[0]
  count[selected]+=1
  by_fold[fold]+=1
  matched.append((fold,selected,str(row.site_key)))
 for fold in range(1,6):
  options=[m for m in matched if m[0]==fold]
  if options:
   options.sort(key=lambda m:hashlib.sha256(m[2].encode()).hexdigest())
   fixed_by_fold[str(fold)]={"anonymous_fixed_hash_source_site_selection":True,
      "original_source_supercell":options[0][1]}
 source_pass=len(matched)>=300 and all(by_fold[f]>=30 for f in range(1,6))
 return {"source_geo_linked_FR_sites_in_both_years_without_crop":sum(count.values()),
  "original_source_FR_sites":len(fr),"per_supercell_number_of_sites":dict(sorted(count.items())),
  "per_preexisting_2018_spatial_fold_source_valid_sites":dict(sorted(by_fold.items())),
  "source_500m_footprint_outside_same_yearpair_utm_cell_by_fold":dict(sorted(failed.items())),
  "source_hash_selected_native_cell_per_fivefold":fixed_by_fold,
  "original_sameUTM_geometry_2018_2022_site_support_GO":source_pass,
  "original_stac_bboxes_do_NOT_verify_true_internal_TIFF_raster_pixels":True,
  "source_geographic_UTM_zone_matches_original_EPSG326xx_required":True,
  "source_antimeridian_stac_bbox_false_france_overlap_rejected":True,
  "original_survey_site_ids_and_geo_coordinates_not_published":True}
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq18","eq22","out"):
  p.add_argument("--"+k,type=pathlib.Path,required=True)
 args=p.parse_args()
 fr,rzmeta=prior.sites(args)
 if len(fr)!=381 or rzmeta.get("status")!="PASS_2018_RZ2018_ALL_ORIGINAL_BENTHIC_SITES_SOURCE_ONLY_NO_ECOLOGY":
  raise RuntimeError("Original France381 source cohort not unchanged")
 original,meta=source_items()
 result=map_sites(fr,original)
 report={"schema":"GESUI_FR381_IO_LULC_V2_2018_22_EEA_RZ_ORIGINAL_SOURCE_SITE_MGRS_MAP_v1",
  "source_real_Copernicus_RZ2018_and_EEA_FRANCE381_methods_recaptured":True,
  "source_2018_RZ_EEA_2022_bio_numeric_NOT_READ":True,
  "original_EEA2018_SHA256":prior.BENTHOS2018_SHA,
  "original_EEA2022_GIT_blob":prior.BENTHOS2022_GITBLOB,
  "yearpair_product":"io-lulc-annual-v02",
  "publisher_original_YEAR":[2018,2022],
  "original_UTM_supercell_STAC_source_items":original,
  "original_spatial_search_source":meta,
  "source_UTM_antimeridian_false_allworld_bbox_item_names_excluded":["01U","60U"],
  "source_france_sites_by_original_MGRS_tiles":result,
  "source_geographical_GO":result["original_sameUTM_geometry_2018_2022_site_support_GO"],
  "source_ANNULUS_native_class_pixel_values_READ":False,
  "source_status":"PASS_ORIGINAL_2018_22_STAC_FRANCE381_SOURCE_GRID_ONLY" if result["original_sameUTM_geometry_2018_2022_site_support_GO"] else "SOURCE_STOP_original_STAC_fewer_than300_FRANCE_sites",
  "not_2018_22_numeric_benthic_eco_analysis":True}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FR381_IO_LULC2018_2022_NATIVE_UTM_TILE_SOURCE",json.dumps({
  "status":report["source_status"],"stac":meta,
  "site_support":result},sort_keys=True),flush=True)
if __name__=="__main__":main()
