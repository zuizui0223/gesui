#!/usr/bin/env python3
"""Authentic 381 France EEA 2018-2022 same-station SOURCE × original CLMS
tree cover density and confidence-layer exact 10m EPSG3035 annual tile coverage.

No 2018/2022 EQR columns loaded, no 2018/2022 raster pixels interpreted,
no station-level coordinates/source-keys copied into public aggregate output.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,pathlib,re,urllib.request,sys
import pandas as pd,numpy as np
from pyproj import Transformer

sys.path.insert(0,"analysis/20261117_RZ2018")
import eea2018_RZ2018_benthic_EQR_geoblock_model_if_fullsource as rz
BASE="https://s3.waw3-1.cloudferro.com/swift/v1/CatalogueCSV/landcover_landuse/tree_cover_density/"
TCD_ID="clms_vlcc_tree-cover-density_europe_10m_yearly_v1"
QA_ID="clms_vlcc_tree-cover-density-confidence-layer_europe_10m_yearly_v1"
TCD_SHA="5a6f0dc9f67ec7e3671f5984e9cf0425988e23dbb59ce1f6a83a18dc8661d6c9"
BENTHOS2018_SHA="2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02"
BENTHOS2022_GITBLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
GRID=re.compile(r"_E(\d+)N(\d+)_03035_",re.I)
def hashfile(path):
 return hashlib.sha256(path.read_bytes()).hexdigest()
def gitblob(path):
 data=path.read_bytes()
 return hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
def get_year_tiles(name):
 url=BASE+name+"/"+name+"_cog.csv"
 req=urllib.request.Request(url,headers={"User-Agent":"gesui official paired 2018-2022 TCD quality inventory QA/20261124"})
 with urllib.request.urlopen(req,timeout=50) as r:
  raw=r.read(12000000)
  if r.status!=200:raise RuntimeError("Original source CLMS catalogue not HTTP200")
 sha=hashlib.sha256(raw).hexdigest()
 if name==TCD_ID and sha!=TCD_SHA:raise RuntimeError("Original 6314 TCD catalogue source SHA drift")
 reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter=";")
 cols=set(reader.fieldnames or [])
 if not {"name","id","s3_path","bbox","checksum_value","content_length"}.issubset(cols):
  raise RuntimeError("CLMS catalogue publisher source schema lost")
 years={"2018":{},"2022":{}}
 others=collections.Counter()
 for row in reader:
  year="2018" if "_S2018_" in row["name"] else ("2022" if "_S2022_" in row["name"] else None)
  if year is None:continue
  match=GRID.search(row["name"])
  if not match:
   others[year]+=1
   continue
  tile=f"E{int(match.group(1))}N{int(match.group(2))}"
  if tile in years[year]:raise RuntimeError("Multiple current publisher tile IDs in same year "+tile)
  years[year][tile]={"product_uuid":row["id"],"name":row["name"],
    "source_size":int(row["content_length"]),"source_hash":row["checksum_value"],
    "original_asset_s3_path":row["s3_path"]}
 return years,{"original_csv_sha256":sha,"source_content_bytes":len(raw),
   "source_EPSG3035_year_counts":{k:len(v) for k,v in years.items()},
   "other_original_projection_counts":dict(others),
   "url":url}
def sites(args):
 if hashfile(args.eq18)!=BENTHOS2018_SHA or gitblob(args.eq22)!=BENTHOS2022_GITBLOB:
  raise RuntimeError("Original 2018 SHA/2022 git blob mismatch")
 land,meta=rz.source_frame(args.rz)
 if len(land)!=978:raise RuntimeError("Original 978 2018 vegetation-mapped station cohort drift")
 oldcols=["site_key","countryCode","procedureClassificationSystem","HYRIV_ID","lat","lon",
    "pre_2018_source_only_spatial_fold","region_3deg_block"]
 futurecols=["site_key","countryCode","procedureClassificationSystem","HYRIV_ID","lat","lon"]
 p=pd.read_csv(args.eq18,usecols=oldcols,dtype={"site_key":str})
 q=pd.read_csv(args.eq22,usecols=futurecols,dtype={"site_key":str})
 if len(p)!=2634 or len(q)!=2475 or p.site_key.duplicated().any() or q.site_key.duplicated().any():
  raise RuntimeError("Original 2018/2022 same source-site population changed")
 match=p.merge(land,on="site_key",how="inner",validate="one_to_one")
 if len(match)!=976:raise RuntimeError("Original 2018 site+RZ site 976 cohort changed")
 both=match.merge(q,on="site_key",how="inner",suffixes=("_2018","_2022"),validate="one_to_one")
 if len(both)!=418:raise RuntimeError("Original both-year RZ matching 418 site frame altered")
 ok=(both.procedureClassificationSystem_2018.fillna("").ne("")&
     both.procedureClassificationSystem_2018.eq(both.procedureClassificationSystem_2022)&
     both.HYRIV_ID_2018.eq(both.HYRIV_ID_2022)&
     both.countryCode_2018.eq(both.countryCode_2022)&
     (np.max(np.abs(both[["lat_2018","lon_2018"]].to_numpy(float)-
              both[["lat_2022","lon_2022"]].to_numpy(float)),axis=1)<=.01))
 selected=both.loc[ok]
 if len(selected)!=409:raise RuntimeError("Original site stability/source-only count 409 changed")
 france=selected.loc[selected.countryCode_2018.eq("FR")].copy()
 if len(france)!=381:raise RuntimeError("Original source France 381 same-site cohort changed")
 folds=france.pre_2018_source_only_spatial_fold.value_counts().to_dict()
 if {int(k):int(v) for k,v in folds.items()}!={1:122,2:58,3:48,4:33,5:120}:
  raise RuntimeError("Original source-only France five folds drift")
 return france,meta
def key_from_coord(x,y):
 return f"E{int(math.floor(x/100000))}N{int(math.floor(y/100000))}"
def ring_needed_tiles(x,y,radius):
 keys=set()
 for dx in (-radius,radius):
  for dy in (-radius,radius):
   keys.add(key_from_coord(x+dx,y+dy))
 return keys
def summarize(france,product,quality):
 transformer=Transformer.from_crs(4326,3035,always_xy=True)
 supported=collections.Counter()
 failed={x:[] for x in ("tcd","tcd_quality")}
 allkeys=set();tiles_per_site=collections.Counter()
 example_by_fold={}
 source_genomes=[]
 for _,s in france.iterrows():
  lon=float(s.lon_2018);lat=float(s.lat_2018)
  X,Y=transformer.transform(lon,lat)
  if not np.isfinite([X,Y]).all():raise RuntimeError("Actual original French source coord projected invalid")
  needed=ring_needed_tiles(X,Y,500.)
  tiles_per_site[len(needed)]+=1
  allkeys.update(needed)
  for field,cat in (("tcd",product),("tcd_quality",quality)):
   present=all(t in cat["2018"] and t in cat["2022"] for t in needed)
   supported[field]+=int(present)
   if not present:
    failed[field].append({"fold":int(s.pre_2018_source_only_spatial_fold),
      "missing_tiles":sorted(t for t in needed if t not in cat["2018"] or t not in cat["2022"])})
  fold=int(s.pre_2018_source_only_spatial_fold)
  source_genomes.append((fold,s.site_key,sorted(needed)))
 # The examples are SHA-selected within PRE-FROZEN 5 geography folds, not
 # selected on taxon/EQR/canopy values or on a valid satellite tile status.
 for fold in (1,2,3,4,5):
  cohort=[q for q in source_genomes if q[0]==fold]
  if not cohort:raise RuntimeError("Original geographical fold has no France station")
  cohort.sort(key=lambda v:hashlib.sha256(v[1].encode()).hexdigest())
  row=cohort[0]
  example_by_fold[str(fold)]={"anonymous_original_stable_source_selection":"lowest_SHA256_site_key",
    "source_France_TCD_100km_tile_keys_required_for_two_rings":row[2]}
 allpairs={
  t:{"tcd_2018_product":product["2018"].get(t,{}).get("name"),
     "tcd_2022_product":product["2022"].get(t,{}).get("name"),
     "quality_2018_product":quality["2018"].get(t,{}).get("name"),
     "quality_2022_product":quality["2022"].get(t,{}).get("name")}
  for t in sorted(allkeys)}
 return {
  "source_real_FR_same_monitoring_sites":len(france),
  "source_riparian_annulus_radii_0to250_and_250to500m":True,
  "source_nearest_raster_tile_bbox_assumes_EPSG3035_100km_grid":True,
  "source_true_site_geometry_only_not_pixel_value":True,
  "native_tree_cover_dens_2018_and_2022_tile_coverage_any_500m_annulus":supported["tcd"],
  "native_TCD_confidence_layer_2018_and_2022_tile_coverage_any_500m_annulus":supported["tcd_quality"],
  "source_tcd_and_quality_100km_tile_union_count":len(allkeys),
  "original_source_tile_keys_needed":sorted(allkeys),
  "number_100km_tiles_per_500m_site":dict(sorted(tiles_per_site.items())),
  "source_missing_coverage":{k:{"site_count":len(v),
    "missing_tiles_example_sanitized":v[:6]} for k,v in failed.items()},
  "source_exact_tile_2018_2022_original_item_names":allpairs,
  "fold_deterministic_source_site_example_tile_keys":example_by_fold,
  "no_true_site_coordinate_or_site_ID_rehosted":True,
  "actual_2018_or_2022_EEA_EQR_values_read":False}
def main():
 p=argparse.ArgumentParser()
 for k in ("rz","eq18","eq22","out"):p.add_argument("--"+k,type=pathlib.Path,required=True)
 a=p.parse_args()
 fr,src=sites(a)
 tcd,tcdqa=get_year_tiles(TCD_ID)
 conf,confqa=get_year_tiles(QA_ID)
 matched=summarize(fr,tcd,conf)
 verdict="PASS_GEO_AND_2YEAR_PUBLIC_CATALOGUES_BUT_RASTER_PIXEL_ACCESS_STOP"
 obj={"schema":"gesui_FR381_2018_22_true_revisited_sites_CLMS_TCD_TCDCL_exact_catalogue_source_map_v1",
   "original_EEA2018_EQR_source_SHA256":BENTHOS2018_SHA,
   "original_EEA2022_EQR_source_git_blob_SHA":BENTHOS2022_GITBLOB,
   "original_RZ2018_source_status":src["status"],
   "TCD_publisher_source":tcdqa,"TCD_confidence_layer_publisher_source":confqa,
   "source_TCD_both_year_original_projected_tile_coverage":matched,
   "original_TCD_annual_pixel_values_accessed":False,
   "original_TCDCL_quality_pixel_values_accessed":False,
   "source_real_year_delta_TCD_or_delta_EQR_calculated":False,
   "original_raster_access_prerequisite":"CDSE OIDC/S3 credentials; anonymous GET of original STAC TIFF returned 403 in Actions 45206093035",
   "source_status":verdict,
   "source_quality_if_pixel_access_not_possible":"SOURCE_STOP_NO_2018_2022_CANOPY_DELTA"};
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(obj,sort_keys=True,indent=2)+"\n")
 print("FR381_TCD2018_2022_TCDCL_SOURCE_SITE_TILE_COUNTS",json.dumps({
  "status":verdict,
  "france_stable_sites":len(fr),
  "TCD_source":tcdqa,"TCDCL_source":confqa,
  "coverage":{k:v for k,v in matched.items() if k not in (
    "source_exact_tile_2018_2022_original_item_names","source_missing_coverage")},
  "failed":matched["source_missing_coverage"]},sort_keys=True),flush=True)
if __name__=="__main__":main()
