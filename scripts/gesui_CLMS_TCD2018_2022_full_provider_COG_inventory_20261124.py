#!/usr/bin/env python3
"""Exact publisher CLMS TCD annual 10m original COG catalogue full inventory.

Complete source 3,059,396-byte official public CSV, semicolon separated.
Collect 2018/2022 tile commonality, literal original BBOX schema, and
source S3 path format without any EQR values or raster pixel measurement.
"""
import argparse,collections,csv,datetime,hashlib,io,json,pathlib,re,urllib.request
URL="https://s3.waw3-1.cloudferro.com/swift/v1/CatalogueCSV/landcover_landuse/tree_cover_density/clms_vlcc_tree-cover-density_europe_10m_yearly_v1/clms_vlcc_tree-cover-density_europe_10m_yearly_v1_cog.csv"
ORIGINAL_SOURCE_SIZE=3059396
ORIGINAL_SOURCE_ETAG="aee6a89c269acc93046033666cf7296e"
GRID=re.compile(r"_E(\d+)N(\d+)_03035_",re.I)
def fetch():
 req=urllib.request.Request(URL,headers={"User-Agent":"gesui original CLMS TCD catalogue audit/20261124"})
 with urllib.request.urlopen(req,timeout=50) as r:
  raw=r.read(12000000)
  head=dict(r.headers)
  if r.status!=200 or len(raw)!=ORIGINAL_SOURCE_SIZE or r.headers.get("ETag","").strip('"')!=ORIGINAL_SOURCE_ETAG:
   raise RuntimeError("CLMS public CSV native source bytes / ETag changed")
  return raw,head
def run(args):
 raw,headers=fetch()
 reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter=";")
 expected=["id","name","content_length","ingestion_date","content_date_start","content_date_end",
           "nominal_date","modification_date","checksum_algorithm","checksum_value","s3_path","bbox"]
 if reader.fieldnames!=expected:raise RuntimeError("Original CLMS catalogue schema changed "+str(reader.fieldnames))
 years=collections.Counter()
 records={str(y):{} for y in (2018,2022)}
 per_year_dup=collections.Counter();bad_bbox=0;extensions=collections.Counter()
 examples={str(y):[] for y in (2018,2022)}
 for row in reader:
  if not row["name"].startswith("CLMS_HRLVLCC_TCD_S"):raise RuntimeError("Unexpected TCD COG source label")
  s=re.search(r"_S(20\d\d)_",row["name"])
  if not s:raise RuntimeError("No original annual source key")
  year=s.group(1);years[year]+=1
  extensions[row["s3_path"].rsplit(".",1)[-1][:12]]+=1
  if not row["bbox"]:bad_bbox+=1
  if year in records:
   m=GRID.search(row["name"])
   if not m:raise RuntimeError("No source/native 3035 COG tile key in annual file")
   key="E"+m.group(1)+"N"+m.group(2)
   if key in records[year]:per_year_dup[year]+=1
   records[year][key]={"original_name":row["name"],
    "original_s3_path":row["s3_path"],
    "original_product_id":row["id"],"original_hash_algorithm":row["checksum_algorithm"],
    "original_content_hash":row["checksum_value"],"original_content_length":row["content_length"],
    "bbox":row["bbox"],"nominal_date":row["nominal_date"],
    "ingestion_date":row["ingestion_date"]}
   if len(examples[year])<4:examples[year].append({
    "original_product":row["name"],"original_size":row["content_length"],
    "bbox_format_example":row["bbox"][:300],
    "s3_path_format_example":row["s3_path"][:300]})
 common=sorted(set(records["2018"])&set(records["2022"]))
 if len(common)<100:raise RuntimeError("CLMS 2018 and 2022 native same-tile overlap <100")
 if any(per_year_dup.values()):raise RuntimeError("CLMS COG duplicate EPSG3035 annual tile ambiguity")
 output={"schema":"gesui_2026_CLMS_original_TCD_10m_2018_2022_same_grid_full_catalogue_source_v1",
  "source_dataset_identifier":"clms_vlcc_tree-cover-density_europe_10m_yearly_v1",
  "source_original_CSV_URL":URL,
  "source_original_CSV_ETAG":ORIGINAL_SOURCE_ETAG,
  "source_original_CSV_bytes":len(raw),
  "source_original_CSV_SHA256":hashlib.sha256(raw).hexdigest(),
  "source_response_row_year_counts":dict(sorted(years.items())),
  "source_all_year_COG_file_count":sum(years.values()),
  "source_EPSG3035_tile_count_2018":len(records["2018"]),
  "source_EPSG3035_tile_count_2022":len(records["2022"]),
  "source_same_TCD_EPSG3035_tile_ID_count_2018_and_2022":len(common),
  "source_same_2018_2022_tile_keys_sorted":common,
  "source_2018_2022_short_native_records":{
    "2018":[records["2018"][x] for x in common[:2]],
    "2022":[records["2022"][x] for x in common[:2]]},
  "original_per_year_preview":examples,
  "source_s3_path_extension_counts":dict(extensions),
  "missing_source_bbox":bad_bbox,
  "source_publisher_actual_2018_and_2022_item_presence":True,
  "original_anonymous_COG_pixel_access_yet_verified":False,
  "source_region_tiled_raster_coverage_at_actual_FR381_sites_verified":False,
  "no_EEA_2018_or_2022_biological_outcomes_opened":True,
  "source_years_requested":[2018,2022],
  "source_status":"PASS_FULL_OFFICIAL_2YEAR_TCD_TILE_INVENTORY_ONLY_STOP_NO_PIXELS",
  "next":"Link exact 381 source-only France sites to 10m COG tiles by authoritative EPSG3035 bbox, then test HTTPS byte-range genuine TIFF/COG header for both years and provider quality layers."}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
 print("CLMS_TCD2018_2022_COMPLETE_ORIGINAL_INVENTORY",json.dumps({
  "source_sha256":output["source_original_CSV_SHA256"],
  "annual_COG_count":output["source_response_row_year_counts"],
  "samegrid":len(common),"years":examples,
  "no_EQR_or_original_pixel_values_read":True},sort_keys=True),flush=True)
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--out",type=pathlib.Path,required=True)
 run(p.parse_args())
