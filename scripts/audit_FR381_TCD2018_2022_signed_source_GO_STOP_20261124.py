#!/usr/bin/env python3
"""Signed ORIGINAL CLMS TCD 2018/2022 access–coverage–quality publication GO/STOP.

Reconciles four separately executed source-only original publisher evidence
artifacts. No biological EQR coefficients or canopy raster values are read.
Any true pixel access must be independently demonstrated; 403 is a source STOP,
not a null result for riparian vegetation or for aquatic organisms.
"""
import argparse,hashlib,json,pathlib
TCD_SHA="5a6f0dc9f67ec7e3671f5984e9cf0425988e23dbb59ce1f6a83a18dc8661d6c9"
TCDCL_SHA="347053a15fe8e3cf5df7fa261daf7ae0cf9dcab4cf1fe5f81237a28f735b4222"
YEARS=("2018","2022")
def sha(p):
 return hashlib.sha256(p.read_bytes()).hexdigest()
def need(ok,why):
 if not ok:raise RuntimeError("TCD SOURCE CONTRACT DRIFT / STOP: "+why)
def main():
 p=argparse.ArgumentParser()
 for key in ("initial","inventory","raster_access","france_sites","out"):
  p.add_argument("--"+key,type=pathlib.Path,required=True)
 a=p.parse_args()
 initial=json.loads(a.initial.read_text())
 inv=json.loads(a.inventory.read_text())
 access=json.loads(a.raster_access.read_text())
 fr=json.loads(a.france_sites.read_text())
 need(initial["schema"]=="gesui_Copernicus_CLMS_TCD2018_TCD2022_France381_catalogue_source_probe_v1","Initial official catalogue source")
 need(inv["schema"]=="gesui_2026_CLMS_original_TCD_10m_2018_2022_same_grid_full_catalogue_source_v1","Official full catalogue")
 need(access["schema"]=="gesui_official_TCD2018_2022_France_real_COG_anonymous_pixel_access_pilot_v1","Actual raster HTTP access")
 need(fr["schema"]=="gesui_FR381_2018_22_true_revisited_sites_CLMS_TCD_TCDCL_exact_catalogue_source_map_v1","Original source matched France381")
 need(inv["source_original_CSV_SHA256"]==TCD_SHA==fr["TCD_publisher_source"]["original_csv_sha256"],"Original TCD CSV source SHA")
 need(fr["TCD_confidence_layer_publisher_source"]["original_csv_sha256"]==TCDCL_SHA,"Original publisher TCDCL QA source SHA")
 need(inv["source_response_row_year_counts"]["2018"]==902 and inv["source_response_row_year_counts"]["2022"]==902,"902 source annual tiles each")
 need(inv["source_same_TCD_EPSG3035_tile_ID_count_2018_and_2022"]==897,"897 paired source EPSG3035 tiles")
 for year in YEARS:
  need(fr["TCD_publisher_source"]["source_EPSG3035_year_counts"][year]==897,"TCD original 897 year "+year)
  need(fr["TCD_confidence_layer_publisher_source"]["source_EPSG3035_year_counts"][year]==897,"TCD confidence quality 897 year "+year)
 need(fr["source_original_2022_EQR_git_blob_SHA"]=="325f199e985d9862ffee5b0edb205c61a7be3833" if "source_original_2022_EQR_git_blob_SHA" in fr else
       fr["original_EEA2022_EQR_source_git_blob_SHA"]=="325f199e985d9862ffee5b0edb205c61a7be3833","Source original 2022 EEA Git blob")
 src=fr["source_TCD_both_year_original_projected_tile_coverage"]
 need(src["source_real_FR_same_monitoring_sites"]==381,"Original same-site FRevents")
 need(src["native_tree_cover_dens_2018_and_2022_tile_coverage_any_500m_annulus"]==381,"Original 381 TCD two-year footprints")
 need(src["native_TCD_confidence_layer_2018_and_2022_tile_coverage_any_500m_annulus"]==381,"Original 381 TCDCL confidence footprints")
 need(src["source_tcd_and_quality_100km_tile_union_count"]==47,"Correct source 47 official EPSG3035 tile labels")
 need(src["number_100km_tiles_per_500m_site"].get("1")==374 and src["number_100km_tiles_per_500m_site"].get("2")==7,"Source original 500m annulus tile boundary counts")
 need(len(src["original_source_tile_keys_needed"])==47 and len(set(src["original_source_tile_keys_needed"]))==47,"Source unique tile keys")
 need(src["actual_2018_or_2022_EEA_EQR_values_read"] is False,"Source-only EQR outcomes must not have been opened")
 need(access["TCD_original_pixel_values_READ"] is False,"Original numeric canopy data should not have been claimed before pixel GO")
 need(access["both_years_real_TIFF_source_pass"] is False,"Source anonymous pixel access unexpectedly PASS; must review!")
 need(access["source_status"]=="STOP_NO_PUBLIC_SOURCE_TIFF_FOR_BOTH_YEARS","Both native COG source headers gate")
 need(access["source_original_STAC_exact_collection_probes"]["clms_vlcc_tree-cover-density_europe_10m_yearly_v1"]["collection"]["code"]==200,"Official STAC collection now confirmed")
 checks=access["source_official_COG_data_byte_range_checks"]
 need(len(checks)==2,"Expected both 2018 & 2022 exact original COG headers")
 need(set(x["year"] for x in checks)==set(YEARS),"Year-specific original COG source checks")
 for item in checks:
  need(item["raw_byte_range"]["code"]==403 and item["first_bytes_are_TIFF"] is False,"Source TIFF anonymously inaccessible 403 "+item["year"])
  need(item["source_asset"]["auth"]==["s3"],"STAC original COG asset authentication flag")
  need(item["source_asset"]["stac_media_type"].startswith("image/tiff"),"Official TIFF data type")
 for year in YEARS:
  need(access["source_original_items_2018_2022"][year]["original_name"].startswith("CLMS_HRLVLCC_TCD_S"+year+"_"),"Actual original year name")
  assets=access["source_original_items_2018_2022"][year]["source_assets"]
  need(any(x["asset_key"]=="product" and x["auth"]==["oidc"] for x in assets),"Original OData archive OIDC flag year "+year)
 report={
  "schema":"GESUI_FR381_TCD2018_2022_ACCESS_AND_SAMPLING_PUBLISHABLE_SOURCE_GO_STOP_v1",
  "original_public_French_tcd_2018_2022_geographic_support":{
    "original_EEA_FR_2018_2022_same_source_sites":381,
    "official_EPSG3035_100km_tile_union":47,
    "2018_TCD_available_matching_product_tiles":381,
    "2022_TCD_available_matching_product_tiles":381,
    "2018_and_2022_annual_confidence_qa_tile_catalogue_coverage":381,
    "real_2018_2022_site_canopy_pixels_ever_downloaded":False,
    "real_TCD_annual_land_cover_change_computed":False
   },
  "native_original_COG_2018_and_2022_anonymous_HTTP":{"2018":403,"2022":403},
  "original_Copernicus_2018_2022_dataset_DOI":"10.2909/e677441e-fb94-431c-b4f9-304f10e4dfd8",
  "fixed_source_definition":"TCD_10m_original_true_riparian_source_0to250_and_250to500m_years_2018_and_2022",
  "result_verdict":{
    "genuine_provider_catalogue_CSV_and_OData":"PASS",
    "exact_original_2018_and_2022_COG_items":"PASS",
    "actual_FR381_both_annual_quality_and_TCD_tile_footprints":"PASS",
    "both_years_anonymous_TIFF_data_access":"STOP_HTTP_403_AUTH_REQUIRED",
    "new_381site_delta_canopy_magnitude":"SOURCE_STOP_NOT_MEASURED",
    "new_land_change_plus_delta_EQR_statistical_model":"SOURCE_STOP_NOT_FITTED",
    "2018_RZ_LCLU_to_2022_yearly_LCLU_class_equivalence":"SOURCE_STOP_NOT_AN_OFFICIAL_2022_RZ_CLASSIFICATION",
    "sewage_specific_biology_or_trophic_subsidy":"SOURCE_STOP_UNMEASURED"
  },
  "source_inventory_signatures":{
    "provider_original_TCD_CSV_SHA256":TCD_SHA,
    "provider_original_TCDCL_CSV_SHA256":TCDCL_SHA,
    "initial_original_HTTP_probe_JSON_SHA256":sha(a.initial),
    "full_2018_and_2022_original_inventory_JSON_SHA256":sha(a.inventory),
    "direct_2018_2022_original_COG_access_JSON_SHA256":sha(a.raster_access),
    "matched_381_sites_and_quality_tile_support_JSON_SHA256":sha(a.france_sites)},
  "source_outcomes_never_opened":True,
  "source_stop_is_not_null_ecological_result":True,
  "next_source_access_requirement":"Authorized Copernicus CDSE S3/OIDC account, native TCD and confidence raster for both 2018 and 2022. All 381 existing exact site IDs and 47 native tiles source-mapped; source-only eligible != quality-pixel GO.",
  "do_not_restore_previous_scientific_positive_by_change_proxy":True
 }
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print("FR381_2018_2022_TCD_ORIGINAL_SOURCE_AUDIT_FINAL",json.dumps(report,sort_keys=True),flush=True)
if __name__=="__main__":main()
