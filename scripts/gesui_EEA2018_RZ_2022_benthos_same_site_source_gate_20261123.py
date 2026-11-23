#!/usr/bin/env python3
"""Original EEA 2018 Riparian Zones station cohort × 2022 same monitoring site SOURCE.

This script deliberately reads ONLY source site/method/river/latlon columns:
DOES NOT READ original 2018 numeric_EQR_2018 or 2022 normalised_EQR.
Reuses frozen Copernicus RZ2018 mapped class-area source, never interpolates.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,pathlib,sys
import numpy as np,pandas as pd
sys.path.insert(0,"analysis/20261117_RZ2018")
import eea2018_RZ2018_benthic_EQR_geoblock_model_if_fullsource as rz

PR14_SOURCE_BLOB="325f199e985d9862ffee5b0edb205c61a7be3833"
PR17_ORIGINAL_2018_SHA="2a1a72d43656aca588a832dfae13860054c01a651e76bf25f984c5a374709b02"
EXPECTED_2018=2634
EXPECTED_2022=2475
EXPECTED_CROSSYEAR=1167
EXPECTED_STABLE_METHOD=1148
SITES_2018=["site_key","countryCode","observedPropertyDeterminandBiologyEQRCode",
 "procedureClassificationSystem","lat","lon","region_3deg_block",
 "pre_2018_source_only_spatial_fold","HYRIV_ID",
 "DIS_AV_CMS","UPLAND_SKM","2018_measured_upstream_P_t_per_year"]
SITES_2022=["site_key","countryCode","observedPropertyDeterminandBiologyEQRCode",
 "procedureClassificationSystem","lat","lon","region_3deg_block",
 "HYRIV_ID","DIS_AV_CMS","UPLAND_SKM","upstream_2022_EEA_P_measured_t_per_yr"]
def sha(path):
 return hashlib.sha256(path.read_bytes()).hexdigest()
def bio_source(path,cols,expected):
 d=pd.read_csv(path,usecols=cols,dtype={"site_key":str})
 if len(d)!=expected or d.site_key.duplicated().any():
  raise RuntimeError("Original source 2018/22 biological unique site population changed")
 if set(d.observedPropertyDeterminandBiologyEQRCode)!={"EEA_13-01-4"}:
  raise RuntimeError("Original benthic BQE source element method code drift")
 return d
def summarize(x):
 return dict(sorted(collections.Counter(x).items()))
def run(args):
 if sha(args.eq18)!=PR17_ORIGINAL_2018_SHA:
  raise RuntimeError("Authentic numeric-2018 publisher SOURCE file original SHA changed")
 # The GH blob SHA of 2022 source must be independently checked in workflow
 # using git hash-object and logged in the output for code-level identity.
 old=bio_source(args.eq18,SITES_2018,EXPECTED_2018)
 future=bio_source(args.eq22,SITES_2022,EXPECTED_2022)
 common=old.merge(future,on="site_key",how="inner",suffixes=("_2018","_2022"),validate="one_to_one")
 if len(common)!=EXPECTED_CROSSYEAR:
  raise RuntimeError("Original EEA2018→22 all-site 1167 source overlap changed")
 same_method=common.procedureClassificationSystem_2018.eq(common.procedureClassificationSystem_2022)
 same_reach=common.HYRIV_ID_2018.eq(common.HYRIV_ID_2022)
 if (int(same_method.sum()),int(same_reach.sum()))!=(EXPECTED_STABLE_METHOD,EXPECTED_CROSSYEAR):
  raise RuntimeError("Original 2018/22 method/reach same-site SOURCE recapture failed")
 if not common.countryCode_2018.eq(common.countryCode_2022).all():
  raise RuntimeError("Country assignment drift for original same site")
 oldcoords=common[["lat_2018","lon_2018"]].to_numpy(float)
 newcoords=common[["lat_2022","lon_2022"]].to_numpy(float)
 diff=np.abs(oldcoords-newcoords)
 small=np.all(diff<=.01,axis=1)
 # Potential source coordinates updated between reporting years: report them
 # and never silently recode site identity from proximity.
 land,meta=rz.source_frame(args.rz)
 if len(land)!=978 or land.site_key.duplicated().any():
  raise RuntimeError("Original real riparian 978 terrestrial eligible stations changed")
 if len(set(land.site_key)&set(old.site_key))!=978:
  raise RuntimeError("Original terrestrial 2018 points absent in 2018 source cohort")
 matched=common.merge(land,on="site_key",how="inner",validate="one_to_one")
 consistent=matched.procedureClassificationSystem_2018.eq(matched.procedureClassificationSystem_2022)&(
   matched.HYRIV_ID_2018.eq(matched.HYRIV_ID_2022))&(
   matched.countryCode_2018.eq(matched.countryCode_2022))&(
   np.max(np.abs(matched[["lat_2018","lon_2018"]].to_numpy(float)-
                     matched[["lat_2022","lon_2022"]].to_numpy(float)),axis=1)<=.01)
 method_n=int(matched.procedureClassificationSystem_2018.eq(matched.procedureClassificationSystem_2022).sum())
 final=matched.loc[consistent].copy()
 exact=len(final)
 year_regions=final.region_3deg_block_2018
 country_counts=summarize(final.countryCode_2018)
 countries_ge10=sum(v>=10 for v in country_counts.values())
 geos=year_regions.nunique()
 fold_counts=summarize(final.pre_2018_source_only_spatial_fold.astype(int))
 quality_for_all=method_n/max(1,len(matched))
 multi=bool(exact>=150 and countries_ge10>=5 and geos>=10 and quality_for_all>=.9)
 france=final.loc[final.countryCode_2018.eq("FR")].copy()
 fr_regions=france.region_3deg_block_2018.value_counts()
 frfold=france.pre_2018_source_only_spatial_fold.value_counts()
 france_go=bool(len(france)>=150 and len(fr_regions)>=5 and
                (frfold>=20).sum()>=3 and quality_for_all>=.9)
 # Original row EQR not opened; source-only landscape primary predictor
 # source ranges and genuine positive/report-absence source only.
 source={
  "schema":"GESUI_EEA2018_RZ2018_TO_EEA2022_BENTHIC_COHORT_SOURCE_ONLY_V1",
  "source_stage_DOES_NOT_READ_2018_or_2022_EQR_numeric":True,
  "2018_source_n":len(old),"2022_source_n":len(future),
  "historical_all_source_exact_revisit_n":len(common),
  "original_same_method_all_source_n":int(same_method.sum()),
  "original_same_HYRIV_all_source_n":int(same_reach.sum()),
  "original_same_site_coords_within01deg_n":int(small.sum()),
  "2018_RZ_true_land_eligible_station_n":len(land),
  "same2018RZ_and_2022_revisit_all_n":len(matched),
  "same2018RZ_2022_method_stable_n":method_n,
  "same2018RZ_2022_method_reach_coords_stable_n":exact,
  "same2018RZ_2022_stable_method_fraction_of_RZ_revisited":float(quality_for_all),
  "country_counts_stable_source":country_counts,
  "independent_2018_3degree_region_source_n":int(geos),
  "original_fixed_2018_geoblock_fold_counts_stable_source":fold_counts,
  "source_three_deg_blocks_france":int(len(fr_regions)),
  "source_france_stable_site_n":len(france),
  "source_france_original_geo_fold_counts":{str(k):int(v) for k,v in frfold.sort_index().items()},
  "multi_country_future_outcome_GO":multi,
  "frozen_france_only_fallback_future_outcome_GO":france_go,
  "source_frozen_acceptance":"MULTICOUNTRY" if multi else ("FRANCE_ONLY" if france_go else "STOP"),
  "2018_reported_source_P_positive_in_stable_RZ_revisit_n":int((final["2018_measured_upstream_P_t_per_year"]>0).sum()),
  "2022_reported_source_P_positive_in_stable_RZ_revisit_n":int((final["upstream_2022_EEA_P_measured_t_per_yr"]>0).sum()),
  "both_years_reported_source_P_positive_n":int(((final["2018_measured_upstream_P_t_per_year"]>0)&(final["upstream_2022_EEA_P_measured_t_per_yr"]>0)).sum()),
  "no_2018_to_2022_true_P_effluent_intervention_effect_identifiable":True,
  "source_original_2018_EQR_file_SHA256":PR17_ORIGINAL_2018_SHA,
  "source_original_2022_EQR_git_blob_SHA":PR14_SOURCE_BLOB,
  "original_RZ2018_source_status":meta["status"],
  "2018_landclass_coverage_is_not_nutrient_or_growth":True,
  "2022_RZ_landcover_update_not_measured":True,
  "2022_independent_raw_biological_source_not_opened":True}
 # We intentionally DO NOT export real EEA site IDs/source lat lon in source JSON.
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(source,indent=2,sort_keys=True)+"\n")
 print("EEA_RZ2018_TO_BENTHIC2022_SOURCE_GATE",json.dumps(source,sort_keys=True),flush=True)
 return source
def main():
 p=argparse.ArgumentParser()
 for key in ("rz","eq18","eq22","out"):
  p.add_argument("--"+key,type=pathlib.Path,required=True)
 run(p.parse_args())
if __name__=="__main__":main()
