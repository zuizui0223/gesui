#!/usr/bin/env python3
"""Source-only correction: EPA NPDES 'Permitted Feature' != 'External Outfall'.

Using EXACT original 2026-11-08 official national ECHO outfall ZIP SHA, match
SUB_TYPE_DESC to existing frozen 568 NWCA sites and original same-2021-date
municipal source records. Zero H_ALL, TN, TP or FQAI outcomes read.
"""
from __future__ import annotations
import argparse,hashlib,json,math,zipfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd

from gesui_NWCA2021_568source_EPA_NPDES_registered_sewage_outfall_geoONLY import (
    original_nwca_source_sites,ndpes_outfall_source,assess_sites,sha as source_sha)

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
PRIOR=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
ZIP_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
MIN_MATCH_NEAR=60
MIN_MATCH_FAR=60
MIN_STATES=10
MIN_STATE_ECO=15

def source_types(path,clean):
    with zipfile.ZipFile(path) as z:
        member=[name for name in z.namelist()
                if name.split("/")[-1].lower()=="npdes_outfalls_layer.csv"]
        if len(member)!=1:raise ValueError("Expected exact one official ECHO NPDES CSV")
        fields=["EXTERNAL_PERMIT_NMBR","PERM_FEATURE_NMBR","LATITUDE83","LONGITUDE83",
                "LATLONG_TYPE","SUB_TYPE_DESC"]
        with z.open(member[0]) as f:
            raw=pd.read_csv(f,usecols=fields,low_memory=False,dtype=str,keep_default_na=False)
    if len(raw)!=866359:raise ValueError("Original 2026-11-08 ECHO outfall row count changed")
    x=raw
    x["PERMIT_ID"]=x.EXTERNAL_PERMIT_NMBR.str.strip().str.upper()
    x["FEATURE_ID"]=x.PERM_FEATURE_NMBR.str.strip().str.upper()
    # Optimize source-type join by only looking up original 58,231 sewage permit
    # point IDs; neither species responses nor source buffers affect source type.
    x=x[x.PERMIT_ID.isin(set(clean.PERMIT_ID))].copy()
    x["LAT"]=pd.to_numeric(x.LATITUDE83,errors="coerce")
    x["LON"]=pd.to_numeric(x.LONGITUDE83,errors="coerce")
    x=x.loc[x.LAT.notna()&x.LON.notna()].copy()
    x["SUB_TYPE_DESC"]=x.SUB_TYPE_DESC.astype("string").str.strip()
    # Exact original official permit + feature + location + coordinate category.
    keys=["PERMIT_ID","FEATURE_ID","LAT","LON","LATLONG_TYPE"]
    type_conflict=x.groupby(keys,dropna=False).SUB_TYPE_DESC.nunique(dropna=False)
    nconf=int((type_conflict>1).sum())
    if nconf:
        raise ValueError(f"{nconf} SAME official permit-feature-coordinate has contradictory source types")
    x=x.drop_duplicates(keys).copy()
    if x.duplicated(keys).any() or clean.duplicated(keys).any():
        raise ValueError("Original source type join cardinality ambiguity")
    both=clean.merge(x[keys+["SUB_TYPE_DESC"]],on=keys,how="left",
                     indicator=True,validate="one_to_one")
    if len(both)!=len(clean) or (both["_merge"]!="both").any():
        mismatched=both.loc[both["_merge"]!="both"]
        raise ValueError(f"Original sewage-filtered ECHO points missing type: {len(mismatched)}")
    both=both.drop(columns="_merge")
    # No recoding unknown types as External Outfall.
    both["STRICT_EXTERNAL_OUTFALL"]=both.SUB_TYPE_DESC.fillna("").str.casefold().eq("external outfall")
    return both,{
        "ECHO_2026_original_complete_national_source_rows":len(raw),
        "source_rows_matching_original_sewage_permit_IDs_before_exact_site_join":len(x),
        "same_original_sewage_feature_geo_rows_exact_source_type_joined":len(both),
        "conflicting_original_ECHO_feature_type_keys":nconf,
        "all_original_sewage_permitted_feature_subtype_counts":dict(Counter(
            both.SUB_TYPE_DESC.fillna("MISSING").astype(str))),
        "external_outfall_exact_permit_feature_source_rows":int(both.STRICT_EXTERNAL_OUTFALL.sum()),
        "external_outfall_distinct_NPDES_administrative_permits":int(both.loc[
            both.STRICT_EXTERNAL_OUTFALL,"PERMIT_ID"].nunique()),
        "unknown_not_external_outfall":int((~both.STRICT_EXTERNAL_OUTFALL).sum()),
        "source_attribute_not_2021_actual_N_or_P_discharge":True,
    }

def geometry_only(root,nwca_dir,npdes_zip):
    original=json.loads((root/PRIOR).read_text(encoding="utf-8"))
    if original["schema"]!="gesui_US_EPA_NWCA2021_568source_NPDES2026_reported_sewage_outfall_GEO_2021permitDates_sourceONLY_v01":
        raise ValueError("Original NPDES 568-site result not frozen")
    if source_sha(npdes_zip)!=ZIP_SHA or original["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]!=ZIP_SHA:
        raise ValueError("Original ECHO ZIP must match pinned Nov 8 2026 SHA")
    if not original["no_original_NWCA_H_ALL_FQAI_or_TNTP_outcomes_accessed_in_source_stage"]:
        raise ValueError("Original pre-biological source data contaminated")
    nwca,siteaudit=original_nwca_source_sites(root,nwca_dir)
    clean,source_stats=ndpes_outfall_source(npdes_zip)
    base_src,base_stats=assess_sites(nwca,clean,"strict")
    expected=original["predeclared_nearest_sewage_outfall_geodesic_source_controls"]["strict"]
    for key in ("site_n","site_near5km_n","site_near10km_n","site_far_geodesic20km_or_more_n",
                "nearest_10km_administrative_unique_NPDES_permits_n",
                "states_with_near_and_far_from_same_state_and_ecoregion",
                "matched_far_sites_from_same_state_and_ecoregion_as_a_near_site_n"):
        if base_stats[key]!=expected[key]:
            raise ValueError(f"Original strict NPDES source not reproducible from exact original ZIP: {key}")
    if base_stats["site_near10km_n"]!=128 or base_stats["site_far_geodesic20km_or_more_n"]!=311:
        raise ValueError("Original 128 near/311 far source cohort changed")
    classified,type_stats=source_types(npdes_zip,clean)
    external=classified.loc[classified.STRICT_EXTERNAL_OUTFALL].copy()
    if len(external)<1000:raise ValueError("National sewage external outfall reporting implausibly sparse")
    new_sites,stats=assess_sites(nwca,external,"strict")
    if len(new_sites)!=568 or len(base_src)!=568:
        raise ValueError("Original fixed NWCA site population changed")
    oldids={x["NWCA2021_source_site_visit_SHA256"]:x for x in base_src}
    newids={x["NWCA2021_source_site_visit_SHA256"]:x for x in new_sites}
    if len(oldids)!=len(newids)!=568 or set(oldids)!=set(newids):
        raise ValueError("Original source NWCA site hashes mismatched after type restriction")
    old_near={i for i,x in oldids.items() if x["outfall_rows_within_10km"]>0}
    new_near={i for i,x in newids.items() if x["outfall_rows_within_10km"]>0}
    if not new_near.issubset(old_near):
        raise ValueError("Restricting prior outfalls created new <=10km sites unexpectedly")
    old_far={i for i,x in oldids.items()
             if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
                x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    new_far={i for i,x in newids.items()
             if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
                x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    if not old_far.issubset(new_far):
        raise ValueError("Restricting original source features shrank far20 control pool")
    nearstateeco={(x["original_state"],x["original_ecoregion"]) for x in new_sites if x["outfall_rows_within_10km"]>0}
    farstateeco={(x["original_state"],x["original_ecoregion"]) for x in new_sites if x["NWCA2021_source_site_visit_SHA256"] in new_far}
    shared=nearstateeco&farstateeco
    near_matched=[x for x in new_sites if x["outfall_rows_within_10km"]>0 and
                  (x["original_state"],x["original_ecoregion"]) in shared]
    far_matched=[x for x in new_sites if x["NWCA2021_source_site_visit_SHA256"] in new_far and
                 (x["original_state"],x["original_ecoregion"]) in shared]
    shared_states={x["original_state"] for x in near_matched}&{x["original_state"] for x in far_matched}
    gate0=bool(stats["site_near10km_n"]>=25 and stats["nearest_10km_administrative_unique_NPDES_permits_n"]>=15 and
        stats["states_near_10km"]>=8 and stats["matched_far_sites_from_same_state_and_ecoregion_as_a_near_site_n"]>=60 and
        stats["states_with_near_and_far_from_same_state_and_ecoregion"]>=6)
    gate1=bool(len(near_matched)>=MIN_MATCH_NEAR and len(far_matched)>=MIN_MATCH_FAR and
        len(shared_states)>=MIN_STATES and len(shared)>=MIN_STATE_ECO)
    return {
        "schema":"gesui_NWCA2021_568_source_2026EPA_NPDES_exact_EXTERNAL_OUTFALL_only_pre_plant_response_gate_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_EPA_NWCA568_SOURCE_only_SHA256":source_sha(root/PRIOR),
        "original_NPDES_external_source_ZIP_SHA256":source_sha(npdes_zip),
        "source_export_2026_11_08_does_not_certify_2021_real_effluent":True,
        "all_568_frozen_nwca_plant_water_locations_source_complete":True,
        "old_source_registered_any_permitted_feature":base_stats,
        "original_regulatory_permit_feature_SUB_TYPE_DESC_QA":type_stats,
        "strict_2021_plausible_exact_EXTERNAL_OUTFALL_geography":stats,
        "original_old_near10_site_n":len(old_near),
        "new_ExternalOutfall_near10_site_n":len(new_near),
        "original_near10_sites_no_longer_near_ANY_external_outfall":len(old_near-new_near),
        "original_old_far20_site_n":len(old_far),
        "new_ExternalOutfall_far20_site_n":len(new_far),
        "new_source_matched_2021_NWCA_sites":{
            "near10_original_site_n":len(near_matched),
            "far20_original_site_n":len(far_matched),
            "n_original_same_state_ecoregion_strata":len(shared),
            "n_distinct_US_states_with_both_original_source_classes":len(shared_states),
            "state_ecoregion_strata":[a+"|"+b for a,b in sorted(shared)],
            "raw_source_eligible_near10_NWCA_site_SHA256s":sorted(new_near),
            "raw_source_eligible_far20_NWCA_site_SHA256s":sorted(new_far),
            "matched_near10_NWCA_site_SHA256s":sorted(x["NWCA2021_source_site_visit_SHA256"] for x in near_matched),
            "matched_far20_NWCA_site_SHA256s":sorted(x["NWCA2021_source_site_visit_SHA256"] for x in far_matched),
        },
        "frozen_before_corrected_source_type_biological_comparison_GATES":{
            "original_2021_source_geometry_minimum_PASSED":gate0,
            "original_same_state_ecoregion_biodiversity_min_60near_60far_10states_15strata_PASSED":gate1,
            "may_open_new_externalOutfall_source_sensitive_2021_biodiversity_effect":bool(gate0 and gate1),
        },
        "all_original_568_site_hashes_and_external_OUTFALL_source_receipts_no_plant_or_nutrient_outcomes":new_sites,
        "all_actual_2021_NP_both_at_nearby_external_outfall_site_count_NOT_IMPROVED_over_previous_26":True,
        "original_registered_geographic_2026_sources_not_hydrological_connections_to_2021_wetlands":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nwca-source-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    v=geometry_only(root,a.nwca_source_dir,a.npdes_zip)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(v,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "old_near10":v["original_old_near10_site_n"],
        "external_near10":v["new_ExternalOutfall_near10_site_n"],
        "old_far20":v["original_old_far20_site_n"],
        "external_far20":v["new_ExternalOutfall_far20_site_n"],
        "original_source_type":v["original_regulatory_permit_feature_SUB_TYPE_DESC_QA"],
        "new_source":{k:v for k,v in v["strict_2021_plausible_exact_EXTERNAL_OUTFALL_geography"].items()
                      if not isinstance(v,(list,dict))},
        "support":{k:v for k,v in v["new_source_matched_2021_NWCA_sites"].items()
                  if not isinstance(v,list)},
        "stop":v["frozen_before_corrected_source_type_biological_comparison_GATES"]
    },indent=2))

if __name__=="__main__":
    main()
