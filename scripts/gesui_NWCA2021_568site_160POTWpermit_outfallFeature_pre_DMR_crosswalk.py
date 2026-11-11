#!/usr/bin/env python3
"""Source-only 568 EPA NWCA wetland sites -> exact municipal NPDES permit+feature.

No wetland Shannon/FQAI/N/P outcomes read. Snapshot ECHO 2026 outfall file
byte-SHA, fixed NWCA site sample and original all-568 outfall distance/source
SHA are re-validated. This is ONLY a candidate point-source crosswalk for
historical 2021 DMR; not a directed river flow exposure assignment.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from gesui_NWCA2021_568source_EPA_NPDES_registered_sewage_outfall_geoONLY import (
    original_nwca_source_sites,ndpes_outfall_source,xyz,distances_km)

DIR=Path("results/us_NWCA2021_plant_TNTP_20261110")
OLD=DIR/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
R=10.

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(2**20),b""):h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--nwca-dir",type=Path,required=True)
    ap.add_argument("--npdes-zip",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    prior=json.loads((root/OLD).read_text())
    if not prior["source_precommitted_GEO_eligibility_gate"]["permitted_2021_outfall_distance_association_source_support_PASS"]:
        raise ValueError("Pre-wetland-biology geospatial source gate was not passed")
    if sha(a.npdes_zip)!=prior["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]:
        raise ValueError("ECHO 2026 national outfall source ZIP has CHANGED")
    sites,_=original_nwca_source_sites(root,a.nwca_dir)
    f,info=ndpes_outfall_source(a.npdes_zip)
    f=f[f.PERMITTED_FEATURE_LATLON&f.PERMIT_DATES_PLAUSIBLE_2021].copy()
    if len(f)!=10865 or f.PERMIT_ID.nunique()!=5176:
        raise ValueError("Original nationwide 2026 strict-outfall source graph changed")
    loc=f[["LAT","LON"]].to_numpy(dtype=float)
    tree=cKDTree(xyz(loc[:,0],loc[:,1]))
    chord=2*math.sin((R/6371.0088)/2)
    previous={r["NWCA2021_source_site_visit_SHA256"]:r for r in prior[
        "all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]}
    if len(previous)!=568:raise ValueError("Old 568 outfall point source records changed")
    pairs=[];rows=[]
    for s in sites.itertuples(index=False):
        point=(float(s.LAT_ANALYS),float(s.LON_ANALYS))
        cand=tree.query_ball_point(xyz(np.array([point[0]]),np.array([point[1]]))[0],r=chord+1e-11)
        distances=distances_km(loc[cand],point) if len(cand) else []
        candidates=[(int(k),float(d)) for k,d in zip(cand,distances) if d<=R+1e-12]
        old=previous[s.source_site_visit_SHA256]
        if len(candidates)!=old["outfall_rows_within_10km"]:
            raise ValueError("Strict source number of outfalls <=10km changed from original 2026 ECHO pre-biology gate")
        permit_hashes={hashlib.sha256(str(f.iloc[k].PERMIT_ID).encode()).hexdigest() for k,d in candidates}
        if sorted(permit_hashes)!=old["independent_10km_NPDES_permits_SHA256"]:
            raise ValueError("Original administrative permit identity changed at source site")
        # Same permit+permitted-feature can have multiple exact coordinates;
        # preserve the nearest plausible source point and source row count.
        selected={}
        for k,d in candidates:
            r=f.iloc[k]
            key=(str(r.PERMIT_ID),str(r.FEATURE_ID))
            if key not in selected or d<selected[key]["geodesic_site_to_permitted_feature_km"]:
                selected[key]={
                    "original_public_NPDES_permit_id":key[0],
                    "original_public_NPDES_feature_number":key[1],
                    "2026_outfall_permit_plausibly_2021":True,
                    "geodesic_site_to_permitted_feature_km":d,
                    "source_Wetland_SITE_ID_VISIT_SHA256":s.source_site_visit_SHA256,
                    "source_original_state":str(s.PSTL_CODE),
                    "source_original_NWCA_ecoregion":str(s.AG_ECO9),
                }
        pairs.extend(selected.values())
        rows.append({
            "source_Wetland_SITE_ID_VISIT_SHA256":s.source_site_visit_SHA256,
            "original_state":str(s.PSTL_CODE),
            "original_ecoregion":str(s.AG_ECO9),
            "strict_source_outfalls10km_n":len(candidates),
            "distinct_permitted_feature_within10km_n":len(selected),
            "strict_source_permits10km_n":len(permit_hashes),
            "candidate_source_only_2021_municipal_registration_NOT_measured_emission":True
        })
    if len(rows)!=568 or len({x["source_Wetland_SITE_ID_VISIT_SHA256"] for x in rows})!=568:
        raise ValueError("Original 568 NWCA first visits not recovered")
    near=[x for x in rows if x["strict_source_outfalls10km_n"]>0]
    if len(near)!=128:
        raise ValueError("The original 128 source-near wetlands changed")
    pairs.sort(key=lambda x:(x["original_public_NPDES_permit_id"],
                             x["original_public_NPDES_feature_number"],
                             x["source_Wetland_SITE_ID_VISIT_SHA256"]))
    permitids={x["original_public_NPDES_permit_id"] for x in pairs}
    if len(permitids)!=160:raise ValueError("Original 160 registered source permits near wetlands changed")
    # Same original 32 comparable source state/ecoregion cells, 240 sites
    nearstrata={(x["original_state"],x["original_ecoregion"]) for x in rows if x["strict_source_outfalls10km_n"]>0}
    farstrata={(x["original_state"],x["original_ecoregion"]) for x in prior[
      "all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
      if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or
         x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
    shared=nearstrata&farstrata
    comparable_near=[x for x in near if (x["original_state"],x["original_ecoregion"]) in shared]
    if len(comparable_near)!=103 or len(shared)!=32:
        raise ValueError("Prior pre-outcome matched 103/32 ecology source changed")
    result={
        "schema":"gesui_NWCA2021_568wetland_160municipalPermit_2021DMR_PREOUTCOME_exactPermitFeature_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_ECHO_original_2026_11_outfalls_zip_SHA256":sha(a.npdes_zip),
        "source_previous_568_site_NPDES_geo_qualifier_SHA256":sha(root/OLD),
        "original_2021_NWCA_wetland_source_n":568,
        "original_near10_source_wetlands":len(near),
        "original_near10_in_32_shared_state_ecoregion_source_wetlands":len(comparable_near),
        "source_shared_state_ecoregion_n":len(shared),
        "original_strict_2026_permitted_feature_nationwide_rows":len(f),
        "original_distinct_NPDES_permit_IDs_within10km_of_original_NWCA_wetlands":len(permitids),
        "original_distinct_NPDES_permit_AND_feature_number_within10km":len(set(
            (x["original_public_NPDES_permit_id"],x["original_public_NPDES_feature_number"]) for x in pairs)),
        "registered_2026_site_to_permit_feature_pairs_geodesic_le10km_without_2021_biology":pairs,
        "original_568_NWCA_source_site_municipal_permit_geometries":rows,
        "SOURCE_CLASS_NPDES_POTW_2026_permitted_feature_2021_possible_date_NOT_actual_release":True,
        "2021_real_effluent_TN_TP_FY2021_FY2022_DMR_VALUES_NOT_YET_READ":True,
        "original_NWCA_2021_H_ALL_FQAI_TNTP_results_NEVER_OPENED":True,
        "EPA_FY2021_and_FY2022_DMR_source_filters_precommitted_00600_00665_50050":True,
        "no_hydrological_wetland_receptor_path_known":True
    }
    return result

if __name__=="__main__":
    main()
