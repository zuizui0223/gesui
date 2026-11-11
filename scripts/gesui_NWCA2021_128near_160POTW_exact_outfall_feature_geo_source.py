#!/usr/bin/env python3
"""Lock EXACT geo-near NPDES permit+outfall feature IDs at 568 NWCA sites.

Original weekly 2026 EPA NPDES zip SHA and original four 2021 NWCA source
CSVs are independently checked. No 2021 DMR, original NWCA H_ALL, FQAI or
TN/TP outcomes read in this stage. This permits a strict later 2021 DMR
feature-ID matching instead of misleading permit-only attribution.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from gesui_NWCA2021_568source_EPA_NPDES_registered_sewage_outfall_geoONLY import (
    original_nwca_source_sites,ndpes_outfall_source,xyz,distances_km)

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
FIX=D/"EPA_NWCA2021_160_SEWAGE_PERMITS_BEFORE_DMR2021_SOURCE.json"
GEO=D/"EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json"
SOURCE_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def textsha(s):return hashlib.sha256(str(s).encode()).hexdigest()
def encode(a,b):
    return textsha(a+"|"+b)

def run(root,nwca_dir,echo):
    perm=json.loads((root/FIX).read_text())
    geo=json.loads((root/GEO).read_text())
    if perm["original_EPA_NWCA2021_568source_2026_NPDES_geo_SHA256"]!=sha(root/GEO):
        raise ValueError("Original preregistered NPDES near160 source changed")
    if sha(echo)!=SOURCE_SHA:
        raise ValueError("Original 2026-11-08 ECHO NPDES outfall ZIP not byte-identical")
    sites,_=original_nwca_source_sites(root,nwca_dir)
    source,source_qc=ndpes_outfall_source(echo)
    strict=source[source.PERMITTED_FEATURE_LATLON & source.PERMIT_DATES_PLAUSIBLE_2021].copy()
    fixed=set(perm["exact_source_frozen_near10_registered_municipal_permit_ID_SHA256"])
    strict=strict[strict.PERMIT_ID.map(textsha).isin(fixed)].copy()
    if len(sites)!=568 or strict.empty:raise ValueError("Expected original study source points absent")
    tree=cKDTree(xyz(strict.LAT.to_numpy(),strict.LON.to_numpy()))
    chord=2*np.sin((10./6371.0088)/2)
    pos=strict[["LAT","LON"]].to_numpy()
    rows=[]
    for record in sites.itertuples(index=False):
        id=record.source_site_visit_SHA256
        center=(float(record.LAT_ANALYS),float(record.LON_ANALYS))
        cand=tree.query_ball_point(xyz(np.array([center[0]]),np.array([center[1]]))[0],r=chord+1e-11)
        if not cand:continue
        distances=distances_km(pos[cand],center)
        for idx,km in zip(cand,distances):
            if km>10.00000001:continue
            s=strict.iloc[idx]
            permit_id=str(s.PERMIT_ID).strip().upper()
            feature_id=str(s.FEATURE_ID).strip().upper()
            rows.append({
                "NWCA2021_site_visit_SHA256":id,
                "source_original_US_state":str(record.PSTL_CODE),
                "source_original_1deg_ecoregion":str(record.AG_ECO9),
                "NPDES_registered_permit_SHA256":textsha(permit_id),
                "NPDES_exact_permit_feature_SHA256":encode(permit_id,feature_id),
                "original_PERM_FEATURE_NMBR_only_SHA256":textsha(feature_id),
                "source_geodesic_outfall_to_wetland_km":float(km),
                "original_2026_official_LATLONG_TYPE_Permitted_Feature":True,
                "original_permit_dates_plausible2021_not_actual_2021_release":True,
            })
    near={r["NWCA2021_site_visit_SHA256"] for r in rows}
    permits={r["NPDES_registered_permit_SHA256"] for r in rows}
    if len(near)!=128 or permits!=fixed:
        raise ValueError(f"Exact outfall-source vs existing pre-DMR original geometry mismatch {len(near)} 128, {len(permits)} 160")
    old=perm["original_site_near10_permits_2021outfallPlausible_from_2026_registry"]
    expected={x["NWCA_site_visit_SHA256"]:set(x["permitted_sewage_NPDES_ID_SHA256"]) for x in old}
    actual={key:{x["NPDES_registered_permit_SHA256"] for x in rows if x["NWCA2021_site_visit_SHA256"]==key} for key in near}
    if expected!=actual:
        raise ValueError("New actual original 2026 feature-level source mapping disagrees with frozen 128 original permit lists")
    features={r["NPDES_exact_permit_feature_SHA256"] for r in rows}
    out={
        "schema":"gesui_NWCA2021_original128geoNear_registered_POTW2026_permitFEATURE_hash_source_only_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "source_original_2026_ECHO_NPDES_ZIP_SHA256":sha(echo),
        "source_NWCA_original568_POTW_geo_result_SHA256":sha(root/GEO),
        "source_original160_permit_DMR_manifest_SHA256":sha(root/FIX),
        "original_568_probability_wetland_sites_with_original_2021_plant_water_source":568,
        "original_128_near10_wetland_site_location_count":len(near),
        "original_160_distinct_administrative_registered_sewage_permits":len(permits),
        "original_near10_exact_permit_feature_ID_count":len(features),
        "original_site_outfall_candidate_feature_pairs":len(rows),
        "original_source_case_site_to_exact_permit_feature_sha256_and_distance":sorted(
            rows,key=lambda x:(x["NWCA2021_site_visit_SHA256"],x["NPDES_exact_permit_feature_SHA256"])),
        "source_only_original_NWCA2021_2026_current_permit_outfalls_NOT_FY2021_release":True,
        "no_2021_actual_DMR_or_wetland_plant_nutrient_outcome_opened":True
    }
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--nwca-source-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();root=Path(__file__).resolve().parents[1]
    out=run(root,a.nwca_source_dir,a.npdes_zip)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in out.items() if k in (
        "original_128_near10_wetland_site_location_count",
        "original_160_distinct_administrative_registered_sewage_permits",
        "original_near10_exact_permit_feature_ID_count","original_site_outfall_candidate_feature_pairs")},indent=2))
if __name__=="__main__":main()
