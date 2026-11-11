#!/usr/bin/env python3
"""2026 EPA permit-year=2016 source-alignment feasibility, H_ALL BLIND.

Frozen posthoc historical permit-date source rule, n50 near/n162 far from
Nov2026 ECHO, not observed 2016 effluent or proof of source water flow.
2011 source gate n12 FAIL: never fit 2011 under historical-date restriction.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd

from gesui_NWCA2011_2016_fixedPOTW2026_near10_far20_source_ONLY import (
    load,sha,assess_sites,ndpes_outfall_source
)
D=Path("results/us_NWCA2021_plant_TNTP_20261110")
S=D/"NWCA2011_2016_ORIGINAL_SITE_VEG_SOURCE_BEFORE_PLANT_OUTCOMES.json"
G=D/"NWCA2011_2016_POTW2026_SOURCE_ONLY_NEAR10_FAR20_2CYCLE.json"
DATE_AUDIT=D/"NWCA2011_2016_EPA_ECHO2026_POTW_ACTUAL_SURVEY_YEAR_PERMIT_SOURCE_QA.json"
WWTP_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"

def run(root,data_dir,zip_file):
    a=json.loads((root/DATE_AUDIT).read_text())
    b=json.loads((root/S).read_text())
    g=json.loads((root/G).read_text())
    if sha(zip_file)!=WWTP_SHA or a["EPA_ECHO_nov2026_NPDES_2021source_snapshot_SHA256"]!=WWTP_SHA:
        raise ValueError("ECHO source changed after historical 2011/2016 permit-date feasibility outcome")
    if a["original_2011_2016_preoutcome_site_source_SHA256"]!=sha(root/S) or a["original_2021plausible_permit_POTW_geo_SHA256"]!=sha(root/G):
        raise ValueError("Prior source frozen 2011/2016 NPDES gate site ids changed")
    if a["year"]["2011"]["historical_year_source_eligibility_comparability_SUPPORT_ONLY"]["posthoc_historical_date_consistency_possible"]:
        raise ValueError("2011 historical outfall source unexpectedly passed")
    if not a["year"]["2016"]["historical_year_source_eligibility_comparability_SUPPORT_ONLY"]["posthoc_historical_date_consistency_possible"]:
        raise ValueError("2016 historical outfall source unexpectedly failed")
    sites,stop=load("2016",root,data_dir,b)
    if stop:raise ValueError("Original 2016 source site frame missing")
    r,stats=ndpes_outfall_source(zip_file)
    effective=pd.to_datetime(r.PERMIT_EFFECTIVE_DATE,errors="coerce")
    expiry=pd.to_datetime(r.PERMIT_EXPIRATION_DATE,errors="coerce")
    term=pd.to_datetime(r.PERMIT_TERMINATION_DATE,errors="coerce")
    eligible=(r.PERMITTED_FEATURE_LATLON &
         effective.le(pd.Timestamp("2016-12-31")) &
         (expiry.isna()|expiry.ge(pd.Timestamp("2016-01-01"))) &
         (term.isna()|term.ge(pd.Timestamp("2016-01-01"))))
    historic=r.loc[eligible].copy()
    if len(historic)!=3341 or int(historic.PERMIT_ID.nunique())!=1792:
        raise ValueError("Previously audited strict 2016-date eligible register source changed")
    receipt,summary=assess_sites(sites,historic,"all_facility_or_feature_proxy")
    nearby=[x for x in receipt if x["outfall_rows_within_10km"]>0]
    far=[x for x in receipt if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or x["nearest_outfall_geodesic_km_within25km_else_null"]>=20]
    match={(x["original_state"],x["original_ecoregion"]) for x in nearby}&{
        (x["original_state"],x["original_ecoregion"]) for x in far}
    output_near=[x for x in nearby if (x["original_state"],x["original_ecoregion"]) in match]
    output_far=[x for x in far if (x["original_state"],x["original_ecoregion"]) in match]
    prior=a["year"]["2016"]["strict_historical_year_plausibility_matched_STATE_ECO_source"]
    if len(output_near)!=prior["near"] or len(output_far)!=prior["far"] or len(match)!=prior["state_ecoregions_with_both_groups"]:
        raise ValueError("2016 historical permit-only source age frame differs from first independent source audit")
    permits={x for r in output_near for x in r["independent_10km_NPDES_permits_SHA256"]}
    if len(permits)!=prior["different_registered_permits_in_near_historical"]:
        raise ValueError("2016 registered sewage permit independent source number differs")
    rows=[]
    for group,vs in (("near10",output_near),("far20",output_far)):
        for r in vs:
            rows.append({
                "original_2016_NWCA_site_visit_hash":r["NWCA2021_source_site_visit_SHA256"],
                "original_state":r["original_state"],
                "original_ecoregion":r["original_ecoregion"],
                "source_original_wetland_class":r["original_wetland_group"],
                "source_original_hydrology_class":r["original_hydro_class"],
                "original_2016_2026_current_permit_asof2016_NEAR10_or_FAR20":group,
                "source_nearest_yearvalid_outfall_km_within25km":r["nearest_outfall_geodesic_km_within25km_else_null"],
                "source_2016_plausible_registered_permit_hashes_10km":r["independent_10km_NPDES_permits_SHA256"],
                "2016_actual_N_P_water_mass_and_hydrological_connection_unconfirmed":True
            })
    rows.sort(key=lambda x:x["original_2016_NWCA_site_visit_hash"])
    return {
        "schema":"gesui_EPA_NWCA2016_2026_ECHO_POTW_permit_effective_asof2016_near10far20_source_BLIND_v01",
        "UTC":datetime.now(timezone.utc).isoformat(),
        "original_2016_veg_SITE_source_SHA256":sha(root/S),
        "original_2016_POTW_2021proxy_source_SHA256":sha(root/G),
        "original_historical_year_permit_age_audit_SHA256":sha(root/DATE_AUDIT),
        "registered_2026_outfall_snapshot_SHA256":WWTP_SHA,
        "2011_historical_permit_source_failure_acknowledged_do_not_fit_old_2011":True,
        "2016_historical_permit_year_source_near10_n":len(output_near),
        "2016_historical_permit_year_source_far20_n":len(output_far),
        "2016_historical_year_source_state_n":len({x["original_state"] for x in output_near}),
        "2016_historical_source_state_ecoregion_n":len(match),
        "2016_historical_near10_unique_regulated_sewage_permits_n":len(permits),
        "original_2016_historical_source_vegmetric_plant_H_ALL_NOT_READ":True,
        "2016_original_source_status_site_hash_receipts":rows,
        "NPDES_regulator_current_permit2026_only_not_verifies_2016_DMR_discharge":True,
        "hist_year_consistency_sensitive_AFTER_original_2021proxy_threecycle_H_ALL_outcome":True
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=run(root,a.source_dir,a.npdes_zip)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in x.items() if k.startswith("2016_historical")},indent=2))
if __name__=="__main__":main()
