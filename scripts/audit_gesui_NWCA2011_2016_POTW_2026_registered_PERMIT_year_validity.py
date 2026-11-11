#!/usr/bin/env python3
"""Retrospective historical PERMIT-DATE plausibility audit on frozen NWCA 2011/16.

The originally fixed 2026 registry restricted to permits plausible in 2021,
which can contain post-2011/16 new permit locations. This is a source DATE
audit AFTER 2011/2016 H_ALL coefficients were seen, not a preregistered
positive treatment validation. No H_ALL, taxa or nutrient values read here.
"""
from __future__ import annotations
import argparse,hashlib,json
from datetime import datetime,timezone
from collections import Counter
from pathlib import Path
import pandas as pd
from gesui_NWCA2011_2016_fixedPOTW2026_near10_far20_source_ONLY import (
    load,sha,assess_sites,ndpes_outfall_source
)

D=Path("results/us_NWCA2021_plant_TNTP_20261110")
SOURCE=D/"NWCA2011_2016_ORIGINAL_SITE_VEG_SOURCE_BEFORE_PLANT_OUTCOMES.json"
ORIGINAL=D/"NWCA2011_2016_POTW2026_SOURCE_ONLY_NEAR10_FAR20_2CYCLE.json"
WWTP_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
def mainrun(root,source_dir,npdes_zip):
    oldsite=json.loads((root/SOURCE).read_text())
    original=json.loads((root/ORIGINAL).read_text())
    if original["original_2011_2016_vegetation_source_only_manifest_SHA256"]!=sha(root/SOURCE):
        raise ValueError("Original historical source SHA has changed")
    if sha(npdes_zip)!=WWTP_SHA:
        raise ValueError("Official Nov8 2026 ECHO source ZIP original SHA changed")
    regs,q=ndpes_outfall_source(npdes_zip)
    if len(regs)!=58231:raise ValueError("Original 2026 registered sewage POTW source cohort changed")
    output={
        "schema":"gesui_NWCA2011_2016_ECHO2026_PERMIT_year_plausibility_postHALL_source_only_QA_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_2011_2016_preoutcome_site_source_SHA256":sha(root/SOURCE),
        "original_2021plausible_permit_POTW_geo_SHA256":sha(root/ORIGINAL),
        "EPA_ECHO_nov2026_NPDES_2021source_snapshot_SHA256":WWTP_SHA,
        "year_specific_permit_definition_ECHO2026_effective_at_or_before_historical_dec31_and_not_expired_terminated_before_jan1":True,
        "AUDIT_AFTER_2011_2016_source_proxy_signs_ALREADY_SEEN":True,
        "2026_register_NOT_proof_of_actual_historical_discharge":True,
        "year":{}
    }
    for y in ("2011","2016"):
        d,stop=load(y,root,source_dir,oldsite)
        if stop or d is None:raise ValueError("Original 2011/2016 source site gate missing")
        r=regs.copy()
        effective=pd.to_datetime(r.PERMIT_EFFECTIVE_DATE,errors="coerce")
        expiry=pd.to_datetime(r.PERMIT_EXPIRATION_DATE,errors="coerce")
        term=pd.to_datetime(r.PERMIT_TERMINATION_DATE,errors="coerce")
        status=(r.PERMITTED_FEATURE_LATLON &
            effective.le(pd.Timestamp(f"{y}-12-31")) &
            (expiry.isna()|expiry.ge(pd.Timestamp(f"{y}-01-01"))) &
            (term.isna()|term.ge(pd.Timestamp(f"{y}-01-01"))))
        historical=r.loc[status].copy()
        rows,geo=assess_sites(d,historical,"all_facility_or_feature_proxy")
        # Historical dataframe has already been STRICT filtered to actual
        # permitted-feature geo and year condition; downstream helper uses
        # only site and Haversine distances, no response.
        near={x["NWCA2021_source_site_visit_SHA256"] for x in rows if x["outfall_rows_within_10km"]>0}
        far={x["NWCA2021_source_site_visit_SHA256"] for x in rows if x["nearest_outfall_geodesic_km_within25km_else_null"] is None or x["nearest_outfall_geodesic_km_within25km_else_null"]>=20}
        strata={(x["original_state"],x["original_ecoregion"]) for x in rows if x["outfall_rows_within_10km"]>0}&{
            (x["original_state"],x["original_ecoregion"]) for x in rows if x["NWCA2021_source_site_visit_SHA256"] in far
        }
        nearMatched={x["NWCA2021_source_site_visit_SHA256"] for x in rows if x["NWCA2021_source_site_visit_SHA256"] in near and (x["original_state"],x["original_ecoregion"]) in strata}
        farMatched={x["NWCA2021_source_site_visit_SHA256"] for x in rows if x["NWCA2021_source_site_visit_SHA256"] in far and (x["original_state"],x["original_ecoregion"]) in strata}
        originalnear={x["site_SHA256"] for x in original["source_2011_2016_years"][y]["original_historical_source_matched_site_date_weight_and_outfall_receipts"] if x["source_group"]=="near10"}
        matching=originalnear & near
        permits={q for x in rows if x["NWCA2021_source_site_visit_SHA256"] in nearMatched for q in x["independent_10km_NPDES_permits_SHA256"]}
        output["year"][y]={
            "original_2021_permit_proxy_2011_16_SOURCE_near10_matched_n":len(originalnear),
            "year_specific_ECHO2026_historical_date_plausible_source_points_n":len(historical),
            "year_specific_ECHO2026_historical_date_plausible_unique_permits_n":int(historical.PERMIT_ID.nunique()),
            "historical_year_plausible_source_sites_near10_n":len(near),
            "historical_year_plausible_source_sites_far20_n":len(far),
            "strict_historical_year_plausibility_matched_STATE_ECO_source":{
                "near":len(nearMatched),"far":len(farMatched),
                "states_with_both_groups":len({x["original_state"] for x in rows if x["NWCA2021_source_site_visit_SHA256"] in nearMatched}),
                "state_ecoregions_with_both_groups":len(strata),
                "different_registered_permits_in_near_historical":len(permits)
            },
            "original_PRE2021_source_near_sites_also_near_year_valid_source_n":len(matching),
            "source_2021_proximate_near_sample_historical_date_overlap_proportion":len(matching)/len(originalnear),
            "historical_year_source_eligibility_comparability_SUPPORT_ONLY":{
                "min_near":50,"min_far":60,"min_states":10,
                "min_geographic_state_eco":20,"min_registered_permits":15,
                "posthoc_historical_date_consistency_possible":(
                    len(nearMatched)>=50 and len(farMatched)>=60 and
                    len(strata)>=20 and len(permits)>=15)
            },
            "actual_2011_2016_WWTP_wastewater_NP_and_wetland_contact_UNMEASURED":True
        }
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--npdes-zip",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    x=mainrun(root,a.source_dir,a.npdes_zip)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(x,sort_keys=True,indent=2)+"\n")
    print(json.dumps(x["year"],indent=2))
if __name__=="__main__":main()
