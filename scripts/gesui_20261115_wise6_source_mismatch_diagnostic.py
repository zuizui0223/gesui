#!/usr/bin/env python3
"""Diagnose how real WISE6 same-EEA-station records fail source-only TNTP gate.

Posthoc SOURCE diagnostic sample: 4 source-reported-positive stations per year,
selected from the already-published 32-site pilot by existence of source WQ
rows (NOT by any chemical value, environmental response or positive sign).
Never treat zero same-year records as zero actual concentration or sewage.
"""
import collections
import json
import os
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0,str(Path(__file__).resolve().parent))
from gesui_20261115_wise6_native_site_np_pilot import read_site_source

TARGET_CODES={"TN":"EEA_31615-01-7","TP":"CAS_7723-14-0"}

def category_counts(rows,station,year):
    country,site,scheme=station.split("|")
    exact=[r for r in rows if str(r.get("monitoringSiteIdentifier"))==site
           and str(r.get("monitoringSiteIdentifierScheme"))==scheme]
    years=collections.Counter(str(r.get("phenomenonTimeReferenceYear")) for r in exact)
    sameyear=[r for r in exact if str(r.get("phenomenonTimeReferenceYear"))==str(year)]
    river=[r for r in sameyear if r.get("parameterWaterBodyCategory")=="RW"]
    water=[r for r in river if r.get("procedureAnalysedMatrix")=="W"]
    counters={}
    for label,records in [("exact_site_any_year",exact),("same_year_all_parameters",sameyear),
                          ("same_year_river",river),("same_year_river_water_matrix",water)]:
        counters[label]=len(records)
    nutrients={}
    for name,code in TARGET_CODES.items():
        target=[r for r in water if r.get("observedPropertyDeterminandCode")==code]
        alltime=[r for r in exact if r.get("observedPropertyDeterminandCode")==code]
        nutrients[name]={"same_year_RWW_exact":len(target),"all_years_exact":len(alltime),
                          "same_year_units":dict(collections.Counter(str(r.get("resultUom")) for r in target)),
                          "same_year_sample_count":[r.get("resultNumberOfSamples") for r in target],
                          "same_year_numeric_present":sum(r.get("resultMeanValue") is not None for r in target),
                          "same_year_below_loq_flags":[r.get("resultQualityMeanBelowLOQ") for r in target]}
    near=[(str(r.get("observedPropertyDeterminandCode")),str(r.get("observedPropertyDeterminandLabel")))
          for r in water if any(k in str(r.get("observedPropertyDeterminandLabel","")).lower()
                                for k in ("nitrogen","phosphorus","nitrate","phosphate","azot","fosfor"))]
    out={"native_site":station,"target_year":year,"row_counts":counters,
         "source_year_counts":dict(years),
         "same_year_river_water_parameter_examples":dict(collections.Counter(
             (str(r.get("observedPropertyDeterminandCode")),str(r.get("observedPropertyDeterminandLabel")))
             for r in water).most_common(16)),
         "nutrient_like_parameter_codes":sorted(set(near))[:30],
         "target_total_TNTP":nutrients}
    return out

def main():
    source=Path(os.environ["ORIGINAL_SELECTED_PILOT_DIR"])
    selected=json.loads((source/"SOURCE_POSITIVE_NATIVE_SITE_SELECTION.json").read_text())
    receipts=json.loads((source/"SOURCE_FETCH_RECEIPTS.json").read_text())
    assert len(selected)==32 and sum(r["upstream_report_positive"] for r in selected)==24
    chosen=[]
    for year in (2018,2022):
        candidates=sorted([r for r in selected if r["year"]==year and r["upstream_report_positive"]
            and len(receipts[r["native_site"]])>=1
            and sum(v.get("n",0) for v in receipts[r["native_site"]])>0],
            key=lambda r:r["native_site"])
        if len(candidates)<4:raise RuntimeError("Not enough original WQ-site accessible positive reports")
        chosen.extend(candidates[:4])
    assert len(chosen)==8
    result={}
    with ThreadPoolExecutor(max_workers=2) as executor:
        futures={executor.submit(read_site_source,q["native_site"]):q for q in chosen}
        for future in as_completed(futures):
            original=futures[future]
            raw=future.result()
            station=original["native_site"]
            info={"reported_upstream_source_positive_original":True,
                "original_H7":original["H7"],"original_period":original["period"],
                "status":raw["status"],"n_fetched_source_rows":len(raw["source_rows"]),
                "source_receipts":raw["receipts"]}
            if raw["status"]=="complete":
                info.update(category_counts(raw["source_rows"],station,original["year"]))
            result[station]=info
            print("WISE6_STAGE",json.dumps({"site":station,"year":original["year"],
                 "status":raw["status"],"original_rows":info["n_fetched_source_rows"],
                 "stages":info.get("row_counts"),
                 "TN":info.get("target_total_TNTP",{}).get("TN"),
                 "TP":info.get("target_total_TNTP",{}).get("TP")},sort_keys=True),flush=True)
    complete=[x for x in result.values() if x["status"]=="complete"]
    summary={
      "schema":"GESUI_WISE6_EEA_SOURCE_SITE_TNTP_MISMATCH_STAGING_V1",
      "selected_original_upstream_report_positive_sites":8,
      "complete_source_site_queries":len(complete),
      "target_year_had_any_wq_records":sum(z.get("row_counts",{}).get("same_year_all_parameters",0)>0 for z in complete),
      "target_year_river_water_records":sum(z.get("row_counts",{}).get("same_year_river_water_matrix",0) for z in complete),
      "target_year_river_water_total_n_rows":sum(z.get("target_total_TNTP",{}).get("TN",{}).get("same_year_RWW_exact",0) for z in complete),
      "target_year_river_water_total_p_rows":sum(z.get("target_total_TNTP",{}).get("TP",{}).get("same_year_RWW_exact",0) for z in complete),
      "posthoc_source_availability_diagnostics_only":True,
      "no_ecological_response_or_wwtp_causal_effect_fitted":True,
      "not_safely_infer_other_124_positive_sites":True,
      "not_false_water_quality_zeros":True
    }
    dest=Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_20261115_wise6_stage_diagnostics"
    dest.mkdir(parents=True,exist_ok=True)
    (dest/"SOURCE_FAILURE_REASON_DIAGNOSTIC.json").write_text(json.dumps({"summary":summary,"sites":result},indent=2,sort_keys=True))
    print("STAGE_DIAGNOSTIC",json.dumps(summary,sort_keys=True),flush=True)
if __name__=="__main__":main()
