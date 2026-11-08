#!/usr/bin/env python3
"""Public EEA 2022 measured-NP Germany/Austria outfall regional index for data coordinators.

Reads ONLY the previously original EEA source-verified public 1,156-outfall
CSV; no butterfly transect positions, observations, counts or habitat data.
Fixed 1° × 1° WGS84 source-only geographic blocks; coarse inventory to
prioritize aggregate protected-transect coverage checks, NOT inferred effects.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

EXPECTED_CSV_SHA256="926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce"
EXPECTED={"DE":630,"AT":526}
DEGREE=1.0

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def tile(v):return math.floor(v/DEGREE)*DEGREE
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    a=p.parse_args()
    if sha(a.source)!=EXPECTED_CSV_SHA256:raise ValueError("Original German/Austrian EEA source manifest changed")
    rows=list(csv.DictReader(a.source.open(encoding="utf8")))
    counts=Counter(r["countryCode"] for r in rows)
    if dict(counts)!=EXPECTED:raise ValueError("Original independently verified country counts changed")
    if len(rows)!=1156:raise ValueError("EEA 2022 source pool unexpectedly changed")
    blocks=defaultdict(list)
    for r in rows:
        country=r["countryCode"]
        la=float(r["discharge_lat_WGS84"]);lo=float(r["discharge_lon_WGS84"])
        n=float(r["original_reported_measured_N_t_year"])
        p=float(r["original_reported_measured_P_t_year"])
        if min(n,p)<=0 or not all(math.isfinite(q) for q in (la,lo,n,p)):
            raise ValueError("Original both-positive source quality violated")
        if not r["source_one_original_FW_outfall"]=="True":raise ValueError("Non FW or multi outlet source included")
        blocks[(country,tile(la),tile(lo))].append((n,p))
    out=[]
    for (country,la,lo),v in sorted(blocks.items()):
        nvals=[t[0] for t in v];pvals=[t[1] for t in v]
        out.append({
          "countryCode":country,"lat_floor":la,"lon_floor":lo,
          "lat_ceiling":la+DEGREE,"lon_ceiling":lo+DEGREE,
          "public_EEA_both_measured_positive_unique_FW_outfalls":len(v),
          "fraction_of_country_original_qualified_sources":len(v)/counts[country],
          "public_EEA_total_reported_N_t_year":sum(nvals),
          "public_EEA_total_reported_P_t_year":sum(pvals),
          "public_EEA_median_per_outfall_N_t_year":statistics.median(nvals),
          "public_EEA_median_per_outfall_P_t_year":statistics.median(pvals),
          "tile_is_not_transect_presence_or_chemistry":True,
        })
    a.out_dir.mkdir(parents=True,exist_ok=True)
    p=a.out_dir/"EEA2022_DE_AT_ORIGINAL_FW_OUTFALLS_IN_ONE_DEGREE_SOURCE_TILES.csv"
    with p.open("w",newline="") as f:
        wr=csv.DictWriter(f,fieldnames=list(out[0]))
        wr.writeheader();wr.writerows(out)
    summary={}
    for country in EXPECTED:
        subset=sorted([r for r in out if r["countryCode"]==country],
                key=lambda r:(-r["public_EEA_both_measured_positive_unique_FW_outfalls"],r["lat_floor"],r["lon_floor"]))
        summary[country]={
          "all_source_1deg_tiles":len(subset),
          "source_original_qualified_fw_outfalls_n":counts[country],
          "source_tiles_ge10_outfalls_n":sum(z["public_EEA_both_measured_positive_unique_FW_outfalls"]>=10 for z in subset),
          "source_count_in_10densest_tiles":sum(z["public_EEA_both_measured_positive_unique_FW_outfalls"] for z in subset[:10]),
          "top10_public_outfall_source_tiles":[{
             "WGS84_lat_floor":z["lat_floor"],"WGS84_lon_floor":z["lon_floor"],
             "original_outfall_count":z["public_EEA_both_measured_positive_unique_FW_outfalls"]
           } for z in subset[:10]]
        }
    receipt={
       "schema":"gesui_2022_original_EEA_DE_AT_public_both_measured_FW_outfalls_1degree_aggregate_source_only_v01",
       "UTC":datetime.now(timezone.utc).isoformat(),
       "source_original_sha256":EXPECTED_CSV_SHA256,
       "source_exact_country_n":EXPECTED,
       "source_spatial_bin_WGS84_degrees":DEGREE,
       "source_tiles_total":len(out),
       "source_tile_CSV_SHA256":sha(p),
       "source_1deg_tile_CSV":p.name,
       "country_source_density_summary":summary,
       "not_independent_country_transect_sample":True,
       "no_butterfly_monitoring_site_location_used":True,
       "no_butterfly_taxa_or_count_outcome_opened":True,
       "no_wastewater_biological_effect_calculated":True,
       "notes":[
          "Tile count is public WWTP original unique freshwater outfall geography, not butterfly sampling effort",
          "Top 10 are highest source count regions independent of any butterfly outcomes",
          "Latitude/longitude 1 degree WGS84 cells have varying km areas by latitude",
          "EIA/EEA reported annual measured N/P load not actual sampled butterfly-substrate nutrient dose",
          "Contact protected transect managers only for aggregate site coverage without disclosing locations"
       ]}
    (a.out_dir/"EEA2022_DE_AT_ORIGINAL_WWTP_SOURCE_1DEG_TILE_QC.json").write_text(
       json.dumps(receipt,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_DE_AT_ORIGINAL_2022_FW_WWTP_SOURCE_ONLY_1DEG_TILE_QC",json.dumps({
       "countries":summary,"total_cells":len(out),"no_butterfly_outcomes":True
    },sort_keys=True),flush=True)
if __name__=="__main__":main()
