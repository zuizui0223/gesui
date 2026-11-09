#!/usr/bin/env python3
"""Original 2022 EEA FW WWTP source origin reaches STRICTLY BETWEEN national WISE6
real 2022 river TP monitoring site pairs, before ANY numeric TP values.

Replays exactly the existing 79 original butterfly-photo-region pair source
classifications to assert pipeline stability, then classifies all new 115
NATIONAL environmental candidate pairs. Original 2022 REST response SHA must
match prior archived source-only report; original 50km, 500m, 8/8 controls
unchanged. No water TP magnitude or butterfly photo outcome accessed.
"""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,geopandas as gpd
import gesui_2022_WISE6_qualifying_EEA_WWTP_exact_between_stations_sourceonly as priorpairs
import gesui_2022_WISE6_TP_original_station_vs_directed_WWTP_sourceonly as link

SOURCE_OLD_SHA="9c07ce9de7f5104b47052c7598a943df405a4b83"
GRAPH_SHA="5e6327c6dd34a81bd7768d840f25c0569d4361b7e078065a97214563bd01b32f"
EXPECTED_SOURCES=4057
EXPECTED_SNAPPED=2967
MAX_KM=50
MIN_POS=8
MIN_ZERO=8
MIN_BLOCKS=4

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def translate(x):
    return {
      "countryCode":x["countryCode"],
      "upstream_station_SHA256":x["upstream_station_SHA256"],
      "downstream_station_SHA256":x["downstream_station_SHA256"],
      "upstream_original_river_reach_id":x["upstream_original_reach_HYRIV_ID"],
      "downstream_original_river_reach_id":x["downstream_original_reach_HYRIV_ID"],
      "upstream_country2degree_block":x["upstream_2degree_region"],
      "downstream_country2degree_block":x["downstream_2degree_region"],
      "full_reach_path_to_downstream_reach_start_km":x["full_reach_path_approx_km"],
      "original_directed_reach_steps":x["directed_hyriv_segments"],
      "NO_numeric_2022_TP_read":True,
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--national-original-pair-stageA",type=Path,required=True)
    p.add_argument("--prior-original-79-pair-stageB",type=Path,required=True)
    p.add_argument("--prior-WWTP-187-sourceonly",type=Path,required=True)
    p.add_argument("--river-parquet",type=Path,required=True)
    p.add_argument("--pinned-2022-UWWTD-script",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    nat=json.loads(a.national_original_pair_stageA.read_text())
    oldpairs=json.loads(a.prior_original_79_pair_stageB.read_text())
    oldsource=json.loads(a.prior_WWTP_187_sourceonly.read_text())
    if not nat["source_only_predeclared_10_pairs_DE3_AT3_fourblocks_gate"]:
        raise ValueError("Entire original national source sampling station pair gate failed")
    if nat["original_thirdparty_HydroRIVERS_EU_SHA256"]!=GRAPH_SHA or sha(a.river_parquet)!=GRAPH_SHA:
        raise ValueError("Original EU river graph is not SHA-pinned mirror")
    if oldpairs["predeclared_stageB_real_WWTP_intervening_vs_sourcefree_pair_gate_pass"]:
        raise ValueError("Original photo-cell 2022 WWTP pair gate unexpectedly passed; historical result changed")
    if gitblob(a.pinned_2022_UWWTD_script)!=SOURCE_OLD_SHA:
        raise ValueError("Original public 2022 EEA UWWTD model source code SHA changed")
    import sys
    sys.path.insert(0,str(a.pinned_2022_UWWTD_script.parent))
    import gesui_eea_nutrient_load_river_crosswalk as eea
    WW,wreceipts=eea.download_layer(eea.PLANT,eea.BOTH_WHERE,eea.PLANT_FIELDS,5492)
    DP,dreceipts=eea.download_layer(eea.OUTFALL,"1=1",eea.OUT_FIELDS,24416)
    if wreceipts+dreceipts!=oldsource["original_source_API_page_response_SHA_receipts"]:
        raise ValueError("EEA 2022 original REST source 31 page hashes changed")
    WW=eea.normalize_key(WW);DP=eea.normalize_key(DP)
    WW=WW[WW.compound_join_key.notna()].copy()
    DP=DP[DP.compound_join_key.notna()].copy()
    nWW=WW.compound_join_key.value_counts()
    nDP=DP.compound_join_key.value_counts()
    WW=WW[WW.compound_join_key.map(nWW).eq(1)].copy()
    WW["outlet_n"]=WW.compound_join_key.map(nDP).fillna(0).astype(int)
    source=WW[WW.outlet_n.eq(1)].merge(DP.drop_duplicates("compound_join_key"),
          on="compound_join_key",how="inner",suffixes=("_WWTP","_DCP"),validate="one_to_one")
    for k in ("dcpLatitude","dcpLongitude","uwwNDischargeMeasured","uwwPDischargeMeasured"):
        source[k]=pd.to_numeric(source[k],errors="coerce")
    source=source[source.dcpLatitude.between(44,64)&source.dcpLongitude.between(-12,26)]
    source=source[source.dcpWaterBodyType.astype("string").str.strip().str.upper().eq("FW")].copy()
    if len(source)!=EXPECTED_SOURCES:raise ValueError("Official 2022 public original unique FW outlets not 4057")
    g=gpd.read_parquet(a.river_parquet).to_crs("EPSG:4326")
    g=g.rename(columns={k:k.upper() for k in g.columns if k!=g.geometry.name})
    g=g.cx[-12:26,44:64].copy()
    for k in ("HYRIV_ID","NEXT_DOWN","LENGTH_KM"):
        g[k]=pd.to_numeric(g[k],errors="coerce")
    g=g[g.HYRIV_ID.notna()].copy()
    graph=priorpairs.topo(g)
    if len(graph)!=199545:raise ValueError("Expected original 199545 directed study network reaches changed")
    source["source_row_id"]=np.arange(len(source),dtype=int)
    source["lat"]=source.dcpLatitude
    source["lon"]=source.dcpLongitude
    snaps,quality=link.snap_unambiguous(source,g,"EEA2022_original_FW")
    if len(snaps)!=EXPECTED_SNAPPED:raise ValueError("Original EEA 2022 source reach snap count changed")
    indexes=[];routes=[]
    for r in snaps.itertuples(index=False):
        reach=int(r.HYRIV_ID)
        massP=float(r.uwwPDischargeMeasured);massN=float(r.uwwNDischargeMeasured)
        if not (massP>0 and massN>0 and math.isfinite(massP) and math.isfinite(massN)):continue
        if reach not in graph:raise ValueError("Original EEA WWTP mapped to river missing topo")
        indexes.append((reach,massP,massN,hashlib.sha256(str(r.compound_join_key).encode()).hexdigest()))
        routes.append(priorpairs.path_from(reach,graph))
    new=[]
    for x in nat["all_original_national_source_only_directed_same_river_station_pair_candidates"]:
        v=translate(x)
        v["real_exact_source_between"]=priorpairs.source_appears_between(v,routes,indexes,graph)
        new.append(v)
    old={ (r["upstream_station_SHA256"],r["downstream_station_SHA256"]):
             r["real_exact_source_between"]
             for r in oldpairs["original_all_sourceonly_stageA_candidate_pairs_with_exact_outlet_interposition"]}
    errors=[]
    for r in new:
        k=(r["upstream_station_SHA256"],r["downstream_station_SHA256"])
        if k not in old:continue
        if old[k]!=r["real_exact_source_between"]:errors.append(k)
    if errors:
        raise ValueError(f"Original earlier 79 photo-cell WWTP source-between assignments failed exact replay on {len(errors)} station pairs")
    if len(old)!=79:raise ValueError("Original 79-site source outcome replay reference unexpectedly changed")
    newsubset={ (r["upstream_station_SHA256"],r["downstream_station_SHA256"]) for r in new}
    if set(old)-newsubset:
        raise ValueError("New national source cohort failed to retain all 79 original photo-region station pairs")
    chosen,qc=priorpairs.pick_strata(new)
    both=chosen["control"]+chosen["intervening_source"]
    c=Counter(r["countryCode"] for r in both)
    blockset=set(r["countryCode"]+"|"+r["upstream_country2degree_block"] for r in both)
    gate=(len(chosen["control"])>=MIN_ZERO and len(chosen["intervening_source"])>=MIN_POS and
           c["DE"]>=1 and c["AT"]>=1 and len(blockset)>=MIN_BLOCKS)
    report={
      "schema":"gesui_2022_WISE6_national_DE_AT_original_EEA_FW_WWTP_between_same_river_stationpairs_sourceonly_STAGEB_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "original_national_2022_TP_sourceonly_station_pair_stageA_sha256":sha(a.national_original_pair_stageA),
      "original_prior_photo_cell_79_pair_WWTP_reach_stageB_sha256":sha(a.prior_original_79_pair_stageB),
      "original_2022_187_public_WWTP_station_sourceonly_sha256":sha(a.prior_WWTP_187_sourceonly),
      "original_EEA_2022_REST_31pages_sha_receipts_REPRODUCED":True,
      "original_pinned_thirdparty_HydroRIVERS_source_SHA256":GRAPH_SHA,
      "NOT_official_river_topology_independently_verified":True,
      "original_EEA2022_Fw_one_unique_outlet_source_count":len(source),
      "original_EEA2022_source_public_snap_500m_n":len(snaps),
      "replayed_all_original_79_photo_cell_2022_WWTP_source_interposition_records_exactly":True,
      "all_national_2022_same_river_source_only_pairs_with_verified_outlet_between_reaches":new,
      "nationwide_source_only_pair_exposure_counts":qc,
      "source_only_selected_disjoint_8_source_and_8_clean_registered_zero_pair_attempt":chosen,
      "selected_disjoint_paired_country_counts":dict(c),
      "selected_disjoint_original_country2degree_region_blocks":len(blockset),
      "original_national_WWTP_between_8_and_clean_zero_8_source_gate_pass":gate,
      "source_registered_outfall_total_P_is_not_waterTP_concentration":True,
      "no_WISE6_TP_numeric_measured_values_or_photo_observer_data_opened":True,
      "2022_national_environmental_population_is_NOT_Lycaenidae_butterfly_response":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,sort_keys=True,ensure_ascii=False,indent=2)+"\n")
    print("GESUI_2022_NATIONAL_RIVER_TP_REAL_2022_EEA_WWTP_BETWEEN_STATION_SOURCE_GATE",
       json.dumps({"national_all_raw_pairs":len(new),
           "source_exposure_counts":qc,
           "source_positive_disjoint_selected":len(chosen["intervening_source"]),
           "no_registered_between_and_no_same_reach_ambiguity_disjoint_controls":len(chosen["control"]),
           "bycountry":dict(c),"country2deg_blocks":len(blockset),
           "gate":gate,
           "all_original_79_EEA_WWTP_source_pair_status_replay_passed":True,
           "no_numeric_TP_values_used":True},sort_keys=True),flush=True)

if __name__=="__main__":main()
