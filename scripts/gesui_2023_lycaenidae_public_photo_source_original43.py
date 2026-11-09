#!/usr/bin/env python3
"""POSTHOC NEW-YEAR independent 2023 Lycaenidae/Pieridae/Papilionidae public
photo-recording source at exact ORIGINAL frozen 43 DE/AT 2022 landcover cells.

Never reads/uses 2023 WWTP source masses (only original 2022); never reads
photograph pixels, observer IDs, butterfly individual counts or behavior.
Fetches only 2023 iNaturalist aggregate 'total_results' photo records and
HTTP response digest after this test's selection/source protocol was frozen.
"""
from __future__ import annotations
import argparse,hashlib,json,time,urllib.parse,urllib.request,urllib.error
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

API="https://api.inaturalist.org/v1/observations"
FAMILIES={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
COUNTRY_PLACES={"DE":7207,"AT":8057}
SEASON=("2023-04-01","2023-10-31")
MIN_PHOTOS=20
MIN_SUPPORT_TOTAL=30
MIN_SUPPORT_DE=18
MIN_SUPPORT_AT=9
MIN_SPATIAL_BLOCKS=12
EXPECTED_PREYEAR={"DE":28,"AT":15}

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def count(row,family,taxon_id):
    la=float(row["lat_floor"]);lo=float(row["lon_floor"])
    params={"taxon_id":taxon_id,"place_id":COUNTRY_PLACES[row["countryCode"]],
      "swlat":la,"swlng":lo,"nelat":la+1,"nelng":lo+1,
      "d1":SEASON[0],"d2":SEASON[1],
      "photos":"true","geo":"true","geoprivacy":"open","taxon_geoprivacy":"open",
      "quality_grade":"research","captive":"false","acc_below":500,"per_page":1}
    url=API+"?"+urllib.parse.urlencode(params)
    source={"original_public_API_url":url,"photo_records_2023":None,"http":None}
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={
             "User-Agent":"gesui-2023-new-year-family-photo-geographic-falsifier/1.0",
             "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=35) as r:
                data=r.read(550000);j=json.loads(data)
                if not isinstance(j.get("total_results"),int):
                    raise ValueError("Official photo query missing total_results")
                source.update({"http":r.status,"original_JSON_sha256":hashlib.sha256(data).hexdigest(),
                               "photo_records_2023":j["total_results"]})
                return source
        except urllib.error.HTTPError as e:
            source.update({"http":e.code,"error":f"HTTP{e.code}"})
            if e.code not in (429,500,502,503) or attempt==3:break
        except Exception as e:
            source.update({"error":type(e).__name__+": "+str(e)})
            if attempt==3:break
        time.sleep(1.6*(attempt+1))
    return source

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-2022-WorldCover-source",type=Path,required=True)
    p.add_argument("--original-2022-family-source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    previous=json.loads(a.original_2022_WorldCover_source.read_text())
    obs2022=json.loads(a.original_2022_family_source.read_text())
    if previous["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED":
        raise ValueError("Original 2022 source was not source-lock complete")
    if previous["original_2022_family_photo_sha256"]!=sha(a.original_2022_family_source):
        raise ValueError("Original 2022 input changed")
    if not obs2022["2022_photo_share_model_gate"]:
        raise ValueError("Original 2022 photo source gate not satisfied")
    cells=previous["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    if len(cells)!=43 or Counter(r["countryCode"] for r in cells)!=Counter(EXPECTED_PREYEAR):
        raise ValueError("Exact original 43-cell 2022 sample not preserved")
    keys=[(r["countryCode"],str(r["lat_floor"]),str(r["lon_floor"])) for r in cells]
    if len(set(keys))!=43:raise ValueError("Duplicate original source/landcover cell key")
    records=[]
    for i,r in enumerate(cells):
        sources={name:count(r,name,taxon_id) for name,taxon_id in FAMILIES.items()}
        counts={f:x["photo_records_2023"] for f,x in sources.items()}
        valid=all(isinstance(x,int) for x in counts.values())
        total=sum(counts.values()) if valid else None
        records.append({
          "countryCode":r["countryCode"],"lat_floor":r["lat_floor"],"lon_floor":r["lon_floor"],
          "original_country_2degree_block":r["country_2deg_block"],
          "2023_Lycaenidae_photo_n":counts["Lycaenidae"],
          "2023_Pieridae_photo_n":counts["Pieridae"],
          "2023_Papilionidae_photo_n":counts["Papilionidae"],
          "2023_three_family_photo_total":total,
          "2023_all_families_query_success":valid,
          "2023_photo_support_ge20":bool(total is not None and total>=MIN_PHOTOS),
          "source_original_2022_reported_P_t_year":r["registered_P_t_year"],
          "public_source_2023_response_receipts":sources,
          "individual_butterfly_abundance_or_adult_puddling_observed":False,
        })
        if i%8==0 or i==len(cells)-1:
            print("GESUI_2023_INDEPENDENT_LYCAENIDAE_PUBLIC_PHOTO_SOURCE_ONLY",
              json.dumps({"progress":i+1,"all43":43,"country":r["countryCode"],
                 "source_2022_P_t_per_year":r["registered_P_t_year"],
                 "2023_family_counts":counts},sort_keys=True),flush=True)
        time.sleep(.7)
    supported=[x for x in records if x["2023_photo_support_ge20"]]
    by=Counter(x["countryCode"] for x in supported)
    blocks=set(x["countryCode"]+"|"+str(x["original_country_2degree_block"]) for x in supported)
    eligible=len(supported)>=MIN_SUPPORT_TOTAL and by.get("DE",0)>=MIN_SUPPORT_DE and by.get("AT",0)>=MIN_SUPPORT_AT and len(blocks)>=MIN_SPATIAL_BLOCKS
    out={
      "schema":"gesui_2023_DE_AT_public_butterfly_family_photo_source_only_original43_WorldCover2021_v01",
      "utc":datetime.now(timezone.utc).isoformat(),
      "POSTHOC_2023_year_check_AFTER_positive_2022_photo_reporting_result":True,
      "original_2022_WorldCover_land_photo_source_sha256":sha(a.original_2022_WorldCover_source),
      "original_2022_EEA_measured_reported_P_source_sha256":obs2022["original_EEA2022_source_manifest_sha256"],
      "original_2022_family_photo_source_sha256":sha(a.original_2022_family_source),
      "exact_source_only_fixed_DE_AT_43_country_1degree_cells":43,
      "source_2023_public_photo_window":list(SEASON),
      "original_family_descendant_taxon_ids":FAMILIES,
      "official_iNat_country_place_ids":COUNTRY_PLACES,
      "min_three_family_photo_records_2023_each_original_cell":MIN_PHOTOS,
      "min_2023_total_supported_cells":MIN_SUPPORT_TOTAL,
      "min_2023_DE":MIN_SUPPORT_DE,
      "min_2023_AT":MIN_SUPPORT_AT,
      "min_2023_country_2degree_blocks":MIN_SPATIAL_BLOCKS,
      "2023_photo_original_API_queries":len(cells)*len(FAMILIES),
      "2023_all_family_API_success_cells":sum(x["2023_all_families_query_success"] for x in records),
      "2023_cells_meeting_photo_support":len(supported),
      "2023_photo_support_by_country":dict(by),
      "2023_photo_support_independent_country_2degree_blocks":len(blocks),
      "2023_outcome_year_source_gate_pass":eligible,
      "source_2023_all_original43_butterfly_photo_count_and_receipts":records,
      "no_2023_WWTP_source_mass_asserted":True,
      "EEA_reported2022_annual_outfall_mass_may_reflect_stable_infrastructure_geography":True,
      "not_individual_butterfly_counts_or_puddling":True,
      "no_photo_pixels_or_individual_observer_identification_downloaded":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2023_PHOTO_SOURCE_RECORD_GATE_FROZEN",
      json.dumps({"original_2022_landcover_cells":43,"supported_2023":len(supported),
         "country":dict(by),"2deg_blocks":len(blocks),"gate":eligible,
         "family_API_success":out["2023_all_family_API_success_cells"]},sort_keys=True),flush=True)

if __name__=="__main__":main()
