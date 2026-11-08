#!/usr/bin/env python3
"""2021 Lycaenidae/Pieridae/Papilionidae *public photo reporting* in 43 original
WorldCover-supported country×1° DE/AT source cells, POSTHOC historical check.

This does NOT test 2021 EEA WWTP exposure. WWTP reported-source data are
from 2022 and may have operated for decades. No photo pixels or observer IDs.
Original 2022 sample was frozen before opening 2021 butterfly family counts.
"""
from __future__ import annotations
import argparse,hashlib,json,math,time,urllib.parse,urllib.request,urllib.error
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path

API="https://api.inaturalist.org/v1/observations"
TAXA={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
COUNTRIES={"DE":7207,"AT":8057}
MIN_2021_FAMILY_PHOTOS_PER_CELL=15
MIN_2021_SUPPORTED_CELLS=30
MIN_2021_DE=18
MIN_2021_AT=9
MIN_2021_COUNTRY_2DEG_BLOCKS=12
EXPECTED_2022={"DE":28,"AT":15}

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()

def request_source(c,taxon):
    lat=float(c["lat_floor"]);lon=float(c["lon_floor"])
    country=c["countryCode"]
    params={
      "taxon_id":taxon,"place_id":COUNTRIES[country],
      "swlat":lat,"swlng":lon,"nelat":lat+1,"nelng":lon+1,
      "d1":"2021-04-01","d2":"2021-10-31",
      "photos":"true","geo":"true","geoprivacy":"open",
      "taxon_geoprivacy":"open","quality_grade":"research",
      "captive":"false","acc_below":500,"per_page":1}
    url=API+"?"+urllib.parse.urlencode(params)
    receipt={"original_public_API_url":url,"http":None,"photo_records_2021":None}
    for attempt in range(4):
        try:
            request=urllib.request.Request(url,headers={
                "User-Agent":"gesui-original_2021_preyear_butterfly_photo_history_posthoc/1.0",
                "Accept":"application/json"})
            with urllib.request.urlopen(request,timeout=35) as r:
                data=r.read(500000)
                obj=json.loads(data)
                if not isinstance(obj.get("total_results"),int):
                    raise ValueError("Official iNaturalist missing total_results")
                receipt.update({"http":r.status,
                      "original_JSON_sha256":hashlib.sha256(data).hexdigest(),
                      "photo_records_2021":int(obj["total_results"])})
                return receipt
        except urllib.error.HTTPError as e:
            receipt.update({"http":e.code,"error":"HTTP"+str(e.code)})
            if e.code not in (429,500,502,503) or attempt==3:break
        except Exception as e:
            receipt["error"]=type(e).__name__+": "+str(e)
            if attempt==3:break
        time.sleep(2*(attempt+1))
    return receipt

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--frozen-landcover-model",type=Path,required=True)
    p.add_argument("--original-2022-photo-source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    past=json.loads(a.frozen_landcover_model.read_text())
    orig=json.loads(a.original_2022_photo_source.read_text())
    if past["status"]!="POSTHOC_LANDCOVER_OUTOFREGION_PREDICTION_AND_PLACEBO_COMPLETED":
        raise ValueError("Prior WorldCover land-adjusted outcome source not final")
    if past["original_2022_family_photo_sha256"]!=sha(a.original_2022_photo_source):
        raise ValueError("Original 2022 family data drift from frozen posthoc source")
    if not orig["2022_photo_share_model_gate"]:
        raise ValueError("Prior original photo model did not pass source/effort quality")
    frozen=past["all_original_same_sample_cells_2022_photo_counts_and_WorldCover_categorical_2021_covariates"]
    if len(frozen)!=43 or Counter(r["countryCode"] for r in frozen)!=Counter(EXPECTED_2022):
        raise ValueError("Original 43 country source cells changed")
    if len(set((r["countryCode"],r["lat_floor"],r["lon_floor"]) for r in frozen))!=len(frozen):
        raise ValueError("Original cell keys not unique")
    rows=[]
    for i,c in enumerate(frozen):
        traces={f:request_source(c,tid) for f,tid in TAXA.items()}
        nums={f:x["photo_records_2021"] for f,x in traces.items()}
        complete=all(isinstance(v,int) for v in nums.values())
        total=sum(nums.values()) if complete else None
        item={
           "countryCode":c["countryCode"],
           "lat_floor":c["lat_floor"],"lon_floor":c["lon_floor"],
           "original_spatial_2degree_block":c["country_2deg_block"],
           "2021_Lycaenidae_photo_n":nums.get("Lycaenidae"),
           "2021_Pieridae_photo_n":nums.get("Pieridae"),
           "2021_Papilionidae_photo_n":nums.get("Papilionidae"),
           "2021_three_family_photo_total":total,
           "original_2021_source_receipts":traces,
           "2021_photo_support_complete":complete,
           "2021_photo_support_min15":bool(total is not None and total>=MIN_2021_FAMILY_PHOTOS_PER_CELL),
           "previous_2022_original_cell_and_source_unchanged":True,
        }
        rows.append(item)
        if i%9==0 or i==len(frozen)-1:
            print("GESUI_2021_HISTORICAL_LYCAENIDAE_PHOTO_REPORTING_SOURCE_QC",
              json.dumps({"queried":i+1,"total":len(frozen),"country":c["countryCode"],
                 "source_city_1deg_lat":c["lat_floor"],"lon":c["lon_floor"],
                 "photos_by_family_2021":nums},sort_keys=True),flush=True)
        time.sleep(.6)
    supported=[x for x in rows if x["2021_photo_support_min15"]]
    by=Counter(x["countryCode"] for x in supported)
    blocks=set(x["original_spatial_2degree_block"] for x in supported)
    eligibility=len(supported)>=MIN_2021_SUPPORTED_CELLS and by.get("DE",0)>=MIN_2021_DE and by.get("AT",0)>=MIN_2021_AT and len(blocks)>=MIN_2021_COUNTRY_2DEG_BLOCKS
    out={
        "schema":"gesui_prior2021_butterfly_photo_composition_source_only_same43_original_2022_WorldCover_sample_v01",
        "date_utc":datetime.now(timezone.utc).isoformat(),
        "POSTHOC_after_positive_2022_source_photo_and_landcover_effect_seen":True,
        "original_2021_WorldCover2022_photo_source_full_SHA256":sha(a.frozen_landcover_model),
        "original_2022_photo_source_full_SHA256":sha(a.original_2022_photo_source),
        "original_2022_DE_AT_WorldCover_fixed_cells":{"DE":28,"AT":15},
        "original_2022_fixed_cell_total":43,
        "original_2021_butterfly_photo_year":2021,
        "2021_photo_window":["2021-04-01","2021-10-31"],
        "original_2021_2022_public_photo_FAMILY_taxon_ids":TAXA,
        "min_prior_2021_family_photos_each":MIN_2021_FAMILY_PHOTOS_PER_CELL,
        "min_prior_2021_supported_countries_1deg_cells":MIN_2021_SUPPORTED_CELLS,
        "min_prior_2021_DE_cells":MIN_2021_DE,
        "min_prior_2021_AT_cells":MIN_2021_AT,
        "min_prior_2021_country2deg_blocks":MIN_2021_COUNTRY_2DEG_BLOCKS,
        "all_2021_API_success":sum(x["2021_photo_support_complete"] for x in rows),
        "all_2021_API_total_queries":len(rows)*len(TAXA),
        "2021_supported_cells_n":len(supported),
        "2021_supported_countries":dict(by),
        "2021_supported_original_country2deg_blocks_n":len(blocks),
        "2021_photo_support_before_historical_comparison_model_gate":eligibility,
        "source_2021_same43_cell_historical_butterfly_photo_records":rows,
        "NO_PHOTO_PIXELS_NO_INDIVIDUAL_OBSERVATION_IDS_OR_PERSONAL_OBSERVERS":True,
        "source_2022_EEA_outfall_is_not_verified_2021_exposure_or_treatment":True,
        "historical_photo_composition_is_not_butterfly_population_history":True,
        "new_2021_photo_data_is_not_independent_preregistered_outcome_test":True,
    }
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_2021_PRIORYEAR_BUTTERFLY_HISTORY_PHOTO_QC_COMPLETE",
      json.dumps({"orig_sites":len(frozen),"photo_2021_supported":len(supported),
         "bycountry":dict(by),"2deg_blocks":len(blocks),"precomparison_gate_pass":eligibility,
         "api_success":out["all_2021_API_success"],"api_queries":out["all_2021_API_total_queries"]},
         sort_keys=True),flush=True)
if __name__=="__main__":main()
