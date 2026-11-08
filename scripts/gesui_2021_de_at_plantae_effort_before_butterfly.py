#!/usr/bin/env python3
"""Freeze 53 EEA source tiles with independent PRIOR YEAR 2021 Plantae photo effort.

Absolutely no 2022 Lycaenidae/Pieridae/Papilionidae queries here.
Source EEA 2022 public one-FW-outfall N/P fields are already SHA-locked.
Year 2021 Plantae photo reporting is a separate observational effort proxy,
not vegetation or nitrate exposure.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time,urllib.error,urllib.parse,urllib.request
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
SOURCE_SHA="926ab310904be048d5d36bbfebd5eabfe0b318279800142f2e2a86bee924f4ce"
SOURCE_TILE_SHA="PLACEHOLDER_SHA_FROM_PREVIOUS_SOURCE_TILE_QC"
API="https://api.inaturalist.org/v1"
SOURCE_COUNTS={"DE":630,"AT":526}
EXPECTED_SOURCE_TILES={"DE":32,"AT":21}
MIN_PLANT=40
MIN_TOTAL=25
MIN_DE=12
MIN_AT=8
MIN_GEO_BLOCKS=10
LATITUDE_DAYS=2021
MAX_QUERIES=70
REQUESTS=0

def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda:f.read(1048576),b""):h.update(chunk)
    return h.hexdigest()

def get(path,params=None):
    global REQUESTS
    REQUESTS+=1
    if REQUESTS>MAX_QUERIES:raise RuntimeError("Public API call bound exceeded")
    url=API+path+("?" + urllib.parse.urlencode(params) if params else "")
    info={"url":url}
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"gesui-2021-plants-before-2022-DE-AT-lycaenid-macro/0.1",
                "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=36) as r:
                raw=r.read(600000)
                j=json.loads(raw)
                info.update({"http":r.status,"response_SHA256":hashlib.sha256(raw).hexdigest(),"response_bytes":len(raw)})
                return info,j
        except urllib.error.HTTPError as e:
            info.update({"http":e.code,"error":str(e)})
            if e.code not in (429,500,502,503) or attempt==2:break
        except Exception as e:
            info["error"]=type(e).__name__+": "+str(e)
            if attempt==2:break
        time.sleep(1.8*(attempt+1))
    return info,None

def resolve_taxon():
    receipt,j=get("/taxa/autocomplete",{"q":"Plantae","per_page":20})
    eligible=[x for x in (j or {}).get("results",[]) if x.get("name")=="Plantae" and x.get("rank")=="kingdom"]
    if len(eligible)!=1:raise ValueError("Original iNat Plantae kingdom not unique")
    return int(eligible[0]["id"]),{"taxon_name":"Plantae","id":eligible[0]["id"],"receipt":receipt}
def place(name):
    receipt,j=get("/places/autocomplete",{"q":name,"per_page":30})
    candidates=[{"name":x.get("name"),"id":x.get("id"),
       "place_type":x.get("place_type"),"place_type_name":x.get("place_type_name")} for x in (j or {}).get("results",[])]
    exact=[x for x in candidates if str(x["name"] or "").casefold()==name.casefold()]
    official=[x for x in exact if (str(x.get("place_type_name") or "").casefold()=="country"
               or x.get("place_type")==12)]
    if len(official)==1:ans=official[0]
    elif len(exact)==1:ans=exact[0]
    else:
        raise ValueError(f"Ambiguous or absent official place {name}, metadata only: {candidates[:10]}")
    return int(ans["id"]),{"official_place":ans,"receipt":receipt,"candidate_count":len(candidates)}

def band(n):
    return "01_04" if n<=4 else "05_19" if n<=19 else "20plus"

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-manifest",required=True,type=Path)
    p.add_argument("--source-tiles",required=True,type=Path)
    p.add_argument("--source-tiles-qc",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    if sha(a.source_manifest)!=SOURCE_SHA:raise ValueError("Original 1156 DE/AT WWTP source manifest SHA drift")
    orig=json.loads(a.source_tiles_qc.read_text())
    if orig["source_original_sha256"]!=SOURCE_SHA:
        raise ValueError("Source-only geographic tiles do not correspond to original EEA manifest")
    if sha(a.source_tiles)!=orig["source_tile_CSV_SHA256"]:
        raise ValueError("Originally frozen 1 degree geographic source tile CSV changed")
    tiles=list(csv.DictReader(a.source_tiles.open()))
    countby=Counter(x["countryCode"] for x in tiles)
    if dict(countby)!=EXPECTED_SOURCE_TILES:
        raise ValueError("Original source-only 53 tiles unexpectedly changed")
    for c in SOURCE_COUNTS:
        if sum(int(x["public_EEA_both_measured_positive_unique_FW_outfalls"]) for x in tiles if x["countryCode"]==c)!=SOURCE_COUNTS[c]:
            raise ValueError("Original country source count allocation changed")
    plants,tmeta=resolve_taxon()
    ids={};places={}
    for c,name in (("DE","Germany"),("AT","Austria")):
        ids[c],places[c]=place(name)
    rows=[]
    for i,t in enumerate(tiles):
        iso=t["countryCode"]
        params={
          "taxon_id":plants,"place_id":ids[iso],"d1":"2021-01-01","d2":"2021-12-31",
          "swlat":float(t["lat_floor"]),"swlng":float(t["lon_floor"]),
          "nelat":float(t["lat_ceiling"]),"nelng":float(t["lon_ceiling"]),
          "photos":"true","geo":"true","geoprivacy":"open","taxon_geoprivacy":"open",
          "quality_grade":"research","captive":"false","acc_below":500,
          "per_page":1}
        rec,j=get("/observations",params)
        total=(int(j["total_results"]) if isinstance(j,dict) and isinstance(j.get("total_results"),int) else None)
        n=int(t["public_EEA_both_measured_positive_unique_FW_outfalls"])
        row={
          **t,"2021_Plantae_photo_count":total,
          "2021_Plantae_API_receipt":rec,
          "2021_plant_effort_eligible":total is not None and total>=MIN_PLANT,
          "original_2022_source_intensity_band":band(n),
          "geographic_2degree_block":f'{math.floor(float(t["lat_floor"])/2)}_{math.floor(float(t["lon_floor"])/2)}',
          "2022_Lycaenidae_Pieridae_Papilionidae_data_queried":False,
        }
        rows.append(row)
        if i%10==0 or i==len(tiles)-1:
            print("GESUI_DE_AT_2021_INDEPENDENT_BOTANICAL_PHOTO_EFFORT",
              json.dumps({"processed":i+1,"total":len(tiles),"country":iso,
                  "reported_2022_qualified_WWTP":n,
                  "plant_2021_photos":total,"http":rec.get("http"),
                  "source_band":band(n)},sort_keys=True),flush=True)
        time.sleep(.55)
    passing=[r for r in rows if r["2021_plant_effort_eligible"]]
    nby=Counter(r["countryCode"] for r in passing)
    bby={c:sorted(set(r["original_2022_source_intensity_band"] for r in passing if r["countryCode"]==c))
         for c in ("DE","AT")}
    geo=set(r["countryCode"]+"_"+r["geographic_2degree_block"] for r in passing)
    allbands=set(r["original_2022_source_intensity_band"] for r in passing)
    eligible=len(passing)>=MIN_TOTAL and nby["DE"]>=MIN_DE and nby["AT"]>=MIN_AT and len(allbands)==3 and all(len(bby[c])>=2 for c in bby) and len(geo)>=MIN_GEO_BLOCKS
    counts_by_band={c:dict(Counter(r["original_2022_source_intensity_band"] for r in passing if r["countryCode"]==c)) for c in ("DE","AT")}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    out={
      "schema":"gesui_DE_AT_original_EEA2022_1degree_tile_2021_Plant_only_effort_freeze_v01",
      "retrieved_UTC":datetime.now(timezone.utc).isoformat(),
      "source_only_original_2022_EEA_manifest_sha256":SOURCE_SHA,
      "source_only_EEA_1degree_tile_QC_sha256":sha(a.source_tiles_qc),
      "source_only_EEA_1degree_tile_csv_sha256":sha(a.source_tiles),
      "original_preselected_source_positive_DE_AT_grid_n":len(tiles),
      "original_preselected_source_positive_by_country":dict(countby),
      "official_iNaturalist_Plantae_kingdom":tmeta,
      "official_iNaturalist_country_places":places,
      "min_2021_Plant_records_per_tile":MIN_PLANT,
      "plant_2021_eligible_DE_AT_source_tiles":len(passing),
      "plant_2021_eligible_by_country":dict(nby),
      "plant_2021_eligible_source_bands_by_country":bby,
      "plant_2021_eligible_source_band_count_by_country":counts_by_band,
      "plant_2021_eligible_distinct_country_2deg_blocks":len(geo),
      "gate_2021_effort_eligibility_before_2022_butterfly_outcome":eligible,
      "number_unavailable_source_queries":sum(r["2021_Plantae_photo_count"] is None for r in rows),
      "source_and_2021_plant_only_full_cell_records":rows,
      "no_2022_butterfly_photo_records_or_2022_butterfly_taxon_used":True,
      "no_actual_sewage_butterfly_effect":True,
      "actual_site_based_adult_puddling_not_observed":True,
      "source_only_EEA_mass_is_reported_WWTP_not_instream_dose":True,
      "inference_target":"spatial photo reporting share across original source-positive landscape cells, not causal sewage effect",
    }
    a.out.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_DE_AT_2021_PLANT_EFFORT_GATE_COMPLETE_NO_2022_BUTTERFLIES",
      json.dumps({"all_source_cells":len(tiles),"passing":len(passing),
         "countries":dict(nby),"bands":counts_by_band,"country_blocks":len(geo),
         "gate":eligible,"api_unavailable":out["number_unavailable_source_queries"]},sort_keys=True),flush=True)

if __name__=="__main__":main()
