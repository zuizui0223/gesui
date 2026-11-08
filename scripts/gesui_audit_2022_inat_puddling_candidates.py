#!/usr/bin/env python3
"""2022 European butterfly public text 'puddling' CANDIDATE micro-audit.

No photo pixels downloaded/rehosted; a free-text match, observation field,
and a photograph are NOT confirmation of puddling behavior.
Queries use validated iNaturalist family descendant taxon IDs and the
same 2022 European bbox as the source-only panel.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,urllib.parse,urllib.request
from pathlib import Path

BASE="https://api.inaturalist.org/v1/observations"
TAXA={"Lycaenidae":47923,"Pieridae":48508,"Papilionidae":47223}
BBOX={"swlat":"44","swlng":"-12","nelat":"64","nelng":"26"}
YEAR="2022"

def query(name,taxid):
    params={**BBOX,"taxon_id":str(taxid),"d1":"2022-01-01","d2":"2022-12-31",
            "photos":"true","geo":"true","geoprivacy":"open",
            "taxon_geoprivacy":"open","q":"puddling","per_page":"100"}
    u=BASE+"?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(u,headers={
        "User-Agent":"gesui-lycaenid-puddling-text-candidate-audit/0.1",
        "Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=40) as r:
        raw=r.read(800000)
        j=json.loads(raw)
        if j.get("total_results",0)>100:
            raise ValueError("Candidate count >100, would need page-scoped exhaustive enumeration")
    found=[]
    for x in j.get("results",[]):
        # Use only explicitly non-obscured open public GPS data, reject otherwise.
        geo=x.get("geojson") or {}
        coords=geo.get("coordinates") if isinstance(geo,dict) else None
        xy_public=(coords and len(coords)>=2 and
                   not x.get("coordinates_obscured") and
                   x.get("geoprivacy") not in ("private","obscured"))
        latitude=float(coords[1]) if xy_public else None
        longitude=float(coords[0]) if xy_public else None
        accuracy=x.get("positional_accuracy")
        fields=x.get("ofvs") or x.get("observation_field_values") or []
        body=" ".join([str(x.get("description") or ""),str(x.get("species_guess") or ""),
                 " ".join(str(f.get("value") or "") for f in fields)])
        photo=x.get("photos") or []
        item={
            "source_observation_url":"https://www.inaturalist.org/observations/"+str(x.get("id")),
            "observation_id":x.get("id"),
            "observed_on":x.get("observed_on"),
            "family_query":name,
            "identified_taxon":(x.get("taxon") or {}).get("name"),
            "identified_taxon_rank":(x.get("taxon") or {}).get("rank"),
            "public_wgs84_latitude":latitude,
            "public_wgs84_longitude":longitude,
            "positional_accuracy_m":accuracy,
            "public_open_precise_location_not_obscured":bool(xy_public),
            "eligible_500m_positional_accuracy":bool(xy_public and isinstance(accuracy,(int,float)) and 0<=accuracy<=500),
            "photos_n":len(photo),
            "photo_first_license":photo[0].get("license_code") if photo else None,
            "has_text_puddling_casefold":("puddling" in body.casefold()),
            "fieldIDs_and_values":[{"id":f.get("observation_field_id"),
                "value_contains_puddling":("puddling" in str(f.get("value") or "").casefold())}
                for f in fields][:12],
            "visual_behavior_confirmed":False,
            "proven_photo_2022_effluent_chemistry":False,
        }
        found.append(item)
    return {"source_url":u,"source_response_sha256":hashlib.sha256(raw).hexdigest(),
            "api_total_results":j.get("total_results"),
            "candidate_source_n":len(found),
            "eligible_500m_public_geo_n":sum(r["eligible_500m_positional_accuracy"] for r in found),
            "IDs":[r["observation_id"] for r in found],
            "candidates":found}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    result={}
    for f,id in TAXA.items():
        item=query(f,id)
        result[f]=item
        print("GESUI_2022_PUDDLING_TEXT_SOURCE_CANDIDATES",json.dumps({
            "family":f,"taxid":id,
            "candidates":item["candidate_source_n"],
            "accurate_open_geo":item["eligible_500m_public_geo_n"],
            "not_photo_confirmed":True},sort_keys=True),flush=True)
    obj={"schema":"gesui_2022_iNat_Europe_puddling_free_text_candidates_unverified_v01",
         "recorded_at":datetime.datetime.now(datetime.timezone.utc).isoformat(),
         "year":2022,"bbox":BBOX,"candidates_by_family":result,
         "visual_confirmation_is_REQUIRED_for_puddling":True,
         "all_obs_pending_visual_review":True,
         "source_EEA_WWTP_and_nutrient_exposure_not_linked":True,
         "no_actual_puddling_rate_or_wwtp_effect_estimated":True,
         "limitations":[
            "Query hits are search-engine responses, not behavior-confirmed photos",
            "Obscured or inaccurate locations excluded from exact WWTP proximity",
            "Photo license and original observer credits govern any reproduction",
            "Observer effort, seasonal detectability, and field presence not controlled",
            "No zero-absence inference from untagged photos",
         ]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    print("GESUI_2022_TEXT_PUDDLING_CANDIDATES_ENUMERATION_SUCCESS",flush=True)
if __name__=="__main__":
    main()
