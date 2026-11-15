#!/usr/bin/env python3
"""Source-authenticated original iNat metadata for THREE possible carrion field500 events.

No image pixel content downloaded or redistributed; no manual photo viewing
or proboscis claim. 2015 source pair same observer/date is not counted as two
independent carcasses. Confirm exact photo licences before any image access.
"""
from __future__ import annotations
import argparse, hashlib, json, time, urllib.parse, urllib.request, collections
from pathlib import Path

CASES=[150764924,150764928,367297297]
ORIGINAL_SOURCE_GIT_SHA="a44fb5b67f3b27b1243b8f1db66a595d4347a2d7"
API="https://api.inaturalist.org/v2/observations"
FIELDS="id,description,observed_on,user.id,taxon.id,taxon.name,taxon.rank,photos.id,photos.url,photos.license_code,photos.attribution,photos.original_dimensions,ofvs.name,ofvs.value,quality_grade,geojson.coordinates"
def blob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def get():
    params={"id":",".join(map(str,CASES)),"per_page":10,"fields":FIELDS}
    url=API+"?"+urllib.parse.urlencode(params)
    e=None
    for i in range(4):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json",
               "User-Agent":"GESUI-field500-carrion-source-metadata-no-photos/1.0"})
            with urllib.request.urlopen(req,timeout=50) as q:
                raw=q.read(2_000_000)
            d=json.loads(raw)
            if not isinstance(d.get("results"),list):raise RuntimeError("Unexpected iNat v2 observations")
            return d,hashlib.sha256(raw).hexdigest()
        except Exception as err:
            e=str(err)
            if i<3:time.sleep(i+1)
    raise RuntimeError("Original 3 case source metadata unreachable: "+e)
def run(a):
    raw=a.baseline.read_bytes()
    if blob(raw)!=ORIGINAL_SOURCE_GIT_SHA:raise ValueError("Original frozen three-case source manifest SHA drift")
    original=json.loads(raw)
    expected=original["carrion_candidate_rows"]
    if len(expected)!=3:raise ValueError("Original three cases changed")
    found,sha=get()
    cases=found["results"]
    if set(x.get("id") for x in cases)!=set(CASES):raise ValueError("Source three original observed IDs not recovered")
    cases=sorted(cases,key=lambda x:CASES.index(x["id"]))
    records=[]; owners=[]; dates=[]
    for x in cases:
        idx=CASES.index(x["id"])
        ofvs=x.get("ofvs") or []
        vals=[str(v.get("value") or "").strip() for v in ofvs if v.get("name")=="Butterfly Puddling Substrate"]
        if len(vals)!=1 or vals[0].casefold()!=expected[idx]["original_value"].casefold():
            raise ValueError("Original source field500 text does not match prior research source")
        user=(x.get("user") or {}).get("id")
        date=x.get("observed_on")
        photos=x.get("photos") or []
        photo_ids=[p.get("id") for p in photos]
        original_photo_ids=expected[idx]["photo_count_metadata"]
        if len(photo_ids)!=original_photo_ids:raise ValueError("Original source photo count changed for case")
        desc=(x.get("description") or "").strip()
        keywords=["feed","proboscis","puddl","liquid","dead","carcass","amphibian","animal","corpse"]
        matched=[term for term in keywords if term in desc.casefold()]
        records.append({"original_observation_url":expected[idx]["original_observation_url"],
            "original_field500_value":vals[0],
            "reported_observation_date":date,
            "source_taxon_name":(x.get("taxon") or {}).get("name"),
            "source_taxon_rank":(x.get("taxon") or {}).get("rank"),
            "original_photo_ids":photo_ids,
            "original_photo_licence_codes":[p.get("license_code") for p in photos],
            "original_photo_dimensions":[p.get("original_dimensions") for p in photos],
            "source_description_available":bool(desc),
            "description_keywords_present":matched,
            "source_grade":x.get("quality_grade"),
            "photo_original_pixels_read":False,
            "substrate_visual_confirmed":False,
            "proboscis_contact_visual_confirmed":False,
            "no_wastewater_source_exposure_tested":True})
        owners.append(user);dates.append(date)
    all_photo_ids=[q for z in records for q in z["original_photo_ids"] if isinstance(q,int)]
    duplicate_photos=[x for x,v in collections.Counter(all_photo_ids).items() if v>1]
    owner_date={(owners[i],dates[i]) for i in range(3)}
    if len(owner_date)!=2:raise ValueError("Original independently audited 2 owner-days changed")
    if not (owners[0]==owners[1] and dates[0]==dates[1]):
        raise ValueError("2015 records changed owner-date identity")
    results={"schema":"GESUI_INAT_3_CARRION_ANNOTATION_ORIGINAL_PHOTO_RIGHTS_DESCRIPTION_SOURCE_GATE_V1",
      "baseline_source_git_blob":ORIGINAL_SOURCE_GIT_SHA,
      "original_source_api_SHA256":sha,
      "original_source_endpoint":API,
      "original_records":records,
      "n_original_carrion_field_text_annotations":3,
      "n_distinct_observer_days":len(owner_date),
      "n_distinct_observation_owners":len(set(owners)),
      "n_independent_carcasses_verified":None,
      "original_photo_ids_reused_among_three":duplicate_photos,
      "no_photographic_pixels_downloaded_or_viewed":True,
      "no_image_recognition_or_manual_proboscis_assessment_performed":True,
      "description_not_evidence_of_actual_ingestion":True,
      "no_original_HNS_or_WWTP_or_nutrient_exposure_tested":True,
      "next_step_photo_permission_and_controlled_visual_review_required":True,
      "source_only_stage":"META_RIGHTS_AND_REPEAT_CASE_REVIEW_NOT_ADULT_FEEDING_CONFIRMATION"}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(results,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_INAT_CARRION_THREE_META_QA="+json.dumps({
      "n_records":3,"owner_day_units":len(owner_date),
      "n_photo_ids":len(all_photo_ids),
      "photo_license_codes":sorted({str(v) for x in records for v in x["original_photo_licence_codes"]}),
      "description_available":sum(x["source_description_available"] for x in records),
      "photo_pixels_seen":False},sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--baseline",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    run(p.parse_args())
