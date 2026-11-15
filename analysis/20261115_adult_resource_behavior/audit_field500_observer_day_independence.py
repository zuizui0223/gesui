#!/usr/bin/env python3
"""Independent event/photograph source audit for iNaturalist field 500.

Answers: are 3 animal-remains labels independent original observer/date/photo events?
How much original substrate annotation data overlap across observers and two
butterfly families? Does source include repeated identical photo IDs?

This file DOES NOT download or inspect actual photograph pixels. No feeding,
proboscis contact, carcass identity, wastewater exposure, habitat chemistry,
abundance or source-specific biological causal effects are estimated.
"""
from __future__ import annotations
import argparse, collections, csv, hashlib, io, json, pathlib, time, urllib.parse, urllib.request
from audit_inat_v2_observer_photo_opportunity import (category, organic_bucket)

ORIGINAL_JSON_BLOB="d23eeb6e1b2fdf9b285e8af6ad5391b366139aa9"
ORIGINAL_REVIEW_BLOB="fb1be650a70283513203f96ebd260ab2619a7d92"
URL="https://api.inaturalist.org/v2/observations"
FIELDS="id,user.id,taxon.id,taxon.ancestor_ids,ofvs.name,ofvs.value,photos.id,observed_on"
FIELD="Butterfly Puddling Substrate"
N=746
PAGE_SIZE=100

def git_blob(x):
    return hashlib.sha1(b"blob "+str(len(x)).encode()+b"\0"+x).hexdigest()
def fetch(page):
    params={"field:"+FIELD:"","per_page":PAGE_SIZE,"page":page,"fields":FIELDS,
            "order_by":"id","order":"asc"}
    url=URL+"?"+urllib.parse.urlencode(params)
    err=None
    for i in range(3):
        try:
            req=urllib.request.Request(url,headers={
              "User-Agent":"gesui-20261115-puddling-observer-event-duplicated-photos-QA/1.0",
              "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=60) as r:
                raw=r.read(5_000_000)
            if len(raw)>=5_000_000:raise RuntimeError("iNat sparse v2 response size exceeded")
            j=json.loads(raw)
            if not isinstance(j.get("results"),list):raise ValueError("Invalid iNat response")
            return j,hashlib.sha256(raw).hexdigest()
        except Exception as e:
            err=e
            if i<2:time.sleep(1+i)
    raise RuntimeError("Failed exact original iNat source retrieval "+str(err))

def normdate(s):
    return (s or "")[:10] if (s or "").count("-")>=2 else None
def main(a):
    old=a.baseline.read_bytes()
    review=a.review.read_bytes()
    if git_blob(old)!=ORIGINAL_JSON_BLOB or git_blob(review)!=ORIGINAL_REVIEW_BLOB:
        raise ValueError("Original frozen independent 746 source or 78 reference photo audit changed")
    original=json.loads(old)["fields"]["substrate500"]
    srcqueue=list(csv.DictReader(io.StringIO(review.decode("utf-8-sig"))))
    if len(srcqueue)!=78 or original["source_total_results"]!=N or not original["fetch_complete"]:
        raise ValueError("Source observation denominators changed")
    fams={k:v["original_id"] for k,v in json.loads(old)["butterfly_family_taxa_from_iNat_authority"].items()}
    rows=[]
    log=[]
    source_total=None
    for page in range(1,12):
        if page>1:time.sleep(1.10)
        data,sha=fetch(page)
        if source_total is None:source_total=int(data.get("total_results",0))
        if source_total!=N:raise ValueError("Original field500 746 source population changed")
        incoming=data["results"]
        if not incoming:break
        log.append({"page":page,"items":len(incoming),"source_response_SHA256":sha})
        for x in incoming:
            v=[z.get("value","") for z in (x.get("ofvs") or []) if z.get("name")==FIELD]
            if len(v)!=1:raise ValueError("Ambiguous missing original iNat field500")
            tx=x.get("taxon") or {}
            taxids=set(tx.get("ancestor_ids") or [])
            if isinstance(tx.get("id"),int):taxids.add(tx["id"])
            names=[f for f,k in fams.items() if k in taxids]
            family=names[0] if len(names)==1 else "OTHER_OR_UNCLASSIFIED"
            value=str(v[0]).strip()
            img=[z.get("id") for z in x.get("photos") or [] if isinstance(z,dict)]
            rows.append({"obs_id":x["id"],"date":normdate(x.get("observed_on")),
                "user":(x.get("user") or {}).get("id"),
                "family":family,"raw":value,"category":category(value),
                "organic":organic_bucket(category(value)),
                "photo_ids":set(z for z in img if isinstance(z,int))})
        if len(rows)>=N:break
        if len(incoming)<PAGE_SIZE:break
    if len(rows)!=N or len({x["obs_id"] for x in rows})!=N:raise ValueError("Missing/duplicated 746 original obs")
    if dict(collections.Counter(x["family"] for x in rows))!=original["family_ancestor_counts_among_selected_tagged_observations"]:
        raise ValueError("Original iNat family source-shape mismatch")
    if set(x["obs_id"] for x in rows if x["organic"] is not None)!={int(z["original_observation_url"].split("/")[-1]) for z in srcqueue}:
        raise ValueError("Original review queue 78 organic ids fail exact validation")
    photomap=collections.defaultdict(set)
    for x in rows:
        for p in x["photo_ids"]:photomap[p].add(x["obs_id"])
    dupimages={id:vv for id,vv in photomap.items() if len(vv)>1}
    dupobs=set().union(*dupimages.values()) if dupimages else set()
    obsbyuser=collections.defaultdict(list)
    for x in rows:obsbyuser[x["user"]].append(x)
    allusers=[len(v) for v in obsbyuser.values()]
    ess=len(rows)**2/sum(n*n for n in allusers)
    carb=[x for x in rows if x["category"]=="animal_carrion_candidate"]
    if len(carb)!=3:raise ValueError("Carrion source candidate count changed")
    daygroups=collections.Counter((x["user"],x["date"]) for x in carb)
    all_daygroups=collections.Counter((x["user"],x["date"]) for x in rows)
    two_fam=[]
    for uid,li in obsbyuser.items():
        ll=[x for x in li if x["family"]=="Lycaenidae"]
        nn=[x for x in li if x["family"]=="Nymphalidae"]
        if not (ll and nn):continue
        la=sum(x["organic"] is not None for x in ll)/len(ll)
        na=sum(x["organic"] is not None for x in nn)/len(nn)
        two_fam.append({"ly_n":len(ll),"nym_n":len(nn),
            "ly_organic_n":sum(x["organic"] is not None for x in ll),
            "nym_organic_n":sum(x["organic"] is not None for x in nn),
            "difference_ly_minus_nym":la-na,
            "informative_has_mineral_and_organic":any(x["organic"] is not None for x in li)
                 and any(x["category"] in ("mud_soil_sand","rock_gravel") for x in li)})
    if len(two_fam)!=18:raise ValueError("Original 18 both-family observer gate changed")
    diffs=[x["difference_ly_minus_nym"] for x in two_fam]
    mean=sum(diffs)/len(diffs)
    loo=[(sum(diffs)-x)/(len(diffs)-1) for x in diffs]
    carrion_samples=[]
    for x in carb:
        others=[z for z in carb if z["obs_id"]!=x["obs_id"]]
        shared_photo=sum(bool(x["photo_ids"]&z["photo_ids"]) for z in others)
        shared_owner_day=sum(x["user"]==z["user"] and x["date"] and x["date"]==z["date"] for z in others)
        linked=[z for z in srcqueue if z["original_observation_url"].endswith("/"+str(x["obs_id"]))]
        if len(linked)!=1:raise ValueError("One original animal candidate not in review queue")
        carrion_samples.append({"original_observation_url":linked[0]["original_observation_url"],
            "original_value":x["raw"],"photo_count_metadata":len(x["photo_ids"]),
            "year":(x["date"] or "")[:4],
            "shares_identical_photos_with_other_two":bool(shared_photo),
            "same_observer_same_day_as_other_two":bool(shared_owner_day),
            "photo_pixel_review_status":"NOT_INSPECTED",
            "proboscis_contact_status":"NOT_INSPECTED"})
    result={"schema":"GESUI_20261115_FIELD500_746_ORIGINAL_PHOTO_AND_OBSERVERDAY_NONINDEPENDENCE_V1",
            "baseline_746_original_JSON_git_blob":ORIGINAL_JSON_BLOB,
            "original_78_photo_queue_git_blob":ORIGINAL_REVIEW_BLOB,
            "sparse_original_API":URL,"n_full_original_field500_observations":len(rows),
            "original_HTTP_response_SHA256_by_page":log,
            "source_unique_observation_owners":len(obsbyuser),
            "inverse_simpson_effective_observation_owners":round(ess,3),
            "source_owner_n_records_top1":max(allusers),
            "n_photo_ids_including_secondary_photos":len(photomap),
            "n_shared_photo_ids_across_distinct_original_observations":len(dupimages),
            "n_observations_with_at_least_one_reused_photo_id":len(dupobs),
            "n_organic_original_tagged_observations":sum(x["organic"] is not None for x in rows),
            "organic_annotated_distinct_observer_days":len({(x["user"],x["date"]) for x in rows if x["organic"] is not None}),
            "all_tagged_distinct_observer_days":len(all_daygroups),
            "carrion_candidate_annotated_observations":len(carb),
            "carrion_candidate_distinct_observers":len({x["user"] for x in carb}),
            "carrion_candidate_distinct_observer_days":len(daygroups),
            "carrion_candidate_rows":carrion_samples,
            "n_original_owners_with_both_Lyca_Nymph":len(two_fam),
            "n_original_owners_with_both_families_and_organic_mineral":sum(x["informative_has_mineral_and_organic"] for x in two_fam),
            "observer_equal_mean_Ly_minus_Nymph_organic_ANNOTATION_share":round(mean,6),
            "observer_equal_leave_one_owner_out_range":[round(min(loo),6),round(max(loo),6)],
            "n_both_fam_observers_positive_ly_minus_nymph":sum(x>0 for x in diffs),
            "n_both_fam_observers_equal_ly_nymph":sum(x==0 for x in diffs),
            "n_both_fam_observers_negative_ly_minus_nymph":sum(x<0 for x in diffs),
            "n_owners_with_any_organic_field_record":sum(any(x["organic"] is not None for x in li) for li in obsbyuser.values()),
            "observation_day_clusters_without_location_not_independent_resource_events":True,
            "photo_ID_reuse_proves_observation_images_reused_not_same_individual":True,
            "all_original_photo_pixels_NOT_DOWNLOADED_OR_VIEWED":True,
            "original_photo_link_healthy_not_tested":True,
            "original_carrion_from_text_not_visually_verified":True,
            "owner_stratified_label_association_NOT_preference_or_wastewater_effect":True,
            "no_personal_id_or_coordinates_in_output":True,
            "no_original_HNS_WWTP_or_nutrient_exposure_assessed":True,
            "status":"EXPLORATORY_SOURCE_DATA_INDEPENDENCE_AND_OBSERVER_QA_ONLY"}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2)+"\n",encoding="utf8")
    print("GESUI_INAT_746_DUPLICATE_AND_OWNERDAY_QA="+json.dumps({
      "carrion_n":len(carb),"carrion_independent_observer_days":len(daygroups),
      "carrion_independent_observers":result["carrion_candidate_distinct_observers"],
      "photo_reused_IDs":len(dupimages),"photo_reused_observations":len(dupobs),
      "n_eff_owners":round(ess,2),"both_families_owners":len(two_fam),
      "observer_equal_annotation_difference":round(mean,3),"observer_LOO":result["observer_equal_leave_one_owner_out_range"],
      "source_only":True},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--baseline",type=pathlib.Path,required=True)
    p.add_argument("--review",type=pathlib.Path,required=True)
    p.add_argument("--output",type=pathlib.Path,required=True)
    main(p.parse_args())
