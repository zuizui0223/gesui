#!/usr/bin/env python3
"""Original iNaturalist field500 746 full-census => dung-only 46 independent source units.

Read original annotation field, owner/day, photo IDs, source photo license and
descriptions. Does NOT view or redistribute image pixels; NO actual proboscis
contact or effluent feeding claims. Original user-generated field's non-use
is not a recorded absence. Keep original 2022 HNS protocol unchanged.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,re,time,urllib.parse,urllib.request
from pathlib import Path
from audit_inat_v2_observer_photo_opportunity import category

BASELINE_GIT_SHA="d23eeb6e1b2fdf9b285e8af6ad5391b366139aa9"
FIVE_FAMILIES=("Lycaenidae","Nymphalidae","Papilionidae","Pieridae","Hesperiidae")
API="https://api.inaturalist.org/v2/observations"
FIELDS="id,user.id,observed_on,taxon.id,taxon.ancestor_ids,taxon.name,taxon.rank,ofvs.name,ofvs.value,photos.id,photos.license_code,description,quality_grade"
PER=100
LABEL="Butterfly Puddling Substrate"
def git_blob(raw):
    return hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
def get(p):
    query={"field:"+LABEL:"","per_page":PER,"page":p,"fields":FIELDS,
           "order_by":"id","order":"asc"}
    req=urllib.request.Request(API+"?"+urllib.parse.urlencode(query),headers={"User-Agent":"gesui-dung46-source-independence-20261115","Accept":"application/json"})
    last=None
    for k in range(4):
        try:
            with urllib.request.urlopen(req,timeout=55) as r:
                body=r.read(5_000_000)
            if len(body)>=5_000_000:raise RuntimeError("iNat sparse source response too large")
            data=json.loads(body)
            if not isinstance(data,dict) or not isinstance(data.get("results"),list):raise RuntimeError("Unexpected API object")
            return data,hashlib.sha256(body).hexdigest()
        except Exception as e:
            last=e
            if k<3:time.sleep(k+1)
    raise RuntimeError("Original source retrieval unavailable "+str(last))
def family(t,fids):
    t=t or {}
    ancestors=set(t.get("ancestor_ids") or [])
    if isinstance(t.get("id"),int):ancestors.add(t["id"])
    matches=[f for f,id in fids.items() if id in ancestors]
    return matches[0] if len(matches)==1 else "OTHER_OR_UNCLASSIFIED"
def main(args):
    baseline=args.baseline.read_bytes()
    if git_blob(baseline)!=BASELINE_GIT_SHA:raise ValueError("original full-746 v2 source SHA changed")
    old=json.loads(baseline)
    families={f:x["original_id"] for f,x in json.loads(args.old_census.read_text())["butterfly_family_taxa_from_iNat_authority"].items()}
    if old["n_original_field500_observations"]!=746 or not old["original_field500_family_census_exactly_reproduced"]:
        raise ValueError("Full original source not proven")
    observed={}; receipts=[]; n=None
    for page in range(1,15):
        if page>1:time.sleep(1)
        data,sha=get(page)
        if n is None:n=int(data.get("total_results",0))
        if n!=746:raise ValueError("Original 746 count drift")
        if not data["results"]:break
        for raw in data["results"]:
            i=raw.get("id")
            if not isinstance(i,int) or i in observed:raise RuntimeError("duplicated API observation")
            values=[x.get("value") for x in raw.get("ofvs") or [] if x.get("name")==LABEL]
            if len(values)!=1:raise RuntimeError("Target field500 missing or duplicated")
            name=str(values[0]).strip()
            photos=raw.get("photos") or []
            obs=raw.get("observed_on") or ""
            observed[i]={"id":i,"family":family(raw.get("taxon"),families),"substrate_label":name,
                "category":category(name),"owner":(raw.get("user") or {}).get("id"),
                "date":obs[:10] if len(obs)>=10 else None,
                "year":obs[:4] if len(obs)>=4 else None,
                "photos":[x.get("id") for x in photos if isinstance(x.get("id"),int)],
                "licences":[x.get("license_code") for x in photos],
                "description":str(raw.get("description") or ""),
                "quality_grade":raw.get("quality_grade")}
        receipts.append({"page":page,"source_sha256":sha,"n":len(data["results"])})
        if len(observed)>=n:break
    if len(observed)!=746:raise ValueError("Incomplete original iNat source census")
    dung=[v for v in observed.values() if v["category"]=="dung_feces"]
    mix=[v for v in observed.values() if v["category"].startswith("mixed_") and "dung_feces" in v["category"]]
    if len(dung)!=46 or len(mix)!=3:raise ValueError("Original revised dung46 + mixed3 changed")
    byfam=collections.Counter(x["family"] for x in dung)
    expected={"Nymphalidae":24,"Lycaenidae":9,"Papilionidae":6,"Hesperiidae":5,"OTHER_OR_UNCLASSIFIED":2}
    if dict(byfam)!=expected:raise ValueError("Original 46 family score/source records changed "+str(byfam))
    by_owner=collections.defaultdict(list)
    for v in dung:by_owner[v["owner"]].append(v)
    photo_obs=collections.defaultdict(list)
    for v in dung:
        for photo_id in v["photos"]:photo_obs[photo_id].append(v["id"])
    dup_photos={k:v for k,v in photo_obs.items() if len(set(v))>1}
    owner_days=collections.defaultdict(list)
    for v in dung:owner_days[(v["owner"],v["date"])].append(v)
    owner_day_fam=collections.defaultdict(list)
    for v in dung:owner_day_fam[(v["owner"],v["date"],v["family"])].append(v)
    byfamily={}
    for fam in FIVE_FAMILIES:
        rows=[x for x in dung if x["family"]==fam]
        owners={x["owner"] for x in rows}
        days={(x["owner"],x["date"]) for x in rows}
        ph=set(p for x in rows for p in x["photos"])
        byfamily[fam]={"annotations":len(rows),"owners":len(owners),
            "owner_days":len(days),"photos_distinct":len(ph),
            "n_owners_with_two_or_more_tagged_dung_obs":sum(sum(x["owner"]==o for x in rows)>=2 for o in owners)}
    cotag={}
    for fam in ("Lycaenidae","Nymphalidae"):
        fam_users=set(x["owner"] for x in dung if x["family"]==fam)
        cotag[fam]=list(fam_users) # transient only
    common_users=set(cotag["Lycaenidae"])&set(cotag["Nymphalidae"])
    photo_licenses=collections.Counter(str(y) for x in dung for y in x["licences"])
    known_keywords=("proboscis","feeding","feed","drinking","puddling","sucking","sipping")
    desc=[v for v in dung if v["description"]]
    desc_keyword=[v for v in desc if any(q in v["description"].lower() for q in known_keywords)]
    # Raw example source URLs only; no images copied or owner identifiers published.
    review=[]
    for r in sorted(dung,key=lambda x:(x["family"],x["id"])):
        review.append({"original_observation_url":"https://www.inaturalist.org/observations/"+str(r["id"]),
            "family":r["family"],"raw_substrate_label":r["substrate_label"],
            "observed_year":r["year"],
            "photo_count":len(r["photos"]),
            "photo_licence_codes":"|".join(sorted(set(str(v) for v in r["licences"]))),
            "description_mentions_feeding_keyword":any(k in r["description"].lower() for k in known_keywords),
            "same_owner_date_dung_annotations":len(owner_days[r["owner"],r["date"]]),
            "same_owner_date_same_family_dung_annotations":len(owner_day_fam[r["owner"],r["date"],r["family"]]),
            "original_photo_pixels_viewed":"NOT_VIEWED",
            "original_substrate_visually_verified":"NOT_REVIEWED",
            "proboscis_contact_visually_verified":"NOT_REVIEWED",
            "two_independent_reviewers_agree":"NOT_REVIEWED",
            "ecological_feeding_event_confirmed":"NOT_REVIEWED"})
    owner_counts=sorted((len(x) for x in by_owner.values()),reverse=True)
    n_eff=len(dung)**2/sum(x*x for x in owner_counts)
    n_multicluster=sum(len(v)>1 for v in owner_days.values())
    q={
      "schema":"GESUI_I_NAT_DUNG46_COMPLETE_746_SOURCE_INDEPENDENCE_NOT_FEEDING_V1",
      "original_v2_746_source_git_blob":BASELINE_GIT_SHA,
      "source_original_iNat_v2_SHA256_pages":receipts,
      "n_all_original_source_field500_records":746,
      "source_dung_text_only_record_count":len(dung),
      "source_mixed_dung_and_mud_text_record_count":len(mix),
      "dung_original_substrate_classification_source_maintained":True,
      "family":byfamily,
      "total_distinct_dung_observation_owners":len(by_owner),
      "n_distinct_dung_owner_days":len(owner_days),
      "n_owner_day_clusters_with_multiple_dung_annotations":n_multicluster,
      "n_dung_records_in_multiple_obs_owner_day_clusters":sum(len(x) for x in owner_days.values() if len(x)>1),
      "dung_inverse_simpson_effective_owners":round(n_eff,4),
      "dung_top_owner_fraction":round(owner_counts[0]/len(dung),4),
      "dung_top5_owner_fraction":round(sum(owner_counts[:5])/len(dung),4),
      "n_dung_annotation_owners_in_both_Lycaenidae_and_Nymphalidae":len(common_users),
      "dung_count_photo_license_codes":dict(photo_licenses),
      "distinct_dung_photo_id_count":len(photo_obs),
      "n_reused_original_photo_ids_across_dung_annotations":len(dup_photos),
      "n_dung_observations_with_reused_photo_ids":len(set(y for arr in dup_photos.values() for y in arr)),
      "n_dung_original_records_with_description":len(desc),
      "n_dung_description_mentions_feeding_keyword":len(desc_keyword),
      "source_descriptions_not_direct_visual_proboscis_confirmation":True,
      "no_photo_pixels_viewed_or_downloaded":True,
      "no_causal_sewage_or_NaN_nutrient_measured":True,
      "no_representative_observation_opportunity_denominator":True,
      "no_independent_proboscis_contact_event_confirmed":True,
      "observed_family_dung_ANNOTATION_COUNTS_NOT_FEEDING_RATES":True,
      "source_decision":"PHOTO_REVIEW_REQUIRED_NONINDEPENDENT_USER_GENERATED_DUNG_TAGGED_OBSERVATIONS"
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(q,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
    args.review.parent.mkdir(parents=True,exist_ok=True)
    with args.review.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(review[0]));writer.writeheader();writer.writerows(review)
    print("GESUI_DUNG46_ORIGINAL_SOURCE_QA="+json.dumps({
      "records":q["source_dung_text_only_record_count"],"dung_owner_days":q["n_distinct_dung_owner_days"],
      "owners":q["total_distinct_dung_observation_owners"],
      "n_eff_owners":q["dung_inverse_simpson_effective_owners"],
      "byfamily":byfamily,"coannotating_owners_Lyca_Nymph":len(common_users),
      "duplicate_photo_ids":len(dup_photos),"descriptions":len(desc),
      "photo_pixels_seen":False},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--baseline",type=Path,required=True)
    p.add_argument("--old-census",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--review",type=Path,required=True)
    main(p.parse_args())
