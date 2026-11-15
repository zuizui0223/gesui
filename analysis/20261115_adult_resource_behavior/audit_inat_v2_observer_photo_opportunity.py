#!/usr/bin/env python3
"""Source-audited iNaturalist v2 sparse complete field-500 records.

One different question from locked wastewater HNS: characterize annotations,
observer dominance, and photo audit eligibility. DO NOT call this documented
substrate use a verified proboscis nutrient uptake, selection coefficient, or
source-linked wastewater effect. Photo links refer to iNat originals, not
redistributed or downloaded photo bytes.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,re,time,urllib.parse,urllib.request
from pathlib import Path
from audit_inat_puddling_substrate_original_obs import substrate as old_substrate

BASE="https://api.inaturalist.org/v2/observations"
FIELDS="id,user.id,taxon.id,taxon.ancestor_ids,taxon.name,taxon.rank,ofvs.id,ofvs.name,ofvs.value,photos.url,observed_on,quality_grade"
QUERY_KEY="field:Butterfly Puddling Substrate"
PUB_FIELD="Butterfly Puddling Substrate"
BASELINE_GIT_BLOB="6fbc3358dc21b466afd2edbcd6a1182bb3797bbf"
EXPECTED_OBS=746
BATCH=100
CAP_BYTES=4_000_000
ORGANIC_CLASSES={"dung_feces","animal_carrion_candidate","plant_fruit_sap"}
SOURCE_EXPOSURE_FIELDS={"latitude","longitude","geojson","location","private_location","place_guess"}
def gitsha(x):
    return hashlib.sha1(b"blob "+str(len(x)).encode()+b"\0"+x).hexdigest()
def get_page(page):
    q={QUERY_KEY:"","per_page":BATCH,"page":page,"fields":FIELDS,
       "order":"asc","order_by":"id"}
    url=BASE+"?"+urllib.parse.urlencode(q)
    last=None
    for k in range(3):
        try:
            req=urllib.request.Request(url,headers={
              "User-Agent":"gesui-field500-observer-photo-source-audit/1.0",
              "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=65) as r:
                raw=r.read(CAP_BYTES+1)
            if len(raw)>CAP_BYTES:raise RuntimeError("Source response exceeded 4MB")
            d=json.loads(raw)
            if not isinstance(d,dict) or not isinstance(d.get("results"),list):
                raise RuntimeError("iNat v2 response not original observations")
            return d,hashlib.sha256(raw).hexdigest()
        except Exception as e:
            last=e
            if k<2:time.sleep(k+1)
    raise RuntimeError("Original v2 field search failed "+str(last))
def family_from(t,family_ids):
    if not isinstance(t,dict):return "OTHER_OR_UNCLASSIFIED"
    ids=set(t.get("ancestor_ids") or [])
    if t.get("id"):ids.add(t["id"])
    v=[f for f,k in family_ids.items() if k in ids]
    return v[0] if len(v)==1 else "OTHER_OR_UNCLASSIFIED"
def category(s):
    s=re.sub(r"\s+"," ",s.casefold().strip())
    # Revised assignment preserves the important organic substrate separation.
    if not s:return "unknown"
    hits=[]
    animal=r"carrion|carcass|corpse|dead\s+(?:fish|animal|bird|frog|body|creature)|animal remains|fish remains|(?:dead|decomposing)\s+organism"
    dung=r"faec|fec|dung|poop|droppings|manure|scat|excrement"
    fruit=r"\b(?:fruit|fruits|banana|sap)\b|rotting\s+(?:fruit|banana)|decaying vegetation"
    salt=r"urine|sweat|perspir|salt(?:ed)?\s+(?:water|soil)|uric"
    mud=r"mud|soil|sand|clay|earth|silt|\bdirt\b|bare ground"
    rock=r"\brock|gravel|stone|pebble"
    water=r"^(?:water|puddle|river|stream|wet(?:\s+ground)?|rainwater)$"
    for name,pat in [("animal_carrion_candidate",animal),("dung_feces",dung),
                     ("plant_fruit_sap",fruit),("urine_sweat_salt",salt),
                     ("mud_soil_sand",mud),("rock_gravel",rock),("water_only",water)]:
        if re.search(pat,s):hits.append(name)
    if len(hits)>1:return "mixed_"+("+".join(hits))
    if hits:return hits[0]
    return "unclassified"
def organic_bucket(cat):
    if cat=="animal_carrion_candidate":return "ANIMAL_REMAINS"
    if cat=="dung_feces":return "DUNG"
    if cat=="plant_fruit_sap":return "PLANT_FRUIT_SAP"
    if cat.startswith("mixed_") and any(x in cat for x in ("animal_carrion","dung_feces","plant_fruit")):
        return "MIXED_WITH_ORGANIC_CANDIDATE"
    return None
def main(args):
    data=args.baseline.read_bytes()
    if gitsha(data)!=BASELINE_GIT_BLOB:
        raise ValueError("Original 746 annotation census SHA changed")
    original=json.loads(data)
    first=original["fields"]["substrate500"]
    if not first["fetch_complete"] or first["fetched_distinct_observations"]!=EXPECTED_OBS:
        raise ValueError("Original v1 complete 746 census no longer locked")
    family_ids={k:v["original_id"] for k,v in original["butterfly_family_taxa_from_iNat_authority"].items()}
    obs={}
    pages=[]
    total=None
    for page in range(1,15):
        if page>1:time.sleep(1.10)
        d,sha=get_page(page)
        if total is None:total=int(d.get("total_results",0))
        if total!=EXPECTED_OBS:raise ValueError("Original annotated source population changed")
        records=d["results"]
        if len(records)==0:break
        accepted=0
        for v in records:
            id=v.get("id")
            if not isinstance(id,int):raise RuntimeError("Observed id not integer")
            if id in obs:raise RuntimeError("Original field pagination duplicated obs")
            values=[str(q.get("value") or "").strip() for q in v.get("ofvs") or []
                if q.get("name")==PUB_FIELD]
            if len(values)!=1:raise RuntimeError("Field 500 not uniquely present in source ID "+str(id))
            user=v.get("user") or {}
            taxon=v.get("taxon") or {}
            photos=v.get("photos") or []
            photo_urls=[p.get("url") for p in photos if isinstance(p,dict) and isinstance(p.get("url"),str)]
            year=str(v.get("observed_on") or "")[:4]
            obs[id]={"id":id,"user":user.get("id"),"family":family_from(taxon,family_ids),
                     "taxon_name":taxon.get("name"),"value":values[0],
                     "old_category":old_substrate(values[0]),
                     "new_category":category(values[0]),
                     "year":year,"has_photo":bool(photo_urls),
                     "n_photos":len(photo_urls),"photo_url":photo_urls[0] if photo_urls else "",
                     "quality_grade":v.get("quality_grade")}
            accepted+=1
        pages.append({"page":page,"source_page_SHA256":sha,"returned":len(records),
                      "matching_valid_field500":accepted})
        if len(obs)>=total:break
        if len(records)<BATCH:break
    if total!=EXPECTED_OBS or len(obs)!=EXPECTED_OBS:
        raise RuntimeError("Full original source not recovered: "+str(len(obs))+"/"+str(total))
    dat=list(obs.values())
    byfamily=collections.Counter(x["family"] for x in dat)
    byold=collections.Counter(x["old_category"] for x in dat)
    if dict(byfamily)!=first["family_ancestor_counts_among_selected_tagged_observations"]:
        raise RuntimeError("Original field500 family counts differ from v1 census: "+str(byfamily))
    if dict(byold)!=first["substrate_text_class_counts_if_field500"]:
        raise RuntimeError("Previous original substrate regex results do not reproduce")

    bynew=collections.Counter(x["new_category"] for x in dat)
    byfamnew=collections.defaultdict(collections.Counter)
    for x in dat:byfamnew[x["family"]][x["new_category"]]+=1
    observers=collections.defaultdict(list)
    for x in dat:
        observers[x["user"]].append(x)
    obs_sorted=sorted(observers.items(),key=lambda pair:len(pair[1]),reverse=True)
    top1_id=obs_sorted[0][0]
    top5_ids={k for k,_ in obs_sorted[:5]}
    def share(rows,focus):
        t=[r for r in rows if r["family"]==focus]
        organic=[r for r in t if organic_bucket(r["new_category"]) is not None]
        return {"annotated_n":len(t),"organic_candidate_n":len(organic),
                "organic_candidate_share_among_tagged":len(organic)/len(t) if t else None}
    samples={}
    for label,rows in [("all",dat),("exclude_top1",[x for x in dat if x["user"]!=top1_id]),
                       ("exclude_top5",[x for x in dat if x["user"] not in top5_ids])]:
        samples[label]={f:share(rows,f) for f in ("Lycaenidae","Nymphalidae","Papilionidae","Pieridae","Hesperiidae")}
    def observers_w_both(f,g):
        return [(k,li) for k,li in observers.items()
                 if any(x["family"]==f for x in li) and any(x["family"]==g for x in li)]
    both=observers_w_both("Lycaenidae","Nymphalidae")
    informative=[(k,li) for k,li in both if
                any(organic_bucket(x["new_category"]) is not None for x in li)
                and any(x["new_category"] in ("mud_soil_sand","rock_gravel") for x in li)]
    def observer_census(fam):
        o=[li for _,li in observers.items() if any(x["family"]==fam for x in li)]
        return {"n_distinct_original_observers_with_tagged_family":len(o),
                "n_observers_with_organic_annotation":sum(any(x["family"]==fam and organic_bucket(x["new_category"]) is not None for x in li) for li in o)}
    # Only public original iNaturalist observation URLs. No photos redistributed.
    review=[]
    for x in dat:
        bucket=organic_bucket(x["new_category"])
        if bucket is None:continue
        review.append({"review_id":len(review)+1,
                       "original_observation_url":"https://www.inaturalist.org/observations/"+str(x["id"]),
                       "family":x["family"],
                       "original_raw_substrate_field500":x["value"],
                       "strict_revised_source_class":x["new_category"],
                       "original_photo_view_url":x["photo_url"],
                       "has_public_photo":x["has_photo"],
                       "photo_count":x["n_photos"],
                       "observed_year":x["year"],
                       "human_substrate_visible":"","human_proboscis_contact_visible":"",
                       "human_actual_feeding_confirmed":"","independent_second_reader_agreement":"",
                       "review_status":"NOT_HUMAN_REVIEWED"})
    preview=[{"original_value":v,"source_frequency":c,"old_bin":old_substrate(v),"revised_bin":category(v)}
             for v,c in collections.Counter(x["value"].casefold().strip() for x in dat).most_common()
             if old_substrate(v)!=category(v)]
    mismatches=collections.Counter(x["value"] for x in dat if x["old_category"]!=x["new_category"])
    results={"schema":"GESUI_INAT_V2_FIELD500_746_PHOTO_REVIEW_AND_OBSERVER_SOURCE_QA_V1",
         "original_v1_full_census_git_blob":BASELINE_GIT_BLOB,
         "v2_endpoint":BASE,"v2_sparse_fields":FIELDS,
         "n_original_field500_observations":len(dat),
         "original_source_reported_count":total,
         "original_field500_family_census_exactly_reproduced":True,
         "original_old_regex_category_census_exactly_reproduced":True,
         "source_HTTP_SHA256_page_receipts":pages,
         "old_raw_substrate_bins":dict(sorted(byold.items())),
         "revised_raw_text_bins":dict(sorted(bynew.items())),
         "family_by_revised_raw_text_bins":{f:dict(sorted(c.items())) for f,c in sorted(byfamnew.items())},
         "n_observation_owner_ids":len(observers),
         "top1_observer_share":len(obs_sorted[0][1])/len(dat),
         "top5_observer_share":sum(len(x) for _,x in obs_sorted[:5])/len(dat),
         "n_Lycaenidae_owners":observer_census("Lycaenidae"),
         "n_Nymphalidae_owners":observer_census("Nymphalidae"),
         "n_observers_with_both_Lycaenidae_and_Nymphalidae_annotations":len(both),
         "n_informative_observers_with_both_families_mud_and_organic_annotations":len(informative),
         "top_owner_sensitivity_descriptive_ONLY":samples,
         "n_original_organic_candidate_annotations_for_photo_review":len(review),
         "n_original_organic_candidate_annotations_with_original_photos":sum(x["has_public_photo"] for x in review),
         "n_original_organic_candidate_annotations_without_original_photos":sum(not x["has_public_photo"] for x in review),
         "original_value_reclassification_source_examples":preview[:70],
         "n_raw_records_reassigned_by_vocabulary":sum(mismatches.values()),
         "unreviewed_original_photo_urls_are_NOT_verified_behavior":True,
         "review_manifest_does_not_include_hidden_coordinates_or_observer_ids":True,
         "source_records_do_not_establish_sewage_or_nutrient_assimilation":True,
         "no_WWTP_or_HNS_joint_model_run":True,
         "no_family_puddling_propensity_estimated":True,
         "source_level_decision":"ANNOTATION_TAXON_OBSERVER_PHOTO_REVIEW_READY_NOT_BIOLOGICAL_SELECTIVITY"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(results,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf8")
    args.review.parent.mkdir(parents=True,exist_ok=True)
    with args.review.open("w",newline="",encoding="utf8") as f:
        fw=csv.DictWriter(f,fieldnames=list(review[0]))
        fw.writeheader();fw.writerows(review)
    print("GESUI_INAT_V2_SOURCE_ANNOTATION_QA="+json.dumps({
      "full":len(dat),"photos":sum(x["has_public_photo"] for x in review),
      "review":len(review),"revised_bins":dict(bynew),
      "both_Lyca_Nymph_observers":len(both),"informative_both":len(informative),
      "top5pct":round(results["top5_observer_share"]*100,2),
      "no_wastewater_model":True},sort_keys=True,ensure_ascii=False),flush=True)
if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--baseline",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    a.add_argument("--review",type=Path,required=True)
    main(a.parse_args())
