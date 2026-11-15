#!/usr/bin/env python3
"""Two 2018 original cow-manure observation clusters: source-only location/photo QA.

Original field500 Lycaenidae observations 16195221 and 16195223 have
license_code=None on their public photos: NEVER download those image pixels.
Papilionidae 18670197 and 18670234 photos are original CC BY-NC; for a
noncommercial source independence comparison, download ONLY authorized medium
frames transiently, compute ORB inlier counts, dispose pixels. No images or
observer IDs/GPS coordinates in public JSON.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,math,pathlib,time,urllib.parse,urllib.request
import cv2
from audit_three_carrion_photo_similarity import loadpix,match,hav,hh

API="https://api.inaturalist.org/v2/observations"
GROUPS={"Lycaenidae":[16195221,16195223],"Papilionidae":[18670197,18670234]}
FIELDS="id,observed_on,user.id,geojson.coordinates,photos.id,photos.url,photos.license_code,ofvs.name,ofvs.value"
def source():
    q={"id":",".join(map(str,sum(GROUPS.values(),[]))),
       "per_page":10,"fields":FIELDS}
    u=API+"?"+urllib.parse.urlencode(q)
    req=urllib.request.Request(u,headers={
        "User-Agent":"GESUI-dung46-original-2018-license-gated-source-cluster-audit",
        "Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=60) as r:
        raw=r.read(2_000_000)
    d=json.loads(raw)
    records=d.get("results") or []
    if len(records)!=4 or {x.get("id") for x in records}!=set(sum(GROUPS.values(),[])):
        raise RuntimeError("Original 2018 four source observations changed")
    return {x["id"]:x for x in records},hashlib.sha256(raw).hexdigest()
def run(args):
    dat,sha=source()
    out={}
    for family,ids in GROUPS.items():
        a,b=[dat[i] for i in ids]
        vals=[[str(f.get("value") or "").strip() for f in row.get("ofvs") or []
            if f.get("name")=="Butterfly Puddling Substrate"] for row in (a,b)]
        if vals!=[["cattle manure"],["cattle manure"]]:
            raise RuntimeError("Original cow manure field name/value changed "+family)
        users=[(q.get("user") or {}).get("id") for q in (a,b)]
        dates=[q.get("observed_on") for q in (a,b)]
        if users[0]!=users[1] or dates[0]!=dates[1]:
            raise RuntimeError("Expected same original owner/date in 2018 source-pair "+family)
        ph=[q.get("photos") or [] for q in (a,b)]
        nphoto=[len(x) for x in ph]
        lic=collections.Counter(str(z.get("license_code")) for row in ph for z in row)
        ids_overlap={x.get("id") for x in ph[0]}&{x.get("id") for x in ph[1]}
        geom=hav(hh(a.get("geojson")),hh(b.get("geojson")))
        s={"source_original_observation_links":[
             "https://www.inaturalist.org/observations/"+str(i) for i in ids],
           "same_original_observation_owner":True,"same_original_observation_day":True,
           "original_public_coordinate_separation_km":round(geom,5) if geom is not None else None,
           "original_photo_counts":nphoto,
           "photo_licence_code_counts":dict(lic),
           "shared_identical_original_photo_ids":len(ids_overlap),
           "both_field500_cattle_manure_text":True,
           "source_photo_pixels_visually_reviewed":False,
           "adult_proboscis_contact_visually_confirmed":False,
           "independent_dung_pile_or_individual_count_NOT_established":True}
        if family=="Lycaenidae":
            s["photo_pixels_downloaded"]=False
            s["reason"]="Original photo licence_code None for both observations, so images not downloaded"
        else:
            # Only original CC BY-NC data used; some photo licence may differ,
            # which loadpix rejects. Compute max homography inlier concordance
            # across photographic angles, not ecological identical event.
            imgs0=[loadpix(z) for z in ph[0]]
            imgs1=[loadpix(z) for z in ph[1]]
            val0=[z for z in imgs0 if z["status"]=="TEMPORARY_SOURCE_ONLY_IMAGE_FEATURES_COMPUTED"]
            val1=[z for z in imgs1 if z["status"]=="TEMPORARY_SOURCE_ONLY_IMAGE_FEATURES_COMPUTED"]
            matches=[]
            for x in val0:
                for y in val1:
                    v=match(x,y)
                    matches.append({"photo_pair_ids":[x["photo_id"],y["photo_id"]],
                                    "inliers":v.get("ORB_RANSAC_5px_homography_inliers",0),
                                    "ratio07_matches":v.get("ORB_Lowe_ratio0p7_matching_keypoint_count",0),
                                    "hamming":v.get("pHash64_hamming")})
            matches.sort(key=lambda z:(z["inliers"],z["ratio07_matches"]),reverse=True)
            s.update({"source_allowed_cc_by_nc_images_processed_in_memory_only":len(val0)+len(val1),
                      "source_photos_rejected_or_unavailable":nphoto[0]+nphoto[1]-len(val0)-len(val1),
                      "highest_original_cross_observation_photo_feature_alignment":matches[:3],
                      "n_photo_pair_comparisons_with_geometric_8plus_inliers":sum(q["inliers"]>=8 for q in matches),
                      "n_original_photo_pixels_committed_to_GitHub":0,
                      "same_cow_manure_pile_identity_not_proven":True})
        out[family]=s
    res={"schema":"GESUI_I_NAT_2018_LYCA_PAPILIO_DUNG_2PAIRS_SOURCE_LICENSE_GATED_GEOGRAPHY_ORB_V1",
         "source_iNaturalist_v2_response_SHA256":sha,
         "family_source_only_pairs":out,
         "all_four_are_original_2018_cattle_manure_annotations":True,
         "n_original_annotated_observations":4,
         "n_original_owner_day_clusters":2,
         "number_of_truly_independent_dung_resources_NOT_KNOWN":True,
         "no_proboscis_or_feeding_visual_confirmation":True,
         "no_sewage_exposure_or_nutrient_concentrations":True,
         "all_original_photo_pixels_absent_from_results":True}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(res,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_I_NAT_2018_COWMANURE_SOURCE_PAIR_QA="+json.dumps({
       "lycaenid_public_coord_distance_km":out["Lycaenidae"]["original_public_coordinate_separation_km"],
       "papilionid_public_coord_distance_km":out["Papilionidae"]["original_public_coordinate_separation_km"],
       "lycaenid_photos_read":False,
       "papilionid_photos_processed":out["Papilionidae"]["source_allowed_cc_by_nc_images_processed_in_memory_only"],
       "best_image_pair":out["Papilionidae"]["highest_original_cross_observation_photo_feature_alignment"][:1],
       "verified_feeding_events":0,
       "same_dung_pile_unproven":True},sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=pathlib.Path)
    run(p.parse_args())
