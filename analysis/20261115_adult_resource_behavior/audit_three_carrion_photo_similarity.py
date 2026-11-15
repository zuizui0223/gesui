#!/usr/bin/env python3
"""Noncommercial source-review of original CC-BY-NC iNaturalist three animal-remains cases.

Source fields, original image licence & original public positions authenticated.
Download original photos ONLY transiently into process RAM, compute hashes/ORB
similarity; NEVER publish photos, image pixels, raw GPS or owner identities.
This cannot identify a carcass or confirm proboscis contact without human
review; do not call any heuristically matching pictures "same feeding event".
"""
from __future__ import annotations
import argparse,hashlib,io,json,math,re,time,urllib.parse,urllib.request
from pathlib import Path
from itertools import combinations
import numpy as np,cv2
from PIL import Image

API="https://api.inaturalist.org/v2/observations"
CASES=(150764924,150764928,367297297)
EXPECTED_BY_OBS={150764924:[260028346],150764928:[260028363,260029204],
                 367297297:[670748495]}
EXPECTED_VALUE={150764924:"animal remains",150764928:"animal remains",
                 367297297:"carcass of amphibian"}
FIELDS="id,observed_on,user.id,geojson.coordinates,photos.id,photos.url,photos.license_code,ofvs.name,ofvs.value,description"
MAX_PHOTO_BYTES=5_000_000
EXPECTED_LICENCES={"cc-by-nc"}
def query(url,expect_json=False):
    last=None
    for try_n in range(3):
        try:
            if not (url.startswith(API) or url.startswith("https://inaturalist-open-data.s3.amazonaws.com/photos/") or url.startswith("https://static.inaturalist.org/photos/")):
                raise ValueError("Unexpected non-approved source: "+url.split("?")[0][:60])
            req=urllib.request.Request(url,headers={
                "User-Agent":"GESUI-photo-original-source-duplicate-CC-BY-NC-no-redistribution/1.0"})
            with urllib.request.urlopen(req,timeout=60) as r:
                mime=r.headers.get("Content-Type","")
                raw=r.read(MAX_PHOTO_BYTES+1)
            if len(raw)>MAX_PHOTO_BYTES:raise RuntimeError("Original image >5MB; STOP")
            if expect_json:
                if not raw.lstrip().startswith(b"{"):raise RuntimeError("Non-JSON response")
                return json.loads(raw),hashlib.sha256(raw).hexdigest()
            if not mime.startswith("image/"):raise RuntimeError("Source photo is not image media: "+mime)
            return raw,hashlib.sha256(raw).hexdigest()
        except Exception as e:
            last=e
            if try_n<2:time.sleep(try_n+1)
    raise RuntimeError(str(last))
def exact_cases():
    params={"id":",".join(str(x) for x in CASES),"per_page":10,"fields":FIELDS}
    data,digest=query(API+"?"+urllib.parse.urlencode(params),expect_json=True)
    arr=data.get("results") or []
    if len(arr)!=3 or {x.get("id") for x in arr}!=set(CASES):
        raise RuntimeError("Original source 3 precise cases not found")
    return {x["id"]:x for x in arr},digest
def hh(g):
    if not isinstance(g,dict):return None
    coords=g.get("coordinates")
    if not isinstance(coords,list) or len(coords)!=2:return None
    lon,lat=coords
    if not isinstance(lon,(int,float)) or not isinstance(lat,(int,float)):
        return None
    if not -180<=lon<=180 or not -90<=lat<=90:return None
    return [float(lat),float(lon)]
def hav(a,b):
    if a is None or b is None:return None
    la,lo=a;lb,lob=b
    dlat=math.radians(la-lb);dlon=math.radians(lo-lob)
    v=math.sin(dlat/2)**2+math.cos(math.radians(la))*math.cos(math.radians(lb))*math.sin(dlon/2)**2
    return 12742*math.asin(min(1,math.sqrt(v)))
def loadpix(photo):
    src=str(photo.get("url") or "")
    ident=photo.get("id")
    if not isinstance(ident,int) or "/"+str(ident)+"/" not in src:
        raise RuntimeError("Original photo ID does not occur in URL")
    lic=photo.get("license_code")
    if lic not in EXPECTED_LICENCES:
        return {"photo_id":ident,"licence":lic,"status":"SOURCE_LICENCE_NOT_AUTHORIZED_FOR_TEMP_IMAGE_QA"}
    # iNaturalist documents alternate sizes; temporary medium is sufficient
    # for PHash and ORB matching and never included in public Actions artifacts.
    medium=re.sub(r"/(?:square|small|thumb|large|original)\.",r"/medium.",src)
    if medium==src:raise RuntimeError("Unknown original photo size URL")
    if not medium.startswith(("https://inaturalist-open-data.s3.amazonaws.com/photos/",
                               "https://static.inaturalist.org/photos/")):
        raise RuntimeError("Original iNat photo host not approved")
    try:
        body,dig=query(medium,expect_json=False)
        im=Image.open(io.BytesIO(body))
        im.verify()
        im=Image.open(io.BytesIO(body)).convert("RGB")
        if im.width*im.height>5_000_000:raise RuntimeError("Image pixel count too large")
        array=np.array(im)
        # In memory only. Dispose body after features computed.
        gray=cv2.cvtColor(array,cv2.COLOR_RGB2GRAY)
        gray32=cv2.resize(gray,(32,32),interpolation=cv2.INTER_AREA).astype(np.float32)
        dct=cv2.dct(gray32)[:8,:8].flatten()
        phash="".join("1" if x>np.median(dct[1:]) else "0" for x in dct)
        orb=cv2.ORB_create(nfeatures=750)
        kp,desc=orb.detectAndCompute(gray,None)
        return {"photo_id":ident,"licence":lic,"status":"TEMPORARY_SOURCE_ONLY_IMAGE_FEATURES_COMPUTED",
                "sha256_original_medium_photo_source":dig,"source_medium_pixel_dimensions":[im.width,im.height],
                "pHash64":phash,"n_ORB_keypoints":len(kp) if kp else 0,
                "_desc":desc,"_kp_xy":np.float32([k.pt for k in kp]) if kp else np.empty((0,2),np.float32)}
    except Exception as e:
        return {"photo_id":ident,"licence":lic,"status":"ORIGINAL_LICENSED_PHOTO_PIXELS_UNAVAILABLE",
                "error_type":type(e).__name__}
def match(a,b):
    if a.get("status")!="TEMPORARY_SOURCE_ONLY_IMAGE_FEATURES_COMPUTED" or b.get("status")!="TEMPORARY_SOURCE_ONLY_IMAGE_FEATURES_COMPUTED":
        return {"status":"NO_BOTH_SOURCE_IMAGE_FEATURES"}
    hd=sum(x!=y for x,y in zip(a["pHash64"],b["pHash64"]))
    k1=a.get("_desc");k2=b.get("_desc")
    goodmatches=[]
    if k1 is not None and k2 is not None and len(k1)>1 and len(k2)>1:
        m=cv2.BFMatcher(cv2.NORM_HAMMING,crossCheck=False).knnMatch(k1,k2,k=2)
        goodmatches=[pair[0] for pair in m if len(pair)==2 and pair[0].distance<0.7*pair[1].distance]
    geometric_inliers=0
    if len(goodmatches)>=8:
        axy=a["_kp_xy"];bxy=b["_kp_xy"]
        src=np.float32([axy[v.queryIdx] for v in goodmatches]).reshape(-1,1,2)
        dst=np.float32([bxy[v.trainIdx] for v in goodmatches]).reshape(-1,1,2)
        homography,mask=cv2.findHomography(src,dst,cv2.RANSAC,5.0)
        if homography is not None and mask is not None:
            geometric_inliers=int(mask.ravel().sum())
    return {"status":"COMPUTABLE_SOURCE_FEATURE_DISTANCE_NOT_SAME_EVENT_PROOF",
            "pHash64_hamming":hd,
            "ORB_Lowe_ratio0p7_matching_keypoint_count":int(len(goodmatches)),
            "ORB_RANSAC_5px_homography_inliers":geometric_inliers,
            "ORB_RANSAC_inlier_fraction":round(geometric_inliers/len(goodmatches),4) if goodmatches else None,
            "automatic_same_butterfly_or_carcass_identity_unproven":True}
def main(args):
    obs,sha=exact_cases()
    images={};rows=[]
    for id in CASES:
        x=obs[id]
        vals=[q.get("value") for q in x.get("ofvs") or [] if q.get("name")=="Butterfly Puddling Substrate"]
        if vals!=[EXPECTED_VALUE[id]]:raise RuntimeError("Original field500 value changed for "+str(id))
        ids=[y.get("id") for y in x.get("photos") or []]
        if ids!=EXPECTED_BY_OBS[id]:raise RuntimeError("Original observed photo IDs changed for "+str(id))
        for item in x.get("photos") or []:
            feat=loadpix(item);images[item["id"]]=feat
        rows.append({"original_observation_url":"https://www.inaturalist.org/observations/"+str(id),
            "original_source_field500_value":vals[0],"reported_observation_day":x.get("observed_on"),
            "reported_photo_ids":ids,
            "original_photo_licences":[q.get("license_code") for q in x.get("photos") or []],
            "source_public_geojson_present":hh(x.get("geojson")) is not None,
            "original_description_is_present":bool(x.get("description")),
            "original_source_substrate_confirmed_visually":False,
            "original_source_proboscis_contact_confirmed_visually":False})
    pair_distance=hav(hh(obs[CASES[0]].get("geojson")),hh(obs[CASES[1]].get("geojson")))
    pairs=[]
    for i,j in combinations([260028346,260028363,260029204],2):
        pairs.append({"original_photo_ids":[i,j],**match(images[i],images[j])})
    # No image features output where they could enable unintended reconstruction
    image_out=[{k:v for k,v in feat.items() if k not in ("_desc","_kp_xy","pHash64")} for feat in images.values()]
    source={"schema":"GESUI_I_NAT_3_ANIMAL_REMAINS_PHOTO_FEATURES_AND_ORIGINAL_SITE_QA_V1",
       "source_api":"https://api.inaturalist.org/v2/observations",
       "original_three_source_response_SHA256":sha,
       "observation_records":rows,"n_observations":3,
       "n_distinct_observer_days":2,
       "same_owner_same_day_2015_public_site_distance_km":round(pair_distance,4) if pair_distance is not None else None,
       "source_original_photo_features":image_out,
       "photo_pair_source_similarity_NO_BEHAVIOR_CONCLUSION":pairs,
       "licensed_noncommercial_research_transient_image_QA_only":True,
       "source_original_photo_copyright":"CC BY-NC; photos never redistributed or committed",
       "no_original_photo_pixel_bytes_stored_in_GitHub_results":True,
       "no_manual_visual_photograph_review":True,
       "no_visible_proboscis_or_carcase_identity_confirmed":True,
       "no_HNS_WWTP_or_measured_nutrient_exposure":True,
       "nonindependence_warning":"Same owner/date does not establish identical carcass; photo similarity is noncausal."
    }
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(source,sort_keys=True,indent=2,ensure_ascii=False)+"\n")
    print("GESUI_CARRION_PHOTO_SOURCE_FEATURE_QA="+json.dumps({
       "n_photo_features_successful":sum(v["status"]=="TEMPORARY_SOURCE_ONLY_IMAGE_FEATURES_COMPUTED" for v in image_out),
       "pair_distance_km":source["same_owner_same_day_2015_public_site_distance_km"],
       "pairs":pairs,"n_image_bytes_committed":0,
       "human_proboscis_confirmed":False},sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    main(p.parse_args())
