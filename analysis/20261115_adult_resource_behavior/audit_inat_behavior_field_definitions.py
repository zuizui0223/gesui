#!/usr/bin/env python3
"""iNaturalist ORIGINAL observation-field schema and specificity audit (no behavioral outcomes).

Surveys the original unregulated user-created observation fields. Does NOT
read butterfly occurrence counts, photo URLs, HNS, mapped outfall/source exposure
or inferred adult feeding labels. Behaviour categories and false positive risks
must be assessed BEFORE any posthoc family ranking/field numerator.
"""
from __future__ import annotations
import argparse,hashlib,json,re,time,urllib.parse,urllib.request
from pathlib import Path

Q=("puddling","mud","carrion","dung","feces","feeding","feeding behavior","behavior","roosting","aggregat")
ROOT="https://www.inaturalist.org/observation_fields.json"
KNOWN_MIXED_FIELD_ID=17428
API_PAGE_LIMIT=75
MAX_RESULTS_PER_TERM=150
MAX_BYTES=2200000
def get(url):
    if not url.startswith(ROOT) and not url.startswith("https://www.inaturalist.org/observation_fields/"):
        raise ValueError("Unexpected resource. Only official original observation_field metadata is allowed")
    last=None
    for attempt in range(3):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-adult-behavior-field-source-only/1.0",
                                                   "Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=50) as stream:
                raw=stream.read(MAX_BYTES+1)
                if len(raw)>MAX_BYTES:raise RuntimeError("More than 2.2MB field metadata")
                media=stream.headers.get("Content-Type","")
            if not raw.lstrip().startswith((b"{",b"[")):
                raise RuntimeError("Non-JSON field metadata "+repr(raw[:40])+" "+media)
            return json.loads(raw),hashlib.sha256(raw).hexdigest()
        except Exception as e:
            last=str(e)
            if attempt<2:time.sleep(1+attempt)
    return None,str(last)
def items(d):
    if isinstance(d,list):return d,None
    if isinstance(d,dict):
        for k in ["results","observation_fields","data"]:
            if isinstance(d.get(k),list):return d[k],d.get("total_results") or d.get("total")
        if isinstance(d.get("id"),int):return [d],1
    raise RuntimeError("Unknown original iNat legacy observation_fields API shape "+str(type(d)))
def short(d):
    return {"id":d.get("id"),"name":d.get("name"),"description":(d.get("description") or "")[:600],
           "datatype":d.get("datatype"),"allowed_values":d.get("allowed_values"),
           "created_at":d.get("created_at"),"updated_at":d.get("updated_at"),
           "not_independently_verified_behavior":True}
def audit(args):
    results={}
    registry={}
    for q in Q:
        out=[]
        receipts=[]
        for page in (1,2):
            params=urllib.parse.urlencode({"q":q,"page":page,"per_page":API_PAGE_LIMIT})
            data,digest=get(ROOT+"?"+params)
            if data is None:
                receipts.append({"page":page,"source_error":digest})
                break
            arr,total=items(data)
            receipts.append({"page":page,"SHA256":digest,"items":len(arr),"reported_total":total})
            for raw in arr:
                if isinstance(raw,dict) and isinstance(raw.get("id"),int):
                    z=short(raw)
                    out.append(z)
                    registry[z["id"]]=z
            if len(arr)<API_PAGE_LIMIT or len(out)>=MAX_RESULTS_PER_TERM:break
        results[q]={"metadata_results_n":len(out),"receipts":receipts,
                    "field_ids":[z["id"] for z in out]}
    mixed=registry.get(KNOWN_MIXED_FIELD_ID)
    # User-generated 17428 definition admits roosting/aggregation, thus
    # an affirmative value cannot establish adult ingestion of any substrate.
    if mixed is None:
        url=f"https://www.inaturalist.org/observation_fields/{KNOWN_MIXED_FIELD_ID}.json"
        raw,digest=get(url)
        if raw is not None:
            data,total=items(raw)
            mixed=short(data[0])
            registry[mixed["id"]]=mixed
    if mixed is not None and "pudd" not in (mixed["name"]+" "+mixed["description"]).lower():
        raise RuntimeError("Original mixed field 17428 no longer has expected naming; needs manual review")
    fields=sorted(registry.values(),key=lambda z:z["id"])
    restricted=[z for z in fields if re.search(r"puddl|mud|carrion|dung|fec|feeding|ingest|drink|water|roost|aggregat",z["name"]+" "+z["description"],re.I)]
    def is_strict(name,desc):
        all=(name+" "+desc).lower()
        if any(k in all for k in ["roost","aggregat","nectar","mate","copulat"]):
            return False
        return bool(re.search(r"mud.?puddl|carrion|corpse|feces|excrement|dung|feeding on dead|mineral feed",all))
    strict=[z for z in restricted if is_strict(z["name"],z["description"])]
    result={"schema":"GESUI_20261115_INAT_ADULT_PUDDLING_SUBSTRATE_ORIGINAL_FIELD_AUDIT_V1",
        "query_terms":list(Q),"source_endpoint":ROOT,"n_unique_fields":len(fields),
        "n_adult_keyword_candidate_fields":len(restricted),
        "strict_name_description_candidate_field_ids":[x["id"] for x in strict],
        "strict_candidates_not_photo_or_value_validated":True,
        "mixed_original_field_17428":mixed,
        "mixed_field_17428_cannot_claim_puddling_ingestion":True,
        "field_17428_might_include_rest_aggregation_and_nectar":True,
        "field_absence_not_true_behavior_absence":True,
        "fields_are_user_created_and_not_curator_standardized":True,
        "request_receipts":results,
        "candidate_fields":restricted,
        "no_species_or_observations_read":True,
        "no_sewage_water_quality_or_map_crossing":True,
        "no_HNS_computed":True,
        "no_adult_resource_family_contrast_run":True,
        "decision":"SOURCE_ONLY_CLASSIFICATION_VALIDITY_NOT_BIOLOGICAL_EFFECT"}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n")
    print("GESUI_INAT_ADULT_FIELD_SCHEMA="+json.dumps({
         "n_unique_fields":len(fields),"n_keyword_fields":len(restricted),
         "strict_id_list":[z["id"] for z in strict],
         "mixed_id_found":mixed is not None,
         "biological_analyses_run":False},ensure_ascii=False,sort_keys=True),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--output",required=True,type=Path)
    audit(p.parse_args())
