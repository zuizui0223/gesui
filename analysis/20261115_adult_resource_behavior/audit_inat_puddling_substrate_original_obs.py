#!/usr/bin/env python3
"""Audit original iNaturalist adult puddling substrate field 500 vs mixed field 17428.

User-generated Observation Field 500 offers named SUBSTRATES, unlike the
aggregation/roost/puddling field 17428. Dataset may be very small and observer
selected. This program does NOT test wastewater exposure, family prevalence,
species abundance, or adult uptake of N/Na/protein. Merely checks original
annotation opportunities by taxonomic family and raw substrate vocabulary.
"""
from __future__ import annotations
import argparse,collections,hashlib,json,re,time,urllib.parse,urllib.request
from pathlib import Path

URL="https://api.inaturalist.org/v1/"
INAT_QUERY={
 "substrate500":{"field_id":500,"name":"Butterfly Puddling Substrate"},
 "mixed17428":{"field_id":17428,"name":"Lepidopterans aggregating/Roosting/mud puddling"},
 "feeding16376":{"field_id":16376,"name":"Is observation one of these special types of feeding?"},
 "carrion18609":{"field_id":18609,"name":"Behavior: Consuming carrion"},
}
FAMILIES=("Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae","Riodinidae")
MAX_FIELD_FETCH=1000
CAP_BY_FIELD={500:1000,17428:30,16376:20,18609:20}
SOURCE_PAGE_SIZE=5
REQUEST_TIMEOUT=55
def request(path,params):
    if path not in ("observations","taxa/autocomplete"):
        raise ValueError("Only original iNaturalist source GET API allowed")
    url=URL+path+"?"+urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"gesui-adult-puddling-source-record-audit/1.0","Accept":"application/json"})
            with urllib.request.urlopen(req,timeout=REQUEST_TIMEOUT) as s:
                raw=s.read(5_000_000)
                if len(raw)>=5_000_000:raise RuntimeError("iNat oversized API response even with 5 records: STOP")
            d=json.loads(raw)
            if not isinstance(d,dict):raise ValueError("Not iNat v1 JSON object")
            return d,hashlib.sha256(raw).hexdigest()
        except Exception:
            if attempt==3:raise
            time.sleep(2+attempt)
def get_fields(v):
    return v.get("ofvs") or v.get("observation_field_values") or []
def exact_values(v,id):
    vals=[]
    for q in get_fields(v):
        if not isinstance(q,dict):continue
        qid=q.get("field_id")
        if qid is None:
            qid=q.get("observation_field",{}).get("id")
        if qid==id:
            vals.append(str(q.get("value") or "").strip())
    return vals

def family_ancestors():
    r={};qa={}
    for name in FAMILIES:
        d,s=request("taxa/autocomplete",{"q":name,"per_page":10})
        matches=[z for z in d.get("results",[]) if z.get("name")==name and z.get("rank")=="family"]
        if len(matches)!=1:raise RuntimeError("Taxonomy family search not unique "+name+" matches "+str(matches[:2]))
        item=matches[0];r[name]=int(item["id"]);qa[name]={"original_id":r[name],"taxonomy_search_response_SHA256":s}
    if len(set(r.values()))!=len(r):raise RuntimeError("Duplicated butterfly family source taxon IDs")
    return r,qa

def substrate(v):
    v=(v or "").casefold().strip()
    if not v:return "EMPTY"
    # Priority: require unambiguous singular category to avoid overinterpretation.
    terms={
      "carrion_or_dead_organism":r"carrion|dead\s+(?:fish|animal|body|bird|frog)|carcass|dead\s+matter|decompos|corpse|rotten|decay|fish\s+(?:body|corpse)",
      "dung_feces":r"dung|faec|fec|excrement|scat|droppings|manure|poop",
      "urine_sweat_salt":r"urine|sweat|perspir|salt(?:ed)?\s+(?:water|soil)|uric",
      "mud_soil_sand":r"mud|soil|sand|clay|earth|silt",
      "rock_gravel":r"rock|gravel|stone|pebble",
      "water_only":r"^(?:water|puddle|river|stream|wet(?:\s+ground)?|rainwater)$"
    }
    hits=[k for k,rx in terms.items() if re.search(rx,v)]
    if len(hits)>1:return "AMBIGUOUS_MULTIPLE_SUBSTRATES"
    if len(hits)==1:return hits[0]
    return "UNCLASSIFIED_TEXT"

def cat_taxon(x, fams):
    t=x.get("taxon") or {}
    ids=set(t.get("ancestor_ids") or [])
    if t.get("id"):ids.add(t.get("id"))
    hits=[name for name,n in fams.items() if n in ids]
    return hits[0] if len(hits)==1 else "OTHER_OR_UNCLASSIFIED"

def audit(args):
    fams,qa=family_ancestors()
    records={}
    fieldresults={}
    for tag,req in INAT_QUERY.items():
        target=req["field_id"]
        key="field:"+req["name"]
        found={}
        pages=[]
        first=True;total=None
        for page in range(1,(CAP_BY_FIELD[target]+SOURCE_PAGE_SIZE-1)//SOURCE_PAGE_SIZE+1):
            # Cap at 1000. No "missing field" means absence of behavior.
            if page>1:time.sleep(1.05)  # follow original iNaturalist ~1 req/sec guidance
            d,digest=request("observations",{key:"","per_page":SOURCE_PAGE_SIZE,"page":page,"verifiable":"any",
                                               "order_by":"id","order":"asc"})
            if total is None:total=int(d.get("total_results",0))
            obs=d.get("results")
            if not isinstance(obs,list):raise RuntimeError("No iNat observation result list")
            if len(obs)==0:break
            valid=0
            for x in obs:
                values=exact_values(x,target)
                if not values:continue
                valid+=1
                if not isinstance(x.get("id"),int):continue
                # Do not serialize original observation/user IDs or exact lat/lon;
                # these are only source QA and aggregated family/field counts.
                found[x["id"]]={"values":values,"family":cat_taxon(x,fams),
                     "observer_id":(x.get("user") or {}).get("id"),
                     "observed_year":(x.get("observed_on") or "")[:4],
                     "taxon_rank":(x.get("taxon") or {}).get("rank"),
                     "has_description":bool(x.get("description")),
                     "has_photo":bool(x.get("photos"))}
            pages.append({"page":page,"returned_observations":len(obs),"matching_explicit_field":valid,
                          "original_HTTP_SHA256":digest})
            if valid<len(obs)//2:
                raise RuntimeError("iNaturalist field filter seems not applied correctly "+tag+": "+str(valid)+"/"+str(len(obs)))
            if len(found)>=CAP_BY_FIELD[target] or page*SOURCE_PAGE_SIZE>=total or len(obs)<SOURCE_PAGE_SIZE:break
        if total is None:raise RuntimeError("Source iNat field query returned no total")
        arr=list(found.values())
        if target==500 and len(found)!=total:
            raise RuntimeError("Field 500 exact original full observation census not reached; cannot report family/substrate census")
        family_n=collections.Counter(x["family"] for x in arr)
        raw_values=collections.Counter(v.casefold().strip() for x in arr for v in x["values"])
        cats=collections.Counter(substrate(v) for x in arr for v in x["values"]) if target==500 else None
        fam_substrate=collections.defaultdict(collections.Counter)
        if target==500:
            for x in arr:
                for v in x["values"]:
                    fam_substrate[x["family"]][substrate(v)]+=1
        obs_counts=collections.Counter(x["observer_id"] for x in arr if x["observer_id"] is not None)
        top1=max(obs_counts.values(),default=0)/len(arr) if arr else None
        top5=sum(v for _,v in obs_counts.most_common(5))/len(arr) if arr else None
        is_complete=(len(found)==total)
        if len(found)>total:raise RuntimeError("Impossible paginated more observations than total")
        fieldresults[tag]={
          "field_id":target,"source_name":req["name"],
          "source_total_results":total,"fetched_distinct_observations":len(found),
          "fetch_complete":is_complete,
          "n_source_matching_field_to_returned":len(found),
          "raw_field_value_counts":dict(raw_values.most_common(90)),
          "substrate_text_class_counts_if_field500":dict(cats) if cats is not None else None,
          "family_x_substrate_count_field500_full_census_if_complete":{
            k:dict(v) for k,v in sorted(fam_substrate.items())} if target==500 else None,
          "total_distinct_unstandardized_raw_value_strings":len(raw_values),
          "raw_field_value_counts_published_top90_not_full_vocabulary":len(raw_values)>90,
          "family_ancestor_counts_among_selected_tagged_observations":dict(family_n),
          "n_distinct_tagger_observers_may_not_equal_observers":len(obs_counts),
          "fraction_top1_observation_owner":top1,
          "fraction_top5_observation_owners":top5,
          "taxon_obs_year_counts":dict(collections.Counter(x["observed_year"] for x in arr)),
          "source_HTTP_SHA256_receipts":pages,
          "not_a_behavioral_absence_rate":True,
          "not_population_sampling_prevalence":True,
        }
        print("GESUI_INAT_FIELD_RESULT="+json.dumps({"tag":tag,"total":total,"fetched":len(found),"complete":is_complete,
          "families":dict(family_n),"substrates":dict(cats) if cats is not None else None,
          "top_one_observer_share":top1},ensure_ascii=False,sort_keys=True),flush=True)
    out={"schema":"GESUI_20261115_ORIGINAL_INAT_ADULT_PUDDLING_FIELD500_TAXON_SOURCE_AUDIT_V1",
         "butterfly_family_taxa_from_iNat_authority":qa,
         "fields":fieldresults,
         "all_categorical_raw_substrates_require_human_picture_validation":True,
         "source_field500_annotated_value_not_verified_active_fluid_intake":True,
         "mixed_field17428_not_specific_to_puddling":True,
         "absence_of_field_not_observed_nonpuddling":True,
         "target_composition_cannot_estimate_family_puddling_propensity":True,
         "single_field_tags_can_cluster_by_photographer_and_project":True,
         "raw_observation_photo_coordinates_and_user_ids_not_republished":True,
         "no_H7_WWTP_or_2022_riparian_source_geography_used":True,
         "no_HNS_or_nutrient_assimilation_or_biological_wastewater_effect":True,
         "analysis_type":"POSTHOC_DESCRIPTIVE_ANNOTATION_SOURCE_AUDIT_NOT_CONFIRMATORY_BEHAVIOUR_PREVALENCE"}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    audit(p.parse_args())
