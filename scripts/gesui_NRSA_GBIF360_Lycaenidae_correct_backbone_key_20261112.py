#!/usr/bin/env python3
"""Repair one invalid taxonomy parameter in frozen 360x3 GBIF source requests.

The 2026-11-11 Lycaenidae taxonKey 47923 was from iNaturalist,
not GBIF. Correct the exact original geometries with GBIF Backbone
family 5473, explicit checklist, preserving original controls and gates.
This is a SOURCE correction, NOT a fitted ecological response.
"""
import collections,concurrent.futures,hashlib,json,os,pathlib,time,urllib.parse,urllib.request,urllib.error
P=pathlib.Path(os.environ["GBIF_360_SOURCE_ARTIFACT"])
OUT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"nrsa_gbif_correct_taxon"
OUT.mkdir(parents=True,exist_ok=True)
origin=json.loads((P/"GBIF_360x3_SOURCE_SUPPORT_REPORT.json").read_text())
rows=json.loads((P/"GBIF_360x3_ORIGINAL_COUNTS_WITH_URL_SHA_ONLY.json").read_text())
sitebytes=(P/"NRSA1819_360_GBIF_PRE_COUNTS_SITE_MANIFEST.json").read_bytes()
sha=lambda b:hashlib.sha256(b).hexdigest()
if sha(sitebytes)!="bff6f7f268232c1cb11296cdf19c377fc58ce4847abc53925ecb0407172f130e":
    raise RuntimeError("Original 360-site source geography changed")
if len(rows)!=1080 or origin["queries_ok"]!=1023 or origin["queries_missing"]!=57:
    raise RuntimeError("Original GBIF first-year source report changed")
old_lyc={r["site_id"]:r for r in rows if r["taxon"]=="Lycaenidae"}
old_plants={r["site_id"]:r for r in rows if r["taxon"]=="Plantae"}
old_animals={r["site_id"]:r for r in rows if r["taxon"]=="Animalia"}
if set(old_lyc)!=set(old_plants) or set(old_lyc)!=set(old_animals) or len(old_lyc)!=360:
    raise RuntimeError("Original site / taxon source crosscheck failed")
ORIGINAL_CHECKLIST="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
def get(url,timeout=60):
    request=urllib.request.Request(url,headers={
      "User-Agent":"gesui-GBIF-taxonomy-correction-study-20261112/1.0",
      "Accept":"application/json"})
    with urllib.request.urlopen(request,timeout=timeout) as response:
        return response.read()
# Taxon identity is verified before any new 360-site occurrence numbers.
taxon_url="https://api.gbif.org/v1/species/5473"
taxon_raw=get(taxon_url)
taxon=json.loads(taxon_raw)
if "lycaenidae" not in str(taxon.get("scientificName","")).lower() or taxon.get("rank")!="FAMILY":
    raise RuntimeError("GBIF taxon 5473 is not verified Lycaenidae FAMILY: "+str(taxon)[:500])
source_provenance={"gbif_verified_key":5473,
  "taxon_url":taxon_url,"taxon_response_sha256":sha(taxon_raw),
  "scientificName":taxon.get("scientificName"),"rank":taxon.get("rank"),
  "original_invalid_inaturalist_key":47923,
  "explicit_gbif_legacy_checklist":ORIGINAL_CHECKLIST}
print("Verified original GBIF taxonomy "+json.dumps(source_provenance,sort_keys=True),flush=True)
# An additional distinct no-geometry precheck differentiates an invalid
# family key from valid GBIF source with no occurrence records.
check_params={"taxonKey":"5473","checklistKey":ORIGINAL_CHECKLIST,
              "year":"2017,2021","hasCoordinate":"true","limit":"0",
              "occurrenceStatus":"PRESENT"}
global_url="https://api.gbif.org/v1/occurrence/search?"+urllib.parse.urlencode(check_params)
global_raw=get(global_url)
global_j=json.loads(global_raw)
if not isinstance(global_j.get("count"),int) or global_j["count"]<=0:
    raise RuntimeError("Correct GBIF Lycaenidae key has no global occurrences — refuse 360 requests")
source_provenance["global_source_count_positive"]=True
source_provenance["global_source_SHA256"]=sha(global_raw)
def create_url(old):
    p=urllib.parse.urlsplit(old)
    q=urllib.parse.parse_qsl(p.query,keep_blank_values=True)
    have=dict(q)
    if have.get("taxonKey")!="47923" or have.get("year")!="2017,2021" \
       or have.get("limit")!="0" or have.get("hasCoordinate")!="true" \
       or have.get("occurrenceStatus")!="PRESENT" or not have.get("geometry"):
        raise RuntimeError("Original frozen GBIF source parameters unexpectedly changed")
    q=[(k,("5473" if k=="taxonKey" else v)) for k,v in q]
    q.append(("checklistKey",ORIGINAL_CHECKLIST))
    return urllib.parse.urlunsplit((p.scheme,p.netloc,p.path,urllib.parse.urlencode(q),p.fragment))
def query_one(source):
    site=source["site_id"]
    url=create_url(source["source_url"])
    last=None
    for k in range(4):
        try:
            raw=get(url,timeout=65)
            j=json.loads(raw)
            if not isinstance(j.get("count"),int) or j["count"]<0:
                raise RuntimeError("Missing/noninteger GBIF official count response")
            return {
                "site_id":site,"taxon":"Lycaenidae",
                "state":source["state"],"eco9":source["eco9"],"huc8":source["huc8"],
                "status":"ok","record_count":j["count"],
                "source_url":url,"response_sha256":sha(raw)}
        except Exception as e:
            last=type(e).__name__+":"+str(e)[:120]
            time.sleep(1+3*k)
    return {
        "site_id":site,"taxon":"Lycaenidae",
        "state":source["state"],"eco9":source["eco9"],"huc8":source["huc8"],
        "status":"missing","record_count":None,
        "source_url":url,"error":last}
repaired=[]
sources=[old_lyc[k] for k in sorted(old_lyc)]
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
    futures=[pool.submit(query_one,s) for s in sources]
    for i,f in enumerate(concurrent.futures.as_completed(futures)):
        repaired.append(f.result())
        if (i+1)%90==0:print("Correct GBIF Lycaenidae original-source requests "+str(i+1)+"/360",flush=True)
repaired.sort(key=lambda r:r["site_id"])
def qcount(vals,p):
    if not vals:return None
    v=sorted(vals)
    return v[min(len(v)-1,round((len(v)-1)*p))]
counts=[r["record_count"] for r in repaired if r["status"]=="ok"]
R={r["site_id"]:r for r in repaired}
comp=[];eligible=[]
for sid in sorted(old_lyc):
    rpl=old_plants[sid];ran=old_animals[sid];rly=R[sid]
    if rpl["status"]=="ok" and ran["status"]=="ok" and rly["status"]=="ok":
        comp.append(sid)
        if rpl["record_count"]>=50 and rly["record_count"]>=20:
            eligible.append(sid)
states={R[x]["state"] for x in eligible}
eco={str(R[x]["eco9"]) for x in eligible}
gate=(len(comp)>=250 and len(eligible)>=100 and len(states)>=15 and len(eco)>=5)
result={
 "schema":"NRSA2018_19_GBIF360_LYCAENIDAE_CORRECT_TAXON_BACKBONE_SOURCE_v1",
 "original_gbif_source_run":43036787547,
 "source_geography_manifest_sha256":"bff6f7f268232c1cb11296cdf19c377fc58ce4847abc53925ecb0407172f130e",
 "corrected_taxon_provenance":source_provenance,
 "original_plants_and_animal_controls_not_requeried":True,
 "new_query_count":len(repaired),"corrected_lycae_ok":len(counts),
 "corrected_lycae_missing":360-len(counts),
 "corrected_lycae_positive":sum(x>0 for x in counts),
 "corrected_lycae_ge20":sum(x>=20 for x in counts),
 "corrected_lycae_quantiles":[qcount(counts,x) for x in [0,.25,.5,.75,.9,.95,.99,1]],
 "three_taxon_complete":len(comp),
 "original_joint_primary_gate_plant_ge50_corrected_lycae_ge20":len(eligible),
 "joint_eligible_states":len(states),"joint_eligible_eco9":len(eco),
 "unchanged_original_gate":{"complete":250,"joint":100,"states":15,"eco9":5},
 "gate_pass":gate,
 "decision":"PASS_TO_POSTHOC_RECORDED_TAXA_ONLY" if gate else "STOP_NO_CROSS_ECOSYSTEM_ASSOCIATION",
 "no_ecological_effect_models_run":True,
 "no_undocumented_taxon_radius_year_threshold_tuning":True
}
# Original controls + repaired taxon source arrays carry URLs and response SHA.
# Preserve all fixed-site metadata in short-lived Actions artifact; no site
# counts or coordinates committed to repo.
(PARENT:=OUT/"NRSA_GBIF360_REPAIRED_LYCAENIDAE_360_RECORDS_URL_SHA.json").write_text(json.dumps(repaired,indent=1))
(OUT/"NRSA_GBIF360_TAXON_KEY_CORRECTION_SOURCE_GATE.json").write_text(json.dumps(result,indent=2))
(OUT/"NRSA_GBIF360_QUALIFIED_SITE_ID_SNAPSHOT_SOURCE_ONLY.json").write_text(json.dumps({
  "sites":eligible,"original_exact_360_source_sha256":result["source_geography_manifest_sha256"],
  "eligible_site_count":len(eligible),"not_response_selected":True
},indent=2))
print("GBIF360_LYCAENIDAE_KEY_REPAIR_SOURCE_RESULT_START")
print(json.dumps(result,sort_keys=True))
print("GBIF360_LYCAENIDAE_KEY_REPAIR_SOURCE_RESULT_END")
