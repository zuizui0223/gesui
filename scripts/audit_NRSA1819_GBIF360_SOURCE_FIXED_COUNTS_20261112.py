#!/usr/bin/env python3
"""Source-only diagnostic of existing, frozen 360-site x 3 GBIF counts.

NO new GBIF requests and no ecological coefficients. NEVER rescue original gate
by lowering observation thresholds. Validates saved URL+SHA rows and counts.
Original 2026-11-11 source run: 43036787547, artifact 13864116088.
"""
import collections,hashlib,json,os,pathlib,statistics,urllib.parse
base=pathlib.Path(os.environ["GBIF_360_SOURCE_ARTIFACT"])
result=json.loads((base/"GBIF_360x3_SOURCE_SUPPORT_REPORT.json").read_text())
rows=json.loads((base/"GBIF_360x3_ORIGINAL_COUNTS_WITH_URL_SHA_ONLY.json").read_text())
manifest=(base/"NRSA1819_360_GBIF_PRE_COUNTS_SITE_MANIFEST.json").read_bytes()
assert hashlib.sha256(manifest).hexdigest()=="bff6f7f268232c1cb11296cdf19c377fc58ce4847abc53925ecb0407172f130e"
assert result["queried_original_sites"]==360 and len(rows)==1080 and result["queries_ok"]==1023
tmap={"Plantae":"6","Lycaenidae":"47923","Animalia":"1"}
seen=set();by=collections.defaultdict(dict);mismatch=[]
for r in rows:
    key=(r["site_id"],r["taxon"])
    assert key not in seen,key
    seen.add(key);by[r["site_id"]][r["taxon"]]=r
    q=urllib.parse.parse_qs(urllib.parse.urlparse(r["source_url"]).query)
    if q.get("taxonKey")!=[tmap[r["taxon"]]] or q.get("year")!=["2017,2021"] \
        or q.get("limit")!=["0"] or q.get("hasCoordinate")!=["true"] \
        or q.get("occurrenceStatus")!=["PRESENT"] or "geometry" not in q:
        mismatch.append(r["taxon"])
assert not mismatch
assert len(by)==360 and all(set(v)==set(tmap) for v in by.values())
def pct(x,q):
    x=sorted(x)
    if not x:return None
    return x[min(len(x)-1,round((len(x)-1)*q))]
summ={}
for label in tmap:
    sub=[r for r in rows if r["taxon"]==label]
    ok=[int(r["record_count"]) for r in sub if r["status"]=="ok"]
    assert len(ok)+sum(r["status"]!="ok" for r in sub)==360
    summ[label]={
        "requests_success":len(ok),"missing":360-len(ok),
        "record_zeros":sum(x==0 for x in ok),
        "record_count_positive":sum(x>0 for x in ok),
        "record_ge_1":sum(x>=1 for x in ok),
        "record_ge_5":sum(x>=5 for x in ok),
        "record_ge_10":sum(x>=10 for x in ok),
        "record_ge_20":sum(x>=20 for x in ok),
        "record_ge_50":sum(x>=50 for x in ok),
        "quantiles":[pct(ok,q) for q in [0,.1,.25,.5,.75,.9,.95,.99,1]],
        "error_categories":dict(collections.Counter(
            (r.get("error") or "missing-without-reason").split(":")[0]
            for r in sub if r["status"]!="ok"))
    }
complete=[]
for sid,rr in by.items():
    if all(r["status"]=="ok" for r in rr.values()):
        complete.append({label:int(rr[label]["record_count"]) for label in tmap})
joint={}
for p in [1,5,10,20,50]:
    for b in [1,5,10,20]:
        joint[f"plantae_ge{p}_lycaenidae_ge{b}"]=sum(z["Plantae"]>=p and z["Lycaenidae"]>=b for z in complete)
out={
 "original_source_gate_pass":result["gate_pass"],
 "original_source_gate_decision":result["decision"],
 "original_source_sites_sha256":"bff6f7f268232c1cb11296cdf19c377fc58ce4847abc53925ecb0407172f130e",
 "site_ids":len(by),"complete_sites":len(complete),
 "raw_GBIF_index_record_count_not_population_or_visitor_rate":True,
 "new_observations_or_ecological_regressions":False,
 "per_taxon":summ,"joint_support_diagnostic_NOT_redefined_primary_gate":joint,
 "original_primary_support_threshold":"Plantae>=50, Lycaenidae>=20 on same 5-km 2017–21 site",
 "same_GBIF_filters_all_taxa":True,
 "interpretation_if_all_butterfly_zeros":"A zero official indexed record count at 5km is not proof of no butterflies or ecological benefit."
}
assert out["complete_sites"]==result["sites_all_three_success"]
assert joint["plantae_ge50_lycaenidae_ge20"]==result["sites_plants50_and_lycae20"]
outpath=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"NRSA_GBIF360_SOURCE_COVERAGE_FAILURE_DIAGNOSTIC_20261112.json"
outpath.write_text(json.dumps(out,indent=2,sort_keys=True))
print("NRSA_GBIF_SOURCE_ONLY_AUDIT_START")
print(json.dumps(out,sort_keys=True))
print("NRSA_GBIF_SOURCE_ONLY_AUDIT_END")
