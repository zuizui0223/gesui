#!/usr/bin/env python3
"""SOURCE-ONLY freeze original 234 EU level-7 basins with real benthic EQR in
both 2018 and 2022, and exact accepted four butterfly family GBIF taxon keys.
No numerical GBIF species richness, butterfly counts or EQR outcome is read.
"""
import collections,hashlib,json,os,pathlib,urllib.parse,urllib.request,time

sourcefile=pathlib.Path(os.environ["HYBAS7_ORIGINAL_NATIVE_RIVER_SOURCE"])
wktfile=pathlib.Path(os.environ["HYBAS7_ORIGINAL_NATIVE_WKT_SOURCE"])
SOURCE_SHA="a396eaed1c0b406eb3bbb7178e04ba91989ab5964c4f0175fc1193673d624d89"
WKT_SHA="47f597e0a798cd6155ee9fb69ce8916ed5629951d046d114feb9f667963310e4"
if hashlib.sha256(sourcefile.read_bytes()).hexdigest()!=SOURCE_SHA or \
   hashlib.sha256(wktfile.read_bytes()).hexdigest()!=WKT_SHA:
    raise RuntimeError("Original raw H7 native river/WKT cohort SHA changed")
sites=json.loads(sourcefile.read_text())["candidate_source_site_period_rows"]
native_wkt=json.loads(wktfile.read_text())["basins"]
bytype=collections.defaultdict(set)
for s in sites:
    h=s["HYBAS_ID_L07"]
    if h is None:continue
    if s["panel"] in ("EEA18","EEA22"):bytype[s["panel"]].add(str(h))
repeated=sorted(bytype["EEA18"] & bytype["EEA22"])
if len(repeated)!=234:raise RuntimeError("Original EEA repeated basin source not exactly 234")
poly={z["HYBAS_ID_L07"]:z for z in native_wkt}
if len(poly)!=1219 or not set(repeated).issubset(poly):
    raise RuntimeError("Native HYBAS_ID polygon 1219 source not valid")
if not all(poly[k]["status"]=="all_original_polygons_WKT_queryable" for k in repeated):
    raise RuntimeError("Some original repeated EQR river basins lack complete valid WKT source")
old="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
expected=["Lycaenidae","Nymphalidae","Pieridae","Papilionidae","Hesperiidae"]
keys={}
for label in expected:
    q=urllib.parse.urlencode({"name":label,"rank":"FAMILY","strict":"true","checklistKey":old})
    url="https://api.gbif.org/v1/species/match?"+q
    err=None
    for t in range(3):
        try:
            req=urllib.request.Request(url,headers={"Accept":"application/json","User-Agent":"gesui-H7-butterfly-family-original-taxonomy-source-20261112/1"})
            with urllib.request.urlopen(req,timeout=50) as r:raw=r.read()
            obj=json.loads(raw)
            name=str(obj.get("scientificName") or "").lower().split(" ")[0]
            rank=obj.get("rank")
            key=obj.get("usageKey")
            if name!=label.lower() or rank!="FAMILY" or not isinstance(key,int):
                raise RuntimeError("GBIF legacy FAMILY mismatch "+label+" "+str(obj)[:350])
            # Cross-check the exact accepted species info at a stable numeric ID.
            species_url="https://api.gbif.org/v1/species/"+str(key)
            with urllib.request.urlopen(urllib.request.Request(species_url,headers={"Accept":"application/json","User-Agent":"gesui-original-L7-family-id/1"}),timeout=45) as r:
                sr=r.read()
            tax=json.loads(sr)
            if tax.get("rank")!="FAMILY" or label.lower() not in str(tax.get("scientificName","")).lower():
                raise RuntimeError("GBIF numeric family taxon key does not resolve to correct original family")
            keys[label]={"taxonKey":key,"rank":"FAMILY","scientificName":tax.get("scientificName"),
             "source_match_URL":url,"source_match_SHA256":hashlib.sha256(raw).hexdigest(),
             "source_species_URL":species_url,"source_species_SHA256":hashlib.sha256(sr).hexdigest()}
            print("Verified legacy butterfly family taxonomy",label,key,flush=True)
            break
        except Exception as e:
            err=e;time.sleep(1+t*2)
    else:raise RuntimeError("Original GBIF family identity unavailable, do not guess "+label+" "+str(err))
if keys["Lycaenidae"]["taxonKey"]!=5473:raise RuntimeError("Legacy GBIF Lycaenidae key must stay 5473")
report={
 "schema":"gesui_234_EU_EEA_both_years_H7_basin_pre_new_species_and_butterfly_counts_v1",
 "source_original_H7_station_manifest_sha256":SOURCE_SHA,
 "source_original_full_native_WKT_manifest_sha256":WKT_SHA,
 "cohort_selection":"all 234 EEA original H7 subbasins with >=1 monitoring station in both independent original years",
 "original_repeated_H7_basin_n":len(repeated),
 "periods":["2017-2019","2020-2022"],
 "source_repeated_basins_with_complete_original_polygons":len(repeated),
 "butterfly_family_original_GBIF_Backbone_checklist":old,
 "butterfly_5_families_source_verified":keys,
 "original_biology_outcomes_NOT_READ":True,
 "new_GBIF_species_facets_or_new_butterfly_counts_NOT_QUERIED":True
}
manifest={"schema":report["schema"],"source_station_manifest_SHA256":SOURCE_SHA,
 "source_WKT_manifest_SHA256":WKT_SHA,
 "verified_butterfly_family_taxa":keys,
 "original_two_year_EEA_H7_subbasins":[poly[x] for x in repeated],
 "fixed_periods":["2017-2019","2020-2022"],
 "source_only_sample_selected_before_new_species_butterfly_results":True}
outdir=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_H7_234_EU_macro_cohort"
outdir.mkdir(parents=True,exist_ok=True)
encoded=(json.dumps(manifest,sort_keys=True,separators=(",",":"))+"\n").encode()
report["fixed_234_H7_original_WKT_and_taxonomy_source_manifest_sha256"]=hashlib.sha256(encoded).hexdigest()
(outdir/"HYBAS7_234_REPEATED_EU_EQR_ORIGINAL_WKT_TAXON_FAMILY_PRE_NEW_OUTCOME.json").write_bytes(encoded)
(outdir/"HYBAS7_234_EU_REPEAT_ORIGINAL_STATION_TAXON_SOURCE_QA.json").write_text(json.dumps(report,indent=2,sort_keys=True))
print("HYBAS7_234_EU_PRE_NEW_BIO_RESPONSE_SOURCE_START")
print(json.dumps(report,sort_keys=True))
print("HYBAS7_234_EU_PRE_NEW_BIO_RESPONSE_SOURCE_END")
