#!/usr/bin/env python3
"""SOURCE-only EPA NRSA2018/19 riparian habitat native-UID match.

No new habitat values printed/fitted. Do not change source support thresholds.
Inputs: immutable BERW+TN/TP UID manifest from run 43034107486.
"""
import collections,csv,hashlib,json,os,pathlib,time,urllib.request
URL="https://www.epa.gov/sites/default/files/2021-04/nrsa_1819_physical_habitat_larger_set_of_metrics_-_data.csv"
SOURCE_SHA="fbf53b8379d4e6564dbb51b3705b02fdd03e4a2caacad06a2314d643c8412e91"
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"nrsa1819_riparian_source"
ROOT.mkdir(exist_ok=True)
src=pathlib.Path(os.environ["NRSA_FROZEN_MANIFEST"]).read_bytes()
if hashlib.sha256(src).hexdigest()!=SOURCE_SHA:raise RuntimeError("Frozen original river site source changed")
m=json.loads(src)["primary_berw_UID_subset"]
M={x["uid"]:x for x in m}
if len(M)!=1221 or len(m)!=1221:raise RuntimeError("Frozen original 1221 UID visits changed")
outpath=ROOT/"EPA_NRSA1819_PHAB_ORIGINAL.csv"
for i in range(4):
    try:
        req=urllib.request.Request(URL,headers={"Accept":"text/csv,*/*","User-Agent":"gesui-riparian-original-physical-habitat-source-audit/1"})
        with urllib.request.urlopen(req,timeout=90) as response,outpath.open("wb") as f:
            while chunk:=response.read(1048576):f.write(chunk)
        if outpath.stat().st_size<1000:raise RuntimeError("Received empty or invalid physical-habitat source")
        break
    except Exception as e:
        outpath.unlink(missing_ok=True)
        if i==3:raise
        time.sleep(1+2*i)
source_hash=hashlib.sha256(outpath.read_bytes()).hexdigest()
with outpath.open(encoding="utf-8-sig",newline="") as f:
    rdr=csv.DictReader(f)
    head=[k.upper().strip() for k in (rdr.fieldnames or [])]
    if "UID" not in head:raise RuntimeError("No unique site-visit UID in original EPA PHab data")
    # Date, sample types, and three physical habitat metrics must be verified
    # from source headers, NOT inferred after seeing numeric measurements.
    fixed=["XGVEG","XPCAN","XFC_OHV","W1H_PIPE","W1_HALL"]
    if any(k not in head for k in fixed):raise RuntimeError("Expected original EPA PHab columns missing: "+str([k for k in fixed if k not in head]))
    all_rows=0;raw_types=collections.Counter();matches=collections.defaultdict(list)
    def val(row,k):
        return str(row.get(k) or "").strip()
    for raw in rdr:
        all_rows+=1
        r={str(k).upper().strip():v for k,v in raw.items() if k is not None}
        typ=val(r,"SAMPLE_TYPE")
        raw_types[typ]+=1
        uid=val(r,"UID")
        if uid not in M:continue
        # Never read habitat magnitudes into source-only output; just
        # support and original independent ID consistency.
        if val(r,"SITE_ID") and val(r,"SITE_ID")!=M[uid]["site_id"]:
            matches[uid].append({"id_conflict":True,"sample_type":typ})
            continue
        fields={}
        for metric in fixed:
            v=val(r,metric)
            fields[metric]=bool(v and v.upper() not in {"NA","N/A","NULL","NAN","MISSING"})
        matches[uid].append({"id_conflict":False,"sample_type":typ,"fields_present":fields})
unique={u:r[0] for u,r in matches.items() if len(r)==1 and not r[0].get("id_conflict")}
joint={u:z for u,z in unique.items() if all(z["fields_present"][c] for c in ["XGVEG","XPCAN","XFC_OHV"])}
states={M[u]["state"] for u in joint}
eco={str(M[u]["eco9"]) for u in joint}
gate=len(joint)>=650 and len(states)>=25 and len(eco)>=6
# All reported qualifying habitat UIDs are unique by construction; duplicate
# original PHab UIDs are omitted; documented below, never chosen by outcome.
result={
 "schema":"NRSA1819_SAME_SITE_RIPARIAN_PHAB_SOURCE_ONLY_v1",
 "not_ecological_or_trophic_result":True,
 "raw_EPA_original_download_url":URL,
 "source_sha256":source_hash,"source_bytes":outpath.stat().st_size,
 "source_columns":head,
 "source_total_rows":all_rows,"source_sample_type_counts":dict(raw_types),
 "frozen_2018_19_river_UID_manifest_sha256":SOURCE_SHA,
 "original_berw_and_tn_tp_source_uid_n":len(M),
 "source_match_UID_keys_any":len(matches),
 "source_match_UID_duplicated":sum(len(z)>1 for z in matches.values()),
 "source_match_UID_id_conflicts":sum(any(z.get("id_conflict") for z in group) for group in matches.values()),
 "source_match_UID_unique_valid":len(unique),
 "source_column_valid_availability":{k:sum(x["fields_present"][k] for x in unique.values()) for k in fixed},
 "source_match_UID_all_3_primary_habitat_nonblank":len(joint),
 "source_match_distinct_states":len(states),
 "source_match_distinct_eco9":len(eco),
 "predeclared_gate_min_uid":650,"predeclared_gate_min_states":25,"predeclared_gate_min_ecoregions":6,
 "gate_pass":gate,
 "decision":"PASS_TO_FROZEN_RIPARIAN_SOURCE_ONLY_UIDS" if gate else "STOP_NO_RIPARIAN_FUNCTIONAL_MODELS"
}
source={
 "schema":"NRSA2018_19_RIPARIAN_PHAB_SOURCE_UIDS_PRE_NUMERIC_v1",
 "source_url":URL,"original_EPA_PHAB_SHA256":source_hash,
 "frozen_berw_river_manifest_SHA256":SOURCE_SHA,
 "eligible_UIDs":sorted(joint),
 "source_sample_types_original":sorted(set(x["sample_type"] for x in joint.values())),
 "no_habitat_numeric_values_examined":True
}
bytestr=(json.dumps(source,sort_keys=True,separators=(",",":"))+"\n").encode()
result["pre_numeric_riparian_UID_manifest_SHA256"]=hashlib.sha256(bytestr).hexdigest()
(ROOT/"NRSA1819_RIPARIAN_UID_SOURCE_ONLY_MANIFEST.json").write_bytes(bytestr)
(ROOT/"NRSA1819_RIPARIAN_QA_SOURCE_GATE_REPORT.json").write_text(json.dumps(result,indent=2))
print("NRSA_RIPARIAN_SOURCE_ONLY_RESULT_START")
print(json.dumps({k:v for k,v in result.items() if k!="source_columns"},sort_keys=True))
print("NRSA_RIPARIAN_SOURCE_ONLY_RESULT_END")
