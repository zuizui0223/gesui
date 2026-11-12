#!/usr/bin/env python3
"""HJA AS015 v2 2022 emergence DATA SOURCE ONLY — pre-outcome gate.

Official EDI PASTA endpoint only; original MD5 MUST match AS015 published
metadata. Freezes 17 stream x interval x taxon keys after checking provenance,
before ecological response quantities are compared.
"""
import csv,hashlib,io,json,math,os,pathlib,re,time,urllib.request
PKG="https://pasta.lternet.edu/package"
SCOPE="knb-lter-and/12548/2"
EXPECTED_MD5="f3e6e400990da883ae0d341731431749"
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_AS015_emergence_original_source"
ROOT.mkdir(parents=True,exist_ok=True)
def fetch(url):
    err=None
    for n in range(4):
        try:
            req=urllib.request.Request(url,headers={"Accept":"text/plain,text/csv,application/octet-stream,*/*",
                "User-Agent":"gesui-academic-lter-as015-20261112/1.0"})
            with urllib.request.urlopen(req,timeout=90) as response:
                buf=response.read()
                if len(buf)==0:raise ValueError("Empty response")
                return buf,response.geturl()
        except Exception as e:
            err=e
            print("EDI source request retry",n+1,type(e).__name__,str(e)[:160],flush=True)
            time.sleep(min(1+3*n,9))
    raise RuntimeError("Official EDI original endpoint unavailable "+url+" "+str(err))
name_url=PKG+"/name/eml/"+SCOPE
names,actual_names_url=fetch(name_url)
txt=names.decode("utf-8-sig","replace")
print("Original EDI entity-name response metadata "+repr(txt[:400]),flush=True)
found=re.findall(r"\b[a-fA-F0-9]{32,40}\b",txt)
if len(set(found))!=1 or "AS01501" not in txt:
    raise RuntimeError("Could not uniquely bind AS01501_v1.csv via original PASTA entity names")
entity=found[0]
url=PKG+"/data/eml/"+SCOPE+"/"+entity
data,actual_data_url=fetch(url)
md5=hashlib.md5(data).hexdigest();sha=hashlib.sha256(data).hexdigest()
if md5!=EXPECTED_MD5:
    raise RuntimeError("AS015 ORIGINAL MD5 mismatch; do not substitute or alter source. expected="+EXPECTED_MD5+" actual="+md5+" prefix="+repr(data[:100]))
csv_path=ROOT/"AS01501_v1_ORIGINAL.csv"
csv_path.write_bytes(data)
text=data.decode("utf-8-sig")
rdr=csv.DictReader(io.StringIO(text,newline=""))
head=rdr.fieldnames or []
expected={"DBCODE","ENTITY","SITE","STREAM_REGIME","SAMPLE_DATE","TIME_PERIOD",
          "TAXA_ORDER","EMERG_ABUNDANCE","EMERG_BIOMASS"}
if not expected.issubset(set(head)):raise RuntimeError("AS015 v2 source schema mismatch: "+str(head))
raw=list(rdr)
if len(raw)!=1564:raise RuntimeError("AS015 published source 1564 rows changed "+str(len(raw)))
def k(row):return (row["SITE"].strip(),row["TIME_PERIOD"].strip(),row["SAMPLE_DATE"].strip(),row["TAXA_ORDER"].strip())
keys=[k(row) for row in raw]
if len(set(keys))!=len(keys):raise RuntimeError("Duplicated original AS015 primary key")
orders={"Ephemeroptera","Plecoptera","Trichoptera","Diptera"}
sites={r["SITE"].strip() for r in raw}
regimes={}
groups={}
numeric_bad={"EMERG_ABUNDANCE":0,"EMERG_BIOMASS":0}
dates=[]
for row in raw:
    site=row["SITE"].strip();regime=row["STREAM_REGIME"].strip()
    if site in regimes and regimes[site]!=regime:raise RuntimeError("Stream regime conflicts in original source")
    regimes[site]=regime
    order=row["TAXA_ORDER"].strip()
    if order not in orders:raise RuntimeError("Unexpected original taxonomic unit "+order)
    p=row["TIME_PERIOD"].strip();date=row["SAMPLE_DATE"].strip()
    dates.append(date)
    g=groups.setdefault((site,p),{"dates":set(),"orders":set()})
    g["dates"].add(date)
    g["orders"].add(order)
    for numeric in numeric_bad:
        value=row[numeric].strip()
        try:
            z=float(value)
            if not math.isfinite(z) or z<0:numeric_bad[numeric]+=1
        except (ValueError,TypeError):numeric_bad[numeric]+=1
complete={g:(len(c["orders"])==4 and len(c["dates"])==1) for g,c in groups.items()}
counts={s:len([g for g in complete if g[0]==s and complete[g]]) for s in sites}
regimecounts={reg:sum(x==reg for x in regimes.values()) for reg in ("I","P")}
gate=bool(len(sites)>=15 and all(v>=18 for v in counts.values())
          and all(v>=6 for v in regimecounts.values())
          and sum(complete.values())>=300
          and all(v<=1564*.05 for v in numeric_bad.values()))
# Source-only output excludes numeric measures. It MUST include source exact
# accession, counts, dates and immutable source key hashing.
manifest={
  "schema":"HJA_AS015_v2_EMERGENCE_PRE_OUTCOME_SAMPLE_MANIFEST_v1",
  "EDI_package":"knb-lter-and.12548.2","source_entity":entity,
  "original_url":url,"original_MD5":md5,"original_SHA256":sha,
  "scope":"Oregon HJ Andrews 2022 17 sites and 23 time periods",
  "source_complete_keys":sorted([list(key) for key in keys if complete[(key[0],key[1])]]),
  "no_emergence_numerical_outcomes_in_manifest":True}
manifest_raw=(json.dumps(manifest,sort_keys=True,separators=(",",":"))+"\n").encode()
msha=hashlib.sha256(manifest_raw).hexdigest()
(ROOT/"AS015_V2_SOURCE_ONLY_MANIFEST.json").write_bytes(manifest_raw)
result={
 "status":"PASS_FREEZE_BEFORE_EMERGENCE_AMOUNTS" if gate else "STOP_NO_BIOMASS_ASSOCIATIONS",
 "original_package":"knb-lter-and.12548.2",
 "original_entity_name_response_url":name_url,
 "original_data_entity_id":entity,
 "original_source_url":url,
 "original_resolved_data_url":actual_data_url,
 "source_bytes":len(data),"source_md5":md5,"source_sha256":sha,
 "original_fields":head,"rows":len(raw),
 "independent_streams":len(sites),"stream_regime_count":regimecounts,
 "site_event_count":len(groups),"complete_4_order_site_period_events":sum(complete.values()),
 "periods_per_site_min":min(counts.values()),"periods_per_site_max":max(counts.values()),
 "taxa_orders":sorted(orders),
 "numeric_invalid_or_missing":numeric_bad,
 "original_date_min":min(dates),"original_date_max":max(dates),
 "source_frozen_key_manifest_SHA256":msha,
 "predeclared_min_streams":15,"predeclared_min_complete_events":300,
 "predeclared_min_periods_per_site":18,
 "numeric_emergence_amounts_not_analyzed":True,
 "n_exported_streams_to_NRSA_EPA_sites":0}
(ROOT/"AS015_V2_SOURCE_GATE_SUMMARY.json").write_text(json.dumps(result,indent=2,sort_keys=True))
print("AS015_SOURCE_GATE_RESULT_START",flush=True)
print(json.dumps(result,sort_keys=True),flush=True)
print("AS015_SOURCE_GATE_RESULT_END",flush=True)
