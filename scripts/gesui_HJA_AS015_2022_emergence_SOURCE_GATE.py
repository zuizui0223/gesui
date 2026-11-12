#!/usr/bin/env python3
"""HJA AS015 v2 2022 emergence DATA SOURCE ONLY — pre-outcome gate.

Official HJA EML-linked original; EDI PASTA endpoint returned HTTP 403.
Original MD5 MUST match AS015 published metadata. Freezes 17 stream x interval x taxon keys after checking provenance,
before ecological response quantities are compared.
"""
import csv,hashlib,io,json,math,os,pathlib,re,time,urllib.request,urllib.parse,base64,xml.etree.ElementTree as ET
EML_URL="https://andrewsforest.oregonstate.edu/sites/default/files/lter/fsdbdata/eml_archive/AS015_2.xml"
CATALOG_URL="https://andrewsforest.oregonstate.edu/data/datacatalog/AS015"
EXPECTED_EML_SHA256="fb4a6eb41374e7d682bb74b77a0841bd802198f8285951ca7c47a86786fc962e"
EXPECTED_MD5="f3e6e400990da883ae0d341731431749"
ROOT=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"gesui_AS015_emergence_original_source"
ROOT.mkdir(parents=True,exist_ok=True)
def fetch(url, tries=2):
    err=None
    for n in range(tries):
        try:
            req=urllib.request.Request(url,headers={"Accept":"text/plain,text/csv,application/octet-stream,application/xml,text/html,*/*",
                "User-Agent":"Mozilla/5.0 (compatible; gesui-hja-source-integrity/1.0)"})
            with urllib.request.urlopen(req,timeout=60) as response:
                buf=response.read()
                if len(buf)==0:raise ValueError("Empty original source response")
                return buf,response.geturl()
        except Exception as e:
            err=e
            print("Official HJA source request:",type(e).__name__,str(e)[:120],flush=True)
            time.sleep(min(1+n*2,3))
    raise RuntimeError("Original HJA source endpoint unavailable "+url+" "+str(err))
eml,eml_resolved=fetch(EML_URL)
eml_sha=hashlib.sha256(eml).hexdigest()
if eml_sha!=EXPECTED_EML_SHA256:raise RuntimeError("Original HJA EML SHA changed; refuse replacement")
root=ET.fromstring(eml)
if root.attrib.get("packageId")!="knb-lter-and.12548.2":
    raise RuntimeError("Official HJA EML packageId mismatch")
def local(el):return el.tag.rsplit("}",1)[-1]
table_candidates=[]
for t in root.iter():
    if local(t)!="dataTable":continue
    names=[(x.text or "").strip() for x in t.iter() if local(x)=="objectName"]
    urls=[(x.text or "").strip() for x in t.iter() if local(x)=="url"]
    if "AS01501_v1.csv" in names:
        table_candidates.append([u for u in urls if u.endswith("/AS01501_v1.csv")])
if len(table_candidates)!=1 or len(table_candidates[0])!=1:
    raise RuntimeError("Original HJA EML has ambiguous/missing AS015 entity download URL")
url=table_candidates[0][0]
if url!="https://andrewsforest.oregonstate.edu/fsdbdata/cache/AS01501_v1.csv":
    raise RuntimeError("Unexpected original authoritative EML-linked CSV URL")
name_url=EML_URL
entity="AS01501_v1.csv"
access_route="EML_LINK_DIRECT"
try:
    data,actual_data_url=fetch(url)
except RuntimeError:
    # HJA's current official catalog distributes exactly the SAME EML URL using
    # its signed /data/datacatalog/download/<encoded-url> handler. Verify that
    # its base64 path component decodes to the EML-published URL.
    page,_=fetch(CATALOG_URL)
    paths=re.findall(rb'href="(https://andrewsforest\.oregonstate\.edu/data/datacatalog/download/[^"]+)"',page)
    proxy=[]
    for path in paths:
        decoded=path.decode("utf-8")
        token=decoded.rsplit("/",1)[-1].split(".",1)[0]
        try:
            source=base64.urlsafe_b64decode(token+"="*(-len(token)%4)).decode("utf8")
            if source==url:proxy.append(decoded)
        except (ValueError, UnicodeError):pass
    if len(set(proxy))!=1:
        raise RuntimeError("Cannot validate unique official HJA signed download against published EML source URL")
    data,actual_data_url=fetch(proxy[0])
    access_route="EML_VERIFIED_OFFICIAL_CATALOG_SIGNED_DOWNLOAD"
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
 "original_EML_SHA256":eml_sha,
 "official_source_access_route":access_route,
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
