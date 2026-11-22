#!/usr/bin/env python3
"""Original Loughborough/figshare four-country fine sediment source metadata gate.

Does NOT download or inspect any biological outcome data. Publisher API only.
Any unknown file, license or field is retained as unknown, never guessed.
"""
from __future__ import annotations
import argparse,datetime,hashlib,json,pathlib,urllib.request,urllib.error,time

ARTICLE=24615570
ENDPOINTS=[
 f"https://api.figshare.com/v2/articles/{ARTICLE}",
 f"https://api.figshare.com/v2/articles/{ARTICLE}/files",
 f"https://api.figshare.com/v2/articles/{ARTICLE}/versions"]
def get(url):
 req=urllib.request.Request(url,headers={
  "User-Agent":"gesui-external-science-source-gate/20261122 (research dataset metadata)",
  "Accept":"application/json"})
 for attempt in range(3):
  try:
   with urllib.request.urlopen(req,timeout=45) as r:
    content=r.read(1500000)
    return {"status":"HTTP_SUCCESS","http_status":r.status,
       "response_sha256":hashlib.sha256(content).hexdigest(),
       "body":json.loads(content)}
  except urllib.error.HTTPError as e:
   return {"status":"HTTP_ERROR","http_status":e.code,"reason":str(e)[:300]}
  except Exception as e:
   if attempt==2:return {"status":"NETWORK_ERROR","reason":type(e).__name__+":"+str(e)[:300]}
   time.sleep(2**attempt)
 raise RuntimeError("unreachable")
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--out",required=True,type=pathlib.Path)
 args=p.parse_args()
 api={}
 for url in ENDPOINTS:api[url]=get(url)
 main=api[ENDPOINTS[0]]
 body=main.get("body",{}) if isinstance(main.get("body"),dict) else {}
 files=body.get("files",[]) if isinstance(body.get("files"),list) else []
 fetched=api[ENDPOINTS[1]].get("body",[])
 if isinstance(fetched,list) and len(fetched)>len(files):files=fetched
 names=[]
 for obj in files:
  if not isinstance(obj,dict):continue
  names.append({k:obj.get(k) for k in (
   "id","name","size","download_url","supplied_md5","is_link_only","is_embargoed","mimetype")
   if k in obj})
 provider=body.get("license")
 if isinstance(provider,dict):
  provider={k:provider.get(k) for k in ("name","url","value") if k in provider}
 output={
  "schema":"gesui_external_fourcountry_6491_sediment_EPT_Loughborough_source_metadata_v1",
  "source_study_title":"Freshwater invertebrate responses to fine sediment stress: A multi-continent perspective",
  "source_publication_DOI":"10.1111/gcb.17084",
  "claimed_original_repository_DOI":"10.17028/rd.lboro.24615570",
  "article_id_probed":ARTICLE,
  "probed_UTC":datetime.datetime.now(datetime.timezone.utc).isoformat(),
  "metadata_endpoints":[{"url":u,"status":r["status"],"http_status":r.get("http_status"),
    "reason":r.get("reason"),"response_sha256":r.get("response_sha256")} for u,r in api.items()],
  "publisher_article":{
   "id":body.get("id"),"doi":body.get("doi"),"title":body.get("title"),
   "published_date":body.get("published_date"),"version":body.get("version"),
   "url_public_html":body.get("url_public_html"),
   "defined_type_name":body.get("defined_type_name"),
   "license":provider,
   "references":body.get("references")},
  "publisher_file_list":names,
  "source_only_discovery_GO":bool(main.get("status")=="HTTP_SUCCESS" and files),
  "national_study_published_sample_counts_reference_only":{
    "Australia":3099,"Brazil":232,"UK":2818,"New Zealand":342},
  "sample_level_taxonomic_abundance_confirmed":False,
  "same_site_measured_deposited_fine_sediment_confirmed":False,
  "individual_site_source_ID_confirmed":False,
  "quantitative_outcomes_read":False,
  "numerical_EPT_replication_permitted":False,
  "reason":"Metadata discovery only; a separate source schema audit must verify individual files, study sampling grain and actual taxon count columns."}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(output,sort_keys=True,indent=2)+"\n")
 print("EXTERNAL_FOURCOUNTRY_SOURCE_META",json.dumps({
   "discovery_GO":output["source_only_discovery_GO"],
   "publisher_article":output["publisher_article"],
   "files":names,"endpoint_receipts":output["metadata_endpoints"]
  },sort_keys=True),flush=True)
if __name__=="__main__":main()
