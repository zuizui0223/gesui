#!/usr/bin/env python3
"""Original Borealis Dataverse SOURCE discovery only, no invented microplastic
individual values, no nutrient flux extrapolation from contaminant particles.

Study Wardlaw 2026 8 municipal WWTPs in Grand River, and related 2025/26
lower-trophic-level study. Discover genuine DOI/datafile identifiers and
check source schema before an actual downstream upstream exposure model.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,time,urllib.error,urllib.parse,urllib.request,zipfile
from pathlib import Path
BASE="https://borealisdata.ca/api"
QUERIES=[
 "Replication data for: Cross-ecosystem transfer of microplastics via emergent aquatic insects to riparian spiders",
 "Replication Data for: Limited evidence of elevated microplastics in lower trophic levels downstream of wastewater treatment plants"]
KNOWN_DOI="doi:10.5683/SP3/TLDY6W"
CROSS_ECOSYSTEM_DOI="doi:10.5683/SP4/U0MN9S"
PUBLIC_UA="gesui-20261117-borealis-site-year-wwtp-pair-source-audit"
LIMIT=2_500_000
def sha(b):return hashlib.sha256(b).hexdigest()
def get(url,limit=LIMIT):
    if not (url.startswith(BASE+"/")):
        raise RuntimeError("Not a trusted original Borealis Dataverse endpoint")
    req=urllib.request.Request(url,headers={"User-Agent":PUBLIC_UA,"Accept":"application/json,*/*"})
    last=None
    for i in range(3):
        try:
            with urllib.request.urlopen(req,timeout=65) as r:
                b=r.read(limit+1)
                if len(b)>limit:raise RuntimeError("Source file exceeds max original scan")
                return b
        except (urllib.error.URLError,urllib.error.HTTPError,TimeoutError) as e:
            last=e
            if isinstance(e,urllib.error.HTTPError) and e.code not in [429,500,502,503,504]:break
            time.sleep(3+3*i)
    raise RuntimeError("Original Borealis API access blocked: "+str(last))
def js(url):return json.loads(get(url).decode("utf-8"))
def csv_headers(raw,name):
    delim="\t" if name.lower().endswith((".tab",".tsv")) else ","
    r=csv.DictReader(io.StringIO(raw.decode("utf-8-sig","replace")),delimiter=delim)
    cols=list(r.fieldnames or [])
    rowcount=0
    uniq={q:set() for q in cols if re.search(r"site|wwtp|outfall|upstream|downstream|locat|year|month|type|stage|taxon|species|spider|insect|sample",q,re.I)}
    for line in r:
        if not line:continue
        rowcount+=1
        for x in uniq:
            value=line.get(x)
            if value and len(uniq[x])<=5000:uniq[x].add(value.strip())
    return {"file":name,"source_row_count":rowcount,"original_fields":cols[:110],
            "category_distinct_counts":{q:len(v) for q,v in uniq.items()},
            "source_WWTP_location_columns":[x for x in uniq if re.search(r"site|wwtp|outfall|upstream|downstream|location",x,re.I)],
            "no_raw_values_rehosted":True}
def search():
    results=[]
    for q in QUERIES:
        url=BASE+"/search?"+urllib.parse.urlencode({"q":q,"type":"dataset","per_page":30})
        data=js(url)
        hits=data.get("data",{}).get("items",[])
        results.append({"query":q,"source_result_total_count":data.get("data",{}).get("total_count"),
            "source_entries":[{"id":x.get("global_id"),"title":x.get("name"),"url":x.get("url"),
              "publicationDate":x.get("published_at"),"entity_type":x.get("type")}
                 for x in hits]})
    return results
def original_data(doi):
    url=BASE+"/datasets/:persistentId/?"+urllib.parse.urlencode({"persistentId":doi})
    j=js(url)
    if j.get("status")!="OK":raise RuntimeError("Original DOI metadata is not published OK")
    v=j["data"].get("latestVersion",{})
    files=v.get("files",[])
    details=[]
    for f in files:
        d=f.get("dataFile",{})
        label=f.get("label") or d.get("filename") or ""
        size=int(d.get("filesize") or 0)
        fileid=d.get("id")
        data={"label":label,"source_file_id":fileid,"published_file_size":size,
              "contentType":d.get("contentType"),"restricted":d.get("restricted"),
              "checksum":d.get("md5") or d.get("checksum",{})}
        if fileid and size>0 and size<=LIMIT and label.lower().endswith((".csv",".tsv",".tab")):
            try:
                raw=get(BASE+"/access/datafile/"+str(fileid),max(LIMIT,size+100))
                data["source_bytes_sha256"]=sha(raw)
                data["verified_source_schema"]=csv_headers(raw,label)
            except Exception as ex:
                data["download_status"]="PUBLIC_SOURCE_ACCESS_STOP"
                data["error"]=str(ex)[:220]
        details.append(data)
    return {"original_DOI":doi,"source_title":v.get("metadataBlocks",{}).get("citation",{}).get("fields",[])[0].get("value") if v.get("metadataBlocks",{}).get("citation",{}).get("fields") else "",
      "source_version_number":v.get("versionNumber"),
      "published_file_n":len(files),"source_file_entries":details,
      "source_wastewater_impact_study_no_raw_ecology_model":True}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    result={"schema":"GESUI_20261117_BOREALIS_ORIGINAL_8_MUNICIPAL_WWTP_AQUATIC_SPIDER_DATALINK_SOURCE_V1",
      "source":"https://borealisdata.ca","known_2025_26_lower_trophic_data_DOI":KNOWN_DOI,
      "source_8_WWTP_study":"https://doi.org/10.1016/j.envpol.2026.129116"}
    for name,fn in (("published_query_results",search),("original_8_WWTP_basal_insect_DOI_source",lambda:original_data(KNOWN_DOI)),
                    ("original_8_WWTP_aquatic_adult_riparian_spider_DOI_source",lambda:original_data(CROSS_ECOSYSTEM_DOI))):
        try:result[name]=fn()
        except Exception as ex:
            result[name]={"status":"STOP_REMOTE_PUBLIC_SOURCE_ACCESS_OR_SCHEMA",
                    "error_type":type(ex).__name__,"message":str(ex)[:350]}
    result["no_terrestrial_nutrient_flux_or_cross_ecosystem_effect_fitted"]=True
    result["same_facilities_across_two_DOIs_NOT_YET_PROVEN"]=True
    result["source_specific_contaminant_transfer_NOT_plant_biomass_or_nutrient_flux"]=True
    (a.out/"BOREALIS_WARDLAW_WWTP8_AQUATIC_INSECT_SPIDER_SOURCE_METADATA.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("BOREALIS_8WWTP_SOURCE="+json.dumps({
        "search":result["published_query_results"].get("status","SEARCH_RESPONSE") if isinstance(result["published_query_results"],dict) else "SEARCH_RETURNED_SOURCE_RECORDS",
        "source_2025_26_DOI":result["original_8_WWTP_basal_insect_DOI_source"].get("original_DOI",
                                  result["original_8_WWTP_basal_insect_DOI_source"].get("status")),
        "files":result["original_8_WWTP_basal_insect_DOI_source"].get("published_file_n"),
        "crossystem_original_DOI":result["original_8_WWTP_aquatic_adult_riparian_spider_DOI_source"].get("original_DOI",result["original_8_WWTP_aquatic_adult_riparian_spider_DOI_source"].get("status")),
        "crossystem_files":result["original_8_WWTP_aquatic_adult_riparian_spider_DOI_source"].get("published_file_n")},ensure_ascii=False),flush=True)
if __name__=="__main__":main()
