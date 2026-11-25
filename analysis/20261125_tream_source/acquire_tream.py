#!/usr/bin/env python3
"""Public TREAM source-only census. No numerical abundance or biodiversity outcomes evaluated."""
import argparse, collections, csv, hashlib, json, pathlib, re, tempfile, urllib.parse, xml.etree.ElementTree as ET
import requests
import traceback
DOI="doi:10.5063/F1NG4P4R"
CN="https://cn.dataone.org/cn/v2/"
MN="https://knb.ecoinformatics.org/knb/d1/mn/v2/"
NAMES=("TREAM_siteLevel.csv","TREAM_siteYEARLevel.csv","TREAM_alltaxa.csv")
WINDOWS=((2007,2009),(2010,2012),(2013,2015),(2016,2018))
DUP={"100000001","100000309","100000002","100000308"}
HOSTS={"knb.ecoinformatics.org","cn.dataone.org","dataone.org"}
def norm(s):return re.sub(r"[^a-z0-9]","",str(s).lower())
def lname(s):return s.rsplit("}",1)[-1]
def nodes(root,name):return [x for x in root.iter() if lname(x.tag)==name]
def texts(root,name):return ["".join(x.itertext()).strip() for x in nodes(root,name)]
def sha(b):return hashlib.sha256(b).hexdigest()
def obj_url(base,pid):return base+"object/"+urllib.parse.quote(pid,safe="")
def get(url,limit=600_000_000):
    assert urllib.parse.urlparse(url).scheme=="https" and urllib.parse.urlparse(url).hostname in HOSTS,url
    with requests.get(url,timeout=(30,180),stream=True,headers={"User-Agent":"Academic TREAM data provenance audit"}) as r:
        r.raise_for_status()
        if urllib.parse.urlparse(r.url).hostname not in HOSTS:raise RuntimeError("Unreviewed redirect host "+str(r.url))
        data=bytearray()
        for block in r.iter_content(1<<20):
            data.extend(block)
            if len(data)>limit:raise RuntimeError("Byte limit exceeded")
        return bytes(data),{"requested_url":url,"resolved_url":r.url,"content_type":r.headers.get("Content-Type"),"bytes":len(data),"sha256":sha(data)}
def entities(x):
    # EML can represent CSVs as dataTable or otherEntity.
    return nodes(x,"dataTable")+nodes(x,"otherEntity")
def validate_eml(b):
    x=ET.fromstring(b)
    if lname(x.tag).lower()!="eml":raise ValueError("Object is not EML: root="+x.tag)
    found={name for t in entities(x) for name in NAMES
           if name in texts(t,"objectName")+texts(t,"entityName")}
    if found!=set(NAMES):
        raise ValueError("Exact published objects not found: "+repr(sorted(found))+
                         "; entity names="+repr([(texts(t,"objectName"),texts(t,"entityName")) for t in entities(x)])[:3000])
    return x
def resolve_eml(out):
    errors=[]
    def attempt(u,label):
        b,rec=get(u,10_000_000)
        (out/("metadata_candidate_"+label+".xml")).write_bytes(b)
        rec["attempt_label"]=label
        x=validate_eml(b)
        return x,b,rec,errors
    for i,u in enumerate([obj_url(MN,DOI),obj_url(CN,DOI)]):
        try:return attempt(u,str(i))
        except Exception as e:errors.append({"url":u,"error":type(e).__name__+": "+str(e)[:1500]})
    q='(id:"'+DOI+'" OR seriesId:"'+DOI+'") AND formatType:METADATA'
    url=CN+"query/solr/?"+urllib.parse.urlencode({"q":q,"fl":"id,seriesId,title,formatType","rows":20,"wt":"json"})
    try:
        b,rec=get(url,5_000_000);docs=json.loads(b)["response"]["docs"]
        (out/"doi_index_response.json").write_bytes(b)
        ids=sorted({d["id"] for d in docs if (d.get("id")==DOI or d.get("seriesId")==DOI) and d.get("formatType")=="METADATA"})
        for i,pid in enumerate(ids):
            try:return attempt(obj_url(CN,pid),"index_"+str(i))
            except Exception as e:errors.append({"id":pid,"error":type(e).__name__+": "+str(e)[:1500]})
    except Exception as e:errors.append({"url":url,"error":type(e).__name__+": "+str(e)[:1500]})
    (out/"metadata_resolution_errors.json").write_text(json.dumps(errors,indent=2))
    raise RuntimeError("Exact DOI EML resolution failed; see saved diagnostics "+repr(errors)[:2500])
def tables(eml):
    results={}
    for table in entities(eml):
        matches=[name for name in NAMES if name in texts(table,"objectName") or name in texts(table,"entityName")]
        if not matches:continue
        assert len(matches)==1 and matches[0] not in results, ("Ambiguous metadata entity",matches)
        attrs=[]
        for a in nodes(table,"attribute"):
            attrs.append({"name":texts(a,"attributeName"),"definition":texts(a,"attributeDefinition"),"missing_codes":texts(a,"code"),"unit":texts(a,"standardUnit")+texts(a,"customUnit")})
        online=[]
        for item in nodes(table,"online"):
            for u in texts(item,"url"):
                if u.startswith("https://") and urllib.parse.urlparse(u).hostname in HOSTS:online.append(u)
                elif u.startswith("http://") and urllib.parse.urlparse(u).hostname in HOSTS:online.append("https://"+u[7:])
        checks=[]
        for a in nodes(table,"authentication"):
            checks.append({"method":a.attrib.get("method",""),"value":"".join(a.itertext()).strip()})
        results[matches[0]]={"urls":list(dict.fromkeys(online)),"attributes":attrs,"checksums":checks,"entity_name":texts(table,"entityName")}
    if set(results)!=set(NAMES):raise RuntimeError("EML table names not exact: "+repr([(texts(t,"objectName"),texts(t,"entityName")) for t in entities(eml)]))
    return results
def receive(name,info,p):
    errors=[]
    for u in info["urls"]:
        try:
            b,rec=get(u)
            if b.lstrip().startswith(b"<"):raise ValueError("XML/HTML rather than CSV")
            for a in info["checksums"]:
                key=norm(a["method"])
                if key in {"md5","sha1","sha256"}:
                    observed=hashlib.new(key,b).hexdigest()
                    assert observed.lower()==a["value"].lower(),("Publisher checksum mismatch",name,key)
            p.write_bytes(b);rec["publisher_checksums"]=info["checksums"];return rec
        except Exception as e:errors.append([u,str(e)[:400]])
    raise RuntimeError("Public data object unavailable "+name+" "+repr(errors))
def pick(headers,aliases,required=True):
    h=[x for x in headers if norm(x) in aliases]
    if len(h)!=1:
        if not required and not h:return None
        raise RuntimeError("SCHEMA_REVIEW_REQUIRED "+repr(sorted(aliases))+" matches "+repr(h))
    return h[0]
def rows(p):
    with p.open("r",encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f)
        assert r.fieldnames and len(set(r.fieldnames))==len(r.fieldnames)
        yield r.fieldnames,r
def clean(v):return str(v or "").strip()
def sitekey(v):return clean(v)
def year(v):
    n=float(clean(v));assert n.is_integer() and 1900<=n<=2100
    return int(n)
def census(work,result):
    p={n:work/n for n in NAMES}
    headers={}
    for n in NAMES:
        for h,r in rows(p[n]):headers[n]=h
    result["headers"]=headers
    result["stage"]="schema"
    id_alias={"siteid","sitecode","sitenumber","siteidentifier","idsite"}
    year_alias={"year","samplingyear","sampleyear"}
    site={}
    for h,r in rows(p[NAMES[0]]):
        sid=pick(h,id_alias);co=pick(h,{"country","countryname","countrycode","countryiso","countryiso2","countryiso3"})
        pr=pick(h,{"projectid","project","projectnumber","projectcode","datasetid","dataid","dataset","idproject"})
        tax=pick(h,{"taxonomicresolution","primarytaxonomicresolution","taxaresolution","taxresolution","taxlevel","primarytaxonomiclevel","taxonomiclevel","taxonomicres"},False)
        result["schema_site"]={"site":sid,"country":co,"project":pr,"taxonomy":tax}
        for a in r:
            k=sitekey(a[sid]);assert k and k not in site
            site[k]={"country":clean(a[co]),"project":clean(a[pr]),"taxonomy":clean(a[tax]) if tax else None}
    assert len(site)==1816
    years=collections.defaultdict(set);siteyears=set()
    for h,r in rows(p[NAMES[1]]):
        sid=pick(h,id_alias);yr=pick(h,year_alias)
        result["schema_siteyear"]={"site":sid,"year":yr}
        for a in r:
            k=sitekey(a[sid]);y=year(a[yr]);assert k in site and (k,y) not in siteyears
            siteyears.add((k,y));years[k].add(y)
    assert set(years)==set(site)
    tax_years=set();order_missing=0;family_missing=0;nrows=0;unknown=set()
    for h,r in rows(p[NAMES[2]]):
        sid=pick(h,id_alias);yr=pick(h,year_alias)
        fam=pick(h,{"family","taxonfamily","familyname"},False);order=pick(h,{"order","taxonorder","ordername"},False)
        result["schema_taxa"]={"site":sid,"year":yr,"family":fam,"order":order}
        for a in r:
            k=sitekey(a[sid]);y=year(a[yr]);nrows+=1
            if (k,y) not in siteyears:unknown.add((k,y))
            tax_years.add((k,y))
            if fam and clean(a[fam]).lower() in {"","na","nan","null","none"}:family_missing+=1
            if order and clean(a[order]).lower() in {"","na","nan","null","none"}:order_missing+=1
    result["source_census"]={"sites":len(site),"siteyears":len(siteyears),"countries":dict(collections.Counter(x["country"] for x in site.values())),"project_count":len({x["project"] for x in site.values()}),"taxa_rows":nrows,"taxa_siteyears":len(tax_years),"taxa_without_siteyear":len(unknown),"siteyears_without_taxa":len(siteyears-tax_years),"missing_family_rows":family_missing if fam else None,"missing_order_rows":order_missing if order else None,"taxonomy_resolution_counts":dict(collections.Counter(str(x["taxonomy"]) for x in site.values())),"known_duplicate_ids_present":sorted(DUP&set(site))}
    assert not unknown
    # Candidate eligibility concerns source calendar support only, not abundance, chemistry or quality.
    eligible=[]
    for k,v in site.items():
        counts=[sum(lo<=y<=hi and (k,y) in tax_years for y in years[k]) for lo,hi in WINDOWS]
        if min(counts)>=2:eligible.append((k,v,counts))
    def subset(x):
        return {"sites":len(x),"countries":dict(collections.Counter(v["country"] for k,v,c in x)),"projects":len({v["project"] for k,v,c in x}),"site_id_sha256":sha(("\n".join(sorted(k for k,v,c in x))+"\n").encode())}
    excluded={"unitedkingdom","uk","gb","gbr","greatbritain","england","scotland","wales","northernireland","france","fr","fra"}
    primary=[x for x in eligible if x[0] not in DUP and norm(x[1]["country"]) not in excluded and x[1]["country"] and x[1]["project"]]
    result["calendar_opportunity_all"]=subset(eligible);result["calendar_opportunity_conservative"]=subset(primary)
    s=result["calendar_opportunity_conservative"]
    result["calendar_threshold_pass"]=s["sites"]>=100 and len(s["countries"])>=5 and s["projects"]>=10 and sum(v>=20 for v in s["countries"].values())>=3
    result["status"]="SOURCE_CENSUS_COMPLETE_REVIEW_TAXONOMY_AND_ENVIRONMENT_BEFORE_FIT"
    with (out_global/"candidate_roster.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["site_id","country","project","reported_taxonomy","n_2007_09","n_2010_12","n_2013_15","n_2016_18","conservative_candidate"])
        ids={k for k,v,c in primary}
        for k,v,c in sorted(eligible):w.writerow([k,v["country"],v["project"],v["taxonomy"],*c,k in ids])
def main(out):
    global out_global
    out_global=out;out.mkdir(parents=True,exist_ok=True)
    result={"doi":DOI,"windows":WINDOWS,"outcomes_evaluated":False,"status":"SOURCE_STOP","source":{}}
    try:
        x,b,rec,errs=resolve_eml(out);result["eml"]=rec;result["direct_resolution_errors"]=errs
        result["titles"]=texts(x,"title");result["rights"]=texts(x,"intellectualRights")
        result["table_metadata"]=tables(x)
        (out/"metadata_schema.json").write_text(json.dumps(result["table_metadata"],indent=2,ensure_ascii=False))
        with tempfile.TemporaryDirectory() as d:
            w=pathlib.Path(d)
            for n in NAMES:
                result["source"][n]=receive(n,result["table_metadata"][n],w/n)
            census(w,result)
    except Exception as e:
        result["status"]="SOURCE_STOP";result["error"]=type(e).__name__+": "+str(e)[:3000]
        result["traceback"]=traceback.format_exc(limit=6)
        print(result["traceback"],flush=True)
    finally:
        (out/"source_receipt.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
        print("TREAM SOURCE",json.dumps({k:v for k,v in result.items() if k not in {"table_metadata","rights"}},ensure_ascii=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args().out)
