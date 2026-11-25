#!/usr/bin/env python3
import argparse,collections,datetime,hashlib,io,json,pathlib,time,urllib.parse,zipfile
import requests
CHECKLIST="d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
BASE="https://api.gbif.org/v1/"
HIST_SHA="4b37c741d1bd28c178b0cef61c9457af0dd30797c173f9384ebe1d1d0036d07c"
INNER_SHA="251f8514684fe8ea62d94ddc8c6eeefb59976471d84d7389b682aff6637144e4"
PILOT=[
("2070012630","2060012310"),("2070012870","2060012730"),("2070012960","2060012960"),
("2070013100","2060013100"),("2070016510","2060016510"),("2070016680","2060016680"),
("2070020830","2060020740"),("2070021220","2060021220"),("2070023010","2060023010"),
("2070023330","2060023330"),("2070026920","2060026920"),("2070026950","2060026950")]
GENERA="Aphaenogaster Anoplolepis Camponotus Cataglyphis Crematogaster Dolichoderus Formica Froggattella Iberoformica Iridomyrmex Lasius Lepisiota Linepithema Myrmica Nylanderia Pheidole Plagiolepis Tapinoma Temnothorax Tetramorium".split()
def sha(b):return hashlib.sha256(b).hexdigest()
def req(path,params,out,label):
    err=[]
    for attempt in range(3):
        try:
            r=requests.get(BASE+path,params=params,timeout=(15,60),headers={"User-Agent":"Academic butterfly-ant catchment occurrence feasibility audit"})
            if r.status_code in (429,500,502,503,504):
                err.append({"status":r.status_code});time.sleep(3*(attempt+1));continue
            r.raise_for_status()
            assert urllib.parse.urlparse(r.url).hostname=="api.gbif.org"
            b=r.content; j=json.loads(b)
            n=sha(b)+".json";(out/"responses"/n).write_bytes(b)
            return j,{"url":r.url,"sha256":sha(b),"response_file":"responses/"+n,"utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"status":"ok","retry_errors":err}
        except Exception as e:
            err.append({"error":type(e).__name__+": "+str(e)[:1000]})
            time.sleep(2*(attempt+1))
    raise RuntimeError(label+" "+repr(err))
def match(name,rank,family,out):
    try:
        p={"name":name,"rank":rank,"kingdom":"Animalia","strict":"true","checklistKey":CHECKLIST}
        if family:p["family"]="Formicidae"
        a,rec=req("species/match",p,out,name)
        assert a.get("matchType")=="EXACT" and a.get("rank")==rank and a.get("canonicalName")==name,(name,a)
        k=a.get("acceptedUsageKey",a.get("usageKey"));assert isinstance(k,int)
        b,rec2=req("species/"+str(k),{},out,name+" accepted")
        assert b.get("rank")==rank and b.get("taxonomicStatus")=="ACCEPTED",b
        assert b.get("datasetKey")==CHECKLIST,b
        if family:assert b.get("familyKey")==family and b.get("family")=="Formicidae",b
        else:assert b.get("canonicalName")=="Formicidae" and b.get("kingdom")=="Animalia",b
        return {"name":name,"status":"ok","key":k,"accepted_name":b.get("canonicalName"),"match":rec,"usage":rec2}
    except Exception as e:return {"name":name,"status":"UNKNOWN","error":type(e).__name__+": "+str(e)[:2000]}
def geos(path):
    b=path.read_bytes();assert sha(b)==HIST_SHA
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        b=z.read("GESUI_H7_234x2x6_SPECIESKEY_MERGED_SOURCE.json");assert sha(b)==INNER_SHA
    raw=json.loads(b);assert len(raw)==2808
    d=collections.defaultdict(list)
    for a in raw:
        h=str(a["HYBAS_ID_L07"])
        if h not in dict(PILOT):continue
        g=set()
        for v in a["source_request_receipts"]:
            if v.get("status")!="source_success" or not v.get("url"):continue
            q=urllib.parse.parse_qs(urllib.parse.urlparse(v["url"]).query)
            assert q.get("checklistKey",[]) in ([],[CHECKLIST]), (h,q.get("checklistKey"))
            # Legacy archived URLs without checklistKey used the default numeric Backbone.
            # Geometry is unchanged; all new ant queries explicitly specify CHECKLIST.
            g.add(q["geometry"][0])
        assert g,(h,a["group"])
        d[h].append((a["period"],a["group"],tuple(sorted(g))))
    assert len(d)==12
    G={}
    for h,h6 in PILOT:
        assert len(d[h])==12 and len(set((p,g) for p,g,_ in d[h]))==12
        u=set(x[2] for x in d[h]);assert len(u)==1
        parts=next(iter(u))
        G[h]=parts
    return G
def audit_facet(j,name):
    facets=[f for f in j.get("facets",[]) if f.get("field","").upper()==name.upper()]
    assert len(facets)==1,(name,j.keys())
    counts=facets[0]["counts"];assert len(counts)<1000
    d={str(v["name"]):int(v["count"]) for v in counts}
    assert len(d)==len(counts) and all(v>0 for v in d.values())
    assert sum(d.values())<=int(j["count"])
    return d
def main(a):
    out=a.out;out.mkdir(parents=True,exist_ok=True);(out/"responses").mkdir(exist_ok=True)
    G=geos(a.historical);tax=[match("Formicidae","FAMILY",None,out)]
    assert tax[0]["status"]=="ok",tax
    k=tax[0]["key"];tax.extend(match(n,"GENUS",k,out) for n in GENERA)
    (out/"taxonomy.json").write_text(json.dumps(tax,indent=2)+"\n")
    names={str(v["key"]):v["name"] for v in tax[1:] if v["status"]=="ok"}
    assert len(names)==sum(v["status"]=="ok" for v in tax[1:]),"synonym collision"
    results=[]
    for h,h6 in PILOT:
        for years in ("2014,2019","2020,2025"):
            parts=G[h]
            if len(parts)!=1:
                results.append({"h":h,"h6":h6,"years":years,"status":"GEOMETRY_SOURCE_STOP","polygon_components":len(parts),"gate":False})
                continue
            q={"geometry":parts[0],"year":years,"limit":0,"hasCoordinate":"true","hasGeospatialIssue":"false","occurrenceStatus":"PRESENT","taxonKey":k,"checklistKey":CHECKLIST,"facet":["genusKey","speciesKey"],"facetLimit":1000}
            row={"h":h,"h6":h6,"years":years,"geometry_sha256":sha(parts[0].encode())}
            try:
                j,rec=req("occurrence/search",q,out,h+" "+years)
                g=audit_facet(j,"GENUS_KEY");s=audit_facet(j,"SPECIES_KEY")
                count=int(j["count"]);ng=sum(g.values());ns=sum(s.values());cov=ng/count if count else 0.
                detected=[names[i] for i in names if g.get(i,0)>0]
                row.update(status="ok",source=rec,total_records=count,genus_records=ng,species_records=ns,genus_coverage=cov,species_coverage=ns/count if count else 0.,genus_keys=g,species_keys=s,detected_partner_genera=sorted(detected),gate=ng>=20 and cov>=.8 and len(detected)>=2)
            except Exception as e:row.update(status="UNKNOWN",error=type(e).__name__+": "+str(e)[:1200],gate=False)
            results.append(row)
            print("ANT PILOT ROW",json.dumps({k:v for k,v in row.items() if k not in ("source","genus_keys","species_keys")}),flush=True)
            time.sleep(.35)
    eligible=[h for h,_ in PILOT if all(r["gate"] for r in results if r["h"]==h)]
    result={"status":"PILOT_GATE_PASS" if len(eligible)>=9 and len(names)==20 else "PILOT_SOURCE_STOP","sites":len(PILOT),"genus_taxonomy_resolved":len(names),"successful_siteperiods":sum(r["status"]=="ok" for r in results),"eligible_sites":len(eligible),"eligible_ids":eligible,"no_ecological_outcomes":True,"no_true_absence_inference":True,"rows":results}
    (out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("ANT PILOT",json.dumps({k:v for k,v in result.items() if k!="rows"}),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--historical",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args())
