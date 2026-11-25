#!/usr/bin/env python3
"""Retrieve the official Zenodo deposit of the identical Dryad DOI. Source only."""
import argparse,hashlib,json,pathlib,urllib.parse,zipfile,io
import requests
ID=7316029
DOI="10.5061/dryad.wm37pvmq1"
PINS={
 "1._Regional_species_lists_of_European_ants.csv":(427349,"49c6bf01f4c5c35a859f3032729ee334b9ee3b8437da053dbd550270caaa99d5","efa00191a912302d72e6a186e9972acc"),
 "2._Binary_ranges_maps_of_European__ants.zip":(3550094,"db79bb9dcb80e681e5a6a821d9686157117a36430fd0e288f81944443fa93a8c","7bc541d14f9af8b2438079be5ef70b25"),
 "README.txt":(5430,"e5859087b64492127b0a046a0b573ba8a6b050f4e7515a4544a44daabe50bcbc","7f30b5967f1cedcf865f0002577415fe")}
def sha(b):return hashlib.sha256(b).hexdigest()
def get(u,limit):
    assert urllib.parse.urlparse(u).hostname=="zenodo.org" and u.startswith("https://"),u
    with requests.get(u,timeout=(20,90),stream=True,allow_redirects=False,headers={"User-Agent":"Academic ant distribution source integrity audit"}) as r:
        r.raise_for_status();assert r.status_code==200,(u,r.status_code)
        b=bytearray()
        for q in r.iter_content(1<<18):
            b.extend(q);assert len(b)<=limit
        return bytes(b)
def main(out):
    out.mkdir(parents=True,exist_ok=True)
    r={"status":"SOURCE_STOP","doi":DOI,"record_id":ID,"received":{},"no_ecological_outcomes":True}
    try:
        b=get("https://zenodo.org/api/records/"+str(ID),4*1024**2)
        (out/"zenodo_record.json").write_bytes(b);j=json.loads(b)
        assert j["id"]==ID and j["metadata"]["doi"].lower()==DOI
        assert set(x["key"] for x in j["files"])==set(PINS)
        r["metadata_sha256"]=sha(b)
        r["license"]=j["metadata"].get("license");r["title"]=j["metadata"]["title"]
        for f in j["files"]:
            name=f["key"];n,h,md=PINS[name]
            assert f["size"]==n and f["checksum"]=="md5:"+md,(name,f)
            u=f["links"].get("content") or f["links"]["self"]
            b=get(u,n+1)
            assert len(b)==n and sha(b)==h and hashlib.md5(b).hexdigest()==md,("not identical to Dryad",name,len(b),sha(b))
            (out/name).write_bytes(b)
            rec={"bytes":n,"sha256":h,"md5":md,"url":u}
            if name.endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(b)) as z:
                    assert z.testzip() is None
                    rec["members"]=[{"name":a.filename,"bytes":a.file_size,"sha256":sha(z.read(a))} for a in z.infolist() if not a.is_dir()]
            elif name.endswith(".csv"):
                rec["first_rows"]=b.decode("utf-8-sig",errors="strict").splitlines()[:3]
            else:
                rec["readme"]=b.decode("utf-8-sig",errors="strict")
            r["received"][name]=rec
        r["status"]="IDENTICAL_DRYAD_ZENODO_SOURCE_PASS"
    except Exception as e:r["error"]=type(e).__name__+": "+str(e)[:1500]
    finally:
        (out/"receipt.json").write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
        print("EUAD SOURCE "+json.dumps({**r,"received":{k:{a:b for a,b in v.items() if a!="members"}|{"members_n":len(v.get("members",[]))} for k,v in r["received"].items()}},ensure_ascii=False),flush=True)
        assert r["status"]=="IDENTICAL_DRYAD_ZENODO_SOURCE_PASS",r["status"]
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args().out)
