#!/usr/bin/env python3
"""Publisher documentation and exact original chemistry header inventory, no outcome analysis."""
import argparse,csv,hashlib,io,json,pathlib,sys,tempfile,zipfile
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"analysis/20261125_chempop_resolution"))
import chemistry_order as prev
src=prev.src
SUPPORT="https://data-package.ceh.ac.uk/sd/23f194b0-1e21-4f14-be9c-6eabdcc1feb7.zip"
def sha(x):return hashlib.sha256(x).hexdigest()
def main(work,out):
    out.mkdir(parents=True,exist_ok=True)
    r=prev.receipt(work/"chemistry.zip","18c764318c2143a40565b8d58be350d91c7b8e3b371c9705038401db8730ff43")
    expected={v["original_file"]:v["source_member_SHA256"] for v in r["original_wq_source_files"].values()}
    assert len(expected)==21
    result={"status":"SOURCE_REVIEW","support_url":SUPPORT,"no_outcomes_read":True}
    with tempfile.TemporaryDirectory() as tmp:
        d=pathlib.Path(tmp)
        try:
            rec=src.retrieve(SUPPORT,d/"support.zip",64*1024*1024)
            with zipfile.ZipFile(d/"support.zip") as z:
                info=[dict(name=i.filename,bytes=i.file_size,sha256=sha(z.read(i))) for i in z.infolist() if not i.is_dir()]
            assert info
            result["support_archive"]=rec
            result["support_members"]=info
            # Save exact publisher archive as a temporary Actions artifact for inspection, not Git history.
            (out/"publisher_support.zip").write_bytes((d/"support.zip").read_bytes())
        except Exception as e:
            result["support_error"]=type(e).__name__+": "+str(e)[:1000]
        raw=src.retrieve(src.EIDC_URL,d/"original.zip",src.SOURCE_DOWNLOAD_LIMIT)
        result["raw_archive_SHA256"]=raw["SHA256"]
        with zipfile.ZipFile(d/"original.zip") as z:
            headers=[]
            for path,expected_sha in sorted(expected.items()):
                b=z.read(path)
                assert sha(b)==expected_sha,path
                # Header only; no chemistry values or biological records inspected.
                header=next(csv.reader(io.StringIO(b.splitlines()[0].decode("utf-8-sig"))))
                headers.append({"member":path,"sha256":expected_sha,"header":header,"temporal_header_candidates":[x for x in header if any(k in x.lower() for k in ["date","time","start","end","period","window","year","month","day"])]})
            result["chemistry_headers"]=headers
    result["status"]="DOCUMENTS_AVAILABLE_REVIEW_REQUIRED" if "support_archive" in result else "SUPPORT_DOCUMENTS_STOP"
    (out/"metadata_receipt.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("CHEMPOP TEMPORAL METADATA "+json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":
    a=argparse.ArgumentParser();a.add_argument("--work",type=pathlib.Path,required=True);a.add_argument("--out",type=pathlib.Path,required=True);v=a.parse_args();main(v.work,v.out)
