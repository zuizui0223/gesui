#!/usr/bin/env python3
"""Read only public Dryad v2 file-name metadata for candidate independent field
emergence/consumer studies after *binary file* API returned HTTP 401.
No scraping access challenge, no unofficial mirror, no ecological measurements.
"""
import hashlib,json,os,pathlib,urllib.request,urllib.parse
panels={
 "A_ON_29":{"doi":"doi:10.5061/dryad.sbcc2frfk","versionid":342267,"expectedfiles":9},
 "B_QC_20":{"doi":"doi:10.5061/dryad.v6wwpzh8z","versionid":430321,"expectedfiles":2}
}
out={}
for label,meta in panels.items():
    u=f"https://datadryad.org/api/v2/versions/{meta['versionid']}/files"
    req=urllib.request.Request(u,headers={"Accept":"application/json","User-Agent":"gesui-2026-source-only-public-file-inventory/1"})
    with urllib.request.urlopen(req,timeout=40) as r:raw=r.read()
    j=json.loads(raw)
    files=j.get("_embedded",{}).get("stash:files")
    if not isinstance(files,list) or len(files)!=meta["expectedfiles"]:
        raise RuntimeError("Original Dryad metadata file list changed "+label)
    entries=[]
    for f in files:
        x={k:f.get(k) for k in ["id","path","size","digest","description","mimeType","status"]}
        x["download_link_from_original_API"]=f.get("_links",{}).get("stash:download")
        entries.append(x)
    out[label]={
       "doi":meta["doi"],"versionid":meta["versionid"],
       "original_public_API_url":u,"original_API_body_sha256":hashlib.sha256(raw).hexdigest(),
       "n_source_files":len(files),"files":entries,
       "binary_download_state_from_prior_source_run":"HTTP401_api_v2_files_download_not_accessible_without_authorization",
       "not_a_raw_file_ingestion":True,
    }
root=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))
path=root/"DRYAD_FIELD_STUDIES_PUBLIC_FILE_INVENTORY_ONLY_20261112.json"
path.write_text(json.dumps(out,indent=2,sort_keys=True))
print("DRYAD_PUBLIC_FILE_INVENTORY_START")
print(json.dumps(out,sort_keys=True))
print("DRYAD_PUBLIC_FILE_INVENTORY_END")
