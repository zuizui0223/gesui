#!/usr/bin/env python3
"""Verify 2020 EEA public WWTP/outfall dataset access without downloading biology.

Probe the 2019–2020 historical official public catalogue URLs and EEA public
folder protocols. Never conflate a public catalogue entry with actual observed
2020 nutrient loads. Record negative HTTP outcomes.
"""
import argparse,base64,hashlib,json,urllib.request,urllib.error
from pathlib import Path
from urllib.parse import quote
from xml.etree import ElementTree as ET
from datetime import datetime,timezone

UA="gesui/2020-original-source-only-archive-access-QC"
SHARE="sptXqwkQr5g7Bp5"  # Public EEA Nextcloud link observed in official catalogue redirects.
DATASETS={
    "public_2020_treatment":"eea_v_4258_100_m_uwwtd-tp_p_2019-2020_v08_r00",
    "public_2020_discharge":"eea_v_4258_100_m_uwwtd-disc-pts_p_2019-2020_v08_r00",
}
REST={
    "public_2020_treatment":
        "https://water.discomap.eea.europa.eu/arcgis/rest/services/UWWTD/SDI_TreatmentPlants_public_ETRS89_v1/MapServer?f=pjson",
    "public_2020_discharge":
        "https://water.discomap.eea.europa.eu/arcgis/rest/services/UWWTD/SDI_DischargePoints_public_ETRS89_v1/MapServer?f=pjson",
}
def run(url, method="GET", user=None):
    h={"User-Agent":UA,"Accept":"*/*"}
    if method=="PROPFIND":
        h["Depth"]="1"
        h["Content-Type"]="application/xml"
    if user:
        h["Authorization"]="Basic "+base64.b64encode((user+":").encode()).decode()
    q=urllib.request.Request(url,method=method,headers=h,
                            data=b'<?xml version="1.0"?><d:propfind xmlns:d="DAV:"><d:allprop/></d:propfind>' if method=="PROPFIND" else None)
    try:
        with urllib.request.urlopen(q,timeout=26) as resp:
            size=120000 if method!="HEAD" else 0
            b=resp.read(size)
            out={"status":resp.status,"final_url":resp.url,
                 "content_type":resp.headers.get("Content-Type"),
                 "content_length":resp.headers.get("Content-Length"),
                 "head":b[:800].decode("utf-8","replace"),
                 "first120KB_sha256":hashlib.sha256(b).hexdigest()}
            if resp.headers.get("Content-Type","").startswith(("application/json","text/plain")):
                try:
                    j=json.loads(b)
                    out["JSON_keys"]=list(j)[:20]
                    out["JSON_error"]=j.get("error")
                    out["GIS_layers"]=[{"id":z.get("id"),"name":z.get("name")} for z in j.get("layers",[])][:30]
                except Exception:pass
            if method=="PROPFIND" and resp.status==207:
                try:
                    root=ET.fromstring(b)
                    ns={"d":"DAV:"}
                    out["webdav_hrefs"]=[x.text for x in root.findall(".//d:href",ns)][:70]
                except Exception as ex:out["webdav_parse_error"]=str(ex)
            return out
    except urllib.error.HTTPError as ex:return {"status":ex.code,"error":str(ex),
                "final_url":ex.url}
    except Exception as ex:return {"error":type(ex).__name__+": "+str(ex)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);args=ap.parse_args()
    probes={}
    for label,ds in DATASETS.items():
        webdav="https://sdi.eea.europa.eu/datashare/public.php/webdav/"+quote(ds)
        share="https://sdi.eea.europa.eu/datashare/s/"+SHARE+"/download?path=/"+quote(ds)
        probes[label]={
            "catalogue_public_path":ds,
            "public_webdav_method_PROPFIND":run(webdav,"PROPFIND",SHARE),
            "public_webdav_no_auth_PROPFIND":run(webdav,"PROPFIND"),
            "public_specific_folder_DOWNLOAD_HEAD":run(share,"HEAD"),
            "EEA_public_MapServer_response":run(REST[label]),
        }
        print("GESUI_PUBLIC_2020_SPATIAL_SOURCE_ACCESS",json.dumps({
            "kind":label,"results":{k:v.get("status",v.get("error")) for k,v in probes[label].items()
                if isinstance(v,dict)}},ensure_ascii=False),flush=True)
    obj={"schema":"gesui_2020_EEA_public_GIS_archive_and_REST_access_gate_v01",
         "created_utc":datetime.now(timezone.utc).isoformat(),
         "official_EEA_public_catalogue_provenance_only":True,
         "no_original_2020_exposure_measurements_claimed":True,
         "no_biological_EQR_opened":True,
         "probes":probes}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(obj,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_PUBLIC_2020_ARCHIVE_ACCESS_QC_DONE",flush=True)
if __name__=="__main__":main()
