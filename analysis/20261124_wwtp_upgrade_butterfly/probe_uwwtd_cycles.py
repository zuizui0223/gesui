#!/usr/bin/env python3
"""Source-only discovery of EEA UWWTD reporting cycles. Never reads biological outcomes."""
import base64, json, re, sys, urllib.request, urllib.error
from pathlib import Path
from xml.etree import ElementTree as ET

OUT = Path(sys.argv[1]); OUT.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "gesui-wwtp-upgrade-source-gate"}
TOKENS = ["sptXqwkQr5g7Bp5", "mGmXg3zcEqKQLTG", "HCxTFz3Ac8xbHX6", "KkqFMjiif9eS7bn"]
DAV = "https://sdi.eea.europa.eu/datashare/public.php/webdav/"
log = {}

def req(url, method="GET", data=None, headers=None, timeout=60):
    h = dict(UA); h.update(headers or {})
    r = urllib.request.Request(url, data=data, headers=h, method=method)
    with urllib.request.urlopen(r, timeout=timeout) as resp:
        return resp.status, resp.read()

def propfind(token, path="", depth="1"):
    auth = "Basic " + base64.b64encode((token + ":").encode()).decode()
    st, body = req(DAV + path, "PROPFIND", headers={"Authorization": auth, "Depth": depth})
    ns = {"d": "DAV:"}
    out = []
    for r in ET.fromstring(body).findall("d:response", ns):
        href = r.find("d:href", ns).text
        size = r.find(".//d:getcontentlength", ns)
        out.append({"href": href, "bytes": int(size.text) if size is not None and size.text else None})
    return out

for t in TOKENS:
    try:
        top = propfind(t)
        entries = []
        for e in top[1:]:
            entries.append(e)
            rel = e["href"].split("/webdav/", 1)[-1]
            if e["bytes"] is None and rel:
                try:
                    entries += [x for x in propfind(t, rel, "infinity")[1:]][:400]
                except Exception as ex:
                    entries.append({"href": rel, "error": repr(ex)[:200]})
        log[t] = entries
    except Exception as ex:
        log[t] = {"error": repr(ex)[:300]}

# GeoNetwork catalogue search for UWWTD records
cat = {}
body = json.dumps({"query": {"query_string": {"query": "UWWTD OR \"Urban Waste Water\""}},
                   "size": 100, "_source": ["resourceTitleObject.default", "link", "resourceTemporalDateRange", "metadataIdentifier"]}).encode()
for url in ["https://sdi.eea.europa.eu/catalogue/srv/api/search/records/_search",
            "https://sdi.eea.europa.eu/catalogue/srv/api/search/records/_search?bucket=metadata"]:
    try:
        st, b = req(url, "POST", body, {"Content-Type": "application/json", "Accept": "application/json"})
        hits = json.loads(b)["hits"]["hits"]
        cat[url] = [{"id": h["_source"].get("metadataIdentifier"),
                     "title": (h["_source"].get("resourceTitleObject") or {}).get("default"),
                     "links": [l.get("urlObject", {}).get("default") for l in (h["_source"].get("link") or [])][:20]}
                    for h in hits]
        break
    except Exception as ex:
        cat[url] = {"error": repr(ex)[:300]}

(OUT / "dav_listing.json").write_text(json.dumps(log, indent=1))
(OUT / "catalogue_uwwtd.json").write_text(json.dumps(cat, indent=1, ensure_ascii=False))
for t, v in log.items():
    print("==", t)
    if isinstance(v, dict): print(v); continue
    for e in v[:80]: print(e.get("bytes"), e["href"])
for u, v in cat.items():
    if isinstance(v, list):
        for h in v: print(h["title"], h["id"])
    else: print(u, v)
