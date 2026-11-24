#!/usr/bin/env bash
# Source-only connectivity and service-directory probe. No biological outcomes.
set -u
OUT="$1"; mkdir -p "$OUT"
for u in https://sdi.eea.europa.eu/ https://sdi.eea.europa.eu/datashare/ https://www.eea.europa.eu/ https://water.discomap.eea.europa.eu/arcgis/rest/services?f=pjson https://discodata.eea.europa.eu/ https://cmshare.eea.europa.eu/ ; do
  echo "HOST $u -> $(curl -s -o /dev/null -w '%{http_code} %{time_total}s' -m 90 "$u")"
done
for d in WISE_UWWTD UWWTD; do
  curl -s -m 90 "https://water.discomap.eea.europa.eu/arcgis/rest/services/$d?f=pjson" -o "$OUT/dir_$d.json"
  echo "DIR $d:"; python3 -c "import json,sys;d=json.load(open(sys.argv[1]));[print(' ',s['name'],s['type']) for s in d.get('services',[])];print(' folders',d.get('folders'))" "$OUT/dir_$d.json" 2>&1 | head -80
done
