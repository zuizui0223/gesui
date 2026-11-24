#!/usr/bin/env bash
# Source-only connectivity probe. No biological outcomes.
set -u
OUT="$1"; mkdir -p "$OUT"
getent ahosts sdi.eea.europa.eu | head -4
for f in -4 -6; do
  echo "HOST sdi $f -> $(curl $f -s -o /dev/null -w '%{http_code} %{time_total}s' --connect-timeout 25 -m 40 https://sdi.eea.europa.eu/webdav/datastore/public/)"
done
echo "HOST datahub -> $(curl -s -o /dev/null -w '%{http_code}' -m 40 'https://www.eea.europa.eu/en/datahub/datahubitem-view/6244937d-1c2c-47f5-bdf1-33ca01ff1715')"
