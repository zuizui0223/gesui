#!/usr/bin/env python3
"""Private ephemeral ID crosswalk from exact official ECHO NPDES outfalls ZIP.

The source-frozen NWCA2021/2026 outfall study retains 160 SHA256(admin permit)
IDs but not the unhashed EPA permit values. This helper reads EXACT historical
source ZIP SHA and resolves only those 160 public permit strings into an
ephemeral /tmp pattern file for rapid streaming 9.15GB FY2021 DMR full CSV
with ripgrep. It NEVER commits or prints raw permit IDs or plant outcomes.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,zipfile
from pathlib import Path

SOURCE=Path("results/us_NWCA2021_plant_TNTP_20261110/EPA_NWCA2021_NPDES_SEWAGE_OUTFALL_568SITE_SOURCE_ONLY_GEOGRAPHY_GATE.json")
ECHO_SHA="f0d1310b1c07085c42527e2ede5cf97ab83d517d50548b7dadb7621a44aa6160"
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for x in iter(lambda:f.read(1024*1024),b""):h.update(x)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--echo-outfall-zip",type=Path,required=True)
    p.add_argument("--output-patterns",type=Path,required=True)
    a=p.parse_args()
    root=Path(__file__).resolve().parents[1]
    source=json.loads((root/SOURCE).read_text())
    if sha(a.echo_outfall_zip)!=ECHO_SHA or source["source_EPA_ECHO_NPDES_2026_live_weekly_original_ZIP_SHA256"]!=ECHO_SHA:
        raise ValueError("Cannot invert source SHA against different 2026 ECHO NPDES point file")
    rows=source["all_original_568_NWCA_source_site_to_sewage_source_distance_receipts_without_plant_or_TNTP_values"]["strict"]
    near=[r for r in rows if r["outfall_rows_within_10km"]>0]
    hashed={k for r in near for k in r["independent_10km_NPDES_permits_SHA256"]}
    if len(near)!=128 or len(hashed)!=160:
        raise ValueError("Original source 128 nearby wetlands / 160 administrative NPDES permits changed")
    with zipfile.ZipFile(a.echo_outfall_zip) as z:
        names=[q for q in z.namelist() if q.split("/")[-1].lower()=="npdes_outfalls_layer.csv"]
        if len(names)!=1:raise ValueError("Unexpected original EPA NPDES outfall member")
        with z.open(names[0]) as raw:
            source_csv=csv.reader(io.TextIOWrapper(raw,encoding="utf-8-sig",errors="replace",newline=""))
            head=[col.upper().strip() for col in next(source_csv)]
            if "EXTERNAL_PERMIT_NMBR" not in head:raise ValueError("Missing original NPDES admin permit key")
            idx=head.index("EXTERNAL_PERMIT_NMBR")
            lookup={}
            for row in source_csv:
                if len(row)<=idx:continue
                permit=row[idx].strip().upper()
                if not permit or permit in lookup:continue
                h=hashlib.sha256(permit.encode()).hexdigest()
                lookup[permit]=permit if h in hashed else None
    found={p for p,hit in lookup.items() if hit}
    if {hashlib.sha256(p.encode()).hexdigest() for p in found}!=hashed or len(found)!=160:
        raise ValueError("Original 160 near-source NPDES public permits cannot be resolved exactly")
    a.output_patterns.parent.mkdir(parents=True,exist_ok=True)
    a.output_patterns.write_text("\n".join(sorted(found))+"\n",encoding="utf-8")
    print(json.dumps({"original_epa_outfall_zip_sha256":ECHO_SHA,
                      "original_160_NEAR_WWTP_NWCA_PERMITS_unhashed_only_in_ephemeral_runner":len(found),
                      "temp_output_stays_private_outside_GitHub_commit":True},indent=2))
if __name__=="__main__":
    main()
