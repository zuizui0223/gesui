#!/usr/bin/env python3
"""Schema-only audit of original Denmark eBMS occurrence extension; no species rows read."""
from __future__ import annotations
import argparse,hashlib,json,urllib.request,zipfile
import xml.etree.ElementTree as ET
from pathlib import Path
SOURCE="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_transects&v=1.4"
EXPECTED="2cc64ec352754223db06a92eef5b01803968eac94bae34057c17f1abdb6da79b"
ROOT="{http://rs.tdwg.org/dwc/text/}"
def main():
    a=argparse.ArgumentParser()
    a.add_argument("--out",required=True,type=Path)
    args=a.parse_args()
    with urllib.request.urlopen(urllib.request.Request(SOURCE,headers={"User-Agent":"gesui-DK-EBMS-schema-only-2026-11"}),timeout=100) as r:
        raw=r.read()
    chk=hashlib.sha256(raw).hexdigest()
    if chk!=EXPECTED:raise ValueError("Original archive changed: "+chk)
    import io
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        xml=z.read("meta.xml")
        root=ET.fromstring(xml)
        cores=root.findall(ROOT+"core");extensions=root.findall(ROOT+"extension")
        if len(cores)!=1:raise ValueError("No single Event core")
        def inspect(node):
            location=node.find(".//"+ROOT+"files/"+ROOT+"location")
            fields={}
            for x in node.findall(ROOT+"field"):
                fields[x.attrib["term"].rsplit("/",1)[-1]]=int(x.attrib["index"])
            ix=node.find(ROOT+"coreid")
            return {"rowtype":node.attrib.get("rowType"),
                    "path":location.text if location is not None else None,
                    "field_terms_to_column_index":fields,
                    "coreid_field_index":ix.attrib.get("index") if ix is not None else None}
        items=[inspect(i) for i in cores+extensions]
        if items[0]["rowtype"]!="http://rs.tdwg.org/dwc/terms/Event":raise ValueError("Not Event source")
        if not any("Occurrence" in x["rowtype"] for x in items[1:]):
            raise ValueError("No source Occurrence extension")
        for x in items:
            if x["path"] not in z.namelist():raise ValueError("Declared extension missing")
        out={"source_SHA256":chk,"source_archive_bytes":len(raw),
             "meta_XML_SHA256":hashlib.sha256(xml).hexdigest(),
             "data_tables_event_core_and_occurrence_extension_schema_only":items,
             "source_file_names":z.namelist(),
             "original_occurrence_row_count_READ":0,
             "original_species_taxon_or_abundance_values_READ":False,
             "already_committed_5km_original_source_gate_expected":True,
             "source_license":"CC-BY-NC 4.0"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
