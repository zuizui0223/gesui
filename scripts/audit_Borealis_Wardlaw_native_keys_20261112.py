#!/usr/bin/env python3
"""Borealis Wardlaw 8-WWTP 2026 source keys, no numerical effect.

Read exactly SHA-pinned authentic public Borealis Dataverse exported-tab records
from independent source-only GitHub Actions #43163206531. The Dataverse
tabular API export differs byte-for-byte from ingested original (original
stored MD5), so source SHA here is *public Dataverse normalized export* SHA.
"""
import collections,csv,hashlib,io,json,os,pathlib,re
BASE=pathlib.Path(os.environ["BOREALIS_SOURCE_ARTIFACT"])
EXPECTED={
 "Borealis1219447__2026-09-11_Wardlaw.READ_ME.txt":"d627b572e97cd7d369cc83e0477b8b641808a48ce90917e51953ee094e6e999a",
 "Borealis1213261__22.23.MP.Invert.Summary.tab":"1b69eb88fe38f601b11ff466d6a9b3177546a9b489a74c61eaf2a64d770dec66",
 "Borealis1213287__22.23.MP.Invert.data.characteristics.tab":"31f14836543456aa418dc21e7a0e3c652f03acf00f6d77e975dd69184b50b71a",
 "Borealis1213230__FTIR_spectra_library_matches":"3078d27faa8c6a9c9dfd2dbe7a6989831b9f96b5c81381d06503f8489c4e3bc0",
 "Borealis1219236__Stable_Isotope_data.tab":"2cc350defbec3210bccaba49bab80be4286b6a916626c315c2f3e6f40181645e"
}
for filename,digest in EXPECTED.items():
 p=BASE/filename
 if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
    raise RuntimeError("Original Borealis normalized source byte digest mismatch "+filename)
q=json.loads((BASE/"BOREALIS_WARDLAW_2026_8WWTP_SOURCE_ONLY_QA.json").read_text())
if q["source_DOI"]!="doi:10.5683/SP4/U0MN9S" or q["dataset_id"]!=1212155:
    raise RuntimeError("Borealis source DOI / dataset ID modified")
def from_tab(filename):
 p=BASE/filename
 s=p.read_text(encoding="utf-8-sig")
 # Only separate using authentic tab-delimited Dataverse export, not based on
 # outcome values. Default read all rows as string and no numerical coercion.
 rdr=csv.DictReader(io.StringIO(s,newline=""),delimiter="\t")
 h=rdr.fieldnames or []
 if not len(h)>2:raise RuntimeError("Not a proper public Dataverse TSV table "+filename+" headers "+repr(h[:10]))
 rows=list(rdr)
 if len(rows)<10:raise RuntimeError("Too few original records "+filename)
 return h,rows
source={}
for fname in EXPECTED:
 if fname.endswith(".txt"):
    text=(BASE/fname).read_text(encoding="utf-8-sig",errors="replace")
    print("README_FIRST_9500_CHARS_BEGIN\n"+text[:9500]+"\nREADME_FIRST_9500_CHARS_END",flush=True)
    source[fname]={"kind":"original_readme","sha256":EXPECTED[fname],"bytes":len(text)}
    continue
 headers,rows=from_tab(fname)
 native=[v for v in headers if re.search(r"(site|station|location|upstream|downstream|taxon|family|species|life|stage|sample|year|month|time|replicate|wwtp|treat|stream|code)",v,re.I)]
 listed={}
 for col in native:
    cnt=collections.Counter((r.get(col) or "").strip() for r in rows)
    listed[col]={"nonblank":sum(v for k,v in cnt.items() if k and k.upper() not in ("NA","N/A","NULL","NAN")),
        "n_unique_nonblank":len([k for k in cnt if k and k.upper() not in ("NA","N/A","NULL","NAN")]),
        "top_categories":cnt.most_common(20) if len(cnt)<=40 else [],
        "example_values":sorted([k for k in cnt if k])[:14] if len(cnt)<=40 else []}
 source[fname]={"n_rows":len(rows),"columns":headers,"candidate_native_fields":listed,
                "native_file_export_sha256":EXPECTED[fname]}
 print("\nBOREALIS_ORIGINAL_TABLE_HEADERS",fname,"rows",len(rows),json.dumps(headers),flush=True)
 print("BOREALIS_ORIGINAL_TABLE_KEYS",fname,json.dumps(listed,sort_keys=True,ensure_ascii=False),flush=True)
out={
 "dataset_DOI":"doi:10.5683/SP4/U0MN9S",
 "Borealis_release":"2026-09-11",
 "source":"official Canadian Dataverse public exported normalized TSV and original README",
 "source_export_SHA256_pins":EXPECTED,
 "original_Dataverse_ingested_file_MD5_is_NOT_same_as_public_export_TAB":True,
 "not_a_published_WWTP_trophic_causality_result":True,
 "no_contaminant_numerical_effects_compared":True,
 "raw_source_grain_audit":source
}
outpath=pathlib.Path(os.environ.get("RUNNER_TEMP","/tmp"))/"BOREALIS_WARDLAW_8WWTP_NATIVE_GRAIN_SOURCE_ONLY_FREEZE_20261112.json"
outpath.write_text(json.dumps(out,indent=2,sort_keys=True,ensure_ascii=False))
print("BOREALIS_WARDLAW_NATIVE_GRAIN_SOURCE_END")
