#!/usr/bin/env python3
"""2018 original French ROSEAU sludge spreading: recover INSEE commune IDs.

Source-only; original public site/year/CSV SHA are locked. Compare embedded
CODE|NAME strings to INSEE COG at Jan 2018, not current merged communes.
A positive STEU-level agricultural sludge mass is NOT a commune-level
application amount; the same station total must never be summed per commune.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,unicodedata,zipfile
from collections import Counter,defaultdict
from decimal import Decimal,InvalidOperation
from datetime import datetime,timezone
from pathlib import Path
import requests

SITE="https://www.assainissement.developpement-durable.gouv.fr"
PAGE=SITE+"/pages/data/basededonneeboue.php"
JS=SITE+"/public/js/exportbasededonneeBoue_min.js"
EXPORT=SITE+"/controllers/GetDataExportBoueHandler.php"
JS_SHA="f4cf1e2fe84e0ccf608ece6c5d363426960e2b730e21900d48265eba4a78f0e5"
ORIG_2018_SHA="bb9660f13e3dc0b7db359cfc414f14e4510db3b51d73ce0baa4bd5cca821f6c2"
COG_2018_ZIP="https://www.insee.fr/fr/statistiques/fichier/3363419/comsimp2018-txt.zip"
FIELD="Commune d'épandage"
SITE_ID="Code SANDRE la station"
SITE_QTY="Quantité épandage agricole (tMS/an)"
PARSER=re.compile(r"^(?P<code>(?:\d{5}|2[AB]\d{3}))\s*\|\s*(?P<name>.+)$",re.I)
MAX_RAWSIZE=30_000_000

def sha(data):return hashlib.sha256(data).hexdigest()

def commune_normalized(text):
    x=unicodedata.normalize("NFKD",str(text).casefold())
    x="".join(c for c in x if not unicodedata.combining(c))
    x=x.replace("œ","oe").replace("æ","ae")
    return "".join(c for c in x if c.isalnum())

def parse_cell(s):
    """One official source cell may contain several CODE|NAME entries."""
    if not s or not s.strip():return [],[]
    toks=[z.strip() for z in s.strip().split(",") if z.strip()]
    if not toks:return [],[s[:120]]
    valid=[];invalid=[]
    for tok in toks:
        match=PARSER.fullmatch(tok)
        if not match:invalid.append(tok[:100]);continue
        code=match.group("code").upper()
        name=match.group("name").strip()
        if not name:invalid.append(tok[:100]);continue
        valid.append((code,name))
    if len({code for code,_ in valid})!=len(valid):
        # one municipality may appear twice, but count once; flag input issue.
        invalid.append("duplicate_CODE_in_original_commune_cell")
    return sorted(set(valid)),invalid

def quantity(s):
    s=(s or "").strip().replace("\u202f","").replace("\u00a0","").replace(" ","").replace(",",".")
    if not re.fullmatch(r"[+]?\d+(?:\.\d+)?",s):return None
    try:
        x=Decimal(s)
        return x if x.is_finite() else None
    except InvalidOperation:return None

def cog2018(zipraw):
    """Parse the actual INSEE 2018 TXT dataset; never infer a modern code."""
    z=zipfile.ZipFile(io.BytesIO(zipraw))
    filenames=[n for n in z.namelist() if n.lower().endswith(".txt") and "comsimp2018" in n.lower()]
    if len(filenames)!=1:raise ValueError("INSEE 2018 ZIP does not include the expected official commune TXT")
    raw=z.read(filenames[0])
    text=raw.decode("iso-8859-1")
    try:
        delim=csv.Sniffer().sniff(text[:6000],delimiters=";,\t").delimiter
    except csv.Error:
        raise ValueError("Official COG2018 historical TXT delimiter could not be determined")
    lines=iter(csv.reader(io.StringIO(text),delimiter=delim))
    first=next(lines)
    # Historic COG 2018 columns are CDC/CHEFLIEU/REG/DEP/COM/...
    if "DEP" in first and "COM" in first:
        idx_dep=first.index("DEP");idx_com=first.index("COM")
        idx_name=first.index("NCC") if "NCC" in first else 9
        allrows=lines
    else:
        idx_dep=3;idx_com=4;idx_name=9
        allrows=iter([first,*list(lines)])
    codes=defaultdict(set)
    for row in allrows:
        if len(row)<=max(idx_dep,idx_com,idx_name):continue
        dep=row[idx_dep].strip().upper()
        com=row[idx_com].strip().upper()
        name=row[idx_name].strip()
        if not dep or not com or not name:continue
        if len(dep)==2 and re.fullmatch(r"\d{3}",com):
            code=dep+com
        elif len(dep)==3 and re.fullmatch(r"\d{2,3}",com) and len(com.lstrip("0"))<=2:
            code=dep+com.lstrip("0").zfill(2)
        else:continue
        if re.fullmatch(r"(?:\d{5}|2[AB]\d{3})",code):
            codes[code].add(commune_normalized(name))
    if len(codes)<30_000:raise ValueError("INSEE 2018 official COG ZIP parsed too few communes")
    return codes,{"official_zip_sha256":sha(zipraw),"inner_txt_sha256":sha(raw),
                  "inner_filename":filenames[0],"code_count":len(codes)}

def inspect(csvraw,cog):
    if sha(csvraw)!=ORIG_2018_SHA:
        raise ValueError("Gov official 2018 public CSV byte digest differs from original archival run")
    s=csvraw.decode("utf-8-sig")
    records=csv.DictReader(io.StringIO(s),delimiter="|")
    original=records.fieldnames or []
    if not all(k in original for k in (SITE_ID,FIELD,SITE_QTY)):
        raise ValueError("Original government ROSEAU source columns missing")
    rows=0
    station_cell=defaultdict(set)
    station_qty=defaultdict(set)
    station_tokens=defaultdict(set)
    failed=Counter()
    all_source_codes=set()
    parsed_slots=0;total_slots=0
    code_name_conflicts=Counter()
    for row in records:
        if not any(v.strip() for v in row.values() if isinstance(v,str)):continue
        rows+=1
        site=(row[SITE_ID] or "").strip()
        if not site:raise ValueError("Source row lacks STEU SANDRE station code")
        cell=(row[FIELD] or "").strip()
        mass=quantity(row[SITE_QTY])
        station_cell[site].add(cell)
        station_qty[site].add(str(mass) if mass is not None else "missing")
        tokens,invalid=parse_cell(cell)
        total_slots+=len(tokens)+len(invalid)
        parsed_slots+=len(tokens)
        for err in invalid:failed[err]+=1
        for code,name in tokens:
            all_source_codes.add(code)
            station_tokens[site].add((code,name))
            if code in cog and commune_normalized(name) not in cog[code]:
                code_name_conflicts[(code,name)]+=1
    if rows!=1479 or len(station_cell)!=55:
        raise ValueError("Original 2018 CSV rows/station count drift")
    if any(len(v)>1 for v in station_cell.values()):
        raise ValueError("Distinct sludge spreading municipality sets within one source station")
    if any(len(v)>1 for v in station_qty.values()):
        raise ValueError("Station agricultural spreading mass differs between repeated lab-assay rows")
    positive_stations={s for s,v in station_qty.items() if v!={"missing"} and Decimal(next(iter(v)))>0}
    positive_pairs={(s,code) for s in positive_stations for code,_ in station_tokens[s]}
    all_pairs={(s,code) for s,toks in station_tokens.items() for code,_ in toks}
    codes_positive={code for s,code in positive_pairs}
    code_found={c for c in all_source_codes if c in cog}
    # COG 2018 name matches can differ for article, hyphen, town rename;
    # explicitly report rather than auto-reassigning a 5-digit ID.
    return {
      "year":2018,
      "row_count_original_repeated_lab_assays":rows,
      "distinct_STEU_station_ids":len(station_cell),
      "commune_cell_parse_total_tokens_including_repeated_assay_rows":total_slots,
      "commune_cell_parse_valid_code_name_tokens_including_repeated_assay_rows":parsed_slots,
      "parse_token_coverage_fraction":parsed_slots/total_slots if total_slots else None,
      "invalid_token_categories":{k:v for k,v in failed.items() if len(k)<100},
      "distinct_embedded_INSEE_code_candidates":len(all_source_codes),
      "distinct_INSEE_code_candidates_in_official_2018_COG":len(code_found),
      "fraction_distinct_codes_in_official_2018_COG":len(code_found)/len(all_source_codes) if all_source_codes else None,
      "codes_with_2018_COG_name_discrepancy_n":len(code_name_conflicts),
      "source_pair_station_INSEE_raw_all":len(all_pairs),
      "station_with_positive_agricultural_spreading_mass_n":len(positive_stations),
      "station_INSEE_code_pairs_with_station_positive_mass_n":len(positive_pairs),
      "distinct_candidate_application_INSEE_codes_at_positive_mass_STEU":len(codes_positive),
      "distinct_2018_COG_valid_application_INSEE_at_positive_mass_STEU":len(codes_positive&code_found),
      "candidate_application_departments_code_prefix_n":len({code[:2] for code in codes_positive}),
      "source_stations_report_total_agricultural_spread_tonnes_NOT_allocated_per_commune":True,
      "source_commune_LIST_does_not_prove_each_municipality_applied_sludge_in_2018":True,
      "source_COG_code_name_discrepancies_require_manual_resolution":True,
      "one_line_per_application_event_assumption_FORBIDDEN":True,
    },sorted(codes_positive&code_found)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    session=requests.Session()
    session.headers.update({"User-Agent":"gesui-French-INSEE2018-COG-sludge-only-source-audit/20261119"})
    result={"schema":"gesui_FR_ROSEAU_2018_original_CSV_embedded_INSEE_vs_official_historic_COG_source_gate_v1",
            "UTC":datetime.now(timezone.utc).isoformat(),
            "no_plant_butterfly_or_EQR_outcomes_read":True}
    try:
        zipr=session.get(COG_2018_ZIP,timeout=(20,65))
        zipr.raise_for_status()
        official,officialreceipt=cog2018(zipr.content)
        result["COG2018_official_INSEE_archive"]=officialreceipt
        page=session.get(PAGE,timeout=(15,40))
        page.raise_for_status()
        js=session.get(JS,timeout=(15,40))
        js.raise_for_status()
        if sha(js.content)!=JS_SHA:raise ValueError("French government exporter JS source drift")
        post=session.post(EXPORT,timeout=(18,110),data={
            "bassin":"ALL","region":"ALL","annee":"2018","format":"csv",
            "tranche":"ALL","departement":"ALL"},
            headers={"Origin":SITE,"Referer":PAGE})
        post.raise_for_status()
        rec,codes=inspect(post.content,official)
        result["source_QA"]=rec
        result["source_2018_CSV_sha256"]=sha(post.content)
        result["status"]="INSEE2018_COG_AND_APPLICATION_CANDIDATES_AUDITED"
        with (a.out/"FR_2018_ROSEAU_PUBLIC_APPLICATION_COMMUNE_CODES_COG_VALIDATED.csv").open("w",newline="",encoding="utf-8") as f:
            w=csv.writer(f);w.writerow(["2018_commune_code_COG_validated","scope"])
            for code in codes:w.writerow([code,"positive_station_mass_and_listed_municipality_NOT_parcel"])
    except Exception as e:
        result["status"]="INSEE_COG_OR_FRENCH_PUBLIC_SOURCE_ACCESS_STOP"
        result["error"]=str(e)[:800]
    (a.out/"FR2018_ROSEAU_APPLICATION_MUNICIPALITIES_INSEE2018_QA.json").write_text(
      json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n")
    print("GESUI_FR2018_COG_SLUDGE_SOURCE",json.dumps({
        "status":result["status"],
        "new_src":result.get("source_QA"),"error":result.get("error")
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
