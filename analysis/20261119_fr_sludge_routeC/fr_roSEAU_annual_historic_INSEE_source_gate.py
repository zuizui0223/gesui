#!/usr/bin/env python3
"""Source-only French sewage sludge candidate spreading communes, 2017–2022.

Each year checked against its actual historic INSEE Code Officiel Géographique.
The original ROSEAU CSV contains repeated laboratory analysis rows, NOT
distinct field applications. No plant, butterfly or EQR outcomes are read.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,zipfile
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path
import requests
from fr_roSEAU_2018_original_INSEE_COG_source_gate import parse_cell,quantity,commune_normalized

BASES=(
    "https://www.assainissement.developpement-durable.gouv.fr",
    "https://assainissement.developpement-durable.gouv.fr",
    "https://www2.assainissement.developpement-durable.gouv.fr",
)
JS="/public/js/exportbasededonneeBoue_min.js"
PAGE="/pages/data/basededonneeboue.php"
HANDLER="/controllers/GetDataExportBoueHandler.php"
JS_SHA="f4cf1e2fe84e0ccf608ece6c5d363426960e2b730e21900d48265eba4a78f0e5"
COG={
  "2017":"https://www.insee.fr/fr/statistiques/fichier/2666684/comsimp2017-txt.zip",
  "2018":"https://www.insee.fr/fr/statistiques/fichier/3363419/comsimp2018-txt.zip",
  "2019":"https://www.insee.fr/fr/statistiques/fichier/3720946/communes-01012019-csv.zip",
  "2020":"https://www.insee.fr/fr/statistiques/fichier/4316069/communes2020-csv.zip",
  "2021":"https://www.insee.fr/fr/statistiques/fichier/5057840/commune2021-csv.zip",
  "2022":"https://www.insee.fr/fr/statistiques/fichier/6051727/commune_2022.csv",
}
SOURCE_SHA={
  "2017":"7ef811828d78ba667fb53f08efc208a8fc472418da734fcf8735abf65477280d",
  "2018":"bb9660f13e3dc0b7db359cfc414f14e4510db3b51d73ce0baa4bd5cca821f6c2",
  "2019":"8738e2c385416a01248b4383dca7d2906af19406a653750b1752a6a13469b645",
  "2020":"9178cb764fc30c0be400087be021f49ec64271fa6c5425b02bc6fd3f46c5e48e",
  "2021":"541aa7ec6459b93d0313e8e8ccb1bafc2dc85d3595ade647e9cced519a6b0807",
  "2022":"b75f7191ca03750958548312867350568398a6454adac5c373f3d2ff91e199f0",
}
COMMUNE_RE=re.compile(r"^(?:\d{5}|2[AB]\d{3})$")

def sha(raw):return hashlib.sha256(raw).hexdigest()

def official_cog(year,payload,min_count=30000):
    if year not in COG:raise ValueError("Year outside fixed COG source period")
    if year in {"2017","2018"}:
        z=zipfile.ZipFile(io.BytesIO(payload))
        files=[x for x in z.namelist() if x.lower().endswith(".txt") and f"comsimp{year}" in x.lower()]
        if len(files)!=1:raise ValueError("Official historic INSEE TXT commune file not unique")
        raw=z.read(files[0])
        txt=raw.decode("iso-8859-1")
        delim=csv.Sniffer().sniff(txt[:6000],delimiters=";,\t").delimiter
        rows=iter(csv.reader(io.StringIO(txt),delimiter=delim))
        first=next(rows)
        if "DEP" in first and "COM" in first:
            dep=first.index("DEP");com=first.index("COM");name=first.index("NCC") if "NCC" in first else 9
        else:
            dep,com,name=3,4,9
            rows=iter([first,*list(rows)])
        codes=defaultdict(set)
        for row in rows:
            if len(row)<=max(dep,com,name):continue
            d,c,n=row[dep].strip().upper(),row[com].strip(),row[name].strip()
            if not d or not c or not n:continue
            if len(d)==2 and len(c)==3:code=d+c
            elif len(d)==3 and len(c.lstrip("0"))<=2:code=d+c.lstrip("0").zfill(2)
            else:continue
            if COMMUNE_RE.fullmatch(code):codes[code].add(commune_normalized(n))
        orig_name=files[0]
    else:
        if year=="2022":raw=payload;orig_name="commune_2022.csv"
        else:
            z=zipfile.ZipFile(io.BytesIO(payload))
            files=[x for x in z.namelist() if x.lower().endswith(".csv") and "commune" in x.lower() and
                   "mvt" not in x.lower() and "histor" not in x.lower()]
            if len(files)!=1:raise ValueError("Expected exactly one official INSEE annual commune CSV: "+repr(files))
            orig_name=files[0];raw=z.read(files[0])
        reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter=",")
        if not {"COM","NCC"}.issubset(reader.fieldnames or []):
            raise ValueError("Official newer COG CSV lacks COM/NCC code-name columns")
        codes=defaultdict(set)
        for row in reader:
            code=str(row["COM"]).strip().upper()
            if row.get("TYPECOM") not in ("COM",None,""):
                continue
            if COMMUNE_RE.fullmatch(code):
                for key in ("NCC","NCCENR","LIBELLE"):
                    if row.get(key):codes[code].add(commune_normalized(row[key]))
    if len(codes)<min_count:
        raise ValueError(f"COG {year} parsed too few real commune codes: {len(codes)}")
    return codes,{"official_URL":COG[year],"archive_sha256":sha(payload),
                   "official_inner_commune_file":orig_name,"inner_communes_sha256":sha(raw),
                   "official_unique_communes":len(codes)}

def retrieve_year(year):
    errs=[]
    for base in BASES:
        try:
            c=requests.Session()
            c.headers.update({"User-Agent":"gesui-France-sludge-Cog-annual-source-only/20261119"})
            page=c.get(base+PAGE,timeout=(12,35));page.raise_for_status()
            if year not in page.text:raise ValueError("Official year not in selectable public source")
            js=c.get(base+JS,timeout=(12,35));js.raise_for_status()
            if sha(js.content)!=JS_SHA:raise ValueError("Official ROSEAU public exporter JS drift")
            export=c.post(base+HANDLER,timeout=(15,110),data={
              "bassin":"ALL","region":"ALL","annee":year,"format":"csv","tranche":"ALL","departement":"ALL"
            },headers={"Origin":base,"Referer":base+PAGE})
            export.raise_for_status()
            raw=export.content
            if sha(raw)!=SOURCE_SHA[year]:
                raise ValueError("ROSEAU raw CSV year source SHA differs from the completed independent six-year audit")
            return raw,base
        except Exception as ex:
            errs.append({"official_host":base,"err":str(ex)[:220]})
    raise RuntimeError("All 3 original government hosts could not supply identical source: "+json.dumps(errs))

def inspect(year,raw,cog):
    if sha(raw)!=SOURCE_SHA[year]:raise ValueError("Year's previously locked original CSV SHA changed")
    f=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")),delimiter="|")
    orig=f.fieldnames or []
    if not {"Code SANDRE la station","Commune d'épandage","Quantité épandage agricole (tMS/an)"}.issubset(orig):
        raise ValueError("Original ROSEAU public column schema changed")
    bysite=defaultdict(set);mass=defaultdict(set)
    raw_rows=0
    errors=Counter()
    for row in f:
        raw_rows+=1
        station=row["Code SANDRE la station"].strip()
        if not station:raise ValueError("Source lab row has no SANDRE station ID")
        area=row["Commune d'épandage"].strip()
        qty=quantity(row["Quantité épandage agricole (tMS/an)"])
        bysite[station].add(area)
        mass[station].add(str(qty) if qty is not None else "missing")
    if any(len(v)>1 for v in bysite.values()) or any(len(v)>1 for v in mass.values()):
        raise ValueError("Station lists or station-level agricultural mass vary among repeated lab rows")
    allcodes=set()
    positivecodes=set()
    strict_positive_codes=set()
    all_valid_codes=set()
    all_cog_conflicts=set()
    parsed_tokens=0;unparsed_tokens=0
    positive_stations=set()
    station_code_pairs=0
    positive_pairs=0
    for station,raw_strings in bysite.items():
        cell=next(iter(raw_strings))
        tokens,invalid=parse_cell(cell)
        parsed_tokens+=len(tokens);unparsed_tokens+=len(invalid)
        for tok in invalid:errors[tok]+=1
        pos=mass[station]!={"missing"} and quantity(next(iter(mass[station]))) is not None and quantity(next(iter(mass[station])))>0
        if pos:positive_stations.add(station)
        for code,name in tokens:
            station_code_pairs+=1
            allcodes.add(code)
            if code in cog:all_valid_codes.add(code)
            if code in cog and commune_normalized(name) not in cog[code]:
                all_cog_conflicts.add((code,name))
            if pos:
                positive_pairs+=1
                positivecodes.add(code)
                if code in cog and commune_normalized(name) in cog[code]:
                    strict_positive_codes.add(code)
    return {
      "original_repeated_lab_rows":raw_rows,
      "distinct_STEU_source_stations":len(bysite),
      "positive_agricultural_spreading_mass_stations":len(positive_stations),
      "all_original_station_commune_code_pairs":station_code_pairs,
      "positive_station_commune_code_pairs":positive_pairs,
      "distinct_embedded_code_candidates":len(allcodes),
      "distinct_codes_in_official_same_year_COG":len(all_valid_codes),
      "distinct_positive_station_listed_codes":len(positivecodes),
      "distinct_positive_codes_in_same_year_COG":len(positivecodes & all_valid_codes),
      "distinct_positive_codes_matching_both_COG_code_and_normalized_name":len(strict_positive_codes),
      "distinct_positive_code_departments_prefix":len({x[:2] for x in strict_positive_codes}),
      "all_station_level_commune_tokens_parsed":parsed_tokens,
      "failed_station_level_commune_tokens":unparsed_tokens,
      "station_level_token_parse_fraction":parsed_tokens/(parsed_tokens+unparsed_tokens) if (parsed_tokens+unparsed_tokens) else None,
      "COG_code_name_conflict_pairs_all_source":len(all_cog_conflicts),
      "source_tokens_unparsed_type_count":dict(Counter(k.split("|")[0] for k in errors)),
      "actual_application_year_per_municipality_not_recorded":True,
      "site_total_tonnes_never_allocated_to_commune":True,
      "original_2018_NUTS2_pointdiffuse_stop_preserved":True,
    },sorted(strict_positive_codes)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--year",choices=sorted(COG),required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    output={"schema":"gesui_France_sludge_2017to22_embedded_INSEE_source_year_specific_COG_v1",
            "UTC":datetime.now(timezone.utc).isoformat(),
            "year":args.year,"no_plant_NUTRI_butterfly_HNS_EQR_used":True}
    try:
        ses=requests.Session()
        ses.headers.update({"User-Agent":"gesui-INSEE-history-sludge-source-only/20261119"})
        res=ses.get(COG[args.year],timeout=(15,70));res.raise_for_status()
        official,rec=official_cog(args.year,res.content)
        output["historic_INSEE_COG_source_receipt"]=rec
        src,host=retrieve_year(args.year)
        output["official_ROSEAU_host"]=host
        output["original_ROSEAU_CSV_SHA256"]=sha(src)
        results,codes=inspect(args.year,src,official)
        output["source_counts"]=results
        output["status"]="OFFICIAL_YEAR_COG_AND_ORIGINAL_SLUDGE_CODE_NAME_VALIDATED"
        with (args.out/f"FR_ROSEAU_{args.year}_COG_VALIDATED_CANDIDATE_COMMUNE_CODES.csv").open("w",newline="") as f:
            wr=csv.writer(f);wr.writerow(["year","INSEE_code","interpretation"])
            for code in codes:
                wr.writerow([args.year,code,"candidate_in_station_spreading_commune_list_and_station_positive_tMS_NOT_field_application"])
    except Exception as ex:
        output["status"]="COG_OR_OFFICIAL_SOURCE_ACCESS_STOP"
        output["error"]=str(ex)[:700]
    (args.out/f"FR_ROSEAU_{args.year}_HISTORIC_INSEE_COG_SOURCE_QA.json").write_text(
      json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n")
    print("GESUI_FR_ROSEAU_HISTORIC_YEAR_INSEE",json.dumps({
      "year":args.year,"status":output["status"],
      "counts":output.get("source_counts"),"error":output.get("error")
    },ensure_ascii=False),flush=True)

if __name__=="__main__":main()
