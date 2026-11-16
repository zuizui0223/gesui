#!/usr/bin/env python3
"""Real same LARVAL HOST PLANT paired NUTRI+MOIST ecological indicator source.

Ferrer-Paris et al 2014 43,645 published butterfly-host records (same source
as original roster, GloBI pinned rehost) joined strict exact ACCEPTED butterfly
binomials and actual EIV Europe 2023 double-scored HOST plant binomial species.
No guessed synonyms, genus imputation or ecological wastewater outcome.

Source availability only; step 2 reconstruct H7 community scores separately
if >=40 accepted butterflies across>=4 families have >=1 double-scored host.
"""
from __future__ import annotations
import argparse,collections,csv,hashlib,io,json,math,re
from pathlib import Path
import requests
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
SCHEMA="GESUI_20261117_FERRERPARIS2014_EIV2023_SAME_LARVAL_HOST_N_MOIST_EXACT_SOURCE_GATE_V1"
HOST_REPO_SHA="2c5b20bd8c0b2839124c71c83ca83369f1d812c9"
HOST_GITBLOB="ea93a5c0fdb17e361a9b6a9637a497a12e90493c"
HOST_URL=("https://raw.githubusercontent.com/globalbioticinteractions/ferrer-paris/"
          +HOST_REPO_SHA+"/20140912_CompiledButterflyHostPlantRecords_JRFP.csv")
HOST_N_ORIGINAL=43645
EIV_URL="https://www.sci.muni.cz/botany/juice/NET/sources/EIV_Europe.csv"
EIV_SHA="ffeb2a5f45c95886625f970734e431098587c65baad7def786d35364a761dbba"
ROSTER=ROOT/"results/20261115_sync/GESUI_303_HOST_NUTRIENT_SCORES.csv"
ROSTER_SHA="00c1a36a9237a3d1a6d77334aa8baa5dc8b329a49925c2e31f033fa2d302ecdb"
SOURCE_FAMILY={"Hesperiidae","Lycaenidae","Nymphalidae","Papilionidae","Pieridae"}
N_MIN=40
FAMILIES_MIN=4
# Original EIV European:
EIV_COLUMNS=["Unnamed: 0","LIGHT","TEMPE","MOIST","REACT","NUTRI","SALIN"]

def digest(s):return hashlib.sha256(s).hexdigest()
def gitblob(s):return hashlib.sha1(("blob "+str(len(s))+"\0").encode()+s).hexdigest()
def name(s):
    if pd.isna(s):return ""
    return re.sub(r"\s+"," ",str(s).strip()).casefold()
def exact_species(s):
    z=s.split()
    return len(z)==2 and all(re.fullmatch("[a-z][a-z-]+",word) for word in z)
def original_roster():
    if digest(ROSTER.read_bytes())!=ROSTER_SHA:
        raise RuntimeError("Original 303 butterfly species source SHA changed")
    d=pd.read_csv(ROSTER,dtype={"butterflyKey":str},encoding="utf-8-sig")
    if len(d)!=303 or d.butterflyKey.duplicated().any() or set(d.family)!=SOURCE_FAMILY:
        raise RuntimeError("Original eligible 303 butterfly data taxonomy or family set altered")
    accepted=(d.taxon_status=="ACCEPTED")&(
        ~d.actual_used_alias_conflict.astype(str).str.lower().isin(["true","yes","1"]))
    if int(accepted.sum())!=290:
        raise RuntimeError("Original strict 290 accepted/alias clean butterflies source changed")
    unique_names=d.loc[accepted,"accepted_name"].map(name)
    if not unique_names.map(exact_species).all() or unique_names.duplicated().any():
        raise RuntimeError("Original accepted taxonomic names not unique valid species binomials")
    return d,accepted
def download_original(url,source_sha=None,source_gitblob=None):
    r=requests.get(url,timeout=(12,130),headers={
        "User-Agent":"gesui-butterfly-LARVALHOST-EIV-20261117-academic/1.0","Accept":"text/csv"})
    if r.status_code!=200:
        raise RuntimeError("Original public source inaccessible HTTP "+str(r.status_code))
    data=r.content
    if source_sha and digest(data)!=source_sha:
        raise RuntimeError("Published exact original EIV source SHA256 changed "+digest(data))
    if source_gitblob and gitblob(data)!=source_gitblob:
        raise RuntimeError("Published original pinned GloBI/Ferrer-Paris source git blob changed "+gitblob(data))
    return data
def host_edges():
    raw=download_original(HOST_URL,source_gitblob=HOST_GITBLOB)
    source=io.StringIO(raw.decode("utf-8-sig"))
    rr=csv.DictReader(source)
    original_head=rr.fieldnames
    def field(valid):
        matching=[c for c in valid if c in original_head]
        if len(matching)!=1:
            raise RuntimeError("GloBI published original Ferrer-Paris schema field mismatch for "+str(valid)+": "+str(original_head))
        return matching[0]
    butterfly=field(["sourceTaxonName","Lepidoptera Name"])
    host=field(["targetTaxonName","Hostplant Name"])
    tab=collections.defaultdict(set)
    all_names=set()
    nrows=0;n_double_species=0
    invalid=collections.Counter()
    for row in rr:
        nrows+=1
        b=name(row[butterfly]);h=name(row[host])
        if not exact_species(b):
            invalid["original_butterfly_not_exact_binomial"]+=1;continue
        all_names.add(b)
        if not exact_species(h):
            invalid["original_host_not_exact_binomial"]+=1;continue
        tab[b].add(h)
        n_double_species+=1
    if nrows!=HOST_N_ORIGINAL:
        raise RuntimeError("Original published Ferrer-Paris larval-host record count changed "+str(nrows))
    return tab,{"source_publisher":"Ferrer-Paris et al. 2014 / original Figshare data, GloBI pin rehost",
        "original_DOI":"10.6084/m9.figshare.1168861.v1",
        "source_rehosting_GitHub_commit_SHA":HOST_REPO_SHA,
        "original_rehosted_source_gitblob_SHA1":HOST_GITBLOB,
        "source_original_2014_url":HOST_URL,
        "source_original_2014_csv_bytes":len(raw),
        "source_original_columns":original_head,
        "source_original_host_source_record_rows":nrows,
        "source_original_butterfly_unique_binomial_names":len(all_names),
        "source_usable_exact_species_binomial_edges_before_ROSTER_match":n_double_species,
        "source_usable_butterfly_bins_w_numeric_hosts":len(tab),
        "rejected_original_malformed_species_concepts":dict(invalid),
        "original_all_source_rows_NOT_REPUBLISHED":True}
def eiv_plant_N_and_F():
    b=download_original(EIV_URL,source_sha=EIV_SHA)
    data=pd.read_csv(io.BytesIO(b),encoding="utf-8-sig",low_memory=False)
    if len(data)!=8908 or data.columns.tolist()!=EIV_COLUMNS:
        raise RuntimeError("Original European 2023 EIV 8908 row 6 trait source changed")
    data["name"]=data["Unnamed: 0"].map(name)
    data=data.loc[data.name.map(exact_species)].copy()
    for field,lo,hi in [("NUTRI",1,9),("MOIST",1,12)]:
        data[field]=pd.to_numeric(data[field],errors="coerce")
        data.loc[~data[field].between(lo,hi),field]=np.nan
    # Some duplicated names with different indicator values are ambiguous;
    # do not use whichever first happened to occur.
    ambiguous=set()
    for field in ("NUTRI","MOIST"):
        vals=data.groupby("name")[field].agg(["min","max"])
        ambiguous|=set(vals.index[vals["max"]-vals["min"]>1e-10])
    filt=data.loc[~data.name.isin(ambiguous)].drop_duplicates("name")
    paired=filt.loc[filt.NUTRI.notna()&filt.MOIST.notna()]
    out={str(z.name):tuple(map(float,(z.NUTRI,z.MOIST))) for z in paired.itertuples()}
    return out,{"original_eiv_source_SHA256":EIV_SHA,
        "original_published_EIV_rows":8908,
        "source_strict_EIV_bino_NF_paired_plants":len(out),
        "source_conflicting_exact_EIV_species_names_excluded":len(ambiguous),
        "source_original_eiv_file_copyright_not_rehosted":True}
def main(a):
    a.output.mkdir(parents=True,exist_ok=True)
    roster,accepted=original_roster()
    edges,source=host_edges()
    trait,source_eiv=eiv_plant_N_and_F()
    rows=[];families=collections.Counter()
    counted=0
    for z in roster.itertuples(index=False):
        bname=name(z.accepted_name)
        orig_hosts=edges.get(bname,set())
        matched=sorted(orig_hosts & trait.keys())
        n=len(matched)
        eligible=(z.taxon_status=="ACCEPTED" and not bool(z.actual_used_alias_conflict) and n>=1)
        if eligible:
            counted+=1;families[z.family]+=1
        nmean=float(np.mean([trait[h][0] for h in matched])) if matched else None
        fmean=float(np.mean([trait[h][1] for h in matched])) if matched else None
        oldN=float(z.host_nutrient_mean) if pd.notna(z.host_nutrient_mean) else None
        oldF=float(z.host_moisture_mean) if pd.notna(z.host_moisture_mean) else None
        rows.append({"butterflyKey":str(z.butterflyKey),
           "original_butterfly_accepted_name":str(z.accepted_name),
           "original_butterfly_family":str(z.family),
           "original_butterfly_taxon_status":str(z.taxon_status),
           "original_used_alias_conflict":str(z.actual_used_alias_conflict),
           "source_exact_butterfly_name_found_in_2014_HOSTS":bool(bname in edges),
           "source_2014_distinct_exact_binomial_host_taxa":len(orig_hosts),
           "source_2014_host_taxa_with_BOTH_EIV_N_and_MOIST":n,
           "source_2014_host_taxa_paired_trait_coverage":n/len(orig_hosts) if orig_hosts else None,
           "source_same_exact_host_NUTRI_mean":nmean,
           "source_same_exact_host_MOIST_mean":fmean,
           "original_20261115_mean_NUTRI_from_fixed_HNS_roster":oldN,
           "original_20261115_mean_MOIST_from_fixed_HNS_roster":oldF,
           "source_original_N_mean_difference_from_existing_roster":nmean-oldN if nmean is not None and oldN is not None else None,
           "source_original_F_mean_difference_from_existing_roster":fmean-oldF if fmean is not None and oldF is not None else None,
           "eligible_strict_same_HOST_N_F_and_accepted_butterfly":bool(eligible),
           "original_HOSTS_reference_is_global_not_real_local_2017_2022_usage":True})
    if len(rows)!=303:
        raise RuntimeError("Derived source butterfly roster rowcount not 303")
    passed=counted>=N_MIN and sum(v>0 for v in families.values())>=FAMILIES_MIN
    outcsv=a.output/"GESUI_ORIGINAL_303_FERRERPARIS_EIV_SAME_HOST_N_F_SOURCE_ONLY.csv"
    with outcsv.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    summary={"schema":SCHEMA,
       "status":"PASS_SAME_LARVAL_HOST_PLANT_N_F_EXACT_SOURCE_FOR_H7_OPTION" if passed
            else "STOP_SAME_LARVAL_HOST_PLANT_N_F_EXACT_SOURCE_4FAMILY40",
       "source_Ferrer_Paris_original":source,
       "source_EIV_original":source_eiv,
       "source_original_butterfly_303_accepted_taxonomy_SHA256":ROSTER_SHA,
       "source_original_303_rows":len(rows),
       "source_valid_strict_exact_butterflies_BOTH_N_F":counted,
       "source_valid_exact_butterflies_by_family":dict(families),
       "source_min40_accepted_bfly_and4_families":bool(passed),
       "source_no_alias_synonym_taxonomy_imputation":True,
       "source_nonliteral_genus_only_or_unknown_HOSTS_rejected":True,
       "source_EIV_non_numeric_x_not_NUTRI_0":True,
       "source_2014_global_HOSTS_not_2017_2022_locally_observed_actual_larval_host_use":True,
       "source_original_existing_roster_means_NOT_REPLACED":True,
       "source_derived_303_original_taxon_trait_CSV_SHA256":digest(outcsv.read_bytes()),
       "no_local_2022_WWTP_water_EQR_plant_butterfly_effect_yet_fitted":True}
    (a.output/"GESUI_ORIGINAL_HOST_N_F_SAME_PLANT_SOURCE_GATE.json").write_text(
       json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("GESUI_TRUE_SAME_HOST_N_F_SOURCE_GATE="+json.dumps({
       "status":summary["status"],
       "original_2014_host_rows":source["source_original_host_source_record_rows"],
       "original_EIV_paired_exact_plant_taxa":source_eiv["source_strict_EIV_bino_NF_paired_plants"],
       "original_303_accepted_trait_N_F_butterfly_score_source_n":counted,
       "source_family_support":dict(families),
       "original_roster_accepted_n":int(accepted.sum()),
       "derived_csv_SHA":summary["source_derived_303_original_taxon_trait_CSV_SHA256"],
       "true_local_HOSTS_usage_or_2022_wastewater_effect_NOT_measured":True
    },ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args())
