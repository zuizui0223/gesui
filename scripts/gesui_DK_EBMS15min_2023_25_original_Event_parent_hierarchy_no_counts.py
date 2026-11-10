#!/usr/bin/env python3
"""Source-only original parentEventID/child Event hierarchy before DK 15-minute count effects."""
from __future__ import annotations
import argparse,hashlib,io,json,math,urllib.request,zipfile,csv
import xml.etree.ElementTree as ET
from collections import defaultdict,Counter
from datetime import datetime,timezone
from pathlib import Path

NS="{http://rs.tdwg.org/dwc/text/}"
SOURCE="results/ebms_denmark_15minute_2023_2025_WWTP_SOURCE_ONLY/ORIGINAL_DK_EBMS15M_2023_25_EEA2022_203_OUTFALL_EVENT_ONLY_QA.json"
EEA="results/ebms_denmark_2023_2024_WWTP_source_only/DENMARK_EBMS2023_2024_EEA2022_WWTP_SOURCEONLY_EVENT_QC.json"
ZIP="https://ipt.brc.ac.uk/archive.do?r=ebms_denmark_15-minute_counts&v=1.3"
PINNED="c32be29200134fd2ee06b7dd34d76e1420bd6a7810574fb8ffb345c1451169cb"

def hash(b):return hashlib.sha256(b).hexdigest()
def num(t):
    try:
        x=float(t)
        return x if math.isfinite(x) else None
    except (TypeError,ValueError):return None
def hav(a,b,c,d):
    r=math.pi/180
    dl=(a-c)*r;dg=(b-d)*r
    v=math.sin(dl/2)**2+math.cos(a*r)*math.cos(c*r)*math.sin(dg/2)**2
    return 12742*math.asin(min(1,math.sqrt(v)))
def download():
    with urllib.request.urlopen(urllib.request.Request(ZIP,headers={"User-Agent":"gesui-ebms15m-parent-only-source/0.1"}),timeout=150) as f:b=f.read()
    if hash(b)!=PINNED:raise ValueError("Original source ZIP changed")
    return b
def run(root):
    sp=root/SOURCE;ss=json.loads(sp.read_text())
    epath=root/EEA;eea=json.loads(epath.read_text())
    assert ss["source_archive_receipt"]["original_ZIP_sha256"]==PINNED
    assert ss["original_Occurrence_extension_READ_in_this_stage"] is False
    outfalls=eea["original_2022_EEA_DK_qualifying_Fw_outfalls"]
    assert len(outfalls)==203
    originals=Counter()
    rootids=defaultdict(set)
    parents=defaultdict(set)
    parentsByWWTP=defaultdict(set)
    strictBANDS=defaultdict(Counter)
    raw=download()
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        xml=ET.fromstring(z.read("meta.xml"))
        core=xml.find(NS+"core")
        if core is None or core.attrib.get("rowType")!="http://rs.tdwg.org/dwc/terms/Event":raise ValueError("No Event core")
        terms={f.attrib["term"].rsplit("/",1)[-1]:int(f.attrib["index"]) for f in core.findall(NS+"field")}
        fieldID=core.find(NS+"id")
        idx=int(fieldID.attrib["index"]) if fieldID is not None else terms.get("eventID")
        if idx is None:raise ValueError("No original Event ID index")
        sourcefile=core.find(".//"+NS+"files/"+NS+"location")
        if sourcefile is None or not sourcefile.text:raise ValueError("No original Event table")
        sep=core.attrib.get("fieldsTerminatedBy",r"\t")
        if sep==r"\t":sep="\t"
        if sep!="\t":raise ValueError("Event delimiter unexpected")
        skip=int(core.attrib.get("ignoreHeaderLines","1"))
        with io.TextIOWrapper(z.open(sourcefile.text),encoding="utf-8",newline="") as f:
            r=csv.reader(f,delimiter="\t")
            for _ in range(skip):next(r,None)
            for row in r:
                originals["original_total_Event_rows"]+=1
                if len(row)<=max(idx,*terms.values()):
                    originals["row_short"]+=1;continue
                get=lambda n:row[terms[n]].strip() if n in terms else ""
                year=get("eventDate")[:4]
                if year not in ("2023","2024","2025"):continue
                latitude=num(get("decimalLatitude"));longitude=num(get("decimalLongitude"));unc=num(get("coordinateUncertaintyInMeters"))
                if latitude is None or longitude is None or unc is None or unc>500 or not(54<=latitude<=58 and 7<=longitude<=16):continue
                yr=int(year)
                id=row[idx].strip();parent=get("parentEventID")
                if not id:raise ValueError("Original Event ID absent")
                group=parent if parent else id
                if parent:originals["child_Event_rows_in_strict_500m_source"]+=1
                else:originals["root_Event_rows_in_strict_500m_source"]+=1
                nearest=min(((hav(latitude,longitude,float(o["lat"]),float(o["lon"])),o["outfall_id_SHA256"]) for o in outfalls),key=lambda a:(a[0],a[1]))
                distance,wwtp=nearest
                band="near2" if distance<=2 else "near5" if distance<=5 else "remote10" if distance>=10 else "middle"
                key=(yr,band)
                strictBANDS[key]["all_event_rows"]+=1
                parents[key].add(group)
                rootids[key].add(id if not parent else "")
                if band in ("near2","near5"):
                    parentsByWWTP[(yr,wwtp)].add(group)
                    if band=="near2":
                        strictBANDS[(yr,"near5cumulative")]["all_event_rows"]+=1
                        parents[(yr,"near5cumulative")].add(group)
                    else:
                        strictBANDS[(yr,"near5cumulative")]["all_event_rows"]+=1
                        parents[(yr,"near5cumulative")].add(group)
    if originals["original_total_Event_rows"]!=44884:raise ValueError("Archive Event row count mismatch")
    result={}
    for yr in (2023,2024,2025):
        d={}
        for band in ("near2","near5","near5cumulative","middle","remote10"):
            d[band]={"original_strict_quality_Event_rows":strictBANDS[(yr,band)]["all_event_rows"],
                     "distinct_immediate_parent_or_self_survey_groups":len(parents[(yr,band)])}
        active={k[1]:len(v) for k,v in parentsByWWTP.items() if k[0]==yr}
        d["near_source_systems_with_ge2_parent_groups"]=sum(c>=2 for c in active.values())
        d["near_source_systems_any_parent_group"]=len(active)
        d["parent_source_independence_gate_prespecified"]=bool(
            d["near5cumulative"]["distinct_immediate_parent_or_self_survey_groups"]>=50
            and d["remote10"]["distinct_immediate_parent_or_self_survey_groups"]>=100
            and d["near_source_systems_with_ge2_parent_groups"]>=10)
        expected=ss["source_grouped_year_and_accuracy_support"]["strict500"][str(yr)]
        assert d["near5cumulative"]["original_strict_quality_Event_rows"]==expected["outfall_near5_cumulative_Events"]
        assert d["remote10"]["original_strict_quality_Event_rows"]==expected["registered_source_remote10_Events"]
        result[str(yr)]=d
    pass2=all(result[str(y)]["parent_source_independence_gate_prespecified"] for y in (2024,2025))
    out={
        "schema":"gesui_DK_2023_25_EBMS15m_original_Event_parent_hierarchy_source_only_v01",
        "utc":datetime.now(timezone.utc).isoformat(),
        "original_15min_2023_2025_Event_source_QA_SHA256":hash(sp.read_bytes()),
        "original_fixed_2022_DK_EEA_source_SHA256":hash(epath.read_bytes()),
        "original_archive_sha256":PINNED,
        "original_Event_only_hierarchy_qc":dict(originals),
        "source_year_parent_group_count_and_WWTP_system_eligibility":result,
        "2024_2025_both_parent_cluster_support_gate_pass":pass2,
        "original_taxonomic_Occurrence_extension_never_opened_here":True,
        "parentEventID_is_not_independent_butterfly_population_or_necessarily_one_15min_observation":True,
        "a_single_parent_group_is_counted_once_for_observer_effort_cluster_support":True,
    }
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);args=p.parse_args()
    root=Path(__file__).resolve().parents[1];out=run(root)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,sort_keys=True,indent=2)+"\n")
    print(json.dumps({"source":out["original_Event_only_hierarchy_qc"],"yearwise":out["source_year_parent_group_count_and_WWTP_system_eligibility"],"gate":out["2024_2025_both_parent_cluster_support_gate_pass"]},indent=2))
if __name__=="__main__":main()
