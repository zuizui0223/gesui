#!/usr/bin/env python3
"""2023–25 GBIF source acquisition; no biological response models."""
import argparse, csv, hashlib, json, os, pathlib, concurrent.futures, collections
import original_facet_task as orig
ROOT=pathlib.Path(__file__).resolve().parents[2]
PANEL=ROOT/"analysis/20261124_lycaenidae_resource_network/gesui_rebuilt_network.csv"
PERIOD="2023-2025"
def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()
def roster():
    assert blob(PANEL)=="3646cebf16691b81aca3d80ab5d2613d47d189fa"
    assert blob(pathlib.Path(orig.__file__))=="38830b600bf4c78e895b603a85a276ba2d551132"
    by=collections.defaultdict(dict)
    for r in csv.DictReader(PANEL.open()):
        assert r["period"] not in by[r["h"]]
        by[r["h"]][r["period"]]=(r["h6"],float(r["minimum_species_coverage"]))
    assert len(by)==232 and sum(map(len,by.values()))==454
    out={}
    for h,v in by.items():
        if set(v)=={"2017-2019","2020-2022"}:
            assert len({x[0] for x in v.values()})==1
            out[h]={"h6":next(iter(v.values()))[0],"prior_quality":all(x[1]>=.8 for x in v.values())}
    assert len(out)==222 and len({x["h6"] for x in out.values()})==89
    assert sum(x["prior_quality"] for x in out.values())==173
    return out
def summary(rows, scope):
    by=collections.defaultdict(dict)
    for r in rows:
        h=str(r["HYBAS_ID_L07"]);g=r["group"]
        assert h in scope and g not in by[h] and r["period"]==PERIOD
        by[h][g]=r
    assert set(by)==set(scope)
    details=[]
    for h,rs in sorted(by.items()):
        assert set(rs)==set(("Plantae",)+orig.FAMILIES)
        complete=all(r["source_status"]=="ok" for r in rs.values())
        p=b=cov=None
        if complete:
            p=rs["Plantae"]["source_facetted_occurrence_count"]
            b=sum(rs[g]["source_facetted_occurrence_count"] for g in orig.FAMILIES)
            cov=min(r["source_facet_coverage"] or 0 for r in rs.values())
        eligible=complete and p>=50 and b>=20
        details.append({"h":h,**scope[h],"complete":complete,"plant_facetted_n":p,"butterfly_facetted_n":b,"minimum_coverage":cov,"rarefaction_support":eligible,"three_period_quality":bool(eligible and cov>=.8 and scope[h]["prior_quality"])})
    return details
def write(path,obj):
    path.write_text(json.dumps(obj,sort_keys=True,indent=2)+"\n")
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mode",choices=["pilot","full","aggregate"],required=True)
    p.add_argument("--shard",type=int,default=0,choices=range(4))
    p.add_argument("--out",type=pathlib.Path,required=True)
    p.add_argument("--inputs",type=pathlib.Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    scope=roster()
    if a.mode=="aggregate":
        rows=[]
        for s in range(4):
            fs=list(a.inputs.rglob(f"third_full_{s}.json"))
            assert len(fs)==1
            rows+=json.loads(fs[0].read_text())
        assert len(rows)==1332
        details=summary(rows,scope)
        write(a.out/"third_period_support.json",details)
        primary=[x for x in details if x["three_period_quality"]]
        eligible=[x for x in details if x["rarefaction_support"]]
        qa={"expected_h7":222,"requested_facets":1332,"complete_h7":sum(x["complete"] for x in details),"eligible_h7":len(eligible),"eligible_h6":len({x["h6"] for x in eligible}),"quality_h7":len(primary),"quality_h6":len({x["h6"] for x in primary}),"note":"Source opportunity only. No host outcome or fitted model. API missing is not biological absence."}
        write(a.out/"third_period_gate.json",qa);print(json.dumps(qa));return
    basins,families=orig.read_frozen()
    by={str(z["HYBAS_ID_L07"]):z for z in basins}
    assert set(scope)<=set(by)
    orig.PERIODS={PERIOD:"2023,2025"}
    if a.mode=="pilot":
        ids=sorted(scope,key=lambda h:hashlib.sha256(("third-period-source-pilot-v1|"+h).encode()).hexdigest())[:8]
    else:
        ids=sorted(h for h in scope if int(hashlib.sha256(h.encode()).hexdigest()[:8],16)%4==a.shard)
    work=[(by[h],PERIOD,g,families) for h in ids for g in ("Plantae",)+orig.FAMILIES]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        rows=list(ex.map(lambda q:orig.task(*q),work))
    rows.sort(key=lambda r:(str(r["HYBAS_ID_L07"]),r["group"]))
    details=summary(rows,{h:scope[h] for h in ids})
    stem=f"third_{a.mode}_{a.shard}"
    write(a.out/(stem+".json"),rows)
    write(a.out/(stem+"_support.json"),details)
    missing=sum(r["source_status"]!="ok" for r in rows)
    print(json.dumps({"mode":a.mode,"n_h7":len(ids),"requested":len(work),"missing":missing,"source_only":True}),flush=True)
    if a.mode=="pilot" and missing:raise SystemExit("SOURCE_STOP: pilot incomplete; no full request")
if __name__=="__main__": main()
