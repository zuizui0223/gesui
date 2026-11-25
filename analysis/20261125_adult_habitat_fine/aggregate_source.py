#!/usr/bin/env python3
"""Combine fixed shard source outputs, never fit butterfly outcomes."""
import argparse,hashlib,json,pathlib,csv,collections
import pandas as pd,numpy as np
CLASSES=[10,20,30,40,50,60,70,80,90,95,100]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(root,out):
 out.mkdir(parents=True,exist_ok=True)
 folders=[root/f"wc-source-{i}" for i in range(4)]
 assert all((p/"receipt.json").is_file() for p in folders),"missing shard"
 roster=json.loads((folders[0]/"tile_roster.json").read_text());assert len(roster)==44
 areas=pd.read_csv(folders[0]/"basin_areas.csv",dtype={"h":str});assert len(areas)==232 and areas.h.is_unique
 source=[];parts=[];input_receipts={}
 for i,p in enumerate(folders):
  assert json.loads((p/"tile_roster.json").read_text())==roster
  assert sha(p/"basin_areas.csv")==sha(folders[0]/"basin_areas.csv")
  r=json.loads((p/"receipt.json").read_text());assert r["shard"]==i
  assert {x["tile"] for x in r["sources"]}=={x["tile"] for x in roster if x["shard"]==i}
  source.extend(r["sources"])
  if (p/"tile_basin_samples.csv").is_file():parts.append(pd.read_csv(p/"tile_basin_samples.csv",dtype={"h":str}))
  input_receipts[str(i)]={x.name:sha(x) for x in p.iterdir() if x.is_file()}
 (out/"tile_source_inventory.json").write_text(json.dumps(source,indent=2)+"\n")
 expected={(h,t["tile"],p) for t in roster for h in t["basins"] for p in (5,15)}
 allok=all(x["status"]=="OK" for x in source)
 r={"status":"SOURCE_STOP","source_tiles":len(source),"ok_tiles":sum(x["status"]=="OK" for x in source),"expected_tile_basin_phase_rows":len(expected),"shard_receipts":input_receipts,"no_butterfly_outcomes":True}
 try:
  f=pd.concat(parts,ignore_index=True);assert not f.duplicated(["h","tile","phase"]).any()
  actual=set(f[["h","tile","phase"]].itertuples(index=False,name=None));assert actual<=expected
  f.to_csv(out/"tile_basin_samples.csv",index=False)
  r["actual_rows"]=len(f)
  missing=sorted(expected-actual);r["missing_rows"]=len(missing)
  (out/"missing_tile_basin_phase.json").write_text(json.dumps(missing))
  assert allok and actual==expected,"missing native source"
  cols=[c for c in f if c not in ["h","tile","phase"]]
  b=f.groupby(["h","phase"],as_index=False)[cols].sum().merge(areas,on="h",validate="many_to_one")
  assert len(b)==464 and b.h.nunique()==232
  b["area_rel_error"]=abs(b.sample_area_km2/b.area_km2-1)
  b["valid_fraction"]=b.valid_area_km2/b.sample_area_km2
  b["neighbor_valid_fraction"]=b.neighbor_valid_area_km2/b.sample_area_km2
  b["eligible"]=(b.area_rel_error<=.02)&(b.valid_fraction>=.95)&(b.sample_n>=1000)
  b["edge_eligible"]=b.eligible&(b.neighbor_valid_fraction>=.95)
  for c in CLASSES:b[f"class_{c}_fraction"]=b[f"class_{c}_km2"]/b.valid_area_km2
  b["open_fraction"]=b.class_30_fraction+b.class_60_fraction
  b["wet_water_fraction"]=b.class_80_fraction+b.class_90_fraction
  b["wet_edge_fraction"]=b.wet_edge_km2/b.neighbor_valid_area_km2
  b["woody_edge_fraction"]=b.woody_edge_km2/b.neighbor_valid_area_km2
  b["wet_edge_conditional_open"]=b.wet_edge_km2/b.open_neighbor_valid_area_km2.replace(0,np.nan)
  b["woody_edge_conditional_open"]=b.woody_edge_km2/b.open_neighbor_valid_area_km2.replace(0,np.nan)
  assert np.allclose(b[[f"class_{c}_fraction" for c in CLASSES]].sum(axis=1),1,atol=1e-9)
  b.to_csv(out/"basin_phase_metrics.csv",index=False)
  g=b.groupby("h"); e=g.eligible.all();ee=g.edge_eligible.all()
  metrics=["open_fraction","wet_water_fraction","class_30_fraction","class_60_fraction","class_80_fraction","class_90_fraction","wet_edge_fraction","woody_edge_fraction","wet_edge_conditional_open","woody_edge_conditional_open"]
  s=g[metrics].mean();d=g[metrics].max()-g[metrics].min()
  s.loc[~e,metrics]=np.nan
  s.loc[~ee,[x for x in metrics if "edge" in x]]=np.nan
  s["eligible"]=e;s["edge_eligible"]=ee
  s.to_csv(out/"basin_habitat_metrics.csv")
  d.to_csv(out/"phase_absolute_difference.csv")
  r.update(status="SOURCE_OVERLAY_COMPLETE",eligible_basins=int(e.sum()),edge_eligible_basins=int(ee.sum()),area_error_max=float(b.area_rel_error.max()),valid_fraction_min=float(b.valid_fraction.min()),sample_n_min=int(b.sample_n.min()),metrics={m:{"n":int(s[m].notna().sum()),"nonzero":int((s[m]>0).sum()),"median":float(s[m].median()) if s[m].notna().any() else None,"phase_difference_median":float(d[m].median()) if d[m].notna().any() else None} for m in metrics})
 except Exception as ex:r["error"]=type(ex).__name__+": "+str(ex)
 finally:
  (out/"receipt.json").write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
  print("WC SOURCE "+json.dumps({k:v for k,v in r.items() if k!="shard_receipts"}),flush=True)
  assert r["status"]=="SOURCE_OVERLAY_COMPLETE",r["status"]
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--root",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);x=p.parse_args();main(x.root,x.out)
