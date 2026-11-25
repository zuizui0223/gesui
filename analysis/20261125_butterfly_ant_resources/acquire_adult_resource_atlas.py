#!/usr/bin/env python3
"""Extract only exact original BasinATLAS v1.0 native level07 DBF; no response models."""
import argparse,collections,csv,hashlib,io,json,pathlib,re,zipfile,urllib.parse
import requests
from dbfread import DBF
URL="https://figshare.com/ndownloader/files/20087237"
MEMBER="BasinATLAS_v10_shp/BasinATLAS_v10_lev07.dbf"
SHA="1c73db5c07cb02ee6d10b963a0fc9fa4321fab55ff43721e8d6cd5c1054fb113"
PANEL="b0babfa63a96f9b7eb5a83122b08b8664edab6e8c6c28d2faf4ff024ca84c39b"
ATLAS="99f42768ac0d215fca5662de449f97fba598bebab9b147288ef99e19db3633b7"
OLD=["crp_pc_sse","crp_pc_use","dis_m3_pyr","glc_pc_s22","glc_pc_u22","hft_ix_s09","urb_pc_sse","urb_pc_use"]
WET=["wet_pc_s%02d"%i for i in range(1,10)]+["wet_pc_sg1","wet_pc_sg2"]
NEW=WET+["pst_pc_sse","for_pc_sse","ele_mt_sav","tmp_dc_syr","pre_mm_syr","pop_ct_ssu","pop_ct_usu","ppd_pk_sav","ari_ix_sav","swc_pc_syr"]
def sha(b):return hashlib.sha256(b).hexdigest()
class RemoteZip(io.RawIOBase):
 def __init__(self):
  self.pos=0;self.transferred=0;self.calls=[];self.s=requests.Session()
  self.url=URL
  b,total=self.get(0,0);self.size=total
 def get(self,a,b):
  assert 0<=a<=b and b-a+1<=40*1024**2
  with self.s.get(self.url,headers={"Range":f"bytes={a}-{b}","Accept-Encoding":"identity"},stream=True,timeout=(20,120)) as r:
   if r.status_code!=206:raise RuntimeError("expected range 206, got "+str(r.status_code)+" "+r.url.split("?")[0])
   m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",r.headers.get("Content-Range",""));assert m and int(m[1])==a and int(m[2])==b,r.headers
   u=urllib.parse.urlparse(r.url);assert u.scheme=="https" and (u.hostname.endswith(".amazonaws.com") or u.hostname.endswith(".figshare.com") or u.hostname=="figshare.com" or u.hostname=="44c1cec316ec56bca76dfd65d2882b20.r2.cloudflarestorage.com"),u.hostname
   v=r.content;assert len(v)==b-a+1;self.calls.append([a,b,len(v)]);self.transferred+=len(v)
   assert self.transferred<=100*1024**2
   return v,int(m[3])
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,o,w=0):
  self.pos=o if w==0 else self.pos+o if w==1 else self.size+o
  assert self.pos>=0
  return self.pos
 def read(self,n=-1):
  if n==0 or self.pos>=self.size:return b""
  end=self.size-1 if n<0 else min(self.size-1,self.pos+n-1)
  v,t=self.get(self.pos,end);assert t==self.size;self.pos=end+1;return v
def main(a):
 out=a.out;out.mkdir(parents=True,exist_ok=True)
 b=a.panel.read_bytes();assert sha(b)==PANEL
 z=zipfile.ZipFile(io.BytesIO(b))
 rows=list(csv.DictReader(io.TextIOWrapper(z.open("EU234_H7_ORIGINAL_STUDY_NATIVE_BASINPERIOD_EQR_PLANT_BUTTERFLY_SOURCE_COHORT.csv"),encoding="utf8")))
 ids={r["HYBAS_ID_L07"] for r in rows};assert len(rows)==454 and len(ids)==232
 old=None
 if a.atlas:
  b=a.atlas.read_bytes();assert sha(b)==ATLAS
  z=zipfile.ZipFile(io.BytesIO(b));old=json.loads(z.read("HYDROATLAS_V1_L07_1219_H7_BASIN_RAW_LANDUSE_ATTRIBUTES_SOURCE_ONLY.json"))["original_1219_H7_covariate_by_real_ID"]
  assert len(old)==1219 and ids<=set(old)
 r=RemoteZip()
 with zipfile.ZipFile(r) as z:
  m=z.getinfo(MEMBER);assert m.file_size==80137382 and m.CRC==int("ef18ec36",16) and m.compress_size==15674749
  data=z.read(m);assert sha(data)==SHA
 p=out/"temporary_original.dbf";p.write_bytes(data)
 try:
  t=DBF(str(p),load=False);cols=t.field_names
  assert "HYBAS_ID" in cols and all(c in cols for c in OLD+WET+["pst_pc_sse"]),cols
  available=[c for c in NEW if c in cols]
  keep=OLD+available;got={};n=0;replay=[]
  for x in t:
   n+=1;k=str(x["HYBAS_ID"])
   replay.append({"h":k,**{c:x[c] for c in OLD}})
   if k in ids:
    assert k not in got
    got[k]={c:x[c] for c in keep}
   if old is not None and k in old:
    assert all(x[c]==old[k][c] for c in OLD),(k,x,old[k])
  assert n==57646 and set(got)==ids and len({v["h"] for v in replay})==n
  with (out/"native_old_field_replay.csv").open("w") as f:
   w=csv.DictWriter(f,fieldnames=["h"]+OLD);w.writeheader();w.writerows(replay)
  fields={c:{"known":sum(got[k][c] is not None and got[k][c]!=-9999 for k in ids),"minimum":min(got[k][c] for k in ids if got[k][c] is not None and got[k][c]!=-9999),"maximum":max(got[k][c] for k in ids if got[k][c] is not None and got[k][c]!=-9999)} for c in available if any(got[k][c] is not None and got[k][c]!=-9999 for k in ids)}
  for c in WET+["pst_pc_sse"]:
   assert fields[c]["known"]==232 and 0<=fields[c]["minimum"]<=fields[c]["maximum"]<=100,(c,fields[c])
  assert all(got[k]["wet_pc_sg2"]<=got[k]["wet_pc_sg1"] for k in ids)
  with (out/"basin_resource_fields.csv").open("w") as f:
   w=csv.DictWriter(f,fieldnames=["h"]+keep);w.writeheader();w.writerows({"h":k,**got[k]} for k in sorted(ids))
  meta={"status":"SOURCE_PASS","n_h7":len(ids),"original_dbf_sha256":sha(data),"original_dbf_crc32":format(m.CRC,"08x"),"n_global_dbf":n,"old_replayed_ids":len(old) if old else 0,"local_prior_1219_replay_required":old is None,"old_replayed_fields":len(OLD),"all_native_fields":cols,"new_fields":fields,"missing_requested_fields":sorted(set(NEW)-set(available)),"transferred_bytes":r.transferred,"ranges":r.calls,"no_ecological_outcomes":True}
  (out/"receipt.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
  print("RESOURCE SOURCE "+json.dumps({k:v for k,v in meta.items() if k!="all_native_fields"}),flush=True)
 finally:p.unlink(missing_ok=True)
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--panel",type=pathlib.Path,required=True);p.add_argument("--atlas",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path,required=True);main(p.parse_args())
