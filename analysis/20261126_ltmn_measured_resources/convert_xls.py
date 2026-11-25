#!/usr/bin/env python3
"""Convert the two original verified legacy XLS tables to UTF-8 CSV without editing the workbooks."""
import argparse,csv,hashlib,json,pathlib,math
import xlrd
SOURCES={
 "butterfly_index":("ff0e3bafed38a9c657635c3dd0cc9a197d9909da35e357c8d0dec6c3131fa8df",["UKBMS_SITENO","UKBMS_SITENAME","YEAR","COMMON_NAME","SITE_INDEX"],31951),
 "butterfly_visits":("c23945feb1568b4b60ac139136c485c234699710eae395a104e0be92a958a84d",["UKBMS_SITENO","UKBMS_SITENAME","DAY","MONTH","YEAR","START_TIME","RECORDER_NAME - deleted","ONLINE_RECORDER_NAME - deleted"],23877),
}
def main(src,out):
 out.mkdir(parents=True,exist_ok=True);r={}
 for key,(expected,columns,n) in SOURCES.items():
  p=src/(key+".xls");b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==expected,key
  w=xlrd.open_workbook(file_contents=b);s=w.sheet_by_index(0)
  assert s.nrows==n and s.row_values(0)==columns,(key,s.nrows,s.row_values(0))
  # Retain exact numeric values (incl. explicit zeros), blank cells and source row order.
  use=range(5) if key=="butterfly_index" else range(6)
  dest=out/(key+".csv")
  with dest.open("w",newline="",encoding="utf-8") as f:
   wr=csv.writer(f);wr.writerow([columns[i] for i in use])
   for i in range(1,s.nrows):
    vals=[]
    for j in use:
     c=s.cell(i,j);v=c.value
     if c.ctype in (xlrd.XL_CELL_EMPTY,xlrd.XL_CELL_BLANK):v=""
     elif c.ctype==xlrd.XL_CELL_ERROR:raise ValueError((key,i,j,"Excel error"))
     vals.append(v)
    wr.writerow(vals)
  r[key]={"source_sha256":expected,"source_rows_including_header":s.nrows,"csv_rows":n-1,"csv_sha256":hashlib.sha256(dest.read_bytes()).hexdigest(),"columns":[columns[i] for i in use]}
 (out/"conversion_receipt.json").write_text(json.dumps(r,indent=2)+"\n");print(json.dumps(r))
if __name__=="__main__":
 a=argparse.ArgumentParser();a.add_argument("--src",type=pathlib.Path,required=True);a.add_argument("--out",type=pathlib.Path,required=True);p=a.parse_args();main(p.src,p.out)
