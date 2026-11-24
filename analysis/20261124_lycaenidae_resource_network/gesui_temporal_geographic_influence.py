"""Post-hoc fixed-loss influence diagnostic; no new model fitting."""
import csv, hashlib, pathlib, sys
from collections import defaultdict
src = pathlib.Path(__file__).with_name("gesui_strict_reference_losses.csv")
raw = src.read_bytes()
assert hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest() == "9cac456111e4660b2d3dd4ecc7b1cf0e8d7421be"
data = defaultdict(dict)
for r in csv.DictReader(raw.decode().splitlines()):
    if r["mode"] != "H6" or r["reference"] != "training_early_exact":
        continue
    key = (r["sample"], r["endpoint"])
    h = r["h6"]
    assert h not in data[key]
    data[key][h] = (r["fold"], float(r["mse0"]), float(r["mse1"]))
def gain(v):
    return 100 * sum(a-b for _,a,b in v) / sum(a for _,a,b in v)
out = csv.writer(sys.stdout)
out.writerow(["sample","plant","consumer","removed","H6","plant_gain","consumer_gain","contrast_pp","plant_improved","consumer_improved","positive_plant_folds"])
for sample,n in [("quality",77),("all",89)]:
    for pe,ce in [("plant_symmetric","consumer_symmetric"),("plant_early_consumer_fixed","consumer_early_plant_fixed")]:
        p,c=data[sample,pe],data[sample,ce]
        assert len(p)==len(c)==n and p.keys()==c.keys()
        assert all(p[h][0]==c[h][0] for h in p)
        ranked=sorted(p,key=lambda h: (-(p[h][1]-p[h][2]), h))
        for k in [0,1,5,10]:
            keep=[h for h in p if h not in set(ranked[:k])]
            pp,cc=[p[h] for h in keep],[c[h] for h in keep]
            gp,gc=gain(pp),gain(cc)
            folds={r[0] for r in pp}
            out.writerow([sample,pe,ce,k,len(keep),gp,gc,gp-gc,sum(a>b for _,a,b in pp),sum(a>b for _,a,b in cc),sum(gain([r for r in pp if r[0]==f])>0 for f in folds)])
