"""Synthetic source-design checks. Real spatial gate runs separately against pinned sources."""
import importlib.util
from pathlib import Path
import math

p=Path(__file__).resolve().parents[1]/"scripts/gesui_river_source_gate.py"
s=importlib.util.spec_from_file_location("source_gate",p)
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

class Reach:
    def __init__(self,nxt):
        self.NEXT_DOWN=nxt
        self.UPLAND_SKM=100.
        self.DIS_AV_CMS=2.

def test_source_river_fractions():
    assert m.fraction_treated(86400,1)==0.5
    assert m.fraction_treated(0,1)==0
    assert math.isnan(m.fraction_treated(1000,0))
    assert math.isnan(m.fraction_treated(-5,5))

def test_directed_topology():
    reaches={10:Reach(11),11:Reach(12),12:Reach(13),13:Reach(0)}
    up={11:[10],12:[11],13:[12]}
    assert m.river_chain(12,up,reaches)==(10,11,12,13)
    assert m.river_chain(11,up,reaches) is None

def test_distance_and_source_only_stratification():
    assert 110 < m.haversine(0,0,0,1) <112
    vals=[]
    for i in range(16):
        vals.append({"wwtp_id":i,"country":"DEU","effluent_fraction_proxy":i/100,
            "outfall_lat":50+i,"outfall_lon":10.})
    got=m.select_preoutcome(vals,max_per_country_stratum=2)
    assert len(got)==4
    assert {x["exposure_stratum"] for x in got}=={"low","high"}
