"""Offline math/source selection tests: all observational outcomes use external GBIF."""
import importlib.util
from pathlib import Path
import numpy as np

p=Path(__file__).resolve().parents[1]/"scripts/gesui_plantae_outfall_pilot.py"
s=importlib.util.spec_from_file_location("pilot",p)
mod=importlib.util.module_from_spec(s)
s.loader.exec_module(mod)

def test_wwtp_fraction_units():
    assert mod.water_fraction(86400,1)==.5
    assert mod.water_fraction(86400,3)==.25
    assert np.isnan(mod.water_fraction(1000,0))

def test_rarefaction_expected():
    assert mod.rarefy_expected([30],30)==1
    assert mod.rarefy_expected([15,15],30)==2
    assert mod.rarefy_expected([2,3],30) is None
    assert 1<mod.rarefy_expected([20,20],30)<2.01

def test_centers_not_overlapping_if_gt_two_radii():
    assert 4.5<mod.haversine(0,0,0,.045)<6
