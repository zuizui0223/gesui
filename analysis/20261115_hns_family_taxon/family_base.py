"""Import the frozen PR94 family decomposition by exact path, not shared module name.

Multiple historical source packages call their scripts analyse.py. Avoid a
circular import when this newer analysis is loaded from sibling sensitivity.py.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

_SOURCE = (Path(__file__).resolve().parents[1] /
           "20261115_hns_family_decomposition" / "analyse.py")
_SPEC = spec_from_file_location("gesui_frozen_pr94_base", _SOURCE)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("Original PR94 analysis source unavailable")
_BASE = module_from_spec(_SPEC)
_SPEC.loader.exec_module(_BASE)

verify_source = _BASE.verify_source
design = _BASE.design
COMMON = _BASE.COMMON
