"""Offline source integrity test suite; no external downloads."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
modpath = ROOT / "scripts" / "verify_migration_integrity.py"
spec = importlib.util.spec_from_file_location("gesui_verify", modpath)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class IntegrityTests(unittest.TestCase):
    def test_original_archive_all_80(self):
        output = mod.verify(ROOT, ROOT / "provenance" / "CHOCHO_PR122_FILE_MANIFEST.json")
        self.assertEqual(output["verified"], 80)
        self.assertTrue(output["all_match"])

    def test_sha_bytes(self):
        self.assertEqual(mod.git_blob_sha(b"hello\n"), "ce013625030ba8dba906f756967f9e9ca394464a")

if __name__ == "__main__":
    unittest.main()
