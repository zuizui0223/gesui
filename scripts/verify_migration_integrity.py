#!/usr/bin/env python3
"""Offline verification of all 80 copied chocho PR122 Git blobs."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

def git_blob_sha(data: bytes) -> str:
    header = b"blob " + str(len(data)).encode("ascii") + b"\0"
    return hashlib.sha1(header + data).hexdigest()

def verify(root: Path, manifest: Path) -> dict:
    root = root.resolve()
    index = json.loads(manifest.read_text(encoding="utf-8"))
    rows = index["files"]
    if index["file_count"] != len(rows) or len(rows) != 80:
        raise ValueError("Expected the fixed 80-file PR122 manifest")
    paths = set()
    for item in rows:
        rel = Path(item["archived_path"])
        if rel.is_absolute() or ".." in rel.parts or not rel.parts or rel.parts[0] != "archive":
            raise ValueError("Unsafe archive path: " + str(rel))
        if str(rel) in paths:
            raise ValueError("Duplicate archive path: " + str(rel))
        paths.add(str(rel))
        path = (root / rel).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError("Missing or escaping archive file: " + str(rel))
        sha = git_blob_sha(path.read_bytes())
        if sha != item["git_blob_sha"]:
            raise ValueError("Source blob mismatch: " + str(rel))
    return {"verified": len(rows), "source_commit": index["source_commit"], "all_match": True}

def main() -> None:
    cli = argparse.ArgumentParser()
    cli.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    cli.add_argument("--manifest", type=Path)
    args = cli.parse_args()
    source = args.manifest or args.root / "provenance" / "CHOCHO_PR122_FILE_MANIFEST.json"
    print(json.dumps(verify(args.root, source), ensure_ascii=False, sort_keys=True))

if __name__ == "__main__":
    main()
