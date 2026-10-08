#!/usr/bin/env python3
"""Zip AffiTox Zenodo information blocks. Record 03 is sharded (<45 GiB)."""

from __future__ import annotations

import json
import os
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("TOXAFFINITY_ROOT", Path(__file__).resolve().parents[2]))
ZENODO = ROOT / "release" / "zenodo"
OUT = Path(os.environ.get("OUT_DIR", ZENODO / "payloads"))
TARGETS = (
    "1g5m", "2v5z", "3eyg", "3jy9", "3lxk", "11ue", "4ase", "4f65",
    "4tz4", "4zau", "5jkv", "5mo4", "6gqj", "6jok", "5lf3", "7kk3",
)
BUCKETS = [
    ("03a", ["3eyg", "4f65", "5lf3", "2v5z"]),
    ("03b", ["3jy9", "3lxk"]),
    ("03c", ["5mo4", "6gqj", "4zau"]),
    ("03d", ["7kk3", "5jkv", "4tz4"]),
    ("03e", ["4ase", "6jok", "11ue", "1g5m"]),
]
LIMIT = 45 * 1024**3


def zip_tree(src: Path, dest: Path, arc_root: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        dest.unlink()
    print(f"zip {src} -> {dest}", flush=True)
    t0 = time.time()
    with zipfile.ZipFile(dest, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for path in src.rglob("*"):
            if path.is_file():
                zf.write(path, Path(arc_root) / path.relative_to(src))
    print(f"  {dest.name} {dest.stat().st_size / 1024**3:.2f} GiB in {time.time()-t0:.0f}s", flush=True)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "03").mkdir(exist_ok=True)
    print(f"root {ROOT}\nout {OUT}\nstarted {time.strftime('%Y-%m-%dT%H:%M:%S')}", flush=True)
    for name in ("01_raw", "02_prepared", "04_merged", "05_posebusters", "06_metrics_figures"):
        zip_tree(ZENODO / name, OUT / f"{name}.zip", name)
    for pdb in TARGETS:
        src = ROOT / "results" / pdb
        dest = OUT / "03" / f"results_{pdb}.zip"
        if src.is_dir():
            zip_tree(src, dest, f"results/{pdb}")
        else:
            print(f"WARNING: missing {src}", file=sys.stderr, flush=True)
    manifest = {"records": []}
    for name, pdbs in BUCKETS:
        files = []
        total = 0
        for pdb in pdbs:
            z = OUT / "03" / f"results_{pdb}.zip"
            size = z.stat().st_size if z.is_file() else 0
            files.append({"pdb": pdb, "zip": f"03/results_{pdb}.zip", "bytes": size})
            total += size
        if total > LIMIT:
            raise SystemExit(f"{name} is {total/1024**3:.1f} GiB > 45 GiB")
        manifest["records"].append({"id": name, "files": files, "bytes": total})
        print(f"{name}: {total/1024**3:.2f} GiB", flush=True)
    (OUT / "03_buckets.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"finished {time.strftime('%Y-%m-%dT%H:%M:%S')}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
