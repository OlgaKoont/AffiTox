#!/usr/bin/env python3
"""Patch editable metadata on published historical AffiTox Zenodo records.

Does not replace files. Token from ZENODO_ACCESS_TOKEN or ~/.zenodo/token.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from urllib.parse import urlencode

_POST = Path(__file__).resolve().parent
if str(_POST) not in sys.path:
    sys.path.insert(0, str(_POST))
from upload_zenodo import API, GITHUB, Client, load_token  # noqa: E402

HISTORICAL_IDS = (20825057, 20825059, 20825061, 20825063, 20825065, 20825067)
NOTE = (
    "HISTORICAL. This is not the AffiTox v1.0.0 data release. "
    "The old per-protein zips include 7awe/2z5x/3mjg and are superseded. "
    f"Current code: {GITHUB}. Do not cite this record as the paper dataset."
)
GITHUB_REL = {
    "identifier": GITHUB,
    "relation": "isDocumentedBy",
    "resource_type": "software",
}


def _json(client: Client, method: str, url: str, payload: dict | None = None) -> dict:
    body = None
    headers = {}
    if payload is not None:
        body = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    _, raw = client.request(method, url, data=body, headers=headers)
    return json.loads(raw) if raw else {}


def list_depositions(client: Client) -> list[dict]:
    out: list[dict] = []
    page = 1
    while True:
        q = urlencode({"size": 100, "page": page, "sort": "mostrecent"})
        batch = _json(client, "GET", f"{API}/deposit/depositions?{q}")
        if not isinstance(batch, list) or not batch:
            break
        out.extend(batch)
        if len(batch) < 100:
            break
        page += 1
        if page > 20:
            break
    return out


def summarize(dep: dict) -> dict:
    meta = dep.get("metadata") or {}
    return {
        "id": dep.get("id"),
        "state": dep.get("state"),
        "submitted": dep.get("submitted"),
        "title": meta.get("title"),
        "doi": (meta.get("prereserve_doi") or {}).get("doi") or meta.get("doi"),
        "access_right": meta.get("access_right"),
        "html": (dep.get("links") or {}).get("html")
        or (dep.get("links") or {}).get("latest_html"),
    }


def merge_related(existing: list, extra: dict) -> list:
    items = list(existing or [])
    key = (extra.get("identifier"), extra.get("relation"))
    if key not in {(i.get("identifier"), i.get("relation")) for i in items}:
        items.append(extra)
    return items


def _editable_metadata(meta: dict) -> dict:
    out = dict(meta)
    for key in ("doi", "prereserve_doi"):
        out.pop(key, None)
    return out


def patch_historical(client: Client, dep_id: int) -> str:
    dep = _json(client, "GET", f"{API}/deposit/depositions/{dep_id}")
    if not dep.get("submitted"):
        return f"{dep_id}: skip (not published, state={dep.get('state')})"
    try:
        _json(client, "POST", f"{API}/deposit/depositions/{dep_id}/actions/edit")
    except Exception as exc:
        if "400" not in str(exc) and "409" not in str(exc):
            raise
    dep = _json(client, "GET", f"{API}/deposit/depositions/{dep_id}")
    meta = _editable_metadata(dict(dep.get("metadata") or {}))
    desc = str(meta.get("description") or "")
    if "HISTORICAL" not in desc:
        meta["description"] = (
            "<p><strong>HISTORICAL — not AffiTox v1.0.0.</strong> "
            "Old per-protein docking zips (panel included 7awe/2z5x/3mjg). "
            "Do not use these files for the paper. "
            f'See <a href="{GITHUB}">{GITHUB}</a>.</p>'
            + desc
        )
    meta["notes"] = NOTE
    meta["related_identifiers"] = merge_related(
        meta.get("related_identifiers") or [], GITHUB_REL
    )
    _json(client, "PUT", f"{API}/deposit/depositions/{dep_id}", {"metadata": meta})
    _json(client, "POST", f"{API}/deposit/depositions/{dep_id}/actions/publish")
    return f"{dep_id}: patched and re-published metadata"


def main() -> int:
    client = Client(load_token())
    listed = [summarize(d) for d in list_depositions(client)]
    out_dir = Path(__file__).resolve().parents[2] / "release" / "zenodo" / "payloads"
    out_dir.mkdir(parents=True, exist_ok=True)
    list_path = out_dir / "zenodo_depositions_list.json"
    list_path.write_text(json.dumps(listed, indent=2) + "\n")
    print(f"listed {len(listed)} depositions -> {list_path}")
    for rec in listed:
        print(
            f"  {rec.get('id')} [{rec.get('state')}] "
            f"submitted={rec.get('submitted')} {rec.get('doi')} | {rec.get('title')}"
        )
    do_patch = "--patch-historical" in sys.argv
    if not do_patch:
        print("pass --patch-historical to update 20825057–067 notes")
        return 0
    for dep_id in HISTORICAL_IDS:
        try:
            print(patch_historical(client, dep_id), flush=True)
        except Exception as exc:
            print(f"{dep_id}: ERROR {exc}", flush=True)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
