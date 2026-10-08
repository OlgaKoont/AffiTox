#!/usr/bin/env python3
"""Publish a cleaned AffiTox Zenodo record 05 (no 20-check DynamicBind CSVs).

23214420 cannot use newversion (`files.enabled`) and its bucket is locked, so
this creates a sibling deposition, then marks 23214420 isObsoletedBy the new DOI
and retargets umbrella 22007993 hasPart.

Token from ZENODO_ACCESS_TOKEN or ~/.zenodo/token. Never prints the token.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_POST = Path(__file__).resolve().parent
ROOT = _POST.resolve().parents[1]
if str(_POST) not in sys.path:
    sys.path.insert(0, str(_POST))
from upload_zenodo import API, GITHUB, Client, load_token  # noqa: E402
from upload_zenodo_blocks import COMMUNITY, create_deposition, metadata_for, upload_file  # noqa: E402
from publish_zenodo_v1 import (  # noqa: E402
    UMBRELLA_DOI,
    UMBRELLA_ID,
    doi_rel,
    editable_metadata,
    enter_edit,
    get_dep,
    merge_related,
    put_and_publish,
)

PARENT_ID = 23214420
ZIP_PATH = ROOT / "release" / "zenodo" / "payloads" / "05_posebusters.zip"
REPORT = ROOT / "release" / "zenodo" / "payloads" / "zenodo_05_newversion.json"
TITLE = "AffiTox: PoseBusters tables (record 05)"


def description_html(parent_doi: str) -> str:
    return (
        "<p>Full-panel PoseBusters CSVs for Boltz-2, DynamicBind, GNINA, and QVina2 "
        "on the canonical 16-target AffiTox panel. PLAPT has no poses and is excluded.</p>"
        "<p>DynamicBind files are <code>posebusters_results_&lt;pdb&gt;_dynamicbind_new.csv</code> "
        "(21 boolean checks, including <code>internal_energy</code>). This deposit replaces "
        f"<a href='https://doi.org/{parent_doi}'>{parent_doi}</a>, which still contains "
        "superseded 20-check <code>posebusters_results_&lt;pdb&gt;_dynamicbind.csv</code> tables.</p>"
        f"<p>Index: <a href='https://doi.org/{UMBRELLA_DOI}'>{UMBRELLA_DOI}</a>. "
        f"Code: <a href='{GITHUB}'>{GITHUB}</a>.</p>"
    )


def restore_parent_if_editing(client: Client) -> None:
    dep = get_dep(client, PARENT_ID)
    if dep.get("state") != "inprogress":
        return
    print(f"republishing locked parent {PARENT_ID} (old zip kept) ...", flush=True)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    meta.pop("communities", None)
    put_and_publish(client, PARENT_ID, meta)
    print("  parent restored", flush=True)


def patch_parent_obsoleted(client: Client, new_doi: str) -> None:
    enter_edit(client, PARENT_ID)
    dep = get_dep(client, PARENT_ID)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    desc = str(meta.get("description") or "")
    banner = (
        "<p><strong>Superseded fileset.</strong> Use "
        f"<a href='https://doi.org/{new_doi}'>{new_doi}</a> for record 05 "
        "(no 20-check DynamicBind CSVs).</p>"
    )
    if new_doi not in desc:
        meta["description"] = banner + desc
    meta["notes"] = (
        f"Superseded AffiTox record 05 zip. Current files: {new_doi}."
    )
    related = merge_related(meta.get("related_identifiers") or [], doi_rel(new_doi, "isObsoletedBy"))
    meta["related_identifiers"] = related
    meta.pop("communities", None)
    print(f"marking {PARENT_ID} isObsoletedBy {new_doi} ...", flush=True)
    put_and_publish(client, PARENT_ID, meta)
    print("  parent metadata published", flush=True)


def patch_umbrella(client: Client, new_doi: str) -> None:
    enter_edit(client, UMBRELLA_ID)
    dep = get_dep(client, UMBRELLA_ID)
    meta = editable_metadata(dict(dep.get("metadata") or {}))
    old = f"10.5281/zenodo.{PARENT_ID}"
    desc = str(meta.get("description") or "").replace(old, new_doi)
    desc = desc.replace(str(PARENT_ID), new_doi.rsplit(".", 1)[-1])
    meta["description"] = desc
    related = []
    replaced = False
    for item in meta.get("related_identifiers") or []:
        ident = str(item.get("identifier") or "")
        if item.get("relation") == "hasPart" and ident in {old, f"https://doi.org/{old}"}:
            related.append(doi_rel(new_doi, "hasPart"))
            replaced = True
            continue
        related.append(item)
    if not replaced:
        related = merge_related(related, doi_rel(new_doi, "hasPart"))
    meta["related_identifiers"] = related
    meta.pop("communities", None)
    print(f"updating umbrella hasPart 05 -> {new_doi} ...", flush=True)
    put_and_publish(client, UMBRELLA_ID, meta)
    print(f"  umbrella published {UMBRELLA_DOI}", flush=True)


def main() -> int:
    if not ZIP_PATH.is_file():
        raise SystemExit(f"missing {ZIP_PATH}; pack 05 first")
    parent_doi = f"10.5281/zenodo.{PARENT_ID}"
    client = Client(load_token())
    restore_parent_if_editing(client)

    meta = metadata_for(TITLE, description_html(parent_doi))
    meta["related_identifiers"] = [
        {
            "identifier": GITHUB,
            "relation": "isDocumentedBy",
            "resource_type": "software",
        },
        doi_rel(UMBRELLA_DOI, "isPartOf"),
        doi_rel(parent_doi, "isNewVersionOf"),
    ]
    meta["notes"] = (
        "Current AffiTox v1.0.0 record 05. Replaces 10.5281/zenodo.23214420 "
        "(old zip still listed 20-check DynamicBind CSVs). "
        f"Community {COMMUNITY}."
    )
    print("creating new deposition ...", flush=True)
    dep = create_deposition(client, meta)
    dep_id = int(dep["id"])
    bucket = (dep.get("links") or {}).get("bucket")
    doi_pre = (dep.get("metadata") or {}).get("prereserve_doi", {}).get("doi")
    print(f"  draft {dep_id} prereserved {doi_pre}", flush=True)
    if not bucket:
        raise SystemExit(f"{dep_id} has no bucket")
    print(
        f"uploading {ZIP_PATH.name} ({ZIP_PATH.stat().st_size / 1024**2:.1f} MiB) ...",
        flush=True,
    )
    upload_file(client, bucket, ZIP_PATH)
    dep = get_dep(client, dep_id)
    if not (dep.get("files") or []):
        raise SystemExit(f"{dep_id} has no files after upload; refuse to publish")
    print(f"publishing {dep_id} ...", flush=True)
    published = put_and_publish(client, dep_id, editable_metadata(dict(dep.get("metadata") or {})))
    pub_meta = published.get("metadata") or {}
    doi = pub_meta.get("doi") or f"10.5281/zenodo.{dep_id}"
    conceptdoi = published.get("conceptdoi") or pub_meta.get("conceptdoi")
    print(f"  published {doi} (concept {conceptdoi})", flush=True)

    patch_parent_obsoleted(client, doi)
    patch_umbrella(client, doi)

    report = {
        "parent_id": PARENT_ID,
        "parent_doi": parent_doi,
        "new_id": dep_id,
        "doi": doi,
        "conceptdoi": conceptdoi,
        "mode": "new_deposition_replace_05",
        "zip": ZIP_PATH.name,
        "zip_bytes": ZIP_PATH.stat().st_size,
        "files": [
            {
                "filename": f.get("filename") or f.get("key"),
                "size": f.get("filesize") or f.get("size"),
            }
            for f in (get_dep(client, dep_id).get("files") or [])
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    print(f"wrote {REPORT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
