"""Provenance for the frozen inputs behind the index.

The M1 contract requires every fact to say *which* source record asserted it
and *when* that record was retrieved.  The graph already carries
``database / record_id / release / field / locator``; this module supplies the
missing ``retrieved_at``.

Ordering of trust, highest first:

1. ``retrieved_at`` recorded in ``data/manifest.json`` when the bytes were
   downloaded (``status: "download"``);
2. the file's own mtime, which is when the bytes landed on this disk
   (``status: "file_mtime"``);
3. nothing, and we say so instead of inventing a date
   (``status: "not_reported"``).
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from typing import Dict, Iterable, Optional, Tuple

UTC = timezone.utc


def iso_utc(timestamp: float) -> str:
    """Render a POSIX timestamp as a second-precision UTC ISO-8601 string."""
    return datetime.fromtimestamp(timestamp, UTC).strftime("%Y-%m-%dT%H:%M:%SZ")


def load_manifest(path: str) -> Dict[str, object]:
    if not os.path.exists(path):
        return {"files": {}}
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def retrieved_at(path: str, record: Optional[Dict[str, object]] = None) -> Tuple[Optional[str], str]:
    """Return ``(retrieved_at, status)`` for one input file.

    ``status`` is one of ``download`` / ``file_mtime`` / ``not_reported``.
    """
    if record and record.get("retrieved_at"):
        return str(record["retrieved_at"]), str(record.get("retrieved_at_source") or "download")
    if os.path.exists(path):
        return iso_utc(os.path.getmtime(path)), "file_mtime"
    return None, "not_reported"


def collect(entries: Iterable[Tuple[str, str]], manifest_path: str, root: str) -> Dict[str, Dict[str, object]]:
    """Build the ``source_key -> provenance`` map stored in the index meta.

    ``entries`` are ``(source_key, absolute_path)`` pairs; ``source_key`` is the
    value the graph uses to name that input (a basename for the three main
    files, ``biomart/<file>`` for the BioMart dumps).
    """
    files = load_manifest(manifest_path).get("files") or {}
    collected: Dict[str, Dict[str, object]] = {}
    for key, path in entries:
        rel = os.path.relpath(path, root).replace(os.sep, "/")
        record = files.get(rel) or {}
        when, status = retrieved_at(path, record)
        sha = record.get("sha256")
        if not sha and os.path.exists(path):
            # derived artefacts (e.g. xref_enrichment.json) are not downloads
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for block in iter(lambda: handle.read(1 << 20), b""):
                    digest.update(block)
            sha = digest.hexdigest()
        collected[key] = {
            "path": rel,
            "sha256": sha,
            "retrieved_at": when,
            "retrieved_at_status": status,
        }
    return collected


def stamp() -> Dict[str, str]:
    """Build-time identity: when this index was produced."""
    now = datetime.now(UTC).replace(microsecond=0)
    return {"built_at": now.strftime("%Y-%m-%dT%H:%M:%SZ")}
