"""Append-only trace of adapter calls.

The contract says every run should be replayable: what was asked, what was
answered, in what order.  This module writes one JSON line per operation when
``ENSEMBL_TRACE`` points at a file; otherwise it costs a single environment
lookup and nothing is written.

    ENSEMBL_TRACE=runs/trace.jsonl python3 contract/conformance.py
"""

from __future__ import annotations

import json
import os
import time
from typing import Any, Dict

ENV = "ENSEMBL_TRACE"


def enabled() -> bool:
    return bool(os.environ.get(ENV))


def _safe(value: Any) -> Any:
    try:
        json.dumps(value)
        return value
    except (TypeError, ValueError):
        return repr(value)[:200]


def record(tool: str, arguments: Dict[str, Any], result: Dict[str, Any],
           started: float, adapter=None) -> None:
    path = os.environ.get(ENV)
    if not path:
        return
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    build_id = None
    try:
        build_id = (adapter.meta() or {}).get("build_id") if adapter is not None else None
    except Exception:  # noqa: BLE001 - tracing must never break a query
        pass
    entry = {
        "ts": time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime()),
        "tool": tool,
        "arguments": _safe(arguments),
        "status": result.get("status"),
        "ok": result.get("ok"),
        "null_status": result.get("null_status"),
        "items": len(result.get("candidates") or result.get("items")
                     or result.get("records") or result.get("walks")
                     or result.get("datasets") or []),
        "elapsed_ms": round((time.perf_counter() - started) * 1000, 2),
        "build_id": build_id,
    }
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry, ensure_ascii=False, sort_keys=True) + "\n")
