#!/usr/bin/env python3
"""Drive the Ensembl slice exactly the way SciDataGym-Demo's M3 drives UniProt.

Each check mirrors a real call site in ``src/scidatagym/m3/engine.py``:

    engine.py:119  adapter.resolve({"text": ..., "type": ...}, constraints)
    engine.py:128  adapter.recall({"text": ..., "target_type": ...}, view, limit)
    engine.py:195  adapter.neighbors(node, [relation], direction, room, cursor)
    engine.py:315  adapter.get_records([node], fields, region)["items"][0]

Usage:
    python3 bridge/smoke.py
"""

from __future__ import annotations

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from bridge.scidatagym_m1 import EnsemblM1  # noqa: E402

GENE = "Ensembl:ENSG00000166913"
PROTEIN = "Ensembl:ENSP00000300161"

checks = []


def check(name, condition, detail=""):
    checks.append({"name": name, "pass": bool(condition), "detail": str(detail)[:160]})
    print(f"[{'PASS' if condition else 'FAIL'}] {name}" + (f"  {detail}" if not condition else ""))
    return bool(condition)


def main() -> int:
    m1 = EnsemblM1()

    # --- 1. describe() -------------------------------------------------
    described = m1.describe()
    check("describe: views are declared, not implied",
          described["views"] == ["identifier", "field"], described["views"])
    check("describe: manifest carries release + build",
          described["manifest"]["release"] == "116" and described["manifest"]["build_id"])
    check("describe: relations have endpoints for M3 to validate",
          described["relations"]["pdb_xref"] == ["Protein", "StructureEntry"],
          described["relations"].get("pdb_xref"))

    # --- 2. resolve() --------------------------------------------------
    resolved = m1.resolve({"text": "P31946"})
    check("resolve: UniProt accession is ambiguous across two proteins",
          resolved["status"] == "ambiguous" and len(resolved["items"]) >= 2,
          resolved["status"])
    check("resolve: items use the adapter's id/type/properties shape",
          {"id", "type", "label", "properties"} <= set(resolved["items"][0]))
    single = m1.resolve({"text": "YWHAB", "type": "Gene"}, {"taxon_id": 9606})
    check("resolve: a named gene resolves to one node",
          single["status"] == "resolved" and single["items"][0]["id"] == GENE,
          single["status"])
    missing = m1.resolve({"text": "ZZZ_NO_SUCH_GENE"})
    check("resolve: unknown symbol is not_found (not a crash)",
          missing["status"] == "not_found" and missing["items"] == [])

    # --- 3. recall() ---------------------------------------------------
    recalled = m1.recall({"text": "YWHAB", "target_type": "Gene"}, "identifier", 5)
    check("recall(identifier): returns ranked hits",
          recalled["items"] and recalled["items"][0]["rank"] == 1)
    text_hits = m1.recall({"text": "YWHAB", "target_type": "Transcript"}, "field", 5)
    check("recall(field): target_type filter is honoured",
          all(hit["node"]["type"] == "Transcript" for hit in text_hits["items"]))
    try:
        m1.recall({"text": "YWHAB"}, "semantic")
        check("recall(semantic): rejected explicitly", False)
    except ValueError:
        check("recall(semantic): rejected explicitly instead of returning junk", True)

    # --- 4. neighbors() with cursor pagination -------------------------
    first = m1.neighbors(GENE, ["has_transcript"], "out", 5, 0)
    second = m1.neighbors(GENE, ["has_transcript"], "out", 5, 5)
    check("neighbors: page size respected", len(first["items"]) == 5, len(first["items"]))
    check("neighbors: next_cursor allows the next page",
          first["next_cursor"] == 5 and second["items"][0]["node"]["id"] != first["items"][0]["node"]["id"])
    check("neighbors: edge carries relation + properties",
          first["items"][0]["edge"]["relation"] == "has_transcript")
    check("neighbors: source block is attached to every hit",
          all(item["source"] for item in first["items"]))
    pdb_page = m1.neighbors(PROTEIN, ["pdb_xref"], "out", 50, 0)
    check("neighbors: protein reaches PDB 4DNK",
          any(item["node"]["id"] == "PDB:4DNK" for item in pdb_page["items"]))
    ghost = m1.neighbors("Ensembl:ENSG00000000000")
    check("neighbors: unknown node returns node_not_found, not an exception",
          ghost["status"] == "node_not_found" and ghost["items"] == [])
    try:
        m1.neighbors(GENE, ["go_annotation"])
        check("neighbors: unsupported relation is rejected", False)
    except ValueError:
        check("neighbors: unsupported relation is rejected before querying", True)

    # --- 5. get_records() with a region --------------------------------
    record = m1.get_records([PROTEIN], ["sequence"], {"start": 230, "end": 246})["items"][0]
    check("get_records: region is translated to our protein_range",
          record["data"].get("sequence") == "WTSENQGDEGDAGEGEN", record["data"].get("sequence"))
    check("get_records: region echo matches the M3 contract",
          record["data"]["region"]["start"] == 230 and record["data"]["region"]["end"] == 246)
    full = m1.get_records([PROTEIN], ["sequence"])["items"][0]
    check("get_records: without region the full sequence is returned",
          len(full["data"]["sequence"]) == 246, len(full["data"].get("sequence") or ""))
    broken = m1.get_records(["Ensembl:ENSG00000000000"], ["sequence"])["items"][0]
    check("get_records: missing id is reported per item",
          broken["status"] == "not_found")

    # --- 6. raw_record / walk -----------------------------------------
    check("raw_record: says not_supported instead of inventing a record",
          m1.raw_record("P31946")["status"] == "not_supported")
    walked = m1.walk([GENE], ["has_transcript"], max_steps=3, budget=10, random_seed=7)
    replay = m1.walk([GENE], ["has_transcript"], max_steps=3, budget=10, random_seed=7)
    check("walk: same seed replays identically", walked == replay)

    failed = [item for item in checks if not item["pass"]]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} bridge checks passed")
    report = os.path.join(ROOT, "runs", "bridge_smoke.json")
    os.makedirs(os.path.dirname(report), exist_ok=True)
    with open(report, "w", encoding="utf-8") as handle:
        json.dump({"checks": checks, "passed": len(checks) - len(failed),
                   "failed": len(failed), "ok": not failed,
                   "manifest": m1.manifest}, handle, indent=2, ensure_ascii=False)
    print(f"-> {report}")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
