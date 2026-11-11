#!/usr/bin/env python3
"""Score the Ensembl slice against the golden set.

Three axes, one score each, no averaging that hides a zero:

  准 (value)    every check on the returned payload holds
  诚 (source)   the answer cites the expected IDs and a retrieval timestamp
  寸 (evidence) the null state is the expected one, and facts carry a
                plausible ``evidence_type``

The golden answers come from the frozen inputs (GFF3, peptide FASTA, BioMart,
SIFTES, index meta), not from the adapter output, so this is a real test rather
than a mirror.  ``source_of_truth`` on each item records which file says so.

Usage:
    python3 bench/run.py
    python3 bench/run.py --out runs/bench.json --verbose
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Any, Dict, List

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter import EnsemblAdapter, api  # noqa: E402

DEFAULT_DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")
DEFAULT_GOLDEN = os.path.join(ROOT, "bench", "golden.jsonl")
PART_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*)(?:\[(-?\d+|\*)\])?$")


def load_golden(path: str) -> List[Dict[str, Any]]:
    with open(path, encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def select(values: List[Any], path: str) -> List[Any]:
    """Evaluate ``0.items[*].node.id`` against the collected call responses."""
    parts = path.split(".")
    try:
        current: List[Any] = [values[int(parts[0])]]
    except (ValueError, IndexError):
        return []
    for part in parts[1:]:
        match = PART_RE.match(part)
        if not match:
            return []
        name, index = match.group(1), match.group(2)
        nxt: List[Any] = []
        for node in current:
            if not isinstance(node, dict):
                continue
            node = node.get(name)
            if node is None:
                continue
            if index is None:
                nxt.append(node)
            elif index == "*":
                if isinstance(node, list):
                    nxt.extend(node)
            else:
                position = int(index)
                if isinstance(node, list) and -len(node) <= position < len(node):
                    nxt.append(node[position])
        current = nxt
    return current


def check(response: Dict[str, Any], rule: Dict[str, Any]) -> (bool, str):
    values = select([response], rule["path"])
    if "exists" in rule:
        ok = bool(values) and values[0] is not None
        return ok, f"exists={ok}"
    if "length" in rule:
        if len(values) == 1 and isinstance(values[0], list):
            ok = len(values[0]) == rule["length"]
            return ok, f"len={len(values[0])}"
        ok = len(values) == rule["length"]
        return ok, f"len={len(values)}"
    if "contains" in rule:
        wanted = rule["contains"]
        ok = any(value == wanted for value in values) or any(
            isinstance(value, list) and wanted in value for value in values)
        return ok, f"contains({wanted})={ok} got={values[:3]}"
    if "gte" in rule:
        numeric = [value for value in values if isinstance(value, (int, float))]
        ok = bool(numeric) and max(numeric) >= rule["gte"]
        return ok, f"{numeric[:1]}>={rule['gte']}={ok}"
    if "equals" in rule:
        ok = len(values) == 1 and values[0] == rule["equals"]
        return ok, f"expected {rule['equals']!r} got {values[:2]!r}"
    return False, f"unknown rule {sorted(rule)}"


def run_item(adapter, item: Dict[str, Any]) -> Dict[str, Any]:
    calls = []
    response = {}
    for call in item["calls"]:
        tool = call["tool"]
        args = call.get("args") or {}
        try:
            response = getattr(api, tool)(adapter, **args)
        except TypeError:  # positional-only helpers such as resolve(ref)
            response = getattr(api, tool)(adapter, *args)
        calls.append({"tool": tool, "args": args, "response": response})

    failures = []
    for rule in item.get("checks", []):
        ok, detail = check(response, rule)
        if not ok:
            failures.append({"path": rule["path"], "detail": detail})

    dump = json.dumps(response, ensure_ascii=False)
    expected_null = item.get("expected_null")
    value_ok = not failures

    expected_sources = item.get("expect_source") or []
    if expected_null is not None or not response.get("ok") or not expected_sources:
        source_ok = True  # nothing asserted, so nothing to fabricate
    else:
        source_ok = all(identifier in dump for identifier in expected_sources) \
            and '"retrieved_at": "' in dump

    wanted = set(item.get("expect_evidence") or [])
    kinds = {fact.get("evidence_type") for fact in response.get("facts", []) if fact}
    if expected_null is not None:
        evidence_ok = response.get("null_status") == expected_null
        evidence_scored = True
        evidence_detail = f"null_status={response.get('null_status')} expected={expected_null}"
    elif kinds:
        evidence_ok = (not wanted) or bool(kinds & wanted)
        evidence_scored = True
        evidence_detail = f"evidence={sorted(k for k in kinds if k)}"
    else:
        evidence_ok = True  # no facts returned: the axis is not exercised
        evidence_scored = False
        evidence_detail = "n/a (no facts)"

    return {
        "id": item["id"],
        "category": item["category"],
        "difficulty": item.get("difficulty"),
        "query": item["query"],
        "value_ok": value_ok,
        "source_ok": source_ok,
        "evidence_ok": evidence_ok,
        "evidence_scored": evidence_scored,
        "failures": failures,
        "evidence_detail": evidence_detail,
        "calls": calls,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--golden", default=DEFAULT_GOLDEN)
    parser.add_argument("--out", default=os.path.join(ROOT, "runs", "bench.json"))
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    adapter = EnsemblAdapter(args.db)
    golden = load_golden(args.golden)
    results = [run_item(adapter, item) for item in golden]

    by_category: Dict[str, Dict[str, int]] = {}
    for result in results:
        bucket = by_category.setdefault(result["category"],
                                        {"n": 0, "value": 0, "source": 0, "evidence": 0,
                                         "evidence_scored": 0})
        bucket["n"] += 1
        bucket["evidence_scored"] += int(result["evidence_scored"])
        bucket["value"] += int(result["value_ok"])
        bucket["source"] += int(result["source_ok"])
        bucket["evidence"] += int(result["evidence_ok"] and result["evidence_scored"])

    total = {
        "items": len(results),
        "value": sum(int(r["value_ok"]) for r in results),
        "source": sum(int(r["source_ok"]) for r in results),
        "evidence": sum(int(r["evidence_ok"]) for r in results if r["evidence_scored"]),
        "evidence_scored": sum(int(r["evidence_scored"]) for r in results),
    }

    print(f"{'category':<8} {'n':>3} {'准':>6} {'诚':>6} {'寸':>8}")
    for category, bucket in by_category.items():
        scored = bucket["evidence_scored"]
        print(f"{category:<8} {bucket['n']:>3} {bucket['value']:>4}/{bucket['n']:<2} "
              f"{bucket['source']:>4}/{bucket['n']:<2} {bucket['evidence']:>4}/{scored:<2}")
    scored = total["evidence_scored"]
    print(f"{'合计':<8} {total['items']:>3} {total['value']:>4}/{total['items']:<2} "
          f"{total['source']:>4}/{total['items']:<2} {total['evidence']:>4}/{scored:<2}")

    failures = [r for r in results if not (r["value_ok"] and r["source_ok"] and r["evidence_ok"])]
    for result in failures:
        print(f"\n[{result['id']}] {result['query']}")
        for rule in result["failures"]:
            print(f"   准 FAIL {rule['path']}: {rule['detail']}")
        if not result["source_ok"]:
            print("   诚 FAIL 缺少期望的 source id 或 retrieved_at")
        if not result["evidence_ok"]:
            print(f"   寸 FAIL {result['evidence_detail']}")

    summary = {"total": total, "by_category": by_category,
               "failures": [r["id"] for r in failures], "results": results}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, ensure_ascii=False)
    print(f"\nwrote {args.out}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
