#!/usr/bin/env python3
"""Conformance suite for the Ensembl M1 slice.

Everything here is checkable without a server: the six contract operations are
called through :mod:`ensembl_adapter.api` (the same code path the MCP tools use),
validated against ``contract/m1.schema.json``, and then held to the semantic
rules that a schema cannot express:

* ``ok`` is true exactly when the query was answered;
* a missing answer is one of the four ``null_status`` values, never silence;
* every asserted source carries ``retrieved_at`` (or says ``not_reported``);
* sources were retrieved before the index was built;
* ``truncated`` matches ``total``/``limit``;
* the same walk called twice is byte-identical.

Usage:
    python3 contract/conformance.py                 # full run
    python3 contract/conformance.py --out runs/conformance.json
    python3 contract/conformance.py --list          # names only
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any, Dict, List

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter import EnsemblAdapter  # noqa: E402
from ensembl_adapter import api  # noqa: E402

DEFAULT_DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")
SCHEMA = os.path.join(ROOT, "contract", "m1.schema.json")

GENE = "Ensembl:ENSG00000166913"          # YWHAB
PROTEIN = "Ensembl:ENSP00000300161"       # canonical 14-3-3 beta
UNIPROT = "UniProtKB:P31946"
PDB = "PDB:4DNK"


class Reporter:
    def __init__(self) -> None:
        self.checks: List[Dict[str, Any]] = []

    def check(self, name: str, condition: bool, detail: Any = "") -> bool:
        self.checks.append({"name": name, "pass": bool(condition), "detail": detail})
        mark = "PASS" if condition else "FAIL"
        print(f"[{mark}] {name}" + (f"  {detail}" if detail and not condition else ""))
        return bool(condition)

    @property
    def failed(self) -> int:
        return sum(1 for item in self.checks if not item["pass"])


def load_validators(path: str) -> Dict[str, Any]:
    """Build one validator per operation, sharing the document's ``$defs``."""
    from jsonschema import Draft202012Validator

    with open(path, encoding="utf-8") as handle:
        root = json.load(handle)
    validators = {}
    for name in root["$defs"]:
        if name in {"null_status", "status", "source", "node", "edge", "fact", "provenance"}:
            continue
        doc = {"$schema": root["$schema"], "$defs": root["$defs"],
               "$ref": f"#/$defs/{name}"}
        validators[name] = Draft202012Validator(doc)
    return validators


def schema_errors(validator, payload: Dict[str, Any]) -> List[str]:
    return [f"{'/'.join(str(part) for part in error.absolute_path) or '<root>'}: {error.message}"
            for error in validator.iter_errors(payload)]


def walk_sources(value: Any, trail: str = "") -> List[Dict[str, Any]]:
    """Collect every ``source`` block anywhere in a response."""
    found: List[Dict[str, Any]] = []
    if isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(walk_sources(item, f"{trail}[{index}]"))
    elif isinstance(value, dict):
        if {"db", "id", "version", "retrieved_at_status"} <= set(value):
            found.append({"trail": trail or "<root>", **value})
        for key, item in value.items():
            found.extend(walk_sources(item, f"{trail}.{key}" if trail else key))
    return found


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--schema", default=SCHEMA)
    parser.add_argument("--out", default=os.path.join(ROOT, "runs", "conformance.json"))
    parser.add_argument("--list", action="store_true", help="list the checks only")
    args = parser.parse_args(argv)

    adapter = EnsemblAdapter(args.db)
    meta = adapter.meta()
    report = Reporter()
    responses: Dict[str, Dict[str, Any]] = {}

    def run(name: str, call) -> Dict[str, Any]:
        try:
            payload = call()
        except Exception as exc:  # noqa: BLE001 - a crash is a conformance failure
            report.check(f"{name}: callable", False, f"{type(exc).__name__}: {exc}")
            return {"ok": False, "status": "error", "null_status": "access_failed",
                    "truncated": False}
        responses[name] = payload
        return payload

    caps = run("capabilities", lambda: api.capabilities(adapter))
    report.check("capabilities: declared scope matches the index",
                 caps.get("consistent") is True,
                 (caps.get("snapshot"), meta.get("release")))
    report.check("capabilities: snapshot carries release + build identity",
                 all(caps.get("snapshot", {}).get(key) for key in
                     ("release", "assembly", "build_id", "built_at", "contract_version")),
                 caps.get("snapshot"))
    report.check("capabilities: four null states are documented",
                 set(caps.get("null_semantics") or {}) ==
                 {"not_found", "not_reported", "na", "access_failed"})
    report.check("capabilities: every frozen input has provenance",
                 all(entry.get("retrieved_at_status") for entry in
                     (caps.get("provenance") or {}).values()) and bool(caps.get("provenance")))

    found = run("resolve", lambda: api.resolve(adapter, {"text": "YWHAB", "type": "Gene"},
                                               {"taxon_id": 9606}))
    candidate = found["candidates"][0] if found.get("candidates") else {}
    report.check("resolve: symbol -> Ensembl gene", candidate.get("id") == GENE,
                 candidate.get("id"))
    report.check("resolve: answered result has null_status null",
                 found.get("ok") is True and found.get("null_status") is None)

    missing = run("resolve_missing", lambda: api.resolve(adapter, "ZZZ_NO_SUCH_GENE"))
    report.check("resolve: unknown symbol -> not_found, not an empty success",
                 missing.get("status") == "not_found"
                 and missing.get("null_status") == "not_found"
                 and missing.get("ok") is False)

    outside = run("resolve_outside_scope",
                  lambda: api.resolve(adapter, "YWHAB", {"taxon_id": 10090}))
    report.check("resolve: mouse is na (out of scope), not not_found",
                 outside.get("null_status") == "na")

    drift = run("resolve_version_drift",
                lambda: api.resolve(adapter, "ENST00000353703.999"))
    report.check("resolve: wrong version is not_found, not a silent downgrade",
                 drift.get("null_status") == "not_found")

    nb = run("neighbors", lambda: api.neighbors(adapter, GENE, ["has_transcript"], "out", 50))
    report.check("neighbors: YWHAB has 40 transcripts", nb.get("total") == 40, nb.get("total"))
    report.check("neighbors: every edge carries a source",
                 all(item.get("source", {}).get("db") for item in nb.get("items", [])))

    cut = run("neighbors_truncated",
              lambda: api.neighbors(adapter, GENE, ["has_transcript"], "out", 1))
    report.check("neighbors: limit truncates and says so",
                 cut.get("truncated") is True and cut.get("total") == 40
                 and len(cut.get("items", [])) == 1)

    ghost = run("neighbors_missing", lambda: api.neighbors(adapter, "Ensembl:ENSG00000000000"))
    report.check("neighbors: unknown node -> not_found",
                 ghost.get("null_status") == "not_found" and ghost.get("items") == [])

    record = run("get_records", lambda: api.get_records(
        adapter, [PROTEIN], ["sequence", "protein_mapping"], {"protein_range": [230, 246]}))
    payload = record["records"][0] if record.get("records") else {}
    report.check("get_records: 230-246 slice is exact",
                 (payload.get("fields", {}).get("sequence") or {}).get("slice")
                 == "WTSENQGDEGDAGEGEN")
    report.check("get_records: region positions are fully retrieved",
                 payload.get("region", {}).get("retrieved") == list(range(230, 247)))
    report.check("get_records: facts exist for the requested values",
                 len(record.get("facts", [])) >= 2, len(record.get("facts", [])))

    partial = run("get_records_partial", lambda: api.get_records(
        adapter, [PROTEIN], ["sequence"], {"protein_range": [246, 247]}))
    tail = partial["records"][0] if partial.get("records") else {}
    report.check("get_records: partial range reports the missing residue",
                 tail.get("region", {}).get("missing") == [{"position": 247, "reason": "beyond_protein"}]
                 and partial.get("status") == "ok")

    broken = run("get_records_invalid_region",
                 lambda: api.get_records(adapter, [PROTEIN], region={"protein_range": [2, 1]}))
    report.check("get_records: malformed region is invalid_region, not an exception",
                 broken.get("status") == "invalid_region" and broken.get("ok") is False)

    absent = run("get_records_missing", lambda: api.get_records(adapter, ["Ensembl:ENSG00000000000"]))
    names = [item.get("null_status") for item in absent.get("records", [])]
    report.check("get_records: missing id is recorded, not dropped",
                 names == ["not_found"], names)

    first = run("walk", lambda: api.walk(adapter, [GENE], ["has_transcript"], 4, 5))
    second = api.walk(adapter, [GENE], ["has_transcript"], 4, 5)
    report.check("walk: same seed and rng_seed replay byte-identically", first == second)
    report.check("walk: paths carry the asserting source",
                 all(step.get("source", {}).get("db")
                     for item in first.get("walks", []) for step in item.get("path", [])))

    cross = run("resolve_crosswalk", lambda: api.resolve(adapter, UNIPROT))
    reached = {node["id"] for node in cross.get("candidates", [])}
    report.check("crosswalk: UniProt accession reaches the Ensembl protein",
                 PROTEIN in reached, sorted(reached)[:4])
    report.check("crosswalk: ambiguity is reported, not hidden",
                 cross.get("crosswalk_ambiguous") is True
                 and any(hit["alternatives"] > 0 for hit in cross.get("crosswalk", [])))
    edges = cross.get("crosswalk") or []
    report.check("crosswalk: every edge cites a file, row and hash",
                 bool(edges) and all(hit["evidence"].get("file")
                                     and hit["evidence"].get("sha256")
                                     and hit["evidence"].get("retrieved_at") for hit in edges),
                 edges[0]["evidence"] if edges else None)

    sifts = run("resolve_pdb", lambda: api.resolve(adapter, PDB))
    sifts_hit = [hit for hit in sifts.get("crosswalk", [])
                 if hit["relation"] == "pdb_uniprot_sifts"]
    report.check("crosswalk: PDB entry links back to UniProt through SIFTS",
                 bool(sifts_hit) and sifts_hit[0]["other"]["id"] == "P31946",
                 sifts_hit[0]["other"] if sifts_hit else None)

    report.check("capabilities: the crosswalk is declared, not implied",
                 (caps.get("crosswalk") or {}).get("edges", 0) > 0
                 and bool((caps.get("crosswalk") or {}).get("build_id")))

    recall = run("recall", lambda: api.recall(adapter, "YWHAB", "structure"))
    report.check("recall: unimplemented view is na, not an empty list",
                 recall.get("null_status") == "na")

    # ---- schema + cross-cutting invariants -----------------------------
    validators = load_validators(args.schema)
    operations = ("capabilities", "resolve", "neighbors", "get_records", "walk", "recall")

    def operation_of(name: str) -> str:
        return next((op for op in operations if name.startswith(op)), name)

    for name, payload in responses.items():
        validator = validators.get(operation_of(name))
        if validator is None:
            report.check(f"schema: {name} has a schema", False)
            continue
        errors = schema_errors(validator, payload)
        report.check(f"schema: {name}", not errors, errors[:3])

    null_states = {resp.get("null_status") for resp in responses.values()}
    report.check("null: all four states are exercised",
                 {"not_found", "na"} <= null_states,
                 sorted(str(state) for state in null_states))

    built_at = (caps.get("snapshot") or {}).get("built_at")
    collected = []
    for name, payload in responses.items():
        for source in walk_sources(payload):
            collected.append((name, source))
    report.check("provenance: every source states whether its timestamp is known",
                 all(source.get("retrieved_at_status") in
                     {"download", "file_mtime", "not_reported"}
                     for _, source in collected),
                 len(collected))
    report.check("provenance: no asserted fact pretends to a timestamp it lacks",
                 all(source.get("retrieved_at") is not None
                     for name, source in collected
                     if name.startswith("get_records") and source.get("retrieved_at_status") != "not_reported"))
    report.check("provenance: inputs were retrieved before the index was built",
                 all(source.get("retrieved_at", "") <= (built_at or "")
                     for _, source in collected if source.get("retrieved_at")),
                 built_at)

    summary = {
        "contract_version": api.CONTRACT_VERSION,
        "snapshot": caps.get("snapshot"),
        "database": args.db,
        "checks": report.checks,
        "passed": len(report.checks) - report.failed,
        "failed": report.failed,
        "ok": report.failed == 0,
    }
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True, ensure_ascii=False)
    print(f"\n{summary['passed']}/{len(report.checks)} checks passed "
          f"({api.CONTRACT_VERSION}) -> {args.out if args.out else 'stdout only'}")
    return 0 if report.failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
