#!/usr/bin/env python3
"""Local console for the Ensembl M1 slice.

A thin HTTP wrapper around :mod:`ensembl_adapter.api` — the same functions the
MCP server exposes — plus a static single-page UI.  Read-only, bound to
127.0.0.1, no third-party packages.

    python3 viz/server.py --port 8765
    open http://127.0.0.1:8765
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import traceback
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import parse_qs, urlparse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter import EnsemblAdapter, api  # noqa: E402

WEB = os.path.join(ROOT, "viz", "web")
DEFAULT_DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")
CONTENT_TYPES = {".html": "text/html; charset=utf-8",
                 ".js": "text/javascript; charset=utf-8",
                 ".css": "text/css; charset=utf-8",
                 ".svg": "image/svg+xml"}


def parse_region(value: str | None):
    """``protein_range:230-246`` -> ``{"protein_range": [230, 246]}``"""
    if not value or ":" not in value:
        return None
    kind, span = value.split(":", 1)
    if "-" not in span:
        return None
    start, end = span.split("-", 1)
    try:
        return {kind: [int(start), int(end)]}
    except ValueError:
        return None


class Handler(BaseHTTPRequestHandler):
    adapter: EnsemblAdapter = None  # set in main()

    def log_message(self, fmt, *args):  # quieter default logging
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    def _send(self, body: bytes, status: int = 200, content_type: str = "application/json"):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _json(self, payload, status: int = 200):
        self._send(json.dumps(payload, ensure_ascii=False).encode("utf-8"), status)

    def do_GET(self):  # noqa: N802 - stdlib naming
        parsed = urlparse(self.path)
        query = {key: values[0] for key, values in parse_qs(parsed.query).items()}
        try:
            if parsed.path.startswith("/api/"):
                return self._json(self._api(parsed.path[5:], query))
            return self._static(parsed.path)
        except Exception as exc:  # noqa: BLE001 - surface the failure to the UI
            traceback.print_exc()
            return self._json({"error": f"{type(exc).__name__}: {exc}"}, 500)

    # ------------------------------------------------------------------
    def _api(self, name: str, query: dict):
        adapter = self.adapter
        if name == "capabilities":
            return api.capabilities(adapter)
        if name == "search":
            ref = {"text": query.get("q", "")}
            if query.get("type"):
                ref["type"] = query["type"]
            constraints = {"taxon_id": int(query["taxon"])} if query.get("taxon") else None
            return api.resolve(adapter, ref, constraints)
        if name == "graph":
            node_id = query.get("id", "")
            record = api.get_records(adapter, [node_id])
            record_row = (record.get("records") or [{}])[0]
            node = record_row.get("node")
            if not node:
                return {"ok": False, "status": "not_found", "null_status": "not_found",
                        "node": None, "items": []}
            result = api.neighbors(adapter, node["id"],
                                   relations=query.get("relations") or None,
                                   direction=query.get("direction", "both"),
                                   limit=int(query.get("limit", 24)))
            hits = adapter.crosswalk_lookup(node["id"])
            external = node["id"].split(":", 1)[0] in {"UniProtKB", "PDB", "HGNC", "CCDS", "MIM"}
            result.update({
                "node": node,
                "facts": record.get("facts", []),
                "crosswalk": hits[:8],
                # many-to-many only matters when the anchor is an external key;
                # one protein mapping to 17 structures is normal, not ambiguity
                "crosswalk_ambiguous": external and any(hit["alternatives"] > 0 for hit in hits),
            })
            return result
        if name == "records":
            ids = [item for item in (query.get("ids") or "").split(",") if item]
            fields = [item for item in (query.get("fields") or "").split(",") if item] or None
            return api.get_records(adapter, ids, fields, parse_region(query.get("region")))
        if name == "crosswalk":
            text = (query.get("q") or "").strip()
            hits = adapter.crosswalk_lookup(text) if text else []
            resolved = api.resolve(adapter, text) if text else None
            return {"query": text, "hits": hits,
                    "crosswalk_ambiguous": bool(resolved and resolved.get("crosswalk_ambiguous")),
                    "candidates": (resolved or {}).get("candidates", [])}
        if name == "health":
            return _health(adapter)
        if name == "targets":
            return _targets(adapter)
        if name == "isoforms":
            return _isoforms(query)
        if name == "bridge":
            return _bridge()
        if name == "bench":
            return {"bench": _read_json(os.path.join(ROOT, "runs", "bench.json")),
                    "conformance": _read_json(os.path.join(ROOT, "runs", "conformance.json")),
                    "verify": _read_json(os.path.join(ROOT, "runs", "verify.json"))}
        return {"error": f"unknown endpoint {name}"}

    def _static(self, path: str):
        name = "index.html" if path in ("/", "") else os.path.basename(path)
        full = os.path.join(WEB, name)
        if not os.path.isfile(full):
            return self._json({"error": "not found", "path": path}, 404)
        extension = os.path.splitext(full)[1]
        with open(full, "rb") as handle:
            self._send(handle.read(), 200, CONTENT_TYPES.get(extension, "application/octet-stream"))


def _targets(adapter) -> dict:
    """A small menu of concrete proteins: canonical per gene, and structure-rich ones."""
    conn = adapter.conn
    genes = []
    rows = conn.execute(
        "SELECT g.stable_id AS gid, g.label AS symbol, COUNT(*) AS n "
        "FROM node g JOIN edge e ON e.src=g.node_id AND e.relation='has_transcript' "
        "WHERE g.kind='gene' AND g.biotype='protein_coding' GROUP BY g.stable_id "
        "ORDER BY n DESC LIMIT 12").fetchall()
    for row in rows:
        transcript = adapter.canonical_transcript(row["gid"])
        if not transcript:
            continue
        protein = adapter.transcript_translation(transcript["stable_id"])
        if not protein:
            continue
        genes.append({
            "value": protein["id"],
            "label": f"{row['symbol']} · canonical {protein['stable_id']} ({protein['length']} aa)",
        })
    structures = []
    for row in conn.execute(
            "SELECT s.stable_id AS pid, COUNT(*) AS n FROM edge e "
            "JOIN node s ON s.node_id=e.src "
            "WHERE e.relation='pdb_xref' AND s.kind='translation' "
            "GROUP BY s.stable_id ORDER BY n DESC LIMIT 10"):
        structures.append({
            "value": f"Ensembl:{row['pid']}",
            "label": f"{row['pid']} · {row['n']} 个 PDB 结构",
        })
    return {"genes": genes, "structures": structures}


def _isoforms(query: dict) -> dict:
    """Every isoform of one gene, and whether it covers the requested window."""
    gene = str(query.get("gene") or "").strip()
    start = int(query.get("start") or 1)
    end = int(query.get("end") or start)
    adapter = Handler.adapter
    resolved = api.resolve(adapter, {"text": gene, "type": "Gene"}, {"taxon_id": 9606})
    if not resolved.get("candidates"):
        return {"gene": gene, "status": resolved.get("status"), "rows": []}
    node = resolved["candidates"][0]
    transcripts = api.neighbors(adapter, node["id"], ["has_transcript"], "out", 1000)
    rows = []
    for item in transcripts.get("items", [])[:60]:
        transcript = item["node"]["stable_id"]
        try:
            coverage = adapter.transcript_coverage(transcript, start, end)
        except Exception:  # noqa: BLE001 - one bad isoform must not kill the table
            continue
        positions = coverage.get("positions") or []
        rows.append({
            "transcript": transcript,
            "version": coverage.get("transcript_version"),
            "protein": coverage.get("protein_id"),
            "protein_length": coverage.get("protein_length"),
            "biotype": coverage.get("biotype"),
            "canonical": bool(coverage.get("is_canonical")),
            "status": coverage.get("status"),
            "retrieved": sum(1 for pos in positions if pos.get("status") == "ok"),
            "expected": max(0, min(end, coverage.get("protein_length") or 0) - start + 1),
        })
    rows.sort(key=lambda row: (row["status"] != "covered", not row["canonical"],
                               -(row["protein_length"] or 0)))
    covered = sum(1 for row in rows if row["status"] == "covered")
    return {"gene": gene, "gene_id": node["id"], "region": [start, end],
            "total": len(rows), "covered": covered, "rows": rows}


def _bridge() -> dict:
    """Run the SciDataGym-Demo M1 surface over this slice, live."""
    from bridge.scidatagym_m1 import EnsemblM1

    m1 = EnsemblM1(Handler.adapter.db_path)
    described = m1.describe()
    resolved = m1.resolve({"text": "P31946"})
    page1 = m1.neighbors("Ensembl:ENSG00000166913", ["has_transcript"], "out", 3, 0)
    page2 = m1.neighbors("Ensembl:ENSG00000166913", ["has_transcript"], "out", 3, 3)
    record = m1.get_records(["Ensembl:ENSP00000300161"], ["sequence"],
                            {"start": 230, "end": 246})["items"][0]
    return {
        "adapter_id": described["adapter_id"],
        "manifest": described["manifest"],
        "views": described["views"],
        "types": described["types"],
        "relations": dict(list(described["relations"].items())[:6]),
        "unsupported": described["unsupported"],
        "smoke": _read_json(os.path.join(ROOT, "runs", "bridge_smoke.json")),
        "demos": {
            "resolve": {
                "status": resolved["status"],
                "items": [{"id": item["id"], "type": item["type"],
                           "label": item["label"]} for item in resolved["items"][:4]],
            },
            "neighbors": {
                "page1": [item["node"]["id"] for item in page1["items"]],
                "page2": [item["node"]["id"] for item in page2["items"]],
                "next_cursor": page1["next_cursor"],
                "relation": page1["items"][0]["edge"]["relation"] if page1["items"] else None,
                "source": page1["items"][0]["source"] if page1["items"] else None,
            },
            "get_records": {
                "id": record["id"],
                "sequence": record["data"].get("sequence"),
                "region": record["data"].get("region"),
            },
            "raw_record": m1.raw_record("P31946"),
        },
    }


def _health(adapter) -> dict:
    """Everything needed to judge the index itself, not one answer from it."""
    conn = adapter.conn

    def scalar(sql: str) -> int:
        row = conn.execute(sql).fetchone()
        return int(row[0]) if row and row[0] is not None else 0

    def covered(relation: str, kind: str) -> int:
        return scalar(
            "SELECT COUNT(DISTINCT e.src) FROM edge e JOIN node s ON s.node_id=e.src "
            f"WHERE e.relation='{relation}' AND s.kind='{kind}'")

    counts = {
        "genes": scalar("SELECT COUNT(*) FROM node WHERE kind='gene'"),
        "transcripts": scalar("SELECT COUNT(*) FROM node WHERE kind='transcript'"),
        "proteins": scalar("SELECT COUNT(*) FROM node WHERE kind='translation'"),
        "exons": scalar("SELECT COUNT(*) FROM node WHERE kind='exon'"),
        "xref_nodes": scalar("SELECT COUNT(*) FROM node WHERE kind='xref'"),
        "nodes": scalar("SELECT COUNT(*) FROM node"),
        "edges": scalar("SELECT COUNT(*) FROM edge"),
    }
    coverage = {
        "proteins_uniprot": covered("uniprot_xref", "translation"),
        "proteins_pdb": covered("pdb_xref", "translation"),
        "genes_hgnc": covered("hgnc_xref", "gene"),
        "genes_entrez": covered("entrez_xref", "gene"),
        "genes_mim": covered("mim_xref", "gene"),
        "transcripts_ccds": covered("ccds_xref", "transcript"),
    }
    edges_by_relation = dict(conn.execute(
        "SELECT relation, COUNT(*) FROM edge GROUP BY relation ORDER BY 2 DESC"))
    stats = _read_json(os.path.join(ROOT, "data", "index", "build_stats.json")) or {}
    return {
        "capabilities": api.capabilities(adapter),
        "counts": counts,
        "coverage": coverage,
        "edges_by_relation": edges_by_relation,
        "validation": (stats.get("counts") or {}).get("validation", {}),
        "bioptypes": (stats.get("counts") or {}).get("gene_biotypes", {}),
        "fingerprint": stats.get("fingerprint"),
        "elapsed_seconds": stats.get("elapsed_seconds"),
        "verify": _read_json(os.path.join(ROOT, "runs", "verify.json")),
        "conformance": _read_json(os.path.join(ROOT, "runs", "conformance.json")),
    }


def _read_json(path: str):
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args(argv)

    Handler.adapter = EnsemblAdapter(args.db)
    # single-threaded on purpose: one read-only SQLite connection, no locking
    server = HTTPServer((args.host, args.port), Handler)
    print(f"SciDataGym · Ensembl console -> http://{args.host}:{args.port}", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
