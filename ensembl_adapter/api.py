"""M1 contract layer for the Ensembl slice.

This module is the single place where the adapter's internal ``status`` values
become the public contract:

* every response carries ``ok`` / ``status`` / ``null_status`` / ``truncated``;
* ``null_status`` is one of ``not_found`` / ``not_reported`` / ``na`` /
  ``access_failed`` and is never faked;
* every fact carries a self-contained ``source`` with ``version``,
  ``retrieved_at`` and ``retrieved_at_status``.

Both the MCP server and the conformance runner call these functions, so the
schema in ``contract/m1.schema.json`` describes exactly what ships.

Layering: ``adapter.py`` reads the graph; this module shapes answers;
``mcp_server.py`` exposes them over MCP.
"""

from __future__ import annotations

import functools
import json
import os
import re
import sqlite3
import time
from typing import Any, Dict, List, Optional, Sequence

from . import trace
from .adapter import TYPE_PRIORITY

CONTRACT_VERSION = "m1.1"

#: internal status -> public null state
NULL_BY_STATUS = {
    "not_found": "not_found",
    "unsupported": "na",
    "error": "access_failed",
}

ANSWERED = ("ok", "partial")

#: which frozen file asserted a cross-reference field
FIELD_TO_SOURCE = {
    "BioMart pdb": "biomart/pdb.tsv",
    "BioMart ccds": "biomart/ccds_mim.tsv",
    "BioMart mim_gene_accession": "biomart/ccds_mim.tsv",
    "BioMart uniprotswissprot": "biomart/uniprot.tsv",
    "BioMart hgnc_id": "biomart/gene_xref.tsv",
    "BioMart entrezgene_id": "biomart/gene_xref.tsv",
}

VERSIONED_STABLE_RE = re.compile(r"(?:Ensembl:)?ENS[A-Z]+\d+\.(\d+)")

MAX_XREF_FACTS = 100
MAX_CROSSWALK = 20
MAX_CANDIDATES = 50


def _traced(tool: str):
    """Record every contract call when tracing is enabled."""
    def decorator(function):
        @functools.wraps(function)
        def wrapper(adapter, *args, **kwargs):
            started = time.perf_counter()
            result = function(adapter, *args, **kwargs)
            if trace.enabled():
                trace.record(tool, {"args": list(args), "kwargs": kwargs}, result,
                             started, adapter)
            return result
        return wrapper
    return decorator


# ----------------------------------------------------------------------
# Envelope
# ----------------------------------------------------------------------
def envelope(result: Dict[str, Any]) -> Dict[str, Any]:
    """Add the contract envelope without hiding the internal status."""
    status = str(result.get("status") or "error")
    out = dict(result)
    out["status"] = status
    out["ok"] = status in ANSWERED
    out["null_status"] = NULL_BY_STATUS.get(status)
    if out["null_status"] == "na":
        out.setdefault("scope_status", "unsupported")
    out.setdefault("truncated", False)
    if status == "ok" and out["null_status"] is not None:  # defensive
        out["null_status"] = None
    return out


def _source_key(field: Optional[str], meta: Dict[str, str]) -> Optional[str]:
    """Map a ``source.field`` label back to the frozen file that asserted it."""
    if not field:
        return None
    if field.startswith("GFF3"):
        return meta.get("source_gff3")
    if field == "FASTA pep":
        return meta.get("source_pep")
    if field in FIELD_TO_SOURCE:
        return FIELD_TO_SOURCE[field]
    if field.startswith("derived:"):
        return meta.get("source_gff3")
    return None


def _normalize_source(value: Any, meta: Dict[str, str], prov: "Provenance") -> Any:
    """Rename the graph's native ``database/record_id/release`` to contract keys.

    Every source block ends up with the same five required fields, including a
    provenance lookup, so a node's origin is as checkable as a fact's.
    """
    if isinstance(value, list):
        return [_normalize_source(item, meta, prov) for item in value]
    if not isinstance(value, dict):
        return value
    out = {key: _normalize_source(item, meta, prov) for key, item in value.items()}
    if "database" in out and "record_id" in out and "db" not in out:
        release = meta.get("release")
        entity_version = out.get("entity_version")
        version = release
        if out.get("field") == "FASTA pep":
            # the FASTA header's "version" is the release the peptides came from
            entity_version = out.get("release")
        out.update({
            "db": out.pop("database"),
            "id": out.pop("record_id"),
            "version": version,
        })
        out.pop("release", None)
        if entity_version is not None:
            out["entity_version"] = entity_version
        stamp = prov.get(_source_key(out.get("field"), meta) or "")
        out.setdefault("retrieved_at", stamp["retrieved_at"])
        out.setdefault("retrieved_at_status", stamp["retrieved_at_status"])
        out.setdefault("assembly", meta.get("assembly"))
        out.setdefault("build_id", meta.get("build_id"))
    return out


def _normalize_payload(value: Any, meta: Dict[str, str], prov: "Provenance") -> Any:
    """Normalize every nested ``source`` block in one response."""
    if isinstance(value, list):
        return [_normalize_payload(item, meta, prov) for item in value]
    if not isinstance(value, dict):
        return value
    out = {key: _normalize_payload(item, meta, prov) for key, item in value.items()}
    if "database" in out and "record_id" in out:
        out = _normalize_source(out, meta, prov)
    return out


# ----------------------------------------------------------------------
# Provenance
# ----------------------------------------------------------------------
class Provenance:
    """Lazy view of ``meta.source_retrieved_at``."""

    def __init__(self, meta: Dict[str, str]):
        self.meta = meta
        try:
            self.sources = json.loads(meta.get("source_retrieved_at") or "{}")
        except json.JSONDecodeError:
            self.sources = {}

    def get(self, key: str) -> Dict[str, Optional[str]]:
        entry = self.sources.get(key) or {}
        return {
            "retrieved_at": entry.get("retrieved_at"),
            "retrieved_at_status": entry.get("retrieved_at_status") or "not_reported",
            "path": entry.get("path"),
            "sha256": entry.get("sha256"),
        }

    def snapshot(self) -> Dict[str, Any]:
        return {
            "build_id": self.meta.get("build_id"),
            "built_at": self.meta.get("built_at"),
            "release": self.meta.get("release"),
            "assembly": self.meta.get("assembly"),
            "chromosome": self.meta.get("chromosome"),
            "region": self.meta.get("chromosome"),
            "schema_version": self.meta.get("schema_version"),
            "contract_version": CONTRACT_VERSION,
        }


def _fact_source(record: Dict[str, Any], prov: Provenance, key: str, meta: Dict[str, str],
                 *, field: Optional[str] = None, locator: Optional[str] = None) -> Dict[str, Any]:
    node = record["node"]
    native = record.get("source") or {}
    stamp = prov.get(key)
    return {
        "db": "Ensembl",
        "id": node.get("stable_id") or native.get("record_id"),
        "version": meta.get("release"),
        "entity_version": str(node.get("version")) if node.get("version") is not None else None,
        "field": field or native.get("field"),
        "locator": locator or native.get("locator"),
        "retrieved_at": stamp["retrieved_at"],
        "retrieved_at_status": stamp["retrieved_at_status"],
        "assembly": meta.get("assembly"),
        "build_id": meta.get("build_id"),
    }


def _fact(record: Dict[str, Any], attribute: str, value: Any, *, value_type: str,
          evidence_type: str, source: Dict[str, Any], scale: str,
          region: Optional[Dict[str, Any]] = None,
          derived_from: Optional[List[Dict[str, Any]]] = None,
          unit: Optional[str] = None) -> Dict[str, Any]:
    namespace, _, identifier = str(record["id"]).partition(":")
    value_block: Dict[str, Any] = {"type": value_type, "data": value}
    if unit:
        value_block["unit"] = unit
    fact = {
        "entity": {"id": identifier, "namespace": namespace or "Ensembl"},
        "attribute": attribute,
        "value": value_block,
        "evidence_type": evidence_type,
        "source": source,
        "scale": scale,
    }
    if region:
        fact["region"] = region
    if derived_from:
        fact["derived_from"] = derived_from
    return fact


# ----------------------------------------------------------------------
# Operations
# ----------------------------------------------------------------------
DATASET_FORMATS = {
    ".gff3.gz": "gff3+gzip",
    ".fa.gz": "fasta+gzip",
    ".fat.gz": "fasta+gzip",
    ".tsv": "tsv",
    ".json": "json",
    ".cif": "mmcif",
    ".sqlite": "sqlite",
}


def _dataset_format(path: str) -> str:
    lowered = path.lower()
    for suffix, label in DATASET_FORMATS.items():
        if lowered.endswith(suffix):
            return label
    return os.path.splitext(lowered)[1].lstrip(".") or "unknown"


def _human_bytes(size: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size} B"


@_traced("datasets")
def datasets(adapter) -> Dict[str, Any]:
    """DatasetRecord for every frozen input and for the derived index itself."""
    meta = adapter.meta()
    prov = Provenance(meta)
    root = os.path.dirname(os.path.abspath(adapter.db_path))
    for _ in range(4):                     # data/index/x.sqlite -> repo root
        if os.path.exists(os.path.join(root, "capabilities.json")):
            break
        root = os.path.dirname(root)
    records: List[Dict[str, Any]] = []

    for key in sorted(prov.sources):
        entry = prov.sources[key]
        path = os.path.join(root, entry.get("path") or "")
        size = os.path.getsize(path) if os.path.exists(path) else None
        records.append({
            "dataset_id": f"source:{key}",
            "name": os.path.basename(entry.get("path") or key),
            "description": "frozen input for the Ensembl graph build",
            "size": size,
            "size_human": _human_bytes(size) if size else None,
            "format": _dataset_format(key),
            "access": {"kind": "download", "fetch": "python3 fetch_data.py",
                       "read_only": True},
            "provenance": {
                "db": "Ensembl", "release": meta.get("release"),
                "path": entry.get("path"), "sha256": entry.get("sha256"),
                "retrieved_at": entry.get("retrieved_at"),
                "retrieved_at_status": entry.get("retrieved_at_status"),
            },
        })

    index_path = os.path.abspath(adapter.db_path)
    index_size = os.path.getsize(index_path) if os.path.exists(index_path) else None
    records.append({
        "dataset_id": f"index:{os.path.basename(index_path)}",
        "name": os.path.basename(index_path),
        "description": "SQLite graph index served by this adapter",
        "size": index_size,
        "size_human": _human_bytes(index_size) if index_size else None,
        "format": "sqlite",
        "access": {"kind": "build", "read_only": True,
                   "fetch": "./run_all.sh", "query": "resolve/neighbors/get_records"},
        "provenance": {
            "db": "Ensembl", "release": meta.get("release"),
            "assembly": meta.get("assembly"), "region": meta.get("chromosome"),
            "build_id": meta.get("build_id"), "built_at": meta.get("built_at"),
            "retrieved_at": meta.get("built_at"),
            "retrieved_at_status": "download" if meta.get("built_at") else "not_reported",
        },
    })

    return {"ok": True, "status": "ok", "null_status": None, "truncated": False,
            "datasets": records, "total": len(records)}


FACT_ATTRIBUTES = {
    "protein_sequence",
    "protein_position_mapping",
    "transcript_sequence",
    "exon_structure",
    "translation",
    "gene_location",
    "cross_reference",
}


@_traced("capabilities")
def capabilities(adapter) -> Dict[str, Any]:
    meta = adapter.meta()
    prov = Provenance(meta)
    declaration = adapter.capabilities()
    declared_snapshot = {
        "release": declaration.get("release"),
        "assembly": declaration.get("assembly"),
        "region": (declaration.get("scope") or {}).get("region"),
    }
    index_snapshot = {
        "release": meta.get("release"),
        "assembly": meta.get("assembly"),
        "region": meta.get("chromosome"),
    }
    return {
        "ok": True,
        "status": "ok",
        "null_status": None,
        "truncated": False,
        "contract_version": CONTRACT_VERSION,
        "database": declaration.get("database", "Ensembl"),
        "release": meta.get("release"),
        "assembly": meta.get("assembly"),
        "schema_version": meta.get("schema_version"),
        "scope": declaration.get("scope"),
        "node_types": declaration.get("node_types", []),
        "relations": declaration.get("relations", []),
        "fields": declaration.get("fields", []),
        "region_types": declaration.get("region_types", []),
        "limitations": declaration.get("limitations", []),
        "namespaces": declaration.get("namespaces") or [
            "Ensembl", "UniProtKB", "PDB", "HGNC", "EntrezGene", "MIM", "CCDS"],
        "key_alignment": declaration.get("key_alignment"),
        "scales": ["gene", "transcript", "exon", "protein"],
        "search_modes": ["identifier", "text"],
        "fact_attributes": sorted(FACT_ATTRIBUTES),
        "dataset_count": len(prov.sources) + 1,
        "trace": {"enabled": trace.enabled(), "env": "ENSEMBL_TRACE"},
        "evidence_types": ["computational", "inferred"],
        "null_semantics": {
            "not_found": "the query ran and nothing in this slice matched",
            "not_reported": "the entity exists but the source asserts no such value",
            "na": "outside the declared scope of this adapter",
            "access_failed": "the underlying read failed; nothing is claimed",
        },
        "crosswalk": {**adapter.crosswalk.stats(), "build_id": adapter.crosswalk.build_id},
        "snapshot": prov.snapshot(),
        "provenance": {key: prov.get(key) for key in sorted(prov.sources)},
        "consistent": declared_snapshot == index_snapshot,
    }


@_traced("resolve")
def resolve(adapter, ref: Any, constraints: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    result = adapter.resolve(ref, constraints)
    raw = ref if isinstance(ref, str) else (ref or {}).get("text") or (ref or {}).get("id") or ""
    version = VERSIONED_STABLE_RE.fullmatch(str(raw).strip())
    if version:
        result["candidates"] = [node for node in result["candidates"]
                                if str(node["version"]) == version.group(1)]
        if result["status"] == "ok" and not result["candidates"]:
            result["status"] = "not_found"
    result = _expand_crosswalk(adapter, result, str(raw).strip())
    result = _normalize_payload(result, adapter.meta(), Provenance(adapter.meta()))
    return envelope(result)


def _expand_crosswalk(adapter, result: Dict[str, Any], text: str) -> Dict[str, Any]:
    """Add authoritative answers from the key-alignment table.

    Resolving ``P31946`` used to stop at the external reference node; with a
    crosswalk it also answers the Ensembl protein the accession maps to, and
    reports how many other proteins share that accession instead of silently
    picking one.
    """
    if result.get("status") not in ANSWERED + ("not_found",) or not text:
        return result
    hits = adapter.crosswalk_lookup(text)
    if not hits:
        return result
    seen = {node["id"] for node in result["candidates"]}
    for hit in hits:
        other = hit["other"]
        if other["namespace"] != "Ensembl":
            continue
        if other["id"] in seen:
            continue
        node = adapter.get(other["id"])
        if not node:
            continue
        seen.add(other["id"])
        node = dict(node)
        node["matched_by"] = f"crosswalk:{hit['relation']}"
        node["crosswalk_alternatives"] = hit["alternatives"]
        result["candidates"].append(node)
    result["candidates"].sort(key=lambda node: (TYPE_PRIORITY.get(node["type"], 99), node["id"]))
    if result["candidates"]:
        result["status"] = "ok"
    if len(result["candidates"]) > MAX_CANDIDATES:
        result["candidates"] = result["candidates"][:MAX_CANDIDATES]
        result["truncated"] = True
    result["crosswalk"] = hits[:MAX_CROSSWALK]
    result["crosswalk_ambiguous"] = any(hit["alternatives"] > 0 for hit in hits)
    return result


@_traced("neighbors")
def neighbors(adapter, node_id: str, relations: Optional[Sequence[str]] = None,
              direction: str = "both", limit: int = 20) -> Dict[str, Any]:
    result = adapter.neighbors(node_id, relations, direction, limit)
    result = _normalize_payload(result, adapter.meta(), Provenance(adapter.meta()))
    return envelope(result)


@_traced("search_entities")
def search_entities(adapter, query: str, k: int = 10) -> List[Dict[str, Any]]:
    """Entity recall behind the platform's search tool: exact first, then text."""
    query = str(query or "").strip()
    if not query:
        return []
    nodes: List[Dict[str, Any]] = []
    seen = set()
    for candidate in api_resolve_candidates(adapter, query, k):
        if candidate["id"] not in seen:
            seen.add(candidate["id"])
            nodes.append(candidate)
    if len(nodes) < k:
        for candidate in adapter.recall(query, "text", k).get("candidates", []):
            if candidate["id"] not in seen:
                seen.add(candidate["id"])
                nodes.append(candidate)
    return nodes[:k]


def api_resolve_candidates(adapter, query: str, k: int) -> List[Dict[str, Any]]:
    result = resolve(adapter, {"text": query}, {"taxon_id": 9606})
    return result.get("candidates", [])[:k]


def validate_region(region: Any) -> Optional[Dict[str, Any]]:
    """Return an error payload when ``region`` is malformed, else ``None``."""
    if region is None:
        return None
    if not isinstance(region, dict) or len(region) != 1:
        return {"status": "invalid_region", "reason": "region must have exactly one key"}
    if not set(region) <= {"protein_range", "genomic_range"}:
        return {"status": "invalid_region", "reason": "unsupported region type"}
    interval = next(iter(region.values()))
    if (not isinstance(interval, list) or len(interval) != 2
            or any(type(pos) is not int for pos in interval)
            or not 1 <= interval[0] <= interval[1]):
        return {"status": "invalid_region", "reason": "interval must be [start, end] with 1 <= start <= end"}
    if "protein_range" in region and interval[1] - interval[0] >= 1000:
        return {"status": "invalid_region", "reason": "protein_range is limited to 1000 residues"}
    return None


@_traced("get_records")
def get_records(adapter, ids: Sequence[str], fields: Optional[Sequence[str]] = None,
                region: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    invalid = validate_region(region)
    if invalid:
        return {**envelope({**invalid, "records": [], "facts": []}), "truncated": False}

    meta = adapter.meta()
    prov = Provenance(meta)
    try:
        result = adapter.get_records(ids, fields, region)
    except (OSError, RuntimeError, sqlite3.OperationalError):
        return envelope({"status": "error", "records": [], "facts": []})

    result = _normalize_payload(result, meta, prov)
    facts: List[Dict[str, Any]] = []
    facts_truncated = False
    for record in result.get("records", []):
        record.update(envelope(record))
        if record.get("status") != "ok":
            continue
        produced, cut = _record_facts(record, meta, prov, region)
        facts.extend(produced)
        facts_truncated = facts_truncated or cut
    return {**envelope(result), "facts": facts, "facts_truncated": facts_truncated}


def _record_facts(record: Dict[str, Any], meta: Dict[str, str], prov: Provenance,
                  region: Optional[Dict[str, Any]]) -> tuple[List[Dict[str, Any]], bool]:
    node = record["node"]
    kind = node.get("kind")
    facts: List[Dict[str, Any]] = []
    scale = {"translation": "protein", "transcript": "transcript",
             "gene": "gene", "exon": "exon"}.get(kind, "gene")
    gff = _fact_source(record, prov, meta.get("source_gff3", ""), meta)
    pep = _fact_source(record, prov, meta.get("source_pep", ""), meta)
    genome = _fact_source(record, prov, meta.get("source_genome", ""), meta)

    if kind == "translation":
        sequence = (record.get("fields") or {}).get("sequence") or {}
        key = "slice" if region and "protein_range" in region else "sequence"
        value = sequence.get(key) or sequence.get("sequence")
        if value:
            fact_region = {"protein_range": list(region["protein_range"])} if key == "slice" else None
            facts.append(_fact(record, "protein_sequence", value, value_type="string",
                               evidence_type="computational", source=pep, scale="protein",
                               region=fact_region, unit="aa"))
        for pos in (record.get("region") or {}).get("positions", []):
            if pos.get("status") != "ok":
                continue
            payload = {key: pos[key] for key in pos if key != "status"}
            facts.append(_fact(
                record, "protein_position_mapping", payload, value_type="object",
                evidence_type="inferred",
                source={**gff, "field": "derived: CDS + genome"},
                scale="protein", unit="aa",
                region={"protein_range": [pos["protein_pos"], pos["protein_pos"]]},
                derived_from=[
                    {"file": meta.get("source_gff3"), "record_id": node.get("stable_id")},
                    {"file": meta.get("source_genome"), "record_id": payload.get("transcript_id")},
                ]))
    elif kind == "transcript":
        sequence = (record.get("fields") or {}).get("sequence") or {}
        if sequence.get("sequence"):
            facts.append(_fact(record, "transcript_sequence", sequence["sequence"],
                               value_type="string", evidence_type="computational",
                               source=genome, scale="transcript"))
        exons = (record.get("fields") or {}).get("exons")
        if isinstance(exons, list) and exons:
            facts.append(_fact(record, "exon_structure", [
                {"id": exon.get("id"), "start": exon.get("start"),
                 "end": exon.get("end"), "strand": exon.get("strand")}
                for exon in exons], value_type="array", evidence_type="computational",
                source=gff, scale="exon"))
        translation = (record.get("fields") or {}).get("translation")
        if isinstance(translation, dict) and translation.get("id"):
            facts.append(_fact(record, "translation", {
                "id": translation["id"], "version": translation.get("version")},
                value_type="object", evidence_type="computational",
                source=gff, scale="protein"))
    elif kind == "gene" and node.get("location"):
        facts.append(_fact(record, "gene_location", node["location"], value_type="object",
                           evidence_type="computational", source=gff, scale="gene", unit="bp"))

    xrefs = (record.get("fields") or {}).get("xrefs") or []
    truncated = len(xrefs) > MAX_XREF_FACTS
    for item in xrefs[:MAX_XREF_FACTS]:
        edge = item.get("edge") or {}
        source = item.get("source") or {}
        key = FIELD_TO_SOURCE.get(source.get("field"), meta.get("source_biomart", ""))
        facts.append(_fact(
            record, "cross_reference",
            {"accession": item["node"]["id"], "relation": edge.get("relation"),
             "role": edge.get("role")},
            value_type="object", evidence_type="computational",
            source=_fact_source(record, prov, key, meta,
                                field=source.get("field"), locator=source.get("locator")),
            scale=scale))
    return facts, truncated


@_traced("walk")
def walk(adapter, seeds: Sequence[str], relations: Optional[Sequence[str]] = None,
         max_steps: int = 6, budget: int = 30, rng_seed: int = 20261106) -> Dict[str, Any]:
    result = adapter.walk(seeds, relations, max_steps, budget, rng_seed)
    result = _normalize_payload(result, adapter.meta(), Provenance(adapter.meta()))
    return envelope(result)


@_traced("recall")
def recall(adapter, query: str, view: str = "exact", limit: int = 20) -> Dict[str, Any]:
    result = adapter.recall(query, view, limit)
    result = _normalize_payload(result, adapter.meta(), Provenance(adapter.meta()))
    return envelope(result)
