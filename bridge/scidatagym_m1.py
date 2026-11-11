"""Expose the Ensembl chr20 slice through SciDataGym-Demo's M1 interface.

SciDataGym-Demo's M3 engine drives an adapter through six methods
(``describe / resolve / recall / neighbors / get_records / raw_record``) and
reads two attributes (``catalog``, ``manifest``).  Our Ensembl adapter answers
the same questions with different names and shapes, so this module is a thin
translation layer — no business logic, no new queries.

    from bridge.scidatagym_m1 import EnsemblM1
    m1 = EnsemblM1("ensembl_proto/data/index/ensembl_chr20_r116.sqlite")
    m1.resolve({"text": "P31946"})
    m1.neighbors("Ensembl:ENSP00000300161", ["pdb_xref"], "out", 5)

What is intentionally *not* translated:

* ``raw_record`` — the Ensembl slice indexes GFF3/FASTA/BioMart, not a single
  source record per entity; it returns ``not_supported`` instead of pretending.
* ``semantic`` view — no LSA index was built here; ``describe()["views"]``
  declares only what exists so M3 can degrade instead of crashing.
"""

from __future__ import annotations

import json
import os
import sys
from typing import Any, Dict, List, Optional, Sequence

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ensembl_adapter import EnsemblAdapter, api  # noqa: E402
from ensembl_adapter.crosswalk import namespace_of  # noqa: E402

DEFAULT_DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")

#: relation -> (source type, target type), taking our own graph schema as truth
RELATION_ENDPOINTS = {
    "has_region": ("Assembly", "Region"),
    "contains": ("Region", "Gene"),
    "has_transcript": ("Gene", "Transcript"),
    "has_exon": ("Transcript", "Exon"),
    "has_translation": ("Transcript", "Protein"),
    "uniprot_xref": ("Protein", "Protein"),
    "pdb_xref": ("Protein", "StructureEntry"),
    "hgnc_xref": ("Gene", "Gene"),
    "entrez_xref": ("Gene", "Gene"),
    "mim_xref": ("Gene", "Gene"),
    "ccds_xref": ("Transcript", "CodingSequence"),
}

TYPE_FIELDS = {
    "Gene": ["stable_id", "version", "biotype", "location", "length", "xrefs", "source"],
    "Transcript": ["stable_id", "version", "biotype", "location", "length", "exons",
                   "translation", "xrefs", "source"],
    "Protein": ["stable_id", "version", "length", "sequence", "xrefs", "source"],
    "Exon": ["stable_id", "version", "location", "length", "source"],
    "StructureEntry": ["stable_id", "xrefs", "source"],
    "CodingSequence": ["stable_id", "source"],
}

STATUS_MAP = {"ok": "resolved", "not_found": "not_found", "partial": "resolved",
              "invalid_region": "not_found", "unsupported": "not_found"}


class EnsemblM1:
    """Implementation of the SciDataGym-Demo M1 surface over ``ensembl_adapter``."""

    adapter_id = "ensembl"

    def __init__(self, db_path: str = DEFAULT_DB, genome_path: Optional[str] = None):
        self.path = db_path
        self._adapter = EnsemblAdapter(db_path, genome_path)
        meta = self._adapter.meta()
        self.manifest = {
            "database": "Ensembl",
            "release": meta.get("release"),
            "assembly": meta.get("assembly"),
            "region": meta.get("chromosome"),
            "build_id": meta.get("build_id"),
            "built_at": meta.get("built_at"),
            "schema_version": meta.get("schema_version"),
            "contract_version": api.CONTRACT_VERSION,
            "taxon_id": 9606,
            "scope_note": "human chromosome 20 only",
        }
        self.catalog = self._build_catalog()

    # ------------------------------------------------------------------
    # contract surface
    # ------------------------------------------------------------------
    def describe(self) -> Dict[str, Any]:
        caps = api.capabilities(self._adapter)
        return {
            "adapter_id": self.adapter_id,
            "schema_version": "1.0",
            "manifest": self.manifest,
            "types": sorted(self.catalog["types"]),
            "relations": {r: list(pair) for r, pair in self.catalog["relations"].items()},
            "views": ["identifier", "field"],
            "unsupported": self.catalog["unsupported_operations"],
            "notes": [
                "graph slice: human chromosome 20, Ensembl release %s" % self.manifest["release"],
                "no semantic (LSA) index in this slice",
                "raw_record is not offered: entities are derived from GFF3/FASTA/BioMart",
                "facts carry db/id/version/retrieved_at per contract m1.1",
            ],
            "capabilities": caps,
        }

    def resolve(self, ref: Any, constraints: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        result = api.resolve(self._adapter, ref, constraints)
        items = [self._node(node) for node in result.get("candidates", [])]
        return {
            "items": items,
            "status": "not_found" if not items else ("resolved" if len(items) == 1 else "ambiguous"),
            "truncated": bool(result.get("truncated")),
            "null_status": result.get("null_status"),
        }

    def recall(self, query: Any, view: str, limit: int = 20) -> Dict[str, Any]:
        """``query`` is ``{text, target_type, constraints}`` like the UniProt adapter."""
        if view not in ("identifier", "field"):
            raise ValueError(f"Unsupported view: {view}")
        query = {"text": query} if isinstance(query, str) else (query or {})
        text = str(query.get("text") or "").strip()
        target = query.get("target_type")
        constraints = query.get("constraints") or {}

        if view == "identifier" or not text:
            resolved = self.resolve({"text": text, "type": target}, constraints)
            items = [{"node": node, "score": 1.0, "rank": index + 1, "view": view}
                     for index, node in enumerate(resolved["items"][:limit])]
            return {"items": items, "status": resolved["status"], "truncated": False}

        result = api.recall(self._adapter, text, "text", limit)
        nodes = [self._node(node) for node in result.get("candidates", [])]
        if target:
            nodes = [node for node in nodes if node["type"] == target]
        items = [{"node": node, "score": 1.0 - index / max(1, len(nodes)),
                  "rank": index + 1, "view": view} for index, node in enumerate(nodes[:limit])]
        return {"items": items, "status": "ok" if items else "not_found", "truncated": False}

    def neighbors(self, node_id: str, relations: Optional[Sequence[str]] = None,
                  direction: str = "out", limit: int = 100, cursor: int = 0) -> Dict[str, Any]:
        if direction not in ("out", "in", "both") or cursor < 0:
            raise ValueError("Invalid direction/cursor")
        known = set(self.catalog["relations"])
        unknown = [relation for relation in (relations or []) if relation not in known]
        if unknown:
            raise ValueError(f"Unsupported relation: {unknown[0]}")

        result = api.neighbors(self._adapter, node_id, relations, direction, 100000)
        if not result.get("ok"):
            return {"items": [], "status": "node_not_found", "truncated": False,
                    "next_cursor": None}
        every = [self._edge_item(item, node_id) for item in result["items"]]
        page = every[cursor:cursor + limit]
        truncated = cursor + limit < len(every)
        return {
            "items": page,
            "status": "ok",
            "truncated": truncated,
            "next_cursor": cursor + limit if truncated else None,
            "total": result.get("total"),
        }

    def get_records(self, ids: Sequence[str], fields: Optional[Sequence[str]] = None,
                    region: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        region_arg = None
        if region:
            start, end = region.get("start"), region.get("end")
            if not isinstance(start, int) or not isinstance(end, int) or start < 1 or end < start:
                raise ValueError("Region must be {start, end} with 1 <= start <= end")
            region_arg = {"protein_range": [start, end]}

        result = api.get_records(self._adapter, list(ids), list(fields) if fields else None,
                                 region_arg)
        items: List[Dict[str, Any]] = []
        for record in result.get("records", []):
            if record.get("status") != "ok":
                items.append({"id": record.get("id"), "status": record.get("status", "not_found"),
                              "null_status": record.get("null_status")})
                continue
            data = dict(record.get("fields") or {})
            sequence_field = data.get("sequence")
            if isinstance(sequence_field, dict):
                # the contract exposes {sequence, length, slice}; M3 expects the string
                data["length"] = sequence_field.get("length", data.get("length"))
                data["sequence"] = (sequence_field.get("slice") if region_arg
                                    else sequence_field.get("sequence"))
            if region_arg:
                data["region"] = {"start": region["start"], "end": region["end"],
                                  "indexing": "one_based_closed",
                                  "sequence_version": record["node"].get("version")}
            data.setdefault("length", record["node"].get("length"))
            data.setdefault("source", record.get("source"))
            items.append({
                "id": record["id"],
                "type": record["type"],
                "data": data,
                "source": record.get("source"),
                "status": "ok",
            })
        return {"items": items,
                "facts": result.get("facts", []),
                "status": "ok" if all(item["status"] == "ok" for item in items) else "partial"}

    def raw_record(self, accession: str) -> Dict[str, Any]:
        return {
            "status": "not_supported",
            "reason": "the Ensembl slice has no single source record per entity",
            "hint": "use get_records(fields=['source']) for the GFF3/FASTA/BioMart locator",
        }

    def walk(self, seeds: Sequence[str], relations: Sequence[str], max_steps: int = 3,
             budget: int = 30, random_seed: int = 0) -> Dict[str, Any]:
        result = api.walk(self._adapter, list(seeds), list(relations), max_steps, budget,
                          random_seed)
        return {"walks": result.get("walks", []), "status": result.get("status"),
                "rng_seed": result.get("rng_seed")}

    # ------------------------------------------------------------------
    # shape translation
    # ------------------------------------------------------------------
    def _node(self, node: Dict[str, Any]) -> Dict[str, Any]:
        native = node.get("source") or {}
        properties = {
            "stable_id": node.get("stable_id"),
            "version": node.get("version"),
            "biotype": node.get("biotype"),
            "length": node.get("length"),
            "location": node.get("location"),
            "label": node.get("label"),
            "source": native,
        }
        if node.get("type") == "StructureEntry":
            properties["structure_id"] = node.get("stable_id")
        return {
            "id": node["id"],
            "type": node.get("type"),
            "label": node.get("label") or node.get("stable_id") or node["id"],
            "taxon_id": 9606 if node.get("type") != "StructureEntry" else None,
            "properties": {key: value for key, value in properties.items() if value is not None},
        }

    def _edge_item(self, item: Dict[str, Any], anchor: str) -> Dict[str, Any]:
        edge = item.get("edge") or {}
        return {
            "node": self._node(item["node"]),
            "edge": {
                "source": edge.get("from"),
                "relation": edge.get("relation"),
                "target": edge.get("to"),
                "properties": {key: value for key, value in edge.items()
                               if key not in ("from", "relation", "to")},
            },
            "source": item.get("source") or {},
            "direction": item.get("direction"),
        }

    def _build_catalog(self) -> Dict[str, Any]:
        caps = api.capabilities(self._adapter)
        relations = {}
        for relation in caps.get("relations", []):
            pair = RELATION_ENDPOINTS.get(relation)
            if pair:
                relations[relation] = pair
        return {
            "version": "1.0",
            "types": TYPE_FIELDS,
            "relations": relations,
            "field_aliases": {},
            "op_aliases": {},
            "taxon_aliases": {"human": 9606, "homo sapiens": 9606},
            "method_aliases": {},
            "notes": ["translated from ensembl capabilities.json by bridge/scidatagym_m1.py"],
            "unsupported_operations": [
                "raw_record", "semantic_view", "go_annotation", "regulation",
                "variation", "non_human_taxa",
            ],
        }


if __name__ == "__main__":
    m1 = EnsemblM1()
    print(json.dumps(m1.describe(), indent=2, ensure_ascii=False)[:900])
