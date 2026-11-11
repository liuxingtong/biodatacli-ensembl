"""MCP tool surface exactly as the platform specifies it.

The platform's required tools are ``describe / resolve / search / get_facts``;
the optional ones are ``search_datasets / get_dataset_metadata / get_access_plan
/ expand``.  Our internal API (``capabilities / resolve / neighbors /
get_records / recall / datasets``) answers the same questions with different
names and payloads, so this module is the translation layer, mirroring what
``bridge/scidatagym_m1.py`` does for SciDataGym-Demo's M1.

Nothing here invents data: every field comes from the graph or from the frozen
input provenance.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Sequence

from . import api
from .crosswalk import namespace_of

LICENSE = ("Ensembl data are produced by EMBL-EBI and freely available; see "
           "https://www.ensembl.org/info/about/legal/disclaimer.html")

#: fields of get_records that produce facts, keyed by the scale they describe
FIELDS_BY_SCALE = {
    "gene": ["attributes", "xrefs"],
    "transcript": ["sequence", "exons", "translation", "xrefs"],
    "exon": ["attributes"],
    "protein": ["sequence", "protein_mapping", "xrefs"],
}


def describe(adapter) -> Dict[str, Any]:
    """Capability declaration in the platform's shape."""
    caps = api.capabilities(adapter)
    return {
        "adapter_id": f"ensembl-r{caps.get('release')}-{caps.get('scope', {}).get('region', '')}",
        "entities": caps.get("node_types", []),
        "scales": caps.get("scales", []),
        "keyspaces": caps.get("namespaces", []),
        "evidence_types": caps.get("evidence_types", []),
        "search_modes": caps.get("search_modes", ["identifier", "text"]),
        "license": LICENSE,
        "release": caps.get("release"),
        "assembly": caps.get("assembly"),
        "snapshot": caps.get("snapshot"),
        "null_semantics": caps.get("null_semantics"),
        "key_alignment": caps.get("key_alignment"),
        "fact_attributes": caps.get("fact_attributes", []),
        "crosswalk": {k: caps["crosswalk"][k] for k in ("edges", "build_id")
                      if k in (caps.get("crosswalk") or {})},
    }


def resolve(adapter, term: str, context: Optional[Dict[str, Any]] = None,
            species: Optional[str] = None) -> Dict[str, Any]:
    """Free text / alias -> canonical key."""
    constraints: Dict[str, Any] = {}
    context = context or {}
    if species and str(species).lower() not in ("human", "homo sapiens", "9606"):
        return {"candidates": [], "status": "unsupported", "null_status": "na",
                "reason": "only human (taxon 9606) is indexed"}
    if context.get("type"):
        constraints["type"] = None  # type is carried inside the ref, not constraints
    ref: Dict[str, Any] = {"text": term}
    if context.get("type"):
        ref["type"] = context["type"]

    result = api.resolve(adapter, ref, {"taxon_id": 9606})
    out = []
    for rank, node in enumerate(result.get("candidates", [])):
        namespace, _, identifier = str(node["id"]).partition(":")
        aliases = [value for value in (node.get("label"), node.get("stable_id"),
                                       node.get("global_id")) if value and value != identifier]
        out.append({
            "canonical_id": identifier,
            "namespace": namespace,
            "type": node.get("type"),
            "score": round(max(0.1, 1.0 - rank * 0.1), 2),
            "aliases": sorted(set(aliases)),
        })
    return {"candidates": out, "status": result.get("status"),
            "null_status": result.get("null_status"),
            "ambiguous": bool(result.get("crosswalk_ambiguous"))}


def search(adapter, query: str, filters: Optional[Dict[str, Any]] = None,
           k: int = 10) -> Dict[str, Any]:
    """Retrieve facts, not just entities: entity + attribute + value + evidence."""
    filters = filters or {}
    attribute = filters.get("attribute")
    scale = filters.get("scale")
    kinds = filters.get("entity_type")

    hits = api.search_entities(adapter, query, k)
    facts: List[Dict[str, Any]] = []
    for rank, node in enumerate(hits):
        if kinds and node.get("type") not in (kinds if isinstance(kinds, list) else [kinds]):
            continue
        record = api.get_records(adapter, [node["id"]],
                                 fields=filters.get("fields"))
        for fact in record.get("facts", []):
            if attribute and fact["attribute"] != attribute:
                continue
            if scale and fact.get("scale") != scale:
                continue
            facts.append({
                "entity": fact["entity"],
                "attribute": fact["attribute"],
                "value": fact["value"],
                "evidence_type": fact["evidence_type"],
                "source": fact["source"],
                "scale": fact.get("scale"),
                "score": round(max(0.1, 1.0 - rank * 0.1), 2),
            })
            if len(facts) >= k:
                break
        if len(facts) >= k:
            break
    return {"facts": facts, "status": "ok" if facts else "not_found",
            "null_status": None if facts else "not_found",
            "total": len(facts)}


def get_facts(adapter, entity_ids: Sequence[str], fields: Optional[Sequence[str]] = None,
              scale: Optional[str] = None) -> Dict[str, Any]:
    """Structured facts for known ids, plus what was and was not covered."""
    record = api.get_records(adapter, list(entity_ids), list(fields) if fields else None)
    facts = record.get("facts", [])
    if scale:
        facts = [fact for fact in facts if fact.get("scale") == scale]

    answered = {fact["entity"]["id"] for fact in facts}
    requested = [str(identifier).split(":")[-1] for identifier in entity_ids]
    missing = [identifier for identifier in requested if identifier not in answered]
    return {
        "facts": facts,
        "coverage": {
            "requested": len(requested),
            "entities_with_facts": len(answered),
            "missing": missing,
            "attributes": sorted({fact["attribute"] for fact in facts}),
            "null_status": record.get("null_status"),
        },
    }


# ----------------------------------------------------------------------
# optional capabilities
# ----------------------------------------------------------------------
def search_datasets(adapter, query: str = "", filters: Optional[Dict[str, Any]] = None,
                    k: int = 10) -> Dict[str, Any]:
    filters = filters or {}
    records = api.datasets(adapter)["datasets"]
    needle = (query or "").lower()

    def matches(record: Dict[str, Any]) -> bool:
        if needle and needle not in (record["name"] + " " + record["description"]).lower():
            return False
        if filters.get("format") and record["format"] != filters["format"]:
            return False
        if filters.get("access") and record["access"]["kind"] != filters["access"]:
            return False
        return True

    hits = [record for record in records if matches(record)][:k]
    return {"datasets": hits, "total": len(hits),
            "status": "ok" if hits else "not_found",
            "null_status": None if hits else "not_found"}


def get_dataset_metadata(adapter, dataset_id: str) -> Dict[str, Any]:
    for record in api.datasets(adapter)["datasets"]:
        if record["dataset_id"] == dataset_id or record["name"] == dataset_id:
            return {"dataset": record, "status": "ok", "null_status": None}
    return {"dataset": None, "status": "not_found", "null_status": "not_found"}


def get_access_plan(adapter, dataset_id: str, target: str = "on_demand") -> Dict[str, Any]:
    if target != "on_demand":
        return {"status": "unsupported", "null_status": "na",
                "reason": "only on_demand access is supported by this adapter"}
    record = get_dataset_metadata(adapter, dataset_id)
    if not record.get("dataset"):
        return {"status": "not_found", "null_status": "not_found"}
    dataset = record["dataset"]
    return {
        "dataset_id": dataset["dataset_id"],
        "fetch_url_or_cmd": dataset["access"]["fetch"],
        "bytes": dataset["size"],
        "format": dataset["format"],
        "license": LICENSE,
        "status": "ok",
        "null_status": None,
        "kind": dataset["access"]["kind"],
    }


def expand(adapter, entity_ids: Sequence[str], key: Sequence[str] | str,
           max_depth: int = 1, k: int = 50) -> Dict[str, Any]:
    """Graph expansion across scales, flattened to from/to/relation edges."""
    relations = [key] if isinstance(key, str) else list(key)
    frontier = list(entity_ids)
    seen_nodes = set(frontier)
    edges: List[Dict[str, Any]] = []
    for _depth in range(max(1, min(max_depth, 3))):
        next_frontier: List[str] = []
        for entity in frontier:
            page = api.neighbors(adapter, entity, relations, "both", k)
            for item in page.get("items", []):
                edge = item.get("edge") or {}
                edges.append({
                    "from": edge.get("from"),
                    "to": edge.get("to"),
                    "relation": edge.get("relation"),
                    "evidence": item.get("source") or {},
                    "score": 1.0,
                })
                other = item["node"]["id"]
                if other not in seen_nodes:
                    seen_nodes.add(other)
                    next_frontier.append(other)
                if len(edges) >= k:
                    break
            if len(edges) >= k:
                break
        if not next_frontier or len(edges) >= k:
            break
        frontier = next_frontier
    return {"edges": edges[:k], "total": len(edges),
            "nodes": sorted(seen_nodes), "status": "ok" if edges else "not_found",
            "null_status": None if edges else "not_found"}
