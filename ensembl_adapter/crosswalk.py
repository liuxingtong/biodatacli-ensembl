"""Load and query the cross-database key alignment table.

The crosswalk is data, not code: ``crosswalk/crosswalk.jsonl.gz`` holds one
edge per line, each carrying the BioMart row (or SIFTS record) that asserted
it, the file's SHA256 and when it was retrieved.  This module only indexes it.

Fan-out is exposed per lookup instead of hidden: a UniProt accession that
matches several chromosome-20 proteins returns every edge and reports
``alternatives`` rather than silently answering with the first hit.
"""

from __future__ import annotations

import gzip
import json
import os
import re
from collections import defaultdict
from typing import Dict, Iterable, List, Optional

UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
ENSEMBL_RE = re.compile(r"^(ENSG|ENST|ENSE|ENSP)([0-9]{6,})(?:\.([0-9]+))?$")
HGNC_RE = re.compile(r"^HGNC:(\d+)$", re.IGNORECASE)
CCDS_RE = re.compile(r"^CCDS\d+", re.IGNORECASE)
MIM_RE = re.compile(r"^\d{6}$")


def namespace_of(text: str) -> Optional[str]:
    """Best-effort namespace for a bare identifier (never guesses from a symbol)."""
    text = str(text).strip()
    if ":" in text:
        prefix, rest = text.split(":", 1)
        upper = prefix.upper()
        if upper in {"ENSEMBL", "ENSG", "ENSP", "ENST", "ENSE"}:
            return "Ensembl" if ENSEMBL_RE.match(rest) else None
        if upper in {"UNIPROTKB", "UNIPROT"}:
            return "UniProtKB"
        if upper in {"PDB", "HGNC", "CCDS", "MIM", "ENTREZGENE", "GO", "REACTOME"}:
            return {"ENTREZGENE": "EntrezGene", "REACTOME": "Reactome"}.get(upper, upper)
        return None
    if UNIPROT_RE.match(text):
        return "UniProtKB"
    if PDB_RE.match(text):
        return "PDB"
    if ENSEMBL_RE.match(text):
        return "Ensembl"
    if HGNC_RE.match(text):
        return "HGNC"
    if CCDS_RE.match(text):
        return "CCDS"
    if text.upper().startswith("R-HSA-"):
        return "Reactome"
    return None


class Crosswalk:
    """In-memory index over the edge list (≈28k edges, ~60 ms to load)."""

    def __init__(self, path: Optional[str] = None):
        self.path = path
        self.edges: List[Dict[str, object]] = []
        self.by_key: Dict[tuple, List[Dict[str, object]]] = defaultdict(list)
        self.by_relation: Dict[str, int] = defaultdict(int)
        self.build_id: Optional[str] = None
        self.release: Optional[str] = None
        if path and os.path.exists(path):
            self._load(path)
            stats_path = os.path.join(os.path.dirname(path), "stats.json")
            if os.path.exists(stats_path):
                with open(stats_path, encoding="utf-8") as handle:
                    stats = json.load(handle)
                self.build_id = stats.get("build_id")
                self.release = stats.get("release")

    def _load(self, path: str) -> None:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                edge = json.loads(line)
                self.edges.append(edge)
                self.by_relation[edge["relation"]] += 1
                for namespace, ident, other in (
                        (edge["from"]["namespace"], edge["from"]["id"], edge["to"]),
                        (edge["to"]["namespace"], edge["to"]["id"], edge["from"])):
                    self.by_key[(namespace, ident)].append({**edge, "other": other,
                                                            "anchor": {"namespace": namespace,
                                                                       "id": ident}})

    def __bool__(self) -> bool:
        return bool(self.edges)

    @property
    def available(self) -> bool:
        return bool(self.edges)

    def lookup(self, namespace: str, ident: str) -> List[Dict[str, object]]:
        """Every edge touching ``namespace:ident``, ambiguity included."""
        ident = str(ident).strip()
        if namespace != "UniProtKB":
            ident = ident.upper() if namespace in ("Ensembl", "CCDS") else ident
        hits = self.by_key.get((namespace, ident))
        if not hits and namespace == "UniProtKB":
            hits = self.by_key.get((namespace, ident.upper()))
        if not hits:
            return []
        counts = defaultdict(int)
        for hit in hits:
            counts[hit["relation"]] += 1
        out = []
        for hit in hits:
            out.append({
                "anchor": hit["anchor"],
                "other": hit["other"],
                "relation": hit["relation"],
                "cardinality": hit["cardinality"],
                "alternatives": counts[hit["relation"]] - 1,
                "evidence": hit["evidence"],
                "via": hit.get("via") or {},
            })
        out.sort(key=lambda item: (item["relation"], item["other"]["namespace"], item["other"]["id"]))
        return out

    def lookup_text(self, text: str) -> List[Dict[str, object]]:
        namespace = namespace_of(text)
        if not namespace:
            return []
        ident = str(text)
        for prefix in ("Ensembl:", "UniProtKB:", "UniProt:", "PDB:"):
            if ident.upper().startswith(prefix.upper()):
                ident = ident[len(prefix):]
                break
        version = None
        match = ENSEMBL_RE.match(ident)
        if match:
            ident, version = match.group(1) + match.group(2), match.group(3)
        hits = self.lookup(namespace, ident)
        if version:
            # a versioned ID is a stronger statement than the unversioned mapping
            for hit in hits:
                hit["requested_version"] = version
        return hits

    def stats(self) -> Dict[str, object]:
        return {
            "available": self.available,
            "path": self.path,
            "edges": len(self.edges),
            "by_relation": dict(sorted(self.by_relation.items())),
            "namespaces": sorted({edge["from"]["namespace"] for edge in self.edges} |
                                 {edge["to"]["namespace"] for edge in self.edges}),
        }
