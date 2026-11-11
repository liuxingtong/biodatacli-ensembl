"""Read-only M1-contract query layer over the Ensembl graph index.

The public surface follows the SciDataGym M1 adapter contract:

    resolve(ref, constraints)
    neighbors(node_id, relations, direction, limit)
    get_records(ids, fields, region)
    walk(seeds, relations, max_steps, budget)
    recall(query, view, limit)
    capabilities()

Every returned neighbour carries the node, the typed edge (with attributes)
and the source record that asserted the relation.
"""

from __future__ import annotations

import json
import os
import random
import re
import sqlite3
from typing import Dict, List, Optional, Sequence

from .crosswalk import Crosswalk
from .fasta import load_genome
from .sequence import reverse_complement, splice, translate

UNIPROT_RE = re.compile(r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")
PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
STABLE_RE = re.compile(r"^(ENSG|ENST|ENSE|ENSP)([0-9]{6,})(?:\.([0-9]+))?$")

TYPE_PRIORITY = {
    "Assembly": 0, "Region": 1, "Gene": 2, "Transcript": 3,
    "Exon": 4, "Protein": 5, "CodingSequence": 6, "StructureEntry": 7,
}

NODE_COLS_SRC = (
    "s.node_id AS src_node_id, s.key AS src_key, s.global_id AS src_global_id, "
    "s.kind AS src_kind, s.node_type AS src_node_type, s.stable_id AS src_stable_id, "
    "s.version AS src_version, s.label AS src_label, s.biotype AS src_biotype, "
    "s.seq_region AS src_seq_region, s.start AS src_start, s.end AS src_end, "
    "s.strand AS src_strand, s.length AS src_length, s.extra AS src_extra, "
    "s.source_db AS src_source_db, s.source_record_id AS src_source_record_id, "
    "s.source_version AS src_source_version, s.source_field AS src_source_field, "
    "s.source_locator AS src_source_locator"
)
NODE_COLS_DST = (
    "d.node_id AS dst_node_id, d.key AS dst_key, d.global_id AS dst_global_id, "
    "d.kind AS dst_kind, d.node_type AS dst_node_type, d.stable_id AS dst_stable_id, "
    "d.version AS dst_version, d.label AS dst_label, d.biotype AS dst_biotype, "
    "d.seq_region AS dst_seq_region, d.start AS dst_start, d.end AS dst_end, "
    "d.strand AS dst_strand, d.length AS dst_length, d.extra AS dst_extra, "
    "d.source_db AS dst_source_db, d.source_record_id AS dst_source_record_id, "
    "d.source_version AS dst_source_version, d.source_field AS dst_source_field, "
    "d.source_locator AS dst_source_locator"
)
EDGE_COLS = (
    "e.relation, e.rank, e.attrs, "
    "e.source_db AS e_source_db, e.source_record_id AS e_source_record_id, "
    "e.source_version AS e_source_version, e.source_field AS e_source_field, "
    "e.source_locator AS e_source_locator"
)


class EnsemblAdapter:
    def __init__(self, db_path: str, genome_path: Optional[str] = None,
                 crosswalk_path: Optional[str] = None):
        self.db_path = db_path
        uri = "file:" + os.path.abspath(db_path) + "?mode=ro"
        self.conn = sqlite3.connect(uri, uri=True)
        self.conn.row_factory = sqlite3.Row
        self._genome: Optional[Dict[str, str]] = None
        self._genome_path = genome_path or self._locate_genome()
        self._capabilities: Optional[Dict[str, object]] = None
        self._crosswalk_path = crosswalk_path or self._locate_crosswalk()
        self._crosswalk: Optional[Crosswalk] = None

    # ------------------------------------------------------------------
    # Low-level helpers
    # ------------------------------------------------------------------
    def _locate_genome(self) -> Optional[str]:
        base = self.meta().get("source_genome")
        if not base:
            return None
        candidates = [
            os.path.join(os.path.dirname(os.path.abspath(self.db_path)), "..", "raw", base),
            os.path.join(os.path.dirname(os.path.abspath(self.db_path)), base),
            os.path.join(os.getcwd(), "ensembl_proto", "data", "raw", base),
            os.path.join(os.getcwd(), "data", "raw", base),
        ]
        for candidate in candidates:
            candidate = os.path.normpath(candidate)
            if os.path.exists(candidate):
                return candidate
        return None

    def _locate_crosswalk(self) -> Optional[str]:
        root = os.path.dirname(os.path.dirname(os.path.abspath(self.db_path)))
        candidates = [
            os.path.join(root, "crosswalk", "crosswalk.jsonl.gz"),
            os.path.join(os.getcwd(), "crosswalk", "crosswalk.jsonl.gz"),
        ]
        for candidate in candidates:
            if os.path.exists(candidate):
                return candidate
        return None

    @property
    def crosswalk(self) -> Crosswalk:
        """Key alignment table; empty (and falsey) when no build is present."""
        if self._crosswalk is None:
            self._crosswalk = Crosswalk(self._crosswalk_path)
        return self._crosswalk

    def crosswalk_lookup(self, text: str):
        return self.crosswalk.lookup_text(text) if self.crosswalk.available else []

    def meta(self) -> Dict[str, str]:
        return {row["key"]: row["value"] for row in self.conn.execute("SELECT key,value FROM meta")}

    def genome(self) -> Dict[str, str]:
        if self._genome is None:
            if not self._genome_path:
                raise RuntimeError("genome FASTA not found; pass genome_path=")
            self._genome = load_genome(self._genome_path)
        return self._genome

    @staticmethod
    def _json(text):
        try:
            return json.loads(text) if text else {}
        except (TypeError, json.JSONDecodeError):
            return {}

    @staticmethod
    def _node_source(row: sqlite3.Row, prefix: str = "") -> Dict[str, object]:
        return {
            "database": row[prefix + "source_db"],
            "record_id": row[prefix + "source_record_id"],
            "release": row[prefix + "source_version"],
            "field": row[prefix + "source_field"],
            "locator": row[prefix + "source_locator"],
        }

    def _row_to_node(self, row: sqlite3.Row, prefix: str = "") -> Dict[str, object]:
        extra = self._json(row[prefix + "extra"])
        seq_region = row[prefix + "seq_region"]
        location = None
        if seq_region is not None:
            location = {
                "seq_region": seq_region,
                "start": row[prefix + "start"],
                "end": row[prefix + "end"],
                "strand": row[prefix + "strand"],
            }
        return {
            "id": row[prefix + "global_id"],
            "type": row[prefix + "node_type"],
            "kind": row[prefix + "kind"],
            "stable_id": row[prefix + "stable_id"],
            "version": row[prefix + "version"],
            "label": row[prefix + "label"],
            "biotype": row[prefix + "biotype"],
            "location": location,
            "length": row[prefix + "length"],
            "attributes": extra,
            "source": self._node_source(row, prefix),
            # flat aliases kept for older call sites / CLI tables
            "key": row[prefix + "key"],
            "global_id": row[prefix + "global_id"],
            "seq_region": seq_region,
            "start": row[prefix + "start"],
            "end": row[prefix + "end"],
            "strand": row[prefix + "strand"],
            "extra": extra,
        }

    @staticmethod
    def _edge_record(row: sqlite3.Row) -> Dict[str, object]:
        edge = {
            "from": row["src_global_id"],
            "relation": row["relation"],
            "to": row["dst_global_id"],
            "rank": row["rank"],
        }
        edge.update(EnsemblAdapter._json(row["attrs"]))
        return edge

    @staticmethod
    def _edge_source(row: sqlite3.Row) -> Dict[str, object]:
        return {
            "database": row["e_source_db"],
            "record_id": row["e_source_record_id"],
            "release": row["e_source_version"],
            "field": row["e_source_field"],
            "locator": row["e_source_locator"],
        }

    @staticmethod
    def _stable_id(value: str) -> str:
        value = str(value).strip()
        if ":" in value:
            value = value.split(":", 1)[1]
        return value.split(".")[0]

    def _node_row(self, node_id: str) -> Optional[sqlite3.Row]:
        if node_id is None:
            return None
        node_id = str(node_id).strip()
        row = self.conn.execute(
            "SELECT * FROM node WHERE global_id=? OR key=?", (node_id, node_id)
        ).fetchone()
        if row:
            return row

        stable = self._stable_id(node_id)
        row = self.conn.execute(
            "SELECT * FROM node WHERE stable_id=? ORDER BY "
            "CASE kind WHEN 'gene' THEN 0 WHEN 'transcript' THEN 1 WHEN 'exon' THEN 2 "
            "WHEN 'translation' THEN 3 ELSE 4 END, node_id",
            (stable,),
        ).fetchone()
        if row:
            return row

        if PDB_RE.match(node_id):
            return self.conn.execute(
                "SELECT * FROM node WHERE global_id=?", (f"PDB:{node_id.upper()}",)
            ).fetchone()
        if UNIPROT_RE.match(node_id):
            return self.conn.execute(
                "SELECT * FROM node WHERE global_id=?", (f"UniProtKB:{node_id}",)
            ).fetchone()
        if node_id.upper().startswith(("ENSG", "ENST", "ENSE", "ENSP")):
            return self.conn.execute(
                "SELECT * FROM node WHERE global_id=?", (f"Ensembl:{node_id.upper()}",)
            ).fetchone()

        return self.conn.execute(
            "SELECT * FROM node WHERE label=? COLLATE NOCASE AND kind IN "
            "('gene','transcript','translation','exon') "
            "ORDER BY CASE kind WHEN 'gene' THEN 0 WHEN 'transcript' THEN 1 "
            "WHEN 'exon' THEN 2 ELSE 3 END, node_id",
            (node_id,),
        ).fetchone()

    def get(self, node_id: str) -> Optional[Dict[str, object]]:
        row = self._node_row(node_id)
        return self._row_to_node(row) if row else None

    def node_by_id(self, node_id: int) -> Optional[Dict[str, object]]:
        row = self.conn.execute("SELECT * FROM node WHERE node_id=?", (node_id,)).fetchone()
        return self._row_to_node(row) if row else None

    def _type_matches(self, node: Dict[str, object], type_filter) -> bool:
        if not type_filter:
            return True
        allowed = [type_filter] if isinstance(type_filter, str) else list(type_filter)
        return node["type"] in allowed

    # ------------------------------------------------------------------
    # M1 contract: resolve / neighbors / get_records / walk / recall
    # ------------------------------------------------------------------
    def resolve(self, ref, constraints: Optional[Dict[str, object]] = None) -> Dict[str, object]:
        if isinstance(ref, str):
            ref = {"text": ref}
        if not isinstance(ref, dict):
            return {"candidates": [], "status": "invalid_ref", "truncated": False}
        text = str(ref.get("text") or ref.get("id") or "").strip()
        type_filter = ref.get("type")
        constraints = constraints or {}

        if not text:
            return {"candidates": [], "status": "not_found", "truncated": False}
        if constraints.get("taxon_id") not in (None, 9606):
            return {"candidates": [], "status": "unsupported",
                    "reason": "only human (taxon_id=9606) is indexed", "truncated": False}
        if constraints.get("release") is not None and str(constraints["release"]) != self.meta().get("release"):
            return {"candidates": [], "status": "unsupported",
                    "reason": f"index release is {self.meta().get('release')}", "truncated": False}
        if constraints.get("assembly") not in (None, self.meta().get("assembly")):
            return {"candidates": [], "status": "unsupported",
                    "reason": f"index assembly is {self.meta().get('assembly')}", "truncated": False}

        candidates: List[Dict[str, object]] = []
        seen = set()

        def add(row, matched_by):
            if row is None:
                return
            node = self._row_to_node(row)
            if node["id"] in seen:
                return
            if not self._type_matches(node, type_filter):
                return
            if constraints.get("biotype") and node["biotype"] != constraints["biotype"]:
                return
            node["matched_by"] = matched_by
            seen.add(node["id"])
            candidates.append(node)

        add(self.conn.execute(
            "SELECT * FROM node WHERE global_id=? OR key=?", (text, text)
        ).fetchone(), "id")

        stable = self._stable_id(text)
        rows = self.conn.execute(
            "SELECT * FROM node WHERE stable_id=? ORDER BY "
            "CASE kind WHEN 'gene' THEN 0 WHEN 'transcript' THEN 1 WHEN 'exon' THEN 2 "
            "WHEN 'translation' THEN 3 ELSE 4 END, node_id",
            (stable,),
        ).fetchall()
        for row in rows:
            add(row, "stable_id")

        if PDB_RE.match(text):
            add(self.conn.execute("SELECT * FROM node WHERE global_id=?",
                                  (f"PDB:{text.upper()}",)).fetchone(), "pdb_id")
        if UNIPROT_RE.match(text):
            add(self.conn.execute("SELECT * FROM node WHERE global_id=?",
                                  (f"UniProtKB:{text}",)).fetchone(), "uniprot_accession")

        rows = self.conn.execute(
            "SELECT * FROM node WHERE label=? COLLATE NOCASE AND kind IN "
            "('gene','transcript','translation','exon')",
            (text,),
        ).fetchall()
        for row in rows:
            add(row, "symbol")

        candidates.sort(key=lambda n: (TYPE_PRIORITY.get(n["type"], 99), n["id"]))
        return {"candidates": candidates,
                "status": "ok" if candidates else "not_found",
                "truncated": False}

    def neighbors(self, node_id: str, relations: Optional[Sequence[str]] = None,
                  direction: str = "both", limit: int = 20) -> Dict[str, object]:
        node = self.get(node_id)
        if not node:
            return {"items": [], "status": "not_found", "truncated": False, "total": 0}
        if direction not in ("out", "in", "both"):
            return {"items": [], "status": "invalid_direction", "truncated": False, "total": 0}

        row = self.conn.execute("SELECT node_id FROM node WHERE key=?", (node["key"],)).fetchone()
        nid = row["node_id"]
        where = []
        params: List[object] = []
        if direction == "out":
            where.append("e.src=?")
            params.append(nid)
        elif direction == "in":
            where.append("e.dst=?")
            params.append(nid)
        else:
            where.append("(e.src=? OR e.dst=?)")
            params.extend([nid, nid])
        if relations:
            rels = [relations] if isinstance(relations, str) else list(relations)
            where.append("e.relation IN (" + ",".join("?" for _ in rels) + ")")
            params.extend(rels)

        sql = (
            "SELECT " + EDGE_COLS + ", " + NODE_COLS_SRC + ", " + NODE_COLS_DST +
            " FROM edge e JOIN node s ON s.node_id=e.src JOIN node d ON d.node_id=e.dst"
            " WHERE " + " AND ".join(where) +
            " ORDER BY e.relation, d.global_id, s.global_id"
        )
        rows = self.conn.execute(sql, params).fetchall()
        items = []
        for r in rows:
            is_out = r["src_global_id"] == node["id"]
            other = self._row_to_node(r, "dst_") if is_out else self._row_to_node(r, "src_")
            items.append({
                "node": other,
                "edge": self._edge_record(r),
                "source": self._edge_source(r),
                "direction": "out" if is_out else "in",
            })
        items.sort(key=lambda item: (item["edge"]["relation"], item["node"]["id"]))
        total = len(items)
        if limit is not None:
            items = items[:max(0, int(limit))]
        return {"items": items, "status": "ok", "total": total,
                "truncated": total > len(items), "limit": limit}

    def _protein_transcript(self, protein_id: str) -> Optional[Dict[str, object]]:
        stable = self._stable_id(protein_id)
        row = self.conn.execute(
            "SELECT n.* FROM edge e JOIN node n ON n.node_id=e.src "
            "WHERE e.dst=(SELECT node_id FROM node WHERE stable_id=? AND kind='translation') "
            "AND e.relation='has_translation' ORDER BY n.node_id LIMIT 1",
            (stable,),
        ).fetchone()
        return self._row_to_node(row) if row else None

    def _default_fields(self, node: Dict[str, object]) -> List[str]:
        node_type = node["type"]
        if node_type == "Protein":
            return ["attributes", "sequence", "xrefs"]
        if node_type == "Transcript":
            return ["attributes", "sequence", "exons", "translation", "xrefs"]
        if node_type == "Gene":
            return ["attributes", "xrefs"]
        return ["attributes"]

    def get_records(self, ids, fields: Optional[Sequence[str]] = None,
                    region: Optional[Dict[str, object]] = None) -> Dict[str, object]:
        if isinstance(ids, str):
            ids = [ids]
        records = []
        for ident in ids:
            node = self.get(ident)
            if not node:
                records.append({"id": str(ident), "status": "not_found"})
                continue
            records.append(self._record(node, fields, region))
        status = "ok" if all(r.get("status") == "ok" for r in records) else "partial"
        return {"records": records, "status": status, "truncated": False}

    def _record(self, node: Dict[str, object], fields: Optional[Sequence[str]],
                region: Optional[Dict[str, object]]) -> Dict[str, object]:
        fields = list(fields) if fields else self._default_fields(node)
        record = {
            "id": node["id"],
            "type": node["type"],
            "node": node,
            "fields": {},
            "region": None,
            "source": node["source"],
            "status": "ok",
        }
        for field in fields:
            if field == "attributes":
                record["fields"]["attributes"] = node["attributes"]
            elif field == "sequence":
                record["fields"]["sequence"] = self._sequence_for(node, region)
            elif field == "exons":
                if node["type"] == "Transcript":
                    record["fields"]["exons"] = self.transcript_exons(node["stable_id"])
                else:
                    record["fields"]["exons"] = {"status": "unsupported"}
            elif field == "translation":
                if node["type"] == "Transcript":
                    record["fields"]["translation"] = self.transcript_translation(node["stable_id"])
                else:
                    record["fields"]["translation"] = {"status": "unsupported"}
            elif field == "xrefs":
                xref_relations = [r for r in self.capabilities().get("relations", []) if r.endswith("_xref")]
                record["fields"]["xrefs"] = self.neighbors(
                    node["id"], relations=xref_relations, direction="both", limit=1000
                )["items"]
            elif field in ("protein_mapping", "residue_mapping", "coverage"):
                record["fields"][field] = self._mapping_for(node, region)
            else:
                record["fields"][field] = {"status": "unsupported",
                                           "reason": f"unknown field {field}"}
        if region:
            record["region"] = self._mapping_for(node, region)
        return record

    def _sequence_for(self, node: Dict[str, object], region) -> Dict[str, object]:
        if node["type"] == "Protein":
            seq = self.protein_seq(node["stable_id"])
            if not seq:
                return {"status": "unsupported", "reason": "no sequence stored for this node"}
            out = {"sequence": seq, "length": len(seq)}
            if region and "protein_range" in region and seq:
                start, end = region["protein_range"]
                out["slice"] = seq[max(0, start - 1):end]
            return out
        if node["type"] == "Transcript":
            seq = self.cdna_seq(node["stable_id"])
            return {"sequence": seq, "length": len(seq)}
        return {"status": "unsupported", "reason": f"sequence not defined for {node['type']}"}

    def _mapping_for(self, node: Dict[str, object], region: Optional[Dict[str, object]]) -> Dict[str, object]:
        if not region:
            return {"status": "unsupported", "reason": "region is required"}
        if "protein_range" not in region:
            if "genomic_range" in region:
                start, end = region["genomic_range"]
                loc = node.get("location") or {}
                if loc and loc.get("seq_region"):
                    overlap = not (loc["end"] < start or loc["start"] > end)
                    return {"type": "genomic_range", "requested": [start, end],
                            "overlap": overlap, "status": "ok"}
                return {"status": "unsupported", "reason": "node has no genomic location"}
            return {"status": "unsupported", "reason": "unsupported region type"}

        start, end = [int(x) for x in region["protein_range"]]
        if node["type"] == "Gene":
            cov = self.gene_coverage(node["stable_id"], start, end)
            cov["type"] = "protein_range"
            cov["requested"] = [start, end]
            return cov
        if node["type"] == "Transcript":
            cov = self.transcript_coverage(node["stable_id"], start, end)
            positions = cov.get("positions", [])
            cov["type"] = "protein_range"
            cov["requested"] = [start, end]
            cov["retrieved"] = [p["protein_pos"] for p in positions if p["status"] == "ok"]
            cov["missing"] = [{"position": p["protein_pos"], "reason": p["status"]}
                              for p in positions if p["status"] != "ok"]
            return cov
        if node["type"] == "Protein":
            tx = self._protein_transcript(node["stable_id"])
            if not tx:
                return {"status": "unsupported", "reason": "no transcript for protein"}
            positions = []
            for pos in range(start, end + 1):
                positions.append(self.map_protein_pos(tx["stable_id"], pos))
            retrieved = [p for p in positions if p["status"] == "ok"]
            missing = [{"position": p["protein_pos"], "reason": p["status"]}
                       for p in positions if p["status"] != "ok"]
            seq = self.protein_seq(node["stable_id"])
            return {
                "type": "protein_range", "requested": [start, end],
                "transcript_id": tx["stable_id"],
                "retrieved": [p["protein_pos"] for p in retrieved],
                "missing": missing,
                "positions": positions,
                "sequence_slice": seq[max(0, start - 1):end],
                "status": "ok" if not missing else "partial",
            }
        return {"status": "unsupported", "reason": f"region mapping not defined for {node['type']}"}

    def walk(self, seeds, relations: Optional[Sequence[str]] = None, max_steps: int = 6,
             budget: int = 30, rng_seed: int = 20261106, steps: Optional[int] = None) -> Dict[str, object]:
        if steps is not None:
            max_steps = steps
        if isinstance(seeds, str):
            seeds = [seeds]
        walks = []
        for index, seed in enumerate(seeds):
            start = self.get(seed)
            if not start:
                walks.append({"seed": str(seed), "status": "not_found",
                              "path": [], "visited": []})
                continue
            rng = random.Random(f"{rng_seed}:{index}:{start['id']}")
            current = start
            visited = [start]
            seen = {start["id"]}
            path = []
            stopped = "steps"
            for _ in range(max_steps):
                nbr = self.neighbors(current["id"], relations=relations,
                                     direction="both", limit=100000)
                items = nbr["items"]
                if not items:
                    stopped = "dead_end"
                    break
                pool = [item for item in items if item["node"]["id"] not in seen]
                if not pool:
                    pool = items
                pool.sort(key=lambda item: (item["node"]["type"], item["node"]["id"]))
                item = rng.choice(pool)
                other = item["node"]
                path.append({
                    "from": current["id"],
                    "relation": item["edge"]["relation"],
                    "to": other["id"],
                    "direction": item["direction"],
                    "source": item["source"],
                })
                current = other
                if other["id"] not in seen:
                    seen.add(other["id"])
                    visited.append(other)
                if len(visited) >= budget:
                    stopped = "budget"
                    break
            walks.append({
                "seed": start["id"],
                "seed_query": str(seed),
                "status": "ok",
                "path": path,
                "visited": visited,
                "visited_by_type": self._count_types(visited),
                "stopped": stopped,
                "max_steps": max_steps,
                "budget": budget,
                "relations": list(relations) if relations else None,
            })
        return {"walks": walks, "status": "ok", "rng_seed": rng_seed}

    def recall(self, query: str, view: str = "exact", limit: int = 20) -> Dict[str, object]:
        if view == "exact":
            result = self.resolve({"text": query})
            result["candidates"] = result.get("candidates", [])[:limit]
            return result
        if view == "text":
            rows = self.conn.execute(
                "SELECT * FROM node WHERE label LIKE ? AND kind IN "
                "('gene','transcript','translation','exon') LIMIT ?",
                (f"%{query}%", limit),
            ).fetchall()
            return {"candidates": [self._row_to_node(r) for r in rows],
                    "status": "ok", "truncated": False}
        return {"candidates": [], "status": "unsupported",
                "reason": f"view {view!r} is not implemented in this M1 slice"}

    def capabilities(self) -> Dict[str, object]:
        if self._capabilities is not None:
            return self._capabilities
        base = os.path.dirname(os.path.dirname(os.path.abspath(self.db_path)))
        candidates = [
            os.path.join(base, "capabilities.json"),
            os.path.join(os.path.dirname(os.path.abspath(self.db_path)), "capabilities.json"),
            os.path.join(os.getcwd(), "capabilities.json"),
        ]
        for path in candidates:
            if os.path.exists(path):
                with open(path, encoding="utf-8") as handle:
                    self._capabilities = json.load(handle)
                return self._capabilities
        self._capabilities = {"database": "Ensembl", "release": self.meta().get("release"),
                              "relations": [], "fields": [], "region_types": []}
        return self._capabilities

    # ------------------------------------------------------------------
    # Gene model helpers
    # ------------------------------------------------------------------
    def _gene_id(self, value: str) -> str:
        value = str(value).strip()
        if value.startswith("Ensembl:"):
            value = value.split(":", 1)[1]
        value = value.split(".")[0]
        if value.startswith("ENSG"):
            return value
        node = self.get(value)
        if node and node["kind"] == "gene":
            return node["stable_id"]
        return value

    def gene_transcripts(self, gene_id: str) -> List[Dict[str, object]]:
        gene_id = self._gene_id(gene_id)
        rows = self.conn.execute(
            "SELECT n.* FROM edge e JOIN node n ON n.node_id=e.dst "
            "WHERE e.src=(SELECT node_id FROM node WHERE key=?) AND e.relation='has_transcript' "
            "ORDER BY n.key",
            (f"gene:{gene_id}",),
        ).fetchall()
        return [self._row_to_node(r) for r in rows]

    def canonical_transcript(self, gene_id: str) -> Optional[Dict[str, object]]:
        for tx in self.gene_transcripts(gene_id):
            if "Ensembl_canonical" in (tx["attributes"].get("tags") or []):
                return tx
        return None

    def transcript_exons(self, tx_id: str) -> List[Dict[str, object]]:
        tx_id = self._stable_id(tx_id)
        rows = self.conn.execute(
            "SELECT ei.rank, ei.phase, ei.end_phase, n.* "
            "FROM exon_instance ei JOIN node n ON n.key='exon:' || ei.exon_id "
            "WHERE ei.transcript_id=? ORDER BY ei.rank",
            (tx_id,),
        ).fetchall()
        out = []
        for row in rows:
            node = self._row_to_node(row)
            node["rank"] = row["rank"]
            node["exon_phase"] = row["phase"]
            node["exon_end_phase"] = row["end_phase"]
            out.append(node)
        return out

    def transcript_translation(self, tx_id: str) -> Optional[Dict[str, object]]:
        tx_id = self._stable_id(tx_id)
        row = self.conn.execute(
            "SELECT n.* FROM edge e JOIN node n ON n.node_id=e.dst "
            "WHERE e.src=(SELECT node_id FROM node WHERE key=?) AND e.relation='has_translation'",
            (f"transcript:{tx_id}",),
        ).fetchone()
        return self._row_to_node(row) if row else None

    def protein_seq(self, protein_id: str) -> Optional[str]:
        protein_id = self._stable_id(protein_id)
        row = self.conn.execute("SELECT seq FROM seq WHERE seq_id=?", (protein_id,)).fetchone()
        return row["seq"] if row else None

    def cdna_seq(self, tx_id: str) -> str:
        tx_id = self._stable_id(tx_id)
        exons = self.transcript_exons(tx_id)
        tx = self.get(f"transcript:{tx_id}")
        if not tx or not exons:
            return ""
        return splice(self.genome()[tx["seq_region"]], exons, tx["strand"])

    def transcript_stats(self, tx_id: str) -> Optional[Dict[str, object]]:
        tx_id = self._stable_id(tx_id)
        row = self.conn.execute("SELECT * FROM transcript_stats WHERE transcript_id=?", (tx_id,)).fetchone()
        return dict(row) if row else None

    # ------------------------------------------------------------------
    # Coordinate mapping
    # ------------------------------------------------------------------
    def _cds_segments(self, tx_id: str) -> List[sqlite3.Row]:
        return list(self.conn.execute(
            "SELECT * FROM cds_segment WHERE transcript_id=? ORDER BY rank", (tx_id,)
        ).fetchall())

    def map_protein_pos(self, tx_id: str, pos: int) -> Dict[str, object]:
        tx_id = self._stable_id(tx_id)
        tx = self.get(f"transcript:{tx_id}")
        if not tx:
            return {"status": "unknown_transcript", "transcript_id": tx_id, "protein_pos": pos}
        protein = self.transcript_translation(tx_id)
        if not protein:
            return {"status": "no_translation", "transcript_id": tx_id, "protein_pos": pos}
        length = protein["length"]
        if pos < 1 or length is None or pos > length:
            return {
                "status": "beyond_protein",
                "transcript_id": tx_id,
                "protein_id": protein["stable_id"],
                "protein_pos": pos,
                "protein_length": length,
            }

        segments = self._cds_segments(tx_id)
        if not segments:
            return {"status": "no_cds", "transcript_id": tx_id, "protein_pos": pos}
        phase0 = segments[0]["phase"]
        cds_len = sum(s["end"] - s["start"] + 1 for s in segments) - phase0
        partial_codon = False
        if phase0 and pos == 1:
            seg = segments[0]
            genomic = []
            for i in range(phase0):
                genomic.append(seg["start"] + i if tx["strand"] > 0 else seg["end"] - i)
            offset = 0
            partial_codon = True
        else:
            offset = (3 * (pos - 2)) if phase0 else (3 * (pos - 1))
            if offset + 3 > cds_len:
                return {
                    "status": "beyond_cds",
                    "transcript_id": tx_id,
                    "protein_id": protein["stable_id"],
                    "protein_pos": pos,
                    "protein_length": length,
                    "cds_length": cds_len,
                }
            genomic = [self._cds_offset_to_genomic(segments, phase0, o)
                       for o in (offset, offset + 1, offset + 2)]

        bases = []
        for g in genomic:
            base = self.genome()[tx["seq_region"]][g - 1].upper()
            if tx["strand"] < 0:
                base = reverse_complement(base)
            bases.append(base)
        codon = "".join(bases)
        residue = "X" if partial_codon else translate(codon)

        exons = self.transcript_exons(tx_id)
        ranks = []
        exon_ids = []
        for g in genomic:
            for exon in exons:
                if exon["start"] <= g <= exon["end"]:
                    ranks.append(exon["rank"])
                    exon_ids.append(exon["id"])
                    break
        unique_ranks = sorted(set(ranks))
        unique_exons = sorted(set(exon_ids))

        seq_pos = pos - 1
        protein_seq = self.protein_seq(protein["stable_id"]) or ""
        expected = protein_seq[seq_pos] if seq_pos < len(protein_seq) else None

        return {
            "status": "ok",
            "transcript_id": tx_id,
            "protein_id": protein["stable_id"],
            "protein_pos": pos,
            "partial_codon": partial_codon,
            "residue": residue,
            "expected_residue": expected,
            "residue_matches": (residue == expected) if expected else None,
            "codon": codon,
            "cds_offset": offset,
            "cds_offset_1based": offset + 1,
            "genomic_positions": genomic,
            "genomic_interval": [min(genomic), max(genomic)],
            "strand": tx["strand"],
            "exon_ranks": unique_ranks,
            "exon_ids": unique_exons,
            "split_codon": len(unique_exons) > 1,
        }

    def _cds_offset_to_genomic(self, segments: Sequence[sqlite3.Row], phase0: int, offset: int):
        acc = 0
        for i, seg in enumerate(segments):
            start = seg["start"]
            end = seg["end"]
            trim = phase0 if i == 0 else 0
            avail = (end - start + 1) - trim
            if offset < acc + avail:
                within = offset - acc
                if seg["strand"] > 0:
                    return start + trim + within
                return end - trim - within
            acc += avail
        raise IndexError(f"CDS offset {offset} outside segments")

    def map_interval(self, tx_id: str, start: int, end: int) -> List[Dict[str, object]]:
        return [self.map_protein_pos(tx_id, pos) for pos in range(start, end + 1)]

    def uniprot_transcripts(self, acc: str) -> List[Dict[str, object]]:
        xref = self.get(f"UniProtKB:{acc}")
        if not xref:
            return []
        items = self.neighbors(xref["id"], relations=["uniprot_xref"],
                               direction="in", limit=1000)["items"]
        out = []
        for item in items:
            protein = item["node"]
            out.append({"protein": protein,
                        "transcript": self._protein_transcript(protein["stable_id"])})
        return out

    def coverage(self, acc: str, start: int, end: int) -> Dict[str, object]:
        entries = []
        for item in self.uniprot_transcripts(acc):
            tx = item["transcript"]
            if not tx:
                continue
            entries.append(self.transcript_coverage(tx["stable_id"], start, end))
        return {"uniprot": acc, "start": start, "end": end, "entries": entries}

    def transcript_coverage(self, tx_id: str, start: int, end: int) -> Dict[str, object]:
        tx_id = self._stable_id(tx_id)
        tx = self.get(f"transcript:{tx_id}")
        protein = self.transcript_translation(tx_id)
        if not tx:
            return {"transcript_id": tx_id, "status": "unknown_transcript", "positions": []}
        if not protein:
            return {
                "transcript_id": tx_id, "transcript_version": tx["version"],
                "protein_id": None, "protein_length": None,
                "is_canonical": "Ensembl_canonical" in (tx["attributes"].get("tags") or []),
                "biotype": tx["biotype"], "status": "no_translation", "positions": [],
            }
        length = protein["length"]
        if start > length:
            status = "beyond_protein"
            positions = []
        else:
            positions = self.map_interval(tx_id, start, min(end, length))
            bad = [p for p in positions if p["status"] != "ok"]
            if bad:
                status = "partial_model"
            elif end > length:
                status = "partial_beyond"
            else:
                status = "covered"
        return {
            "transcript_id": tx_id,
            "transcript_version": tx["version"],
            "protein_id": protein["stable_id"],
            "protein_length": length,
            "is_canonical": "Ensembl_canonical" in (tx["attributes"].get("tags") or []),
            "biotype": tx["biotype"],
            "status": status,
            "positions": positions,
        }

    def gene_coverage(self, gene_id: str, start: int, end: int) -> Dict[str, object]:
        gene_id = self._gene_id(gene_id)
        entries = [self.transcript_coverage(tx["stable_id"], start, end)
                   for tx in self.gene_transcripts(gene_id)]
        return {"gene_id": gene_id, "start": start, "end": end, "entries": entries}

    # ------------------------------------------------------------------
    # Graph / region / stats
    # ------------------------------------------------------------------
    def neighborhood(self, seed: str, depth: int = 2, budget: int = 100) -> Dict[str, object]:
        start = self.get(seed)
        if not start:
            return {"seed": seed, "status": "not_found", "nodes": []}
        frontier = [start]
        seen = {start["id"]}
        nodes = [start]
        edges = []
        stopped = "depth"
        for _ in range(depth):
            next_frontier = []
            for node in frontier:
                for item in self.neighbors(node["id"], limit=1000)["items"]:
                    other = item["node"]
                    edges.append({"from": node["id"], "relation": item["edge"]["relation"],
                                  "to": other["id"]})
                    if other["id"] in seen:
                        continue
                    seen.add(other["id"])
                    nodes.append(other)
                    next_frontier.append(other)
                    if len(nodes) >= budget:
                        stopped = "budget"
                        break
                if stopped == "budget":
                    break
            if stopped == "budget":
                break
            frontier = sorted(next_frontier, key=lambda n: (n["type"], n["id"]))
            if not frontier:
                stopped = "exhausted"
                break
        return {"seed": seed, "status": "ok", "depth": depth, "budget": budget,
                "stopped": stopped, "nodes": nodes, "edges": edges}

    def region_genes(self, chrom: str, start: int, end: int) -> List[Dict[str, object]]:
        rows = self.conn.execute(
            "SELECT * FROM node WHERE kind='gene' AND seq_region=? AND start<=? AND end>=? "
            "ORDER BY start, key",
            (chrom, end, start),
        ).fetchall()
        return [self._row_to_node(r) for r in rows]

    @staticmethod
    def _count_types(nodes: Sequence[Dict[str, object]]) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        for node in nodes:
            counts[node["type"]] = counts.get(node["type"], 0) + 1
        return dict(sorted(counts.items()))

    def stats(self) -> Dict[str, object]:
        meta = self.meta()
        nodes = dict(self.conn.execute("SELECT kind, COUNT(*) FROM node GROUP BY kind ORDER BY kind"))
        node_types = dict(self.conn.execute(
            "SELECT node_type, COUNT(*) FROM node GROUP BY node_type ORDER BY node_type"))
        relations = dict(self.conn.execute(
            "SELECT relation, COUNT(*) FROM edge GROUP BY relation ORDER BY relation"))
        return {"meta": meta, "nodes": nodes, "node_types": node_types, "relations": relations}
