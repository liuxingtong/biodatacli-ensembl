#!/usr/bin/env python3
"""Build the cross-database key alignment table (the "crosswalk").

Principle: every adapter output must be expressible in authoritative IDs, and
the mapping between those IDs is itself data with provenance.  This script
turns the frozen BioMart exports plus the SIFTS enrichment into one
deterministic edge list:

    ENSP ↔ UniProtKB        protein_uniprot
    ENSP ↔ PDB              protein_pdb
    ENSG ↔ HGNC             gene_hgnc
    ENSG ↔ EntrezGene       gene_entrez
    ENSG ↔ MIM              gene_mim
    ENST ↔ CCDS             transcript_ccds
    PDB  ↔ UniProtKB        pdb_uniprot_sifts   (SIFTS, not BioMart)

Links are many-to-many in the real data (789 UniProt accessions map to more
than one chromosome-20 protein).  Rather than silently picking the first hit,
each edge records how many alternatives share its source, so callers can see
the ambiguity.

Output:
    crosswalk/crosswalk.jsonl.gz   one edge per line, sorted
    crosswalk/stats.json           counts, ambiguity and build identity
"""

from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter.provenance import collect as collect_provenance  # noqa: E402

RAW = os.path.join(ROOT, "data", "raw")
OUT_DIR = os.path.join(ROOT, "crosswalk")
DEFAULT_OUT = os.path.join(OUT_DIR, "crosswalk.jsonl.gz")
DEFAULT_STATS = os.path.join(OUT_DIR, "stats.json")


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


BIOMART_DIR = os.path.join(RAW, "biomart")


def read_tsv(rel: str):
    """Yield ``(line_number, row)`` with the header on line 1."""
    path = os.path.join(RAW, rel)
    if rel.startswith("biomart/"):
        path = os.path.join(BIOMART_DIR, rel.split("/", 1)[1])
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for number, line in enumerate(handle, start=2):
            values = line.rstrip("\n").split("\t")
            if len(values) < len(header):
                values += [""] * (len(header) - len(values))
            yield number, dict(zip(header, values))


def load_enrichment() -> dict:
    with open(os.path.join(RAW, "xref_enrichment.json"), encoding="utf-8") as handle:
        return json.load(handle)


def cardinality(edges, source_key, target_key) -> dict:
    """How many targets per source, and sources per target."""
    forward = collections.defaultdict(set)
    backward = collections.defaultdict(set)
    for edge in edges:
        forward[edge[source_key]].add(edge[target_key])
        backward[edge[target_key]].add(edge[source_key])
    return {"forward": forward, "backward": backward}


def classify(n_out: int, n_in: int) -> str:
    if n_out == 1 and n_in == 1:
        return "1:1"
    if n_out == 1:
        return "1:n"
    if n_in == 1:
        return "n:1"
    return "n:m"


def build(manifest_path: str, biomart_dir: str = os.path.join(RAW, "biomart"),
          out_path: str = DEFAULT_OUT) -> dict:
    provenance = collect_provenance(
        [(f"biomart/{name}", os.path.join(biomart_dir, name))
         for name in ("uniprot.tsv", "pdb.tsv", "gene_xref.tsv", "ccds_mim.tsv",
                      "go.tsv", "reactome.tsv")]
        + [("xref_enrichment.json", os.path.join(RAW, "xref_enrichment.json")),
           ("kegg/pathway.txt", os.path.join(RAW, "kegg", "pathway.txt")),
           ("kegg/enzyme.txt", os.path.join(RAW, "kegg", "enzyme.txt"))],
        manifest_path, ROOT,
    )
    release = "116"
    edges = []

    def ev(key: str, row: int, field: str) -> dict:
        stamp = provenance[key]
        return {
            "file": stamp["path"], "row": row, "field": field,
            "sha256": stamp["sha256"], "release": release,
            "retrieved_at": stamp["retrieved_at"],
            "retrieved_at_status": stamp["retrieved_at_status"],
        }

    # --- raw pairs, one list per relation -----------------------------
    pairs = collections.defaultdict(list)   # relation -> [(from, to, evidence, via)]

    for number, row in read_tsv("biomart/uniprot.tsv"):
        protein, accession = row["Protein stable ID"], row["UniProtKB/Swiss-Prot ID"]
        if protein and accession:
            pairs["protein_uniprot"].append((
                ("Ensembl", protein), ("UniProtKB", accession),
                ev("biomart/uniprot.tsv", number, "UniProtKB/Swiss-Prot ID"),
                {"gene": row["Gene stable ID"], "transcript": row["Transcript stable ID"]}))

    for number, row in read_tsv("biomart/pdb.tsv"):
        protein, pdb = row["Protein stable ID"], row["PDB ID"]
        if protein and pdb:
            pairs["protein_pdb"].append((
                ("Ensembl", protein), ("PDB", pdb.upper()),
                ev("biomart/pdb.tsv", number, "PDB ID"),
                {"gene": row["Gene stable ID"], "transcript": row["Transcript stable ID"]}))

    for number, row in read_tsv("biomart/gene_xref.tsv"):
        gene = row["Gene stable ID"]
        if not gene:
            continue
        if row["HGNC ID"]:
            pairs["gene_hgnc"].append((("Ensembl", gene), ("HGNC", row["HGNC ID"]),
                                       ev("biomart/gene_xref.tsv", number, "HGNC ID"), {}))
        if row["NCBI gene (formerly Entrezgene) ID"]:
            pairs["gene_entrez"].append((("Ensembl", gene), ("EntrezGene", row["NCBI gene (formerly Entrezgene) ID"]),
                                         ev("biomart/gene_xref.tsv", number, "NCBI gene (formerly Entrezgene) ID"), {}))

    for number, row in read_tsv("biomart/ccds_mim.tsv"):
        gene, transcript = row["Gene stable ID"], row["Transcript stable ID"]
        if transcript and row["CCDS ID"]:
            pairs["transcript_ccds"].append((("Ensembl", transcript), ("CCDS", row["CCDS ID"]),
                                             ev("biomart/ccds_mim.tsv", number, "CCDS ID"),
                                             {"gene": gene}))
        if gene and row["MIM gene accession"]:
            pairs["gene_mim"].append((("Ensembl", gene), ("MIM", row["MIM gene accession"]),
                                      ev("biomart/ccds_mim.tsv", number, "MIM gene accession"), {}))

    for number, row in read_tsv("biomart/go.tsv"):
        gene, go_id = row["Gene stable ID"], row["GO term accession"]
        if gene and go_id:
            pairs["gene_go"].append((
                ("Ensembl", gene), ("GO", go_id),
                ev("biomart/go.tsv", number, "GO term accession"),
                {"go_name": row["GO term name"], "go_domain": row["GO domain"]}))

    for number, row in read_tsv("biomart/reactome.tsv"):
        gene, pathway = row["Gene stable ID"], row["Reactome ID"]
        if gene and pathway:
            pairs["gene_reactome"].append((
                ("Ensembl", gene), ("Reactome", pathway),
                ev("biomart/reactome.tsv", number, "Reactome ID"), {}))

    # --- KEGG: Ensembl gene -> pathway / EC, bridged on Entrez ID ------
    entrez_to_gene = {}
    for _number, row in read_tsv("biomart/gene_xref.tsv"):
        gene_id = row["Gene stable ID"]
        entrez = row["NCBI gene (formerly Entrezgene) ID"]
        if gene_id and entrez:
            entrez_to_gene[entrez] = gene_id

    def read_kegg(rel: str, target_prefix: str, relation: str):
        path = os.path.join(RAW, rel)
        if not os.path.exists(path):
            return
        key = rel.replace("data/raw/", "").replace("raw/", "")
        with open(path, encoding="utf-8") as handle:
            for number, line in enumerate(handle, start=1):
                parts = line.split()
                if len(parts) != 2 or not parts[0].startswith("hsa:"):
                    continue
                entrez = parts[0].split(":", 1)[1]
                gene_id = entrez_to_gene.get(entrez)
                if not gene_id or not parts[1].startswith(target_prefix):
                    continue
                target = parts[1].split(":", 1)[1]
                pairs[relation].append((
                    ("Ensembl", gene_id), (("KEGG", target) if relation == "gene_kegg_pathway"
                                           else ("EC", target)),
                    ev(key, number, "KEGG link"), {"entrez": entrez}))

    read_kegg("kegg/pathway.txt", "path:", "gene_kegg_pathway")
    read_kegg("kegg/enzyme.txt", "ec:", "gene_ec")

    # --- SIFTS: PDB entry <-> UniProt accession -----------------------
    enrichment = load_enrichment()
    cif_sha = provenance["xref_enrichment.json"]["sha256"]
    for pdb_id, record in sorted((enrichment.get("pdb") or {}).items()):
        if not record.get("uniprot"):
            continue
        pairs["pdb_uniprot_sifts"].append((
            ("PDB", pdb_id.upper()), ("UniProtKB", record["uniprot"]),
            {"file": "data/raw/xref_enrichment.json", "row": None,
             "field": "SIFTS", "sha256": cif_sha, "release": "2026-11-01",
             "retrieved_at": provenance["xref_enrichment.json"]["retrieved_at"],
             "retrieved_at_status": provenance["xref_enrichment.json"]["retrieved_at_status"],
             "note": record.get("evidence")},
            {"chains": record.get("chains"), "protein_range": record.get("protein_range")}))

    # BioMart repeats a gene once per transcript; keep the first evidence row
    # per (from, to) pair so the edge list is a set, not a join artefact.
    deduped = {}
    for relation, rows in pairs.items():
        seen = {}
        for src, dst, evidence, via in rows:
            seen.setdefault((src, dst), (evidence, via))
        deduped[relation] = [(src, dst, evidence, via)
                             for (src, dst), (evidence, via) in seen.items()]
    pairs = deduped

    # --- emit with explicit cardinality -------------------------------
    relation_stats = {}
    for relation, rows in sorted(pairs.items()):
        forward = collections.defaultdict(set)
        backward = collections.defaultdict(set)
        for src, dst, _ev, _via in rows:
            forward[src].add(dst)
            backward[dst].add(src)
        ambiguous = 0
        for src, dst, evidence, via in rows:
            n_out, n_in = len(forward[src]), len(backward[dst])
            direction_ambiguous = n_out > 1 or n_in > 1
            ambiguous += int(direction_ambiguous)
            edges.append({
                "from": {"namespace": src[0], "id": src[1]},
                "to": {"namespace": dst[0], "id": dst[1]},
                "relation": relation,
                "cardinality": classify(n_out, n_in),
                "source_fanout": n_out,
                "target_fanout": n_in,
                "evidence": evidence,
                "via": via or {},
            })
        relation_stats[relation] = {
            "edges": len(rows),
            "ambiguous_edges": ambiguous,
            "distinct_sources": len(forward),
            "distinct_targets": len(backward),
        }

    edges.sort(key=lambda e: (e["from"]["namespace"], e["from"]["id"],
                              e["relation"], e["to"]["namespace"], e["to"]["id"]))

    digest = hashlib.sha256()
    for key in sorted(provenance):
        digest.update(f"{key}:{provenance[key]['sha256']}\n".encode())
    build_id = digest.hexdigest()[:16]

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with gzip.open(out_path, "wt", encoding="utf-8", compresslevel=9) as handle:
        for edge in edges:
            handle.write(json.dumps(edge, sort_keys=True, ensure_ascii=False) + "\n")

    return {
        "build_id": build_id,
        "release": release,
        "edge_count": len(edges),
        "relations": relation_stats,
        "namespaces": sorted({e["from"]["namespace"] for e in edges} |
                             {e["to"]["namespace"] for e in edges}),
        "provenance": provenance,
        "output": os.path.relpath(out_path, ROOT),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=os.path.join(ROOT, "data", "manifest.json"))
    parser.add_argument("--stats", default=DEFAULT_STATS)
    parser.add_argument("--biomart", default=os.path.join(RAW, "biomart"))
    parser.add_argument("--out", default=DEFAULT_OUT)
    args = parser.parse_args(argv)

    global BIOMART_DIR
    BIOMART_DIR = args.biomart
    stats = build(args.manifest, args.biomart, args.out)
    with open(args.stats, "w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2, sort_keys=True, ensure_ascii=False)
    print(json.dumps({key: stats[key] for key in ("build_id", "edge_count", "namespaces")},
                     indent=2))
    for relation, counts in stats["relations"].items():
        print(f"  {relation:<20} {counts['edges']:>6} edges  "
              f"{counts['ambiguous_edges']:>5} ambiguous  "
              f"{counts['distinct_sources']:>6} sources -> {counts['distinct_targets']} targets")
    print(f"wrote {args.stats}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
