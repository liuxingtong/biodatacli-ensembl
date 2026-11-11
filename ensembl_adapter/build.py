"""Build a deterministic, M1-contract SQLite graph for one Ensembl chromosome.

The graph is intentionally close to the contract described in the
SciDataGym M1 section:

  * nodes and edges are the graph
  * every edge has a typed relation, direction (implicit in src/dst), attributes
  * every node and edge carries a source (database, record id, version, field,
    locator) so a neighbour can be traced back to the raw GFF3/FASTA/BioMart row
  * external references (UniProtKB, PDB, HGNC, ...) are ordinary nodes; the
    Ensembl graph records why it believes the link exists

Extra tables (cds_segment, exon_instance, seq, transcript_stats) are the
substrate for ``get_records``; they are not part of the graph contract.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sqlite3
import time
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

from .fasta import load_genome, load_peptides
from .gff import GffData, load_gff3
from .provenance import collect as collect_provenance, stamp as provenance_stamp
from .sequence import reconstruct_cds, reconstruct_cdna, translate

SCHEMA = """
PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA temp_store = MEMORY;

CREATE TABLE meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE node (
    node_id          INTEGER PRIMARY KEY,
    key              TEXT NOT NULL UNIQUE,
    global_id        TEXT NOT NULL UNIQUE,
    kind             TEXT NOT NULL,
    node_type        TEXT NOT NULL,
    stable_id        TEXT NOT NULL,
    version          INTEGER,
    label            TEXT,
    biotype          TEXT,
    seq_region       TEXT,
    start            INTEGER,
    end              INTEGER,
    strand           INTEGER,
    length           INTEGER,
    extra            TEXT,
    source_db        TEXT,
    source_record_id TEXT,
    source_version   TEXT,
    source_field     TEXT,
    source_locator   TEXT
);
CREATE INDEX node_kind   ON node(kind);
CREATE INDEX node_type   ON node(node_type);
CREATE INDEX node_stable ON node(stable_id);
CREATE INDEX node_pos    ON node(seq_region, start, end);

CREATE TABLE edge (
    edge_id          INTEGER PRIMARY KEY,
    src              INTEGER NOT NULL REFERENCES node(node_id),
    dst              INTEGER NOT NULL REFERENCES node(node_id),
    relation         TEXT NOT NULL,
    rank             INTEGER NOT NULL DEFAULT 0,
    attrs            TEXT,
    source_db        TEXT,
    source_record_id TEXT,
    source_version   TEXT,
    source_field     TEXT,
    source_locator   TEXT,
    UNIQUE(src, dst, relation, rank)
);
CREATE INDEX edge_src ON edge(src, relation);
CREATE INDEX edge_dst ON edge(dst, relation);

CREATE TABLE cds_segment (
    id            INTEGER PRIMARY KEY,
    transcript_id TEXT NOT NULL,
    protein_id    TEXT NOT NULL,
    rank          INTEGER NOT NULL,
    seq_region    TEXT NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        INTEGER NOT NULL,
    phase         INTEGER NOT NULL
);
CREATE INDEX cds_tx ON cds_segment(transcript_id);

CREATE TABLE exon_instance (
    id            INTEGER PRIMARY KEY,
    transcript_id TEXT NOT NULL,
    exon_id       TEXT NOT NULL,
    rank          INTEGER NOT NULL,
    start         INTEGER NOT NULL,
    end           INTEGER NOT NULL,
    strand        INTEGER NOT NULL,
    phase         INTEGER,
    end_phase     INTEGER
);
CREATE INDEX exon_tx ON exon_instance(transcript_id);

CREATE TABLE transcript_stats (
    transcript_id  TEXT PRIMARY KEY,
    cds_length     INTEGER,
    protein_length INTEGER,
    protein_md5    TEXT,
    cdna_length    INTEGER,
    cdna_md5       TEXT,
    validation     TEXT NOT NULL
);

CREATE TABLE seq (
    seq_id TEXT PRIMARY KEY,
    kind   TEXT NOT NULL,
    seq    TEXT NOT NULL,
    md5    TEXT NOT NULL,
    length INTEGER NOT NULL
);
"""

XREF_PREFIX = {
    "UniProt": "UniProtKB",
    "PDB": "PDB",
    "HGNC": "HGNC",
    "EntrezGene": "EntrezGene",
    "MIM": "MIM",
    "CCDS": "CCDS",
}
XREF_TYPE = {
    "UniProt": "Protein",
    "PDB": "StructureEntry",
    "HGNC": "Gene",
    "EntrezGene": "Gene",
    "MIM": "Gene",
    "CCDS": "CodingSequence",
}
XREF_RELATION = {
    "UniProt": "uniprot_xref",
    "PDB": "pdb_xref",
    "HGNC": "hgnc_xref",
    "EntrezGene": "entrez_xref",
    "MIM": "mim_xref",
    "CCDS": "ccds_xref",
}


def sha256_file(path, chunk=1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def gff_metadata(path) -> Dict[str, str]:
    from .gff import open_text

    out: Dict[str, str] = {}
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith("#!"):
                if not line.startswith("#"):
                    break
                continue
            body = line[2:].strip()
            if " " in body:
                key, value = body.split(" ", 1)
                out[key] = value
    return out


def read_tsv(path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for line_no, row in enumerate(reader, start=2):
            row["_line"] = str(line_no)
            rows.append(row)
    return rows


def load_enrichment(path: Optional[str]) -> Dict[str, Dict[str, Dict[str, object]]]:
    out: Dict[str, Dict[str, Dict[str, object]]] = {"pdb": {}, "uniprot": {}}
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        for key in ("pdb", "uniprot"):
            out[key].update(data.get(key, {}))
    return out


def load_xrefs(biomart_dir: str, enrichment: Dict[str, Dict[str, Dict[str, object]]]) -> List[Dict[str, object]]:
    records: List[Dict[str, object]] = []

    path = os.path.join(biomart_dir, "uniprot.tsv")
    if os.path.exists(path):
        for row in read_tsv(path):
            acc = (row.get("UniProtKB/Swiss-Prot ID") or "").strip()
            pid = (row.get("Protein stable ID") or "").strip()
            if not acc or not pid:
                continue
            attrs = {"external_id": acc}
            attrs.update(enrichment["uniprot"].get(acc, {}))
            records.append({
                "source": "UniProt", "acc": acc, "entity_id": pid, "entity_kind": "translation",
                "field": "uniprotswissprot", "file": "uniprot.tsv", "line": int(row["_line"]),
                "attrs": attrs,
            })

    path = os.path.join(biomart_dir, "pdb.tsv")
    if os.path.exists(path):
        for row in read_tsv(path):
            acc = (row.get("PDB ID") or "").strip()
            pid = (row.get("Protein stable ID") or "").strip()
            if not acc or not pid:
                continue
            attrs = {"external_id": acc}
            attrs.update(enrichment["pdb"].get(acc, {}))
            records.append({
                "source": "PDB", "acc": acc, "entity_id": pid, "entity_kind": "translation",
                "field": "pdb", "file": "pdb.tsv", "line": int(row["_line"]),
                "attrs": attrs,
            })

    path = os.path.join(biomart_dir, "gene_xref.tsv")
    if os.path.exists(path):
        for row in read_tsv(path):
            gid = (row.get("Gene stable ID") or "").strip()
            if not gid:
                continue
            hgnc = (row.get("HGNC ID") or "").strip()
            entrez = (row.get("NCBI gene (formerly Entrezgene) ID") or "").strip()
            if hgnc:
                records.append({
                    "source": "HGNC", "acc": hgnc, "entity_id": gid, "entity_kind": "gene",
                    "field": "hgnc_id", "file": "gene_xref.tsv", "line": int(row["_line"]),
                    "attrs": {"external_id": hgnc},
                })
            if entrez:
                records.append({
                    "source": "EntrezGene", "acc": entrez, "entity_id": gid, "entity_kind": "gene",
                    "field": "entrezgene_id", "file": "gene_xref.tsv", "line": int(row["_line"]),
                    "attrs": {"external_id": entrez},
                })

    path = os.path.join(biomart_dir, "ccds_mim.tsv")
    if os.path.exists(path):
        for row in read_tsv(path):
            gid = (row.get("Gene stable ID") or "").strip()
            tid = (row.get("Transcript stable ID") or "").strip()
            ccds = (row.get("CCDS ID") or "").strip()
            mim = (row.get("MIM gene accession") or "").strip()
            if ccds and tid:
                version = ccds.rsplit(".", 1)[1] if "." in ccds else None
                records.append({
                    "source": "CCDS", "acc": ccds, "entity_id": tid, "entity_kind": "transcript",
                    "field": "ccds", "file": "ccds_mim.tsv", "line": int(row["_line"]),
                    "attrs": {"external_id": ccds, "ccds_version": version},
                })
            if mim and gid:
                records.append({
                    "source": "MIM", "acc": mim, "entity_id": gid, "entity_kind": "gene",
                    "field": "mim_gene_accession", "file": "ccds_mim.tsv", "line": int(row["_line"]),
                    "attrs": {"external_id": mim},
                })

    records.sort(key=lambda r: (r["source"], r["entity_id"], r["acc"], r["line"]))
    return records


def reconstruct_transcript(genome: str, gff: GffData, tid: str, pep):
    tx = gff.transcripts[tid]
    exons = gff.tx_exons.get(tid, [])
    cds = gff.tx_cds.get(tid, [])

    cdna = reconstruct_cdna(genome, exons, tx.strand) if exons else ""
    cds_seq = ""
    phase = 0
    reconstructed = ""
    if cds:
        cds_seq, phase = reconstruct_cds(genome, cds, tx.strand)
        translated = translate(cds_seq)
        if translated.endswith("*"):
            translated = translated[:-1]
        if phase:
            translated = "X" + translated
        reconstructed = translated

    validation = "non_coding"
    if cds:
        if pep is None:
            validation = "no_peptide_fasta"
        elif reconstructed == pep.seq:
            validation = "exact"
        elif pep.seq.startswith("M") and reconstructed[1:] == pep.seq[1:]:
            validation = "alt_start"
        elif reconstructed.replace("*", "X") == pep.seq.replace("*", "X"):
            validation = "stop_as_X"
        else:
            validation = "mismatch"

    return {
        "tx": tx,
        "cdna": cdna,
        "cds": cds_seq,
        "phase": phase,
        "reconstructed": reconstructed,
        "pep": pep,
        "validation": validation,
    }


def _provenance_meta(args, gff_base: str, pep_base: str, release: str) -> Dict[str, str]:
    """Record when each input was retrieved and give this build an identity."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    entries = [
        (gff_base, args.gff3),
        (os.path.basename(args.genome), args.genome),
        (pep_base, args.pep),
    ]
    if args.enrichment:
        entries.append((os.path.basename(args.enrichment), args.enrichment))
    if os.path.isdir(args.biomart):
        for name in sorted(os.listdir(args.biomart)):
            if name.endswith(".tsv"):
                entries.append((f"biomart/{name}", os.path.join(args.biomart, name)))
    sources = collect_provenance(entries, os.path.join(root, "data", "manifest.json"), root)

    digest = hashlib.sha256()
    for key in sorted(sources):
        digest.update(f"{key}:{sources[key].get('sha256') or sources[key]['path']}\n".encode())
    digest.update(f"release={release}\n".encode())
    return {
        "build_id": digest.hexdigest()[:16],
        **provenance_stamp(),
        "source_retrieved_at": json.dumps(sources, sort_keys=True, ensure_ascii=False),
    }


def build(args) -> Dict[str, object]:
    started = time.time()
    gff = load_gff3(args.gff3)
    genome = load_genome(args.genome)
    peptides = load_peptides(args.pep)
    enrichment = load_enrichment(args.enrichment)
    xref_records = load_xrefs(args.biomart, enrichment)

    chrom = args.chromosome
    if chrom not in genome:
        raise SystemExit(f"chromosome {chrom!r} not found in {args.genome}: {list(genome)}")
    genome_seq = genome[chrom]

    gff_base = os.path.basename(args.gff3)
    pep_base = os.path.basename(args.pep)
    release = str(args.release)

    genes = {gid: g for gid, g in gff.genes.items() if g.seq_region == chrom}
    transcripts = {tid: t for tid, t in gff.transcripts.items()
                   if t.seq_region == chrom and t.gene_id in genes}
    exons = {eid: e for eid, e in gff.exons.items() if e.seq_region == chrom}

    tx_stats: Dict[str, Dict[str, object]] = {}
    for tid in transcripts:
        cds = gff.tx_cds.get(tid, [])
        pep = peptides.get(cds[0].protein_id) if cds else None
        rec = reconstruct_transcript(genome_seq, gff, tid, pep)
        tx_stats[tid] = {
            "cds_length": len(rec["cds"]),
            "protein_length": len(pep.seq) if pep else (len(rec["reconstructed"]) if cds else None),
            "protein_md5": hashlib.md5(pep.seq.encode()).hexdigest() if pep else None,
            "cdna_length": len(rec["cdna"]),
            "cdna_md5": hashlib.md5(rec["cdna"].encode()).hexdigest() if rec["cdna"] else None,
            "validation": rec["validation"],
        }

    pid_to_txs: Dict[str, List[str]] = defaultdict(list)
    for tid in transcripts:
        cds = gff.tx_cds.get(tid, [])
        if cds:
            pid_to_txs[cds[0].protein_id].append(tid)
    protein_ids: Set[str] = set(pid_to_txs)

    # ------------------------------------------------------------------
    # Nodes
    # ------------------------------------------------------------------
    nodes: Dict[str, Dict[str, object]] = {}

    def add_node(key, global_id, kind, node_type, stable_id, *,
                 source_db, source_record_id, source_version, source_field, source_locator,
                 version=None, label=None, biotype=None, seq_region=None, start=None, end=None,
                 strand=None, length=None, extra=None):
        if key in nodes:
            return
        nodes[key] = {
            "key": key, "global_id": global_id, "kind": kind, "node_type": node_type,
            "stable_id": stable_id, "version": version, "label": label, "biotype": biotype,
            "seq_region": seq_region, "start": start, "end": end, "strand": strand,
            "length": length, "extra": extra or {},
            "source_db": source_db, "source_record_id": source_record_id,
            "source_version": source_version, "source_field": source_field,
            "source_locator": source_locator,
        }

    assembly_key = f"assembly:{args.assembly}"
    region_key = f"region:{chrom}"
    add_node(
        assembly_key, f"Ensembl:{args.assembly}", "assembly", "Assembly", args.assembly,
        source_db="Ensembl", source_record_id=args.assembly, source_version=release,
        source_field="GFF3 header", source_locator=f"{gff_base}:1",
        label=args.assembly,
    )
    add_node(
        region_key, f"Ensembl:{chrom}", "region", "Region", chrom,
        source_db="Ensembl", source_record_id=chrom, source_version=release,
        source_field="GFF3 header", source_locator=f"{gff_base}:1",
        label=chrom, seq_region=chrom, start=1, end=len(genome_seq), length=len(genome_seq),
    )

    for gid, gene in genes.items():
        add_node(
            f"gene:{gid}", f"Ensembl:{gid}", "gene", "Gene", gid,
            source_db="Ensembl", source_record_id=gid, source_version=release,
            source_field=f"GFF3 {gene.feature_type}", source_locator=f"{gff_base}:{gene.line}",
            version=gene.version, label=gene.name, biotype=gene.biotype,
            seq_region=gene.seq_region, start=gene.start, end=gene.end,
            strand=gene.strand, length=gene.end - gene.start + 1,
            extra={"description": gene.description, "source": gene.source,
                   "feature_type": gene.feature_type},
        )

    for tid, tx in transcripts.items():
        add_node(
            f"transcript:{tid}", f"Ensembl:{tid}", "transcript", "Transcript", tid,
            source_db="Ensembl", source_record_id=tid, source_version=release,
            source_field=f"GFF3 {tx.feature_type}", source_locator=f"{gff_base}:{tx.line}",
            version=tx.version, label=tx.name, biotype=tx.biotype,
            seq_region=tx.seq_region, start=tx.start, end=tx.end,
            strand=tx.strand, length=tx.end - tx.start + 1,
            extra={"tags": tx.tags, "tsl": tx.tsl, "ccds": tx.ccds,
                   "feature_type": tx.feature_type, "gene_id": tx.gene_id},
        )

    for eid, exon in exons.items():
        add_node(
            f"exon:{eid}", f"Ensembl:{eid}", "exon", "Exon", eid,
            source_db="Ensembl", source_record_id=eid, source_version=release,
            source_field="GFF3 exon", source_locator=f"{gff_base}:{exon.line}",
            version=exon.version, label=eid, seq_region=exon.seq_region,
            start=exon.start, end=exon.end, strand=exon.strand,
            length=exon.end - exon.start + 1,
            extra={"rank": exon.rank, "phase": exon.phase,
                   "end_phase": exon.end_phase, "constitutive": exon.constitutive},
        )

    for pid in sorted(protein_ids):
        pep = peptides.get(pid)
        tx_ids = pid_to_txs[pid]
        tx0 = gff.transcripts[tx_ids[0]] if tx_ids else None
        if pep is not None:
            source_field = "FASTA pep"
            locator = f"{pep_base}:header {pep.header_line}"
            source_version = str(pep.version) if pep.version is not None else release
        else:
            cds0 = gff.tx_cds[tx_ids[0]][0]
            source_field = "GFF3 CDS"
            locator = f"{gff_base}:{cds0.line}"
            source_version = release
        add_node(
            f"translation:{pid}", f"Ensembl:{pid}", "translation", "Protein", pid,
            source_db="Ensembl", source_record_id=pid, source_version=source_version,
            source_field=source_field, source_locator=locator,
            version=pep.version if pep else (gff.tx_cds[tx_ids[0]][0].protein_version if tx_ids else None),
            label=(pep.gene_symbol if pep and pep.gene_symbol else (tx0.name if tx0 else pid)),
            biotype=(pep.biotype if pep else (tx0.biotype if tx0 else "")),
            length=len(pep.seq) if pep else None,
            extra={"gene_id": tx0.gene_id if tx0 else "", "transcript_id": tx0.transcript_id if tx0 else ""},
        )

    # ------------------------------------------------------------------
    # Cross-reference nodes
    # ------------------------------------------------------------------
    def xref_key(source: str, acc: str) -> str:
        return f"xref:{source}:{acc}"

    def add_xref(record: Dict[str, object]) -> Optional[str]:
        source = record["source"]
        acc = record["acc"]
        entity_id = record["entity_id"]
        key = xref_key(source, acc)
        global_id = f"{XREF_PREFIX[source]}:{acc}"
        locator = f"biomart/{record['file']}:row {record['line']}"
        add_node(
            key, global_id, "xref", XREF_TYPE[source], acc,
            source_db="Ensembl", source_record_id=entity_id, source_version=release,
            source_field=f"BioMart {record['field']}", source_locator=locator,
            label=acc,
            extra={"source": source, "external_id": acc, "entity_id": entity_id},
        )
        return key

    entity_prefix = {"gene": "gene:", "transcript": "transcript:", "translation": "translation:"}
    entity_sets = {"gene": genes, "transcript": transcripts, "translation": protein_ids}

    # ------------------------------------------------------------------
    # Edges
    # ------------------------------------------------------------------
    edges: Dict[Tuple[str, str, str, int], Dict[str, object]] = {}

    def add_edge(src, dst, relation, rank=0, attrs=None, *,
                 source_db, source_record_id, source_version, source_field, source_locator):
        key = (src, dst, relation, rank)
        if key in edges:
            return
        edges[key] = {
            "src": src, "dst": dst, "relation": relation, "rank": rank,
            "attrs": attrs or {},
            "source_db": source_db, "source_record_id": source_record_id,
            "source_version": source_version, "source_field": source_field,
            "source_locator": source_locator,
        }

    add_edge(
        assembly_key, region_key, "has_region", 0, {"chromosome": chrom},
        source_db="Ensembl", source_record_id=args.assembly, source_version=release,
        source_field="GFF3 header", source_locator=f"{gff_base}:1",
    )

    for gid, gene in genes.items():
        add_edge(
            region_key, f"gene:{gid}", "contains", 0, {"gene_version": gene.version},
            source_db="Ensembl", source_record_id=gid, source_version=release,
            source_field=f"GFF3 {gene.feature_type}", source_locator=f"{gff_base}:{gene.line}",
        )

    for tid, tx in transcripts.items():
        add_edge(
            f"gene:{tx.gene_id}", f"transcript:{tid}", "has_transcript", 0,
            {"transcript_version": tx.version, "biotype": tx.biotype},
            source_db="Ensembl", source_record_id=tid, source_version=release,
            source_field=f"GFF3 {tx.feature_type}", source_locator=f"{gff_base}:{tx.line}",
        )
        for rank, exon in enumerate(gff.tx_exons.get(tid, []), start=1):
            add_edge(
                f"transcript:{tid}", f"exon:{exon.exon_id}", "has_exon", rank,
                {"rank": exon.rank or rank, "exon_version": exon.version},
                source_db="Ensembl", source_record_id=exon.exon_id, source_version=release,
                source_field="GFF3 exon", source_locator=f"{gff_base}:{exon.line}",
            )
        cds = gff.tx_cds.get(tid, [])
        if cds:
            pid = cds[0].protein_id
            pep = peptides.get(pid)
            locator = f"{gff_base}:{cds[0].line}" + (f"; {pep_base}:header {pep.header_line}" if pep else "")
            add_edge(
                f"transcript:{tid}", f"translation:{pid}", "has_translation", 0,
                {"protein_version": pep.version if pep else cds[0].protein_version,
                 "cds_length": tx_stats[tid]["cds_length"]},
                source_db="Ensembl", source_record_id=pid, source_version=release,
                source_field="GFF3 CDS + FASTA pep" if pep else "GFF3 CDS",
                source_locator=locator,
            )

    # Cross-reference edges
    for record in xref_records:
        eid = record["entity_id"]
        kind = record["entity_kind"]
        if eid not in entity_sets[kind]:
            continue
        dst = add_xref(record)
        if dst is None:
            continue
        src = entity_prefix[kind] + eid
        attrs = dict(record["attrs"])
        attrs.setdefault("external_id", record["acc"])
        add_edge(
            src, dst, XREF_RELATION[record["source"]], 0, attrs,
            source_db="Ensembl", source_record_id=eid, source_version=release,
            source_field=f"BioMart {record['field']}",
            source_locator=f"biomart/{record['file']}:row {record['line']}",
        )

    # ------------------------------------------------------------------
    # Write database
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.db)), exist_ok=True)
    if os.path.exists(args.db):
        os.remove(args.db)
    conn = sqlite3.connect(args.db)
    conn.executescript(SCHEMA)

    meta = {
        "format": "ensembl-graph-m1",
        "schema_version": "m1-contract-v1",
        "release": release,
        "assembly": args.assembly,
        "chromosome": chrom,
        "gene_build_last_updated": gff_metadata(args.gff3).get("genebuild-last-updated", ""),
        "genome_build_accession": gff_metadata(args.gff3).get("genome-build-accession", ""),
        "genome_length": str(len(genome_seq)),
        "source_gff3": gff_base,
        "source_genome": os.path.basename(args.genome),
        "source_pep": pep_base,
        "source_biomart": os.path.basename(args.biomart.rstrip("/")),
        "capabilities": os.path.basename(args.capabilities) if args.capabilities else "capabilities.json",
    }
    meta.update(_provenance_meta(args, gff_base, pep_base, release))
    for key, value in meta.items():
        conn.execute("INSERT INTO meta(key,value) VALUES (?,?)", (key, value))

    node_keys = sorted(nodes)
    key_to_id = {key: i for i, key in enumerate(node_keys, start=1)}
    conn.executemany(
        "INSERT INTO node(node_id,key,global_id,kind,node_type,stable_id,version,label,biotype,"
        "seq_region,start,end,strand,length,extra,source_db,source_record_id,source_version,"
        "source_field,source_locator) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
        [
            (
                key_to_id[key], key, n["global_id"], n["kind"], n["node_type"], n["stable_id"],
                n["version"], n["label"], n["biotype"], n["seq_region"], n["start"], n["end"],
                n["strand"], n["length"], json.dumps(n["extra"], sort_keys=True, ensure_ascii=False),
                n["source_db"], n["source_record_id"], n["source_version"],
                n["source_field"], n["source_locator"],
            )
            for key, n in ((k, nodes[k]) for k in node_keys)
        ],
    )

    ordered_edges = sorted(edges)
    conn.executemany(
        "INSERT INTO edge(src,dst,relation,rank,attrs,source_db,source_record_id,source_version,"
        "source_field,source_locator) VALUES (?,?,?,?,?,?,?,?,?,?)",
        [
            (
                key_to_id[e["src"]], key_to_id[e["dst"]], e["relation"], e["rank"],
                json.dumps(e["attrs"], sort_keys=True, ensure_ascii=False) if e["attrs"] else None,
                e["source_db"], e["source_record_id"], e["source_version"],
                e["source_field"], e["source_locator"],
            )
            for e in (edges[k] for k in ordered_edges)
        ],
    )

    cds_rows = []
    exon_rows = []
    for tid in sorted(transcripts):
        for seg in gff.tx_cds.get(tid, []):
            cds_rows.append((tid, seg.protein_id, seg.rank, seg.seq_region,
                             seg.start, seg.end, seg.strand, seg.phase))
        for rank, exon in enumerate(gff.tx_exons.get(tid, []), start=1):
            exon_rows.append((tid, exon.exon_id, rank, exon.start, exon.end, exon.strand,
                              exon.phase, exon.end_phase))
    conn.executemany(
        "INSERT INTO cds_segment(transcript_id,protein_id,rank,seq_region,start,end,strand,phase)"
        " VALUES (?,?,?,?,?,?,?,?)",
        cds_rows,
    )
    conn.executemany(
        "INSERT INTO exon_instance(transcript_id,exon_id,rank,start,end,strand,phase,end_phase)"
        " VALUES (?,?,?,?,?,?,?,?)",
        exon_rows,
    )

    seq_rows = []
    for pid in sorted(protein_ids):
        pep = peptides.get(pid)
        if pep is None:
            continue
        seq_rows.append((pid, "protein", pep.seq, hashlib.md5(pep.seq.encode()).hexdigest(), len(pep.seq)))
    conn.executemany("INSERT INTO seq(seq_id,kind,seq,md5,length) VALUES (?,?,?,?,?)", seq_rows)

    conn.executemany(
        "INSERT INTO transcript_stats(transcript_id,cds_length,protein_length,protein_md5,"
        "cdna_length,cdna_md5,validation) VALUES (?,?,?,?,?,?,?)",
        [
            (tid, s["cds_length"], s["protein_length"], s["protein_md5"],
             s["cdna_length"], s["cdna_md5"], s["validation"])
            for tid, s in sorted(tx_stats.items())
        ],
    )
    conn.commit()

    # ------------------------------------------------------------------
    # Stats / fingerprint
    # ------------------------------------------------------------------
    node_by_kind = dict(conn.execute("SELECT kind, COUNT(*) FROM node GROUP BY kind ORDER BY kind"))
    node_by_type = dict(conn.execute("SELECT node_type, COUNT(*) FROM node GROUP BY node_type ORDER BY node_type"))
    edge_by_relation = dict(conn.execute(
        "SELECT relation, COUNT(*) FROM edge GROUP BY relation ORDER BY relation"))
    validation_counts = dict(conn.execute(
        "SELECT validation, COUNT(*) FROM transcript_stats GROUP BY validation ORDER BY validation"))
    xref_by_source = dict(conn.execute(
        "SELECT json_extract(extra,'$.source'), COUNT(*) FROM node WHERE kind='xref' GROUP BY 1 ORDER BY 1"))
    gene_biotypes = dict(conn.execute(
        "SELECT biotype, COUNT(*) FROM node WHERE kind='gene' GROUP BY biotype ORDER BY 2 DESC, 1"))

    counts = {
        "genes": len(genes),
        "transcripts": len(transcripts),
        "exons_unique": len(exons),
        "exon_instances": len(exon_rows),
        "translations": len(protein_ids),
        "nodes": len(node_keys),
        "edges": len(ordered_edges),
        "nodes_by_kind": node_by_kind,
        "nodes_by_type": node_by_type,
        "edges_by_relation": edge_by_relation,
        "xref_by_source": xref_by_source,
        "gene_biotypes": gene_biotypes,
        "validation": validation_counts,
    }

    fingerprint = db_fingerprint(conn)
    conn.close()

    result = {
        "release": release,
        "assembly": args.assembly,
        "chromosome": chrom,
        "schema_version": "m1-contract-v1",
        "counts": counts,
        "fingerprint": fingerprint,
        "inputs": {
            "gff3": {"path": args.gff3, "sha256": sha256_file(args.gff3)},
            "genome": {"path": args.genome, "sha256": sha256_file(args.genome)},
            "pep": {"path": args.pep, "sha256": sha256_file(args.pep)},
            "biomart_dir": args.biomart,
            "enrichment": args.enrichment,
        },
        "elapsed_seconds": round(time.time() - started, 2),
    }
    with open(args.stats, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, ensure_ascii=False)
    return result


def db_fingerprint(conn: sqlite3.Connection) -> str:
    """Hash logical table contents, independent of page layout."""
    digest = hashlib.sha256()
    tables = [
        "meta", "node", "edge", "cds_segment", "exon_instance", "transcript_stats", "seq",
    ]
    for table in tables:
        digest.update(("TABLE:" + table + "\n").encode())
        cursor = conn.execute(f"SELECT * FROM {table} ORDER BY rowid")
        cols = [d[0] for d in cursor.description]
        digest.update(("|".join(cols) + "\n").encode())
        for row in cursor:
            digest.update(json.dumps(tuple(row), sort_keys=True, ensure_ascii=False).encode())
            digest.update(b"\n")
    return digest.hexdigest()


def main(argv=None):
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--pep", required=True)
    parser.add_argument("--biomart", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--enrichment", default=os.path.join(base, "data", "raw", "xref_enrichment.json"))
    parser.add_argument("--capabilities", default=os.path.join(base, "capabilities.json"))
    parser.add_argument("--chromosome", default="20")
    parser.add_argument("--release", default="116")
    parser.add_argument("--assembly", default="GRCh38.p14")
    args = parser.parse_args(argv)
    result = build(args)
    print(json.dumps(result["counts"], indent=2, sort_keys=True))
    print("fingerprint", result["fingerprint"])
    print("elapsed", result["elapsed_seconds"], "s")


if __name__ == "__main__":
    main()
