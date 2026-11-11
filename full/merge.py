#!/usr/bin/env python3
"""Merge the per-chromosome SQLite graphs into one database.

Node ids are per-part integers, so merging renumbers them and rewrites the edges
through the shared ``key`` (``gene:ENSG…``, ``xref:PDB:4DNK``, …).  That also
de-duplicates entities which appear on several chromosomes — a PDB entry
referenced by proteins from chr3 and chr20 becomes one node with two edges,
which is what a whole-genome graph should look like.

    python3 full/merge.py
    # -> full/index/ensembl_grch38_r116.sqlite + stats.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter.build import SCHEMA, db_fingerprint  # noqa: E402

PARTS = os.path.join(ROOT, "full", "parts")
OUT_DIR = os.path.join(ROOT, "full", "index")
WANTED = [str(n) for n in range(1, 23)] + ["X", "Y", "MT"]

NODE_COLS = ("key", "global_id", "kind", "node_type", "stable_id", "version", "label",
             "biotype", "seq_region", "start", "end", "strand", "length", "extra",
             "source_db", "source_record_id", "source_version", "source_field", "source_locator")
EDGE_COLS = ("src", "dst", "relation", "rank", "attrs", "source_db", "source_record_id",
             "source_version", "source_field", "source_locator")


def part_files() -> list[str]:
    found = []
    for chrom in WANTED:
        path = os.path.join(PARTS, f"chr{chrom}.sqlite")
        if os.path.exists(path):
            found.append(path)
    return found


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts", default=PARTS)
    parser.add_argument("--out", default=os.path.join(OUT_DIR, "ensembl_grch38_r116.sqlite"))
    args = parser.parse_args(argv)

    paths = part_files()
    if not paths:
        raise SystemExit("no parts found: run full/build_all.py first")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    if os.path.exists(args.out):
        os.remove(args.out)

    started = time.time()
    final = sqlite3.connect(args.out)
    final.executescript(SCHEMA)

    next_id = 1
    key_to_id: dict[str, int] = {}
    seen_as_primary = 0

    for path in paths:
        part = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        part.row_factory = sqlite3.Row
        local: dict[int, str] = {}

        rows = part.execute("SELECT node_id, " + ", ".join(NODE_COLS) + " FROM node").fetchall()
        payload = []
        for row in rows:
            key = row["key"]
            local[row["node_id"]] = key
            if key in key_to_id:
                seen_as_primary += 1
                continue
            key_to_id[key] = next_id
            payload.append((next_id,) + tuple(row[column] for column in NODE_COLS))
            next_id += 1
        final.executemany(
            "INSERT INTO node(node_id," + ",".join(NODE_COLS) + ") VALUES ("
            + ",".join("?" * (len(NODE_COLS) + 1)) + ")", payload)

        edges = []
        for row in part.execute("SELECT " + ", ".join(EDGE_COLS) + " FROM edge"):
            src = key_to_id.get(local.get(row["src"]))
            dst = key_to_id.get(local.get(row["dst"]))
            if src is None or dst is None:
                continue
            edges.append((src, dst) + tuple(row[column] for column in EDGE_COLS[2:]))
        final.executemany(
            "INSERT OR IGNORE INTO edge(" + ",".join(EDGE_COLS) + ") VALUES ("
            + ",".join("?" * len(EDGE_COLS)) + ")", edges)

        for table, columns in (
                ("cds_segment", "transcript_id,protein_id,rank,seq_region,start,end,strand,phase"),
                ("exon_instance", "transcript_id,exon_id,rank,start,end,strand,phase,end_phase"),
                ("transcript_stats", "transcript_id,cds_length,protein_length,protein_md5,"
                                     "cdna_length,cdna_md5,validation"),
                ("seq", "seq_id,kind,seq,md5,length")):
            data = part.execute(f"SELECT {columns} FROM {table}").fetchall()
            final.executemany(
                f"INSERT OR IGNORE INTO {table}({columns}) VALUES ("
                + ",".join("?" * len(columns.split(","))) + ")",
                [tuple(row) for row in data])
        part.close()
        final.commit()
        print(f"  merged {os.path.basename(path)}: nodes={len(rows):,} "
              f"(shared keys so far: {seen_as_primary:,})")

    # --- metadata -----------------------------------------------------
    first = paths[0]
    meta = dict(sqlite3.connect(f"file:{first}?mode=ro", uri=True)
                .execute("SELECT key,value FROM meta").fetchall())
    chromosomes = [os.path.basename(path)[3:-7] for path in paths]
    total_bases = 0
    for path in paths:
        base = os.path.basename(path)
        stats_path = os.path.join(os.path.dirname(path), base.replace(".sqlite", ".stats.json"))
        if os.path.exists(stats_path):
            total_bases += int(json.load(open(stats_path, encoding="utf-8"))["counts"].get("nodes", 0))
    digest = hashlib.sha256()
    for path in paths:
        stats_path = os.path.join(os.path.dirname(path),
                                  os.path.basename(path).replace(".sqlite", ".stats.json"))
        if os.path.exists(stats_path):
            digest.update(json.load(open(stats_path, encoding="utf-8"))["fingerprint"].encode())
    # the merged graph is built from the whole-genome inputs, not one part's split files
    raw = os.path.join(ROOT, "full", "raw")
    provenance = {}
    for key, name in (("gff3.gz", "gff3.gz"), ("genome.fa.gz", "genome.fa.gz"),
                      ("pep.fa.gz", "pep.fa.gz")):
        path = os.path.join(raw, name)
        if os.path.exists(path):
            provenance[key] = {
                "path": os.path.relpath(path, ROOT),
                "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest()
                if os.path.getsize(path) < 200_000_000 else None,
                "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                              time.gmtime(os.path.getmtime(path))),
                "retrieved_at_status": "download",
            }
    meta.update({
        "source_gff3": "gff3.gz",
        "source_genome": "genome.fa.gz",
        "source_pep": "pep.fa.gz",
        "source_retrieved_at": json.dumps(provenance, sort_keys=True),
        "chromosome": ",".join(chromosomes),
        "build_id": digest.hexdigest()[:16],
        "genome_length": str(final.execute(
            "SELECT SUM(length) FROM node WHERE kind='region'").fetchone()[0] or 0),
        "format": "ensembl-graph-m1-full",
    })
    final.execute("DELETE FROM meta")
    final.executemany("INSERT INTO meta(key,value) VALUES (?,?)", sorted(meta.items()))
    final.commit()

    counts = {
        "genes": final.execute("SELECT COUNT(*) FROM node WHERE kind='gene'").fetchone()[0],
        "transcripts": final.execute("SELECT COUNT(*) FROM node WHERE kind='transcript'").fetchone()[0],
        "translations": final.execute("SELECT COUNT(*) FROM node WHERE kind='translation'").fetchone()[0],
        "exons": final.execute("SELECT COUNT(*) FROM node WHERE kind='exon'").fetchone()[0],
        "xref_nodes": final.execute("SELECT COUNT(*) FROM node WHERE kind='xref'").fetchone()[0],
        "nodes": final.execute("SELECT COUNT(*) FROM node").fetchone()[0],
        "edges": final.execute("SELECT COUNT(*) FROM edge").fetchone()[0],
    }
    fingerprint = db_fingerprint(final)
    elapsed = round(time.time() - started, 1)
    stats = {"chromosomes": chromosomes, "counts": counts, "fingerprint": fingerprint,
             "elapsed_seconds": elapsed, "shared_keys_deduped": seen_as_primary,
             "path": os.path.relpath(args.out, ROOT)}
    with open(os.path.join(os.path.dirname(args.out), "stats.json"), "w", encoding="utf-8") as handle:
        json.dump(stats, handle, indent=2, sort_keys=True)
    print(json.dumps(stats, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
