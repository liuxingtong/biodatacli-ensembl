#!/usr/bin/env python3
"""Split the full Ensembl GFF3/FASTA into one file per chromosome.

A whole-genome build in one process would hold ~3 GB of sequence plus every
parsed feature at once.  Splitting first lets the existing ``build()`` run per
chromosome (a few hundred MB each), and a merge step stitches the SQLite files
together afterwards.

    python3 full/split.py
    # -> full/split/gff/chr20.gff3.gz, full/split/fa/chr20.fa.gz, ...
"""

from __future__ import annotations

import gzip
import os
import sys
from collections import OrderedDict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "full", "raw")
OUT = os.path.join(ROOT, "full", "split")

WANTED = [str(n) for n in range(1, 23)] + ["X", "Y", "MT"]


def split_gff(source: str, target_dir: str) -> dict:
    os.makedirs(target_dir, exist_ok=True)
    handles: "OrderedDict[str, gzip.GzipFile]" = OrderedDict()
    headers: list[str] = []
    counts: dict[str, int] = {}
    try:
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    headers.append(line)
                    continue
                chrom = line.split("\t", 1)[0]
                if chrom not in WANTED:
                    continue
                if chrom not in handles:
                    path = os.path.join(target_dir, f"chr{chrom}.gff3.gz")
                    stream = gzip.open(path, "wt", encoding="utf-8", compresslevel=6)
                    stream.writelines(headers)
                    handles[chrom] = stream
                    counts[chrom] = 0
                handles[chrom].write(line)
                counts[chrom] += 1
    finally:
        for stream in handles.values():
            stream.close()
    return counts


def split_fasta(source: str, target_dir: str) -> dict:
    os.makedirs(target_dir, exist_ok=True)
    handles: dict[str, gzip.GzipFile] = {}
    counts: dict[str, int] = {}
    chrom = None
    try:
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith(">"):
                    chrom = line[1:].split()[0]
                    if chrom not in WANTED:
                        chrom = None
                        continue
                    if chrom not in handles:
                        path = os.path.join(target_dir, f"chr{chrom}.fa.gz")
                        handles[chrom] = gzip.open(path, "wt", encoding="utf-8", compresslevel=6)
                        counts[chrom] = 0
                    handles[chrom].write(line)
                    counts[chrom] += len(line.strip())
                elif chrom is not None:
                    handles[chrom].write(line)
                    counts[chrom] += len(line.strip())
    finally:
        for stream in handles.values():
            stream.close()
    return counts


def main() -> int:
    gff = os.path.join(RAW, "gff3.gz")
    fa = os.path.join(RAW, "genome.fa.gz")
    for path in (gff, fa):
        if not os.path.exists(path):
            raise SystemExit(f"missing {path}")
    print("splitting GFF3 …")
    gff_counts = split_gff(gff, os.path.join(OUT, "gff"))
    print(f"  {len(gff_counts)} chromosomes, {sum(gff_counts.values()):,} feature lines")
    print("splitting FASTA …")
    fa_counts = split_fasta(fa, os.path.join(OUT, "fa"))
    total = sum(fa_counts.values())
    print(f"  {len(fa_counts)} sequences, {total:,} bases")
    missing = [chrom for chrom in WANTED if chrom not in gff_counts]
    if missing:
        print(f"  note: no GFF records for {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
