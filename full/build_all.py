#!/usr/bin/env python3
"""Build one SQLite graph per chromosome, then merge them.

``build.build()`` itself is unchanged: each run sees one chromosome's GFF3/FASTA
and writes a complete little database.  Merging renumbers nodes (the graph uses
integer keys internally) and de-duplicates shared external references such as
``PDB:4DNK``, which many proteins on different chromosomes point at.

    python3 full/build_all.py --chromosomes 21,22     # pilot
    python3 full/build_all.py                          # everything
    python3 full/merge.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter.build import build  # noqa: E402

SPLIT = os.path.join(ROOT, "full", "split")
PARTS = os.path.join(ROOT, "full", "parts")
FULL = os.path.join(ROOT, "full")
WANTED = [str(n) for n in range(1, 23)] + ["X", "Y", "MT"]


def chromosomes() -> list[str]:
    found = []
    for chrom in WANTED:
        gff = os.path.join(SPLIT, "gff", f"chr{chrom}.gff3.gz")
        fa = os.path.join(SPLIT, "fa", f"chr{chrom}.fa.gz")
        if os.path.exists(gff) and os.path.exists(fa):
            found.append(chrom)
    return found


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chromosomes", default="", help="comma list, default: all available")
    parser.add_argument("--pep", default=os.path.join(FULL, "raw", "pep.fa.gz"))
    parser.add_argument("--biomart", default=os.path.join(FULL, "biomart"))
    parser.add_argument("--enrichment", default=os.path.join(ROOT, "data", "raw", "xref_enrichment.json"))
    parser.add_argument("--capabilities", default=os.path.join(ROOT, "capabilities.json"))
    parser.add_argument("--release", default="116")
    parser.add_argument("--assembly", default="GRCh38.p14")
    args = parser.parse_args(argv)

    targets = [c.strip() for c in args.chromosomes.split(",") if c.strip()] or chromosomes()
    if not targets:
        raise SystemExit("nothing to build: run full/split.py first")
    os.makedirs(PARTS, exist_ok=True)

    started = time.time()
    summary = []
    for index, chrom in enumerate(targets, start=1):
        gff = os.path.join(SPLIT, "gff", f"chr{chrom}.gff3.gz")
        fa = os.path.join(SPLIT, "fa", f"chr{chrom}.fa.gz")
        db = os.path.join(PARTS, f"chr{chrom}.sqlite")
        stats = os.path.join(PARTS, f"chr{chrom}.stats.json")
        part_started = time.time()
        print(f"[{index}/{len(targets)}] chr{chrom} …", flush=True)
        result = build(SimpleNamespace(
            gff3=gff, genome=fa, pep=args.pep, biomart=args.biomart,
            db=db, stats=stats, enrichment=args.enrichment,
            capabilities=args.capabilities, chromosome=chrom,
            release=args.release, assembly=args.assembly,
        ))
        elapsed = time.time() - part_started
        size = os.path.getsize(db) / 2**20
        print(f"    genes={result['counts']['genes']:,} nodes={result['counts']['nodes']:,} "
              f"edges={result['counts']['edges']:,} {size:.0f} MiB {elapsed:.0f}s", flush=True)
        summary.append({"chromosome": chrom, "genes": result["counts"]["genes"],
                        "nodes": result["counts"]["nodes"], "edges": result["counts"]["edges"],
                        "seconds": round(elapsed, 1), "db": db, "stats": stats})

    payload = {"parts": summary, "chromosomes": targets,
               "seconds": round(time.time() - started, 1)}
    with open(os.path.join(PARTS, "build_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    print(f"\ndone in {payload['seconds']}s -> {os.path.join(PARTS, 'build_summary.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
