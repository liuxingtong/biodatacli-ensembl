#!/usr/bin/env python3
"""Fetch whole-genome BioMart exports (same attributes, no chromosome filter).

The chr20 slice used four queries pinned to ``chromosome_name=20``.  A
whole-genome build needs the same columns for every chromosome, so this strips
that one filter line and re-runs the query against the same frozen archive.
"""

from __future__ import annotations

import os
import sys
import time
import urllib.parse
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
QUERIES = os.path.join(ROOT, "queries")
OUT = os.path.join(ROOT, "full", "biomart")
BIOMART = "https://jun2026.archive.ensembl.org/biomart/martservice"

TARGETS = {
    "pdb": "biomart_pdb.xml",
    "uniprot": "biomart_uniprot.xml",
    "gene_xref": "biomart_gene_xref.xml",
    "ccds_mim": "biomart_ccds_mim.xml",
}


def strip_filter(xml: str) -> str:
    lines = [line for line in xml.splitlines(keepends=True)
             if 'name="chromosome_name"' not in line]
    return "".join(lines)


def fetch(xml: str, destination: str, retries: int = 3) -> None:
    url = BIOMART + "?" + urllib.parse.urlencode({"query": xml})
    partial = destination + ".part"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=600) as response:
                total = 0
                with open(partial, "wb") as handle:
                    while True:
                        block = response.read(1 << 20)
                        if not block:
                            break
                        handle.write(block)
                        total += len(block)
            os.replace(partial, destination)
            print(f"  {os.path.basename(destination)}: {total:,} bytes")
            return
        except Exception as exc:  # noqa: BLE001 - network flakiness
            print(f"  attempt {attempt} failed: {exc}", file=sys.stderr)
            if attempt == retries:
                raise
            time.sleep(2 * attempt)


def main() -> int:
    os.makedirs(OUT, exist_ok=True)
    for name, query_file in TARGETS.items():
        path = os.path.join(QUERIES, query_file)
        if not os.path.exists(path):
            print(f"missing query {path}", file=sys.stderr)
            continue
        xml = strip_filter(open(path, encoding="utf-8").read())
        print(f"fetching {name} …", flush=True)
        fetch(xml, os.path.join(OUT, f"{name}.tsv"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
