#!/usr/bin/env python3
"""Fetch the release-pinned Ensembl inputs used by the chr20 graph.

Everything is pinned to Ensembl release 116 (GRCh38.p14):

  * GFF3, genome FASTA and peptide FASTA for human chromosome 20
  * the release-115 GFF3 for the same chromosome (version-drift check)
  * four BioMart TSV exports (UniProt, PDB, HGNC/Entrez, CCDS/MIM)
  * 4DNK mmCIF + SIFTS XML (cross-database position check)

Usage:
    python3 fetch_data.py            # download missing files, verify checksums
    python3 fetch_data.py --check    # verify only
    python3 fetch_data.py --force    # re-download everything

The checksums are recorded in data/manifest.json.  The manifest is an input
artifact: it is not used by the graph builder except for provenance.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ensembl_adapter.provenance import iso_utc  # noqa: E402

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "data", "raw")
MANIFEST = os.path.join(ROOT, "data", "manifest.json")

FTP = "https://ftp.ensembl.org/pub/release-116"
FTP_115 = "https://ftp.ensembl.org/pub/release-115"
BIOMART = "https://jun2026.archive.ensembl.org/biomart/martservice"

FILES = [
    (
        os.path.join(RAW, "Homo_sapiens.GRCh38.116.chromosome.20.gff3.gz"),
        f"{FTP}/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.chromosome.20.gff3.gz",
    ),
    (
        os.path.join(RAW, "Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz"),
        f"{FTP}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz",
    ),
    (
        os.path.join(RAW, "Homo_sapiens.GRCh38.pep.all.fa.gz"),
        f"{FTP}/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz",
    ),
    (
        os.path.join(RAW, "Homo_sapiens.GRCh38.115.chromosome.20.gff3.gz"),
        f"{FTP_115}/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.chromosome.20.gff3.gz",
    ),
    (
        os.path.join(RAW, "4DNK.cif"),
        "https://files.rcsb.org/download/4DNK.cif",
    ),
    (
        os.path.join(RAW, "4dnk.xml.gz"),
        "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/4dnk.xml.gz",
    ),
    (
        os.path.join(RAW, "uniprot_P31946.fasta"),
        "https://rest.uniprot.org/uniprotkb/P31946.fasta",
    ),
]

BIOMART_QUERIES = [
    ("biomart_uniprot.xml", "uniprot.tsv"),
    ("biomart_gene_xref.xml", "gene_xref.tsv"),
    ("biomart_pdb.xml", "pdb.tsv"),
    ("biomart_ccds_mim.xml", "ccds_mim.tsv"),
    ("biomart_go.xml", "go.tsv"),
    ("biomart_reactome.xml", "reactome.tsv"),
]


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, dest: str, retries: int = 3) -> None:
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    partial = dest + ".part"
    for attempt in range(1, retries + 1):
        try:
            print(f"  GET {url}")
            req = urllib.request.Request(url, headers={"User-Agent": "biodata-ensembl-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=180) as response, open(partial, "wb") as out:
                total = response.headers.get("Content-Length")
                read = 0
                while True:
                    block = response.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
                    read += len(block)
            os.replace(partial, dest)
            print(f"  -> {dest} ({read} bytes)")
            return
        except Exception as exc:  # noqa: BLE001 - retry all network failures
            print(f"  attempt {attempt}/{retries} failed: {exc}", file=sys.stderr)
            if attempt == retries:
                raise
            time.sleep(2 * attempt)


def download_biomart(query_file: str, dest: str) -> None:
    query = open(query_file, encoding="utf-8").read()
    url = BIOMART + "?" + urllib.parse.urlencode({"query": query})
    download(url, dest)


def stamp(record: dict, path: str, *, downloaded: bool) -> dict:
    """Attach when these bytes were retrieved; never invent a date."""
    if downloaded:
        record["retrieved_at"] = iso_utc(time.time())
        record["retrieved_at_source"] = "download"
    elif os.path.exists(path):
        record.setdefault("retrieved_at", iso_utc(os.path.getmtime(path)))
        record.setdefault("retrieved_at_source", "file_mtime")
    return record


def load_manifest() -> dict:
    if os.path.exists(MANIFEST):
        with open(MANIFEST, encoding="utf-8") as handle:
            return json.load(handle)
    return {"release": "116", "files": {}}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify only; never download")
    parser.add_argument("--force", action="store_true", help="re-download even if present")
    args = parser.parse_args(argv)

    manifest = load_manifest()
    files = manifest.setdefault("files", {})
    failed = 0

    for dest, url in FILES:
        rel = os.path.relpath(dest, ROOT)
        record = files.get(rel)
        if os.path.exists(dest) and not args.force:
            actual = sha256(dest)
            if record and record.get("sha256") == actual:
                print(f"[ok] {rel} ({os.path.getsize(dest)} bytes)")
                continue
            if record:
                print(f"[FAIL] {rel} checksum mismatch", file=sys.stderr)
                failed += 1
                continue
            print(f"[new] {rel} present, recording checksum")
            files[rel] = stamp({"url": url, "sha256": actual,
                                "size": os.path.getsize(dest)}, dest, downloaded=False)
            continue
        if args.check:
            print(f"[missing] {rel}", file=sys.stderr)
            failed += 1
            continue
        download(url, dest)
        files[rel] = stamp({"url": url, "sha256": sha256(dest),
                            "size": os.path.getsize(dest)}, dest, downloaded=True)

    os.makedirs(os.path.join(RAW, "biomart"), exist_ok=True)
    for query_name, out_name in BIOMART_QUERIES:
        query_file = os.path.join(ROOT, "queries", query_name)
        dest = os.path.join(RAW, "biomart", out_name)
        rel = os.path.relpath(dest, ROOT)
        record = files.get(rel)
        if os.path.exists(dest) and not args.force:
            actual = sha256(dest)
            if record and record.get("sha256") == actual:
                print(f"[ok] {rel} ({os.path.getsize(dest)} bytes)")
                continue
            if record:
                print(f"[FAIL] {rel} checksum mismatch", file=sys.stderr)
                failed += 1
                continue
            print(f"[new] {rel} present, recording checksum")
            files[rel] = stamp({"url": BIOMART, "query": query_name,
                                "sha256": actual, "size": os.path.getsize(dest)},
                               dest, downloaded=False)
            continue
        if args.check:
            print(f"[missing] {rel}", file=sys.stderr)
            failed += 1
            continue
        download_biomart(query_file, dest)
        files[rel] = stamp({"url": BIOMART, "query": query_name,
                            "sha256": sha256(dest), "size": os.path.getsize(dest)},
                           dest, downloaded=True)

    # Files that predate per-file timestamps get their mtime, labelled as such.
    for rel, record in files.items():
        path = os.path.join(ROOT, rel)
        if not record.get("retrieved_at") and os.path.exists(path):
            record["retrieved_at"] = iso_utc(os.path.getmtime(path))
            record["retrieved_at_source"] = "file_mtime"

    manifest["release"] = "116"
    manifest["assembly"] = "GRCh38.p14"
    manifest["updated_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    with open(MANIFEST, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    print(f"manifest -> {MANIFEST}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
