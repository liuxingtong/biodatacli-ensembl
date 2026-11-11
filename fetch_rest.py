#!/usr/bin/env python3
"""Snapshot the small set of Ensembl REST responses used for verification.

The graph build never calls REST.  These JSON files are frozen inputs so the
verification step stays reproducible after rest.ensembl.org moves to release
117.  Re-run with --force only when intentionally refreshing the snapshot.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
import urllib.error
import urllib.request

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "snapshots", "rest")
BASE = "https://rest.ensembl.org"

ENDPOINTS = {
    "info_software": "/info/software",
    "info_data": "/info/data",
    "info_assembly_human": "/info/assembly/homo_sapiens",
    "lookup_symbol_ywhab": "/lookup/symbol/homo_sapiens/YWHAB",
    "lookup_ywhab_expand": "/lookup/id/ENSG00000166913?expand=1",
    "xrefs_ensp_p31946": "/xrefs/id/ENSP00000300161",
    "sequence_cdna_canonical": "/sequence/id/ENST00000353703?type=cdna",
    "sequence_protein_canonical": "/sequence/id/ENSP00000300161?type=protein",
    "sequence_cdna_alt": "/sequence/id/ENST00000372839?type=cdna",
    "sequence_protein_short": "/sequence/id/ENST00000428262?type=protein",
    "overlap_ywhab": "/overlap/region/homo_sapiens/20:44885676-44908715?feature=gene;feature=transcript;feature=exon;feature=cds",
    "error_versioned_id": "/lookup/id/ENSG00000166913.15",
    "homology_ywhab": "/homology/id/homo_sapiens/ENSG00000166913?type=orthologues",
    "overlap_variation_ywhab": "/overlap/region/human/20:44885676-44908715?feature=variation",
    "overlap_regulatory_ywhab": "/overlap/region/human/20:44885676-44908715?feature=regulatory",
    "vep_rs1042522": "/vep/human/id/rs1042522",
}


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def fetch(path: str):
    url = BASE + path
    req = urllib.request.Request(
        url,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "biodata-ensembl-snapshot/1.0",
        },
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as response:
            body = response.read().decode("utf-8")
            return response.status, dict(response.headers), body
    except urllib.error.HTTPError as exc:
        return exc.code, dict(exc.headers), exc.read().decode("utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)

    manifest_path = os.path.join(OUT, "manifest.json")
    manifest = {"base": BASE, "endpoints": {}}
    if os.path.exists(manifest_path) and not args.force:
        with open(manifest_path, encoding="utf-8") as handle:
            manifest = json.load(handle)

    release = None
    for name, path in ENDPOINTS.items():
        file_name = name + ".json"
        dest = os.path.join(OUT, file_name)
        if os.path.exists(dest) and not args.force:
            text = open(dest, encoding="utf-8").read()
            print(f"[ok] {file_name}")
            if name == "info_software":
                try:
                    release = json.loads(text).get("release")
                except json.JSONDecodeError:
                    pass
            continue
        status, headers, text = fetch(path)
        with open(dest, "w", encoding="utf-8") as handle:
            handle.write(text)
        manifest["endpoints"][name] = {
            "url": BASE + path,
            "status": status,
            "sha256": sha256_text(text),
            "bytes": len(text.encode("utf-8")),
            "content_type": headers.get("Content-Type", ""),
        }
        print(f"[{status}] {file_name} ({len(text)} bytes)")
        if name == "info_software":
            try:
                release = json.loads(text).get("release")
            except json.JSONDecodeError:
                pass

    manifest["release"] = release
    manifest["fetched_at"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    with open(manifest_path, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
    print(f"release={release} manifest={manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
