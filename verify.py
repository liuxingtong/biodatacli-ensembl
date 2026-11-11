#!/usr/bin/env python3
"""Offline verification of the Ensembl chr20 graph.

Checks, in order:
  1. frozen REST snapshot is still release 116
  2. gene/transcript versions and coordinates agree with REST
  3. canonical protein and cDNA sequences agree byte-for-byte with REST
  4. every CDS translates to the released peptide (or is a documented exception)
  5. 4DNK/SIFTS numbering and C-alpha coverage agree with the local files
  6. graph counts and the logical fingerprint match build_stats.json
  7. BioMart snapshots contain the P31946 / 4DNK edges
  8. versioned stable IDs are rejected by REST (the snapshot records this)

No network access is required.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ensembl_adapter import EnsemblAdapter
from ensembl_adapter.build import db_fingerprint
from ensembl_adapter.crossref import load_ca_positions, load_sifts

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "data", "raw")
REST = os.path.join(ROOT, "snapshots", "rest")


class Checker:
    def __init__(self):
        self.checks = []

    def check(self, name, condition, detail=""):
        status = "PASS" if condition else "FAIL"
        self.checks.append({"name": name, "status": status, "detail": detail})
        print(f"[{status}] {name}" + (f"  {detail}" if detail else ""))
        return condition

    def finish(self, out_path):
        failed = [c for c in self.checks if c["status"] == "FAIL"]
        result = {"checks": self.checks, "passed": len(self.checks) - len(failed),
                  "failed": len(failed)}
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
        print(f"\n{result['passed']} passed, {result['failed']} failed -> {out_path}")
        return 1 if failed else 0


def read_json(name):
    with open(os.path.join(REST, name + ".json"), encoding="utf-8") as handle:
        return json.load(handle)


def md5(text):
    return hashlib.md5(text.encode("utf-8")).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite"))
    parser.add_argument("--stats", default=os.path.join(ROOT, "data", "index", "build_stats.json"))
    parser.add_argument("--out", default=os.path.join(ROOT, "runs", "verify.json"))
    args = parser.parse_args(argv)

    c = Checker()
    adapter = EnsemblAdapter(args.db)
    meta = adapter.meta()
    stats = json.load(open(args.stats, encoding="utf-8"))

    # 1. release pinned
    software = read_json("info_software")
    data_info = read_json("info_data")
    c.check("snapshot release is 116",
            software.get("release") == 116 and 116 in data_info.get("releases", []),
            f"software={software.get('release')} data={data_info.get('releases')}")
    c.check("graph meta release/assembly",
            meta["release"] == "116" and meta["assembly"] == "GRCh38.p14",
            f"{meta['release']} {meta['assembly']}")

    assembly = read_json("info_assembly_human")
    c.check("assembly identity",
            assembly.get("assembly_name") == "GRCh38.p14"
            and assembly.get("assembly_accession") == "GCA_000001405.29",
            f"{assembly.get('assembly_name')} {assembly.get('assembly_accession')}")

    # 2. gene and transcripts
    lookup = read_json("lookup_ywhab_expand")
    gene = adapter.get("gene:ENSG00000166913")
    c.check("YWHAB gene coordinates/version",
            gene is not None
            and gene["version"] == lookup.get("version")
            and gene["start"] == lookup.get("start")
            and gene["end"] == lookup.get("end")
            and gene["strand"] == lookup.get("strand"),
            f"version={gene['version'] if gene else None} {gene['start'] if gene else None}-{gene['end'] if gene else None}")
    canonical = adapter.canonical_transcript(gene["stable_id"])
    c.check("canonical transcript",
            canonical is not None and lookup.get("canonical_transcript", "").split(".")[0] == canonical["stable_id"],
            f"{canonical['stable_id']}.{canonical['version']} vs REST {lookup.get('canonical_transcript')}")

    rest_tx = {t["id"]: t for t in lookup.get("Transcript", [])}
    graph_tx = {t["stable_id"]: t for t in adapter.gene_transcripts(gene["stable_id"])}
    c.check("transcript set",
            set(rest_tx) == set(graph_tx),
            f"REST={len(rest_tx)} graph={len(graph_tx)} diff={sorted(set(rest_tx) ^ set(graph_tx))[:5]}")
    version_mismatch = [tid for tid in rest_tx if tid in graph_tx
                        and rest_tx[tid].get("version") != graph_tx[tid]["version"]]
    c.check("transcript versions", not version_mismatch, f"mismatch={version_mismatch[:5]}")

    # 3. sequences
    protein_rest = read_json("sequence_protein_canonical")
    protein = adapter.transcript_translation(canonical["stable_id"])
    protein_seq = adapter.protein_seq(protein["stable_id"])
    c.check("canonical protein identity",
            protein_rest.get("seq") == protein_seq and len(protein_seq) == 246,
            f"{len(protein_seq)} aa, REST version={protein_rest.get('version')}")

    cdna_rest = read_json("sequence_cdna_canonical")
    cdna = adapter.cdna_seq(canonical["stable_id"])
    c.check("canonical cDNA reconstruction",
            cdna_rest.get("seq") == cdna,
            f"len={len(cdna)} md5={md5(cdna)}")

    tx_stats = adapter.transcript_stats(canonical["stable_id"])
    c.check("cDNA md5 stored",
            tx_stats["cdna_md5"] == md5(cdna) and tx_stats["cdna_length"] == len(cdna),
            f"{tx_stats['cdna_length']} bp")

    alt_cdna = read_json("sequence_cdna_alt")
    c.check("second isoform cDNA reconstruction",
            alt_cdna.get("seq") == adapter.cdna_seq("ENST00000372839"),
            f"len={len(adapter.cdna_seq('ENST00000372839'))}")

    # 4. proteome-wide translation check
    validation = stats["counts"]["validation"]
    total_coding = sum(v for k, v in validation.items() if k != "non_coding")
    c.check("all CDS translations accounted for",
            total_coding == 9585 and validation.get("mismatch", 0) == 0,
            f"coding={total_coding} exact={validation.get('exact')} alt_start={validation.get('alt_start')} "
            f"stop_as_X={validation.get('stop_as_X')} mismatch={validation.get('mismatch', 0)}")

    # 5. 4DNK / SIFTS
    sifts = load_sifts(os.path.join(RAW, "4dnk.xml.gz"), "P31946")
    ca = load_ca_positions(os.path.join(RAW, "4DNK.cif"))
    mapping = sifts.get("A", {})
    observed_auth = [p for p in range(230, 247) if mapping.get(p, {}).get("pdb_auth") is not None]
    auth_ok = all(mapping[p]["pdb_auth"] == p for p in observed_auth)
    label_ok = all(mapping.get(p, {}).get("pdb_label") == p + 1 for p in range(230, 247))
    c.check("SIFTS UniProt->PDB numbering", auth_ok and label_ok and len(observed_auth) == 4,
            f"observed auth=UniProt ({len(observed_auth)}), label=UniProt+1")
    a_obs = sum(1 for p in range(230, 247) if mapping.get(p, {}).get("pdb_auth") in ca.get("A", {}))
    b_map = sifts.get("B", {})
    b_obs = sum(1 for p in range(230, 247) if b_map.get(p, {}).get("pdb_auth") in ca.get("B", {}))
    c.check("4DNK C-alpha coverage", a_obs == 4 and b_obs == 5, f"A={a_obs}/17 B={b_obs}/17")

    # 6. graph counts / fingerprint
    counts = stats["counts"]
    c.check("node/edge counts",
            counts["nodes"] == 67001 and counts["edges"] == 173437,
            f"nodes={counts['nodes']} edges={counts['edges']}")
    c.check("logical fingerprint",
            db_fingerprint(adapter.conn) == stats["fingerprint"],
            stats["fingerprint"][:16])

    # 7. BioMart edges
    with open(os.path.join(RAW, "biomart", "uniprot.tsv"), encoding="utf-8") as handle:
        uniprot_rows = list(csv.DictReader(handle, delimiter="\t"))
    with open(os.path.join(RAW, "biomart", "pdb.tsv"), encoding="utf-8") as handle:
        pdb_rows = list(csv.DictReader(handle, delimiter="\t"))
    c.check("BioMart P31946 -> ENSP00000300161",
            any(r["Protein stable ID"] == "ENSP00000300161"
                and r["UniProtKB/Swiss-Prot ID"] == "P31946" for r in uniprot_rows))
    c.check("BioMart 4DNK -> ENSP00000300161",
            any(r["Protein stable ID"] == "ENSP00000300161" and r["PDB ID"] == "4DNK" for r in pdb_rows))
    c.check("graph has P31946 and 4DNK xrefs",
            adapter.get("xref:UniProt:P31946") is not None and adapter.get("xref:PDB:4DNK") is not None)

    # 8. M1 contract smoke checks
    c.check("graph schema is m1-contract-v1",
            meta.get("schema_version") == "m1-contract-v1", meta.get("schema_version", ""))
    pdb_items = adapter.neighbors("Ensembl:ENSP00000300161", relations=["pdb_xref"],
                                  direction="out", limit=50)["items"]
    four = [item for item in pdb_items if item["node"]["id"] == "PDB:4DNK"]
    c.check("PDB:4DNK neighbour has chains/range/source",
            bool(four) and four[0]["edge"].get("chains") == ["A", "B"]
            and four[0]["edge"].get("protein_range") == [1, 246]
            and bool(four[0]["source"].get("locator")),
            f"items={len(pdb_items)}")
    records = adapter.get_records(["Ensembl:ENSP00000300161"],
                                  fields=["protein_mapping"],
                                  region={"protein_range": [230, 246]})["records"][0]
    c.check("get_records 230-246 mapping",
            records["region"]["retrieved"] == list(range(230, 247))
            and not records["region"]["missing"],
            f"retrieved={len(records['region']['retrieved'])} missing={len(records['region']['missing'])}")
    capped = adapter.neighbors("Ensembl:ENSG00000166913", relations=["has_transcript"],
                               direction="out", limit=20)
    c.check("neighbors limit/truncated",
            capped["total"] == 40 and len(capped["items"]) == 20 and capped["truncated"],
            f"total={capped['total']} returned={len(capped['items'])}")

    # 9. versioned-ID behaviour (snapshot)
    error = read_json("error_versioned_id")
    manifest = json.load(open(os.path.join(REST, "manifest.json"), encoding="utf-8"))
    status = manifest["endpoints"].get("error_versioned_id", {}).get("status")
    c.check("REST rejects versioned stable IDs",
            status == 400 and "not found" in json.dumps(error).lower(),
            f"status={status} body={json.dumps(error)[:80]}")

    return c.finish(args.out)


if __name__ == "__main__":
    raise SystemExit(main())
