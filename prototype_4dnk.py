#!/usr/bin/env python3
"""Minimal prototype: UniProt P31946 230-246 -> Ensembl -> PDB 4DNK.

This is deliberately a single-gene, single-structure walk-through, mirroring
the 4DNK/SIFTS prototype from the PDB phase.  It answers four questions with
one reproducible script:

  1. Which Ensembl gene/transcript/protein is P31946?
  2. Where do residues 230-246 live in the genome (codon, exon)?
  3. What are the corresponding 4DNK author/label residue numbers?
  4. Which positions actually have a C-alpha in chains A and B?
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ensembl_adapter import EnsemblAdapter
from ensembl_adapter.crossref import load_ca_positions, load_sifts, sifts_metadata
from ensembl_adapter.fasta import iter_fasta

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "data", "raw")


def read_single_fasta(path):
    for header, seq in iter_fasta(path):
        return header, seq
    return "", ""


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite"))
    parser.add_argument("--cif", default=os.path.join(RAW, "4DNK.cif"))
    parser.add_argument("--sifts", default=os.path.join(RAW, "4dnk.xml.gz"))
    parser.add_argument("--uniprot-fasta", default=os.path.join(RAW, "uniprot_P31946.fasta"))
    parser.add_argument("--out", default=os.path.join(ROOT, "runs", "prototype_4dnk.json"))
    parser.add_argument("--start", type=int, default=230)
    parser.add_argument("--end", type=int, default=246)
    args = parser.parse_args(argv)

    adapter = EnsemblAdapter(args.db)
    gene = adapter.resolve({"text": "YWHAB", "type": "Gene"}, {"taxon_id": 9606})["candidates"][0]
    canonical = adapter.canonical_transcript(gene["stable_id"])
    protein = adapter.transcript_translation(canonical["stable_id"])
    protein_seq = adapter.protein_seq(protein["stable_id"])
    _, uniprot_seq = read_single_fasta(args.uniprot_fasta)

    sifts = load_sifts(args.sifts, "P31946")
    chains = sorted(sifts)
    ca = load_ca_positions(args.cif)
    pdb_chains = sorted(ca)

    rows = []
    for pos in range(args.start, args.end + 1):
        mapping = adapter.map_protein_pos(canonical["stable_id"], pos)
        per_chain = {}
        for chain in pdb_chains:
            s = sifts.get(chain, {}).get(pos, {})
            auth = s.get("pdb_auth")
            per_chain[chain] = {
                "auth": auth,
                "label": s.get("pdb_label"),
                "ca": bool(auth is not None and auth in ca.get(chain, {})),
            }
        row = {
            "uniprot_pos": pos,
            "aa": protein_seq[pos - 1] if pos <= len(protein_seq) else None,
            "codon": mapping.get("codon"),
            "genomic": mapping.get("genomic_interval"),
            "exon_rank": (mapping.get("exon_ranks") or [None])[0],
            "exon_id": (mapping.get("exon_ids") or [None])[0],
            "pdb": per_chain,
        }
        rows.append(row)

    observed = {chain: [r["uniprot_pos"] for r in rows if r["pdb"][chain]["ca"]] for chain in pdb_chains}
    missing = {chain: [r["uniprot_pos"] for r in rows if not r["pdb"][chain]["ca"]] for chain in pdb_chains}

    isoform = adapter.gene_coverage(gene["stable_id"], args.start, args.end)
    iso_counts = Counter(entry["status"] for entry in isoform["entries"])
    covered = [e for e in isoform["entries"] if e["status"] == "covered"]
    beyond = [e for e in isoform["entries"] if e["status"] == "beyond_protein"]
    partial = [e for e in isoform["entries"] if e["status"] == "partial_beyond"]
    no_protein = [e for e in isoform["entries"] if e["status"] == "no_translation"]

    summary = {
        "gene": gene["stable_id"],
        "gene_global_id": gene["id"],
        "gene_version": gene["version"],
        "gene_label": gene["label"],
        "gene_location": f"{gene['seq_region']}:{gene['start']}-{gene['end']}:{'+' if gene['strand'] > 0 else '-'}",
        "transcripts_total": len(isoform["entries"]),
        "canonical_transcript": canonical["stable_id"],
        "canonical_global_id": canonical["id"],
        "canonical_version": canonical["version"],
        "canonical_exons": len(adapter.transcript_exons(canonical["stable_id"])),
        "protein": protein["stable_id"],
        "protein_global_id": protein["id"],
        "protein_version": protein["version"],
        "protein_length": protein["length"],
        "uniprot_identical": protein_seq == uniprot_seq,
        "window": [args.start, args.end],
        "window_covered_by_canonical": sum(1 for r in rows if r["codon"]),
        "pdb_chains": pdb_chains,
        "pdb_observed_counts": {chain: len(observed[chain]) for chain in pdb_chains},
        "pdb_missing": missing,
        "isoform_status": dict(sorted(iso_counts.items())),
        "isoform_covered": len(covered),
        "isoform_beyond": len(beyond),
        "isoform_partial_beyond": len(partial),
        "isoform_no_translation": len(no_protein),
    }

    result = {
        "release": adapter.meta()["release"],
        "assembly": adapter.meta()["assembly"],
        "sifts": sifts_metadata(args.sifts),
        "summary": summary,
        "rows": rows,
        "isoforms": [
            {
                "transcript_id": e["transcript_id"],
                "version": e["transcript_version"],
                "protein_length": e["protein_length"],
                "status": e["status"],
                "is_canonical": e["is_canonical"],
            }
            for e in isoform["entries"]
        ],
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, ensure_ascii=False)

    print(f"# P31946 {args.start}-{args.end} -> Ensembl release {result['release']} -> 4DNK")
    print()
    print(f"gene      {gene['stable_id']}.{gene['version']}  {gene['label']}  {summary['gene_location']}")
    print(f"canonical {canonical['stable_id']}.{canonical['version']}  {len(adapter.transcript_exons(canonical['stable_id']))} exons")
    print(f"protein   {protein['stable_id']}.{protein['version']}  {protein['length']} aa  identical_to_UniProt={protein_seq == uniprot_seq}")
    print(f"4DNK      chains {pdb_chains}  model 1")
    print()
    header = "| pos | aa | codon | genomic (chr20) | exon | auth A | auth B | label | CA A | CA B |"
    print(header)
    print("|---:|:--:|:-----:|-----------------|-----:|-------:|-------:|------:|:----:|:----:|")
    for row in rows:
        genomic = "-".join(str(x) for x in row["genomic"]) if row["genomic"] else "-"
        auth_a = row["pdb"].get("A", {}).get("auth")
        auth_b = row["pdb"].get("B", {}).get("auth")
        auth_a = "-" if auth_a is None else auth_a
        auth_b = "-" if auth_b is None else auth_b
        label_a = row["pdb"].get("A", {}).get("label")
        label_a = "-" if label_a is None else label_a
        print(
            f"| {row['uniprot_pos']} | {row['aa']} | {row['codon'] or '-'} | {genomic} | "
            f"{row['exon_rank'] or '-'} | {auth_a} | {auth_b} | {label_a} | "
            f"{'Y' if row['pdb'].get('A', {}).get('ca') else '-'} | "
            f"{'Y' if row['pdb'].get('B', {}).get('ca') else '-'} |"
        )
    print()
    print(f"canonical model covers {summary['window_covered_by_canonical']}/{len(rows)} positions")
    print(f"4DNK C-alpha observed: " + ", ".join(f"{k} {v}/{len(rows)}" for k, v in summary["pdb_observed_counts"].items()))
    for chain in pdb_chains:
        print(f"  missing {chain}: {missing[chain]}")
    print()
    print(f"YWHAB isoforms for this window: {summary['isoform_status']}")
    print(f"  covered={summary['isoform_covered']}  beyond_protein={summary['isoform_beyond']}  "
          f"partial_beyond={summary['isoform_partial_beyond']}  no_translation={summary['isoform_no_translation']}")
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
