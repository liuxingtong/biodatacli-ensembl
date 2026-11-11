#!/usr/bin/env python3
"""M1 acceptance walk-through: P31946 -> Ensembl -> 4DNK.

This is the Ensembl-side equivalent of the PDF's P31946 example.  It uses only
the common M1 calls:

    resolve(ref, constraints)
    neighbors(node_id, relations, direction, limit)
    get_records(ids, fields, region)
    walk(seeds, relations, max_steps, budget)

The final JSON record has the same node / edge / source shape as the PDF
example, with Ensembl as the asserting source and 4DNK as the external node.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ensembl_adapter import EnsemblAdapter

ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--out", default=os.path.join(ROOT, "runs", "acceptance_m1.json"))
    args = parser.parse_args(argv)

    adapter = EnsemblAdapter(args.db)
    checks = []

    def check(name, condition, detail=""):
        checks.append({"name": name, "pass": bool(condition), "detail": detail})
        print(f"[{'PASS' if condition else 'FAIL'}] {name}  {detail}")
        return condition

    # 1. resolve the UniProt accession to a real node
    resolved = adapter.resolve({"text": "P31946", "type": "Protein"}, {"taxon_id": 9606})
    check("resolve P31946", resolved["status"] == "ok" and resolved["candidates"],
          resolved["candidates"][0]["id"] if resolved["candidates"] else resolved["status"])
    uniprot = resolved["candidates"][0]

    # 2. walk from the UniProt node back into Ensembl
    back = adapter.neighbors(uniprot["id"], relations=["uniprot_xref"], direction="in", limit=20)
    proteins = [item for item in back["items"] if item["node"]["type"] == "Protein"]
    check("UniProtKB:P31946 -> Ensembl proteins", len(proteins) >= 1,
          f"{len(proteins)} proteins")
    canonical = next((p for p in proteins if p["node"]["stable_id"] == "ENSP00000300161"), proteins[0])

    # 3. the Ensembl protein -> PDB structure edge
    pdb_edges = adapter.neighbors(canonical["node"]["id"], relations=["pdb_xref"],
                                  direction="out", limit=50)
    four = next((item for item in pdb_edges["items"] if item["node"]["id"] == "PDB:4DNK"), None)
    check("ENSP00000300161 -> PDB:4DNK", four is not None,
          f"{len(pdb_edges['items'])} PDB edges")
    if four:
        check("4DNK edge carries chains and range",
              four["edge"].get("chains") == ["A", "B"]
              and four["edge"].get("protein_range") == [1, 246],
              json.dumps({k: four["edge"].get(k) for k in ("chains", "protein_range")}))

    # 4. get_records: sequence + per-residue mapping for 230-246
    records = adapter.get_records(
        [canonical["node"]["id"]],
        fields=["sequence", "protein_mapping", "xrefs"],
        region={"protein_range": [230, 246]},
    )
    rec = records["records"][0]
    retrieved = rec["region"]["retrieved"] if rec.get("region") else []
    missing = rec["region"]["missing"] if rec.get("region") else []
    check("get_records 230-246", retrieved == list(range(230, 247)) and not missing,
          f"retrieved={len(retrieved)} missing={len(missing)}")
    seq_slice = rec["fields"]["sequence"].get("slice")
    check("protein slice", seq_slice == "WTSENQGDEGDAGEGEN", seq_slice or "")
    first = rec["region"]["positions"][0]
    check("residue 230 maps to chr20:44,906,385-387",
          first["genomic_interval"] == [44906385, 44906387] and first["codon"] == "TGG",
          f"{first['genomic_interval']} {first['codon']}")

    # 5. transcript view: exons + translation
    tx_rec = adapter.get_records(["Ensembl:ENST00000353703"],
                                 fields=["exons", "translation", "sequence"])["records"][0]
    check("transcript exons/translation",
          len(tx_rec["fields"]["exons"]) == 6
          and tx_rec["fields"]["translation"]["stable_id"] == "ENSP00000300161",
          f"exons={len(tx_rec['fields']['exons'])}")

    # 6. bounded random walk, reproducible and budget-limited
    relations = ["uniprot_xref", "pdb_xref", "has_translation", "has_exon"]
    walk_a = adapter.walk(["UniProtKB:P31946"], relations=relations,
                          max_steps=8, budget=10, rng_seed=20261106)
    walk_b = adapter.walk(["UniProtKB:P31946"], relations=relations,
                          max_steps=8, budget=10, rng_seed=20261106)
    walk = walk_a["walks"][0]
    check("walk reproducible", walk_a == walk_b, walk["stopped"])
    check("walk budget respected", len(walk["visited"]) <= 10, f"visited={len(walk['visited'])}")

    # PDF-shaped record
    pdf_record = None
    if four:
        pdf_record = {
            "node": {"id": four["node"]["id"], "type": four["node"]["type"]},
            "edge": four["edge"],
            "source": four["source"],
        }

    result = {
        "release": adapter.meta()["release"],
        "assembly": adapter.meta()["assembly"],
        "schema_version": adapter.meta()["schema_version"],
        "checks": checks,
        "passed": sum(1 for c in checks if c["pass"]),
        "failed": sum(1 for c in checks if not c["pass"]),
        "resolved": uniprot,
        "canonical_protein": canonical,
        "pdb_xref": pdf_record,
        "get_records": rec,
        "transcript_record": tx_rec,
        "walk": walk,
    }
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, ensure_ascii=False)

    if pdf_record:
        print("\nPDF-shaped neighbour record:")
        print(json.dumps(pdf_record, indent=2, sort_keys=True, ensure_ascii=False))
    print(f"\n{result['passed']} passed, {result['failed']} failed -> {args.out}")
    return 0 if result["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
