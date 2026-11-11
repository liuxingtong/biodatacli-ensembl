#!/usr/bin/env python3
"""Build the small cross-reference enrichment file used by the graph.

The BioMart PDB export gives only the PDB ID, not the chain list or the mapped
UniProt range.  For the M1 acceptance case (P31946 / 4DNK) we derive those from
the frozen SIFTS XML and the mmCIF file, and record the UniProt entry version
from the frozen REST xref snapshot.

Output: data/raw/xref_enrichment.json
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ensembl_adapter.crossref import load_ca_positions, load_sifts

ROOT = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(ROOT, "data", "raw")
OUT = os.path.join(RAW, "xref_enrichment.json")


def main():
    sifts = load_sifts(os.path.join(RAW, "4dnk.xml.gz"), "P31946")
    ca = load_ca_positions(os.path.join(RAW, "4DNK.cif"))
    chains = sorted(sifts)
    positions = sorted({pos for chain in sifts.values() for pos in chain})
    protein_range = [min(positions), max(positions)]

    snapshot = json.load(open(os.path.join(ROOT, "snapshots", "rest", "xrefs_ensp_p31946.json"),
                              encoding="utf-8"))
    uniprot_meta = next(
        (x for x in snapshot
         if x.get("dbname") == "Uniprot/SWISSPROT" and x.get("primary_id") == "P31946"),
        {},
    )

    window = range(230, 247)
    observed = {
        chain: [pos for pos in window
                if sifts.get(chain, {}).get(pos, {}).get("pdb_auth") in ca.get(chain, {})]
        for chain in chains
    }

    enrichment = {
        "pdb": {
            "4DNK": {
                "uniprot": "P31946",
                "chains": chains,
                "protein_range": protein_range,
                "observed_window": observed,
                "evidence": "SIFTS 2026-11-01 + 4DNK.cif model 1",
            }
        },
        "uniprot": {
            "P31946": {
                "entry_version": uniprot_meta.get("version"),
                "ensembl_range": [uniprot_meta.get("ensembl_start"), uniprot_meta.get("ensembl_end")],
                "xref_range": [uniprot_meta.get("xref_start"), uniprot_meta.get("xref_end")],
                "evidence": "snapshots/rest/xrefs_ensp_p31946.json",
            }
        },
    }
    with open(OUT, "w", encoding="utf-8") as handle:
        json.dump(enrichment, handle, indent=2, sort_keys=True)
    print(json.dumps(enrichment, indent=2, sort_keys=True))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
