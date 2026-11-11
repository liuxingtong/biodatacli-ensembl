"""Tests for the offline Ensembl adapter.

Run from ensembl_proto/:
    python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter import EnsemblAdapter
from ensembl_adapter.build import build, db_fingerprint
from ensembl_adapter.sequence import reverse_complement


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


class SyntheticBuildTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="ensembl_test_")
        self.genome = "ACGT" * 50
        cds_plus = "ATGAAACCCGGGTTTAAATAG"
        cds_minus = reverse_complement(cds_plus)
        genome = list(self.genome)
        genome[10:31] = list(cds_plus)       # 11-31, 1-based
        genome[104:125] = list(cds_minus)    # 105-125, 1-based
        self.genome = "".join(genome)

        write(os.path.join(self.tmp, "genome.fa"), ">1\n" + self.genome + "\n")
        gff = "\n".join([
            "##gff-version 3",
            "##sequence-region 1 1 200",
            "1\ttest\tgene\t5\t45\t.\t+\t.\tID=gene:G1;Name=AAA;biotype=protein_coding;gene_id=G1;version=1",
            "1\ttest\tmRNA\t5\t45\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;Name=AAA-201;biotype=protein_coding;transcript_id=T1;tag=Ensembl_canonical;version=1",
            "1\ttest\texon\t5\t45\t.\t+\t.\tParent=transcript:T1;exon_id=E1;rank=1;version=1",
            "1\ttest\tCDS\t11\t31\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1;protein_id=P1;version=1",
            "1\ttest\tgene\t100\t130\t.\t-\t.\tID=gene:G2;Name=BBB;biotype=protein_coding;gene_id=G2;version=1",
            "1\ttest\tmRNA\t100\t130\t.\t-\t.\tID=transcript:T2;Parent=gene:G2;Name=BBB-201;biotype=protein_coding;transcript_id=T2;tag=Ensembl_canonical;version=1",
            "1\ttest\texon\t100\t130\t.\t-\t.\tParent=transcript:T2;exon_id=E2;rank=1;version=1",
            "1\ttest\tCDS\t105\t125\t.\t-\t0\tID=CDS:P2;Parent=transcript:T2;protein_id=P2;version=1",
        ]) + "\n"
        write(os.path.join(self.tmp, "test.gff3"), gff)
        pep = (
            ">P1.1 pep chromosome:test:1:5:45:1 gene:G1.1 transcript:T1.1 gene_symbol:AAA\nMKPGFK\n"
            ">P2.1 pep chromosome:test:1:100:130:-1 gene:G2.1 transcript:T2.1 gene_symbol:BBB\nMKPGFK\n"
        )
        write(os.path.join(self.tmp, "pep.fa"), pep)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def build(self, name):
        db = os.path.join(self.tmp, name + ".sqlite")
        stats = os.path.join(self.tmp, name + ".json")
        args = SimpleNamespace(
            gff3=os.path.join(self.tmp, "test.gff3"),
            genome=os.path.join(self.tmp, "genome.fa"),
            pep=os.path.join(self.tmp, "pep.fa"),
            biomart=os.path.join(self.tmp, "biomart"),
            enrichment=None, capabilities=None,
            db=db, stats=stats, chromosome="1", release="test", assembly="test",
        )
        result = build(args)
        return db, result

    def test_graph_and_translation(self):
        db, result = self.build("a")
        self.assertEqual(result["counts"]["genes"], 2)
        self.assertEqual(result["counts"]["transcripts"], 2)
        self.assertEqual(result["counts"]["exons_unique"], 2)
        self.assertEqual(result["counts"]["translations"], 2)
        self.assertEqual(result["counts"]["validation"].get("exact"), 2)

        adapter = EnsemblAdapter(db)
        gene = adapter.resolve({"text": "AAA", "type": "Gene"})["candidates"][0]
        self.assertEqual(gene["stable_id"], "G1")
        tx = adapter.canonical_transcript("G1")
        self.assertEqual(tx["stable_id"], "T1")
        self.assertEqual(adapter.protein_seq("P1"), "MKPGFK")

        first = adapter.map_protein_pos("T1", 1)
        self.assertEqual(first["codon"], "ATG")
        self.assertEqual(first["genomic_positions"], [11, 12, 13])
        last = adapter.map_protein_pos("T1", 6)
        self.assertEqual(last["codon"], "AAA")
        self.assertEqual(last["genomic_interval"], [26, 28])

        minus_first = adapter.map_protein_pos("T2", 1)
        self.assertEqual(minus_first["codon"], "ATG")
        self.assertEqual(minus_first["genomic_positions"], [125, 124, 123])
        self.assertEqual(minus_first["genomic_interval"], [123, 125])
        self.assertIn("ATG", adapter.cdna_seq("T2"))

    def test_deterministic_fingerprint(self):
        db_a, result_a = self.build("a")
        db_b, result_b = self.build("b")
        self.assertEqual(result_a["fingerprint"], result_b["fingerprint"])
        a = EnsemblAdapter(db_a)
        b = EnsemblAdapter(db_b)
        self.assertEqual(db_fingerprint(a.conn), db_fingerprint(b.conn))
        self.assertEqual(a.walk("T1", max_steps=4, budget=10, rng_seed=7),
                         a.walk("T1", max_steps=4, budget=10, rng_seed=7))


@unittest.skipUnless(
    os.path.exists(os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")),
    "chr20 index not built",
)
class RealIndexTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.adapter = EnsemblAdapter(os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite"))

    def test_ywhab_canonical(self):
        gene = self.adapter.resolve({"text": "YWHAB", "type": "Gene"})["candidates"][0]
        self.assertEqual(gene["stable_id"], "ENSG00000166913")
        self.assertEqual(gene["version"], 15)
        tx = self.adapter.canonical_transcript(gene["stable_id"])
        self.assertEqual(tx["stable_id"], "ENST00000353703")
        protein = self.adapter.transcript_translation(tx["stable_id"])
        self.assertEqual(protein["stable_id"], "ENSP00000300161")
        self.assertEqual(protein["length"], 246)

    def test_window_230_246(self):
        result = self.adapter.gene_coverage("ENSG00000166913", 230, 246)
        counts = {}
        for entry in result["entries"]:
            counts[entry["status"]] = counts.get(entry["status"], 0) + 1
        self.assertEqual(counts, {"covered": 23, "beyond_protein": 12,
                                  "partial_beyond": 3, "no_translation": 2})
        first = self.adapter.map_protein_pos("ENST00000353703", 230)
        self.assertEqual(first["genomic_interval"], [44906385, 44906387])
        self.assertEqual(first["codon"], "TGG")

    def test_phase_gt_zero_first_residue(self):
        # ENST00000413737 has a 5'-incomplete CDS (phase=1); Ensembl emits X.
        tx = "ENST00000413737"
        protein = self.adapter.transcript_translation(tx)
        seq = self.adapter.protein_seq(protein["stable_id"])
        first = self.adapter.map_protein_pos(tx, 1)
        self.assertEqual(first["residue"], "X")
        self.assertTrue(first["partial_codon"])
        second = self.adapter.map_protein_pos(tx, 2)
        self.assertEqual(second["residue"], seq[1])

    def test_pdb_anchor(self):
        self.assertIsNotNone(self.adapter.get("xref:PDB:4DNK"))
        self.assertIsNotNone(self.adapter.get("xref:UniProt:P31946"))
        neighbors = self.adapter.neighbors("PDB:4DNK", relations=["pdb_xref"], direction="in")
        self.assertTrue(any(n["node"]["id"] == "Ensembl:ENSP00000300161"
                            for n in neighbors["items"]))

    def test_walk_budget_and_reproducibility(self):
        first = self.adapter.walk("4DNK", relations=["pdb_xref", "has_translation", "has_exon"],
                                  max_steps=20, budget=8, rng_seed=11)
        walk = first["walks"][0]
        self.assertLessEqual(len(walk["visited"]), 8)
        self.assertEqual(walk["stopped"], "budget")
        self.assertEqual(first, self.adapter.walk("4DNK", relations=["pdb_xref", "has_translation", "has_exon"],
                                                  max_steps=20, budget=8, rng_seed=11))


if __name__ == "__main__":
    unittest.main()
