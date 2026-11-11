"""Sequence reconstruction from genome + GFF3 (Ensembl semantics).

Ensembl's GFF3 CDS phase is the number of bases to skip at the start of a
CDS segment (in transcript orientation) to reach the first complete codon.
For a 5'-incomplete CDS the released peptide starts with ``X``; the same
convention is reproduced here.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {}
_i = 0
for _a in _BASES:
    for _b in _BASES:
        for _c in _BASES:
            CODON_TABLE[_a + _b + _c] = _AAS[_i]
            _i += 1

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def translate(seq: str) -> str:
    seq = seq.upper()
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def splice(genome: str, segments: Sequence, strand: int, offset: int = 0) -> str:
    """Concatenate genomic segments in transcript orientation.

    ``segments`` must already be ordered in transcript orientation.  Each
    segment needs ``start`` and ``end`` (1-based, inclusive) attributes.
    """
    parts: List[str] = []
    for seg in segments:
        start = seg["start"] if isinstance(seg, dict) else seg.start
        end = seg["end"] if isinstance(seg, dict) else seg.end
        chunk = genome[start - 1:end]
        if strand < 0:
            chunk = reverse_complement(chunk)
        parts.append(chunk)
    joined = "".join(parts)
    return joined[offset:] if offset else joined


def reconstruct_cds(genome: str, segments: Sequence, strand: int):
    """Return (cds_sequence, first_phase).

    The first segment's phase bases are trimmed; the caller prepends ``X``
    to the peptide when the phase is non-zero.
    """
    if not segments:
        return "", 0
    first_phase = int(getattr(segments[0], "phase", 0) or 0)
    seq = splice(genome, segments, strand, offset=first_phase)
    return seq, first_phase


def reconstruct_cdna(genome: str, exons: Sequence, strand: int) -> str:
    return splice(genome, exons, strand)
