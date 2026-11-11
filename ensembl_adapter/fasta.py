"""Small FASTA helpers for Ensembl genome / peptide files."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from typing import Dict, Iterator, Optional, Tuple

_ID_RE = re.compile(r"^(\S+?)(?:\.(\d+))?\s")


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def iter_fasta_records(path) -> Iterator[Tuple[str, str, int]]:
    header: Optional[str] = None
    header_line = 0
    chunks = []
    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks), header_line
                header = line[1:].strip()
                header_line = line_no
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks), header_line


def iter_fasta(path) -> Iterator[Tuple[str, str]]:
    for header, seq, _line in iter_fasta_records(path):
        yield header, seq


def split_id_version(token: str) -> Tuple[str, Optional[int]]:
    if "." in token:
        base, _, ver = token.rpartition(".")
        if ver.isdigit():
            return base, int(ver)
    return token, None


def load_genome(path) -> Dict[str, str]:
    """Load all records.  The chr20 file has a single record."""
    out: Dict[str, str] = {}
    for header, seq in iter_fasta(path):
        name = header.split()[0]
        out[name] = seq.upper()
    return out


@dataclass
class Peptide:
    protein_id: str
    version: Optional[int]
    seq: str
    gene_id: str
    transcript_id: str
    gene_symbol: str
    seq_region: str
    start: Optional[int]
    end: Optional[int]
    strand: Optional[int]
    biotype: str
    header_line: int = 0

    @property
    def length(self) -> int:
        return len(self.seq)


def _token(header: str, key: str) -> str:
    match = re.search(r"(?:^|\s)" + re.escape(key) + r":(\S+)", header)
    return match.group(1) if match else ""


def parse_peptide(header: str, seq: str, header_line: int = 0) -> Peptide:
    first = header.split()[0]
    protein_id, version = split_id_version(first)
    gene_raw = _token(header, "gene")
    tx_raw = _token(header, "transcript")
    gene_id, _ = split_id_version(gene_raw) if gene_raw else ("", None)
    tx_id, _ = split_id_version(tx_raw) if tx_raw else ("", None)
    chrom_raw = _token(header, "chromosome")
    region, start, end, strand = "", None, None, None
    if chrom_raw:
        parts = chrom_raw.split(":")
        if len(parts) >= 5:
            region = parts[1]
            try:
                start, end = int(parts[2]), int(parts[3])
                strand = 1 if parts[4] == "1" else -1
            except ValueError:
                pass
    return Peptide(
        protein_id=protein_id,
        version=version,
        seq=seq.upper(),
        gene_id=gene_id,
        transcript_id=tx_id,
        gene_symbol=_token(header, "gene_symbol"),
        seq_region=region,
        start=start,
        end=end,
        strand=strand,
        biotype=_token(header, "transcript_biotype") or _token(header, "gene_biotype"),
        header_line=header_line,
    )


def load_peptides(path) -> Dict[str, Peptide]:
    out: Dict[str, Peptide] = {}
    for header, seq, header_line in iter_fasta_records(path):
        pep = parse_peptide(header, seq, header_line=header_line)
        out[pep.protein_id] = pep
    return out
