"""Minimal Ensembl GFF3 reader.

Only the records needed for the gene-model graph are kept:
gene / transcript-like / exon / CDS.  The reader deliberately does not
hard-code transcript feature types: Ensembl uses ``mRNA`` for most coding
transcripts but ``lnc_RNA``, ``transcript``, ``miRNA``, ``snRNA`` ... for
others.  Any record carrying a ``transcript_id`` attribute is a transcript.
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional
from urllib.parse import unquote

GENE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
TRANSCRIPT_LEVEL_TYPES = {
    "mRNA",
    "transcript",
    "lnc_RNA",
    "pseudogenic_transcript",
    "unconfirmed_transcript",
    "miRNA",
    "snRNA",
    "snoRNA",
    "rRNA",
    "scRNA",
}


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_attributes(text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in text.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = unquote(value)
    return attrs


def _strip_prefix(value: str, prefix: str) -> str:
    if value.startswith(prefix):
        return value[len(prefix):]
    return value


def _parent_ids(attrs: Dict[str, str]) -> List[str]:
    raw = attrs.get("Parent", "")
    out: List[str] = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            item = item.split(":", 1)[1]
        out.append(item)
    return out


def _version(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(str(value).strip())
    except ValueError:
        return None


@dataclass
class Gene:
    gene_id: str
    name: str
    biotype: str
    seq_region: str
    start: int
    end: int
    strand: int
    version: Optional[int]
    source: str
    description: str = ""
    feature_type: str = "gene"
    line: int = 0


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    name: str
    biotype: str
    seq_region: str
    start: int
    end: int
    strand: int
    version: Optional[int]
    source: str
    feature_type: str
    tags: List[str] = field(default_factory=list)
    tsl: str = ""
    ccds: str = ""
    line: int = 0


@dataclass
class Exon:
    exon_id: str
    seq_region: str
    start: int
    end: int
    strand: int
    version: Optional[int]
    rank: int
    phase: Optional[int] = None
    end_phase: Optional[int] = None
    constitutive: Optional[int] = None
    line: int = 0


@dataclass
class CdsSegment:
    transcript_id: str
    protein_id: str
    protein_version: Optional[int]
    seq_region: str
    start: int
    end: int
    strand: int
    phase: int
    rank: int = 0
    line: int = 0


@dataclass
class GffData:
    genes: Dict[str, Gene]
    transcripts: Dict[str, Transcript]
    exons: Dict[str, Exon]  # unique ENSE records
    tx_exons: Dict[str, List[Exon]]  # transcript -> ordered exon instances
    tx_cds: Dict[str, List[CdsSegment]]
    feature_counts: Dict[str, int]


def load_gff3(path) -> GffData:
    genes: Dict[str, Gene] = {}
    transcripts: Dict[str, Transcript] = {}
    exons: Dict[str, Exon] = {}
    tx_exons_raw: Dict[str, List[Exon]] = defaultdict(list)
    tx_cds_raw: Dict[str, List[CdsSegment]] = defaultdict(list)
    feature_counts: Dict[str, int] = defaultdict(int)

    with open_text(path) as handle:
        for line_no, line in enumerate(handle, 1):
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            seq_region, source, ftype, start, end, _score, strand, phase, attr_text = cols
            attrs = parse_attributes(attr_text)
            feature_counts[ftype] += 1
            strand_i = 1 if strand == "+" else (-1 if strand == "-" else 0)
            start_i, end_i = int(start), int(end)

            if ftype in GENE_TYPES:
                gid = attrs.get("gene_id") or _strip_prefix(attrs.get("ID", ""), "gene:")
                if not gid:
                    continue
                genes[gid] = Gene(
                    gene_id=gid,
                    name=attrs.get("Name", gid),
                    biotype=attrs.get("biotype", ""),
                    seq_region=seq_region,
                    start=start_i,
                    end=end_i,
                    strand=strand_i,
                    version=_version(attrs.get("version")),
                    source=source,
                    description=attrs.get("description", ""),
                    feature_type=ftype,
                    line=line_no,
                )
                continue

            if "transcript_id" in attrs:
                tid = attrs["transcript_id"]
                parents = _parent_ids(attrs)
                gid = parents[0] if parents else ""
                tags = [t for t in attrs.get("tag", "").split(",") if t]
                transcripts[tid] = Transcript(
                    transcript_id=tid,
                    gene_id=gid,
                    name=attrs.get("Name", tid),
                    biotype=attrs.get("biotype", ""),
                    seq_region=seq_region,
                    start=start_i,
                    end=end_i,
                    strand=strand_i,
                    version=_version(attrs.get("version")),
                    source=source,
                    feature_type=ftype,
                    tags=tags,
                    tsl=attrs.get("transcript_support_level", ""),
                    ccds=attrs.get("ccdsid", ""),
                    line=line_no,
                )
                continue

            if ftype == "exon":
                eid = attrs.get("exon_id") or attrs.get("Name", "")
                if not eid:
                    continue
                rank = _version(attrs.get("rank")) or 0
                exon = Exon(
                    exon_id=eid,
                    seq_region=seq_region,
                    start=start_i,
                    end=end_i,
                    strand=strand_i,
                    version=_version(attrs.get("version")),
                    rank=rank,
                    phase=_version(attrs.get("ensembl_phase")),
                    end_phase=_version(attrs.get("ensembl_end_phase")),
                    constitutive=_version(attrs.get("constitutive")),
                    line=line_no,
                )
                exons.setdefault(eid, exon)
                for tid in _parent_ids(attrs):
                    tx_exons_raw[tid].append(exon)
                continue

            if ftype == "CDS":
                pid = attrs.get("protein_id")
                if not pid:
                    continue
                for tid in _parent_ids(attrs):
                    tx_cds_raw[tid].append(
                        CdsSegment(
                            transcript_id=tid,
                            protein_id=pid,
                            protein_version=_version(attrs.get("version")),
                            seq_region=seq_region,
                            start=start_i,
                            end=end_i,
                            strand=strand_i,
                            phase=_version(phase) or 0,
                            line=line_no,
                        )
                    )

    # Order exon instances and CDS segments in transcript orientation.
    tx_exons: Dict[str, List[Exon]] = {}
    for tid, records in tx_exons_raw.items():
        tx = transcripts.get(tid)
        strand = tx.strand if tx else (records[0].strand if records else 1)
        records = sorted(records, key=lambda e: e.rank)
        # Some GFF3 writers omit rank; fall back to genomic order.
        if all(e.rank == 0 for e in records):
            records.sort(key=lambda e: (e.start, e.end), reverse=(strand < 0))
        tx_exons[tid] = records

    tx_cds: Dict[str, List[CdsSegment]] = {}
    for tid, records in tx_cds_raw.items():
        tx = transcripts.get(tid)
        strand = tx.strand if tx else (records[0].strand if records else 1)
        records.sort(key=lambda c: (c.start, c.end), reverse=(strand < 0))
        for i, seg in enumerate(records, start=1):
            seg.rank = i
        tx_cds[tid] = records

    return GffData(
        genes=genes,
        transcripts=transcripts,
        exons=exons,
        tx_exons=tx_exons,
        tx_cds=tx_cds,
        feature_counts=dict(feature_counts),
    )
