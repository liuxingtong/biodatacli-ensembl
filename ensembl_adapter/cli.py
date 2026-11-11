"""Command-line interface for the offline Ensembl M1 adapter."""

from __future__ import annotations

import argparse
import json
import os
import sys

from . import EnsemblAdapter, api

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")


def _adapter(args):
    return EnsemblAdapter(args.db, genome_path=args.genome)


def _first(adapter, text, node_type=None):
    ref = {"text": text}
    if node_type:
        ref["type"] = node_type
    result = adapter.resolve(ref, {"taxon_id": 9606})
    if not result.get("candidates"):
        raise SystemExit(f"not found: {text}")
    return result["candidates"][0]


def cmd_stats(args):
    print(json.dumps(_adapter(args).stats(), indent=2, sort_keys=True))


def cmd_capabilities(args):
    """The contract view clients see, not just the static declaration."""
    print(json.dumps(api.capabilities(_adapter(args)), indent=2, sort_keys=True))


def cmd_resolve(args):
    adapter = _adapter(args)
    result = adapter.resolve({"text": args.query, "type": args.type}, {"taxon_id": 9606})
    if result["status"] != "ok":
        print(f"{result['status']}: {args.query}", file=sys.stderr)
        return 1
    for row in result["candidates"]:
        loc = row.get("location") or {}
        print(f"{row['type']:14s} {row['id']:36s} v={str(row['version']):>4s} "
              f"{loc.get('seq_region', '')}:{loc.get('start', '')}-{loc.get('end', '')} "
              f"{row['label'] or ''}  [{row['matched_by']}]")
    return 0


def cmd_transcripts(args):
    adapter = _adapter(args)
    gene = _first(adapter, args.gene, "Gene")
    rows = adapter.gene_transcripts(gene["stable_id"])
    print(f"{gene['stable_id']} {gene['label']} transcripts={len(rows)}")
    print(f"{'transcript':18s} {'ver':>4s} {'len':>6s} {'biotype':24s} canonical")
    for tx in rows:
        protein = adapter.transcript_translation(tx["stable_id"])
        print(f"{tx['stable_id']:18s} {str(tx['version']):>4s} "
              f"{str(protein['length'] if protein else '-'):>6s} {tx['biotype']:24s} "
              f"{'yes' if 'Ensembl_canonical' in tx['attributes'].get('tags', []) else ''}")
    return 0


def cmd_map(args):
    adapter = _adapter(args)
    if args.transcript:
        tx_id = _first(adapter, args.transcript, "Transcript")["stable_id"]
    else:
        gene = _first(adapter, args.gene, "Gene")
        tx_id = adapter.canonical_transcript(gene["stable_id"])["stable_id"]
    print(f"{'pos':>4s} {'aa':>2s} {'codon':6s} {'genomic':>19s} {'exon':>5s}")
    for pos in range(args.start, args.end + 1):
        row = adapter.map_protein_pos(tx_id, pos)
        if row["status"] != "ok":
            print(f"{pos:>4d} {row['status']}")
            continue
        genomic = f"{row['genomic_interval'][0]}-{row['genomic_interval'][1]}"
        print(f"{pos:>4d} {row.get('residue', '-'):>2s} {row['codon']:6s} "
              f"{genomic:>19s} {str(row['exon_ranks']):>5s}")
    return 0


def cmd_coverage(args):
    adapter = _adapter(args)
    if args.uniprot:
        result = adapter.coverage(args.uniprot, args.start, args.end)
    else:
        gene = _first(adapter, args.gene, "Gene")
        result = adapter.gene_coverage(gene["stable_id"], args.start, args.end)
    counts = {}
    for entry in result["entries"]:
        counts[entry["status"]] = counts.get(entry["status"], 0) + 1
    print(f"window {args.start}-{args.end} entries={len(result['entries'])} "
          f"{json.dumps(counts, sort_keys=True)}")
    for entry in result["entries"]:
        mark = "*" if entry.get("is_canonical") else " "
        print(f" {mark}{entry['transcript_id']:18s} v{str(entry.get('transcript_version')):>3s} "
              f"L={str(entry.get('protein_length')):>4s} {entry['status']}")
    return 0


def cmd_records(args):
    adapter = _adapter(args)
    fields = args.fields.split(",") if args.fields else None
    region = None
    if args.start is not None and args.end is not None:
        region = {"protein_range": [args.start, args.end]}
    result = adapter.get_records(args.ids.split(","), fields=fields, region=region)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "ok" else 1


def cmd_walk(args):
    adapter = _adapter(args)
    relations = args.relations.split(",") if args.relations else None
    result = adapter.walk(args.seed, relations=relations, max_steps=args.steps,
                          budget=args.budget, rng_seed=args.rng_seed)
    walk = result["walks"][0]
    print(f"seed {walk.get('seed_query')} -> {walk.get('seed')}  stopped={walk.get('stopped')} "
          f"visited={len(walk.get('visited', []))} {walk.get('visited_by_type')}")
    for step in walk.get("path", []):
        print(f"  {step['from']} --{step['relation']}--> {step['to']}")
    return 0 if result.get("status") == "ok" else 1


def cmd_neighborhood(args):
    adapter = _adapter(args)
    result = adapter.neighborhood(args.seed, depth=args.depth, budget=args.budget)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("status") == "ok" else 1


def cmd_region(args):
    adapter = _adapter(args)
    rows = adapter.region_genes(args.chrom, args.start, args.end)
    for gene in rows:
        print(f"{gene['stable_id']:18s} {gene['start']:>9d}-{gene['end']:<9d} "
              f"{gene['label']:16s} {gene['biotype']}")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=DEFAULT_DB)
    parser.add_argument("--genome", default=None)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("stats"); p.set_defaults(func=cmd_stats)
    p = sub.add_parser("capabilities"); p.set_defaults(func=cmd_capabilities)

    p = sub.add_parser("resolve")
    p.add_argument("query")
    p.add_argument("--type", default=None)
    p.set_defaults(func=cmd_resolve)

    p = sub.add_parser("transcripts"); p.add_argument("gene"); p.set_defaults(func=cmd_transcripts)

    p = sub.add_parser("map")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene"); group.add_argument("--transcript")
    p.add_argument("--start", type=int, required=True); p.add_argument("--end", type=int, required=True)
    p.set_defaults(func=cmd_map)

    p = sub.add_parser("coverage")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene"); group.add_argument("--uniprot")
    p.add_argument("--start", type=int, required=True); p.add_argument("--end", type=int, required=True)
    p.set_defaults(func=cmd_coverage)

    p = sub.add_parser("records")
    p.add_argument("--ids", required=True, help="comma-separated contract IDs")
    p.add_argument("--fields", default=None, help="comma-separated fields")
    p.add_argument("--start", type=int, default=None); p.add_argument("--end", type=int, default=None)
    p.set_defaults(func=cmd_records)

    p = sub.add_parser("walk")
    p.add_argument("--seed", default="4DNK"); p.add_argument("--steps", type=int, default=10)
    p.add_argument("--budget", type=int, default=30); p.add_argument("--rng-seed", type=int, default=20261106)
    p.add_argument("--relations", default=None, help="comma-separated relation whitelist")
    p.set_defaults(func=cmd_walk)

    p = sub.add_parser("neighborhood")
    p.add_argument("--seed", default="YWHAB"); p.add_argument("--depth", type=int, default=2)
    p.add_argument("--budget", type=int, default=40)
    p.set_defaults(func=cmd_neighborhood)

    p = sub.add_parser("region")
    p.add_argument("--chrom", default="20"); p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.set_defaults(func=cmd_region)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
