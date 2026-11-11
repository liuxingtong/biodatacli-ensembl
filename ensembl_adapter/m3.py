"""Local stand-in for the M3 seed step.

M3 (the future natural-language layer) is simulated here by a deterministic
seed resolver: the caller supplies a symbol / stable ID / UniProt / PDB ID,
the adapter resolves it to one graph node, and a seeded random walk explores
the neighbourhood.  Two runs with the same seed must produce byte-identical
JSON.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ensembl_adapter import EnsemblAdapter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", default=os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite"))
    parser.add_argument("--seed", default="4DNK", help="symbol, ENSG/ENST/ENSP, UniProt or PDB ID")
    parser.add_argument("--steps", type=int, default=10)
    parser.add_argument("--budget", type=int, default=30)
    parser.add_argument("--rng-seed", type=int, default=20261106)
    parser.add_argument("--relations", default=None, help="comma-separated relation whitelist")
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)

    adapter = EnsemblAdapter(args.db)
    relations = args.relations.split(",") if args.relations else None
    first = adapter.walk(args.seed, relations=relations, max_steps=args.steps,
                         budget=args.budget, rng_seed=args.rng_seed)
    second = adapter.walk(args.seed, relations=relations, max_steps=args.steps,
                          budget=args.budget, rng_seed=args.rng_seed)
    reproducible = first == second
    walk = first["walks"][0]

    result = {
        "release": adapter.meta()["release"],
        "assembly": adapter.meta()["assembly"],
        "simulated_m3": True,
        "seed_query": args.seed,
        "rng_seed": args.rng_seed,
        "steps": args.steps,
        "budget": args.budget,
        "reproducible": reproducible,
        "relations": relations,
        "walk": first,
    }
    out = args.out or os.path.join(ROOT, "runs", f"m3_{args.seed.replace(':', '_')}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, ensure_ascii=False)

    print(f"seed {args.seed!r} -> {walk.get('seed', 'not found')}")
    print(f"stopped={walk.get('stopped')} visited={len(walk.get('visited', []))} "
          f"by_type={walk.get('visited_by_type')}")
    for step in walk.get("path", []):
        print(f"  {step['from']} --{step['relation']}--> {step['to']}")
    print(f"reproducible={reproducible} out={out}")
    return 0 if reproducible else 1


if __name__ == "__main__":
    raise SystemExit(main())
