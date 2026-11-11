"""MCP surface for the Ensembl M1 slice.

The tools here are deliberately thin: argument validation lives in this file,
while the contract shape (envelope, null semantics, facts, provenance) lives in
:mod:`ensembl_adapter.api`.  The conformance runner talks to the same module, so
what passes locally is what an MCP session returns.
"""

import argparse
from contextlib import closing, contextmanager
import json
from pathlib import Path
from typing import Annotated, Any, Literal

from mcp.server.fastmcp import FastMCP
from pydantic import Field

from . import api
from .adapter import EnsemblAdapter

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data/index/ensembl_chr20_r116.sqlite"
CONTRACT = api.CONTRACT_VERSION
Limit = Annotated[int, Field(ge=1, le=1000)]


def create_server(db_path=DEFAULT_DB, genome_path=None, *, transport="stdio", port=8000):
    mcp = FastMCP("ensembl-chr20", host="127.0.0.1", port=port,
                  stateless_http=True, json_response=True)
    declaration = json.loads((ROOT / "capabilities.json").read_text())

    # One request-scoped read-only connection; refuse an index that drifted
    # away from the declared snapshot.
    @contextmanager
    def adapter():
        instance = EnsemblAdapter(str(db_path), genome_path)
        with closing(instance.conn):
            meta = instance.meta()
            if (meta.get("release"), meta.get("assembly"), meta.get("chromosome")) != (
                    declaration["release"], declaration["assembly"], declaration["scope"]["region"]):
                raise ValueError("Index does not match the declared Ensembl snapshot")
            yield instance

    @mcp.tool(structured_output=True)
    def capabilities() -> dict[str, Any]:
        """Describe this service's coverage, provenance and null semantics."""
        with adapter() as a:
            return {**api.capabilities(a), "tool_contract": CONTRACT, "transport": transport}

    @mcp.tool(structured_output=True)
    def resolve(ref: str | dict, constraints: dict | None = None) -> dict[str, Any]:
        """Resolve a name or native ID within the declared Ensembl slice."""
        with adapter() as a:
            return api.resolve(a, ref, constraints)

    @mcp.tool(structured_output=True)
    def neighbors(node_id: str, relations: list[str] | None = None,
                  direction: Literal["in", "out", "both"] = "both",
                  limit: Limit = 20) -> dict[str, Any]:
        """Look up typed adjacent objects and the source of each edge."""
        with adapter() as a:
            return api.neighbors(a, node_id, relations, direction, limit)

    @mcp.tool(structured_output=True)
    def get_records(ids: Annotated[list[str], Field(min_length=1, max_length=50)],
                    fields: list[str] | None = None,
                    region: dict | None = None) -> dict[str, Any]:
        """Read records and emit self-sourced facts for what was requested."""
        with adapter() as a:
            return api.get_records(a, ids, fields, region)

    @mcp.tool(structured_output=True)
    def walk(seeds: Annotated[list[str], Field(min_length=1, max_length=20)],
             relations: list[str] | None = None,
             max_steps: Annotated[int, Field(ge=0, le=100)] = 6,
             budget: Limit = 30, rng_seed: int = 20261106) -> dict[str, Any]:
        """Walk typed relations with deterministic sampling and a per-seed budget."""
        with adapter() as a:
            return api.walk(a, seeds, relations, max_steps, budget, rng_seed)

    @mcp.tool(structured_output=True)
    def recall(query: str, view: str = "exact", limit: Limit = 20) -> dict[str, Any]:
        """Recall exact or text candidates in the local index."""
        with adapter() as a:
            return api.recall(a, query, view, limit)

    return mcp


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=DEFAULT_DB)
    parser.add_argument("--genome", type=Path)
    parser.add_argument("--transport", choices=["stdio", "streamable-http"], default="stdio")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    create_server(args.db, args.genome, transport=args.transport,
                  port=args.port).run(transport=args.transport)
