"""MCP server exposing the platform's required and optional tools.

Required: describe / resolve / search / get_facts
Optional: search_datasets / get_dataset_metadata / get_access_plan / expand

All business logic lives in :mod:`ensembl_adapter.mcp_contract`; this file only
adds transport and argument validation.
"""

import argparse
from contextlib import closing, contextmanager
import json
from pathlib import Path
from typing import Annotated, Any

from mcp.server.fastmcp import FastMCP
from pydantic import Field

from . import api, mcp_contract as platform
from .adapter import EnsemblAdapter

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data/index/ensembl_chr20_r116.sqlite"
Limit = Annotated[int, Field(ge=1, le=200)]


def create_server(db_path=DEFAULT_DB, genome_path=None, *, transport="stdio", port=8000):
    mcp = FastMCP("ensembl", host="127.0.0.1", port=port,
                  stateless_http=True, json_response=True)
    declaration = json.loads((ROOT / "capabilities.json").read_text())

    @contextmanager
    def adapter():
        instance = EnsemblAdapter(str(db_path), genome_path)
        with closing(instance.conn):
            meta = instance.meta()
            if (meta.get("release"), meta.get("assembly")) != (
                    declaration["release"], declaration["assembly"]):
                raise ValueError("Index does not match the declared Ensembl snapshot")
            yield instance

    @mcp.tool(structured_output=True)
    def describe() -> dict[str, Any]:
        """Capability declaration: entities, scales, keyspaces, evidence types."""
        with adapter() as a:
            return platform.describe(a)

    @mcp.tool(structured_output=True)
    def resolve(term: str, context: dict | None = None,
                species: str | None = "human") -> dict[str, Any]:
        """Free text or alias to canonical key."""
        with adapter() as a:
            return platform.resolve(a, term, context, species)

    @mcp.tool(structured_output=True)
    def search(query: str, filters: dict | None = None, k: Limit = 10) -> dict[str, Any]:
        """Retrieve facts matching a query, with evidence and source."""
        with adapter() as a:
            return platform.search(a, query, filters, k)

    @mcp.tool(structured_output=True)
    def get_facts(entity_ids: list[str], fields: list[str] | None = None,
                  scale: str | None = None) -> dict[str, Any]:
        """Structured facts for known ids, plus coverage of the request."""
        with adapter() as a:
            return platform.get_facts(a, entity_ids, fields, scale)

    @mcp.tool(structured_output=True)
    def search_datasets(query: str = "", filters: dict | None = None,
                        k: Limit = 10) -> dict[str, Any]:
        """Dataset catalogue search."""
        with adapter() as a:
            return platform.search_datasets(a, query, filters, k)

    @mcp.tool(structured_output=True)
    def get_dataset_metadata(dataset_id: str) -> dict[str, Any]:
        """Full metadata for one dataset."""
        with adapter() as a:
            return platform.get_dataset_metadata(a, dataset_id)

    @mcp.tool(structured_output=True)
    def get_access_plan(dataset_id: str, target: str = "on_demand") -> dict[str, Any]:
        """How to fetch a dataset on demand."""
        with adapter() as a:
            return platform.get_access_plan(a, dataset_id, target)

    @mcp.tool(structured_output=True)
    def expand(entity_ids: list[str], key: list[str] | str,
               max_depth: Annotated[int, Field(ge=1, le=3)] = 1,
               k: Limit = 50) -> dict[str, Any]:
        """Graph expansion across scales."""
        with adapter() as a:
            return platform.expand(a, entity_ids, key, max_depth, k)

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
