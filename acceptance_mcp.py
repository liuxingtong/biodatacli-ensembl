"""End-to-end check through a real MCP session.

Asserts the platform tool surface (4 required + 4 optional) and the shape of
what each returns. Run against a local index:

    python3 acceptance_mcp.py --db data/index/ensembl_chr20_r116.sqlite
    python3 acceptance_mcp.py --url http://127.0.0.1:8000/mcp   # already running
"""

import argparse
import asyncio
from contextlib import asynccontextmanager
import json
from pathlib import Path
import sys

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.client.streamable_http import streamable_http_client

ROOT = Path(__file__).resolve().parent
TOOLS = {"describe", "resolve", "search", "get_facts",
         "search_datasets", "get_dataset_metadata", "get_access_plan", "expand"}


@asynccontextmanager
async def session_streams(url, db):
    if url:
        async with streamable_http_client(url) as (reader, writer, _):
            yield reader, writer
        return
    server = StdioServerParameters(
        command=sys.executable,
        args=["-m", "ensembl_adapter.mcp_server", "--db", str(db)],
        cwd=str(ROOT),
    )
    async with stdio_client(server) as streams:
        yield streams


async def acceptance(url=None, db=ROOT / "data/index/ensembl_chr20_r116.sqlite"):
    trace = []
    async with session_streams(url, db) as (reader, writer):
        async with ClientSession(reader, writer) as session:
            initialized = await session.initialize()
            tools = {tool.name for tool in (await session.list_tools()).tools}
            assert tools == TOOLS, sorted(tools ^ TOOLS)

            async def call(name, **arguments):
                result = await session.call_tool(name, arguments)
                assert not result.isError, (name, result.content)
                payload = result.structuredContent
                assert isinstance(payload, dict), name
                trace.append({"tool": name, "arguments": arguments, "response": payload})
                return payload

            scope = await call("describe")
            assert scope["adapter_id"].startswith("ensembl"), scope["adapter_id"]
            assert {"entities", "scales", "keyspaces", "evidence_types",
                    "search_modes", "license"} <= set(scope)

            resolved = await call("resolve", term="P31946")
            assert resolved["status"] == "ok" and len(resolved["candidates"]) >= 2
            first = resolved["candidates"][0]
            assert {"canonical_id", "namespace", "score", "aliases"} <= set(first)

            found = await call("search", query="YWHAB", filters={"attribute": "gene_location"}, k=3)
            assert found["facts"], found
            assert {"entity", "attribute", "value", "evidence_type", "source", "score"} \
                <= set(found["facts"][0])

            facts = await call("get_facts", entity_ids=["Ensembl:ENSP00000300161"],
                               fields=["sequence"], scale="protein")
            assert facts["facts"] and "missing" in facts["coverage"]

            catalogue = await call("search_datasets", query="")
            assert catalogue["total"] >= 2
            dataset_id = catalogue["datasets"][0]["dataset_id"]

            metadata = await call("get_dataset_metadata", dataset_id=dataset_id)
            assert metadata["dataset"]["provenance"]["db"] == "Ensembl"

            plan = await call("get_access_plan", dataset_id=dataset_id)
            assert {"fetch_url_or_cmd", "bytes", "format", "license"} <= set(plan)

            expanded = await call("expand", entity_ids=["Ensembl:ENSG00000166913"],
                                  key="has_transcript", max_depth=1, k=5)
            assert expanded["edges"]
            assert {"from", "to", "relation", "evidence"} <= set(expanded["edges"][0])

            return {"protocol": initialized.protocolVersion, "snapshot": scope.get("snapshot"),
                    "calls": trace}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", type=Path, default=ROOT / "data/index/ensembl_chr20_r116.sqlite")
    parser.add_argument("--url", help="Test an already-running streamable HTTP endpoint")
    parser.add_argument("--out", type=Path, help="Optional JSON replay trace")
    args = parser.parse_args()
    result = asyncio.run(acceptance(args.url, args.db))
    if args.out:
        args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"PASS: 8 platform tools, {len(result['calls'])} calls")
