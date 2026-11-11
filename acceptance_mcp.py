"""Run the Ensembl-only golden queries through an actual MCP stdio session."""

import argparse
import asyncio
from contextlib import asynccontextmanager
import json
from pathlib import Path
import sys

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from mcp.client.streamable_http import streamable_http_client


@asynccontextmanager
async def session_streams(url):
    if url:
        async with streamable_http_client(url) as (reader, writer, _):
            yield reader, writer
        return
    server = StdioServerParameters(
        command=sys.executable, args=["-m", "ensembl_adapter.mcp_server"],
        cwd=str(Path(__file__).resolve().parent),
    )
    async with stdio_client(server) as streams:
        yield streams


async def acceptance(url=None):
    trace = []
    async with session_streams(url) as (reader, writer):
        async with ClientSession(reader, writer) as session:
            initialized = await session.initialize()
            tools = {tool.name: tool for tool in (await session.list_tools()).tools}
            assert set(tools) == {"capabilities", "resolve", "neighbors", "get_records", "walk", "recall"}
            assert "node_id" in tools["neighbors"].inputSchema["properties"]

            async def call(name, **arguments):
                result = await session.call_tool(name, arguments)
                assert not result.isError, (name, result.content)
                payload = result.structuredContent
                assert isinstance(payload, dict), name
                trace.append({"tool": name, "arguments": arguments, "response": payload})
                return payload

            scope = await call("capabilities")
            assert scope["scope"]["region"] == "20"
            gene, = (await call("resolve", ref={"text": "YWHAB", "type": "Gene"},
                                constraints={"taxon_id": 9606}))["candidates"]
            assert gene["id"] == "Ensembl:ENSG00000166913"
            neighbors = await call("neighbors", node_id=gene["id"], relations=["has_transcript"],
                                   direction="out", limit=50)
            assert neighbors["total"] == 40 and not neighbors["truncated"]
            transcript, = [item["node"] for item in neighbors["items"]
                           if "Ensembl_canonical" in item["node"]["attributes"].get("tags", [])]
            assert (transcript["id"], transcript["version"]) == ("Ensembl:ENST00000353703", 9)
            model, = (await call("get_records", ids=[transcript["id"]],
                                 fields=["exons", "translation"]))["records"]
            assert len(model["fields"]["exons"]) == 6
            protein_id = model["fields"]["translation"]["id"]
            assert protein_id == "Ensembl:ENSP00000300161"
            args = dict(ids=[protein_id], fields=["sequence", "protein_mapping"],
                        region={"protein_range": [230, 246]})
            response = await call("get_records", **args)
            assert response == await call("get_records", **args)
            protein, = response["records"]
            assert protein["fields"]["sequence"]["slice"] == "WTSENQGDEGDAGEGEN"
            assert protein["region"]["retrieved"] == list(range(230, 247))
            fact = next(f for f in response["facts"] if f["attribute"] == "protein_position_mapping")
            assert fact["value"]["data"]["codon"] == "TGG"
            assert fact["value"]["data"]["genomic_interval"] == [44906385, 44906387]
            assert fact["source"]["version"] == "116"
            assert fact["source"]["retrieved_at"].endswith("Z"), "facts must carry retrieval time"
            assert fact["source"]["retrieved_at_status"] in {"download", "file_mtime"}
            assert protein["source"]["entity_version"] == "4"
            assert protein["source"]["release"] == "116"
            assert (await call("resolve", ref="ZZZ_NO_SUCH_GENE"))["null_status"] == "not_found"
            assert (await call("resolve", ref="YWHAB",
                               constraints={"taxon_id": 10090}))["scope_status"] == "unsupported"
            assert (await call("recall", query="YWHAB", view="structure"))["null_status"] == "na"
            assert (await call("resolve", ref="ENST00000353703.999"))["null_status"] == "not_found"
            truncated = await call("neighbors", node_id=gene["id"],
                                   relations=["has_transcript"], direction="out", limit=1)
            assert truncated["truncated"] and truncated["total"] == 40
            partial, = (await call("get_records", ids=[protein_id], fields=["sequence"],
                                   region={"protein_range": [246, 247]}))["records"]
            assert partial["region"]["retrieved"] == [246]
            assert partial["region"]["missing"] == [{"position": 247, "reason": "beyond_protein"}]
            invalid = await call("get_records", ids=[protein_id], region={"protein_range": [2, 1]})
            assert invalid["status"] == "invalid_region"
            walk_args = dict(seeds=[gene["id"]], relations=["has_transcript"], budget=2)
            assert await call("walk", **walk_args) == await call("walk", **walk_args)
            return {"protocol": initialized.protocolVersion, "snapshot": scope["snapshot"],
                    "calls": trace}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="Optional JSON replay trace")
    parser.add_argument("--url", help="Test an already-running Streamable HTTP MCP endpoint")
    args = parser.parse_args()
    result = asyncio.run(acceptance(args.url))
    if args.out:
        args.out.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print(f"PASS: 6 MCP tools, {len(result['calls'])} calls, Ensembl-only YWHAB replay")
