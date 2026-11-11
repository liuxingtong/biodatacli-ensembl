"""Contract tests: the rules the schema cannot express.

These run against the built ch20 index and exercise the same ``api`` layer the
MCP server exposes, so a regression here is a regression in what clients see.
"""

import json
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from ensembl_adapter import EnsemblAdapter  # noqa: E402
from ensembl_adapter import api  # noqa: E402

DB = os.path.join(ROOT, "data", "index", "ensembl_chr20_r116.sqlite")
GENE = "Ensembl:ENSG00000166913"
PROTEIN = "Ensembl:ENSP00000300161"


def collect_sources(value, found=None):
    if found is None:
        found = []
    if isinstance(value, list):
        for item in value:
            collect_sources(item, found)
    elif isinstance(value, dict):
        if {"db", "id", "version", "retrieved_at_status"} <= set(value):
            found.append(value)
        for item in value.values():
            collect_sources(item, found)
    return found


@unittest.skipUnless(os.path.exists(DB), "chr20 index not built")
class ContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.adapter = EnsemblAdapter(DB)

    def test_capabilities_match_the_index(self):
        caps = api.capabilities(self.adapter)
        self.assertTrue(caps["consistent"])
        self.assertEqual(caps["snapshot"]["release"], self.adapter.meta()["release"])
        self.assertTrue(caps["snapshot"]["build_id"])
        self.assertTrue(caps["snapshot"]["built_at"].endswith("Z"))
        self.assertEqual(set(caps["null_semantics"]),
                         {"not_found", "not_reported", "na", "access_failed"})

    def test_every_frozen_input_has_provenance(self):
        caps = api.capabilities(self.adapter)
        self.assertTrue(caps["provenance"])
        for key, entry in caps["provenance"].items():
            self.assertIn(entry["retrieved_at_status"], {"download", "file_mtime", "not_reported"}, key)
            self.assertRegex(entry["retrieved_at"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$", key)

    def test_facts_carry_retrieved_at(self):
        """The old index said retrieved_at=None for every fact; it must not come back."""
        result = api.get_records(self.adapter, [PROTEIN], ["sequence"],
                                 {"protein_range": [230, 246]})
        self.assertGreaterEqual(len(result["facts"]), 2)
        for fact in result["facts"]:
            self.assertIsNotNone(fact["source"]["retrieved_at"], fact["attribute"])
            self.assertEqual(fact["source"]["retrieved_at_status"], "file_mtime")
            self.assertEqual(fact["entity"]["namespace"], "Ensembl")
            self.assertIn(fact["evidence_type"], {"experimental", "computational",
                                                  "literature", "inferred"})
        attributes = {fact["attribute"] for fact in result["facts"]}
        self.assertEqual(attributes, {"protein_sequence", "protein_position_mapping"})

    def test_record_facts_carry_provenance_too(self):
        result = api.get_records(self.adapter, [GENE], ["xrefs"])
        facts = [fact for fact in result["facts"] if fact["attribute"] == "cross_reference"]
        self.assertTrue(facts)
        for fact in facts[:20]:
            self.assertTrue(fact["source"]["field"].startswith("BioMart"), fact["source"])
            self.assertIsNotNone(fact["source"]["retrieved_at"])

    def test_null_states_are_distinct(self):
        self.assertEqual(api.resolve(self.adapter, "ZZZ_NO_SUCH_GENE")["null_status"], "not_found")
        self.assertEqual(api.resolve(self.adapter, "YWHAB", {"taxon_id": 10090})["null_status"], "na")
        self.assertEqual(api.recall(self.adapter, "YWHAB", "structure")["null_status"], "na")
        self.assertIsNone(api.resolve(self.adapter, "YWHAB")["null_status"])

    def test_failed_query_is_not_an_empty_success(self):
        missing = api.resolve(self.adapter, "ZZZ_NO_SUCH_GENE")
        self.assertFalse(missing["ok"])
        self.assertEqual(missing["candidates"], [])
        ghost = api.neighbors(self.adapter, "Ensembl:ENSG00000000000")
        self.assertFalse(ghost["ok"])
        self.assertEqual(ghost["items"], [])

    def test_truncation_is_explicit(self):
        cut = api.neighbors(self.adapter, GENE, ["has_transcript"], "out", 1)
        self.assertTrue(cut["truncated"])
        self.assertEqual(cut["total"], 40)
        self.assertEqual(len(cut["items"]), 1)

    def test_invalid_region_is_reported_not_raised(self):
        result = api.get_records(self.adapter, [PROTEIN], region={"protein_range": [2, 1]})
        self.assertEqual(result["status"], "invalid_region")
        self.assertFalse(result["ok"])
        self.assertEqual(result["records"], [])

    def test_walk_replays(self):
        first = api.walk(self.adapter, [GENE], ["has_transcript"], 4, 5)
        second = api.walk(self.adapter, [GENE], ["has_transcript"], 4, 5)
        self.assertEqual(first, second)

    def test_sources_share_one_shape(self):
        result = api.neighbors(self.adapter, GENE, ["has_transcript"], "out", 3)
        for item in result["items"]:
            for source in collect_sources(item):
                self.assertIn("db", source)
                self.assertIn("version", source)
                self.assertIn("retrieved_at_status", source)

    def test_crosswalk_reaches_the_ensembl_protein(self):
        result = api.resolve(self.adapter, "P31946")
        reached = {node["id"]: node for node in result["candidates"]}
        self.assertIn("Ensembl:ENSP00000300161", reached)
        self.assertEqual(reached["Ensembl:ENSP00000300161"]["matched_by"],
                         "crosswalk:protein_uniprot")
        self.assertTrue(result["crosswalk_ambiguous"])

    def test_crosswalk_edges_carry_provenance(self):
        hits = self.adapter.crosswalk_lookup("4DNK")
        self.assertTrue(hits)
        for hit in hits:
            self.assertTrue(hit["evidence"]["file"])
            self.assertTrue(hit["evidence"]["sha256"])
            self.assertIn(hit["evidence"]["retrieved_at_status"], {"download", "file_mtime"})
            self.assertIn(hit["cardinality"], {"1:1", "1:n", "n:1", "n:m"})

    def test_crosswalk_stats_are_declared(self):
        stats = api.capabilities(self.adapter)["crosswalk"]
        self.assertTrue(stats["available"])
        self.assertGreater(stats["edges"], 0)
        self.assertTrue(stats["build_id"])
        self.assertIn("protein_uniprot", stats["by_relation"])

    def test_facts_declare_units(self):
        result = api.get_records(self.adapter, [PROTEIN], ["sequence"],
                                 {"protein_range": [230, 246]})
        units = {fact["attribute"]: fact["value"].get("unit") for fact in result["facts"]}
        self.assertEqual(units.get("protein_sequence"), "aa")
        self.assertEqual(units.get("protein_position_mapping"), "aa")

    def test_dataset_records(self):
        result = api.datasets(self.adapter)
        self.assertGreaterEqual(result["total"], 2)
        for record in result["datasets"]:
            self.assertTrue(record["dataset_id"])
            self.assertIsInstance(record["size"], int)
            self.assertIn(record["access"]["kind"], {"download", "build", "api"})
            self.assertTrue(record["provenance"]["db"])
            self.assertTrue(record["provenance"]["retrieved_at"])
        index = [r for r in result["datasets"] if r["access"]["kind"] == "build"]
        self.assertEqual(len(index), 1)
        self.assertEqual(index[0]["format"], "sqlite")

    def test_trace_records_calls(self):
        import tempfile
        path = os.path.join(tempfile.mkdtemp(), "trace.jsonl")
        os.environ["ENSEMBL_TRACE"] = path
        try:
            api.resolve(self.adapter, "YWHAB")
            api.neighbors(self.adapter, GENE, ["has_transcript"], "out", 2)
        finally:
            os.environ.pop("ENSEMBL_TRACE", None)
        lines = [json.loads(line) for line in open(path, encoding="utf-8")]
        self.assertEqual([entry["tool"] for entry in lines], ["resolve", "neighbors"])
        self.assertTrue(all(entry["build_id"] for entry in lines))
        self.assertTrue(all(entry["elapsed_ms"] >= 0 for entry in lines))

    def test_platform_tool_shapes(self):
        from ensembl_adapter import mcp_contract as platform
        described = platform.describe(self.adapter)
        self.assertTrue({"adapter_id", "entities", "scales", "keyspaces",
                         "evidence_types", "search_modes", "license"} <= set(described))

        resolved = platform.resolve(self.adapter, "P31946")
        first = resolved["candidates"][0]
        self.assertTrue({"canonical_id", "namespace", "score", "aliases"} <= set(first))
        self.assertTrue(resolved["ambiguous"])

        found = platform.search(self.adapter, "YWHAB", {"attribute": "gene_location"}, 5)
        self.assertTrue(found["facts"])
        self.assertTrue({"entity", "attribute", "value", "evidence_type",
                         "source", "score"} <= set(found["facts"][0]))

        facts = platform.get_facts(self.adapter, [PROTEIN], ["sequence"], "protein")
        self.assertTrue(facts["facts"])
        self.assertIn("missing", facts["coverage"])

    def test_access_plan_and_expand(self):
        from ensembl_adapter import mcp_contract as platform
        plan = platform.get_access_plan(self.adapter, "index:ensembl_chr20_r116.sqlite")
        self.assertEqual(plan["status"], "ok")
        self.assertTrue(plan["fetch_url_or_cmd"] and plan["format"] and plan["license"])
        self.assertIsInstance(plan["bytes"], int)

        expanded = platform.expand(self.adapter, [GENE], "has_transcript", 1, 3)
        self.assertEqual(expanded["status"], "ok")
        self.assertTrue({"from", "to", "relation", "evidence"} <= set(expanded["edges"][0]))

    def test_responses_validate_against_the_schema(self):
        from jsonschema import Draft202012Validator
        with open(os.path.join(ROOT, "contract", "m1.schema.json"), encoding="utf-8") as handle:
            root = json.load(handle)
        response = api.get_records(self.adapter, [PROTEIN], ["sequence"],
                                   {"protein_range": [230, 246]})
        validator = Draft202012Validator(
            {"$schema": root["$schema"], "$defs": root["$defs"], "$ref": "#/$defs/get_records"})
        errors = list(validator.iter_errors(response))
        self.assertEqual(errors, [], [error.message for error in errors])


if __name__ == "__main__":
    unittest.main()
