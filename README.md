# ensembl_proto — offline Ensembl M1 adapter

A dependency-free Python implementation of the SciDataGym **M1 single-database
adapter** for Ensembl release 116 (GRCh38.p14, human chromosome 20).

It builds a deterministic SQLite graph from release-pinned GFF3/FASTA/BioMart
files, keeps per-node and per-edge provenance, and exposes the common M1 calls:

```
resolve(ref, constraints)
neighbors(node_id, relations, direction, limit)
get_records(ids, fields, region)
walk(seeds, relations, max_steps, budget)
recall(query, view, limit)
capabilities()
```

The familiarisation note is at [`../ENSEMBL_RESEARCH_2026-11-06.md`](../ENSEMBL_RESEARCH_2026-11-06.md).

No third-party packages are required (Python 3.9+ standard library only).

## Quick start

```bash
cd ensembl_proto

# 1. download the pinned inputs (resume-safe; existing files are verified)
python3 fetch_data.py

# 2. derive the 4DNK/P31946 chain + range enrichment from the frozen SIFTS file,
#    then build the cross-database key table
python3 make_enrichment.py
python3 crosswalk/build.py

# 3. build the chr20 index
python3 -m ensembl_adapter.build \
  --gff3 data/raw/Homo_sapiens.GRCh38.116.chromosome.20.gff3.gz \
  --genome data/raw/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz \
  --pep data/raw/Homo_sapiens.GRCh38.pep.all.fa.gz \
  --biomart data/raw/biomart \
  --enrichment data/raw/xref_enrichment.json \
  --capabilities capabilities.json \
  --db data/index/ensembl_chr20_r116.sqlite \
  --stats data/index/build_stats.json

# 4. minimal prototype and M1 acceptance
python3 prototype_4dnk.py
python3 acceptance_m1.py

# 5. contract conformance, benchmark, offline checks and tests
python3 contract/conformance.py
python3 bench/run.py
python3 verify.py
python3 -m unittest discover -s tests -v
```

Or run everything:

```bash
./run_all.sh
```

## M1 contract

Contracts IDs are cross-database strings:

```
Ensembl:ENSG00000166913   Gene
Ensembl:ENST00000353703   Transcript
Ensembl:ENSE00001458785   Exon
Ensembl:ENSP00000300161   Protein (Ensembl translation)
UniProtKB:P31946          Protein (external reference)
PDB:4DNK                  StructureEntry (external reference)
```

Example: the PDF-style P31946 → 4DNK walk.

```python
from ensembl_adapter import EnsemblAdapter

a = EnsemblAdapter("data/index/ensembl_chr20_r116.sqlite")

# 1. resolve the external accession
a.resolve({"text": "P31946", "type": "Protein"}, {"taxon_id": 9606})
# -> candidates: [UniProtKB:P31946]

# 2. cross into Ensembl and then to the structure
a.neighbors("UniProtKB:P31946", relations=["uniprot_xref"], direction="in", limit=20)
a.neighbors("Ensembl:ENSP00000300161", relations=["pdb_xref"], direction="out", limit=50)
```

One neighbour item has the PDF's `node / edge / source` shape:

```json
{
  "node": {"id": "PDB:4DNK", "type": "StructureEntry"},
  "edge": {
    "from": "Ensembl:ENSP00000300161",
    "relation": "pdb_xref",
    "to": "PDB:4DNK",
    "chains": ["A", "B"],
    "protein_range": [1, 246],
    "observed_window": {"A": [230, 231, 232, 233], "B": [230, 231, 232, 233, 234]}
  },
  "source": {
    "database": "Ensembl",
    "record_id": "ENSP00000300161",
    "release": "116",
    "field": "BioMart pdb",
    "locator": "biomart/pdb.tsv:row 9524"
  }
}
```

`get_records` returns sequence and per-residue mapping, including what was
retrieved and what is missing:

```python
a.get_records(
    ["Ensembl:ENSP00000300161"],
    fields=["sequence", "protein_mapping", "xrefs"],
    region={"protein_range": [230, 246]},
)
# region.retrieved = 230..246 (17/17)
# region.missing   = []
# fields.sequence.slice = "WTSENQGDEGDAGEGEN"
```

`walk` takes a relation whitelist, samples uniformly within the allowed
transitions, and reports visited nodes, path and stop reason:

```python
a.walk(
    ["UniProtKB:P31946"],
    relations=["uniprot_xref", "pdb_xref", "has_translation", "has_exon"],
    max_steps=8,
    budget=10,
    rng_seed=20261106,
)
```

## MCP tools (platform contract)

The server exposes exactly the tools the platform requires, plus the optional
ones:

| tool | purpose |
|---|---|
| `describe` | adapter_id, entities, scales, keyspaces, evidence_types, search_modes, license |
| `resolve` | free text / alias → `{canonical_id, namespace, score, aliases[]}` |
| `search` | query → facts (`entity, attribute, value, evidence_type, source, score`) |
| `get_facts` | ids → `{facts[], coverage{requested, entities_with_facts, missing}}` |
| `search_datasets` | dataset catalogue (optional) |
| `get_dataset_metadata` | full dataset record (optional) |
| `get_access_plan` | `{fetch_url_or_cmd, bytes, format, license}` (optional) |
| `expand` | `from/to/relation/evidence` edges across scales (optional) |

`ensembl_adapter/mcp_contract.py` holds the logic; `mcp_server.py` only adds
transport. The internal API (`capabilities / neighbors / get_records / walk /
recall / datasets`) stays available for the bench, bridge and console.

## Contract layer (m1.1)

`adapter.py` answers questions about the graph; `api.py` turns those answers
into the shipped contract. Both the MCP server and the conformance runner call
`api.py`, so a check that passes here describes what a client actually receives.

Every response carries the same envelope:

```json
{"ok": true, "status": "ok", "null_status": null, "truncated": false}
```

`null_status` is the only way to say "I do not know", and it has four values:

| null_status | meaning | example |
|---|---|---|
| `not_found` | the query ran, nothing in this slice matched | `resolve("ZZZ_NO_SUCH_GENE")` |
| `not_reported` | the entity exists, the source asserts no such value | a fact with no retrieval time |
| `na` | outside the declared scope | mouse `taxon_id=10090`, `recall(view="structure")` |
| `access_failed` | the underlying read failed; nothing is claimed | unreadable index |

A missing answer is never an empty success: `ok` is false whenever the question
was not answered, and `candidates`/`items` being empty never implies "no such
thing exists".

Provenance is attached to every fact and to every node/edge source block:

```json
{
  "entity": {"id": "ENSP00000300161", "namespace": "Ensembl"},
  "attribute": "protein_position_mapping",
  "value": {"type": "object", "data": {"protein_pos": 230, "codon": "TGG", "...": "..."}},
  "evidence_type": "inferred",
  "source": {
    "db": "Ensembl", "id": "ENSP00000300161", "version": "116",
    "entity_version": "4", "field": "derived: CDS + genome",
    "retrieved_at": "2026-11-06T12:22:48Z",
    "retrieved_at_status": "file_mtime",
    "assembly": "GRCh38.p14", "build_id": "9b6e009147221c77"
  },
  "scale": "protein"
}
```

`retrieved_at_status` says how the timestamp was obtained, highest trust first:
`download` (recorded by `fetch_data.py` when the bytes arrived), `file_mtime`
(the bytes' own mtime, used for files that predate per-file stamping), and
`not_reported` (unknown — never a fabricated date).

Values that have a unit declare it (`protein_sequence` in `aa`, `gene_location`
in `bp`), matching the `value{type,unit?}` shape of the contract.

`datasets()` returns a `DatasetRecord` for every frozen input plus the index
itself: `dataset_id`, `name`, `description`, `size`, `format`,
`access{kind,fetch,read_only}` and `provenance{db,release,sha256,retrieved_at}`.

Every contract call can be traced by pointing `ENSEMBL_TRACE` at a file; each
line records tool, arguments, status, item count, elapsed time and `build_id`.

Facts are emitted for `protein_sequence`, `protein_position_mapping`,
`transcript_sequence`, `exon_structure`, `translation`, `gene_location` and
`cross_reference`; `evidence_type` is `computational` for values the source
asserts and `inferred` for values this adapter derives.

`capabilities()` reports what can be asked, plus `snapshot`
(`release`/`assembly`/`build_id`/`built_at`) and `provenance` (per input file).
A `consistent: false` there means the index drifted away from
`capabilities.json`.

Machine-readable shapes live in [`contract/m1.schema.json`](contract/m1.schema.json)
(JSON Schema 2020-12); the runnable checks live in
[`contract/conformance.py`](contract/conformance.py):

```bash
python3 contract/conformance.py            # 70 checks, exit 1 on failure
python3 contract/conformance.py --out runs/conformance.json
```

The suite covers what a schema cannot: exact 230–246 recovery, truncation
accounting, replay determinism, the four null states, timestamps present and
predating the build, and every response validating against the schema.

## Cross-database keys (crosswalk)

`crosswalk/crosswalk.jsonl.gz` is the key-alignment table: 20,087 edges over
nine namespaces (Ensembl, UniProtKB, PDB, HGNC, EntrezGene, MIM, CCDS), built
from the frozen BioMart exports plus SIFTS. Each edge carries the row that
asserted it, the file's SHA256, the release and when it was retrieved.

```
protein_uniprot        ENSP ⇄ UniProtKB     1,235 edges
protein_pdb            ENSP ⇄ PDB          13,444
gene_go                ENSG ⇄ GO           9,316
gene_reactome          ENSG ⇄ Reactome     3,604
transcript_ccds        ENST ⇄ CCDS          3,135
gene_hgnc              ENSG ⇄ HGNC          1,066
gene_entrez            ENSG ⇄ EntrezGene      748
gene_kegg_pathway      ENSG ⇄ KEGG            516
gene_mim               ENSG ⇄ MIM             458
gene_ec                ENSG ⇄ EC               86
pdb_uniprot_sifts      PDB  ⇄ UniProtKB         1
```

KEGG and EC edges are bridged through the Entrez ID (KEGG's gene key is the
Entrez number), and each edge records that bridge in `via`.

Authoritative keys covered: ENSG, UniProtKB, PDB, GO, Reactome, KEGG, EC (7 of
the 9 named in the spec). CL and UBERON are not covered: they need ontology
files plus expression-based mapping, not a BioMart join.

Fan-out is data, not a bug: 1,012 of the UniProt edges share their accession
with another chromosome-20 protein, and 2,763 PDB entries map to more than one
protein. `resolve` therefore reports `crosswalk_ambiguous` and gives each hit
an `alternatives` count instead of silently answering with the first one:

```
resolve("P31946")
  Ensembl:ENSP00000300161  matched_by=crosswalk:protein_uniprot  alternatives=1
  Ensembl:ENSP00000361930  matched_by=crosswalk:protein_uniprot  alternatives=1
  UniProtKB:P31946         matched_by=stable_id
  crosswalk_ambiguous=true
```

The table is declared in `capabilities().crosswalk` with its edge count and
`build_id`, so a client can tell whether an answer came from the graph, from
the crosswalk, or from both.

## Bridging into SciDataGym-Demo (M1 shim)

`bridge/scidatagym_m1.py` exposes this slice through the M1 interface the demo's
M3 engine already calls. No changes on either side — it is shape translation
only.

| demo's M1 | here | note |
|---|---|---|
| `describe()` | `api.capabilities()` | adds the `views` / `catalog` / `manifest` that M3 reads |
| `resolve(ref, constraints)` | `api.resolve()` | `ok/not_found` → `resolved/ambiguous/not_found` |
| `recall(query, view, limit)` | `api.recall()` | only `identifier` and `field`; `semantic` is rejected, not faked |
| `neighbors(id, rel, dir, limit, cursor)` | `api.neighbors()` | `cursor` slicing lives in the shim |
| `get_records(ids, fields, region{start,end})` | `api.get_records()` | `region` ↔ `protein_range` |
| `raw_record(accession)` | — | returns `not_supported`; the slice has no one-record-per-entity source |
| `walk(..., random_seed)` | `api.walk(rng_seed)` | same determinism |

```bash
python3 bridge/smoke.py     # 23 checks, organised by the four engine.py call sites
```

The smoke run mirrors real call sites (`engine.py:119/128/195/315`): ambiguous
`resolve("P31946")`, cursor paging over 40 transcripts, a 230–246 region, and a
missing id. It also asserts the negative paths — unsupported relation, unknown
node, unknown view.

**What does not transfer:** the demo's catalog is `Protein / GO / Keyword`, this
slice is `Gene / Transcript / Exon / Protein`; the relation vocabularies differ.
The mechanism plugs in; the C1–C3 questions themselves are UniProt questions and
would need Ensembl equivalents.

## Whole-genome build (all 25 chromosomes)

```bash
python3 full/fetch_parallel.py    # 1.0 GB from Ensembl, 24 parallel range streams
python3 full/split.py             # GFF3/FASTA -> one file per chromosome
python3 full/build_all.py         # 25 builds, one SQLite each
python3 full/merge.py             # renumber nodes, dedupe shared keys
```

Why not one build: `build()` holds the whole genome (3.1 GB of sequence plus every
parsed feature) in memory; 16 GB is not enough. Splitting keeps each build at a
few hundred MB, and `build()` itself is unchanged.

Measured on this machine (16 GB RAM, 8 cores):

| stage | time | output |
|---|---|---|
| download 1.0 GB | ~55 min (10–15 KB/s per stream, 24 streams ≈ 300 KB/s) | `full/raw/` |
| split 10.7M GFF lines + 3.09 Gb sequence | 7m11s | `full/split/` |
| 25 per-chromosome builds | 7m27s | `full/parts/` (3.2 GB) |
| merge | 2m36s | `full/index/ensembl_grch38_r116.sqlite` (3.1 GB) |

Merged graph: **2,622,629 nodes / 6,782,632 edges, 3.26 GB** — 78,733 genes,
644,292 transcripts, 370,886 proteins, 1,319,589 exons, 209,103 external
reference nodes. 27,467 duplicate keys collapsed when merging (shared assembly
node plus every PDB/UniProt entry referenced from more than one chromosome):
46,147 UniProt edges and 359,037 PDB edges survive the merge. Counts are the GFF3's own feature counts:
chr20 alone has 16,324 distinct transcript features, all with unique IDs.

The whole-genome BioMart exports (`full/fetch_biomart.py`) feed the rebuild, so
UniProt/PDB xref edges are present. GO / Reactome / KEGG / EC are not graph
edges; they come from the crosswalk, and the whole-genome one is built next to
the index:

```bash
python3 crosswalk/build.py --biomart full/biomart   --out full/index/crosswalk.jsonl.gz --stats full/index/crosswalk_stats.json
```

1,123,188 edges (GO 289k, PDB 394k, Reactome 146k, CCDS 115k, UniProt 51k, KEGG
20k, EC 3.9k). The adapter prefers a crosswalk that sits next to the database.
Loading it costs ~35 s and ~1.1M in-memory edges; candidate lists are capped at
50 with `truncated: true`.

## Visual console

```bash
python3 viz/server.py --port 8765
open http://127.0.0.1:8765
```

A read-only local console over the same `api.py` the MCP server exposes. Four
questions, four tabs:

| tab | what it shows | why it matters |
|---|---|---|
| 图谱探索 | neighbourhood graph, node card, provenance block, facts table | 看得见一条结论是怎么被 3–4 条边和它自己的 source 撑起来的 |
| 序列与坐标 | full protein with the requested window highlighted, per-residue codon and genomic interval | 坐标映射不是口头声明，逐残基可核 |
| 跨库对齐 | crosswalk edges, cardinality, alternatives, evidence file/row/sha256 | 一个 accession 对应多个对象时，歧义当场暴露 |
| 数据库体检 | verdict, snapshot identity, input provenance, coverage rates, CDS→peptide consistency, 24 reconciliation checks, 47 contract checks by group, graph scale, declared limits | 判断的是**数据本身**够不够格，而不是某次查询答得漂不漂亮 |

Deep links: `?q=YWHAB` opens a search, `#health` / `#seq` / `#cross` / `#eval`
open a tab directly. Screenshots of each view are in `viz/screenshots/`.

## Benchmark: 准 / 诚 / 寸

`bench/golden.jsonl` holds 32 queries in six categories (单实体, 跨尺度,
数据集, 模糊, 不存在, 陷阱). Each item names the frozen file that justifies its
answer (`source_of_truth`), so the gold is not a copy of adapter output:

| axis | question it asks | how it is scored |
|---|---|---|
| 准 | is the value right? | every check on the payload holds |
| 诚 | is it sourced? | expected IDs appear and a `retrieved_at` is present |
| 寸 | is the uncertainty right? | `null_status` equals the expected state; facts carry a plausible `evidence_type` |

```bash
python3 bench/run.py            # 32/32 准, 32/32 诚, 13/13 寸 (scored subset)
python3 bench/run.py --out runs/bench.json
```

The 寸 axis only applies to items that return facts or must return null (13 of
32); the rest are marked n/a rather than counted as passes.

**Scope of this suite.** Every item is derived from the frozen inputs and
carries its `source_of_truth`; there is no human annotation. It is a regression
suite with a fixed rubric (准/诚/寸): it catches breakage and keeps the scoring
convention stable. It is not an independent accuracy claim — the same snapshot
writes the questions and answers them, so a green run proves internal
consistency, nothing more.

## What is in the graph

```
assembly --has_region--> region --contains--> gene
gene --has_transcript--> transcript --has_exon--> exon (rank)
transcript --has_translation--> protein
gene / transcript / protein --<typed xref>--> external node
```

Typed relations:

```
has_region, contains, has_transcript, has_exon, has_translation,
uniprot_xref, pdb_xref, hgnc_xref, entrez_xref, mim_xref, ccds_xref
```

Every node and edge stores `source_db`, `source_record_id`, `source_version`,
`source_field` and `source_locator` (GFF3 line, FASTA header line, or BioMart
row). The extra tables (`cds_segment`, `exon_instance`, `seq`,
`transcript_stats`) back `get_records`; they are not part of the graph
contract.

Chr20 build: **67,001 nodes / 173,437 edges / 91 MiB**.

## CLI

```bash
python3 -m ensembl_adapter.cli stats
python3 -m ensembl_adapter.cli capabilities
python3 -m ensembl_adapter.cli resolve YWHAB --type Gene
python3 -m ensembl_adapter.cli transcripts YWHAB
python3 -m ensembl_adapter.cli map --gene YWHAB --start 230 --end 246
python3 -m ensembl_adapter.cli coverage --gene YWHAB --start 230 --end 246
python3 -m ensembl_adapter.cli records --ids Ensembl:ENSP00000300161 \
  --fields sequence,protein_mapping --start 230 --end 246
python3 -m ensembl_adapter.cli walk --seed 4DNK --steps 10 --budget 30 \
  --relations pdb_xref,has_translation,has_exon
```

## Layout

| Path | Purpose |
|---|---|
| `fetch_data.py` | fixed-URL downloader + SHA256 manifest |
| `fetch_rest.py` | freezes the REST responses used by verification |
| `make_enrichment.py` | derives 4DNK chains/range from frozen SIFTS |
| `capabilities.json` | M1 capability declaration (scope, relations, fields) |
| `ensembl_adapter/gff.py` | GFF3 reader with line numbers |
| `ensembl_adapter/fasta.py` | genome/peptide FASTA readers with header lines |
| `ensembl_adapter/sequence.py` | splice, reverse-complement, translate |
| `ensembl_adapter/build.py` | SQLite M1 graph builder |
| `ensembl_adapter/adapter.py` | graph queries + gene-model helpers |
| `ensembl_adapter/api.py` | M1 contract layer: envelope, null states, facts, provenance |
| `ensembl_adapter/provenance.py` | per-input `retrieved_at` and build stamping |
| `contract/m1.schema.json` | JSON Schema for the six responses |
| `contract/conformance.py` | 47-check contract suite |
| `crosswalk/build.py` | builds the cross-database key table |
| `ensembl_adapter/crosswalk.py` | crosswalk loader + ambiguity-aware lookup |
| `bench/golden.jsonl` | 32 golden queries, six categories |
| `bench/run.py` | 准/诚/寸 scorer |
| `viz/server.py` | read-only HTTP console (stdlib only) |
| `viz/web/` | single-page UI, no build step |
| `viz/screenshots/` | rendered views of the console |
| `bridge/scidatagym_m1.py` | M1 shim for SciDataGym-Demo's M3 |
| `bridge/smoke.py` | 23 checks mirroring engine.py call sites |
| `ensembl_adapter/crossref.py` | mmCIF CA + SIFTS UniProt mapping |
| `ensembl_adapter/m3.py` | deterministic M3-seed simulation |
| `verify.py` | 24 offline checks |
| `acceptance_m1.py` | 10-check P31946 → Ensembl → 4DNK acceptance |
| `prototype_4dnk.py` | 230–246 residue table and 4DNK overlay |
| `snapshots/rest/` | frozen REST JSON (release 116) |
| `runs/` | prototype/drift/verification/acceptance JSON |
| `tests/` | synthetic, real-index and contract tests |

## Version pinning

* Ensembl FTP: `release-116`
* Assembly: `GRCh38.p14` / `GCA_000001405.29`
* BioMart: `https://jun2026.archive.ensembl.org/biomart/martservice`
* REST snapshot: release 116 (`rest.ensembl.org` at snapshot time)

The graph builder never calls the network. Re-running `fetch_data.py` after a
new Ensembl release will **not** silently pick it up: URLs are fixed to
release 116. To move to another release, edit `fetch_data.py`, rebuild, and
refresh `fetch_rest.py` snapshots.

## Determinism

Three independent builds produce byte-identical SQLite files and the same
logical fingerprint:

```
SQLite SHA256  3d3df60344376aa21ec024fbeb37c750dcb35a2d270f51f7e3fea5bc87c7882c
fingerprint    883f71c192d332e8af5ad5c9b3fe8c0e0936faf3759b406a2249a51d6255a53d
```

See `runs/reproducibility.json`.

## Out of scope in this slice

* only human chromosome 20 (declared in `capabilities.json`)
* no patch / alt / PAR_Y / variation / regulation / compara nodes
* PDB chains and range are enriched only for 4DNK
* structure coordinates are read by a PDB adapter, not this one

## Ensembl MCP 本地试运行

MCP 进程复用现有 SQLite 和 chr20 基因组 FASTA，原有 Python Adapter 可继续独立使用。本轮用 MCP Python SDK 1.26.0 验证 stdio 和本机 Streamable HTTP。六个工具的工具名和参数来自现有 M1 入口，合同标识为 `ensembl-mcp-draft-v1`。

```sh
uv run --python 3.12 --with mcp==1.26.0 python acceptance_mcp.py --out runs/acceptance_mcp.json
```

平台本地对接时也可以启动 HTTP 进程，并在另一终端验收。

```sh
uv run --python 3.12 --with mcp==1.26.0 python -m ensembl_adapter.mcp_server --transport streamable-http --port 8000
uv run --python 3.12 --with mcp==1.26.0 python acceptance_mcp.py --url http://127.0.0.1:8000/mcp
```

自定义快照位置时通过 `--db` 和 `--genome` 指定。服务检查 SQLite 声明的发布版、组装和染色体，HTTP 监听本机回环地址。平台注册时先调用 `capabilities`，再用 `resolve → neighbors → get_records` 做 YWHAB 联调。结果里的序列和坐标映射分别生成事实，坐标映射标注为计算得到。现有源快照清单只记了清单更新时间，缺少文件的实际抓取时刻，因此事实的 `retrieved_at` 返回 null 并注明原因。

[验收脚本](acceptance_mcp.py)包含首批可复跑用例，覆盖标准 ID、转录本选择、区间映射、来源、无匹配、超范围、截断和固定种子游走。[验收 trace](runs/acceptance_mcp.json)记录每一步工具输入和返回值。当前实现使用草案工具签名，事实输出限定蛋白序列和位置映射。正式工具表、完整 Fact/DatasetRecord schema、独立留出黄金集与实际 M4 注册待对齐。
