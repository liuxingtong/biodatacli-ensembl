# 接入说明（给平台/主库负责人）

本仓库实现一个数据库 = 一个 adapter = 一个独立 MCP 服务。平台只通过 8 个工具调用，
不碰库内业务逻辑。

## 一、启动

```bash
pip install -r requirements.txt          # 契约层与测试只需要 jsonschema
pip install -r requirements-mcp.txt      # MCP 服务；注意 mcp 需 <2，v2 改了 FastMCP 的导入路径

# 数据：下载 + 建图 + 全量检查（约 1 小时，主要耗在下载）
python3 fetch_data.py            # 8 个冻结输入，约 48 MB
python3 make_enrichment.py
python3 -m ensembl_adapter.build --gff3 data/raw/Homo_sapiens.GRCh38.116.chromosome.20.gff3.gz \
  --genome data/raw/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz \
  --pep data/raw/Homo_sapiens.GRCh38.pep.all.fa.gz \
  --biomart data/raw/biomart --enrichment data/raw/xref_enrichment.json \
  --capabilities capabilities.json \
  --db data/index/ensembl_chr20_r116.sqlite --stats data/index/build_stats.json
python3 crosswalk/build.py
./run_all.sh                     # 契约 70 + bench 32 + 桥接 23 + 单测 26

# 服务
python3 -m ensembl_adapter.mcp_server --db data/index/ensembl_chr20_r116.sqlite
```

## 二、工具与返回

必选四个：`describe / resolve / search / get_facts`；
可选四个：`search_datasets / get_dataset_metadata / get_access_plan / expand`。

| 工具 | 入参 | 出参要点 |
|---|---|---|
| describe | `{}` | adapter_id、entities、scales、keyspaces、evidence_types、search_modes、license |
| resolve | `{term, context?, species?}` | `[{canonical_id, namespace, score, aliases[]}]`，附 `ambiguous` |
| search | `{query, filters?, k}` | `[{entity, attribute, value, evidence_type, source, score}]` |
| get_facts | `{entity_ids[], fields?, scale?}` | `{facts[], coverage{requested, entities_with_facts, missing}}` |
| search_datasets | `{query, filters?, k}` | 数据集目录 |
| get_dataset_metadata | `{dataset_id}` | 完整 DatasetRecord |
| get_access_plan | `{dataset_id, target}` | `{fetch_url_or_cmd, bytes, format, license}` |
| expand | `{entity_ids[], key, max_depth?, k}` | `[{from, to, relation, evidence, score}]` |

统一约定：

- 每条 fact 自带 `source{db,id,version,retrieved_at,retrieved_at_status,field,locator}`；
- 查不到时返回 `null_status`，取值 `not_found / not_reported / na / access_failed`，不编造；
- 版本：`contract_version = m1.1`，索引里带 `build_id` 与 `built_at`；
- `ENSEMBL_TRACE=<file>` 时每个调用落一行 JSONL（工具、参数、状态、耗时、build_id）。

## 三、能力边界（先说清楚，避免误用）

| 项 | 现状 |
|---|---|
| 物种/范围 | 人类 GRCh38；chr20 切片（默认）或全基因组（`full/`） |
| 证据类型 | 只产出 `computational`、`inferred`；无实验/文献级断言 |
| 检索模式 | `identifier`、`text`；**没有语义（LSA/向量）检索** |
| 标准键 | ENSG、UniProtKB、PDB、GO、Reactome、KEGG、EC；CL、UBERON 未覆盖 |
| GO/Reactome/KEGG/EC | 在 crosswalk 里，不是图边；图边只有 UniProt/PDB/HGNC/Entrez/MIM/CCDS |
| `raw_record` | 不提供（切片无“一个实体一条原始记录”） |
| 变异/调控/比较基因组 | 不提供 |
| 全量库 | 3.26 GB，2,622,629 节点；crosswalk 首次加载约 35 s |

## 四、验收方式

```bash
python3 contract/conformance.py     # 70 项：schema、null 四态、溯源、截断、幂等
python3 bench/run.py                # 32 题，准/诚/寸
python3 bridge/smoke.py             # 23 项，按 SciDataGym-Demo 的 M3 调用方式
python3 -m unittest discover -s tests   # 26 项
python3 acceptance_mcp.py --db <构建好的索引>   # 走真实 MCP 会话，校验 8 个工具
# 已在本机验证：PASS: 8 platform tools, 8 calls（python 3.12 + mcp 1.x）
```

bench 的 32 题全部由冻结数据推导（每题记录 `source_of_truth`），定位是回归与口径固定，
不作为对外准确率背书。

## 五、已知待办

1. `search` 依赖实体召回后再展开 facts，没有独立的事实级索引；大结果集需要分页。
2. 全量库的 crosswalk 加载 35 s、占内存约 1 GB，首次查询明显慢。
3. 没有 CI；检查靠 `./run_all.sh` 本地跑。
