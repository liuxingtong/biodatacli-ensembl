# SciDataGym 进展 — 2026-11-10

**一句话**：连上了服务器、把 SciDataGym-Demo 的 M1 接口接到我们的 Ensembl 图上（23 条桥接检查全过），并把可视化控制台补到 6 个页签——现在"数据够不够格、接口接不接得上、答案对不对"都能在一个页面里看到。

---

## 一、服务器接入（卡了两天的那个）

之前一直失败，今天定位到根因：**`jump.pjlab.org.cn` 不是 Linux 登录节点，是 JumpServer 堡垒机**。纯公钥永远过不去，必须密码登录后选资产。

```
本机
 └─① 43.142.81.116:41237   Go 网关，只放行 tunnel，公钥 ~/.ssh/id_ed25519_pjlab 已授权
      └─ tunnel --target target-a525cc35-…
           └─② jump.pjlab.org.cn   JumpServer 堡垒机，账号 gaozhangyang + 密码
                └─③ 10.140.84.67   host-10-140-84-67, Ubuntu 24.04
                     /mnt/data3 · /mnt/data6 · /mnt/data1/gaozhangyang
```

做了三件事：

1. **封装执行**：`~/bin/pjrun '命令'`，一次 10–13 秒，自动过密码、过菜单、剥 ANSI。
2. **只读取回源码**：服务器上 `tar czf` + `base64`（必须用换行包裹，`-w0` 会卡死），sha256 双向核对一致，落到 `scidatagym-demo/`（65 KB）。**服务器上的代码一个字节没动**。
3. **写进 skill**：`biodatacli` 里记了链路、密码、公钥、`pjrun` 用法、传文件方法，下次直接可用。

拉回来的东西：`src/scidatagym/{m1,m2,m3}`、`tests`、`examples/cases.json`、`reports/`（C1–C3 真实输出）、`代码阅读指南.md`。

---

## 二、M1 桥接：能不能接上 → 能，23/23

`ensembl_proto/bridge/scidatagym_m1.py`，**两端代码零改动**，只做形状翻译：

| 他们的 M1 | 我们的实现 | 状态 |
|---|---|---|
| `describe()` | `api.capabilities()` | 已映射（补了 M3 依赖的 `views`/`catalog`/`manifest`） |
| `resolve(ref, constraints)` | `api.resolve()` | 已映射，`ok/not_found` → `resolved/ambiguous/not_found` |
| `recall(query, view, limit)` | `api.recall()` | 部分：只有 `identifier`/`field`，`semantic` 显式拒绝 |
| `neighbors(…, cursor)` | `api.neighbors()` | 已映射，cursor 切片在桥接层 |
| `get_records(…, region{start,end})` | `api.get_records()` | 已映射，`region` ↔ `protein_range` |
| `raw_record()` | — | **不提供**，返回 `not_supported` |
| `walk(…, random_seed)` | `api.walk(rng_seed)` | 已映射 |

`bridge/smoke.py` 的 23 条检查**照 `engine.py` 的四个真实调用点写**（119/128/195/315），包括负路径：不支持的关系查询前抛错、未知节点返回 `node_not_found`、未知视图拒绝而不是返回空列表。

**为什么是桥接不是直接改成他们的**：他们接口没有地方放 `null_status` 四态、`facts` 逐条溯源、`truncated+total`；直接统一会把"未报道"和"没查到"压成同一个值，溯源信息直接丢失。长期应该收敛到 `m1.v2` 超集（方法名用他们的，字段保留我们的），然后删掉 bridge——**这是下一步的第一优先**。

**接不上的不是接口，是语义**：他们 catalog 是 `Protein/GO/Keyword`，我们是 `Gene/Transcript/Exon/Protein`；C1–C3 的题干本身是 UniProt 题，换库要另写题干。

---

## 三、可视化控制台：6 个页签

`python3 viz/server.py --port 8765` → http://127.0.0.1:8765

| 页签 | 回答什么问题 |
|---|---|
| 图谱探索 | 一条结论是怎么被几条边和它自己的 source 撑起来的；蛋白节点可一键跳到序列页 |
| 序列与坐标 | 全长蛋白高亮目标区间 + 逐残基密码子/基因组坐标；**按基因比全部异构体**；**快速选择下拉框**（3 组 61 项） |
| 跨库对齐 | 一个 accession 对应多个对象时，歧义当场暴露，附证据文件/行号/sha256 |
| M1 桥接 | 接口对照 + 现场调用 + 23 条检查结果 |
| 数据库体检 | 合格判定、快照身份、来源与取回时间、覆盖率、翻译一致性、24 项对账、47 项契约分组、已知边界 |
| 评测面板 | 准/诚/寸 32 条、分类得分、逐条结果 |

今天新补的三处：

- **默认加载**：打开页面不再是空画布（之前没查询就是一片空，这是 bug）。
- **异构体视图**：YWHAB 的 40 个转录本里 23 个覆盖 230–246，canonical 标出，长度 246 vs 287 —— 同一个"230–246"在不同异构体上不是同一段，只看一个蛋白会把这个歧义藏起来。
- **下拉框**：`YWHAB 的 40 个异构体` / `有实验结构的蛋白`（按 PDB 数排序）/ `基因的 canonical 蛋白`，共 61 项，全部来自真实索引。

---

## 四、验证总账（今天结束状态）

```
47/47  契约检查          contract/conformance.py
23/23  M3 桥接检查        bridge/smoke.py
32/32  准  32/32 诚  13/13 寸   bench/run.py
24/24  与 REST 快照对账   verify.py
21/21  单元测试
```

一条 `./run_all.sh` 全跑通。

---

## 五、还没做的（按优先级）

1. **`m1.v2` 契约收敛**：把两边的字段合并成一份，双方实现，删掉 bridge。
2. **golden set 独立标注**：现在 32 条是照同一份快照写的，首轮全绿属预期；"不存在/陷阱"两类需要人工精标 + hold-out。
3. **用 Ensembl 题跑一遍他们的 M3**：机制已通，缺题干（例如"YWHAB 哪个异构体覆盖 230–246 且有实验结构"）。
4. **跨库反查**：跨库对齐页目前只能外部 → Ensembl，缺 Ensembl → 外部。
5. **M2/M3 的 capabilities 分发**：他们代码里 M3 直接读 `catalog`，越界请求靠 `unsupported_operations`；这块要在 v2 契约里定死。

---

## 六、截图

| 文件 | 内容 |
|---|---|
| `viz/screenshots/explore.png` | 图谱探索（基因邻域） |
| `viz/screenshots/explore-protein.png` | 蛋白节点 + "在序列页打开" |
| `viz/screenshots/sequence.png` | 序列 + 下拉框 + 异构体表 |
| `viz/screenshots/crosswalk.png` | 跨库对齐与歧义 |
| `viz/screenshots/bridge.png` | M1 桥接 23/23 |
| `viz/screenshots/health.png` | 数据库体检 |
| `viz/screenshots/bench.png` | 评测面板 |
