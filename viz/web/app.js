/* SciDataGym · Ensembl console — vanilla JS, no build step. */

const $ = (sel) => document.querySelector(sel);
const esc = (value) => String(value ?? "").replace(/[&<>"']/g, (c) => (
  { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

const TYPE_COLOR = {
  Gene: "#4cc2ff", Transcript: "#7ef0c0", Protein: "#c9a3ff",
  Exon: "#ffd166", StructureEntry: "#ff9f6b", CodingSequence: "#9ad0ff",
  Assembly: "#8ea0c8", Region: "#8ea0c8",
};

async function api(path, params = {}) {
  const query = new URLSearchParams(
    Object.entries(params).filter(([, v]) => v !== undefined && v !== null && v !== ""));
  const response = await fetch(`/api/${path}?${query}`);
  const payload = await response.json();
  if (!response.ok) throw new Error(payload.error || response.statusText);
  return payload;
}

/* ------------------------------------------------------------------ */
/* boot                                                                */
/* ------------------------------------------------------------------ */
let caps = null;

async function boot() {
  try {
    caps = await api("capabilities");
    const snap = caps.snapshot || {};
    const regions = String((caps.scope || {}).region || "").split(",").filter(Boolean);
    const scopeLabel = regions.length > 1
      ? `· Ensembl ${regions.length} 条染色体`
      : `· Ensembl chr${regions[0] || "?"}`;
    const title = $("#scope-title");
    if (title) title.textContent = scopeLabel;
    $("#meta").textContent =
      `release ${caps.release} · ${caps.assembly} · chr${(caps.scope || {}).region} · ` +
      `build ${String(snap.build_id || "").slice(0, 8)} · contract ${caps.contract_version}`;
    $("#legend").innerHTML = Object.entries(TYPE_COLOR)
      .filter(([type]) => (caps.node_types || []).includes(type))
      .map(([type, color]) => `<span><i style="background:${color}"></i>${type}</span>`).join("")
      + `<span><i style="background:#8ea0c8"></i>共 ${(caps.crosswalk || {}).edges || 0} 条跨库边</span>`;
  } catch (error) {
    $("#meta").textContent = `服务不可用：${error.message}`;
  }
}

/* ------------------------------------------------------------------ */
/* tabs + forms                                                        */
/* ------------------------------------------------------------------ */
const loaded = { eval: false, seq: false, cross: false, health: false, bridge: false };

document.querySelectorAll("#tabs button").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll("#tabs button").forEach((item) => item.classList.remove("active"));
    document.querySelectorAll(".view").forEach((view) => view.classList.remove("active"));
    button.classList.add("active");
    const name = button.dataset.view;
    $(`#view-${name}`).classList.add("active");
    history.replaceState(null, "", `#${name}`);
    if (name === "eval" && !loaded.eval) { loaded.eval = true; renderEval(); }
    if (name === "health" && !loaded.health) { loaded.health = true; renderHealth(); }
    if (name === "bridge" && !loaded.bridge) { loaded.bridge = true; renderBridge(); }
    if (name === "seq" && !loaded.seq) {
      loaded.seq = true;
      loadSequence("Ensembl:ENSP00000300161", 230, 246);
      loadIsoforms("YWHAB", 230, 246);
    }
    if (name === "cross" && !loaded.cross) { loaded.cross = true; loadCrosswalk("P31946"); }
  });
});

$("#search").addEventListener("submit", (event) => {
  event.preventDefault();
  doSearch($("#q").value.trim(), $("#type").value);
});

function bindSeeds() {
  document.querySelectorAll("[data-seed]").forEach((button) => {
    button.addEventListener("click", () => {
      $("#q").value = button.dataset.seed;
      doSearch(button.dataset.seed, "");
    });
  });
}
bindSeeds();

$("#seq-pick").addEventListener("change", (event) => {
  const value = event.target.value;
  if (!value) return;
  $("#seq-id").value = value;
  loadSequence(value, Number($("#seq-start").value), Number($("#seq-end").value));
});

$("#seq-form").addEventListener("submit", (event) => {
  event.preventDefault();
  loadSequence($("#seq-id").value.trim(), Number($("#seq-start").value), Number($("#seq-end").value));
});
$("#seq-demo").addEventListener("click", () => {
  $("#seq-start").value = 1; $("#seq-end").value = 40;
  loadSequence($("#seq-id").value.trim(), 1, 40);
});
$("#iso-form").addEventListener("submit", (event) => {
  event.preventDefault();
  loadIsoforms($("#iso-gene").value.trim(), Number($("#iso-start").value), Number($("#iso-end").value));
});
$("#cross-form").addEventListener("submit", (event) => {
  event.preventDefault();
  loadCrosswalk($("#cross-q").value.trim());
});

/* ------------------------------------------------------------------ */
/* explore                                                             */
/* ------------------------------------------------------------------ */
async function doSearch(query, type) {
  if (!query) return;
  const url = new URL(location.href);
  url.searchParams.set("q", query);
  if (type) url.searchParams.set("type", type); else url.searchParams.delete("type");
  history.replaceState(null, "", url);
  showLoading("检索中…");
  try {
    const result = await api("search", { q: query, type: type || "" });
    if (!result.candidates || !result.candidates.length) {
      renderMiss(query, result);
      return;
    }
    const preferred = result.candidates.find((node) => node.namespace !== "UniProtKB") || result.candidates[0];
    await loadNode(preferred.id);
  } catch (error) {
    renderError(error);
  }
}

function renderMiss(query, result) {
  clearGraph();
  const status = result.null_status || result.status;
  const explain = {
    not_found: "检索执行了，但这个切片里没有匹配的实体。",
    na: "该查询超出本适配器声明的范围（例如非人类物种、未实现的检索视图）。",
    access_failed: "底层读取失败，因此不作任何断言。",
    not_reported: "实体存在，但来源没有给出这个值。",
  }[status] || "未返回结果。";
  $("#detail").innerHTML = `
    <div class="empty">
      <h3>没有答案</h3>
      <p><code>${esc(query)}</code> → <span class="badge bad">${esc(status)}</span></p>
      <p>${esc(explain)}</p>
      <p class="hint">ok = ${esc(result.ok)} · status = ${esc(result.status)}</p>
    </div>`;
}

async function loadNode(id) {
  showLoading(`展开 ${id} …`);
  try {
    const graph = await api("graph", { id, limit: 40 });
    if (!graph.node) { renderMiss(id, graph); return; }
    $("#q").value = id;
    drawGraph(graph);
    renderDetail(graph);
  } catch (error) {
    renderError(error);
  }
}

function clearGraph() {
  $("#graph").innerHTML = `<text x="450" y="310" text-anchor="middle" fill="#5a6b91"
    font-size="13">在右上角搜索一个实体，或点下面的示例</text>`;
}
function showLoading(text) { $("#detail").innerHTML = `<div class="empty"><h3>${esc(text)}</h3></div>`; }
function renderError(error) {
  $("#detail").innerHTML = `<div class="empty"><h3 class="bad-txt">请求失败</h3><p>${esc(error.message)}</p></div>`;
}

/* Force layout: small enough to be readable, deterministic enough to trust. */
function layout(nodes, links, width, height) {
  const margin = 70;
  const cx = width / 2, cy = height / 2;
  nodes.forEach((node, index) => {
    if (node.root) { node.x = cx; node.y = cy; return; }
    const angle = (index / Math.max(1, nodes.length - 1)) * Math.PI * 2;
    const radius = 170 + (index % 3) * 40;
    node.x = cx + Math.cos(angle) * radius;
    node.y = cy + Math.sin(angle) * radius * 0.78;
  });

  const byId = new Map(nodes.map((node) => [node.id, node]));
  for (let step = 0; step < 520; step += 1) {
    const alpha = 1 - step / 520;
    for (let i = 0; i < nodes.length; i += 1) {
      for (let j = i + 1; j < nodes.length; j += 1) {
        const a = nodes[i], b = nodes[j];
        let dx = b.x - a.x, dy = b.y - a.y;
        let dist = Math.hypot(dx, dy) || 0.01;
        const force = (5400 * alpha) / (dist * dist);
        dx = (dx / dist) * force; dy = (dy / dist) * force;
        a.x -= dx; a.y -= dy; b.x += dx; b.y += dy;
      }
    }
    links.forEach((link) => {
      const a = byId.get(link.source), b = byId.get(link.target);
      if (!a || !b) return;
      const dx = b.x - a.x, dy = b.y - a.y;
      const dist = Math.hypot(dx, dy) || 0.01;
      const target = 165;
      const force = ((dist - target) / dist) * 0.045 * alpha;
      a.x += dx * force; a.y += dy * force;
      b.x -= dx * force; b.y -= dy * force;
    });
    nodes.forEach((node) => {
      if (node.root) { node.x = cx; node.y = cy; return; }
      node.x += (cx - node.x) * 0.006 * alpha;
      node.y += (cy - node.y) * 0.006 * alpha;
      node.x = Math.max(margin, Math.min(width - margin, node.x));
      node.y = Math.max(margin, Math.min(height - margin, node.y));
    });
  }
  return { byId };
}

function drawGraph(graph) {
  const svg = $("#graph");
  const width = 900, height = 620;
  const root = { id: graph.node.id, label: shortLabel(graph.node), type: graph.node.type, root: true };
  const nodes = [root];
  const links = [];
  const seen = new Set([root.id]);
  graph.items.forEach((item) => {
    const node = item.node;
    if (!seen.has(node.id)) {
      seen.add(node.id);
      nodes.push({ id: node.id, label: shortLabel(node), type: node.type });
    }
    links.push({
      source: item.edge.from, target: item.edge.to,
      relation: item.edge.relation,
      direction: item.direction,
    });
  });

  const { byId } = layout(nodes, links, width, height);
  const edgeMarkup = links.map((link) => {
    const a = byId.get(link.source), b = byId.get(link.target);
    if (!a || !b) return "";
    return `<line class="edge" x1="${a.x.toFixed(1)}" y1="${a.y.toFixed(1)}"
              x2="${b.x.toFixed(1)}" y2="${b.y.toFixed(1)}"></line>`;
  }).join("");

  const nodeMarkup = nodes.map((node) => {
    const color = TYPE_COLOR[node.type] || "#8ea0c8";
    const r = node.root ? 11 : 7;
    const anchor = node.x > width * 0.62 ? "end" : "start";
    const dx = node.x > width * 0.62 ? -14 : 14;
    return `<g class="node ${node.root ? "root" : ""}" data-id="${esc(node.id)}"
               transform="translate(${node.x.toFixed(1)},${node.y.toFixed(1)})">
              <circle r="${r}" fill="${color}" fill-opacity="${node.root ? 0.95 : 0.75}"
                      stroke="${color}" stroke-opacity=".9"></circle>
              <text x="${dx}" y="4" text-anchor="${anchor}">${esc(node.label)}</text>
            </g>`;
  }).join("");

  svg.innerHTML = `<g class="edges">${edgeMarkup}</g><g class="nodes">${nodeMarkup}</g>`;
  svg.querySelectorAll(".node").forEach((group) => {
    group.addEventListener("click", () => loadNode(group.dataset.id));
  });
  $("#graph-hint").textContent =
    `${nodes.length - 1} 个邻居 · ${links.length} 条边 · 点击节点继续展开`;
}

function shortLabel(node) {
  const id = String(node.id || "");
  const stable = node.stable_id || id.split(":").pop();
  if (node.type === "StructureEntry") return stable;
  if (node.type === "Protein") return `ENSP…${stable.slice(-6)}`;
  if (node.type === "Transcript") return `ENST…${stable.slice(-6)}`;
  if (node.type === "Gene" && node.label) return node.label;
  return stable;
}

async function renderDetail(graph) {
  const node = graph.node;
  const source = node.source || {};
  const items = graph.items || [];

  const location = node.location
    ? `chr${node.location.seq_region}:${node.location.start.toLocaleString()}–${node.location.end.toLocaleString()} (${node.location.strand > 0 ? "+" : "−"})`
    : "—";

  const rows = [
    ["stable ID", node.stable_id ? `${node.stable_id}.${node.version ?? ""}` : "—"],
    ["biotype", node.biotype || "—"],
    ["location", location],
    ["length", node.length ?? "—"],
  ].map(([k, v]) => `<tr><td>${k}</td><td>${esc(v)}</td></tr>`).join("");

  const badgeClass = source.retrieved_at_status === "download" ? "ok"
    : source.retrieved_at_status === "file_mtime" ? "warn" : "neutral";

  const ambiguous = graph.crosswalk_ambiguous
    ? `<div class="note" style="margin:12px 0 0">该外部标识符对应多个 Ensembl 对象，已全部保留：
        ${(graph.crosswalk || []).filter((hit) => hit.alternatives > 0)
          .slice(0, 4).map((hit) => `<code>${esc(hit.other.namespace)}:${esc(hit.other.id)}</code>`).join("、")}
        （见「跨库对齐」页）</div>`
    : "";

  const neighbors = items.length ? items.map((item) => `
    <li data-id="${esc(item.node.id)}">
      <span class="id">${esc(item.node.id)}</span>
      <span class="rel">${esc(item.edge.relation)} · ${esc(item.direction)}</span>
    </li>`).join("") : `<li class="hint">该节点没有可展开的邻居</li>`;

  const seqAction = node.type === "Protein" && String(node.id).startsWith("Ensembl:")
    ? `<button class="link" id="open-seq" data-id="${esc(node.id)}">在序列页打开 →</button>`
    : "";

  $("#detail").innerHTML = `
    <h3>${esc(node.id)}<span class="type-chip">${esc(node.type)}</span></h3>
    ${seqAction}
    ${ambiguous}
    <table class="kv">${rows}</table>
    <div class="src">
      <div class="row"><span class="k">source</span>
        <span>${esc(source.db || source.database)}:${esc(source.id || source.record_id)}@${esc(source.version || source.release)}</span></div>
      <div class="row"><span class="k">field</span><span>${esc(source.field || "—")}</span></div>
      <div class="row"><span class="k">locator</span><span>${esc(source.locator || "—")}</span></div>
      <div class="row"><span class="k">retrieved</span>
        <span>${esc(source.retrieved_at || "—")}
          <span class="badge ${badgeClass}">${esc(source.retrieved_at_status || "unknown")}</span></span></div>
    </div>
    <div id="facts-block"></div>
    <h4 class="hint" style="margin:16px 0 4px">邻居（${items.length}${graph.truncated ? "，已截断" : ""}）</h4>
    <ul class="nbr">${neighbors}</ul>
    <div class="chips">${sideNote(graph)}</div>`;

  document.querySelectorAll(".nbr li[data-id]").forEach((row) => {
    row.addEventListener("click", () => loadNode(row.dataset.id));
  });

  const opener = $("#open-seq");
  if (opener) {
    opener.addEventListener("click", () => {
      document.querySelector('#tabs button[data-view="seq"]').click();
      $("#seq-id").value = opener.dataset.id;
      loadSequence(opener.dataset.id, Number($("#seq-start").value), Number($("#seq-end").value));
    });
  }

  renderFacts(graph.facts || []);
}

function sideNote(graph) {
  if (graph.truncated) return `总邻居 ${graph.total}，已按 limit 截断`;
  if (graph.crosswalk_ambiguous) return "跨库命中存在歧义，见「跨库对齐」页";
  return "";
}

function renderFacts(all) {
  const box = $("#facts-block");
  if (!box) return;
  const facts = (all || []).slice(0, 6);
  if (!facts.length) {
    box.innerHTML = `<div class="chips" style="margin-top:12px">该节点没有可展示的断言</div>`;
    return;
  }
  {
    const total = (all || []).length;
    box.innerHTML = `
      <h4 class="hint" style="margin:16px 0 4px">断言（${total} 条，显示前 ${facts.length}）</h4>
      <table class="data">
        <thead><tr><th>attribute</th><th>evidence</th><th>source</th><th>retrieved</th></tr></thead>
        <tbody>${facts.map((fact) => `
          <tr>
            <td>${esc(fact.attribute)}</td>
            <td><span class="badge ${fact.evidence_type === "inferred" ? "warn" : "ok"}">${esc(fact.evidence_type)}</span></td>
            <td>${esc(fact.source.id)}@${esc(fact.source.version)}</td>
            <td>${esc((fact.source.retrieved_at || "").slice(0, 19))}</td>
          </tr>`).join("")}</tbody>
      </table>`;
  }
}

/* ------------------------------------------------------------------ */
/* sequence                                                            */
/* ------------------------------------------------------------------ */
let pickerStatic = null;   // {structures, genes}
let pickerIso = null;      // {label, items}

function optgroup(label, items) {
  if (!items || !items.length) return "";
  return `<optgroup label="${esc(label)}">` + items.map((item) =>
    `<option value="${esc(item.value)}">${esc(item.label)}</option>`).join("") + `</optgroup>`;
}

function renderPicker() {
  const select = $("#seq-pick");
  if (!select) return;
  select.innerHTML = `<option value="">— 从列表选一个 —</option>`
    + (pickerIso ? optgroup(pickerIso.label, pickerIso.items) : "")
    + (pickerStatic ? optgroup("有实验结构的蛋白", pickerStatic.structures)
        + optgroup("基因的 canonical 蛋白", pickerStatic.genes) : "");
}

async function fillTargets() {
  try {
    pickerStatic = await api("targets");
  } catch (error) {
    pickerStatic = null;
  }
  renderPicker();
  renderSeeds();
}

function renderSeeds() {
  const boxes = [$("#seeds"), $("#quick-seeds")].filter(Boolean);
  if (!boxes.length) return;
  const picks = [
    { value: "YWHAB", label: "YWHAB" },
    { value: "P31946", label: "P31946" },
    { value: "4DNK", label: "4DNK" },
    ...((pickerStatic && pickerStatic.genes) || []).slice(3, 7).map((item) => ({
      value: item.value, label: item.label.split(" · ")[0],
    })),
  ];
  const markup = picks.map((pick, index) =>
    `<button class="link" data-seed="${esc(pick.value)}">${esc(pick.label)}</button>`
    + (index < picks.length - 1 ? "、" : "")).join("");
  boxes.forEach((box) => { box.innerHTML = markup; });
  bindSeeds();
}

async function loadSequence(id, start, end) {
  const box = $("#seq-out");
  box.innerHTML = `<div class="empty"><h3>取数中…</h3></div>`;
  try {
    const region = start && end ? `protein_range:${start}-${end}` : "";
    const payload = await api("records", {
      ids: id, fields: "sequence,protein_mapping", region,
    });
    const record = (payload.records || [])[0];
    if (!record || record.status !== "ok") {
      box.innerHTML = `<div class="empty"><h3>没有取到记录</h3>
        <p><code>${esc(id)}</code> → <span class="badge bad">${esc(payload.null_status || payload.status)}</span></p></div>`;
      return;
    }
    const sequence = (record.fields || {}).sequence || {};
    const facts = payload.facts || [];
    const mapping = record.region || {};
    const full = sequence.sequence || "";
    const [from, to] = mapping.requested || [start, end];
    const before = full.slice(0, Math.max(0, from - 1));
    const middle = full.slice(Math.max(0, from - 1), to);
    const after = full.slice(to);
    const positions = mapping.positions || [];
    const missing = mapping.missing || [];

    box.innerHTML = `
      <div class="seq-meta">
        <span>蛋白 <b>${esc(record.id)}</b></span>
        <span>全长 <b>${full.length}</b> aa</span>
        <span>取回 <b>${positions.filter((p) => p.status === "ok").length}</b> / ${positions.length}</span>
        <span>缺失 <b class="${missing.length ? "bad-txt" : "good-txt"}">${missing.length}</b></span>
        <span>证据 <b>${facts.length}</b> 条</span>
      </div>
      <div class="seq-wrap">${esc(before)}<span class="hl">${esc(middle)}</span>${esc(after)}</div>
      <table class="data">
        <thead><tr><th>protein pos</th><th>codon</th><th>genomic interval</th><th>transcript</th></tr></thead>
        <tbody>${positions.map((pos) => `
          <tr>
            <td>${esc(pos.protein_pos)}</td>
            <td>${esc(pos.codon || "—")}</td>
            <td>${(pos.genomic_interval || []).join("–") || "—"}</td>
            <td>${esc(pos.transcript_id || mapping.transcript_id || "")}</td>
          </tr>`).join("")}</tbody>
      </table>`;
  } catch (error) {
    box.innerHTML = `<div class="empty"><h3 class="bad-txt">请求失败</h3><p>${esc(error.message)}</p></div>`;
  }
}


/* ------------------------------------------------------------------ */
/* isoforms: every transcript of a gene against one window            */
/* ------------------------------------------------------------------ */
async function loadIsoforms(gene, start, end) {
  const box = $("#iso-out");
  if (!box) return;
  box.innerHTML = `<div class="empty"><h3>比对中…</h3></div>`;
  try {
    const d = await api("isoforms", { gene, start, end });
    if (!d.rows || !d.rows.length) {
      box.innerHTML = `<div class="empty"><h3>没有找到 ${esc(gene)} 的转录本</h3>
        <p><code>${esc(d.status || "not_found")}</code></p></div>`;
      return;
    }
    const badge = (status) => {
      const cls = status === "covered" ? "ok" : status === "partial_model" ? "warn" : "bad";
      return `<span class="badge ${cls}">${esc(status)}</span>`;
    };
    const rows = d.rows.map((row) => `
      <tr>
        <td>${esc(row.transcript)}${row.canonical ? ' <span class="badge ok">canonical</span>' : ""}</td>
        <td>${esc(row.protein || "—")}</td>
        <td>${row.protein_length ?? "—"}</td>
        <td>${row.retrieved}/${row.expected}</td>
        <td>${badge(row.status)}</td>
      </tr>`).join("");
    pickerIso = {
      label: `${d.gene} 的 ${d.rows.length} 个异构体`,
      items: d.rows.filter((row) => row.protein).map((row) => ({
        value: `Ensembl:${row.protein}`,
        label: `${row.transcript}${row.canonical ? " (canonical)" : ""} · ${row.protein} · ${row.protein_length} aa`,
      })),
    };
    renderPicker();

    box.innerHTML = `
      <div class="panel">
        <div class="panel-head">
          <h2>${esc(d.gene_id)} · ${d.rows.length} 个转录本 · ${d.covered} 个覆盖 ${d.region[0]}–${d.region[1]}</h2>
          <span class="chips">蛋白长度不同 → 同一个位点在不同异构体上含义不同</span>
        </div>
        <table class="data">
          <thead><tr><th>转录本</th><th>蛋白</th><th>长度</th><th>取回/应得</th><th>状态</th></tr></thead>
          <tbody>${rows}</tbody>
        </table>
      </div>`;
  } catch (error) {
    box.innerHTML = `<div class="empty"><h3 class="bad-txt">比对失败</h3><p>${esc(error.message)}</p></div>`;
  }
}

/* ------------------------------------------------------------------ */
/* crosswalk                                                           */
/* ------------------------------------------------------------------ */
async function loadCrosswalk(query) {
  const box = $("#cross-out");
  box.innerHTML = `<div class="empty"><h3>查询中…</h3></div>`;
  try {
    const payload = await api("crosswalk", { q: query });
    const hits = payload.hits || [];
    if (!hits.length) {
      box.innerHTML = `<div class="empty"><h3>没有跨库命中</h3>
        <p><code>${esc(query)}</code> 不在对照表中（或不是可识别的标识符）。</p>
        <p class="hint">对照表覆盖 Ensembl / UniProtKB / PDB / HGNC / EntrezGene / MIM / CCDS。</p></div>`;
      return;
    }
    const ambiguous = payload.crosswalk_ambiguous;
    const banner = ambiguous
      ? `<div class="note">该标识符对应多个对象（fan-out > 1）。这里全部列出，不替你做选择——这正是「歧义显式化」。</div>`
      : `<div class="chips">一对一映射，无歧义。</div>`;
    box.innerHTML = `
      <h3 style="font-family:var(--mono)">${esc(query)}</h3>
      ${banner}
      <table class="data">
        <thead><tr><th>边</th><th>关系</th><th>基数</th><th>备选</th><th>证据（文件 · 行 · 取回）</th></tr></thead>
        <tbody>${hits.map((hit) => {
          const from = `${hit.anchor.namespace}:${hit.anchor.id}`;
          const to = `${hit.other.namespace}:${hit.other.id}`;
          const ev = hit.evidence || {};
          const short = (ev.sha256 || "").slice(0, 8);
          return `<tr>
            <td>${esc(from)} → ${esc(to)}</td>
            <td>${esc(hit.relation)}</td>
            <td>${esc(hit.cardinality)}</td>
            <td>${hit.alternatives > 0 ? `<span class="badge warn">+${hit.alternatives}</span>` : "—"}</td>
            <td>${esc((ev.file || "").split("/").pop())}${ev.row ? `:${ev.row}` : ""}
              · ${esc(short)} · ${esc((ev.retrieved_at || "").slice(0, 10))}</td>
          </tr>`;
        }).join("")}</tbody>
      </table>`;
  } catch (error) {
    box.innerHTML = `<div class="empty"><h3 class="bad-txt">请求失败</h3><p>${esc(error.message)}</p></div>`;
  }
}



/* ------------------------------------------------------------------ */
/* bridge: the SciDataGym-Demo M1 surface over this slice             */
/* ------------------------------------------------------------------ */
const BRIDGE_MAP = [
  ["describe()", "api.capabilities()", "已映射", "补了 views / catalog / manifest 三个 M3 依赖字段"],
  ["resolve(ref, constraints)", "api.resolve()", "已映射", "status 转成 resolved / ambiguous / not_found"],
  ["recall(query, view, limit)", "api.recall()", "部分", "视图只有 identifier / field；semantic 显式拒绝"],
  ["neighbors(id, rel, dir, limit, cursor)", "api.neighbors()", "已映射", "cursor 在桥接层切片实现"],
  ["get_records(ids, fields, region{start,end})", "api.get_records()", "已映射", "region 双向翻译成 protein_range"],
  ["raw_record(accession)", "—", "不提供", "返回 not_supported，不编造原文"],
  ["walk(seeds, rel, steps, budget, seed)", "api.walk()", "已映射", "random_seed → rng_seed"],
  ["catalog / manifest 属性", "capabilities / snapshot", "已补齐", "从 capabilities.json 翻译"],
];

function renderBridgeRow([theirs, ours, state, note]) {
  const cls = state === "已映射" || state === "已补齐" ? "ok"
    : state === "部分" ? "warn" : "bad";
  return `<tr>
    <td style="font-family:var(--mono)">${esc(theirs)}</td>
    <td style="font-family:var(--mono)">${esc(ours)}</td>
    <td><span class="badge ${cls}">${esc(state)}</span></td>
    <td style="font-family:inherit;white-space:normal">${esc(note)}</td>
  </tr>`;
}

async function renderBridge() {
  const box = $("#bridge-out");
  box.innerHTML = `<div class="empty"><h3>连接 M3 …</h3></div>`;
  let d;
  try {
    d = await api("bridge");
  } catch (error) {
    box.innerHTML = `<div class="empty"><h3 class="bad-txt">桥接层不可用</h3><p>${esc(error.message)}</p></div>`;
    return;
  }
  const smoke = d.smoke || {};
  const checks = smoke.checks || [];
  const passed = checks.filter((c) => c.pass).length;
  const demos = d.demos || {};
  const neighbors = demos.neighbors || {};

  box.innerHTML = `
    <div class="verdict">
      <div>
        <div class="kpi-num good" style="font-size:22px">M3 可以直接驱动这个 Ensembl 切片</div>
        <div class="chips" style="margin-top:4px">
          桥接层 <code>bridge/scidatagym_m1.py</code>：不改两端代码，只做形状翻译；下面 ${checks.length} 条检查按他们 M3 的调用方式逐条跑
        </div>
      </div>
      <div class="verdict-badges">
        <span class="badge ${passed === checks.length ? "ok" : "bad"}">桥接 ${passed}/${checks.length}</span>
        <span class="badge warn">视图 ${esc((d.views || []).join(" / "))}</span>
      </div>
    </div>

    <div class="panel">
      <div class="panel-head"><h2>接口对照</h2>
        <span class="chips">左：SciDataGym-Demo 的 M1 · 右：我们的实现</span></div>
      <table class="data">
        <thead><tr><th>他们的方法</th><th>我们的方法</th><th>状态</th><th>说明</th></tr></thead>
        <tbody>${BRIDGE_MAP.map(renderBridgeRow).join("")}</tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>现场跑一遍（真实调用，不是截图）</h2>
        <span class="chips">release ${esc((d.manifest || {}).release)} · build ${esc(((d.manifest || {}).build_id || "").slice(0, 8))}</span></div>
      <table class="data">
        <thead><tr><th>调用</th><th>结果</th></tr></thead>
        <tbody>
          <tr><td>resolve("P31946")</td>
              <td><span class="badge ${demos.resolve.status === "ambiguous" ? "warn" : "ok"}">${esc(demos.resolve.status)}</span>
                  ${(demos.resolve.items || []).map((i) => esc(i.id)).join(" · ")}</td></tr>
          <tr><td>neighbors(YWHAB, has_transcript, 3)</td>
              <td>${esc((neighbors.page1 || []).join(", "))} <span class="chips">next_cursor=${esc(neighbors.next_cursor)}</span></td></tr>
          <tr><td>neighbors(cursor=3)</td>
              <td>${esc((neighbors.page2 || []).join(", "))}</td></tr>
          <tr><td>get_records(ENSP00000300161, 230–246)</td>
              <td><code>${esc(demos.get_records.sequence)}</code> ·
                  region ${esc(JSON.stringify(demos.get_records.region))}</td></tr>
          <tr><td>raw_record("P31946")</td>
              <td><span class="badge bad">${esc(demos.raw_record.status)}</span> ${esc(demos.raw_record.reason || "")}</td></tr>
        </tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>逐条桥接检查（${passed}/${checks.length}）</h2>
        <span class="chips">来自 bridge/smoke.py，按 engine.py 的四个调用点组织</span></div>
      <table class="data">
        <thead><tr><th>检查</th><th>结果</th></tr></thead>
        <tbody>${checks.map((c) => `<tr>
          <td style="font-family:inherit">${esc(c.name)}</td>
          <td><span class="badge ${c.pass ? "ok" : "bad"}">${c.pass ? "PASS" : "FAIL"}</span></td>
        </tr>`).join("")}</tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>没有假装能做的部分</h2></div>
      <ul class="limit-list">
        <li><b>semantic 视图</b>：这版没有 LSA 索引，<code>describe().views</code> 只声明 <code>identifier / field</code>，M3 会降级而不是崩。</li>
        <li><b>raw_record</b>：Ensembl 切片没有“一个实体对应一条原始 DAT 记录”的概念，返回 <code>not_supported</code>。</li>
        <li><b>语义不通用</b>：他们 catalog 是 Protein/GO/Keyword，我们是 Gene/Transcript/Exon/Protein，关系词表也不同；能对接的是机制，C1–C3 的题本身换不过来。</li>
        <li><b>范围</b>：只有 human chr20，越界一律 <code>na</code>。</li>
      </ul>
    </div>`;
}

/* ------------------------------------------------------------------ */
/* health: judge the index itself, not one answer                     */
/* ------------------------------------------------------------------ */
async function renderHealth() {
  const box = $("#health-out");
  box.innerHTML = `<div class="empty"><h3>体检中…</h3></div>`;
  const d = await api("health");
  const caps = d.capabilities || {};
  const snap = caps.snapshot || {};
  const prov = caps.provenance || {};
  const conformance = d.conformance || {};
  const verify = d.verify || {};
  const counts = d.counts || {};
  const coverage = d.coverage || {};
  const validation = d.validation || {};

  const conformanceChecks = conformance.checks || [];
  const conformancePassed = conformanceChecks.filter((c) => c.pass).length;
  const verifyChecks = verify.checks || [];
  const verifyPassed = verifyChecks.filter((c) => (c.status || "").toUpperCase() === "PASS").length;
  const provEntries = Object.entries(prov);
  const provStamped = provEntries.filter(([, v]) => v.retrieved_at).length;

  const healthy = conformance.ok === true && verify.failed === 0 && caps.consistent === true;

  const pct = (part, whole) => (whole ? Math.round((part * 1000) / whole) / 10 : 0);
  const bar = (part, whole, label, value) => {
    const ratio = pct(part, whole);
    return `<tr>
      <td style="font-family:inherit">${esc(label)}</td>
      <td>${value}</td>
      <td>${ratio}%</td>
      <td style="width:240px"><div class="bar"><i style="width:${ratio}%"></i></div></td>
    </tr>`;
  };

  const provRows = provEntries.map(([key, entry]) => {
    const status = entry.retrieved_at_status || "not_reported";
    const cls = status === "download" ? "ok" : status === "file_mtime" ? "warn" : "bad";
    return `<tr>
      <td style="font-family:var(--mono)">${esc((entry.path || key).split("/").pop())}</td>
      <td>${esc((entry.sha256 || "").slice(0, 12) || "—")}</td>
      <td>${esc((entry.retrieved_at || "—").slice(0, 19))}</td>
      <td><span class="badge ${cls}">${esc(status)}</span></td>
    </tr>`;
  }).join("");

  const verifyRows = verifyChecks.map((check) => {
    const ok = (check.status || "").toUpperCase() === "PASS";
    return `<tr>
      <td style="font-family:inherit">${esc(check.name)}</td>
      <td><span class="badge ${ok ? "ok" : "bad"}">${ok ? "PASS" : "FAIL"}</span></td>
      <td style="font-family:var(--mono);color:var(--muted)">${esc(String(check.detail || "").slice(0, 70))}</td>
    </tr>`;
  }).join("");

  const groups = {};
  conformanceChecks.forEach((check) => {
    const name = String(check.name || "");
    const group = name.includes(":") ? name.split(":")[0] : "其他";
    groups[group] = groups[group] || { passed: 0, total: 0 };
    groups[group].total += 1;
    if (check.pass) groups[group].passed += 1;
  });
  const groupRows = Object.entries(groups).map(([group, bucket]) => {
    const ok = bucket.passed === bucket.total;
    return `<tr>
      <td style="font-family:inherit">${esc(group)}</td>
      <td>${bucket.passed}/${bucket.total}</td>
      <td><span class="badge ${ok ? "ok" : "bad"}">${ok ? "全通过" : "有失败"}</span></td>
    </tr>`;
  }).join("");

  const relationRows = Object.entries(d.edges_by_relation || {}).slice(0, 8).map(([relation, n]) => {
    const ratio = pct(n, counts.edges);
    return `<tr>
      <td>${esc(relation)}</td>
      <td>${(n || 0).toLocaleString()}</td>
      <td style="width:220px"><div class="bar"><i style="width:${Math.min(100, ratio * 6)}%"></i></div></td>
    </tr>`;
  }).join("");

  const limitations = (caps.limitations || []).map((item) => `<li>${esc(item)}</li>`).join("");

  box.innerHTML = `
    <div class="verdict ${healthy ? "" : "warn"}">
      <div>
        <div class="kpi-num ${healthy ? "good" : "bad-txt"}" style="font-size:22px">
          ${healthy ? "这套数据当前可判为「合格」" : "体检发现问题，暂不可判为合格"}
        </div>
        <div class="chips" style="margin-top:4px">
          判据：契约检查全过 + 与 REST 快照对账全过 + 声明与实际一致，三条同时成立
        </div>
      </div>
      <div class="verdict-badges">
        <span class="badge ${conformance.ok ? "ok" : "bad"}">契约 ${conformancePassed}/${conformanceChecks.length}</span>
        <span class="badge ${verify.failed === 0 ? "ok" : "bad"}">对账 ${verifyPassed}/${verifyChecks.length}</span>
        <span class="badge ${caps.consistent ? "ok" : "bad"}">声明一致</span>
      </div>
    </div>

    <div class="kpis">
      <div class="kpi"><div class="num">${(counts.genes || 0).toLocaleString()}</div><div class="lbl">基因</div></div>
      <div class="kpi"><div class="num">${(counts.transcripts || 0).toLocaleString()}</div><div class="lbl">转录本</div></div>
      <div class="kpi"><div class="num">${(counts.proteins || 0).toLocaleString()}</div><div class="lbl">蛋白</div></div>
      <div class="kpi"><div class="num">${(counts.edges || 0).toLocaleString()}</div><div class="lbl">图边</div></div>
      <div class="kpi"><div class="num">${((caps.crosswalk || {}).edges || 0).toLocaleString()}</div><div class="lbl">跨库映射</div></div>
    </div>

    <div class="panel">
      <div class="panel-head"><h2>快照身份</h2>
        <span class="chips">${esc(caps.contract_version)} · 构建耗时 ${esc(d.elapsed_seconds ?? "—")}s</span></div>
      <table class="kv">
        <tr><td>release / assembly</td><td><code>${esc(caps.release)}</code> · ${esc(caps.assembly)}</td></tr>
        <tr><td>覆盖范围</td><td>chr${esc((caps.scope || {}).region)} · taxon ${esc((caps.scope || {}).taxon_id)} · ${esc((caps.scope || {}).description || "")}</td></tr>
        <tr><td>build_id / built_at</td><td><code>${esc(snap.build_id)}</code> · ${esc(snap.built_at)}</td></tr>
        <tr><td>内容指纹</td><td><code>${esc((d.fingerprint || "").slice(0, 16))}</code>（verify.py 会重算并与 build_stats 比对）</td></tr>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>数据来源与取回时间</h2>
        <span class="chips">${provStamped}/${provEntries.length} 个输入有明确取回时间</span></div>
      <table class="data">
        <thead><tr><th>输入文件</th><th>sha256</th><th>retrieved_at</th><th>可信度</th></tr></thead>
        <tbody>${provRows}</tbody>
      </table>
      <p class="chips" style="margin:10px 0 0">
        <code>download</code> = 下载时记录；<code>file_mtime</code> = 用文件自身时间（历史文件）；
        <code>not_reported</code> = 未知，不编造。
      </p>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>覆盖度：有多少实体真的连到了外部库</h2></div>
      <table class="data">
        <thead><tr><th>指标</th><th>数量</th><th>占比</th><th></th></tr></thead>
        <tbody>
          ${bar(coverage.proteins_uniprot, counts.proteins, "蛋白 → UniProtKB", `${coverage.proteins_uniprot}/${counts.proteins}`)}
          ${bar(coverage.proteins_pdb, counts.proteins, "蛋白 → PDB 结构", `${coverage.proteins_pdb}/${counts.proteins}`)}
          ${bar(coverage.transcripts_ccds, counts.transcripts, "转录本 → CCDS", `${coverage.transcripts_ccds}/${counts.transcripts}`)}
          ${bar(coverage.genes_hgnc, counts.genes, "基因 → HGNC", `${coverage.genes_hgnc}/${counts.genes}`)}
          ${bar(coverage.genes_entrez, counts.genes, "基因 → EntrezGene", `${coverage.genes_entrez}/${counts.genes}`)}
          ${bar(coverage.genes_mim, counts.genes, "基因 → MIM", `${coverage.genes_mim}/${counts.genes}`)}
        </tbody>
      </table>
      <p class="chips" style="margin:10px 0 0">覆盖度是「这个库值不值得用」的第一道判断：低于某个阈值时，缺的不是检索能力，而是数据。</p>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>翻译一致性（CDS → 蛋白）</h2>
        <span class="chips">共 ${(counts.transcripts || 0).toLocaleString()} 个转录本</span></div>
      <table class="data">
        <thead><tr><th>结果</th><th>数量</th><th>说明</th></tr></thead>
        <tbody>
          <tr><td>exact</td><td>${(validation.exact || 0).toLocaleString()}</td><td style="font-family:inherit">本地重建 CDS 与发布肽段逐字节一致</td></tr>
          <tr><td>alt_start</td><td>${(validation.alt_start || 0).toLocaleString()}</td><td style="font-family:inherit">注释起始密码子非标准，属已知例外</td></tr>
          <tr><td>stop_as_X</td><td>${(validation.stop_as_X || 0).toLocaleString()}</td><td style="font-family:inherit">以终止符结尾、被记为 X 的个例</td></tr>
          <tr><td>non_coding</td><td>${(validation.non_coding || 0).toLocaleString()}</td><td style="font-family:inherit">非编码转录本，无 CDS 可对</td></tr>
        </tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>与外部快照对账（${verifyPassed}/${verifyChecks.length}）</h2>
        <span class="chips">verify.py：序列、坐标、版本、指纹逐项比对冻结的 REST 响应</span></div>
      <table class="data">
        <thead><tr><th>检查项</th><th>结果</th><th>细节</th></tr></thead>
        <tbody>${verifyRows}</tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>契约检查分组（${conformancePassed}/${conformanceChecks.length}）</h2>
        <span class="chips">${esc(conformance.contract_version || "")}</span></div>
      <table class="data">
        <thead><tr><th>分组</th><th>通过</th><th>状态</th></tr></thead>
        <tbody>${groupRows}</tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>图规模与关系分布</h2>
        <span class="chips">节点 ${(counts.nodes || 0).toLocaleString()} · 边 ${(counts.edges || 0).toLocaleString()}</span></div>
      <table class="data">
        <thead><tr><th>关系</th><th>边数</th><th></th></tr></thead>
        <tbody>${relationRows}</tbody>
      </table>
    </div>

    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>已知边界（写进 capabilities，不藏着）</h2></div>
      <ul class="limit-list">${limitations}</ul>
      <p class="chips" style="margin:10px 0 0">
        这些不是 bug，而是这份切片的声明范围：越界查询会返回 <code>na</code>，而不是编一个答案。
      </p>
    </div>`;
}

/* ------------------------------------------------------------------ */
/* eval                                                                */
/* ------------------------------------------------------------------ */
async function renderEval() {
  const box = $("#eval-out");
  box.innerHTML = `<div class="empty"><h3>加载中…</h3></div>`;
  const payload = await api("bench");
  const bench = payload.bench || { total: {}, by_category: {}, results: [] };
  const conformance = payload.conformance || {};
  const verify = payload.verify || {};
  const total = bench.total || {};

  const kpi = (value, label, cls = "") =>
    `<div class="kpi"><div class="num ${cls}">${esc(value)}</div><div class="lbl">${esc(label)}</div></div>`;

  const categories = Object.entries(bench.by_category || {}).map(([name, bucket]) => {
    const ratio = bucket.n ? bucket.value / bucket.n : 0;
    return `<tr>
      <td>${esc(name)}</td>
      <td>${bucket.value}/${bucket.n}</td>
      <td>${bucket.source}/${bucket.n}</td>
      <td>${bucket.evidence_scored ? `${bucket.evidence}/${bucket.evidence_scored}` : "n/a"}</td>
      <td style="width:220px"><div class="bar"><i style="width:${(ratio * 100).toFixed(0)}%"></i></div></td>
    </tr>`;
  }).join("");

  const items = (bench.results || []).map((result) => `
    <tr>
      <td>${esc(result.id)}</td>
      <td>${esc(result.category)}</td>
      <td style="font-family:inherit">${esc(result.query)}</td>
      <td>${result.value_ok ? '<span class="good-txt">✓</span>' : '<span class="bad-txt">✗</span>'}</td>
      <td>${result.source_ok ? '<span class="good-txt">✓</span>' : '<span class="bad-txt">✗</span>'}</td>
      <td>${result.evidence_scored
        ? (result.evidence_ok ? '<span class="good-txt">✓</span>' : '<span class="bad-txt">✗</span>')
        : '<span class="hint">n/a</span>'}</td>
    </tr>`).join("");

  box.innerHTML = `
    <div class="kpis">
      ${kpi(`${total.value ?? 0}/${total.items ?? 0}`, "准 · 值对不对", "good")}
      ${kpi(`${total.source ?? 0}/${total.items ?? 0}`, "诚 · 有没有来源", "good")}
      ${kpi(`${total.evidence ?? 0}/${total.evidence_scored ?? 0}`, "寸 · 空值与证据类型", "good")}
      ${kpi(`${conformance.passed ?? 0}/${(conformance.checks || []).length}`, "契约检查", "good")}
      ${kpi(`${verify.passed ?? 0}/${(verify.checks || []).length}`, "与 REST 快照对账", "")}
    </div>
    <div class="panel">
      <div class="panel-head"><h2>分类得分</h2>
        <span class="chips">寸 只统计真正考到的条目，其余标 n/a</span></div>
      <table class="data">
        <thead><tr><th>类别</th><th>准</th><th>诚</th><th>寸</th><th>进度</th></tr></thead>
        <tbody>${categories}</tbody>
      </table>
    </div>
    <div class="panel" style="margin-top:16px">
      <div class="panel-head"><h2>逐条结果（${bench.results ? bench.results.length : 0} 条）</h2></div>
      <table class="data">
        <thead><tr><th>ID</th><th>类别</th><th>问题</th><th>准</th><th>诚</th><th>寸</th></tr></thead>
        <tbody>${items}</tbody>
      </table>
    </div>
    <div class="note">
      本套题全部由冻结数据推导（每题带 <code>source_of_truth</code>），不含人工标注。
      用途是回归保护与口径锁定：改代码后能立刻发现哪里坏了，并保证"准/诚/寸"的判定标准不变。
      它不对外部准确率背书——同一份数据既出题又答题，只能证明内部一致。
    </div>`;
}

/* ------------------------------------------------------------------ */
/* boot: ?q=YWHAB deep-links straight into a search                     */
/* ------------------------------------------------------------------ */
(async function start() {
  await boot();
  fillTargets();
  const params = new URLSearchParams(location.search);
  const query = params.get("q") || (location.hash ? null : "YWHAB");
  if (query) {
    $("#q").value = query;
    $("#type").value = params.get("type") || "";
    doSearch(query, params.get("type") || "");
  }
  const hash = (location.hash || "").replace("#", "");
  if (hash && ["seq", "cross", "eval", "health", "bridge"].includes(hash)) {
    document.querySelector(`#tabs button[data-view="${hash}"]`).click();
  }
})();
