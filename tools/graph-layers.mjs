// 依存グラフを **層状のボックス図** として書き出す(自己完結 HTML、外部ライブラリ無し)。
//
// ★このレイアウトが答える問い: 「どこから手を付ければよいか」
//   右端 = [IUTchIII] Corollary 3.12。左へ行くほど「先に要るもの」。
//   **左端の層は依存を持たない**ので、そこから始められる。
//
// ★循環の扱い(ここが設計の核心)
//   循環があると層に並べられないので、**強連結成分(SCC)に潰す**。
//   潰した結果が凝縮 DAG で、これは必ず層に並ぶ。
//   ★そして SCC は「まとめてしか扱えない塊」なので、**そのままグルーピングの単位**になる。
//
//   ★★訂正の記録(2026-08-15): 当初この欄には
//   「824 節点 → 最大 SCC 262 節点(IUTchI/II/III/IV にまたがる)。
//    IUT 本体は1つの循環であり、『IUTchI から順に』という順序付けは存在しない」
//   と書いていた。**誤りである。**
//   その 262 節点の循環は、辺の定義が「本文が名指しした」だったことの副作用だった。
//   辺を「証明が依拠している」に定め直した実測(ResearchPaper/cycle-analysis.md):
//     節点 824 → 653 / 最大 SCC 262 → **17** / 循環する成分 30 → 12 / 層 19 → **55**
//   ★残る最大 SCC 17 は IUTchIII + IUTchIV のみで、**IUTchI・IUTchII は循環から外れる。**
//   すなわち **IUTchI → II → III → IV の順序付けは存在する。**
//
// ★ボックスの決め方
//   ・サイズ>1 の SCC は、それ自体で1ボックス(分割できないので)
//   ・残りの単独節点は (層, 論文) でまとめて1ボックス
//     = 「依存関係が少ないものはグルーピング」
//
// ★分母の限界: 原文側は pdftotext 由来で目視していない。番号の無い依存を数え落とし、
//   「the discussion of ...」型の案内を数え過ぎる。**下界でも上界でもない。**
//
// 使い方: node tools/graph-layers.mjs [出力先.html]

import { readFileSync, writeFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const LEAN = join(ROOT, 'lean', 'ABC3');
const OUT = process.argv[2] ?? join(ROOT, 'dependency-graph.html');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const FILE_OF = JSON.parse(readFileSync(join(ROOT, 'tools', '_fileof.json'), 'utf8'));

// ── 原文側の抽出 ──────────────────────────────────────────
const loaded = new Map();
function load(tag) {
  if (loaded.has(tag)) return loaded.get(tag);
  const f = FILE_OF[tag]; const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) { loaded.set(tag, null); return null; }
  const lines = readFileSync(p, 'utf8').split(/\r?\n/);
  const pageOf = new Array(lines.length).fill(0);
  { let pg = 0; for (let i = 0; i < lines.length; i++) {
      const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  for (let i = 0; i < lines.length; i++) { const m = lines[i].match(declRe); if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] }); }
  for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;
  const decls = new Map(); for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  const o = { decls, lines }; loaded.set(tag, o); return o;
}
// ★辺の意味(2026-08-15 に定めた): 「原文の証明が実際に依拠しているもの」。
//   次は辺ではない —— 前方参照(cf. …, below)/ 解説への案内(the discussion of …)。
//   実測(ResearchPaper/cycle-analysis.md): この2つを落とすだけで
//   最大 SCC が 262 → 16 になり、**論文をまたぐ循環が消える**。
//   ★Remark を種別で落とすことはしない。[IUTchIII] Remark 3.9.5, (i) は
//    正則包を**定義**しており、我々は実際にそれに依拠している。
const NOT_A_DEP = /\bbelow\b|the discussion of|the discussion surrounding|\bsee also\b/i;

function outs(tag, name) {
  const P = load(tag); if (!P) return null;
  const d = P.decls.get(name); if (!d) return null;
  const body = P.lines.slice(d.line, d.end).join('\n');
  const ctx = (i, len) => body.slice(Math.max(0, i - 90), i + len + 60).replace(/\s+/g, ' ');
  const es = new Map(); const spans = [];
  const keep = (to, c) => {
    const soft = NOT_A_DEP.test(c);
    const cur = es.get(to);
    // 同じ先が複数回出るなら、1回でも「依拠」の文脈があれば辺として残す
    if (cur === undefined || (cur && !soft)) es.set(to, soft);
  };
  const xre = new RegExp(`\\[([A-Za-z]+)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(xre)) { keep(`${m[1]} / ${m[2]} ${m[3]}`, ctx(m.index, m[0].length)); spans.push([m.index, m.index + m[0].length]); }
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to !== name && P.decls.has(to)) keep(`${tag} / ${to}`, ctx(m.index, m[0].length));
  }
  return { edges: [...es].filter(([, soft]) => !soft).map(([k]) => k), page: d.page,
           dropped: [...es].filter(([, soft]) => soft).length };
}

const ROOTK = 'IUTchIII / Corollary 3.12';
const adj = new Map(), page = new Map();
{
  const q = [ROOTK], seen = new Set();
  while (q.length) {
    const k = q.shift(); if (seen.has(k)) continue; seen.add(k);
    const [t, n] = k.split(' / ');
    const r = outs(t, n);
    adj.set(k, r ? r.edges : []); page.set(k, r ? r.page : 0);
    for (const e of (r ? r.edges : [])) if (!seen.has(e)) q.push(e);
  }
  for (const k of seen) if (!adj.has(k)) { adj.set(k, []); page.set(k, 0); }
}

// ── Tarjan SCC(反復版) ────────────────────────────────────
const index = new Map(), low = new Map(), onstk = new Set(), stk = [], comp = new Map();
let idx = 0, nc = 0;
for (const s of adj.keys()) {
  if (index.has(s)) continue;
  const work = [[s, 0]];
  while (work.length) {
    const fr = work[work.length - 1], v = fr[0];
    if (fr[1] === 0) { index.set(v, idx); low.set(v, idx); idx++; stk.push(v); onstk.add(v); }
    let rec = false;
    const es = adj.get(v) ?? [];
    for (let i = fr[1]; i < es.length; i++) {
      const w = es[i]; if (!adj.has(w)) continue;
      if (!index.has(w)) { fr[1] = i + 1; work.push([w, 0]); rec = true; break; }
      else if (onstk.has(w)) low.set(v, Math.min(low.get(v), index.get(w)));
    }
    if (rec) continue;
    if (low.get(v) === index.get(v)) { const c = nc++; while (true) { const w = stk.pop(); onstk.delete(w); comp.set(w, c); if (w === v) break; } }
    work.pop();
    if (work.length) { const p = work[work.length - 1][0]; low.set(p, Math.min(low.get(p), low.get(v))); }
  }
}
const members = new Map();
for (const [k, c] of comp) { if (!members.has(c)) members.set(c, []); members.get(c).push(k); }

// 凝縮 DAG と層(左 0 = 依存なし)
const cadj = new Map();
for (const [v, es] of adj) {
  const cv = comp.get(v); if (!cadj.has(cv)) cadj.set(cv, new Set());
  for (const w of es) { const cw = comp.get(w); if (cw !== undefined && cw !== cv) cadj.get(cv).add(cw); }
}
for (const c of members.keys()) if (!cadj.has(c)) cadj.set(c, new Set());
const memo = new Map();
const layerOf = (c, st = new Set()) => {
  if (memo.has(c)) return memo.get(c);
  if (st.has(c)) return 0;
  st.add(c); let d = 0;
  for (const w of cadj.get(c) ?? []) d = Math.max(d, layerOf(w, st) + 1);
  st.delete(c); memo.set(c, d); return d;
};

// ── 我々の位置(3層) ───────────────────────────────────────
function walk(d, a = []) { for (const f of readdirSync(d)) { const p = join(d, f); if (statSync(p).isDirectory()) walk(p, a); else if (p.endsWith('.lean')) a.push(p); } return a; }
const ours = new Map();
{
  const SRC_RE = /paper\s*:=\s*"([^"]*)"[\s\S]{0,400}?item\s*:=\s*"([^"]*)"/g;
  const EDGE_RE = /\.otherPaper\s+"\[?([A-Za-z]+)\]?"\s+"([^"]*)"/g;
  const ITEM_RE = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`);
  const rank = { named: 0, skeleton: 1, landed: 2 };
  const put = (k, kind, file) => { const c = ours.get(k); if (!c || rank[kind] > rank[c.kind]) ours.set(k, { kind, file }); };
  for (const f of walk(LEAN)) {
    const rel = f.slice(ROOT.length + 1).replace(/\\/g, '/'), t = readFileSync(f, 'utf8');
    const landed = /\.inProject|\.inMathlib/.test(t);
    for (const m of t.matchAll(SRC_RE)) { const im = ITEM_RE.exec(m[2]); if (im) put(`${m[1]} / ${im[1]} ${im[2]}`, landed ? 'landed' : 'skeleton', rel); }
    for (const m of t.matchAll(EDGE_RE)) { const im = ITEM_RE.exec(m[2]); if (im) put(`${m[1]} / ${im[1]} ${im[2]}`, 'named', rel); }
  }
}

// ── ボックス化: SCC(>1) はそのまま / 単独節点は (層, 論文) でまとめる ──
const boxes = new Map();  // id -> {layer, label, tags, items[], scc}
const boxOfNode = new Map();
for (const [c, mem] of members) {
  const L = layerOf(c);
  if (mem.length > 1) {
    const id = `scc${c}`;
    const tg = new Map(); for (const m of mem) { const t = m.split(' / ')[0]; tg.set(t, (tg.get(t) ?? 0) + 1); }
    boxes.set(id, { layer: L, tags: [...tg].sort((a, b) => b[1] - a[1]), items: mem, scc: true });
    for (const m of mem) boxOfNode.set(m, id);
  } else {
    const k = mem[0], t = k.split(' / ')[0];
    const id = `L${L}:${t}`;
    if (!boxes.has(id)) boxes.set(id, { layer: L, tags: [[t, 0]], items: [], scc: false });
    boxes.get(id).items.push(k);
    boxes.get(id).tags[0][1]++;
    boxOfNode.set(k, id);
  }
}
// ボックス間の辺
const bedge = new Map();
for (const [v, es] of adj) {
  const a = boxOfNode.get(v); if (!a) continue;
  for (const w of es) { const b = boxOfNode.get(w); if (!b || a === b) continue; const k = `${a} ${b}`; bedge.set(k, (bedge.get(k) ?? 0) + 1); }
}

// ── 配置 ─────────────────────────────────────────────────
const BW = 150, GAPX = 58, GAPY = 10, PADT = 62;
const maxL = Math.max(...[...boxes.values()].map((b) => b.layer));
const byLayer = new Map();
for (const [id, b] of boxes) { if (!byLayer.has(b.layer)) byLayer.set(b.layer, []); byLayer.get(b.layer).push(id); }
let maxH = 0;
for (const [L, ids] of byLayer) {
  ids.sort((a, b) => boxes.get(b).items.length - boxes.get(a).items.length);
  let y = PADT;
  for (const id of ids) {
    const b = boxes.get(id);
    b.h = Math.max(28, Math.min(170, 24 + b.items.length * 1.9));
    b.w = BW;
    // ★層 0(依存なし)が最も左、根(最大層)が最も右。
    //   以前 (maxL - b.layer) と書いており **左右が反転していた**(2026-08-15 訂正)。
    b.x = b.layer * (BW + GAPX);
    b.y = y; y += b.h + GAPY;
  }
  maxH = Math.max(maxH, y);
}
// x を反転(層 0 = 左)。上で (maxL - layer) にしているので根が x 最大 = 右。
const B = [...boxes.entries()].map(([id, b]) => ({
  id, layer: b.layer, x: Math.round(b.x), y: Math.round(b.y), w: b.w, h: Math.round(b.h),
  tags: b.tags, n: b.items.length, scc: b.scc ? 1 : 0,
  ours: b.items.filter((k) => ours.has(k)).length,
  landed: b.items.filter((k) => ours.get(k)?.kind === 'landed').length,
  root: b.items.includes(ROOTK) ? 1 : 0,
  items: b.items.slice(0, 400).map((k) => ({ k, p: page.get(k) ?? 0, o: ours.get(k)?.kind ?? null })),
}));
const bidx = new Map(B.map((b, i) => [b.id, i]));
const BE = [...bedge.entries()].map(([k, w]) => { const [a, b] = k.split(' '); return [bidx.get(a), bidx.get(b), w]; })
  .filter(([a, b]) => a !== undefined && b !== undefined);

const layerStats = [...byLayer.keys()].sort((a, b) => a - b).map((L) => ({
  L, boxes: byLayer.get(L).length,
  nodes: byLayer.get(L).reduce((s, id) => s + boxes.get(id).items.length, 0),
  ours: byLayer.get(L).reduce((s, id) => s + boxes.get(id).items.filter((k) => ours.has(k)).length, 0),
}));

const data = {
  B, BE, layerStats, maxL,
  stats: { nodes: adj.size, sccs: nc, bigSccs: [...members.values()].filter((m) => m.length > 1).length,
    boxes: B.length, layers: maxL + 1, ours: [...ours.keys()].filter((k) => adj.has(k)).length },
  width: (maxL + 1) * (BW + GAPX), height: maxH + 40,
};
const html = readFileSync(join(ROOT, 'tools', 'graph-layers.template.html'), 'utf8').replace('/*__DATA__*/', JSON.stringify(data));
writeFileSync(OUT, html, 'utf8');
console.log(`書き出し: ${OUT}`);
console.log(`  節点 ${data.stats.nodes} / SCC ${data.stats.sccs}(サイズ>1 は ${data.stats.bigSccs})`);
console.log(`  ボックス ${data.stats.boxes} / 層 ${data.stats.layers}(左 0 = 依存なし、右 ${maxL} = 根)`);
console.log(`  我々が触れている節点 ${data.stats.ours}`);
