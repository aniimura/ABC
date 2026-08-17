// 依存グラフを自己完結の HTML(外部ライブラリ無し)に書き出す。
//
// ★何を描くか
//   下地  = 原文側のグラフ(tools/full-graph.mjs と同じ抽出。824 節点)
//   重ね  = **我々が書き出した節点**(Lean の `.src` / `.needs` から)
//   すなわち「地形」と「我々がどこまで来たか」を1枚に重ねる。被覆率が目で見える。
//
// ★分母の限界(HTML 側にも出す)
//   原文側は pdftotext 由来で目視していない。番号の無い依存を数え落とし、
//   「the discussion of ...」型の案内を数え過ぎる。**下界でも上界でもない。**
//
// 使い方: node tools/graph-html.mjs [出力先.html]

import { readFileSync, writeFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const LEAN = join(ROOT, 'lean', 'ABC3');
const OUT = process.argv[2] ?? join(ROOT, 'dependency-graph.html');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';

const FILE_OF = {
  IUTchI: 'Inter-universal Teichmuller Theory I',
  IUTchII: 'Inter-universal Teichmuller Theory II',
  IUTchIII: 'Inter-universal Teichmuller Theory III',
  IUTchIV: 'Inter-universal Teichmuller Theory IV',
  AbsTopI: 'Topics in Absolute Anabelian Geometry I',
  AbsTopII: 'Topics in Absolute Anabelian Geometry II',
  AbsTopIII: 'Topics in Absolute Anabelian Geometry III',
  EtTh: 'The Etale Theta Function and its Frobenioid-theoretic Manifestations',
  FrdI: 'The Geometry of Frobenioids I',
  FrdII: 'The Geometry of Frobenioids II',
  SemiAnbd: 'Semi-graphs of Anabelioids',
  pGC: 'A Version of the Grothendieck Conjecture for p-adic Local Fields',
  GenEll: 'Arithmetic Elliptic Curves in General Position',
  AbsAnab: 'Absolute Anabelian Geometry',
  CanLift: 'Canonical Liftings',
  Config: 'The Algebraic and Anabelian Geometry of Configuration Spaces',
  CombGC: 'Combinatorial Grothendieck Conjecture',
  CombCusp: 'Combinatorial Cuspidalization',
  AbsSect: 'Galois Sections',
  MT: 'An Introduction to p-adic Teichmuller Theory',
  HASurII: 'A Survey of the Hodge-Arakelov Theory of Elliptic Curves II',
  NodNon: 'Combinatorial Anabelian Topics I',
};

// ── 原文側の抽出(full-graph.mjs と同じ規則) ─────────────────────────
const loaded = new Map();
function load(tag) {
  if (loaded.has(tag)) return loaded.get(tag);
  const f = FILE_OF[tag];
  const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) { loaded.set(tag, null); return null; }
  const lines = readFileSync(p, 'utf8').split(/\r?\n/);
  const pageOf = new Array(lines.length).fill(0);
  { let pg = 0; for (let i = 0; i < lines.length; i++) {
      const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(declRe);
    if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] });
  }
  for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;
  const decls = new Map();
  for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  const o = { decls, lines };
  loaded.set(tag, o); return o;
}
function outEdges(tag, name) {
  const P = load(tag); if (!P) return null;
  const d = P.decls.get(name); if (!d) return null;
  const body = P.lines.slice(d.line, d.end).join('\n');
  const es = []; const spans = [];
  // ★タグに数字を許すこと(2026-08-15 修正)。FrdI は [Mzk15] の形で引用するので、
  //   [A-Za-z]+ だと masking し損ね、直後の同一論文内スキャンが偽の辺を作る。
  //   実測: FrdI に偽の辺 Definition 1.2 → Definition 3.1 が生まれ、偽の循環に見えていた。
  const xre = new RegExp(`\\[([A-Za-z][A-Za-z0-9]*)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(xre)) { es.push([m[1], `${m[2]} ${m[3]}`]); spans.push([m.index, m.index + m[0].length]); }
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to !== name && P.decls.has(to)) es.push([tag, to]);
  }
  return { edges: es, page: d.page };
}

const rootTag = 'IUTchIII', rootName = 'Corollary 3.12';
const key = (t, n) => `${t} / ${n}`;
const nodes = new Map();   // key -> {tag,name,depth,page,external}
const edges = [];
{
  const q = [[rootTag, rootName, 0]];
  const seen = new Set();
  while (q.length) {
    const [t, n, dp] = q.shift();
    const k = key(t, n);
    if (seen.has(k)) continue;
    seen.add(k);
    const r = outEdges(t, n);
    nodes.set(k, { tag: t, name: n, depth: dp, page: r?.page ?? 0, external: !FILE_OF[t] });
    if (!r) continue;
    for (const [et, en] of r.edges) {
      const ek = key(et, en);
      edges.push([k, ek]);
      if (!seen.has(ek) && !q.some(([a, b]) => a === et && b === en)) q.push([et, en, dp + 1]);
    }
  }
}
// 到達先が節点として作られていない辺(外部文献など)も節点にする
for (const [, to] of edges) if (!nodes.has(to)) {
  const [t, n] = to.split(' / ');
  nodes.set(to, { tag: t, name: n, depth: 99, page: 0, external: !FILE_OF[t] });
}

// ── 我々の側: Lean の `.src` から「どの項目を張ったか」を拾う ────────
function walk(d, a = []) {
  for (const f of readdirSync(d)) {
    const p = join(d, f);
    if (statSync(p).isDirectory()) walk(p, a); else if (p.endsWith('.lean')) a.push(p);
  }
  return a;
}
// 3層に分ける。層が違えば意味が違うので、混ぜて1つの率にしない。
//   landed   … `.needs` の LeanStatus が `inProject`/`inMathlib`。実物に着地している
//   skeleton … `.src` を持つ節点を張った(statement は書いた)
//   named    … `.needs` の辺の先として**名前だけ**知っている(まだ張っていない)
const ours = new Map();  // key -> {file, kind}
{
  const SRC_RE = /paper\s*:=\s*"([^"]*)"[\s\S]{0,400}?item\s*:=\s*"([^"]*)"/g;
  const EDGE_RE = /\.otherPaper\s+"\[?([A-Za-z]+)\]?"\s+"([^"]*)"/g;
  const ITEM_RE = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`);
  const rank = { named: 0, skeleton: 1, landed: 2 };
  const put = (k, kind, file) => {
    const cur = ours.get(k);
    if (!cur || rank[kind] > rank[cur.kind]) ours.set(k, { kind, file });
  };
  for (const f of walk(LEAN)) {
    const rel = f.slice(ROOT.length + 1).replace(/\\/g, '/');
    const t = readFileSync(f, 'utf8');
    // ① 張った節点
    for (const m of t.matchAll(SRC_RE)) {
      const im = ITEM_RE.exec(m[2]);
      if (im) put(key(m[1], `${im[1]} ${im[2]}`), 'skeleton', rel);
    }
    // ② 辺の先として名前を知っている / 着地している
    for (const m of t.matchAll(EDGE_RE)) {
      const im = ITEM_RE.exec(m[2]);
      if (im) put(key(m[1], `${im[1]} ${im[2]}`), 'named', rel);
    }
  }
  // ③ 着地: `.needs` に inProject/inMathlib を持つ節点(ファイル単位で近似)
  for (const f of walk(LEAN)) {
    const t = readFileSync(f, 'utf8');
    if (!/\.inProject|\.inMathlib/.test(t)) continue;
    const rel = f.slice(ROOT.length + 1).replace(/\\/g, '/');
    for (const m of t.matchAll(SRC_RE)) {
      const im = ITEM_RE.exec(m[2]);
      if (im) put(key(m[1], `${im[1]} ${im[2]}`), 'landed', rel);
    }
  }
}

// ── Interface の義務(構造体)を節点として足す ──────────────────────
// ★なぜ足すか: グラフの節点語彙は**原文の項目番号**なので、
//   `GenEll Definition 1.1` に 9 本の obligation がぶら下がっていても
//   **1 節点に潰れて見えなかった**(2026-08-17 実測)。
//   check.mjs の「Interface 実装待ち」と同じ単位を図にも出す。
// ★被覆率は汚さない —— 義務節点は `ob:1` を持ち、既存の統計から除く。
{
  const witness = new Set();
  for (const f of walk(LEAN)) {
    const t = readFileSync(f, 'utf8');
    for (const m of t.matchAll(/^theorem\s+([A-Za-z_][\w'.]*)\.nonvacuous/gm))
      witness.add(m[1].split('.').pop());
  }
  for (const f of walk(LEAN)) {
    const rel = f.slice(ROOT.length + 1).replace(/\\/g, '/');
    if (!/\/Interface\//.test(rel)) continue;
    const t = readFileSync(f, 'utf8');
    for (const sm of t.matchAll(/^structure\s+([A-Za-z_][A-Za-z0-9_']*)\s/gm)) {
      const nm = sm[1];
      // ★`.src` を持たない構造体もある(`.waiting` だけ)——錨を付けずに数える
      const sm2 = new RegExp(`${nm}\\.src\\s*:[^=]*:=\\s*\\{([\\s\\S]{0,400}?)\\}`).exec(t);
      const pm = sm2 ? /paper\s*:=\s*"([^"]*)"/.exec(sm2[1]) : null;
      const im = sm2 ? /item\s*:=\s*"([^"]*)"/.exec(sm2[1]) : null;
      const km = (pm && im) ? new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`).exec(im[1]) : null;
      const anchorKey = (km && pm) ? key(pm[1], `${km[1]} ${km[2]}`) : null;
      const anc = anchorKey && nodes.has(anchorKey) ? nodes.get(anchorKey) : null;
      const k = `義務 / ${nm}`;
      if (nodes.has(k)) continue;
      nodes.set(k, {
        tag: 'Interface', name: nm,
        depth: anc ? (anc.depth === 99 ? 99 : anc.depth + 1) : 99,
        page: anc ? anc.page : 0, external: false,
        oblig: true, done: witness.has(nm), file: rel, item: im ? im[1] : null,
      });
      if (anc) edges.push([anchorKey, k]);
    }
  }
}

// ── 配置: 深さで同心円、角度は論文ごとにまとめる ─────────────────
const byTag = new Map();
for (const [k, v] of nodes) { if (!byTag.has(v.tag)) byTag.set(v.tag, []); byTag.get(v.tag).push(k); }
const tags = [...byTag.keys()].sort((a, b) => byTag.get(b).length - byTag.get(a).length);
const tagAngle = new Map();
{
  const total = [...byTag.values()].reduce((s, a) => s + a.length, 0);
  let acc = 0;
  for (const t of tags) {
    const span = (byTag.get(t).length / total) * Math.PI * 2;
    tagAngle.set(t, [acc, acc + span]);
    acc += span;
  }
}
const maxDepth = Math.max(...[...nodes.values()].map((v) => (v.depth === 99 ? 0 : v.depth)));
for (const t of tags) {
  const list = byTag.get(t).sort((a, b) => nodes.get(a).depth - nodes.get(b).depth);
  const [a0, a1] = tagAngle.get(t);
  list.forEach((k, i) => {
    const v = nodes.get(k);
    const dep = v.depth === 99 ? maxDepth : v.depth;
    const r = 120 + dep * 105 + ((i * 37) % 55);
    const ang = a0 + ((i + 0.5) / list.length) * (a1 - a0);
    v.x = Math.cos(ang) * r;
    v.y = Math.sin(ang) * r;
  });
}

// ── 出力 ────────────────────────────────────────────────
const idx = new Map([...nodes.keys()].map((k, i) => [k, i]));
const N = [...nodes.entries()].map(([k, v]) => ({
  k, t: v.tag, n: v.name, d: v.depth === 99 ? maxDepth : v.depth, p: v.page,
  x: Math.round(v.x), y: Math.round(v.y),
  o: v.oblig ? (v.done ? '義務(埋まった)' : '義務(実装待ち)')
             : (ours.has(k) ? ours.get(k).kind : null),
  f: v.oblig ? v.file : (ours.has(k) ? ours.get(k).file : null),
  e: v.external ? 1 : 0,
  ob: v.oblig ? (v.done ? 2 : 1) : 0,
  it: v.oblig ? v.item : null,
}));
const E = edges.filter(([a, b]) => idx.has(a) && idx.has(b)).map(([a, b]) => [idx.get(a), idx.get(b)]);

// ★義務節点は被覆率の分母・分子の**どちらにも入れない**(語彙が違うから)
const NP = N.filter((n) => !n.ob);
const stats = {
  nodes: NP.length, edges: E.filter(([a, b]) => !N[a].ob && !N[b].ob).length, maxDepth,
  ours: NP.filter((n) => n.o).length,
  landed: NP.filter((n) => n.o === 'landed').length,
  skeleton: NP.filter((n) => n.o === 'skeleton').length,
  named: NP.filter((n) => n.o === 'named').length,
  papers: tags.filter((t) => t !== 'Interface').length,
  external: NP.filter((n) => n.e).length,
  oblig: N.filter((n) => n.ob).length,
  obligDone: N.filter((n) => n.ob === 2).length,
};

const html = readFileSync(join(ROOT, 'tools', 'graph-html.template.html'), 'utf8')
  .replace('/*__DATA__*/', JSON.stringify({ N, E, stats, tags }));
writeFileSync(OUT, html, 'utf8');
console.log(`書き出し: ${OUT}`);
console.log(`  節点 ${stats.nodes} / 辺 ${stats.edges} / 最大深さ ${stats.maxDepth} / 論文 ${stats.papers} / 外部文献 ${stats.external}`);
console.log(`  我々: 着地 ${stats.landed} / 張った ${stats.skeleton} / 名前だけ ${stats.named}`);
console.log(`  Interface の義務: ${stats.obligDone} / ${stats.oblig}(★被覆率には入れていない——語彙が違う)`);
