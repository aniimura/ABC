// 論文をまたぐ辺と同一論文内の辺を **両方** 辿った推移閉包を測る。
//
// ★この道具の位置づけ(先に書く)
//   `check.mjs` の依存グラフは我々が **`.needs` に書き出した** 辺だけを見る。
//   こちらは **原文の .txt から機械抽出した** 辺を見る。前者は我々の作業量、
//   後者は本来必要な量であり、**その差が残りの仕事**である。
//   ★`pdftotext` 経由なので目視していない。事実2により項目名は壊れうる。
//   **逐語には使えない。規模の当たりを付けるためだけのもの。下界である。**
//
// 使い方: node tools/full-graph.mjs [<タグ> <項目>]   既定は IUTchIII / "Corollary 3.12"

import { readFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';

// 引用タグ → 0_Source のファイル名(拡張子なし)。
// ★papers.json に登記されていない論文も辺の先に出るので、ここで補う。
//   登記の要否は check.mjs の D26 が別途落とす。
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
  // ★以下は略号の慣行から当てたもので、**目視で確認していない**。
  //   誤った写像は誤った辺を生むので、この行の下の数字は上の12件より弱い。
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

const loaded = new Map(); // tag -> { decls: Map(name -> {line,end,page}), lines }

function load(tag) {
  if (loaded.has(tag)) return loaded.get(tag);
  const f = FILE_OF[tag];
  const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) { loaded.set(tag, null); return null; }

  const lines = readFileSync(p, 'utf8').split(/\r?\n/);   // ★CRLF を落とす
  const pageOf = new Array(lines.length).fill(0);
  { let pg = 0; for (let i = 0; i < lines.length; i++) {
      const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }

  // ★宣言行 = 行頭 + 番号 + ピリオド + (行末 または 括弧つき表題)。
  //   緩めると継続行("Definition 3.1. Let")を宣言と誤検出して本体が切れる(実測)。
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(declRe);
    if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] });
  }
  for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;

  const decls = new Map();
  for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  const out = { decls, lines };
  loaded.set(tag, out);
  return out;
}

// 節点 = "TAG / Kind N.M"
function outEdges(tag, name) {
  const P = load(tag);
  if (!P) return { edges: [], unknownPaper: true };
  const d = P.decls.get(name);
  if (!d) return { edges: [], undeclared: true };
  const body = P.lines.slice(d.line, d.end).join('\n');

  const edges = [];
  // ① 論文をまたぐ辺: "[Tag], Kind N.M" / "[Tag] Kind N.M"
  // ★タグに数字を許すこと(2026-08-15 修正)。IUT 系は `[IUTchIII]` の形だが、
  //   **FrdI は `[Mzk15]` の形で引用する**。`[A-Za-z]+` だとこの引用が masking されず、
  //   直後の ② が `Definition 3.1` を**同一論文内の参照として拾ってしまう**。
  //   実測: FrdI に偽の辺 `Definition 1.2 → Definition 3.1` が生まれ、
  //   §1 の 8 項目 + §3 の 3 項目が **11 節点の偽の循環**に見えていた。
  const xre = new RegExp(`\\[([A-Za-z][A-Za-z0-9]*)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  const spans = [];
  for (const m of body.matchAll(xre)) {
    edges.push({ tag: m[1], name: `${m[2]} ${m[3]}` });
    spans.push([m.index, m.index + m[0].length]);
  }
  // ② 同一論文内の辺: ① に食われていない "Kind N.M"
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to === name) continue;
    if (P.decls.has(to)) edges.push({ tag, name: to });
  }
  return { edges };
}

const rootTag = process.argv[2] ?? 'IUTchIII';
const rootName = process.argv[3] ?? 'Corollary 3.12';
const key = (t, n) => `${t} / ${n}`;

const depth = new Map([[key(rootTag, rootName), 0]]);
const seen = new Set();
const perPaper = new Map();
const unknownPapers = new Set();
const undeclared = new Set();
let edgeCount = 0;
const queue = [{ tag: rootTag, name: rootName }];

while (queue.length) {
  const cur = queue.shift();
  const k = key(cur.tag, cur.name);
  if (seen.has(k)) continue;
  seen.add(k);
  perPaper.set(cur.tag, (perPaper.get(cur.tag) ?? 0) + 1);

  const r = outEdges(cur.tag, cur.name);
  if (r.unknownPaper) { unknownPapers.add(cur.tag); continue; }
  if (r.undeclared) { undeclared.add(k); continue; }
  for (const e of r.edges) {
    edgeCount++;
    const ek = key(e.tag, e.name);
    if (!depth.has(ek)) { depth.set(ek, depth.get(k) + 1); queue.push(e); }
  }
}

console.log(`根: ${key(rootTag, rootName)}`);
console.log(`到達した節点: ${seen.size} 件 / 辺 ${edgeCount} 本 / 最大深さ ${Math.max(...[...seen].map((k) => depth.get(k) ?? 0))}`);
console.log('論文別:');
for (const [t, n] of [...perPaper].sort((a, b) => b[1] - a[1])) {
  console.log(`  ${String(n).padStart(4)}  ${t}${FILE_OF[t] ? '' : '  ← 0_Source に対応不明'}`);
}
if (unknownPapers.size) console.log(`★.txt が見つからない論文: ${[...unknownPapers].join(', ')}`);
console.log(`★節点として辿れなかった参照(その論文に宣言が無い): ${undeclared.size} 件`);
for (const u of [...undeclared].slice(0, 8)) console.log(`     ${u}`);
