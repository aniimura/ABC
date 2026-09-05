#!/usr/bin/env node
/**
 * 依存グラフ —— **ノード = Lean ファイル**(2026-09-03、第 1454)
 *
 * ★当初の意図(設計者の言葉):
 *   ① 論文の主張から `Skeleton` を作る
 *   ② **Lean ファイル単位のグラフノード**が現れる
 *   ③ ノードがどの論文(または理論)に属するか判別できる
 *   ④ 形式化作業中にノードを足しながら進める
 *   ⑤ グラフで進捗を管理する
 *   ⑥ **ファイル数とノード数が一致**し、Explorer 上でも確認できる
 *
 * ★★これまでのグラフ(`tools/check.mjs` の中)は ② を満たしていなかった。
 *   節点を「条なしの原典項目」にしていたので:
 *     - `Lemma 3.5` は **626 宣言 / 135 ファイル**あるのに **1 節点**にしか見えない
 *     - 条つき `.src`(`Proposition 3.4(段 A)`)が親の項目へ丸められ、
 *       「段 A が Proposition 3.4 を引く」が**自己ループ 14 件**として現れていた
 *     - 理論(`GaloisRep` 等、983 ファイル)は**そもそも節点にならない**
 *
 * ☆本ツールはノードをファイルにする。ファイル数 = ノード数は**定義から**成り立つ(⑥)。
 *   辺は `import`(Lean が保証するので循環しない)。所属はディレクトリから読む(③)。
 *   `.src` の項目・`sorry` の有無は**ノードの属性**として持つ。
 *
 * 使い方:
 *   node tools/graph.mjs              # 要約
 *   node tools/graph.mjs --owner GenEll   # 所属で絞る
 *   node tools/graph.mjs --from Skeleton/GenEll/Section3.lean  # そのノードの下流
 *   node tools/graph.mjs --sorry      # sorry に依存している頂点を出す
 *   node tools/graph.mjs --violations # 「1 ファイル = 1 ノード」の違反
 *   node tools/graph.mjs --json
 */

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const SRC = join(ROOT, 'lean', 'ABC3');
const PAPERS = JSON.parse(readFileSync(join(ROOT, 'ResearchPaper', 'papers.json'), 'utf8')).papers;
/** ディレクトリ名が Lean の慣習(大文字始まり)で論文タグと違うもの。 */
const DIR_ALIAS = { PGC: 'pGC' };

const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

const KIND = 'Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example';
const ITEM_BARE = new RegExp(`^\\s*(?:${KIND})\\s+\\d+(?:\\.\\d+)+\\s*$`);
const SRC_RE = /\.src[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?item\s*:=\s*"([^"]*)"/g;
const SORRY_RE = /(^|[^\w.])sorry([^\w]|$)/;
const THEORY_RE = /def\s+theory\s*:\s*(?:ABC3\.Meta\.)?Theory/;

/** コメントと**文字列リテラル**を潰す。
 *  ★文字列を潰さないと、`.needs` の中の「… sorry 0 …」という**説明文**を
 *    本物の `sorry` と取り違える(2026-09-03 実測: 47 件中 8 件が誤検出だった)。 */
const stripComments = (s) => s
  .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
  .replace(/--[^\n]*/g, ' ')
  .replace(/"(?:[^"\\\n]|\\.)*"/g, '""');

function walk(d, acc = []) {
  for (const n of readdirSync(d)) {
    const p = join(d, n);
    if (statSync(p).isDirectory()) walk(p, acc);
    else if (p.endsWith('.lean')) acc.push(p);
  }
  return acc;
}

/** ノード。1 ファイル 1 個。 */
const nodes = new Map();          // モジュール名 → node
const byPath = new Map();         // 相対パス → node
const theories = new Set();

for (const f of walk(SRC)) {
  const rel = relative(SRC, f).replace(/\\/g, '/');
  const parts = rel.split('/');
  const bucket = parts[0];
  const dir = parts.length > 2 ? parts[1] : null;
  const mod = 'ABC3.' + rel.slice(0, -5).replace(/\//g, '.');
  const text = readFileSync(f, 'utf8');
  if (THEORY_RE.test(text) && dir) theories.add(dir);
  const items = [];
  SRC_RE.lastIndex = 0;
  let m;
  while ((m = SRC_RE.exec(text)) !== null) items.push({ paper: m[1], item: m[2] });
  const node = {
    mod, rel, bucket, dir, text,
    owner: dir ? (DIR_ALIAS[dir] ?? dir) : '(バケツ直下)',
    ownerKind: dir ? (PAPERS[DIR_ALIAS[dir] ?? dir] ? 'paper' : 'theory') : 'root',
    items,
    bare: [...new Set(items.filter((x) => ITEM_BARE.test(x.item)).map((x) => `[${x.paper}] ${x.item.trim()}`))],
    allItems: [...new Set(items.map((x) => `[${x.paper}] ${x.item.trim()}`))],
    hasSorry: stripComments(text).split('\n').some((l) => SORRY_RE.test(l)),
    imports: (text.match(/^import\s+(ABC3\.[\w.]+)/gm) ?? []).map((s) => s.replace(/^import\s+/, '')),
  };
  nodes.set(mod, node);
  byPath.set(rel, node);
}
for (const n of nodes.values()) n.imports = n.imports.filter((i) => nodes.has(i));

/* ─────────────────────────────────────────────────────────────────────────
 * ★`.needs` に書かれた**数学的依存**を辺として拾う(2026-09-05)
 *
 * ★なぜ要るか(実測)
 *   `import` 依存と数学的依存は一致しない。実例:
 *     `Skeleton/PGC/Section1Cor13.lean` の `inertia_recoverable` は
 *     `.needs` に「Proposition 1.2 を (L, H) に適用して q_L を得る段」と書いてあり、
 *     実際 `Found/PGC/InertiaTransport.lean::inertia_recoverable_of_prop12` で
 *     Prop 1.2 に完全に帰着している。しかし Prop 1.2 のある
 *     `Skeleton/PGC/Section1.lean` を **import していない**(型として要らないから)。
 *   その結果 `tools/frontier.mjs` はこのノードを **着手可能(下流 27、第 1 位)**と
 *   出していた。agent を割いても上流の sorry に当たって止まる。
 *
 * ★`Meta/Claim.lean` の `otherPaper` の docstring が予告していた誤読が、
 *   実際に起きている——同一論文内の依存を `.derivation`(指す先を持たない自由文字列)で
 *   書くと辺にならない。だからここでは **自由文の中の項目番号も拾う**。
 *
 * ★これは「疑い」であって証明ではない。文字列照合なので取りこぼしも誤検出もある。
 *   `--math-edges` で全件を印字できるようにしてあるので、疑わしければ目で見ること。
 * ───────────────────────────────────────────────────────────────────────── */

/** 項目名を「Kind N.N」だけに正規化する(`Proposition 3.4(段 A)` → `Proposition 3.4`)。 */
const bareItem = (it) => {
  const m = it.match(new RegExp(`(?:${KIND})\\s+\\d+(?:\\.\\d+)+`));
  return m ? m[0].replace(/\s+/g, ' ') : null;
};

/** `.needs … :=` の直後の `[ … ]` を括弧対応で切り出す(文字列の中の括弧は無視)。
 *  ★正規表現で `\][\s\S]*?` と書くと閉じ括弧が同じ行にある形を取りこぼす。 */
function needsBlocks(text) {
  const out = [];
  const head = /\.needs\s*:\s*List\s+(?:ABC3\.Meta\.)?ProofObligation\s*:=/g;
  let h;
  while ((h = head.exec(text)) !== null) {
    let i = head.lastIndex;
    while (i < text.length && /\s/.test(text[i])) i++;
    if (text[i] !== '[') continue;
    let depth = 0, inStr = false, j = i;
    for (; j < text.length; j++) {
      const c = text[j];
      if (inStr) { if (c === '\\') j++; else if (c === '"') inStr = false; continue; }
      if (c === '"') { inStr = true; continue; }
      if (c === '[') depth++;
      else if (c === ']') { depth--; if (depth === 0) { j++; break; } }
    }
    out.push(text.slice(i, j));
  }
  return out;
}

const NEEDS_OTHER = /\.otherPaper\s+"((?:[^"\\]|\\.)*)"\s+"((?:[^"\\]|\\.)*)"/g;
const NEEDS_FREE = /\.(derivation|implicitStep|folklore)\s+((?:\s*"(?:[^"\\]|\\.)*"\s*(?:\+\+)?)+)/g;
const ITEM_IN_TEXT = new RegExp(`(?:${KIND})\\s*\\d+(?:\\.\\d+)+`, 'g');

for (const n of nodes.values()) {
  const refs = [];
  // 自由文の項目番号は「どの論文の」かが書いていない。自分の `.src` が名乗る論文を候補にする。
  const own = [...new Set(n.items.map((x) => x.paper))];
  const cand = own.length ? own : [n.owner];
  for (const body of needsBlocks(n.text)) {
    let m;
    NEEDS_OTHER.lastIndex = 0;
    while ((m = NEEDS_OTHER.exec(body)) !== null) {
      const b = bareItem(m[2]);
      if (b) refs.push({ paper: m[1], item: b, kind: 'otherPaper' });
    }
    NEEDS_FREE.lastIndex = 0;
    while ((m = NEEDS_FREE.exec(body)) !== null) {
      const kind = m[1];
      ITEM_IN_TEXT.lastIndex = 0;
      let t;
      while ((t = ITEM_IN_TEXT.exec(m[2])) !== null) {
        const it = t[0].replace(/\s+/g, ' ');
        for (const p of cand) refs.push({ paper: p, item: it, kind });
      }
    }
  }
  n.needsRefs = refs;
}

/** `paper|項目` → その項目を `.src` で名乗るモジュール群。 */
const declaredBy = new Map();
for (const n of nodes.values()) {
  for (const it of n.items) {
    const b = bareItem(it.item);
    if (!b) continue;
    const k = `${it.paper}|${b}`;
    if (!declaredBy.has(k)) declaredBy.set(k, new Set());
    declaredBy.get(k).add(n.mod);
  }
}

/** 推移的な import の上流(自分を含まない)。 */
const upCache = new Map();
function upstreamOf(mod) {
  if (upCache.has(mod)) return upCache.get(mod);
  const seen = new Set();
  const st = [...(nodes.get(mod)?.imports ?? [])];
  while (st.length) {
    const m = st.pop();
    if (seen.has(m)) continue;
    seen.add(m);
    for (const w of nodes.get(m)?.imports ?? []) if (!seen.has(w)) st.push(w);
  }
  upCache.set(mod, seen);
  return seen;
}

/** ★import に現れない数学的依存。
 *
 *  `needsUncovered` … `.needs` が指す項目のうち、**どの提供元も import していない**もの。
 *                     これは衛生の問題であって、必ずしも塞がってはいない。
 *  `mathEdges`      … そのうち **提供元が全部 sorry を持つ**もの。
 *                     ★これだけが「数学的に塞がっている」——項目を実際に確立している
 *                     ファイルがどこにも無いという意味だから。
 *                     提供元が 1 つでも sorry 無しなら、中身は在る(import していないだけ)。
 *
 *  ★この区別が要る理由(実測): `Skeleton/GenEll/Section2.lean` の
 *    `.implicitStep` は `Proposition 1.4` を指し、その項目を `.src` で名乗るファイルは
 *    **136 本**ある(「1 ファイル = 1 ノード」の違反側)。全部を辺にすると 175 辺の雑音になる。
 *    そのうち sorry を持つものは 0 本なので、塞がってはいない。 */
for (const n of nodes.values()) {
  const up = upstreamOf(n.mod);
  const unc = [];
  const seen = new Set();
  for (const r of n.needsRefs) {
    const key = `${r.paper}|${r.item}`;
    if (seen.has(key)) continue;
    const tgt = declaredBy.get(key);
    if (!tgt) continue;                        // 木にまだ無い項目(別論文・未着手)
    const mods = [...tgt].filter((m) => m !== n.mod);
    if (!mods.length || mods.some((m) => up.has(m))) continue;   // import で覆えている
    seen.add(key);
    unc.push({ ...r, mods });
  }
  n.needsUncovered = unc;
  n.mathEdges = unc.filter((u) => u.mods.every((m) => nodes.get(m).hasSorry))
    .map((u) => ({ mods: u.mods, via: `${u.kind}: ${u.paper} ${u.item}` }));
}

/** n から import で辿れる集合(n を含む)。 */
function closure(start) {
  const seen = new Set([start.mod]);
  const st = [start.mod];
  while (st.length) {
    for (const c of nodes.get(st.pop()).imports) if (!seen.has(c)) { seen.add(c); st.push(c); }
  }
  return seen;
}

if (flag('--json')) {
  console.log(JSON.stringify({
    nodes: [...nodes.values()].map((n) => ({
      mod: n.mod, rel: n.rel, bucket: n.bucket, owner: n.owner, ownerKind: n.ownerKind,
      items: n.allItems, bare: n.bare, hasSorry: n.hasSorry, imports: n.imports,
      mathEdges: n.mathEdges,
    })),
  }, null, 1));
  process.exit(0);
}

// ── ★`.needs` が指しているのに import していない辺 ────────────────────────
if (flag('--math-edges')) {
  const all = flag('--all');
  const key = all ? 'needsUncovered' : 'mathEdges';
  const rows = [...nodes.values()].filter((n) => n[key].length);
  const cnt = rows.reduce((a, n) => a + n[key].length, 0);
  console.log('★`.needs` に書かれているが `import` に現れない依存');
  console.log(all
    ? `  ${cnt} 件 / ${rows.length} ノード（--all: 提供元が sorry 無しのものも出す＝衛生の問題）`
    : `  ${cnt} 件 / ${rows.length} ノード（★提供元が全部 sorry ＝ 数学的に塞がっている。--all で全部）`);
  console.log('  ☆文字列照合による**疑い**である。辺として採るかは目で見て決めること。\n');
  for (const n of rows.sort((a, b) => (b.hasSorry - a.hasSorry) || a.rel.localeCompare(b.rel))) {
    console.log(`  ${n.rel}${n.hasSorry ? '  ★sorry' : ''}`);
    for (const e of n[key]) {
      const ms = (e.mods ?? []).map((m) => nodes.get(m));
      const head = ms.slice(0, 3).map((t) => `${t.rel}${t.hasSorry ? ' ★sorry' : ''}`).join(', ');
      console.log(`      [${e.via ?? `${e.kind}: ${e.paper} ${e.item}`}]`);
      console.log(`      → ${head}${ms.length > 3 ? ` … 他 ${ms.length - 3} 本` : ''}`);
    }
  }
  process.exit(0);
}

// ── 個別の照会 ─────────────────────────────────────────────────────────────
const from = opt('--from');
if (from) {
  const n = byPath.get(from.replace(/\\/g, '/'));
  if (!n) { console.error(`そのノードが無い: ${from}`); process.exit(1); }
  const cl = closure(n);
  console.log(`★${n.rel}  (${n.owner} / ${n.ownerKind})`);
  console.log(`  直接の import : ${n.imports.length}`);
  console.log(`  推移閉包       : ${cl.size} ノード(自分を含む)`);
  const bad = [...cl].map((m) => nodes.get(m)).filter((x) => x.hasSorry);
  console.log(`  うち sorry を含む: ${bad.length}${bad.length ? '\n     ' + bad.map((x) => x.rel).join('\n     ') : ''}`);
  process.exit(0);
}

// ── 「1 ファイル = 1 ノード」の違反 ────────────────────────────────────────
if (flag('--violations')) {
  const v = [...nodes.values()].filter((n) => n.allItems.length >= 2)
    .sort((a, b) => b.allItems.length - a.allItems.length);
  console.log(`★「1 ファイル = 1 ノード」の違反: ${v.length} / ${nodes.size} ファイル`);
  console.log('  ☆原典の節を写したファイル(Section1.lean 等)は設計どおりである。');
  console.log('  ★作業中に足した節点が 1 ファイルに積み上がったものが、割る候補である。\n');
  for (const n of v.slice(0, 25)) {
    console.log(`  ${String(n.allItems.length).padStart(3)} 項目  ${n.rel}`);
  }
  if (v.length > 25) console.log(`  … 他 ${v.length - 25} ファイル`);
  process.exit(0);
}

// ── sorry に依存している頂点 ───────────────────────────────────────────────
if (flag('--sorry')) {
  const sorries = [...nodes.values()].filter((n) => n.hasSorry);
  console.log(`★sorry を含むノード: ${sorries.length}`);
  for (const s of sorries) console.log(`   ${s.rel}  (${s.owner})`);
  // ★`import` の到達は**上からの評価**である——ファイルを import しても、
  //   その中の `sorry` を使っているとは限らない。真偽は `#print axioms` が持つ。
  //   ここで出すのは「疑いのある範囲」であって「壊れている集合」ではない。
  console.log('\n★sorry を含むノードを(推移的に)import しているもの —— 所属別');
  console.log('  ☆import 到達は上からの評価である。実際にその sorry を使っているかは');
  console.log('    `#print axioms` でしか分からない(例: alpha_in_modl_image は誰も消費していない)。');
  const per = new Map();
  for (const n of nodes.values()) {
    const cl = closure(n);
    if (!sorries.some((s) => s.mod !== n.mod && cl.has(s.mod))) continue;
    const k = `${n.ownerKind === 'paper' ? '論文' : '理論'} ${n.owner}`;
    per.set(k, (per.get(k) ?? 0) + 1);
  }
  for (const [k, v] of [...per].sort((a, b) => b[1] - a[1])) {
    console.log(`    ${k.padEnd(22)}${String(v).padStart(6)} ノード`);
  }
  console.log('\n★`.src` を持つノードのうち、sorry を含むファイルを直接 import しているもの:');
  for (const n of [...nodes.values()].filter((x) => x.allItems.length)) {
    const hit = n.imports.filter((i) => nodes.get(i).hasSorry);
    if (hit.length) {
      console.log(`    ${n.rel}`);
      console.log(`      → ${hit.map((h) => nodes.get(h).rel).join(', ')}`);
    }
  }
  process.exit(0);
}

// ── 表示用の群団化(論文 × セクション) ─────────────────────────────────────
//
// ★内部のグラフは **1 ファイル = 1 ノード**で持つ(そのまま扱うと 1,814 節点で読めない)。
//   表示するときだけ、**論文とそのセクション**で群団にまとめる。
//   セクションは `.src` の item の番号の先頭から取る(`Lemma 3.5` → §3)。
//   ★`sectionId` は項目ごと(`genell-lemma-3-5`、172 種)なので、節の単位には使えない。
const SECT_RE = new RegExp(`^\\s*(?:${KIND})\\s+(\\d+)`);

/** ノードの群団。`.src` があれば `[論文] §N`、無ければ所属で束ねる。 */
function clusterOf(n) {
  // ★`ABC3/Found.lean` などバケツの根は「import をまとめるだけ」のファイルである。
  //   群団に混ぜると、根 → 全群団の辺が最上位に出てしまい表示が読めなくなる。
  if (!n.dir) return '(根: import をまとめるファイル)';
  const keys = n.items
    .map((x) => { const m = SECT_RE.exec(x.item); return m ? `[${x.paper}] §${m[1]}` : null; })
    .filter(Boolean);
  if (!keys.length) return `${n.ownerKind === 'theory' ? '理論' : '論文'} ${n.owner} (項目なし)`;
  // 最も多く現れた群団を主とする(ノード数の合計が総数と一致するようにするため)
  const c = new Map();
  for (const k of keys) c.set(k, (c.get(k) ?? 0) + 1);
  return [...c].sort((a, b) => b[1] - a[1])[0][0];
}

if (flag('--clusters') || flag('--display')) {
  const cl = new Map();
  for (const n of nodes.values()) {
    const k = clusterOf(n);
    const r = cl.get(k) ?? { n: 0, sorry: 0, bare: new Set(), files: [] };
    r.n++; r.sorry += n.hasSorry ? 1 : 0; n.bare.forEach((b) => r.bare.add(b)); r.files.push(n);
    cl.set(k, r);
  }
  console.log('★表示用の依存グラフ —— 論文 × セクションで群団化');
  console.log(`  内部は 1 ファイル = 1 ノード(${nodes.size} 節点)。ここではそれを ${cl.size} 群団にまとめる。\n`);
  console.log(`  ${'群団'.padEnd(26)}${'ノード'.padStart(7)}${'sorry'.padStart(7)}  条なし .src`);
  for (const [k, r] of [...cl].sort((a, b) => b[1].n - a[1].n)) {
    const b = [...r.bare].map((x) => x.replace(/^\[[^\]]*\] /, '')).sort().join(' ');
    console.log(`  ${k.padEnd(26)}${String(r.n).padStart(7)}${String(r.sorry).padStart(7)}  ${b.slice(0, 60)}`);
  }
  // 群団のあいだの辺
  const ce = new Map();
  const ROOTC = '(根: import をまとめるファイル)';
  for (const n of nodes.values()) {
    const a = clusterOf(n);
    if (a === ROOTC) continue;              // ★根からの辺は依存ではなく取りまとめである
    for (const i of n.imports) {
      const b = clusterOf(nodes.get(i));
      if (a === b || b === ROOTC) continue;
      const k = `${a} → ${b}`;
      ce.set(k, (ce.get(k) ?? 0) + 1);
    }
  }
  console.log(`\n  群団のあいだの辺: ${ce.size} 種 / ${[...ce.values()].reduce((s, v) => s + v, 0)} 本(上位 15)`);
  for (const [k, v] of [...ce].sort((a, b) => b[1] - a[1]).slice(0, 15)) {
    console.log(`    ${String(v).padStart(5)}  ${k}`);
  }
  process.exit(0);
}

// ── 要約 ───────────────────────────────────────────────────────────────────
const owner = opt('--owner');
const sel = [...nodes.values()].filter((n) => !owner || n.owner === owner);
const edges = sel.reduce((s, n) => s + n.imports.length, 0);

console.log('★依存グラフ(ノード = Lean ファイル)');
console.log(`  ノード ${sel.length} / ファイル ${sel.length}  ← ★定義から一致する`);
console.log(`  辺(import) ${edges} 本  ★Lean が保証するので循環は無い`);

const byOwner = new Map();
for (const n of sel) {
  const k = `${n.ownerKind === 'paper' ? '論文' : n.ownerKind === 'theory' ? '理論' : '—'} ${n.owner}`;
  const r = byOwner.get(k) ?? { n: 0, bare: 0, items: 0, sorry: 0 };
  r.n++; r.bare += n.bare.length; r.items += n.allItems.length; r.sorry += n.hasSorry ? 1 : 0;
  byOwner.set(k, r);
}
console.log('\n  所属別(★③ どの論文・理論に属するか):');
console.log(`    ${'所属'.padEnd(22)}${'ノード'.padStart(7)}${'条なし'.padStart(7)}${'項目'.padStart(7)}${'sorry'.padStart(7)}`);
for (const [k, r] of [...byOwner].sort((a, b) => b[1].n - a[1].n)) {
  console.log(`    ${k.padEnd(22)}${String(r.n).padStart(7)}${String(r.bare).padStart(7)}${String(r.items).padStart(7)}${String(r.sorry).padStart(7)}`);
}

const vio = sel.filter((n) => n.allItems.length >= 2).length;
console.log(`\n  ★「1 ファイル = 1 ノード」: 守れているのは ${sel.length - vio} / ${sel.length}` +
            ` (${(100 * (sel.length - vio) / sel.length).toFixed(0)}%)  —— 内訳は --violations`);
console.log(`  ★sorry を含むノード: ${sel.filter((n) => n.hasSorry).length}  —— 影響は --sorry`);
console.log(`  ☆登記済みの理論: ${[...theories].join(', ')}`);
