#!/usr/bin/env node
/**
 * 持ち場（brief）—— 1 ノードぶんの文脈を、sub-agent に渡せる形に切り出す。
 *
 * ★この道具の位置づけ（`ResearchPaper/orchestration.md` §3）
 *   数学全体を 1 つのコンテキストに入れることはできない。かわりに
 *   **依存グラフを外部の真実**とし、agent には「そのノードだけ」を渡す。
 *   ここが「グラフ」と「agent」の唯一の受け渡し口である。
 *
 * ★何を渡すか（渡さないか）
 *   渡す:   sorry の宣言そのもの / その docstring（原文の逐語と監査の経緯が入っている）
 *           / `.src`（原典の位置）/ `.needs`（証明が要求するもの）
 *           / 直接の依存と被依存 / 上流に残る sorry / 下流への影響
 *   渡さない: 木の全体像。**agent は自分の持ち場の外を知らなくてよい**。
 *
 * ★2026-09-07（メタ backlog M34）**まだ `.lean` が無いノードの持ち場**を足した。
 *   それまでの入口は `--node` / `--mod` の 2 つだけで、**どちらも既存ファイルを要求する**。
 *   ところが 2026-09-06 に Lean が着地した 15 session のうち **13 が新規ファイル**で、
 *   その持ち場は人が手で書いていた。★**報告された brief の誤り 5 件は 5 件とも
 *   「まだファイルが無いノード」のもの**で、共通の原因は
 *   **本体が原文を自分で読まず、agent の報告から段取りを書いたこと**である。
 *   ⇒ `--paper <鍵> --item "<項目名>"` で、**原典の逐語と在庫を機械で差し込む**。
 *   ★この入口は**段取りを書かない**。原文に無いことを brief に足すと、
 *     いまの手書きの誤りを機械化するだけになる。並べるところまでで止める。
 *
 * 使い方:
 *   node tools/brief.mjs --node Skeleton/PGC/Section1.lean
 *   node tools/brief.mjs --mod ABC3.Skeleton.PGC.Section1
 *   node tools/brief.mjs --node ... --json
 *   node tools/brief.mjs --paper Yoshida08 --item "Proposition 6.6"   ← まだ無いノード
 *   node tools/brief.mjs --paper Yoshida08 --id prop-6-6              ← id で直に引く
 *   node tools/brief.mjs --paper Yoshida08 --list                     ← 項目の一覧
 */

import { execFileSync } from 'node:child_process';
import { readFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));

const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

const relArg = opt('--node');
const modArg = opt('--mod');
const paperArg = opt('--paper');

// ★新しい入口は既存の 2 つと排他。ここで分岐して抜ける（以下の既存経路は無傷）。
if (paperArg) { paperBrief(paperArg, opt('--item'), opt('--id'), flag('--list'), flag('--json')); process.exit(0); }

if (!relArg && !modArg) {
  console.error('使い方: node tools/brief.mjs --node <Skeleton/…/X.lean> | --mod <ABC3.…>');
  console.error('        node tools/brief.mjs --paper <鍵> --item "<項目名>"   … まだ .lean が無いノード');
  console.error('        node tools/brief.mjs --paper <鍵> --list              … その論文の項目一覧');
  process.exit(2);
}

const graph = JSON.parse(execFileSync('node', [join(ROOT, 'tools', 'graph.mjs'), '--json'],
  { encoding: 'utf8', maxBuffer: 1 << 28 }));
const nodes = graph.nodes;
const byMod = new Map(nodes.map((n) => [n.mod, n]));

const norm = (s) => s.replace(/\\/g, '/');
const node = modArg
  ? byMod.get(modArg)
  : nodes.find((n) => norm(n.rel) === norm(relArg));
if (!node) { console.error(`ノードが見つからない: ${relArg ?? modArg}`); process.exit(1); }

/* ---------- 依存・被依存 ---------- */
/* ★`import` 依存と数学的依存は一致しない。`graph.mjs` が `.needs` から拾った
 *   `mathEdges` を import と同じ扱いにする（`tools/frontier.mjs` と同じ規則）。
 *   ★これが無いと、数学的に塞がっているノードの持ち場を agent に渡してしまう。 */
const mathUp = new Map();
const mathVia = new Map();
for (const n of nodes) {
  const ms = [...new Set((n.mathEdges ?? []).flatMap((e) => e.mods))]
    .filter((m) => byMod.has(m) && m !== n.mod);
  if (ms.length) mathUp.set(n.mod, ms);
  for (const e of n.mathEdges ?? []) for (const m of e.mods) mathVia.set(`${n.mod}→${m}`, e.via);
}

const rdeps = new Map();
const addR = (from, to) => {
  if (!rdeps.has(from)) rdeps.set(from, []);
  if (!rdeps.get(from).includes(to)) rdeps.get(from).push(to);
};
for (const n of nodes) {
  for (const im of n.imports) { if (byMod.has(im)) addR(im, n.mod); }
  for (const m of mathUp.get(n.mod) ?? []) addR(m, n.mod);
}
const reach = (start, edgesOf) => {
  const seen = new Set(); const st = [...(edgesOf(start) ?? [])];
  while (st.length) { const m = st.pop(); if (seen.has(m)) continue; seen.add(m);
    for (const w of edgesOf(m) ?? []) if (!seen.has(w)) st.push(w); }
  return seen;
};
const downOf = (m) => rdeps.get(m) ?? [];
const importUpOf = (m) => (byMod.get(m)?.imports ?? []).filter((x) => byMod.has(x));
const upOf = (m) => [...importUpOf(m), ...(mathUp.get(m) ?? [])];
const sorryMods = new Set(nodes.filter((n) => n.hasSorry).map((n) => n.mod));

const down = reach(node.mod, downOf);
const up = reach(node.mod, upOf);
const upImport = reach(node.mod, importUpOf);
const blockers = [...up].filter((m) => sorryMods.has(m)).sort();
let dsItems = 0; for (const m of down) dsItems += (byMod.get(m)?.items.length ?? 0);

/* ---------- Lean ファイルの中身 ---------- */
const path = join(ROOT, 'lean', 'ABC3', node.rel);
const text = readFileSync(path, 'utf8');

/** コメントと文字列を空白で潰す（sorry の誤検出を避ける。graph.mjs と同じ手筋）。 */
const strip = (s) => s
  .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
  .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length))
  .replace(/"(?:[^"\\\n]|\\.)*"/g, (m) => ' '.repeat(m.length));
const bare = strip(text);

const DECL = /^(?:@\[[^\]]*\]\s*)?(?:noncomputable\s+|private\s+|protected\s+|scoped\s+)*(theorem|lemma|def|instance|abbrev|structure)\s+([A-Za-z_À-￿][\w.'!?À-￿]*)/gm;
const decls = [];
for (let m; (m = DECL.exec(bare));) decls.push({ kind: m[1], name: m[2], start: m.index });
for (let i = 0; i < decls.length; i++) {
  decls[i].end = i + 1 < decls.length ? decls[i + 1].start : text.length;
}

const SORRY = /(^|[^\w.])sorry([^\w]|$)/;
const withSorry = decls.filter((d) => SORRY.test(bare.slice(d.start, d.end)));

/** 直前の `/-- … -/` を拾う（原文の逐語・監査の経緯はここにある）。 */
function docstringBefore(start) {
  const head = text.slice(0, start);
  const close = head.lastIndexOf('-/');
  if (close < 0) return null;
  const openA = head.lastIndexOf('/--', close);
  if (openA < 0) return null;
  const between = head.slice(close + 2).trim();
  // 宣言との間に他の宣言やコードが挟まっていたら、その docstring ではない
  if (between && !/^(?:@\[[^\]]*\]\s*)*$/.test(between)) return null;
  return text.slice(openA + 3, close).trim();
}

/** `NAME.src` / `NAME.needs` の本体をそのまま取る。 */
function bodyOf(name) {
  const d = decls.find((x) => x.name === name);
  return d ? text.slice(d.start, d.end).trimEnd() : null;
}

const goals = withSorry.map((d) => ({
  kind: d.kind,
  name: d.name,
  statement: text.slice(d.start, d.end).trimEnd(),
  doc: docstringBefore(d.start),
  src: bodyOf(`${d.name}.src`),
  needs: bodyOf(`${d.name}.needs`),
}));

const result = {
  node: node.rel,
  mod: node.mod,
  owner: node.owner,
  ownerKind: node.ownerKind,
  bucket: node.bucket,
  items: node.items,
  status: blockers.length ? 'blocked' : (node.hasSorry ? 'startable' : 'done'),
  impact: { downstream: down.size, downstreamItems: dsItems },
  blockers: blockers.map((m) => ({
    rel: byMod.get(m)?.rel ?? m,
    // ★import には現れず、`.needs` から拾って初めて見えた依存
    viaNeeds: upImport.has(m) ? null : (mathVia.get(`${node.mod}→${m}`) ?? '（推移的）'),
  })),
  imports: node.imports.filter((m) => byMod.has(m)).map((m) => ({
    rel: byMod.get(m).rel, hasSorry: byMod.get(m).hasSorry,
  })),
  usedByDirect: downOf(node.mod).map((m) => byMod.get(m).rel),
  goals,
  file: `lean/ABC3/${node.rel}`,
};

if (flag('--json')) { console.log(JSON.stringify(result, null, 1)); process.exit(0); }

/* ---------- 人と agent が読む形 ---------- */
const L = [];
L.push(`# 持ち場: ${result.node}`);
L.push('');
L.push(`所属: ${result.owner}（${result.ownerKind === 'paper' ? '論文' : '理論'}） / バケツ: ${result.bucket}`);
L.push(`状態: **${result.status}**`);
L.push(`影響: これが解けると下流 **${result.impact.downstream} ノード** / **${result.impact.downstreamItems} 項目** が sorry から解放される`);
L.push(`ファイル: \`${result.file}\``);
L.push('');

if (result.blockers.length) {
  L.push('## ★上流に残る sorry（これらが先）');
  for (const b of result.blockers) {
    L.push(`- \`${b.rel}\`${b.viaNeeds ? `  ★**import には現れない依存**（\`.needs\` の ${b.viaNeeds}）` : ''}`);
  }
  L.push('');
  L.push('☆このノードには**まだ着手できない**。上流を先に片付けること。');
  if (result.blockers.some((b) => b.viaNeeds)) {
    L.push('★import していないので Lean は何も言わないが、**数学的には塞がっている**。');
    L.push('  この形は 2026-09-05 に `tools/graph.mjs --math-edges` で機械的に拾えるようにした。');
  }
  L.push('');
}

L.push('## 目標（このノードの sorry）');
if (!goals.length) L.push('（sorry 無し）');
for (const g of goals) {
  L.push('');
  L.push(`### \`${g.name}\``);
  if (g.doc) {
    L.push('');
    L.push('<details><summary>docstring（原文の逐語と監査の経緯）</summary>');
    L.push('');
    L.push('```');
    L.push(g.doc);
    L.push('```');
    L.push('');
    L.push('</details>');
  }
  L.push('');
  L.push('```lean');
  L.push(g.statement);
  L.push('```');
  if (g.src) { L.push(''); L.push('```lean'); L.push(g.src); L.push('```'); }
  if (g.needs) {
    L.push('');
    L.push('**証明が要求するもの（`.needs`、下界）**');
    L.push('```lean'); L.push(g.needs); L.push('```');
  }
}
L.push('');
L.push('## 依存（このノードが import しているもの）');
for (const im of result.imports) L.push(`- \`${im.rel}\`${im.hasSorry ? '  ★sorry あり' : ''}`);
if (!result.imports.length) L.push('（木の外＝Mathlib のみ）');
L.push('');
L.push('## 被依存（このノードを直接 import しているもの）');
for (const u of result.usedByDirect) L.push(`- \`${u}\``);
if (!result.usedByDirect.length) L.push('（無し＝葉）');
L.push('');
L.push('## 許された手順（この順で）');
L.push('1. **在庫**を引く（既にあるものを書き直さない）:');
L.push('   `node tools/decl-index.mjs` → `.cache/decl-index.txt` を結論のリテラルで grep');
L.push('   Mathlib は `.cache/mathlib-index.txt`。★`Unknown constant` は「無い」ではなく');
L.push('   「import していない」ことが多い（`tools/lean-idioms.md` #68）');
L.push('   ★**「mathlib に無い」と書く前に** `node tools/absent-recheck.mjs --try \'<正規表現>\'`');
L.push('   で件数を測り、`.absent` に ``re:`<正規表現>`→<件数>`` と測定日を書く（G11）。');
L.push('   索引には**無名 instance も入っている**（`instance : CompactSpace Gal(K/k)` のような');
L.push('   名前を持たないものが mathlib の instance の 58.8%。名前で grep しても出ない）');
L.push('2. **証明は MCP `lean_check` で通す**（0.01〜1 秒）。ファイルに書くのはその後。');
L.push('3. 書いたら `lake build <対象モジュール>` のみ。全体ビルドは最後に 1 回。');
L.push('4. 配管で詰まったら `tools/lean-idioms.md` を引く。新しい失敗形なら 1 行足す。');
L.push('5. 逸脱（原典の読み替え・前提の追加）をしたら docstring に必ず記録する。');
console.log(L.join('\n'));

/* ══════════════════════════════════════════════════════════════════════════
 *  ここから下: **まだ `.lean` が無いノード**の持ち場（`--paper` / `--item`）
 *
 *  ★出すもの（この順。優先度そのもの）
 *    1. 原文の逐語（`1_Structured` の `.verbatim`）と、あれば `.reading`
 *       ——★**最重要**。本体が原文を読まずに段取りを書くのを防ぐのが目的。
 *    2. `.src` の雛形。`pdfPage` / `item` / `sectionId` は **HTML の属性から取る**
 *       ——手で書くとずれる（G1 が項目名の一致を見るようになった。実測 3 件）。
 *    3. 在庫。**同じ節で既に立っているノードの宣言を「型」ごと**並べる
 *       ——`span_smul_uniformizer` は Cor 6.3 用に作られたが型は Lemma 6.5 の入力だった。
 *       名前（＝作られた目的）で引くと落ちるので、隣接項目の**型を全部見せる**。
 *    4. 依存（本文が名指しする項目）と被依存（この項目を名指しする項目）。
 *    5. 省略の合図（`hedge-index.mjs --paper … --item …` をそのまま埋める）。
 *
 *  ★出さないもの: **段取り**。原文に無いことをここで書くと、いまの手書きの誤りを
 *    機械化するだけになる。証明の方針は原文（1）を読んだ人が決める。
 * ══════════════════════════════════════════════════════════════════════════ */

/* ★以下は**関数宣言だけ**にしてある（`const` を使わない）。この節の入口は
 *   ファイル先頭で呼ばれるので、`const` で書くと既存の top-level コードより
 *   後ろに置いた宣言が TDZ に落ちる。関数宣言は巻き上がるので安全である。 */

/** 構造化 HTML の置き場と、論文の登記簿。 */
function structDir() { return join(ROOT, 'ResearchPaper', '1_Structured'); }
function papersJson() { return join(ROOT, 'ResearchPaper', 'papers.json'); }

/** 装飾クラス。★`check.mjs` の `DECOR_CLASSES` と**同じ語彙**でなければならない
 *  （`1_Structured/README.md` §3。1 つの体系・2 つの表現）。 */
function decorClasses() { return 'ul1|ul2|bar|hat|tilde|dot1|dot2|bb|scr|frak|prime'; }

function entTable() {
  return {
  amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ',
  ldquo: '“', rdquo: '”', lsquo: '‘', rsquo: '’',
  minus: '−', times: '×', ge: '≥', le: '≤', ne: '≠', equiv: '≡', cong: '≅',
  rarr: '→', larr: '←', harr: '↔', rArr: '⇒', lArr: '⇐', hArr: '⇔', mapsto: '↦',
  sube: '⊆', supe: '⊇', sub: '⊂', sup: '⊃', isin: '∈', notin: '∉', ni: '∋',
  cap: '∩', cup: '∪', empty: '∅', sum: '∑', prod: '∏', part: '∂', infin: '∞',
  or: '∨', and: '∧', not: '¬', forall: '∀', exist: '∃', radic: '√', prop: '∝',
  sect: '§', hellip: '…', mdash: '—', ndash: '–', middot: '·', bull: '•',
  prime: '′', Prime: '″', deg: '°', plusmn: '±', divide: '÷', sdot: '⋅',
  alpha: 'α', beta: 'β', gamma: 'γ', delta: 'δ', epsilon: 'ε', zeta: 'ζ',
  eta: 'η', theta: 'θ', iota: 'ι', kappa: 'κ', lambda: 'λ', mu: 'μ', nu: 'ν',
  xi: 'ξ', pi: 'π', rho: 'ρ', sigma: 'σ', tau: 'τ', upsilon: 'υ', phi: 'φ',
  chi: 'χ', psi: 'ψ', omega: 'ω',
  Gamma: 'Γ', Delta: 'Δ', Theta: 'Θ', Lambda: 'Λ', Xi: 'Ξ', Pi: 'Π',
  Sigma: 'Σ', Phi: 'Φ', Psi: 'Ψ', Omega: 'Ω',
  };
}
function decodeEnt(s) {
  const ent = entTable();
  return s
    .replace(/&#x([0-9a-fA-F]+);/g, (_, h) => String.fromCodePoint(parseInt(h, 16)))
    .replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(parseInt(d, 10)))
    .replace(/&([a-zA-Z]+);/g, (m, n) => (n in ent ? ent[n] : m));
}

/**
 * HTML の断片を平文にする。
 *
 * ★`applyTxt` が 2 つの用途を分ける。ここが**手で書くと必ずずれる**ところである:
 *   - `false`（**読む用**）: 要素の中身をそのまま使う。原文どおり ⟨σ⟩ が出る。
 *   - `true` （**貼る用**）: `data-txt` の値に置き換える。すなわち
 *     **`pdftotext` に見えている形**。`check.mjs` の引用照合（`原文 (鍵 p.N):` の
 *     `>` 行）は Lean の逐語を **PDF テキストと**照合するので、
 *     `data-txt=""`（＝抽出に出ない記号。生成部分群の ⟨ ⟩ 等）を残したまま
 *     docstring に貼ると**落ちる**。貼る用は必ず `applyTxt = true`。
 *
 * 装飾は `X[bb]` の角括弧記法で出す。これは `check.mjs` の
 * `leanQuoteProjection` が落とす記法そのもの（下付き `_`・上付き `^` も落ちる）
 * なので、貼る用の出力は**そのまま引用照合を通る**。
 *
 * ★`braces`（読む用だけ true）は 2 文字以上の添字を `{}` で括る。
 *   括らないと `i<sub>j−1</sub>` が `i_j−1` になり、**`i_{j−1}` と `i_j − 1` が
 *   区別できない**——Proposition 6.6 (i) は前者で、後者だと自明な主張になる。
 *   貼る用では括らない。`leanQuoteProjection` は `_` `^` は落とすが
 *   `{` `}` は落とさないので、括ると引用照合が通らなくなる。
 */
function renderText(html, applyTxt, braces = false) {
  let s = html;
  const reTxt = new RegExp('<(\\w+)\\b[^>]*\\bdata-txt="([^"]*)"[^>]*>([^<]*)</\\1>');
  const reDecor = new RegExp(`<span class="(${decorClasses()})">([^<]*)</span>`);
  const script = (mark, x) => (braces && x.length > 1 ? `${mark}{${x}}` : `${mark}${x}`);
  for (let guard = 0; guard < 2000; guard++) {
    const before = s;
    s = s.replace(reTxt, (_m, _t, txt, inner) => (applyTxt ? txt : inner));
    s = s.replace(reDecor, (_m, cls, inner) => `${inner}[${cls}]`);
    s = s.replace(/<sub>([^<]*)<\/sub>/, (_m, x) => script('_', x));
    s = s.replace(/<sup>([^<]*)<\/sup>/, (_m, x) => script('^', x));
    s = s.replace(/<code>([^<]*)<\/code>/, (_m, x) => `\`${x}\``);
    if (s === before) break;
  }
  return decodeEnt(s.replace(/<br\s*\/?>/g, '\n').replace(/<[^>]+>/g, ''))
    .replace(/[ \t]+/g, ' ')
    .replace(/\n\s*/g, '\n')
    .trim();
}

/** section を属性・verbatim・reading・hdr に割る（`check.mjs` の `extractSections` の拡張）。 */
function sectionsOf(html) {
  const out = [];
  const re = /<section\b([^>]*)>([\s\S]*?)<\/section>/g;
  for (let m; (m = re.exec(html));) {
    const attrs = {};
    for (let a, ar = /([\w-]+)="([^"]*)"/g; (a = ar.exec(m[1]));) attrs[a[1]] = a[2];
    const g = (re2) => { const x = re2.exec(m[2]); return x ? x[1] : null; };
    out.push({
      attrs,
      verbatim: g(/<div class="verbatim">([\s\S]*?)<\/div>/),
      reading: g(/<div class="reading">([\s\S]*?)<\/div>/),
      hdr: g(/<p class="hdr">([\s\S]*?)<\/p>/),
    });
  }
  return out;
}

function walkFiles(dir, acc = []) {
  if (!existsSync(dir)) return acc;
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walkFiles(p, acc); else acc.push(p);
  }
  return acc;
}

/** 項目名から「種別 + 番号」だけを取る。★`check.mjs` の `itemKeyOf` と同じ規則
 *  （G1 が同じ鍵で一致を見るので、ここがずれると brief と門で答えが変わる）。 */
function itemKeyOf(s) {
  const m = /(Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example|Claim|Assertion|Step|Section|Chapter)\s*([0-9]+(?:[.-][0-9]+)*)/
    .exec(s ?? '');
  return m ? `${m[1]} ${m[2].replace(/-/g, '.')}` : null;
}

/** 論文タグ → その論文の section（構造化されている分だけ）。 */
function collectPaper(paperKey) {
  const secs = [];
  const files = walkFiles(structDir())
    .filter((p) => p.endsWith('.html'))
    .filter((p) => !p.includes('.legacy.'));
  for (const f of files) {
    const html = readFileSync(f, 'utf8');
    if (!html.includes(`data-paper="${paperKey}"`)) continue;
    const h1 = /<h1[^>]*>([\s\S]*?)<\/h1>/.exec(html);
    for (const s of sectionsOf(html)) {
      if (!/\bstatement\b/.test(s.attrs['class'] ?? '')) continue;
      if (s.attrs['data-paper'] !== paperKey) continue;
      secs.push({ ...s, file: relative(ROOT, f).replace(/\\/g, '/'),
        sectionTitle: h1 ? renderText(h1[1], false) : null });
    }
  }
  return secs;
}

/** 構造化されている論文の鍵（`--paper` に何が渡せるか）。 */
function structuredPaperKeys() {
  const keys = new Set();
  for (const f of walkFiles(structDir()).filter((p) => p.endsWith('.html'))) {
    for (let m, re = /data-paper="([^"]+)"/g; (m = re.exec(readFileSync(f, 'utf8')));) keys.add(m[1]);
  }
  return [...keys].sort();
}

/** 木の中で「その項目」を名乗っているノードを引く（`graph.mjs` の `items` が正典）。
 *
 *  ★鍵は **`論文タグ|種別+番号`**。番号だけで引いてはならない——
 *  `hedge-index.mjs` の実測どおり、**見出し 129 件のうち 25 件（19%）を 2 本以上の
 *  論文が名乗る**（`Theorem 4.2` は FrdI / CorrHyp / Falt1 / pGC の 4 本）。
 *  実際、番号だけで引いた最初の版は Yoshida08 の `Lemma 6.5` に **FrdI の
 *  `Lemma65.jean`（六指数定理）**を在庫として出した。この取り違えは
 *  M34 が減らそうとしている誤りそのものなので、鍵に論文を入れる。
 *  `graph.mjs` の items は 4549 件すべてが `[タグ] 項目` の形（実測）。 */
function itemsInTree() {
  const g = JSON.parse(execFileSync('node', [join(ROOT, 'tools', 'graph.mjs'), '--json'],
    { encoding: 'utf8', maxBuffer: 1 << 28 }));
  const map = new Map(); // `タグ|種別+番号` → [{rel, hasSorry, raw}]
  for (const n of g.nodes) {
    for (const it of n.items ?? []) {
      const tag = /^\[([^\]]+)\]/.exec(it);
      const k = itemKeyOf(it);
      if (!tag || !k) continue;
      const key = `${tag[1]}|${k}`;
      if (!map.has(key)) map.set(key, []);
      map.get(key).push({ rel: n.rel.replace(/\\/g, '/'), hasSorry: n.hasSorry, raw: it });
    }
  }
  return map;
}

/** 宣言の**結論**だけを取る（束縛子を捨てる）。
 *
 *  ★`.cache/decl-index.txt` の statement は完全だが、先頭の 150 字はほぼ
 *  `{A B : Type*} [CommRing A] …` の束縛子で、**結論が切り落とされる**。
 *  在庫は「作られた目的（名前）」ではなく**型**で引くもので、効くのは結論の形である
 *  （`span_smul_uniformizer` は Cor 6.3 用に作られたが、結論
 *  `maximalIdeal B = Ideal.span {τ • α}` が Lemma 6.5 の入力だった）。 */
function conclusionOf(stmt) {
  let s = (stmt ?? '')
    .replace(/^@\[[^\]]*\]\s*/, '')
    .replace(/^(?:noncomputable|private|protected|scoped|nonrec)\s+/g, '')
    .replace(/^(?:theorem|lemma|def|instance|abbrev|structure)\s+\S*\s*/, '');
  const OPEN = { '(': ')', '[': ']', '{': '}', '⦃': '⦄' };
  for (;;) {
    s = s.replace(/^\s+/, '');
    const close = OPEN[s[0]];
    if (!close) break;
    let depth = 0; let i = 0;
    for (; i < s.length; i++) {
      if (OPEN[s[i]]) depth++;
      else if (s[i] === close || s[i] === ')' || s[i] === ']' || s[i] === '}' || s[i] === '⦄') {
        depth--; if (depth === 0) break;
      }
    }
    if (i >= s.length) return (stmt ?? '').trim(); // 対応が取れない → 丸ごと返す
    s = s.slice(i + 1);
  }
  s = s.replace(/^\s*/, '');
  return s.startsWith(':') ? s.slice(1).trim() : (stmt ?? '').trim();
}

/** ノードの `.lean` から `NAME.src` の中身を拾う（どの宣言がその項目を担いでいるか）。 */
function srcDeclsOf(rel) {
  const p = join(ROOT, 'lean', 'ABC3', rel);
  if (!existsSync(p)) return [];
  const t = readFileSync(p, 'utf8');
  const out = [];
  const re = /def\s+([\w.'!?À-￿]+)\.src\s*:\s*[\w.]*Source\s*:=\s*\{([^}]*)\}/g;
  for (let m; (m = re.exec(t));) {
    const f = (k) => { const x = new RegExp(`${k}\\s*:=\\s*"?([^",}]*)"?`).exec(m[2]); return x ? x[1].trim() : null; };
    out.push({ decl: m[1], item: f('item'), sectionId: f('sectionId'), pdfPage: f('pdfPage') });
  }
  return out;
}

/** ノードの `.lean` から theorem/lemma の**署名**を取り、結論だけにする。
 *
 *  ★`.cache/decl-index.txt` を使わない理由は 2 つ、どちらも実測:
 *    1) 索引は statement を **約 400 字で切る**。束縛子が長い宣言では
 *       `IsPGroup p (` のように**結論が落ちる**——在庫として役に立たない。
 *    2) 索引は作った時刻で固まる。当日できたノード（`ConjugateSumValuation.lean`）は
 *       **索引に無い**。持ち場は「いま木にあるもの」を見せなければならない。
 *  読むのは隣接ノード数本だけなので、木を全部読むわけではない。 */
function declConclusionsOf(rel) {
  const p = join(ROOT, 'lean', 'ABC3', rel);
  if (!existsSync(p)) return null;
  const text = readFileSync(p, 'utf8');
  const bare = text
    .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
    .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length))
    .replace(/"(?:[^"\\\n]|\\.)*"/g, (m) => ' '.repeat(m.length));
  const re = /^(?:@\[[^\]]*\]\s*)?(?:noncomputable\s+|private\s+|protected\s+|scoped\s+|nonrec\s+)*(theorem|lemma)\s+([A-Za-z_À-￿][\w.'!?À-￿]*)/gm;
  const out = [];
  for (let m; (m = re.exec(bare));) {
    if (/\.(src|needs|nonvacuous|waiting|record|loadBearing|negControl)$/.test(m[2])) continue;
    // 署名の終わり = 深さ 0 の最初の `:=`
    let depth = 0; let end = -1;
    for (let i = m.index; i < bare.length; i++) {
      const c = bare[i];
      if ('([{⦃⟨'.includes(c)) depth++;
      else if (')]}⦄⟩'.includes(c)) depth--;
      else if (depth === 0 && c === ':' && bare[i + 1] === '=') { end = i; break; }
      else if (depth === 0 && c === '\n' && /^\s*(?:@\[|theorem |lemma |def |instance |\/--)/.test(bare.slice(i + 1, i + 40))) break;
    }
    if (end < 0) continue;
    const sig = text.slice(m.index, end).replace(/\s+/g, ' ').trim();
    out.push({ name: m[2], line: text.slice(0, m.index).split('\n').length, concl: conclusionOf(sig) });
  }
  return out;
}

/** `.cache/decl-index.txt`（無ければ null）。行は `kind\t完全名\tfile:line\tstatement`。 */
function declIndex() {
  const p = join(ROOT, '.cache', 'decl-index.txt');
  if (!existsSync(p)) return null;
  const rows = readFileSync(p, 'utf8').split('\n').filter(Boolean).map((l) => {
    const [kind, full, loc, stmt] = l.split('\t');
    return { kind: (kind ?? '').trim(), full, loc: loc ?? '', stmt: stmt ?? '' };
  });
  return { rows, mtime: statSync(p).mtime };
}

/** 索引を引く語から落とす散文語。★数学の語（ramification / uniformizer …）は残す。 */
function stopWords() {
  return new Set(`about above after again against all also and any are because been before
being below between both but can cannot did does doing down during each few for from further
had has have having here hence how into its itself just let more most must not now off once
only other ought our out over own same shall should since some such than that the their them
then there these they this those through thus too under until upon very was were what when
where which while who whom why will with would you your set get put obtain follows following
have write given gives take taken thus note assume suppose respectively above below case cases
first second third form forms map maps such that with which where there exists exist unique
proof lemma theorem proposition corollary definition remark section chapter equation`
    .split(/\s+/).filter(Boolean));
}

function paperBrief(paperKey, itemArg, idArg, listOnly, asJson) {
  /* ── 論文の同定。★構造化されていなければ**素直に失敗する**（劣化出力を出さない）。 */
  const reg = JSON.parse(readFileSync(papersJson(), 'utf8')).papers;
  const meta = reg[paperKey];
  if (!meta) {
    console.error(`papers.json に鍵 "${paperKey}" が無い。`);
    console.error(`使える鍵（papers.json、${Object.keys(reg).length} 本）: ${Object.keys(reg).sort().join(' ')}`);
    process.exit(2);
  }
  const secs = collectPaper(paperKey);
  if (!secs.length) {
    console.error(`★"${paperKey}"（${meta.title}）は **1_Structured/ に無い**。`);
    console.error('  この入口は構造化された逐語の上でしか動かない（原文を機械で差し込むのが目的なので、');
    console.error('  構造化前に劣化した brief を出すことはしない）。先に 1_Structured を作ること。');
    console.error(`  いま構造化されている鍵: ${structuredPaperKeys().join(' ')}`);
    process.exit(3);
  }
  const files = [...new Set(secs.map((s) => s.file))].sort();

  if (listOnly || (!itemArg && !idArg)) {
    console.log(`# [${paperKey}] ${meta.title} —— 1_Structured の項目 ${secs.length} 件`);
    console.log(`ファイル: ${files.join(' / ')}`);
    console.log('');
    for (const s of secs) {
      console.log(`  ${(s.attrs['id'] ?? '').padEnd(14)} p.${String(s.attrs['data-pdf-page'] ?? '?').padStart(3)}  ${s.attrs['data-item'] ?? ''}`);
    }
    if (!listOnly) {
      console.log('');
      console.log('★--item "<項目名>" か --id <id> で 1 件の持ち場を出す。');
      process.exit(2);
    }
    return;
  }

  /* ── 項目の同定。id が最優先、次に完全一致、次に「種別+番号」の鍵。 */
  let hits;
  let how;
  if (idArg) { hits = secs.filter((s) => s.attrs['id'] === idArg); how = `id="${idArg}"`; }
  else {
    const want = itemArg.trim().toLowerCase();
    hits = secs.filter((s) => (s.attrs['data-item'] ?? '').trim().toLowerCase() === want);
    how = '完全一致';
    if (!hits.length) {
      const k = itemKeyOf(itemArg);
      if (k) { hits = secs.filter((s) => itemKeyOf(s.attrs['data-item']) === k); how = `種別+番号 "${k}"（check.mjs G1 と同じ鍵）`; }
    }
    if (!hits.length) { hits = secs.filter((s) => (s.attrs['data-item'] ?? '').toLowerCase().includes(want)); how = '部分一致'; }
  }
  if (!hits.length) {
    console.error(`★[${paperKey}] に "${idArg ?? itemArg}" は無い。`);
    console.error('  ある項目（--list と同じ）:');
    for (const s of secs) console.error(`    ${(s.attrs['id'] ?? '').padEnd(14)} ${s.attrs['data-item'] ?? ''}`);
    process.exit(4);
  }
  if (hits.length > 1) {
    console.error(`★"${idArg ?? itemArg}" が ${hits.length} 件に当たる。--id で 1 つに絞ること:`);
    for (const s of hits) console.error(`    --id ${s.attrs['id']}   ${s.attrs['data-item']}`);
    process.exit(4);
  }
  const sec = hits[0];
  const id = sec.attrs['id'];
  const item = sec.attrs['data-item'] ?? '';
  const page = sec.attrs['data-pdf-page'] ?? null;
  const key = itemKeyOf(item);

  /* ── 直前の setup（記号の定義はそこにある）。同じファイルの中だけ見る。 */
  const inFile = secs.filter((s) => s.file === sec.file);
  const at = inFile.findIndex((s) => s.attrs['id'] === id);
  let setup = null;
  for (let i = at - 1; i >= 0; i--) if (/\bsetup\b/.test(inFile[i].attrs['class'] ?? '')) { setup = inFile[i]; break; }

  /* ── 本文が名指しする項目（＝依存の候補）と、この項目を名指しする項目（＝被依存）。
   *    ★ここは**原文の literal**しか見ない。推測はしない。 */
  const REF = /(Theorem|Proposition|Corollary|Definition|Lemma|Remark|Claim)\s*([0-9]+(?:\.[0-9]+)+)/g;
  const bodyText = `${renderText(sec.verbatim ?? '', false, true)}\n${renderText(sec.reading ?? '', false)}`;
  const refs = [...new Set([...bodyText.matchAll(REF)].map((m) => `${m[1]} ${m[2]}`))].filter((r) => r !== key);
  const citedBy = secs.filter((s) => {
    if (s.attrs['id'] === id) return false;
    const t = `${renderText(s.verbatim ?? '', false, true)}\n${renderText(s.reading ?? '', false, true)}`;
    return key && [...t.matchAll(REF)].some((m) => `${m[1]} ${m[2]}` === key);
  }).map((s) => s.attrs['data-item']);

  /* ── 在庫: 木の中の項目 → ノード → そこにある宣言の**結論**。
   *    ★`nodeOf` は必ず論文タグつきで引く（番号だけだと別の論文を掴む）。 */
  const tree = itemsInTree();
  const nodeOf = (k) => tree.get(`${paperKey}|${k}`) ?? [];
  const selfNodes = nodeOf(key);
  const refNodes = refs.map((r) => ({ item: r, nodes: nodeOf(r) }));
  // 同じ節（同じ HTML ファイル）で、この項目より**前**にある項目 = 結果が使えるもの
  const priorItems = inFile.slice(0, at)
    .map((s) => ({ label: s.attrs['data-item'] ?? '', key: itemKeyOf(s.attrs['data-item']) }))
    .filter((x) => x.key);
  const sameFileItems = [...new Set(inFile.map((s) => itemKeyOf(s.attrs["data-item"])).filter(Boolean))];
  const neighbourRels = [...new Set(sameFileItems.flatMap((k) => nodeOf(k).map((n) => n.rel)))];

  const di = declIndex();
  const neighbourDecls = neighbourRels.map((rel) => ({
    rel,
    items: [...new Set(sameFileItems.filter((k) => nodeOf(k).some((n) => n.rel === rel)))],
    decls: declConclusionsOf(rel) ?? [],
  }));

  /* ── 在庫（語で引く）。★逐語（＋直前の setup。どちらも原文）の語を宣言**名**に当てる。
   *    よくある語は落とす。★これは**弱い**照会である——記号中心の項目では 0 に近い。
   *    そのときは上の「同じ節のノードの結論」を見ること（正直にそう出す）。 */
  let wordHits = [];
  let usedWords = [];
  if (di) {
    const stop = stopWords();
    const source = `${renderText(sec.verbatim ?? '', false, true)} ${renderText(setup?.verbatim ?? '', false, true)}`;
    const words = [...new Set(source.toLowerCase().match(/[a-z]{5,}/g) ?? [])].filter((w) => !stop.has(w));
    // ★名前の**最後の成分**だけを見る。名前空間（`ABC3.Found.GenEll.…`）まで含めると
    //   `extension` の 1 語で無関係の論文の宣言が大量に釣れる（実測）。
    //   台帳の付随宣言（`.src` / `.needs` …）も落とす——在庫ではない。
    const names = di.rows
      .filter((r) => (r.kind === 'theorem' || r.kind === 'lemma')
        && !/\.(src|needs|nonvacuous|waiting|record|loadBearing|negControl)$/.test(r.full ?? ''))
      .map((r) => ({ r, n: (r.full ?? '').split('.').pop().toLowerCase() }));
    // ★語は**珍しさで重みを付ける**。素の件数で並べると、`finite`(217 件) と
    //   `extension`(20 件) の 2 語一致だけで GenEll の楕円曲線の補題が上位を占めた（実測）。
    //   150 件を超える語（`finite` / `local` …）は落とし、残りを 1/件数 で足す。
    const freq = new Map(words.map((w) => [w, names.filter((x) => x.n.includes(w)).length]));
    usedWords = words.filter((w) => freq.get(w) > 0 && freq.get(w) <= 150);
    const score = new Map();
    const nmatch = new Map();
    for (const w of usedWords) {
      for (const x of names) if (x.n.includes(w)) {
        score.set(x.r, (score.get(x.r) ?? 0) + 1 / freq.get(w));
        nmatch.set(x.r, (nmatch.get(x.r) ?? 0) + 1);
      }
    }
    wordHits = [...score.entries()]
      .filter(([r]) => nmatch.get(r) >= Math.min(2, usedWords.length))  // 1 語だけの一致は雑音
      .sort((a, b) => b[1] - a[1] || a[0].full.localeCompare(b[0].full)).slice(0, 8)
      .map(([r]) => [r, nmatch.get(r)]);
  }

  /* ── 省略の合図。0_Source は gitignore 下なので、無ければ静かに落とす。 */
  let hedge = null;
  try {
    hedge = execFileSync('node', [join(ROOT, 'tools', 'hedge-index.mjs'), '--paper', paperKey, '--item', key ?? item],
      { encoding: 'utf8', maxBuffer: 1 << 24 });
  } catch { hedge = null; }

  const srcTemplate =
    `{ paper := "${paperKey}", pdfPage := ${page}, item := "${key ?? item}", sectionId := "${id}" }`;

  if (asJson) {
    console.log(JSON.stringify({
      paper: paperKey, title: meta.title, sectionId: id, item, itemKey: key, pdfPage: page ? Number(page) : null,
      file: sec.file, sectionTitle: sec.sectionTitle, matchedBy: how,
      verbatim: renderText(sec.verbatim ?? '', false, true),
      verbatimForDocstring: renderText(sec.verbatim ?? '', true),
      reading: sec.reading ? renderText(sec.reading, false, true) : null,
      setup: setup ? { id: setup.attrs['id'], item: setup.attrs['data-item'], verbatim: renderText(setup.verbatim ?? '', false, true) } : null,
      srcTemplate, refs: refNodes, citedBy, selfNodes,
      inventory: neighbourDecls, priorItems, wordHits: wordHits.map(([r, n]) => ({ ...r, matched: n })),
      hedge,
    }, null, 1));
    return;
  }

  /* ── 人と agent が読む形 ─────────────────────────────────────── */
  const O = [];
  O.push(`# 持ち場（まだ .lean が無い）: [${paperKey}] ${item}`);
  O.push('');
  O.push(`原典: **${meta.title}**（${[meta.author, meta.year].filter(Boolean).join(', ')}） / 物理 p.${page}`
    + ` / 記法の危険度 **${meta.notationRisk ?? '未記載'}**`);
  if (meta.notationRisk === 'high') {
    O.push('★**記法の危険度 high** ——`.txt` 抽出で対象の同定が壊れる。locator は PDF 目視で確かめること');
    O.push('（`ResearchPaper/papers.json` の `notationNotes` に何が壊れるかが書いてある）。');
  }
  O.push(`構造化: \`${sec.file}\` の \`#${id}\`${sec.sectionTitle ? `（${sec.sectionTitle}）` : ''}`);
  O.push(`同定: ${how}`);
  if (selfNodes.length) {
    O.push('');
    O.push(`★**この項目は既に木にある**: ${selfNodes.map((n) => `\`${n.rel}\`${n.hasSorry ? '（sorry あり）' : ''}`).join(' / ')}`);
    O.push(`  ⇒ 新規ファイルではない。\`node tools/brief.mjs --node ${selfNodes[0].rel}\` の方が情報が多い。`);
  }
  O.push('');
  O.push('## 1. ★原文の逐語（これを読まずに段取りを書かないこと）');
  O.push('');
  O.push('```');
  O.push(renderText(sec.verbatim ?? '（.verbatim が無い）', false, true));
  O.push('```');
  if (sec.reading) {
    O.push('');
    O.push('**読み**（★原典ではない。我々が構造化のときに付けた注）');
    O.push('');
    for (const l of renderText(sec.reading, false, true).split('\n')) O.push(`> ${l}`);
  }
  if (setup) {
    O.push('');
    O.push(`<details><summary>直前の setup（記号の定義。\`#${setup.attrs['id']}\` = ${setup.attrs['data-item']}）</summary>`);
    O.push('');
    O.push('```');
    O.push(renderText(setup.verbatim ?? '', false, true));
    O.push('```');
    O.push('');
    O.push('</details>');
  }
  O.push('');
  O.push('## 2. `.src` の雛形（★属性から機械で取った。手で書き直さないこと）');
  O.push('');
  O.push('```lean');
  O.push(`def <宣言名>.src : ABC3.Meta.Source :=`);
  O.push(`  ${srcTemplate}`);
  O.push('```');
  O.push('');
  O.push(`★\`pdfPage\` と \`item\` は HTML の \`data-pdf-page\` / \`data-item\` から取っている`);
  O.push(`（\`data-item\` は "${item}"、G1 は「種別+番号」= "${key}" で照合する）。`);
  O.push('');
  O.push('**docstring に貼る逐語**（★`data-txt` を適用済み＝`pdftotext` に見えている形。');
  O.push('この形でないと `check.mjs` の引用照合が落ちる）:');
  O.push('');
  O.push('```lean');
  O.push(`原文 (${paperKey} p.${page}):`);
  for (const l of renderText(sec.verbatim ?? '', true).split('\n')) O.push(`> ${l}`);
  O.push('```');
  O.push('');
  O.push('## 3. 在庫 —— 同じ節で**既に立っている**ノードの結論（`.lean` を直に読んだ）');
  O.push('');
  O.push('★名前ではなく**結論の形**を見ること。`span_smul_uniformizer` は Corollary 6.3 用に');
  O.push('作られたが、結論 `maximalIdeal B = Ideal.span {τ • α}` が Lemma 6.5 の入力だった。');
  if (!neighbourDecls.length) {
    O.push('');
    O.push(`（[${paperKey}] のこの節を名乗るノードは木にまだ無い＝この節の最初の 1 本）`);
  }
  for (const nb of neighbourDecls) {
    O.push('');
    O.push(`### \`${nb.rel}\`（${nb.items.join(' / ')}） — theorem/lemma ${nb.decls.length} 件`);
    O.push('');
    for (const d of nb.decls.slice(0, 20)) O.push(`- \`${d.name}\` : \`${d.concl.slice(0, 150)}\``);
    if (nb.decls.length > 20) {
      O.push(`- …他 ${nb.decls.length - 20} 件（\`grep -n "^theorem\\|^lemma" lean/ABC3/${nb.rel}\`）`);
    }
  }
  O.push('');
  if (!di) O.push('（`.cache/decl-index.txt` が無い。木全体を語で引くには `node tools/decl-index.mjs`）');
  else {
    O.push(`### 逐語の語で索引の**名前**を引いた結果（★弱い照会。語: ${usedWords.join(' ') || 'なし'}）`);
    O.push(`索引の作成: ${di.mtime.toISOString().slice(0, 16).replace('T', ' ')}（古ければ \`node tools/decl-index.mjs\`）`);
    O.push('');
    if (!wordHits.length) {
      O.push('（0 件。★逐語が記号中心の項目では語では引けない——これは「在庫が無い」ではない。');
      O.push('  上の「同じ節のノードの結論」と `.cache/mathlib-index.txt` を見ること）');
    } else {
      for (const [r, n] of wordHits) O.push(`- (${n}語一致) \`${r.full}\` — \`${r.loc}\``);
    }
  }
  O.push('');
  O.push('## 4. 依存（★原文が名指ししているものだけ。推測は足していない）');
  O.push('');
  if (!refNodes.length) O.push('（本文・読みが他の項目を名指ししていない）');
  for (const r of refNodes) {
    const st = r.nodes.length
      ? r.nodes.map((n) => `\`${n.rel}\`${n.hasSorry ? ' ★sorry あり' : ' ✓'}`).join(' / ')
      : '**木に無い**（着手前に立てるか、依存として記録する）';
    O.push(`- ${r.item} → ${st}`);
  }
  O.push('');
  O.push('**同じ節でこの項目より前にある項目**（原文の依存表ではない。記号と結果の出所の候補）');
  if (!priorItems.length) O.push('（無し＝節の先頭）');
  for (const it2 of priorItems) {
    const ns = nodeOf(it2.key);
    O.push(`- ${it2.label} → ${ns.length ? ns.map((n) => `\`${n.rel}\`${n.hasSorry ? ' ★sorry あり' : ' ✓'}`).join(' / ') : '木に無い'}`);
  }
  O.push('');
  O.push('**被依存**（この項目を名指ししている原典の項目）');
  if (!citedBy.length) O.push('（構造化済みの範囲では無し）');
  for (const c of citedBy) O.push(`- ${c}`);
  O.push('');
  O.push('## 5. 省略の合図（`hedge-index.mjs`）');
  O.push('');
  if (!hedge) O.push('（測れなかった。`0_Source` の `.txt` が無いか、論文の鍵が hedge 側に無い）');
  else { O.push('```'); O.push(hedge.trimEnd()); O.push('```'); }
  O.push('');
  O.push('## 許された手順（この順で）');
  O.push('0. ★**上の 1 を読む。** 原文に書いていない段取りをここに求めないこと');
  O.push('   ——この入口は逐語と在庫を並べるだけで、証明の方針は書かない（それが誤りの元だった）。');
  O.push('1. **在庫**を引く（3 の型を見る。名前で引くと落ちる）。');
  O.push('   Mathlib は `.cache/mathlib-index.txt`。★`Unknown constant` は「無い」ではなく');
  O.push('   「import していない」ことが多い（`tools/lean-idioms.md` #68）。');
  O.push('   ★**「mathlib に無い」と書く前に** `node tools/absent-recheck.mjs --try \'<正規表現>\'`。');
  O.push('2. **証明は MCP `lean_check` で通す**（0.01〜1 秒）。ファイルに書くのはその後。');
  O.push('3. 書いたら `lake build <対象モジュール>` のみ。全体ビルドは最後に 1 回。');
  O.push('4. 配管で詰まったら `tools/lean-idioms.md` を引く。新しい失敗形なら 1 行足す。');
  O.push('5. 逸脱（原典の読み替え・前提の追加）をしたら docstring に必ず記録する。');
  console.log(O.join('\n'));
}
