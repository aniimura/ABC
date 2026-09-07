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
/* ★名前つき実体の表は `tools/entities.mjs` に 1 本化した(2026-09-07 メタ第 15 回、M53)。
 *   ★ここが**開けていなかった**もの: 生きた `1_Structured` 69 本で使われる 53 種のうち
 *   `&otimes;`(4 回)・`&eacute;`(3 回)・`&ouml;`(1 回)の **3 種・延べ 8 回**。
 *   `&otimes;`(テンソル積)は数学の記号なので、体裁ではなく**内容の穴**だった。
 *   ★import 束縛はモジュール本体より先に初期化されるので、下の「関数宣言だけ」の
 *     規約(TDZ 回避)とは衝突しない。 */
import { ENTITIES, decodeEntities } from './entities.mjs';
/* ★見出しの形の判定（M59 / M65 / M71 / M76）は `tools/heading-shape.mjs` に 1 本化した
 *   (2026-09-07 メタ第 20 回、M80)。★理由: 見張りが `--audit-proof` の側にしか無く
 *   **ゲートで落ちなかった**ため、`check.mjs --selftest` からも同じ実装を叩けるようにした。
 *   ★`check.mjs` に写すのではなく共有モジュールにしたのは、写すと「写した方を測る」ことになるから。 */
import { headingRe, headingLineShape, itemHeadingRes, pickHeading } from './heading-shape.mjs';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));

const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

const relArg = opt('--node');
const modArg = opt('--mod');
const paperArg = opt('--paper');

// ★`--audit-proof`: `proofParagraphOf`（1b の Proof 段落抽出）の**無音の失敗**を数える。
//   メタ第 14 回（backlog M47-1）。★**出力するだけで、brief の挙動は 1 行も変えない。**
//   ここに相乗りさせた理由: 測る対象が `proofParagraphOf` という **brief.mjs の内部関数**で、
//   別ファイルに写すと「写した方を測っている」ことになるため（`_*.mjs` を増やさない、も兼ねる）。
if (flag('--audit-proof')) { auditProof(flag('--json')); process.exit(0); }

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

/** 名前つき実体の表。★実体は `tools/entities.mjs`(`check.mjs` と共有の 113 種)。
 *  ★ここで**別の表を持たない**こと —— 90 種の写しを持っていた間、
 *  生きた HTML の `&otimes;` `&eacute;` `&ouml;` が延べ 8 回**開かれずに出ていた**。 */
function entTable() { return ENTITIES; }
/** ★名前は数字を含みうる(`sup1` `sup2` `sup3` が legacy に実在する)。
 *  旧実装の `/&([a-zA-Z]+);/` はそれらを 1 件も開かなかった。共有の実装に寄せる。 */
function decodeEnt(s) { return decodeEntities(s); }


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

/**
 * ★**PyMuPDF が大きな演算子を 1 文字に潰す**ときの凡例（2026-09-07 に Yoshida08 で実測）。
 *
 * Y9 の実装者の要請で足した ——「`Q τ∈H(στ(π′)−π′)` の `Q` が `∏` だと**推測**できたから
 * 読めた。凡例が付いていれば推測が要らない」。
 * ★`pdftotext` 側の壊れ方とは**別物**である（あちらは `Σ` を落とす）。
 *
 * ★**単独の大文字は本当にその記号のこともある**（Yoshida は `L/K` の `L` を使う）ので、
 * 曖昧なものは曖昧と書く。
 *
 * ★★**`S`（⋃）と `L`（⊕）は 2026-09-07 に落とした**（メタ第 17 回、backlog M60 / M64）。
 * 標本ではなく**全数**で数えた —— `--audit-proof` が Proof を取れた 160 件を母数に、
 * `S` が出る **18 件**・`L` が出る **15 件**の一致箇所を 1 つずつ原文で読んだ:
 *
 * | 記号 | 凡例が出た項目 | ★本当にその記号 | 偽陽性 |
 * |---|---|---|---|
 * | `S` | 18 | ★**0**（Sylow 部分群 / スキーム `S` / 有限集合 `S`） | 18/18 |
 * | `L` | 15 | ★**0**（体 `L/K` / 直線束 `L`。[BK] の 1 件だけ大演算子だが正体は `Σ` で `⊕` ではない） | 15/15 |
 * | `T` | 7 | **2**（[Yoshida08] `cor-5-12` / `lemma-5-13` の `⋂m≥1`） | 5/7 |
 *
 * ⇒ `S` / `L` は**当たったためしが無い**ので落とす。★**`T` は残す**——
 * 全数 7 件中 2 件（29%）が本物で、これは捨てるには高い
 * （★メタ第 16 回は標本 12 件から「≒2 箇所」と**外挿**していたが、
 *   実数を数えると「7 件中 2 件」= 3 倍以上の当たり率だった）。
 * ★**却下した案**: 「行末に単独で立ち、次の行が添字で始まるものだけ拾う」に絞ると
 * `lemma-5-13` の本物の `T` を落とす（本物 2 → 1）。**絞り込みは退行**なので採らなかった。
 */
function glyphLegend() {
  return [
    { tok: 'Q', means: '∏', re: /(^|\s)Q(\s|$)/m, note: '' },
    { tok: 'P', means: 'Σ', re: /(^|\s)P(\s|$)/m, note: '' },
    { tok: 'T', means: '⋂', re: /(^|\s)T(\s|$)/m, note: '★曖昧（全数 7 件中 本物は 2 件。残りはスキーム / 有限集合の名前 `T`）' },
    { tok: 'b + 大文字', means: 'ハット（`bK` = `K̂`）', re: /(^|\s)b\s*[A-Z]/m, note: '★上付き・下付きのハットは `b` として残る。本文サイズでは丸ごと落ちる' },
    { tok: '∼ = （2 行に割れる）', means: '≅', re: /∼/m, note: '' },
    { tok: '̸=', means: '≠', re: /̸=/m, note: '★`pdftotext` 側では斜線が落ちて `=` に見える' },
  ];
}

/* ★`headingRe`（M51）は `tools/heading-shape.mjs` へ移した（メタ第 20 回 M80）。
 *   ここには**実装を残さない** —— 2 つある状態を作らないため。 */

/**
 * ★**原典の Proof 段落を `0_Source/<file>.txt` から出す。**
 *
 * 2026-09-07 に足した。Y7a と Y7b の実装者が**独立に**「brief に Proof が無く、
 * `.txt` を直読して段落を取ったのが決定的だった」と報告したため
 * （Y7b は 981–1030 行を直読して 12 段の段取りをそのまま得た）。
 *
 * ★**これは発見的な切り出しである。** 当てにならない場合があるので、
 * 必ず**行番号を添えて**出し、読み手が `.txt` を直読できるようにする。
 * ★`0_Source` は gitignore 下なので、無ければ静かに落とす（`hedge` と同じ作法）。
 * ★抽出器は **PyMuPDF**（`tools/source-text.py`）——`pdftotext` とは壊れ方が違う
 * （`≠` は `̸=` で残る／ハットは本文サイズで落ち上付きでは `b` になる）。
 */

/* ★`headingLineShape`（M59）は `tools/heading-shape.mjs` へ移した（メタ第 20 回 M80）。
 *   ★動機と実測（誤った境界 36/36 を落とし、失う真の見出しは [Tate] の 1 個だけ）は移した先に全部ある。
 *   ここには**実装を残さない** —— 2 つある状態を作らないため。 */

function proofParagraphOf(meta, key, sec, item) {
  if (!meta?.file) return null;
  // ★原典に見出しが無い単位（我々が切り出した setup / 暗黙の定義）は、
  //   「見つからない」ではなく「原典に無い」と言わせる（2026-09-07 の実測で 68 件中 17 件がこれ）。
  const implicit = /\b(setup|remark)\b/.test(sec?.attrs?.['class'] ?? '')
    || (sec?.attrs?.['data-implicit'] ?? '') === 'true';
  const p = join(ROOT, 'ResearchPaper', '0_Source', `${meta.file}.txt`);
  if (!existsSync(p)) return { err: 'missing', want: `ResearchPaper/0_Source/${meta.file}.txt` };
  const lines = readFileSync(p, 'utf8').split('\n');
  const HEAD = headingRe();
  // 見出しの同定。`Proposition 6.6` は `Proposition 6.6.` とも `Proposition 6.6 (Sen [14]).` とも出る。
  // ★`itemKeyOf`（check.mjs G1 と共有）は数字番号しか拾わないので、
  //   `Theorem A` のような**英字番号**はここで `data-item` から拾い直す。
  let [kind, num] = (key ?? '').split(/\s+/);
  // ★★2026-09-07（メタ第 12 回 M39）**`kind` が見出し語かを検査していなかった。**
  //   `itemKeyOf`（check.mjs G1 と共有）は `Section` / `Chapter` / `Step` / `Assertion` も拾うので、
  //   `key = "Section 2"` が [pGC] の `.txt` 138 行目の節見出し
  //   「Section 2: Higher Ramification Groups」に当たり、
  //   ★**その手前 13 行が 9 件の項目すべての「証明」として出ていた**（誤報 20 件）。
  //   ⇒ 見出し語でなければ「原典に見出しを持たない単位」として扱う（それが実態）。
  const HEADWORD = /^(Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example|Claim|Fact|Exercise)$/;
  if (kind && !HEADWORD.test(kind)) { kind = null; num = null; }
  if (!kind || !num) {
    const m = /(Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example|Claim)\s+([A-Z])(?![a-zA-Z0-9])/.exec(item ?? '');
    if (m) { kind = m[1]; num = m[2]; } else return { err: 'nokey', implicit };
  }
  // ★M51 / M62 / M64（番号先行の見出し・`(?![a-z])` で閉じる理由）と
  // ★正規表現 2 本の作り方（`want` / `wantDot` と、末尾ドットの落とし穴）は
  //   `tools/heading-shape.mjs` の `itemHeadingRes` に移した（メタ第 20 回 M80）。
  const { want, wantDot } = itemHeadingRes(kind, num);
  // ★★★★★M71（折り返した相互参照を見出しと誤認する壊れ方）と M76（その見張り）の
  //   実装は `tools/heading-shape.mjs` の `pickHeading` に移した（メタ第 20 回 M80）。
  //   ★**「行頭 + 見出しの形」で全部拾い、最初を採る**。`hits` が 2 本以上なら破れうる持ち場。
  //   ★この規則は `check.mjs --selftest` の `txtCases` がゲートで毎回叩く。
  // ★探す順。**下に行くほど信頼度が落ちる**。★最後の段は従来の挙動そのままで、
  //   「形の判定でどれも通らなかったときに退行しない」ための保険である（実測で 337 件がここに残る）。
  const picked = pickHeading(lines, kind, num);   // 1) 行頭に立ち、かつ**見出しの形**
  const shapedAtLines = picked.hits;
  let at = picked.at > 0 ? picked.at - 1 : -1;
  let forcedRunin = false;
  if (at < 0) {                                     // 2) **字下げされた**行頭（2 段組の走査）
    // [Tate] "   Theorem 1. We have…" / [BK] "  Lemma 3.8.1. Let V be…"。
    // ★M66 は「[Tate] の真の見出しは .txt に 1 行も残っていない」と書いたが、★**それは誤りだった**
    //   —— 555 行に字下げで立っている（`^` に当たらないので見えていなかっただけ）。
    // ★段組の走査は本文が左右で混ざるので、`runin` を立てて**信頼度を落として**出す
    //   （`runin` は境界の判定も行中形に切り替えるので、切り出しが 200 行に膨らむのを防ぐ）。
    for (let i = 0; i < lines.length; i++) {
      const t = lines[i].replace(/^\s+/, '');
      const m = wantDot.exec(t);
      if (m && headingLineShape(t, m)) { at = i; forcedRunin = true; break; }
    }
  }
  if (at < 0) at = lines.findIndex((l) => want.test(l));   // 3) 従来どおり（形を問わない）
  const shapedAt = (i) => { if (i < 0) return null; const t = (lines[i] ?? '').replace(/^\s+/, ''); const m = wantDot.exec(t); return !!(m && headingLineShape(t, m)); };
  // ★退避: 見出しが**走り込み式**（行中に出る）論文がある。実測 2026-09-07:
  //   行頭に立つ = BK CorrHyp EtTh GenEll LocProP NCBelyi pGC SemiAnbd Stacks Tate Yoshida08、
  //   行中に走り込む = MilneCFT(154 件) Sharifi(973 件)。
  //   ★このときは切り出しの信頼度が落ちるので `runin` を立てて出力側で断る。
  let runin = forcedRunin;
  if (at < 0) {
    const anywhere = new RegExp(`(?:^|\\s)${kind}\\s+${num.replace(/\./g, '\\.')}(?!\\.?[0-9])`);
    at = lines.findIndex((l) => anywhere.test(l));
    runin = at >= 0;
  }
  if (at < 0) return { err: 'noheading', key: `${kind} ${num}`, implicit };
  // 次の見出しまでの範囲でだけ `Proof` を探す（他の項目の証明を掴まないため）。
  const nextHead = runin
    ? new RegExp(`(?:^|\\s)(Theorem|Proposition|Corollary|Definition|Lemma|Remark|Example|Claim|Fact|Exercise)\\s+[0-9]+(?:\\.[0-9]+)*\\.`)
    : HEAD;
  // ★★**「行頭に見出しが立っている」だけでは境界にならない。** 2026-09-07 に実測:
  //   Yoshida `Corollary 4.9` の証明の中で、折り返しのせいで
  //   `Proposition 4.8. Lemma 4.6 shows bKm` という行ができ、**見出しと誤認されて
  //   証明が 6 行のうち 2 行で切れた**。★引用や折り返しは同じ見出し語を**再出現**させる。
  //   ⇒ **各見出し鍵の「最初の出現」だけを境界として採る**（項目はファイル内で
  //   昇順に一度だけ立つ、という原典の性質を使う）。
  const boundary = new Set();
  {
    const seen = new Set();
    for (let i = 0; i < lines.length; i++) {
      const m = nextHead.exec(lines[i]);
      if (!m) continue;
      // ★★2026-09-07（メタ第 16 回 M59）**行頭に立った「見出しに見える語」の大半は
      //   折り返した相互参照である。** `.txt` は折り返しが固定なので
      //   `… uniformizer of K by` の次の行が `Proposition 4.4(ii), which shows v(β) = qi.`
      //   になり、これが境界と誤認されて **文の途中で証明が切れていた**（実測 54 件中 36 件）。
      if (!headingLineShape(lines[i], m)) continue;
      const k = lines[i].slice(m.index, m.index + m[0].length).trim();
      if (seen.has(k)) continue;      // 2 度目以降は引用か折り返し
      seen.add(k);
      boundary.add(i);
    }
  }
  const isBoundary = (i) => boundary.has(i);
  let stop = lines.length;
  for (let i = at + 1; i < lines.length; i++) if (isBoundary(i)) { stop = i; break; }
  let ps = -1;
  for (let i = at; i < stop; i++) if (/^\s*Proof\b/.test(lines[i])) { ps = i; break; }
  // ★`Proof` も走り込むことがある（`… . Proof. Let …`）。行頭で取れなければ行中を探す。
  if (ps < 0) for (let i = at; i < stop; i++) if (/(?:^|[.\s])(Proof|PROOF)\b/.test(lines[i])) { ps = i; runin = true; break; }
  if (ps < 0) {
    // ★**論証が主張の「手前」にある論文がある。** [pGC] §1 で実測(2026-09-07):
    //   Mochizuki は `Proof:` ブロックを置かず、地の文で導いてから
    //   `Proposition 1.2:` と宣言する。⇒ 直前の地の文こそが証明である。
    //   ★40 行を上限に、直前の見出しは越えずに拾う。
    let lead0 = Math.max(0, at - 40);
    for (let i = at - 1; i >= lead0; i--) if (isBoundary(i)) { lead0 = i + 1; break; }
    const lead = lines.slice(lead0, at).join('\n').replace(/^\s+|\s+$/g, '');
    return {
      err: 'noproof', key: `${kind} ${num}`, headingLine: at + 1, stopLine: stop,
      headShaped: shapedAt(at), headText: (lines[at] ?? '').trim().slice(0, 90),
      shapedHits: shapedAtLines,
      lead: lead || null, leadFrom: lead0 + 1, leadTo: at,
      rel: `ResearchPaper/0_Source/${meta.file}.txt`,
    };
  }
  // 終端: 証明終わりの記号の行（単独でも行末でもよい）／次の見出し／行数の上限（下記 CAP）。
  // ★**記号は論文ごとに違う**（2026-09-07 実測）: Yoshida08 は `□`、
  //   ★Mochizuki は `⃝`（U+20DD 囲み丸。`… assertion (iv). ⃝` のように行末に付く）。
  //   ⇒ GenEll の 16 件が全部「終端なし」だったのはこれが理由だった。
  const END = new RegExp('(\u25A1|\u2293\u2294|\u2294\u2293|\u25A0|\u220E|\u20DD|\u25CB|\u25EF'
    + '|Q\\.\\s?E\\.\\s?D\\.)\\s*$');
  let pe = ps;
  /* ★★上限（メタ第 19 回 M72。★**200 → 800**）
   *
   * 200 だった理由は「暴走を止める」だけで、★**原典の証明の長さを測って決めた数字ではなかった**。
   * 2026-09-07 に 360 件の全数で測った結果:
   *   · 上限に当たっていたのは **6 件**。真の長さは 351 / 274 / 740×4 行。
   *   · 上限を外す（10 万行）と 6 件とも終端記号に当たって `ok` になり、**他の 354 件は 1 件も動かない**。
   *   · ★**800 と 10 万行は今日の 48 本では同じ結果**（1b 総行数 13052 で一致）。
   *     ⇒ 800 は「暴走を止める」役目を保ったまま、今日の原典を**取り足りなく**しない最小級の数字。
   * ★**なぜ 400 のような中間値を採らないか**（実測）: 400 では [IUTchIII] の 4 件（740 行）が
   *   依然 `trunc` のまま **各 +200 行だけ太る** ＝ 1b が +1025 行増えて誰も救われない。
   *   ★中間値は「高くつくのに直らない」ので**支配されている**。
   * ★★決め手（`cor-3-12-step-xi` の実測）: この項目は `Step (xi-e), (xi-f)` を名指しているが、
   *   その本体は **12468–12506 行**にあり、200 行の打ち切り（12061 行）より **407 行うしろ**である。
   *   ⇒ ★**その項目の主題が 1 行も載っていなかった。** 鍵語の数でも 0 対 17。
   * ★`BRIEF_PROOF_CAP` で上書きできる（測り直し用。★本番では立てない）。 */
  const CAP = Number(process.env.BRIEF_PROOF_CAP) || 800;
  const cap = Math.min(ps + CAP, lines.length);
  for (let i = ps; i < cap; i++) {
    pe = i;
    if (i > ps && isBoundary(i)) { pe = i - 1; break; }
    if (END.test(lines[i])) break;
  }
  const text = lines.slice(ps, pe + 1).join('\n').replace(/\s+$/, '');
  const truncated = pe >= cap - 1 && !END.test(lines[pe] ?? '');
  /* ★★上限で切ったとき「**残りが何行あるか**」を数える（メタ第 19 回 M72）。
   *
   * なぜ本文ではなく行数だけを出すのか:
   *   終端を**探す**のは `.txt` を走査するだけで済む（本文を brief に載せる必要が無い）。
   *   ⇒ 1 行で「取り足りない」を**数字で**言える。読み手が「200 行で終わり」と誤解しない。
   * ★実測（2026-09-07、360 件）: 上限に当たっているのは 6 件で、真の長さは
   *   [FrdI] Thm 3.4 = 351 行 / Thm 4.2 = 274 行 / [IUTchIII] Corollary 3.12 = **740 行**（4 件が共有）。
   *   ★`cor-3-12-step-xi` が指す **Step (xi) は 12339 行から**始まる ＝ 上限 200 行の
   *   打ち切り（12061 行）より **278 行うしろ**。★**その項目の本体が 1 行も載っていない。**
   *   ⇒ 「あと何行あるか」を言わないと、読み手はそこを読みに行く理由を持てない。 */
  let restTo = null, restEndFound = false;
  if (truncated) {
    for (let i = cap; i < lines.length; i++) {
      if (isBoundary(i)) { restTo = i; break; }               // 直前の行までが証明
      if (END.test(lines[i])) { restTo = i + 1; restEndFound = true; break; }
    }
    if (restTo === null) restTo = lines.length;               // 終端も見出しも無いまま EOF
  }
  /* ★★見張り（メタ第 21 回 M84。M81 が「直さないと決めた」ことの**上界だけ**を数える）
   *
   * 埋め込み補題（embedded lemma）: 外側の証明の中に `Lemma` / `Claim` が立つと、
   * 境界の規則が**正しく働いた結果として**外側の証明が切れる。実測（M77 / M81）:
   * `[CorrHyp] corrhyp-thm-5-3` は真の長さ 151 行のうち **17 行（11%）**しか出ていない。
   * ★これは `trunc` とは**別の壊れ方**で、★`noend` の中にしか出ない。
   *
   * ★**直さない**（M81 の実測: V1「終端まで走る」は 1 件助けて **21 件を壊す**、
   *   V2 も CorrHyp を直せない）。★**数えるだけ**にする。M76 と同じ「上界を数える」作法。
   *
   * 条件（M84 のとおり 3 つ）:
   *   (a) `noend`（終端記号に当たらず**境界で**切った）
   *   (b) 止めた境界の行が `Lemma` / `Claim` の見出し
   *   (c) ★**その先に終端記号がある**（＝終端記号を持つ論文である。
   *       [Falt1] / [Stacks] / [NCBelyi] / [Tate] は `.txt` に 1 つも無いので落ちる）
   * ★★上界の意味: ここに出た件数だけが V3（M81）の得になりうる。★今日は 1 件。
   *   ★2 件目が出たら V3 の損得が変わる ⇒ そのときに人が判断する。
   * ★出す `embedEndAt` は「境界の先で**最初に**出る終端記号」であって、
   *   ★**外側の証明の終端とは限らない**（CorrHyp では埋め込み補題自身の ⃝ が先に出る）。
   *   上界なので、これでよい。★推測して外側を当てにいかない。 */
  let embedKind = null; let embedEndAt = null;
  if (!truncated && !END.test(lines[pe] ?? '')) {
    const b = lines[pe + 1] ?? '';
    const m = /(?:^|\s)(Lemma|Claim)\s+\d/.exec(b) || /^\s*\d+(?:\.\d+)*\.\s*(Lemma|Claim)\b/.exec(b);
    if (m) {
      embedKind = m[1];
      for (let i = pe + 2; i < lines.length; i++) {
        if (END.test(lines[i])) { embedEndAt = i + 1; break; }
      }
      if (embedEndAt === null) embedKind = null;              // (c) を満たさない
    }
  }
  return {
    embedKind, embedEndAt, embedBoundaryAt: embedKind ? pe + 2 : null,
    text,
    fromLine: ps + 1,
    toLine: pe + 1,
    headingLine: at + 1,
    headShaped: shapedAt(at), headText: (lines[at] ?? '').trim().slice(0, 90),
    shapedHits: shapedAtLines,
    runin,
    endFound: END.test(lines[pe] ?? ''),
    truncated, cap: CAP,
    restTo, restEndFound,
    restLines: restTo === null ? 0 : restTo - (pe + 1),
    rel: `ResearchPaper/0_Source/${meta.file}.txt`,
  };
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

  /* ── ★原典の Proof 段落（2026-09-07 に足した。実装者 2 名が独立に「これが決定打」と報告）。 */
  const proof = proofParagraphOf(meta, key, sec, item);

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
      hedge, proof,
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
  /* ── ★★照合は**項目単位**である。原典の 1 項目が (i)(ii)(iii) に分かれていると
   *    「一部だけ埋まっている」場合に誤解を招く（2026-09-07、Y24 の報告）:
   *    `cor-6-13` で「既に木にある」と出たが、埋まっていたのは (i)(ii) だけで
   *    ★**(iii) は無かった**。実装者は「半分ミスリード」と評した。
   *    ⇒ 逐語に部分番号が見えるときは、それを名指しして断る。 */
  const partRe = /\((i{1,3}|iv|v|vi{1,3}|[1-9])\)/g;
  const partsFound = [...new Set([...(renderText(sec.verbatim ?? '', false, true).match(partRe) ?? [])])];
  if (selfNodes.length) {
    O.push('');
    O.push(`★**この項目は既に木にある**: ${selfNodes.map((n) => `\`${n.rel}\`${n.hasSorry ? '（sorry あり）' : ''}`).join(' / ')}`);
    O.push(`  ⇒ 新規ファイルではない。\`node tools/brief.mjs --node ${selfNodes[0].rel}\` の方が情報が多い。`);
    if (partsFound.length >= 2) {
      O.push('');
      O.push(`★★**ただしこの照合は「項目」単位である。** 逐語に部分番号 ${partsFound.join(' ')} が見える`);
      O.push('  ので、★**一部だけが埋まっている可能性がある**（2026-09-07 に `cor-6-13` で実際に起きた');
      O.push('  ——(i)(ii) は埋まっていたが (iii) は無く、実装者は「半分ミスリード」と評した）。');
      O.push('  ★**上のファイルを直読して、どの部分が埋まっているかを自分で確かめること。**');
    }
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
  O.push('## 1b. ★原典の Proof 段落（`0_Source` の `.txt` から機械で切り出した）');
  O.push('');
  if (!proof) {
    O.push('（この論文の `0_Source` の入口が引けない）');
  } else if (proof.err === 'missing') {
    O.push(`（\`${proof.want}\` が無い。\`0_Source\` は gitignore 下なので手元にしか置かれない）`);
  } else if (proof.err === 'nokey' || (proof.err === 'noheading' && proof.implicit)) {
    O.push('★**この項目は原典に見出しを持たない**——`setup` / `remark` / 暗黙の定義として');
    O.push('  **我々が切り出した単位**なので、Proof 段落が無いのは当然である（欠落ではない）。');
    O.push('  記号の出所を知りたければ、この節の**前後の項目**の Proof を見ること。');
  } else if (proof.err === 'noheading') {
    O.push(`★**\`.txt\` の行頭にも行中にも "${proof.key}" という見出しが見つからなかった。**`);
    O.push('  抽出器が見出しを潰しているか、番号の書式が違う。`.txt` を直に grep すること。');
  } else if (proof.err === 'noproof') {
    O.push(`★**原典はこの項目に \`Proof\` を付けていない**（見出しは \`.txt\` の ${proof.headingLine} 行目、`);
    O.push(`  次の見出しまで ${proof.stopLine - proof.headingLine} 行）。定義・設定・系ではこれが普通である。`);
    O.push('  ★**「証明が無い」＝「自明」ではない。** §5 の省略の合図を見ること。');
    if (proof.lead) {
      O.push('');
      O.push('★★**論証が主張の「手前」にある論文がある。** [pGC] は `Proof:` ブロックを置かず、');
      O.push('地の文で導いてから `Proposition 1.2:` と宣言する（2026-09-07 に実測）。');
      O.push(`直前の地の文（\`${proof.rel}\` の **${proof.leadFrom}–${proof.leadTo} 行**、上限 40 行）:`);
      O.push('');
      O.push('```');
      O.push(proof.lead);
      O.push('```');
      O.push('');
      O.push('★**これが証明とは限らない**——ただの導入かもしれない。行番号から `.txt` を直読して');
      O.push('  どこから論証が始まっているかを自分で決めること。');
    }
  } else {
    O.push(`出所: \`${proof.rel}\` の **${proof.fromLine}–${proof.toLine} 行**`
      + `（見出しは ${proof.headingLine} 行目）${proof.truncated ? ` ★**終端記号に届かず ${proof.cap} 行で打ち切った**` : ''}`);
    if (proof.truncated && proof.restLines > 0) {
      // ★M72（メタ第 19 回）: 「200 行で終わり」と読まれないよう、**残りを数字で**言う。
      O.push('');
      O.push(`★★**この証明はまだ ${proof.restLines} 行つづく**`
        + `（\`${proof.rel}\` の **${proof.toLine + 1}–${proof.restTo} 行**`
        + `${proof.restEndFound ? '、そこで終端記号に当たる' : '、終端記号は見つからない'}）。`);
      O.push(`  ★**上に出したのは全 ${proof.restTo - proof.fromLine + 1} 行中の先頭 ${proof.toLine - proof.fromLine + 1} 行`
        + `（${Math.round((proof.toLine - proof.fromLine + 1) * 100 / (proof.restTo - proof.fromLine + 1))}%）にすぎない。**`);
      O.push('  ★この項目の根拠が後半にあるなら、上の行番号で `.txt` を直読すること。');
    }
    if (proof.runin) {
      O.push('');
      O.push('★★**見出しか `Proof` が行頭に立たず、行中で拾った**（走り込み式の組版）。');
      O.push('  この切り出しは**前後を巻き込んでいる可能性が高い**。行番号から `.txt` を直読すること。');
    } else if (!proof.endFound && !proof.truncated) {
      // ★`truncated` のときは「次の見出しの手前で切った」は**嘘**（切ったのは行数の上限）。
      //   M72 の注記を上で出しているので、ここは重ねない（メタ第 19 回で直した）。
      O.push('');
      O.push('★**終端記号（`□` 等）に当たらず、次の見出しの手前で切った。**');
      O.push('  証明の末尾が落ちているか、次の項目を巻き込んでいるかもしれない。');
    }
    O.push('');
    O.push('```');
    O.push(proof.text);
    O.push('```');
    O.push('');
    // ★凡例（この Proof の中に**実際に出ている**ものだけ挙げる。noise を出さない）
    const hits = glyphLegend().filter((g) => g.re.test(proof.text));
    if (hits.length) {
      O.push('★**この抽出器（PyMuPDF）の字形の潰れ方**——上の引用に実際に出ているものだけ:');
      O.push('');
      O.push('| 見えている形 | 原典 | |');
      O.push('|---|---|---|');
      for (const g of hits) O.push(`| \`${g.tok}\` | ${g.means} | ${g.note} |`);
      O.push('');
      O.push('★**単独の大文字は本当にその文字のこともある**（Yoshida は体の名前に `L` を使う）。');
      O.push('  曖昧な行は行番号から `.txt` を直読して決めること。');
      O.push('');
    }
    O.push('★**これは発見的な切り出しである**（行頭の見出し／`Proof`／証明終端記号で区切っただけ）。');
    O.push('原文の意図が取れないと感じたら**上の行番号から `.txt` を直読すること**——');
    O.push('2026-09-07 に Y7a・Y7b の実装者が独立に「直読が決定打だった」と報告している。');
    O.push('★`.txt` には `===== [page N] =====` が**行として挿入される**ので、');
    O.push('  ページ境界で 1 つの文が割れる（Y9 が「by Lemma / 5.11.」で踏んだ）。');
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

  /* ── ★★「既に木にある」を**末尾にも**出す（2026-09-07）。
   *    本体が `sed -n '/## 1\./,/…/p'` で出力を絞ったせいで冒頭の警告を切り落とし、
   *    ★**既に埋まっている項目（Cor 6.3）を持ち場として配ってしまった**。
   *    実装者が 5 分で気づいて内容を組み替えたが、道具の側で塞げる失敗形である。
   *    ⇒ **冒頭と末尾の両方に出す。** どちらか片方だけを読んでも当たる。
   *    ★`--json` には `selfNodes` として既に入っているので変更しない。 */
  if (selfNodes.length) {
    O.push('');
    O.push('---');
    O.push('');
    O.push('# ★★★警告（冒頭にも出している）: この項目は**既に木にある**');
    O.push('');
    for (const n of selfNodes) O.push(`- \`lean/ABC3/${n.rel}\`${n.hasSorry ? ' ★**sorry あり**' : ' ✓ sorry なし'}`);
    O.push('');
    O.push('★**新規ファイルの持ち場として配られたなら、その前提が誤っている。**');
    O.push('  まず `node tools/brief.mjs --node ' + selfNodes[0].rel + '` を読み、');
    O.push('  **何が既にあり何が無いかを測ってから**進むこと。');
    O.push('  ★**そのまま新規に書き直すと、同じ数学を 2 度形式化することになる。**');
    if (partsFound.length >= 2) {
      O.push('');
      O.push(`★★**ただし照合は「項目」単位である。** 逐語に部分番号 ${partsFound.join(' ')} が見えるので、`);
      O.push('  ★**一部だけが埋まっている可能性がある**（2026-09-07 に `cor-6-13` で実際に起きた）。');
      O.push('  ★**「既に在る」を鵜呑みにせず、どの部分が埋まっているかを直読で確かめること。**');
    }
  }
  console.log(O.join('\n'));
}

/**
 * ★`proofParagraphOf` の**無音の失敗**を数える（メタ第 14 回、backlog M47-1）。
 *
 * 動機: 1b の Proof 段落は「16 人の実装者が役に立ったと報告した」機能だが、
 * **切り出しに失敗しても brief はその節を黙って落とすだけ**である。
 * 失敗の件数も、失敗が集中する論文も、誰も数えたことがなかった。
 *
 * ★結果の読み方（`ok` でも当たっているとは限らない。信頼度の印を分けて数える）:
 *   ok            … Proof ブロックを取れた
 *     /runin      … 見出しか `Proof` が**行中**に走り込んでいた（切り出しの信頼度が落ちる）
 *     /noend      … 証明終わりの記号が見つからず、**次の見出しか行数の上限**で切った
 *     /trunc      … 行数の上限で切った（＝**取り足りない**）。★メタ第 19 回（M72）に
 *                   上限を 200 → 800 にしたので、★**今日の 48 本では 0 件**である。
 *                   （200 のときの 6 件の真の長さ: [FrdI] Thm 3.4 = 351 / Thm 4.2 = 274 /
 *                     [IUTchIII] Corollary 3.12 = 740 行を 4 件が共有）。
 *                   ★1 件でも立ったら「上限が原典に追い付いていない」合図。
 *                   brief 本体には「まだ N 行つづく（X–Y 行）」の 3 行が出る。
 *   noproof/lead  … `Proof` が無く、直前の地の文を代わりに出した（[pGC] 型）
 *   noproof/nolead… `Proof` が無く、代わりも無い（★brief に 1b が出ない）
 *   noheading     … `.txt` に見出しが見つからない（★無音の失敗の本体）
 *   nokey         … `data-item` が見出し語でない（setup など。★これは正常）
 *   missing       … `0_Source/<file>.txt` が無い（gitignore 下。★環境の話）
 */
function auditProof(asJson) {
  const reg = JSON.parse(readFileSync(papersJson(), 'utf8')).papers;
  const rows = [];
  for (const pk of structuredPaperKeys()) {
    const meta = reg[pk];
    for (const sec of collectPaper(pk)) {
      const item = sec.attrs['data-item'] ?? '';
      const r = meta ? proofParagraphOf(meta, itemKeyOf(item), sec, item) : null;
      let outcome;
      if (!meta) outcome = 'nopaper';
      else if (r === null) outcome = 'nofile';
      else if (r.err === 'noproof') outcome = r.lead ? 'noproof/lead' : 'noproof/nolead';
      else if (r.err) outcome = r.err;
      else {
        outcome = 'ok';
        if (r.truncated) outcome += '/trunc';
        else if (!r.endFound) outcome += '/noend';
        if (r.runin) outcome += '/runin';
      }
      // ★`glyphLegend` の被覆も同時に測る（メタ第 14 回、backlog M47-3）:
      //   「実際に出た字だけを出す」表が、何件の持ち場で何行になるか。
      const gl = r && r.text ? glyphLegend().filter((g) => g.re.test(r.text)).map((g) => g.tok) : [];
      // ★診断欄（メタ第 17 回 M65）: `headingLine` が**本当に見出しの行か**を後から機械で検算するため。
      //   `.txt` の行を読み直せば `headingLineShape`（M59）に掛けられる。出力の文面は変えない。
      rows.push({ paper: pk, id: sec.attrs['id'] ?? '', item, outcome,
        lines: r && r.text ? r.text.split('\n').length : 0, glyphs: gl,
        headingLine: (r && r.headingLine) ?? null, rel: (r && r.rel) ?? null,
        // ★診断欄（メタ第 18 回 M71）: **どこからどこまでを出したか**を後から機械で突き合わせるため。
        //   `headingLine` だけでは「見出しは動かずに本文の範囲だけ動いた」退行を見逃す。
        //   ★出力（brief 本体）の文面は変えない。`--audit-proof --json` の診断欄だけが増える。
        from: (r && (r.fromLine ?? r.leadFrom)) ?? null,
        to: (r && (r.toLine ?? r.leadTo)) ?? null,
        leadLines: r && r.lead ? r.lead.split('\n').length : 0,
        headShaped: (r && r.headShaped) ?? null, headText: (r && r.headText) ?? null,
        // ★見張り（メタ第 19 回）: 形が合う行が**何箇所あったか**。2 箇所以上なら
        //   「最初を採る」規則が破れうる持ち場（＝人が見るべき上界）。
        shapedHits: (r && r.shapedHits) ?? [],
        restLines: (r && r.restLines) ?? 0, restTo: (r && r.restTo) ?? null,
        // ★見張り（メタ第 21 回 M84）: 埋め込み補題で切れた**上界**。
        embedKind: (r && r.embedKind) ?? null,
        embedBoundaryAt: (r && r.embedBoundaryAt) ?? null,
        embedEndAt: (r && r.embedEndAt) ?? null });
    }
  }
  if (asJson) { console.log(JSON.stringify(rows, null, 1)); return; }
  const tally = new Map();
  for (const r of rows) tally.set(r.outcome, (tally.get(r.outcome) ?? 0) + 1);
  console.log(`# proofParagraphOf の被覆 —— 構造化された項目 ${rows.length} 件`);
  console.log('');
  for (const [k, v] of [...tally].sort((a, b) => b[1] - a[1])) {
    console.log(`  ${String(v).padStart(4)}  ${(v * 100 / rows.length).toFixed(1).padStart(5)}%  ${k}`);
  }
  console.log('');
  console.log('## 論文ごと（★ok 以外が多い論文＝切り出しが効いていない論文）');
  const byPaper = new Map();
  for (const r of rows) {
    if (!byPaper.has(r.paper)) byPaper.set(r.paper, new Map());
    const m = byPaper.get(r.paper);
    m.set(r.outcome, (m.get(r.outcome) ?? 0) + 1);
  }
  const okOf = (m) => [...m].filter(([k]) => k.startsWith('ok')).reduce((a, [, v]) => a + v, 0);
  const rowsP = [...byPaper].map(([p, m]) => {
    const n = [...m.values()].reduce((a, b) => a + b, 0);
    return { p, n, ok: okOf(m), m };
  }).sort((a, b) => (a.ok / a.n) - (b.ok / b.n) || b.n - a.n);
  for (const { p, n, ok, m } of rowsP) {
    const det = [...m].sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k}=${v}`).join(' ');
    console.log(`  ${p.padEnd(12)} ${String(ok).padStart(3)}/${String(n).padEnd(4)} ok  ${det}`);
  }
  console.log('');
  console.log('## glyphLegend の被覆（Proof を取れた持ち場だけが母数）');
  const withText = rows.filter((r) => r.outcome.startsWith('ok'));
  const withGl = withText.filter((r) => r.glyphs.length);
  console.log(`  Proof を取れた ${withText.length} 件 / うち凡例が出る ${withGl.length} 件` +
    `（${(withGl.length * 100 / Math.max(1, withText.length)).toFixed(0)}%）`);
  const gt = new Map();
  for (const r of withText) for (const g of r.glyphs) gt.set(g, (gt.get(g) ?? 0) + 1);
  for (const [g, n] of [...gt].sort((a, b) => b[1] - a[1])) console.log(`    ${String(n).padStart(4)}  ${g}`);
  console.log('');
  console.log('## ★noheading（`.txt` に見出しが見つからない＝無音の失敗）の明細');
  const bad = rows.filter((r) => r.outcome === 'noheading');
  for (const r of bad) console.log(`  [${r.paper}] ${r.id.padEnd(20)} ${r.item}`);
  if (!bad.length) console.log('  （0 件）');

  /* ★★見張り 1（メタ第 19 回）: M71 が自分から申告した「破れる形」。
   *   「形が合う最初の行」を採る規則は、★**同じ鍵が見出しの形で 2 箇所以上に立つとき**にだけ破れうる。
   *   ⇒ その持ち場を**厳密に**数えて明細を出す。★増えたら人が見る。 */
  console.log('');
  console.log('## ★★見張り: 同じ鍵が「見出しの形」で 2 箇所以上に立つ持ち場');
  console.log('   （M71 の危険側。「最初を採る」規則が破れうるのはここだけ ＝ **上界**）');
  const amb = rows.filter((r) => (r.shapedHits ?? []).length >= 2);
  const ambPapers = new Set(amb.map((r) => r.paper));
  console.log(`  ${amb.length} 件 / ${rows.length} 件（${ambPapers.size} 本の論文）`);
  const seen = new Set();
  for (const r of amb) {
    const k = `${r.paper}|${r.shapedHits.join(',')}`;
    if (seen.has(k)) continue;                    // 同じ見出しを共有する項目は 1 行にまとめる
    seen.add(k);
    const others = r.shapedHits.filter((n) => n !== r.headingLine);
    console.log(`  [${r.paper}] ${r.item.padEnd(24)} 採用 L${r.headingLine} / ほか L${others.join(' L')}`);
  }
  if (!amb.length) console.log('  （0 件）');

  /* ★★見張り 2（メタ第 19 回 M72）: 行数の上限に当たった持ち場と、その「残り」。
   *   ★上限を 800 にしたので今日は 0 件。1 件でも立ったら上限が原典に追い付いていない。 */
  console.log('');
  console.log('## ★★見張り: 行数の上限で切れている持ち場（`ok/trunc`）と残りの行数');
  const tr = rows.filter((r) => r.outcome.includes('trunc'));
  for (const r of tr) {
    console.log(`  [${r.paper}] ${r.id.padEnd(22)} 出した ${r.lines} 行 / ★残り ${r.restLines} 行`
      + `（〜${r.restTo} 行。全 ${r.lines + r.restLines} 行）`);
  }
  if (!tr.length) console.log('  （0 件）');

  /* ★★見張り 3（メタ第 21 回 M84）: 埋め込み補題で外側の証明が切れている持ち場の**上界**。
   *   ★M81 が「直さない」と決めた壊れ方を、**数えるところだけ**残す。 */
  console.log('');
  console.log('## ★★見張り: 埋め込み補題で切れた持ち場（`ok/noend` の中。★**上界**）');
  console.log('   （M77/M81。境界が `Lemma`/`Claim` で、かつその先に終端記号がある ＝ V3 の得になりうる件数）');
  const emb = rows.filter((r) => r.embedKind && r.outcome.includes('noend'));
  console.log(`  ${emb.length} 件 / ${rows.length} 件`);
  for (const r of emb) {
    console.log(`  [${r.paper}] ${r.id.padEnd(22)} 出した ${r.lines} 行（L${r.from}-${r.to}）`
      + ` / 境界 L${r.embedBoundaryAt} は ${r.embedKind} / 先の終端記号 L${r.embedEndAt}`);
  }
  if (!emb.length) console.log('  （0 件）');
}
