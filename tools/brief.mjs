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
 * 使い方:
 *   node tools/brief.mjs --node Skeleton/PGC/Section1.lean
 *   node tools/brief.mjs --mod ABC3.Skeleton.PGC.Section1
 *   node tools/brief.mjs --node ... --json
 */

import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));

const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

const relArg = opt('--node');
const modArg = opt('--mod');
if (!relArg && !modArg) {
  console.error('使い方: node tools/brief.mjs --node <Skeleton/…/X.lean> | --mod <ABC3.…>');
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
const rdeps = new Map();
for (const n of nodes) for (const im of n.imports) {
  if (!byMod.has(im)) continue;
  if (!rdeps.has(im)) rdeps.set(im, []);
  rdeps.get(im).push(n.mod);
}
const reach = (start, edgesOf) => {
  const seen = new Set(); const st = [...(edgesOf(start) ?? [])];
  while (st.length) { const m = st.pop(); if (seen.has(m)) continue; seen.add(m);
    for (const w of edgesOf(m) ?? []) if (!seen.has(w)) st.push(w); }
  return seen;
};
const downOf = (m) => rdeps.get(m) ?? [];
const upOf = (m) => (byMod.get(m)?.imports ?? []).filter((x) => byMod.has(x));
const sorryMods = new Set(nodes.filter((n) => n.hasSorry).map((n) => n.mod));

const down = reach(node.mod, downOf);
const up = reach(node.mod, upOf);
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
  blockers: blockers.map((m) => byMod.get(m)?.rel ?? m),
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
  for (const b of result.blockers) L.push(`- \`${b}\``);
  L.push('');
  L.push('☆このノードには**まだ着手できない**。上流を先に片付けること。');
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
L.push('2. **証明は MCP `lean_check` で通す**（0.01〜1 秒）。ファイルに書くのはその後。');
L.push('3. 書いたら `lake build <対象モジュール>` のみ。全体ビルドは最後に 1 回。');
L.push('4. 配管で詰まったら `tools/lean-idioms.md` を引く。新しい失敗形なら 1 行足す。');
L.push('5. 逸脱（原典の読み替え・前提の追加）をしたら docstring に必ず記録する。');
console.log(L.join('\n'));
