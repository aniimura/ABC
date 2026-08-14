// 同一論文内の依存グラフを .txt から機械構築する(親の測定用、使い捨て)
// ★これは pdftotext 経由なので目視していない。規模の当たりを付けるためだけに使う。
import { readFileSync } from 'node:fs';

const path = process.argv[2];
const root = process.argv[3] ?? 'Corollary 3.12';
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';

const txt = readFileSync(path, 'utf8');

// ── 1. 各行に物理ページを振る
// ★CRLF を落とす。落とさないと行末の \r で `$` が効かず、宣言が 0 件になる(実測)。
const lines = txt.split(/\r?\n/);
const pageOf = new Array(lines.length).fill(0);
{
  let pg = 0;
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(/^===== \[page (\d+)\]/);
    if (m) pg = Number(m[1]);
    pageOf[i] = pg;
  }
}

// ── 2. 宣言(行頭の `Kind N.M.`)を拾い、その項目の行範囲を決める
// ★`\.` の後ろに空白か行末を要求する。これが無いと貪欲な後戻りで
//   "Remark 4.6.2," が「Remark 4.6 の宣言」として誤検出され、本体が途中で切れる
//   (2026-08-15 実測。IUTchII Cor 4.10 が到達 1 件になっていた原因)。
// ★宣言行の条件は「行頭 + 番号 + ピリオド + (行末 または 括弧つき表題)」。
//   これを緩めると継続行を宣言と誤検出して本体が途中で切れる。実測した誤検出2種:
//     "Remark 4.6.2, ..."   → 貪欲な後戻りで "Remark 4.6" の宣言に見える
//     "Definition 3.1. Let" → ピリオド+空白だけを条件にすると宣言に見える
//   前者は IUTchII Cor 4.10 の到達を 1 件に、後者も同じ箇所を 4 行に切っていた。
const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
const decls = [];
for (let i = 0; i < lines.length; i++) {
  const m = lines[i].match(declRe);
  if (m) decls.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] });
}
// 同名の重複は最初の宣言を採る
const declOf = new Map();
for (const d of decls) if (!declOf.has(d.name)) declOf.set(d.name, d);
// 範囲 = 次の宣言まで
for (let i = 0; i < decls.length; i++) decls[i].end = i + 1 < decls.length ? decls[i + 1].line : lines.length;

// ── 3. 各項目の本文から、同一論文内の参照を拾う
const refRe = new RegExp(`(?<!\\]\\s?)(?<!\\],\\s)(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
const edges = new Map(); // name -> Set(name)
for (const d of decls) {
  if (edges.has(d.name)) continue;              // 最初の宣言のぶんだけ
  const body = lines.slice(d.line, d.end).join('\n');
  const out = new Set();
  for (const m of body.matchAll(refRe)) {
    const to = `${m[1]} ${m[2]}`;
    if (to !== d.name && declOf.has(to)) out.add(to);
  }
  edges.set(d.name, out);
}

// ── 4. root から推移閉包
const seen = new Set();
const depth = new Map([[root, 0]]);
const queue = [root];
const missing = new Set();
while (queue.length) {
  const cur = queue.shift();
  if (seen.has(cur)) continue;
  seen.add(cur);
  const outs = edges.get(cur);
  if (!outs) { missing.add(cur); continue; }
  for (const to of outs) {
    if (!depth.has(to)) { depth.set(to, depth.get(cur) + 1); queue.push(to); }
  }
}

// ── 5. 出力
const byDepth = new Map();
for (const n of seen) {
  const d = depth.get(n) ?? -1;
  if (!byDepth.has(d)) byDepth.set(d, []);
  byDepth.get(d).push(n);
}
console.log(`root: ${root}`);
console.log(`論文内の番号付き項目: ${declOf.size} 件`);
console.log(`root から到達可能: ${seen.size} 件 / 最大深さ ${Math.max(...depth.values())}`);
for (const d of [...byDepth.keys()].sort((a, b) => a - b)) {
  const ns = byDepth.get(d);
  console.log(`  深さ ${d}: ${ns.length} 件  ${ns.slice(0, 6).join(' / ')}${ns.length > 6 ? ' …' : ''}`);
}
const totalEdges = [...seen].reduce((s, n) => s + (edges.get(n)?.size ?? 0), 0);
console.log(`到達部分の辺: ${totalEdges} 本 / 平均出次数 ${(totalEdges / seen.size).toFixed(1)}`);
