#!/usr/bin/env node
// 論文の番号付き項目(Kind N.M[.K])を §ごとに数える(新しい論文のゴール設定用)。
// ★pdftotext(.txt)由来であり、目視していない。「桁を見るための数」であって、
//   着手可能性や逐語の正しさは保証しない(事実2、PLAN.md)。
// 宣言規則は tools/intra-graph.mjs と同じ:
//   行頭 + Kind + 空白 + N(.M)+ + ピリオド + (行末 または 開き括弧)
// これを緩めると継続行(例: "Definition 3.1. Let X ...")を宣言と誤検出する。
//
// 使い方:
//   node tools/paper-items.mjs <papers.json のタグ>            # 例: LocProP
//   node tools/paper-items.mjs --path "<0_Source の .txt への相対/絶対パス>"
import { readFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const args = process.argv.slice(2);

let path;
let tag = null;
const pathIdx = args.indexOf('--path');
if (pathIdx >= 0) {
  path = args[pathIdx + 1];
} else {
  tag = args[0];
  if (!tag) {
    console.error('使い方: node tools/paper-items.mjs <papers.json のタグ>  または  --path <.txt>');
    process.exit(1);
  }
  const papers = JSON.parse(readFileSync(join(ROOT, 'ResearchPaper', 'papers.json'), 'utf8')).papers;
  const entry = papers[tag];
  if (!entry) {
    console.error(`papers.json に無いタグ: ${tag}`);
    process.exit(1);
  }
  path = join(ROOT, 'ResearchPaper', '0_Source', `${entry.file}.txt`);
}
if (!existsSync(path)) {
  console.error(`.txt が無い(0_Source は再配布しないので各自 pdftotext -layout で作る): ${path}`);
  process.exit(1);
}

const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example|Claim|Fact|Convention';
const txt = readFileSync(path, 'utf8');
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

const sectionAt = [];
for (let i = 0; i < lines.length; i++) {
  const m = lines[i].match(/^Section ([0-9a]+):\s*(.*)$/);
  if (m) sectionAt.push({ line: i, page: pageOf[i], num: m[1], title: m[2].trim() });
}

const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+([0-9a]+(?:\\.[0-9]+)+)\\.[ \\t]*(?:$|\\()`);
const decls = [];
for (let i = 0; i < lines.length; i++) {
  const m = lines[i].match(declRe);
  if (m) decls.push({ kind: m[1], num: m[2], line: i, page: pageOf[i] });
}
const seen = new Set();
const uniq = [];
for (const d of decls) {
  const key = `${d.kind} ${d.num}`;
  if (seen.has(key)) continue;
  seen.add(key);
  uniq.push(d);
}

const bySec = new Map();
for (const d of uniq) {
  const sec = d.num.split('.')[0];
  if (!bySec.has(sec)) bySec.set(sec, []);
  bySec.get(sec).push(d);
}

console.log(`★${tag ?? path}`);
console.log(`論文内の番号付き宣言(重複除去): ${uniq.length} 件`);
if (sectionAt.length) console.log(`Section 見出し: ${sectionAt.map((s) => s.num).join(', ')}`);
console.log('');

const secNums = [...bySec.keys()].sort((a, b) => {
  const na = parseFloat(a.replace(/^a/, '')), nb = parseFloat(b.replace(/^a/, ''));
  if (a.startsWith('a') !== b.startsWith('a')) return a.startsWith('a') ? 1 : -1;
  return na - nb;
});
const goalParts = [];
for (const sec of secNums) {
  const items = bySec.get(sec).sort((a, b) => a.line - b.line);
  const hdr = sectionAt.find((s) => s.num === sec);
  goalParts.push(`§${sec} 0/${items.length}`);
  console.log(`§${sec}${hdr ? ` "${hdr.title}"` : ''}  —— ${items.length} 件  (p.${items[0].page}-${items[items.length - 1].page})`);
  for (const it of items) console.log(`    ${it.kind} ${it.num}  (p.${it.page})`);
}
console.log(`\n★ゴール表記案: ${tag ?? ''} ${goalParts.join(', ')} —— 合計 0/${uniq.length}`);
