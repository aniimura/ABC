// 残っている [FrdI] の項目のうち「葉」(依存がすべて実装済みのもの)を測る。
//
// ★動機(2026-08-18): 進め方は「葉から形式化する」だが、
//   **どれが葉かは物理ページ順では分からない**。実際、この測定で
//   `Theorem 5.2`(§5 の入口・依存 0)が葉であり、しかも `Proposition 5.3`・
//   `Example 6.1`・`Example 6.3`・`Theorem 6.2` がそこで詰まっていると分かった。
//
// ★これは pdftotext 経由の機械抽出なので**当たりを付けるためだけ**に使う。
//   原文が番号で挙げない依存(§0 の語彙・「immediate」の段)は数えられない。
//
// 使い方: node tools/frdi-leaves.mjs

import { readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, '..');
const TXT = path.join(REPO, 'ResearchPaper', '0_Source',
  'The Geometry of Frobenioids I.txt');
const NEED = path.join(REPO, 'ResearchPaper', 'frdi-needed.json');

const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const HEAD = new RegExp(`^\\s*(${KIND})\\s+(\\d+\\.\\d+(?:\\.\\d+)?)\\.`);
const CITE = new RegExp(`\\b(${KIND})\\s+(\\d+\\.\\d+(?:\\.\\d+)?)`, 'g');

const lines = readFileSync(TXT, 'utf8').split(/\r?\n/);

// 見出しの位置 → 各項目の本文(次の見出しまで)
const heads = [];
lines.forEach((ln, i) => {
  const m = ln.match(HEAD);
  if (m) heads.push([i, `${m[1]} ${m[2]}`]);
});
const body = new Map();
heads.forEach(([i, name], j) => {
  const end = j + 1 < heads.length ? heads[j + 1][0] : lines.length;
  const text = lines.slice(i, end).join('\n');
  body.set(name, (body.get(name) ?? '') + '\n' + text);
});

const needed = JSON.parse(readFileSync(NEED, 'utf8')).needed;
const done = new Map(needed.map((x) => [x.item, x.done]));
const todo = needed.filter((x) => !x.done);

const rows = todo.map((x) => {
  const deps = new Set();
  for (const m of (body.get(x.item) ?? '').matchAll(CITE)) {
    const d = `${m[1]} ${m[2]}`;
    if (d !== x.item && done.has(d)) deps.add(d);
  }
  const undone = [...deps].filter((d) => !done.get(d)).sort();
  return { ...x, deps: deps.size, undone };
});

rows.sort((a, b) => a.undone.length - b.undone.length || a.page - b.page);

console.log(`★[FrdI] の残り ${todo.length} 件 —— 未実装の依存が少ない順`);
console.log('  (依存は原文が番号で挙げたもののみ。必要分に入っているものだけ数える)');
for (const r of rows) {
  const mark = r.undone.length === 0 ? '★葉' : '  ';
  console.log(
    `${mark} §${r.section} p.${String(r.page).padStart(3)} ${r.item.padEnd(18)}` +
      ` 未実装 ${r.undone.length}/${r.deps}` +
      (r.undone.length ? `  ← ${r.undone.join(', ')}` : ''),
  );
}

// 「その葉を開けると何が動くか」——各葉を実装したときに未実装依存が減る項目の数
console.log('\n★葉を1つ開けたときに「未実装依存が減る」項目の数(=波及)');
const leaves = rows.filter((r) => r.undone.length === 0);
const impact = leaves
  .map((l) => ({
    item: l.item,
    n: rows.filter((r) => r.undone.includes(l.item)).length,
  }))
  .sort((a, b) => b.n - a.n);
for (const i of impact) console.log(`   ${i.item.padEnd(18)} ${i.n} 件`);
