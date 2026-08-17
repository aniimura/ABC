// 「壁」から到達不能になっている [FrdI] の項目を数える。
//
// ★動機(2026-08-18): 残り 24 件のうち、いくつかは**在庫が無い**ために
//   どうやっても閉じない。その伝播を機械で測る。
//
// 壁(2026-08-18 時点、`tools/index-html.mjs` の逸脱表と `Gap/` の記録に対応):
//   - Definition 2.8  …… mathlib に pro-l 群が無い(位相的有限生成副有限アーベル群の分解)
//   - Lemma 6.5       …… mathlib に six exponentials theorem が無い
//   - Proposition 4.4 …… (ii) の `otriBase`。model / birat-Frobenius-normalized では閉じたが一般は未
//
// ★これは pdftotext 経由の依存抽出に乗っているので、当たりを付けるためだけに使う。
//
// 使い方: node tools/frdi-blocked.mjs

import { readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, '..');
const TXT = path.join(REPO, 'ResearchPaper', '0_Source',
  'The Geometry of Frobenioids I.txt');
const NEED = path.join(REPO, 'ResearchPaper', 'frdi-needed.json');

const WALLS = new Map([
  ['Definition 2.8', 'mathlib に pro-l 群が無い'],
  ['Lemma 6.5', 'mathlib に six exponentials theorem が無い'],
  ['Proposition 4.4', '(ii) の otriBase(一般の 𝒞)'],
]);

const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const HEAD = new RegExp(`^\\s*(${KIND})\\s+(\\d+\\.\\d+(?:\\.\\d+)?)\\.`);
const CITE = new RegExp(`\\b(${KIND})\\s+(\\d+\\.\\d+(?:\\.\\d+)?)`, 'g');

const lines = readFileSync(TXT, 'utf8').split(/\r?\n/);
const heads = [];
lines.forEach((ln, i) => {
  const m = ln.match(HEAD);
  if (m) heads.push([i, `${m[1]} ${m[2]}`]);
});
const body = new Map();
heads.forEach(([i, name], j) => {
  const end = j + 1 < heads.length ? heads[j + 1][0] : lines.length;
  body.set(name, (body.get(name) ?? '') + '\n' + lines.slice(i, end).join('\n'));
});

const needed = JSON.parse(readFileSync(NEED, 'utf8')).needed;
const known = new Set(needed.map((x) => x.item));
const deps = new Map(
  needed.map((x) => {
    const d = new Set();
    for (const m of (body.get(x.item) ?? '').matchAll(CITE)) {
      const n = `${m[1]} ${m[2]}`;
      if (n !== x.item && known.has(n)) d.add(n);
    }
    return [x.item, d];
  }),
);

// 壁から上流へ伝播(「壁に依存する項目」は閉じない)
const blocked = new Map();
for (const [w, why] of WALLS) blocked.set(w, `★壁: ${why}`);
let changed = true;
while (changed) {
  changed = false;
  for (const x of needed) {
    // ★済んだ項目は「閉じた」ことが証拠なので、壁の伝播に乗せない
    if (x.done || blocked.has(x.item)) continue;
    for (const d of deps.get(x.item)) {
      if (blocked.has(d)) {
        blocked.set(x.item, `${d} を経由`);
        changed = true;
        break;
      }
    }
  }
}

const todo = needed.filter((x) => !x.done);
const bad = todo.filter((x) => blocked.has(x.item));
const ok = todo.filter((x) => !blocked.has(x.item));

console.log(`★[FrdI] 残り ${todo.length} 件のうち`);
console.log(`  ★★壁に阻まれている: ${bad.length} 件`);
for (const x of bad) {
  console.log(`     §${x.section} p.${String(x.page).padStart(3)} ${x.item.padEnd(18)} ${blocked.get(x.item)}`);
}
console.log(`  ★手を動かせば閉じる: ${ok.length} 件`);
for (const x of ok) {
  console.log(`     §${x.section} p.${String(x.page).padStart(3)} ${x.item}`);
}

// 節ごとの到達可能な上限
const bySec = new Map();
for (const x of needed) {
  const s = bySec.get(x.section) ?? { total: 0, done: 0, blocked: 0 };
  s.total += 1;
  if (x.done) s.done += 1;
  else if (blocked.has(x.item)) s.blocked += 1;
  bySec.set(x.section, s);
}
console.log('\n★節ごとの「到達可能な上限」(= 全体 − 壁に阻まれた数)');
for (const [sec, s] of [...bySec].sort()) {
  console.log(
    `  §${sec}  現在 ${s.done}/${s.total}   到達可能な上限 ${s.total - s.blocked}/${s.total}` +
      (s.blocked ? `   (壁 ${s.blocked} 件)` : ''),
  );
}
