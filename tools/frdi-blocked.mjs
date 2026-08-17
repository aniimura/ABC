// **まだ閉じていない 3 つの塊**の下流に何件あるかを数える。
//
// ★動機(2026-08-18): 残り 24 件のうち、いくつかは 3 つの未着手の塊の下流にある。
//   その伝播を機械で測る。
//
// ★★**訂正(2026-08-18 夕)**: この道具は当初「壁」と呼び、
//   節ごとの「到達可能な上限」を印字していた。**その読み方は捨てる。**
//   CLAUDE.md の姿勢——「工数の山を『壁』と呼ばない。既知数学の person-years は
//   壁でなく道」——に従い、3 つとも `ResearchPaper/frdi-decomposition.json` の
//   **チェーン**(内部の小目標の DAG)に割った。
//   ★葉と層は `node tools/frdi-newleaves.mjs` が印字する。
//   **この道具が印字するのは「そのチェーンを閉じると何件動くか」であって、
//   到達不能の証明ではない。**
//
// 3 つの塊(`ResearchPaper/frdi-decomposition.json` のチェーンに対応):
//   - Definition 2.8  …… 副有限アーベル群の pro-l 分解        (チェーン prol、葉 4)
//   - Lemma 6.5       …… six exponentials theorem            (チェーン sixexp、葉 3)
//   - Proposition 4.4 …… (ii) の `otriBase`(一般の 𝒞)       (チェーン otricomm、葉 2 は済)
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
  ['Definition 2.8', 'チェーン prol(pro-l 分解)'],
  ['Lemma 6.5', 'チェーン sixexp(six exponentials)'],
  ['Proposition 4.4', 'チェーン otricomm((ii) の otriBase)'],
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
console.log(`  ★★3 つのチェーンの下流にある: ${bad.length} 件`);
for (const x of bad) {
  console.log(`     §${x.section} p.${String(x.page).padStart(3)} ${x.item.padEnd(18)} ${blocked.get(x.item)}`);
}
console.log(`  ★チェーンに依らず単独で閉じる: ${ok.length} 件`);
for (const x of ok) {
  console.log(`     §${x.section} p.${String(x.page).padStart(3)} ${x.item}`);
}

// 節ごとの「いまチェーン待ちの数」
const bySec = new Map();
for (const x of needed) {
  const s = bySec.get(x.section) ?? { total: 0, done: 0, blocked: 0 };
  s.total += 1;
  if (x.done) s.done += 1;
  else if (blocked.has(x.item)) s.blocked += 1;
  bySec.set(x.section, s);
}
console.log('\n★節ごとの「いまチェーン待ちの数」');
console.log('  ★★これは到達不能の証明ではない。3 つのチェーンはいずれも既知数学であり、');
console.log('     葉と層は `node tools/frdi-newleaves.mjs` が印字する。');
for (const [sec, s] of [...bySec].sort()) {
  console.log(
    `  §${sec}  現在 ${s.done}/${s.total}   チェーンに依らず届く ${s.total - s.blocked}/${s.total}` +
      (s.blocked ? `   (チェーン待ち ${s.blocked} 件)` : ''),
  );
}
