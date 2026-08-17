// 「壁」を割って出てきた**新しい葉**を数える。
//
// ★動機(2026-08-18)
//   前日まで 3 つの塊を『壁』と呼んで到達不能と報告していた。CLAUDE.md の
//   「進め方」は、必要なものが出たらスケルトンを足してグラフを更新し、
//   **新しい葉から**形式化することで「工数の大きな塊を壁のように認識しない」
//   ことを目的としている。これはその手続きを、原典の項目より**下の粒度**で
//   実行するための道具である。
//
// ★入力  ResearchPaper/frdi-decomposition.json(我々の設計。原典にも mathlib にも無い)
// ★出力  チェーンごとの層、いま着手できる葉、そのチェーンが開ける [FrdI] 項目数
//
// ★限界: 分解の粒度と依存の向きは**作業のための仮説**である。実装すると変わる。
//   「葉」は依存が無いという意味であって、易しいという意味ではない(size を見ること)。
//
// 使い方: node tools/frdi-newleaves.mjs [--json]

import { readFileSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, '..');
const DEC = path.join(REPO, 'ResearchPaper', 'frdi-decomposition.json');

const doc = JSON.parse(readFileSync(DEC, 'utf8'));

/** 依存の推移閉包で層番号を出す。循環があれば名指しして落とす。 */
function layers(nodes) {
  const byId = new Map(nodes.map((n) => [n.id, n]));
  for (const n of nodes) {
    for (const d of n.deps) {
      if (!byId.has(d)) throw new Error(`未知の依存: ${n.id} -> ${d}`);
    }
  }
  const lay = new Map();
  let changed = true;
  let guard = nodes.length + 2;
  while (changed && guard-- > 0) {
    changed = false;
    for (const n of nodes) {
      const ds = n.deps.map((d) => lay.get(d));
      if (ds.some((x) => x === undefined)) continue;
      const v = ds.length ? Math.max(...ds) + 1 : 0;
      if (lay.get(n.id) !== v) { lay.set(n.id, v); changed = true; }
    }
  }
  const missing = nodes.filter((n) => !lay.has(n.id));
  if (missing.length) throw new Error(`循環がある: ${missing.map((n) => n.id).join(', ')}`);
  return lay;
}

/** 未着手か(done / inMathlib / 半済 は着手対象ではない) */
const isTodo = (n) => n.status === 'todo';

const chains = doc.chains.map((c) => {
  const lay = layers(c.nodes);
  const done = new Set(c.nodes.filter((n) => !isTodo(n)).map((n) => n.id));
  // ★「いま着手できる」= 未着手で、依存がすべて済んでいるもの(層 0 とは限らない)
  const ready = c.nodes.filter((n) => isTodo(n) && n.deps.every((d) => done.has(d)));
  return { c, lay, ready, todo: c.nodes.filter(isTodo) };
});

if (process.argv.includes('--json')) {
  console.log(JSON.stringify(
    chains.map(({ c, lay, ready, todo }) => ({
      id: c.id, title: c.title, unblocks: c.unblocks.length,
      todo: todo.length,
      ready: ready.map((n) => ({ id: n.id, title: n.title, size: n.size })),
      nodes: c.nodes.map((n) => ({ id: n.id, layer: lay.get(n.id), status: n.status })),
    })), null, 2));
  process.exit(0);
}

console.log(`★[FrdI] 塊を割って出た新しい葉(${doc.measuredAt} の設計)\n`);

// ★開ける項目数の多い順(= 律速の順)
chains.sort((a, b) => b.c.unblocks.length - a.c.unblocks.length);

for (const { c, lay, ready, todo } of chains) {
  const depth = Math.max(...c.nodes.map((n) => lay.get(n.id))) + 1;
  console.log(`──────────────────────────────────────────────────────────`);
  console.log(`★チェーン ${c.id} —— ${c.title}`);
  console.log(`   奉じる項目: [${c.serves.paper}] ${c.serves.item}(物理 p.${c.serves.page})`);
  console.log(`   ★これが閉じると [FrdI] の ${c.unblocks.length} 件が動く: ${c.unblocks.join(' / ')}`);
  console.log(`   節点 ${c.nodes.length} 個・層 ${depth} 段・未着手 ${todo.length} 個\n`);
  for (let L = 0; L < depth; L++) {
    const here = c.nodes.filter((n) => lay.get(n.id) === L);
    if (!here.length) continue;
    console.log(`   層${L}`);
    for (const n of here) {
      const mark = n.status === 'done' ? '済  '
        : n.status === 'inMathlib' ? '在庫'
        : n.status === '半済' ? '半済'
        : ready.includes(n) ? '★葉' : '    ';
      console.log(`     ${mark} ${n.id.padEnd(20)} ${n.title}`);
      if (isTodo(n)) console.log(`          大きさ: ${n.size}`);
    }
  }
  console.log('');
}

console.log(`──────────────────────────────────────────────────────────`);
console.log('★★いま着手できる葉(依存がすべて済んでいる未着手の節点)\n');
const all = chains.flatMap(({ c, ready }) => ready.map((n) => ({ c, n })));
for (const { c, n } of all) {
  console.log(`  [${c.id}] ${n.id.padEnd(20)} ${n.title}`);
  console.log(`           大きさ ${n.size} / このチェーンは ${c.unblocks.length} 件を開ける`);
}
console.log(`\n  合計 ${all.length} 葉。`);
console.log('  ★「葉」は依存が無いという意味であって、易しいという意味ではない。');
console.log('  ★分解そのものが我々の設計であり、実装すると変わりうる(限界は JSON の ★限界)。');
