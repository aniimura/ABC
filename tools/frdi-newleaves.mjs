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

/** ★依存の推移閉包で層番号を出す。**鎖をまたぐ依存**も解決する
 *  (2026-08-20: 因子論の展開で `cartier` → `weil` のような鎖跨ぎが出たため。
 *   それまでは鎖の中だけで解決しており、跨いだ依存は「未知の依存」で落ちていた)。 */
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

/** ★全鎖の節点を 1 つのグラフとして層を出す。id の重複は名指して落とす。 */
const ALL_NODES = doc.chains.flatMap((c) => c.nodes.map((n) => ({ ...n, chain: c.id })));
{
  const seen = new Map();
  for (const n of ALL_NODES) {
    if (seen.has(n.id)) throw new Error(`節点 id が重複している: ${n.id}`);
    seen.set(n.id, n.chain);
  }
}
const GLOBAL_LAY = layers(ALL_NODES);
const GLOBAL_BY_ID = new Map(ALL_NODES.map((n) => [n.id, n]));

/** 未着手か(done / inMathlib / 半済 は着手対象ではない) */
const isTodo = (n) => n.status === 'todo';

const chains = doc.chains.map((c) => {
  const lay = new Map(c.nodes.map((n) => [n.id, GLOBAL_LAY.get(n.id)]));
  // ★「済んでいる」判定は**全体**で行う(鎖跨ぎの依存があるため)
  const ready = c.nodes.filter((n) => isTodo(n)
    && n.deps.every((d) => !isTodo(GLOBAL_BY_ID.get(d))));
  return { c, lay, ready, todo: c.nodes.filter(isTodo) };
});

/** ★この鎖が開ける [FrdI] 項目(`unblocks` が無い鎖は `serves` を使う) */
const unblocksOf = (c) => c.unblocks || [c.serves && c.serves.item].filter(Boolean);

if (process.argv.includes('--json')) {
  console.log(JSON.stringify(
    chains.map(({ c, lay, ready, todo }) => ({
      id: c.id, title: c.title, unblocks: unblocksOf(c).length,
      todo: todo.length,
      ready: ready.map((n) => ({ id: n.id, title: n.title, size: n.size })),
      nodes: c.nodes.map((n) => ({ id: n.id, layer: lay.get(n.id), status: n.status })),
    })), null, 2));
  process.exit(0);
}

console.log(`★[FrdI] 塊を割って出た新しい葉(${doc.measuredAt} の設計)\n`);

// ★開ける項目数の多い順(= 律速の順)
chains.sort((a, b) => unblocksOf(b.c).length - unblocksOf(a.c).length);

for (const { c, lay, ready, todo } of chains) {
  const lo = Math.min(...c.nodes.map((n) => lay.get(n.id)));
  const depth = Math.max(...c.nodes.map((n) => lay.get(n.id))) + 1;
  console.log(`──────────────────────────────────────────────────────────`);
  console.log(`★チェーン ${c.id} —— ${c.title}`);
  console.log(`   奉じる項目: [${c.serves.paper}] ${c.serves.item}(物理 p.${c.serves.page})`);
  console.log(`   ★これが閉じると [FrdI] の ${unblocksOf(c).length} 件が動く: ${unblocksOf(c).join(' / ')}`);
  console.log(`   節点 ${c.nodes.length} 個・層 ${lo}–${depth - 1}(全体の層番号)・未着手 ${todo.length} 個\n`);
  for (let L = lo; L < depth; L++) {
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
  console.log(`           大きさ ${n.size} / このチェーンは ${unblocksOf(c).length} 件を開ける`);
}
console.log(`\n  合計 ${all.length} 葉。`);
console.log('  ★「葉」は依存が無いという意味であって、易しいという意味ではない。');
console.log('  ★分解そのものが我々の設計であり、実装すると変わりうる(限界は JSON の ★限界)。');
