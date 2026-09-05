#!/usr/bin/env node
/**
 * 前線 —— 「次にどのノードをやると一番効くか」を依存グラフから計算する。
 *
 * ★この道具の位置づけ
 *   `tools/graph.mjs` は「グラフがどうなっているか」を答える。
 *   こちらは **「だから次に何をすべきか」** を答える。
 *   Orchestrator（本体セッション）が sub-agent に仕事を配る前に、必ずここを引く。
 *
 * ★何を測るか（`ResearchPaper/orchestration.md` §2）
 *   sorry を持つノード v について:
 *     downstream  … v を推移的に import しているノード数（v が解けると解放される数）
 *     dsItems     … その下流が持つ `.src` 項目の総数（＝原典の主張いくつ分か）
 *     blockers    … v が推移的に import している **他の** sorry ノード
 *     startable   … blockers が空。すなわち **今すぐ着手できる**
 *
 *   ★`startable` でない節点に人手（agent）を割いても、上流の sorry に当たって
 *     止まる。だから配る順序は「startable のうち downstream が大きいもの」から。
 *
 * ★同時に起動する agent は **5 個まで**（`ResearchPaper/orchestration.md` §0）。
 *   だから既定の出力も 5 件で切る。足りなければ**次の波**にする。
 *
 * 使い方:
 *   node tools/frontier.mjs                 # 着手可能なものを効果の大きい順に（上位 5 件）
 *   node tools/frontier.mjs --all           # 着手不可のものも出す（blockers 付き）
 *   node tools/frontier.mjs --owner pGC     # 所属で絞る
 *   node tools/frontier.mjs --limit 0       # 件数の上限を外す（俯瞰したいときだけ）
 *   node tools/frontier.mjs --json          # Orchestrator が食う形
 */

import { execFileSync } from 'node:child_process';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));

const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

/** `graph.mjs --json` を唯一の真実として使う（ここで木を再走査しない）。 */
function loadGraph() {
  const out = execFileSync('node', [join(ROOT, 'tools', 'graph.mjs'), '--json'], {
    encoding: 'utf8', maxBuffer: 1 << 28,
  });
  return JSON.parse(out);
}

const { nodes } = loadGraph();
const byMod = new Map(nodes.map((n) => [n.mod, n]));

/** 逆辺: mod → それを直接 import しているノードの mod 一覧。 */
const rdeps = new Map();
for (const n of nodes) {
  for (const im of n.imports) {
    if (!byMod.has(im)) continue;            // Mathlib 等、木の外は辿らない
    if (!rdeps.has(im)) rdeps.set(im, []);
    rdeps.get(im).push(n.mod);
  }
}

/** v から辺を辿って到達できる集合（v 自身は含めない）。 */
function reach(start, edgesOf) {
  const seen = new Set();
  const stack = [...(edgesOf(start) ?? [])];
  while (stack.length) {
    const m = stack.pop();
    if (seen.has(m)) continue;
    seen.add(m);
    for (const w of edgesOf(m) ?? []) if (!seen.has(w)) stack.push(w);
  }
  return seen;
}

const downOf = (m) => rdeps.get(m) ?? [];
const upOf = (m) => (byMod.get(m)?.imports ?? []).filter((x) => byMod.has(x));

const sorryMods = new Set(nodes.filter((n) => n.hasSorry).map((n) => n.mod));

const rows = [];
for (const n of nodes) {
  if (!n.hasSorry) continue;
  const down = reach(n.mod, downOf);
  const up = reach(n.mod, upOf);
  const blockers = [...up].filter((m) => sorryMods.has(m)).sort();
  let dsItems = 0;
  for (const m of down) dsItems += (byMod.get(m)?.items.length ?? 0);
  rows.push({
    mod: n.mod,
    rel: n.rel,
    owner: n.owner,
    ownerKind: n.ownerKind,
    bucket: n.bucket,
    items: n.items,
    downstream: down.size,
    dsItems,
    blockers,
    startable: blockers.length === 0,
  });
}

const ownerFilter = opt('--owner');
let sel = ownerFilter ? rows.filter((r) => r.owner === ownerFilter) : rows;
if (!flag('--all')) sel = sel.filter((r) => r.startable);

sel.sort((a, b) =>
  (b.startable - a.startable) || (b.downstream - a.downstream) ||
  (b.dsItems - a.dsItems) || a.rel.localeCompare(b.rel));

/** ★既定 5 件。1 波で配る持ち場の上限（`ResearchPaper/orchestration.md` §0）。
 *  上限を外したいときだけ `--limit 0`。 */
const shown = sel.length;
const limit = Number(opt('--limit') ?? 5);
if (limit > 0) sel = sel.slice(0, limit);

if (flag('--json')) {
  console.log(JSON.stringify({ generated: new Date().toISOString(), frontier: sel }, null, 1));
  process.exit(0);
}

const total = rows.length;
const startable = rows.filter((r) => r.startable).length;
console.log('★前線（sorry を持つノード）');
console.log(`  sorry ノード ${total} / うち着手可能 ${startable}`);
if (ownerFilter) console.log(`  所属で絞り込み: ${ownerFilter}`);
console.log(`  ${flag('--all') ? '★--all: 着手不可も表示' : '（着手不可を隠している。--all で全部）'}`);
console.log();
console.log('  下流  項目  ノード');
for (const r of sel) {
  const mark = r.startable ? ' ' : '×';
  console.log(`${mark} ${String(r.downstream).padStart(5)} ${String(r.dsItems).padStart(5)}  ${r.rel}  (${r.owner})`);
  if (r.items.length) console.log(`                  項目: ${r.items.slice(0, 4).join(' / ')}${r.items.length > 4 ? ' …' : ''}`);
  if (!r.startable) {
    console.log(`                  ★上流の sorry ${r.blockers.length} 件で止まる:`);
    for (const b of r.blockers.slice(0, 3)) console.log(`                    ${byMod.get(b)?.rel ?? b}`);
    if (r.blockers.length > 3) console.log(`                    … 他 ${r.blockers.length - 3} 件`);
  }
}
console.log();
if (limit > 0 && shown > limit) {
  console.log(`  … 他 ${shown - limit} 件（★既定は 5 件まで。全部見るなら --limit 0）`);
  console.log();
}
console.log('☆次の一手は「下流」が大きい startable なノード。');
console.log('  そのノードの持ち場を sub-agent に渡すには node tools/brief.mjs --node <rel>');
console.log('★同時に起動する agent は 5 個まで。足りなければ次の波にする。');
console.log('  前線が 5 件も無いときに agent を増やしても、上流の sorry に当たって止まる。');
