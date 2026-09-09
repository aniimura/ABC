#!/usr/bin/env node
// 最終目標 abcConjecture_holds から見て、Skeleton が「過不足無く」繋がっているかを測る。
//   node tools/goal-chain.mjs             要約
//   node tools/goal-chain.mjs --run       Lean 計測器を回してから要約（★lake を呼ぶ）
//   --unreached / --orphans / --sorries / --vacuity   各一覧を全部出す
import { readFileSync, existsSync } from 'node:fs';
import { execFileSync } from 'node:child_process';

const args = process.argv.slice(2);
const has = (f) => args.includes(f);
if (has('--run')) {
  execFileSync('node', ['tools/leanfile.mjs', 'lean/ABC3/Check/GoalChain.lean'], { stdio: 'inherit' });
}
const path = ['lean/.cache/goal-chain.json', '.cache/goal-chain.json'].find(existsSync);
if (!path) { console.error('goal-chain.json が無い。`node tools/goal-chain.mjs --run`。'); process.exit(1); }
const d = JSON.parse(readFileSync(path, 'utf8'));

const paperOf = (n) => (n.match(/^ABC3\.Skeleton\.([^.]+)\./) || [, '?'])[1];
const tally = (list) => {
  const m = new Map();
  for (const n of list) m.set(paperOf(n), (m.get(paperOf(n)) || 0) + 1);
  return m;
};
const reachedSkel = d.skeletonAll.length - d.skeletonUnreached.length;
const emptyIf = d.interfaces.filter((x) => x.inhabitants === 0);
const noClosed = d.interfaces.filter((x) => x.inhabitants > 0 && x.closed === 0);

if (has('--unreached')) { d.skeletonUnreached.forEach((n) => console.log(n)); process.exit(0); }
if (has('--orphans')) { d.skeletonOrphans.forEach((n) => console.log(n)); process.exit(0); }
if (has('--sorries')) { d.sorryLeaves.forEach((n) => console.log(n)); process.exit(0); }
const isTrivialWitness = (n) => n.startsWith('ABC3.Check.');
if (has('--vacuity')) {
  for (const x of d.interfaces)
    console.log(`${String(x.inhabitants).padStart(3)} ${String(x.closed).padStart(3)}  ${x.name}`
      + (x.closedNames?.length ? `
        住人: ${x.closedNames.join(', ')}` : ''));
  process.exit(0);
}

console.log(`目標: ${d.goal}    （ABC3 の宣言 ${d.abc3Constants} 件を走査）`);
console.log('');
console.log('## ★不足 —— 目標が現在載っている sorry（ここが切れ目）');
if (d.sorryLeaves.length === 0) console.log('  （無し）');
for (const n of d.sorryLeaves) console.log(`  ${n}`);
console.log('');
console.log(`## ★到達 —— 目標から型で辿れる Skeleton の主張: ${reachedSkel} / ${d.skeletonAll.length}`);
const un = tally(d.skeletonUnreached), or = tally(d.skeletonOrphans), all = tally(d.skeletonAll);
console.log('');
console.log('## ★過 —— 2 通りに測る');
console.log('  「未到達」= 目標から辿れない（★切れ目の向こう側も含むので「不要」ではない）');
console.log('  「孤児」  = ABC3 のどの宣言からも参照されていない（★切れ目に依存しない指標）');
console.log('');
console.log('  論文            全体   未到達   孤児');
for (const [p, tot] of [...all.entries()].sort((a, b) => b[1] - a[1]))
  console.log(`  ${p.padEnd(14)} ${String(tot).padStart(4)}   ${String(un.get(p) || 0).padStart(4)}   ${String(or.get(p) || 0).padStart(4)}`);
console.log(`  ${'合計'.padEnd(13)} ${String(d.skeletonAll.length).padStart(4)}   ${String(d.skeletonUnreached.length).padStart(4)}   ${String(d.skeletonOrphans.length).padStart(4)}`);
console.log('');
console.log('## ★空虚の疑い —— Interface 構造体の住人');
console.log(`  構造体 ${d.interfaces.length} 件 / 住人 0 が ${emptyIf.length} 件 / 住人はあるが閉じた項が 0 のものが ${noClosed.length} 件`);
const trivialOnly = d.interfaces.filter((x) => x.closed > 0 && (x.closedNames || []).every(isTrivialWitness));
console.log(`  ★閉じた住人が Check/（自明・退化の反例置き場）にしか無いもの: ${trivialOnly.length} 件`);
for (const x of emptyIf) console.log(`  住人 0        ${x.name}`);
for (const x of noClosed) console.log(`  閉じた項 0    ${x.name}  (条件付き ${x.inhabitants})`);
for (const x of trivialOnly) console.log(`  Check のみ    ${x.name}  ← ${x.closedNames.join(', ')}`);
