// 循環が「理論の循環」なのか「辺の定義の副作用」なのかを測る。
//
// 我々の辺は「項目の本文が Kind N.M を名指しした」である。これには
//   (a) 証明が実際に依拠している        ← 真の依存
//   (b) 前方参照("cf. …, below")        ← 依存ではない
//   (c) 解説への案内("the discussion of …") ← 依存ではない
//   (d) Remark からの言及(注釈)          ← 依存とは限らない
// が混ざっている。(b)(c)(d) を落としたときに循環が解けるなら、
// **循環は辺の定義の副作用であって理論の循環ではない。**
//
// 使い方: node tools/cycle-probe.mjs

import { readFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const FILE_OF = JSON.parse(readFileSync(join(ROOT, 'tools', '_fileof.json'), 'utf8'));

const loaded = new Map();
function load(tag) {
  if (loaded.has(tag)) return loaded.get(tag);
  const f = FILE_OF[tag]; const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) { loaded.set(tag, null); return null; }
  const lines = readFileSync(p, 'utf8').split(/\r?\n/);
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  for (let i = 0; i < lines.length; i++) { const m = lines[i].match(declRe); if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i }); }
  for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;
  const decls = new Map(); for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  const o = { decls, lines }; loaded.set(tag, o); return o;
}

/** 辺を、文脈つきで返す。 */
function edgesOf(tag, name) {
  const P = load(tag); if (!P) return null;
  const d = P.decls.get(name); if (!d) return null;
  const body = P.lines.slice(d.line, d.end).join('\n');
  const es = []; const spans = [];
  const ctx = (i, len) => body.slice(Math.max(0, i - 90), i + len + 60).replace(/\s+/g, ' ');
  const xre = new RegExp(`\\[([A-Za-z]+)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(xre)) {
    es.push({ to: `${m[1]} / ${m[2]} ${m[3]}`, ctx: ctx(m.index, m[0].length), cross: true });
    spans.push([m.index, m.index + m[0].length]);
  }
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to !== name && P.decls.has(to)) es.push({ to: `${tag} / ${to}`, ctx: ctx(m.index, m[0].length), cross: false });
  }
  // 同じ先が複数回出るときは、最も「依存らしい」文脈を残す(below/discussion が無いもの)
  const best = new Map();
  for (const e of es) {
    const cur = best.get(e.to);
    const soft = /below|the discussion of|cf\.|see |Remark/i.test(e.ctx);
    if (!cur || (cur.soft && !soft)) best.set(e.to, { ...e, soft });
  }
  return [...best.values()];
}

const ROOTK = 'IUTchIII / Corollary 3.12';
const all = new Map();          // v -> [{to, ctx, cross, soft}]
{
  const q = [ROOTK], seen = new Set();
  while (q.length) {
    const k = q.shift(); if (seen.has(k)) continue; seen.add(k);
    const [t, n] = k.split(' / ');
    const es = edgesOf(t, n) ?? [];
    all.set(k, es);
    for (const e of es) if (!seen.has(e.to)) q.push(e.to);
  }
  for (const k of seen) if (!all.has(k)) all.set(k, []);
}

function scc(filter) {
  const adj = new Map();
  for (const [v, es] of all) adj.set(v, es.filter((e) => all.has(e.to) && filter(v, e)).map((e) => e.to));
  const index = new Map(), low = new Map(), onstk = new Set(), stk = [], comp = new Map();
  let idx = 0, nc = 0;
  for (const s of adj.keys()) {
    if (index.has(s)) continue;
    const work = [[s, 0]];
    while (work.length) {
      const fr = work[work.length - 1], v = fr[0];
      if (fr[1] === 0) { index.set(v, idx); low.set(v, idx); idx++; stk.push(v); onstk.add(v); }
      let rec = false;
      const es = adj.get(v) ?? [];
      for (let i = fr[1]; i < es.length; i++) {
        const w = es[i]; if (!adj.has(w)) continue;
        if (!index.has(w)) { fr[1] = i + 1; work.push([w, 0]); rec = true; break; }
        else if (onstk.has(w)) low.set(v, Math.min(low.get(v), index.get(w)));
      }
      if (rec) continue;
      if (low.get(v) === index.get(v)) { const c = nc++; while (true) { const w = stk.pop(); onstk.delete(w); comp.set(w, c); if (w === v) break; } }
      work.pop();
      if (work.length) { const p = work[work.length - 1][0]; low.set(p, Math.min(low.get(p), low.get(v))); }
    }
  }
  const size = new Map();
  for (const [, c] of comp) size.set(c, (size.get(c) ?? 0) + 1);
  const big = [...size.values()].filter((s) => s > 1).sort((a, b) => b - a);
  return { nc, big, max: big[0] ?? 1, inBig: big.reduce((a, b) => a + b, 0), comp, size };
}

const base = scc(() => true);
console.log(`基準: 節点 ${all.size} / SCC ${base.nc} / 最大 ${base.max} / 循環に属する節点 ${base.inBig}`);

// 最大 SCC の中の辺を、文脈で分類
const bigC = [...base.size.entries()].sort((a, b) => b[1] - a[1])[0][0];
const inBig = new Set([...base.comp].filter(([, c]) => c === bigC).map(([k]) => k));
let tot = 0; const cat = { below: 0, discussion: 0, cf: 0, fromRemark: 0, toRemark: 0, plain: 0 };
const samples = [];
for (const v of inBig) for (const e of all.get(v) ?? []) {
  if (!inBig.has(e.to)) continue;
  tot++;
  const c = e.ctx;
  if (/\bbelow\b/i.test(c)) cat.below++;
  else if (/the discussion of|the discussion surrounding/i.test(c)) cat.discussion++;
  else if (/\bcf\.\s*$|\bcf\.\s*\[|\[cf\./i.test(c)) cat.cf++;
  else if (v.includes('/ Remark')) cat.fromRemark++;
  else if (e.to.includes('/ Remark')) cat.toRemark++;
  else { cat.plain++; if (samples.length < 6) samples.push(`${v}  →  ${e.to}\n      …${c.slice(60, 190)}…`); }
}
console.log(`\n最大 SCC(${inBig.size} 節点)の内部の辺 ${tot} 本の文脈:`);
for (const [k, v] of Object.entries(cat)) console.log(`  ${k.padEnd(12)} ${String(v).padStart(4)}  ${(v / tot * 100).toFixed(1)}%`);
console.log(`\n"plain"(どの目印も無い = 依存の可能性が高い)の例:`);
for (const s of samples) console.log('  ' + s);


// 残った循環の中身を出す
{
  const f=(v,e)=>!/below|the discussion of/i.test(e.ctx) && !v.includes('/ Remark') && !e.to.includes('/ Remark');
  const r=scc(f);
  const bySize=[...r.size.entries()].filter(([,n])=>n>1).sort((a,b)=>b[1]-a[1]);
  console.log('\n★注釈・前方参照・解説案内を落としたあとに残る循環:');
  for(const [c,n] of bySize.slice(0,6)){
    const mem=[...r.comp].filter(([,cc])=>cc===c).map(([k])=>k);
    console.log('  SCC(' + n + ' 節点): ' + mem.join('  |  '));
  }
  console.log('  残る循環の総数: ' + bySize.length + ' 個 / 循環に属する節点 ' + r.inBig);
}
// 辺を落としたときに循環が解けるか
const tests = [
  ['below を落とす', (v, e) => !/\bbelow\b/i.test(e.ctx)],
  ['discussion を落とす', (v, e) => !/the discussion of|the discussion surrounding/i.test(e.ctx)],
  ['Remark からの辺を落とす', (v) => !v.includes('/ Remark')],
  ['Remark への辺を落とす', (v, e) => !e.to.includes('/ Remark')],
  ['Remark を両方向で落とす', (v, e) => !v.includes('/ Remark') && !e.to.includes('/ Remark')],
  ['below + discussion を落とす', (v, e) => !/\bbelow\b|the discussion of/i.test(e.ctx)],
  ['below + discussion + Remark 両方向', (v, e) => !/\bbelow\b|the discussion of/i.test(e.ctx) && !v.includes('/ Remark') && !e.to.includes('/ Remark')],
];
console.log(`\n辺を落としたときの循環:`);
console.log(`  ${'条件'.padEnd(34)} 最大SCC  循環内の節点  SCC数`);
console.log(`  ${'(基準)'.padEnd(36)} ${String(base.max).padStart(5)} ${String(base.inBig).padStart(11)} ${String(base.nc).padStart(7)}`);
for (const [name, f] of tests) {
  const r = scc(f);
  console.log(`  ${name.padEnd(36)} ${String(r.max).padStart(5)} ${String(r.inBig).padStart(11)} ${String(r.nc).padStart(7)}`);
}
