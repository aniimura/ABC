// 使い捨て: [FrdI] 同一論文内の強連結成分を測り、**その内部の辺を1本ずつ文脈つきで出す**。
//
// ★目的は「7節点の循環が本物か、辺の抽出の副産物か」の判定である。
//   我々は既に2回、辺の定義の副作用で誤った構造を見ている(262節点の大循環 / 深さ3)。
//   `cycle-probe.mjs` は IUTchIII を根とした到達部分しか見ないので、
//   FrdI 内部を直接見るためにこれを書いた。
//
// 使い方: node tools/_frdi-scc.mjs [<0_Source のファイル名(拡張子なし)>] [<焦点にする項目>]

import { readFileSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const NAME = process.argv[2] ?? 'The Geometry of Frobenioids I';
const FOCUS = process.argv[3] ?? 'Definition 1.3';
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';

const lines = readFileSync(join(ROOT, 'ResearchPaper', '0_Source', `${NAME}.txt`), 'utf8').split(/\r?\n/);
const pageOf = new Array(lines.length).fill(0);
{ let pg = 0; for (let i = 0; i < lines.length; i++) { const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }

// 宣言行 = 行頭 + 番号 + ピリオド + (行末 または 括弧つき表題)。intra-graph.mjs と同じ条件。
const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
const list = [];
for (let i = 0; i < lines.length; i++) { const m = lines[i].match(declRe); if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] }); }
for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;
const decls = new Map();
for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);

// 辺(同一論文内)。`[Tag], Kind N.M` に食われている箇所は除く。
const adj = new Map();      // name -> Map(to -> [{ctx, line, page}])
for (const d of decls.values()) {
  const body = lines.slice(d.line, d.end).join('\n');
  const spans = [];
  // ★TAGRE=old は既存の道具と同じ `[A-Za-z]+`(数字を許さない)。
  //   FrdI は `[Mzk15]` の形で引用するので、old だと **論文をまたぐ引用が
  //   masking されず、同一論文内の参照として二重に拾われる**。
  const TAG = process.env.TAGRE === 'old' ? '[A-Za-z]+' : '[A-Za-z][A-Za-z0-9]*';
  const xre = new RegExp(`\\[(${TAG})\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(xre)) spans.push([m.index, m.index + m[0].length]);
  const out = new Map();
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to === d.name || !decls.has(to)) continue;
    const ctx = body.slice(Math.max(0, m.index - 160), m.index + m[0].length + 120).replace(/\s+/g, ' ');
    const ln = d.line + body.slice(0, m.index).split('\n').length - 1;
    if (!out.has(to)) out.set(to, []);
    out.get(to).push({ ctx, line: ln + 1, page: pageOf[ln] });
  }
  adj.set(d.name, out);
}

// Tarjan
function scc(keep) {
  const g = new Map();
  for (const [v, outs] of adj) g.set(v, [...outs.keys()].filter((w) => adj.has(w) && keep(v, w, outs.get(w))));
  const index = new Map(), low = new Map(), onstk = new Set(), stk = [], comp = new Map();
  let idx = 0, nc = 0;
  for (const s of g.keys()) {
    if (index.has(s)) continue;
    const work = [[s, 0]];
    while (work.length) {
      const fr = work[work.length - 1], v = fr[0];
      if (fr[1] === 0) { index.set(v, idx); low.set(v, idx); idx++; stk.push(v); onstk.add(v); }
      let rec = false;
      const es = g.get(v) ?? [];
      for (let i = fr[1]; i < es.length; i++) {
        const w = es[i];
        if (!index.has(w)) { fr[1] = i + 1; work.push([w, 0]); rec = true; break; }
        else if (onstk.has(w)) low.set(v, Math.min(low.get(v), index.get(w)));
      }
      if (rec) continue;
      if (low.get(v) === index.get(v)) { const c = nc++; for (;;) { const w = stk.pop(); onstk.delete(w); comp.set(w, c); if (w === v) break; } }
      work.pop();
      if (work.length) { const p = work[work.length - 1][0]; low.set(p, Math.min(low.get(p), low.get(v))); }
    }
  }
  const size = new Map();
  for (const [, c] of comp) size.set(c, (size.get(c) ?? 0) + 1);
  return { comp, size, nc };
}

const base = scc(() => true);
const bigs = [...base.size.entries()].filter(([, n]) => n > 1).sort((a, b) => b[1] - a[1]);
console.log(`${NAME}`);
console.log(`番号付き項目 ${decls.size} 件 / SCC ${base.nc} 個 / サイズ>1 の SCC ${bigs.length} 個`);
for (const [c, n] of bigs) {
  const mem = [...base.comp].filter(([, cc]) => cc === c).map(([k]) => k)
    .sort((a, b) => decls.get(a).line - decls.get(b).line);
  console.log(`  SCC(${n}): ${mem.map((m) => `${m}(p.${decls.get(m).page})`).join('  |  ')}`);
}

// ★NODES=... を与えたら、その集合の誘導部分グラフ(出辺・入辺)を出す。
if (process.env.NODES) {
  const set = new Set(process.env.NODES.split(',').map((s) => s.trim()));
  const order = [...set].sort((a, b) => (decls.get(a)?.line ?? 0) - (decls.get(b)?.line ?? 0));
  console.log(`\n===== 誘導部分グラフ(${set.size} 節点) =====`);
  for (const v of order) {
    if (!decls.has(v)) { console.log(`  ${v}  ← 宣言が見つからない`); continue; }
    const outs = [...adj.get(v)].filter(([w]) => set.has(w)).map(([w, o]) => `${w}(${o.length})`);
    const ins = order.filter((u) => u !== v && set.has(u) && adj.get(u)?.has(v)).map((u) => `${u}(${adj.get(u).get(v).length})`);
    console.log(`  ${v} (p.${decls.get(v).page})`);
    console.log(`      出→ ${outs.length ? outs.join(', ') : '(無し)'}`);
    console.log(`      入← ${ins.length ? ins.join(', ') : '(無し)'}`);
  }
  // 位相順序(Kahn)。残ったものがあれば循環。
  const indeg = new Map(order.map((v) => [v, 0]));
  for (const v of order) for (const [w] of adj.get(v) ?? []) if (set.has(w) && w !== v) indeg.set(w, indeg.get(w) + 1);
  const q = order.filter((v) => indeg.get(v) === 0); const topo = [];
  while (q.length) { const v = q.shift(); topo.push(v); for (const [w] of adj.get(v) ?? []) { if (!set.has(w) || w === v) continue; indeg.set(w, indeg.get(w) - 1); if (indeg.get(w) === 0) q.push(w); } }
  console.log(`\n  位相順序(先に作るべき順): ${topo.join('  →  ')}`);
  console.log(`  ★循環に残った節点: ${order.length - topo.length} 件`);
}

// 焦点の SCC の内部の辺を1本ずつ
const fc = base.comp.get(FOCUS);
if (fc === undefined) { console.log(`\n★${FOCUS} は宣言として見つからない`); process.exit(0); }
const mem = new Set([...base.comp].filter(([, c]) => c === fc).map(([k]) => k));
console.log(`\n===== ${FOCUS} を含む SCC(${mem.size} 節点)の内部の辺 =====`);
let n = 0;
for (const v of [...mem].sort((a, b) => decls.get(a).line - decls.get(b).line)) {
  for (const [w, occ] of adj.get(v)) {
    if (!mem.has(w)) continue;
    n++;
    console.log(`\n[${n}] ${v} (p.${decls.get(v).page})  →  ${w} (p.${decls.get(w).page})   出現 ${occ.length} 回`);
    for (const o of occ.slice(0, 4)) console.log(`     p.${o.page} L${o.line}: …${o.ctx}…`);
  }
}
console.log(`\n内部の辺 ${n} 本`);
