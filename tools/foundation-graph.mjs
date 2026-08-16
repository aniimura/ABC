// ★**基礎理論**を節点に持つグラフを書き出す(自己完結 HTML、外部ライブラリ無し)。
//
// ★この道具が答える問い: **「北極星に着くには、どの理論を何本作る必要があるのか」**
//
// ★既存のグラフとの違い(ここが要点):
//   `tools/full-graph.mjs` / `graph-layers.mjs` / `graph-html.mjs` の節点は
//   **論文の番号付き項目**である。だから「その項目を書くのに何の理論が要るか」は**出ない**。
//   実際 [GenEll] の 24 件は、グラフ上は 24 個の節点だが、
//   その下には **Arakelov 理論**と **l 捩れ Galois 表現**という
//   「mathlib に 1 件も無い古典理論 2 本」が横たわっている。
//   ★**その層を節点にしたのが本図である。**
//
// ★入力:
//   ResearchPaper/foundations.json    — 基礎理論の登記(status・探索範囲・公開プロジェクト)
//   ResearchPaper/genell-needed.json  — [GenEll] の需要と完了(node tools/genell-progress.mjs --json)
//   ResearchPaper/frdi-needed.json    — [FrdI] の需要と完了(node tools/frdi-progress.mjs --json)
//
// ★限界(先に書く):
//   1. `foundations.json` は**人手の登記**である。機械が「この項目にはこの理論が要る」と
//      判定しているのではない。★したがって**抜けがありうる**——`neededBy` は下界である。
//   2. status は測定値だが**古びる**。`measuredAt` を必ず見ること。
//   3. 公開プロジェクトの「在る」は「使える」ではない。`lean-ecosystem.json` の実測では
//      drop-in で使えたものは **0 件**だった。
//
// 使い方: node tools/foundation-graph.mjs [出力先.html]

import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const OUT = process.argv[2] ?? join(ROOT, 'foundation-graph.html');

const readJson = (p) => (existsSync(p) ? JSON.parse(readFileSync(p, 'utf8')) : null);
const F = readJson(join(ROOT, 'ResearchPaper', 'foundations.json'));
if (!F) { console.error('★ResearchPaper/foundations.json が無い。'); process.exit(1); }
const genell = readJson(join(ROOT, 'ResearchPaper', 'genell-needed.json'));
const frdi = readJson(join(ROOT, 'ResearchPaper', 'frdi-needed.json'));

const esc = (s) => String(s)
  .replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');

const STATUS = {
  inMathlib: { label: 'mathlib にある', cls: 'ok' },
  partial: { label: '近いものはあるが枠組みが違う', cls: 'partial' },
  inProject: { label: '公開プロジェクトに在る/計画', cls: 'project' },
  building: { label: '我々が構築中', cls: 'building' },
  absent: { label: '見つからなかった', cls: 'absent' },
};

// ── 進捗の数(あれば)
const progress = [];
if (genell) progress.push(['[GenEll]', genell.done, genell.total]);
if (frdi) progress.push(['[FrdI]', frdi.done, frdi.total]);

// ── 基礎理論ごとに、支えている項目数を数える
const rows = F.foundations.map((f) => {
  const items = (f.neededBy ?? []).flatMap((n) => n.items.map((i) => `${n.paper} ${i}`));
  return { ...f, items };
});

const byStatus = {};
for (const r of rows) byStatus[r.status] = (byStatus[r.status] ?? 0) + 1;

const card = (r) => `
<div class="node ${STATUS[r.status]?.cls ?? 'absent'}">
  <div class="name">${esc(r.name)}</div>
  <div class="status">${esc(STATUS[r.status]?.label ?? r.status)}</div>
  <div class="why">${esc(r.why)}</div>
  <div class="items"><b>支えている項目 ${r.items.length} 件</b>: ${r.items.map(esc).join(' / ')}</div>
  ${(r.publicProjects ?? []).length ? `<div class="proj"><b>公開プロジェクト</b>${r.publicProjects.map((p) => `
    <div class="p"><a href="${esc(p.url)}">${esc(p.name)}</a> — ${esc(p.scope)}
      <div class="caveat">${esc(p.caveat)}</div></div>`).join('')}</div>` : ''}
  <details><summary>探索範囲(この「無い」の根拠)</summary>
    <ul>${(r.searched ?? []).map((s) => `<li>${esc(s)}</li>`).join('')}</ul>
    ${r.note ? `<p class="note">${esc(r.note)}</p>` : ''}
  </details>
</div>`;

const html = `<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<title>基礎理論グラフ — 北極星に着くのに何本の理論が要るか</title>
<style>
 body { font-family: Georgia, serif; max-width: 70em; margin: 2em auto; line-height: 1.6; color: #1a1a1a; padding: 0 1em; }
 h1 { font-size: 1.3em; } h2 { font-size: 1.1em; margin-top: 2em; }
 code { background: #eee; padding: 0 .2em; }
 .warn { color: #8a2f00; }
 table { border-collapse: collapse; margin: 1em 0; width: 100%; }
 td, th { border: 1px solid #ccc; padding: .35em .6em; font-size: .95em; vertical-align: top; }
 .lanes { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 1em; margin: 1.5em 0; }
 .lane h3 { font-size: .95em; margin: 0 0 .5em; padding-bottom: .3em; border-bottom: 2px solid #333; }
 .node { border: 1px solid #bbb; border-left-width: 6px; padding: .6em .8em; margin-bottom: .8em; font-size: .9em; background: #fafafa; }
 .node .name { font-weight: bold; }
 .node .status { font-size: .85em; margin: .2em 0; }
 .node .why { font-size: .85em; color: #444; }
 .node .items { font-size: .8em; color: #555; margin-top: .4em; }
 .node .proj { font-size: .82em; margin-top: .5em; background: #eef4ff; padding: .4em; }
 .node .caveat { color: #8a2f00; font-size: .95em; }
 .node details { margin-top: .4em; font-size: .8em; }
 .node .note { color: #8a2f00; }
 .absent { border-left-color: #b33; } .partial { border-left-color: #d90; }
 .project { border-left-color: #36c; } .building { border-left-color: #693; }
 .ok { border-left-color: #093; }
 .legend span { display: inline-block; margin-right: 1em; font-size: .9em; }
 .legend i { display: inline-block; width: 1em; height: .8em; margin-right: .3em; vertical-align: middle; border: 1px solid #999; }
 .legend i.absent { background: #b33; } .legend i.partial { background: #d90; }
 .legend i.project { background: #36c; } .legend i.building { background: #693; }
 .legend i.ok { background: #093; }
 .star { color: #8a2f00; font-weight: bold; }
</style>
</head>
<body>

<h1>基礎理論グラフ —— 北極星に着くのに、どの理論を何本作る必要があるのか</h1>

<p>生成: <code>node tools/foundation-graph.mjs</code>。登記は <code>ResearchPaper/foundations.json</code>(測定日 <b>${esc(F.measuredAt)}</b>)。
toolchain <code>${esc(F.ourToolchain)}</code>。</p>

<p class="warn"><span class="star">★この図が既存のグラフと違う点</span>:
<code>tools/full-graph.mjs</code> などの節点は<b>論文の番号付き項目</b>である。
だから「その項目を書くのに何の理論が要るか」は<b>出ない</b>。
実際 [GenEll] の 24 件はグラフ上 24 個の節点だが、その下には
<b>mathlib に 1 件も無い古典理論</b>が横たわっている。<b>その層を節点にしたのが本図である。</b></p>

<h2>進捗(項目の数)</h2>
<table>
<tr><th>論文</th><th>完了 / 需要</th><th>器具</th></tr>
${progress.map(([n, d, t]) => `<tr><td>${esc(n)}</td><td><b>${d} / ${t}</b></td><td><code>node tools/${n === '[GenEll]' ? 'genell' : 'frdi'}-progress.mjs</code></td></tr>`).join('')}
</table>

<h2>基礎理論の状態</h2>

<p class="legend">
${Object.entries(STATUS).map(([, v]) => `<span><i class="${v.cls}"></i>${esc(v.label)}</span>`).join('')}
</p>

<table>
<tr><th>状態</th><th>件数</th></tr>
${Object.entries(byStatus).map(([k, v]) => `<tr><td>${esc(STATUS[k]?.label ?? k)}</td><td>${v}</td></tr>`).join('')}
</table>

<div class="lanes">
  <div class="lane">
    <h3>① 見つからなかった / 枠組みが違う</h3>
    ${rows.filter((r) => r.status === 'absent' || r.status === 'partial').map(card).join('')}
  </div>
  <div class="lane">
    <h3>② 公開プロジェクトに在る/計画</h3>
    ${rows.filter((r) => r.status === 'inProject').map(card).join('')}
  </div>
  <div class="lane">
    <h3>③ 我々が構築中 / mathlib にある</h3>
    ${rows.filter((r) => r.status === 'building' || r.status === 'inMathlib').map(card).join('')}
  </div>
</div>

<h2 class="warn">★この図の限界(先に書く)</h2>
<ol>
${(F._comment ?? []).filter((c) => c.startsWith('★') || c.includes('限界')).map((c) => `<li>${esc(c)}</li>`).join('')}
<li><b>これは人手の登記である。</b>機械が「この項目にはこの理論が要る」と判定しているのではない。
    したがって <code>neededBy</code> は<b>下界</b>であり、抜けがありうる。</li>
<li><b>status は測定値だが古びる。</b><code>measuredAt</code> を必ず見ること。</li>
<li><b>公開プロジェクトの「在る」は「使える」ではない。</b>
    <code>ResearchPaper/lean-ecosystem.json</code> の実測では、drop-in で使えたものは <b>0 件</b>だった。</li>
</ol>

</body>
</html>
`;

writeFileSync(OUT, html);
console.log(`書き出し: ${OUT}`);
console.log(`★基礎理論 ${rows.length} 件 / 状態別: ${Object.entries(byStatus).map(([k, v]) => `${STATUS[k]?.label ?? k} ${v}`).join(' / ')}`);
for (const r of rows) {
  console.log(`  [${(STATUS[r.status]?.label ?? r.status).padEnd(24)}] ${r.name} — 支えている項目 ${r.items.length} 件`);
}
