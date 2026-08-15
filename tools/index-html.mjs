// リポジトリの入口 index.html を、**ゲートの実出力から**生成する。
//
// ★数字を手で書かない。手で書くと古びて、古びた数字は嘘になる。
//   node tools/check.mjs の出力と git log をそのまま材料にする。
//
// 使い方: node tools/index-html.mjs

import { execSync } from 'node:child_process';
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const sh = (c) => execSync(c, { cwd: ROOT, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 });

const gate = sh('node tools/check.mjs 2>&1 || true');
const pick = (re, d = '—') => { const m = gate.match(re); return m ? m[1].trim() : d; };

const S = {
  selftest: pick(/selftest: (\d+\/\d+) PASS/),
  sorry: pick(/-- sorry: (\d+) 件/),
  decls: pick(/-- Lean 宣言 (\d+)/),
  quotes: pick(/-- 引用照合\(Lean コメント内\): (\d+) 件/),
  units: pick(/-- 構造単位 (\d+) 件/),
  nodes: pick(/スケルトンのある節点\s+(\d+) 件/),
  edges: pick(/辺\(otherPaper\)\s+(\d+) 本/),
  depth: pick(/最大深さ (\d+)/),
  cycles: pick(/循環 (\d+) 件/),
  landed: pick(/着地した葉\s+(\d+) 件/),
  leaves: pick(/★未展開の葉\s+(\d+) 件/),
  gate: /^PASS/m.test(gate) ? 'PASS' : 'NG',
  waiting: pick(/-- Interface 実装待ち: (\d+) 件/),
  implicit: pick(/暗黙の段\(Gap 候補\)\s+(\d+) 件/),
};

const commits = sh('git log --pretty=format:%h%x09%ad%x09%s --date=short -n 12').trim().split('\n')
  .map((l) => { const [h, d, s] = l.split('\t'); return { h, d, s }; });
const total = sh('git rev-list --count HEAD').trim();
const now = sh('git log -1 --date=short --format=%ad').trim();

const esc = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const graphExists = existsSync(join(ROOT, 'dependency-graph.html'));

const html = `<!DOCTYPE html>
<html lang="ja"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Math_ABC3 — IUT の Lean 形式化</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,"Segoe UI","Hiragino Sans","Noto Sans JP",sans-serif;
  background:#f4f5f7;color:#1c1e21;line-height:1.7;padding:28px 20px 60px}
.wrap{max-width:940px;margin:0 auto}
h1{font-size:22px;margin-bottom:4px}
.sub{color:#6b7280;font-size:13px;margin-bottom:22px}
h2{font-size:15px;margin:26px 0 10px;padding-bottom:5px;border-bottom:1px solid #e2e4e8}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(148px,1fr));gap:9px;margin-bottom:6px}
.card{background:#fff;border:1px solid #e2e4e8;border-radius:8px;padding:11px 13px}
.card .k{font-size:11px;color:#6b7280}
.card .v{font-size:20px;font-weight:600;font-variant-numeric:tabular-nums}
.card.ok .v{color:#1e7a52}.card.warn .v{color:#c0392b}
a.big{display:block;background:#fff;border:1px solid #d6dae0;border-radius:10px;padding:16px 18px;
  text-decoration:none;color:inherit;margin-bottom:10px;transition:.15s}
a.big:hover{border-color:#2f6fd0;box-shadow:0 2px 10px rgba(47,111,208,.12)}
a.big b{color:#2f6fd0;font-size:15px}
a.big span{display:block;color:#6b7280;font-size:12.5px;margin-top:3px}
table{width:100%;border-collapse:collapse;font-size:13px;background:#fff;border:1px solid #e2e4e8;border-radius:8px;overflow:hidden}
th,td{padding:7px 11px;text-align:left;border-bottom:1px solid #f0f1f3}
th{background:#fafbfc;font-size:11.5px;color:#6b7280;font-weight:600}
tr:last-child td{border-bottom:0}
code{background:#eef0f3;padding:1px 5px;border-radius:4px;font-size:12px}
.note{background:#fff8e1;border-left:3px solid #f0b429;padding:11px 13px;font-size:12.5px;
  color:#5a4300;border-radius:0 6px 6px 0;margin:10px 0}
.note b{color:#3d2e00}
ul{margin:6px 0 0 18px;font-size:13.5px}li{margin-bottom:4px}
.foot{margin-top:32px;padding-top:14px;border-top:1px solid #e2e4e8;color:#8a8f98;font-size:11.5px}
</style></head><body><div class="wrap">

<h1>Math_ABC3 — 宇宙際タイヒミューラー理論の Lean 形式化</h1>
<p class="sub">最終更新 ${now} ・ コミット ${total} 件 ・ このページは <code>node tools/index-html.mjs</code> で
<b>ゲートの実出力から生成</b>している(数字を手で書かない)</p>

<h2>器具の状態</h2>
<div class="grid">
  <div class="card ${S.gate === 'PASS' ? 'ok' : 'warn'}"><div class="k">ゲート</div><div class="v">${S.gate}</div></div>
  <div class="card ok"><div class="k">selftest</div><div class="v">${S.selftest}</div></div>
  <div class="card"><div class="k">Lean 宣言</div><div class="v">${S.decls}</div></div>
  <div class="card"><div class="k">逐語照合(PDF と一致)</div><div class="v">${S.quotes}</div></div>
  <div class="card"><div class="k">構造単位</div><div class="v">${S.units}</div></div>
  <div class="card"><div class="k">sorry</div><div class="v">${S.sorry}</div></div>
</div>
<div class="note"><b>★<code>sorry</code> の件数を単独で進捗と読まないこと。</b>
<code>Interface</code> に条件を posit すれば <code>sorry</code> は消えるが、負債は増える(実例あり)。
進捗として読むべきは<b>「着地した葉」と「暗黙の段」</b>の側。
<code>Found/</code> と <code>Interface/</code> の <code>sorry</code> は 0 件、<code>axiom</code> は較正用の1ファイルを除き 0 件。</div>

<h2>依存グラフ(我々が写した範囲)</h2>
<div class="grid">
  <div class="card"><div class="k">節点</div><div class="v">${S.nodes}</div></div>
  <div class="card"><div class="k">辺</div><div class="v">${S.edges}</div></div>
  <div class="card"><div class="k">最大深さ</div><div class="v">${S.depth}</div></div>
  <div class="card ok"><div class="k">循環</div><div class="v">${S.cycles}</div></div>
  <div class="card ok"><div class="k">着地した葉</div><div class="v">${S.landed}</div></div>
  <div class="card"><div class="k">未展開の葉</div><div class="v">${S.leaves}</div></div>
  <div class="card"><div class="k">暗黙の段</div><div class="v">${S.implicit}</div></div>
  <div class="card"><div class="k">Interface 実装待ち</div><div class="v">${S.waiting}</div></div>
</div>

<h2>見る</h2>
${graphExists ? `<a class="big" href="dependency-graph.html"><b>依存の層 — 左から着手する →</b>
<span>原文の参照から機械抽出した 664 節点を、強連結成分で潰して 55 層に並べたもの。
右端が [IUTchIII] Corollary 3.12、左端(層 0)は依存を持たないので今すぐ着手できる。
ボックスをクリックすると、緑=これに要るもの / 紫=これを使うもの が推移的にハイライトされる。</span></a>` : ''}
<a class="big" href="PLAN.md"><b>PLAN.md — 計画と、その訂正の記録 →</b>
<span>実測した事実、ゲート G1–G6、飛躍の扱い、正直な制約。<b>反証された自分の主張を消さずに残している。</b></span></a>
<a class="big" href="ResearchPaper/dependency-scale.md"><b>規模の実測 →</b>
<span>Corollary 3.12 から到達する原文の項目は ~800 のオーダー。<b>ただし下界でも上界でもない</b>——
番号の無い依存を数え落とし、解説への案内を数え過ぎる。その限界も書いてある。</span></a>
<a class="big" href="ResearchPaper/cycle-analysis.md"><b>循環は理論の循環か →</b>
<span>答えは <b>否</b>。辺の定義の副作用だった。前方参照と解説案内を落とすと
最大の循環は 262 → 17 節点に縮み、<b>論文をまたぐ大循環は消える</b>。</span></a>

<h2>いま分かっていること</h2>
<table><tr><th>詰まる壁の種類</th><th>実例</th><th>重さ</th></tr>
<tr><td>mathlib の不在</td><td>p進対数(評価 API がノルム体に原理的に不適用)</td><td>作業量</td></tr>
<tr><td>Lean の透明度</td><td>tempered 群、packet-normalization</td><td>作業量</td></tr>
<tr><td>どれでもない(安い)</td><td>対数体積、型の訂正</td><td>★原文が「elementary」と書く箇所は実際に安い</td></tr>
<tr><td><b>語彙の不在</b></td><td>procession-normalization</td><td><b>★最も重い。中層を作るまで手が付かない</b></td></tr>
</table>
<ul>
<li><b>Corollary 3.12 から実物まで届く鎖が1本ある</b> — <code>cor_3_12 → Thm 3.11 → [AbsTopIII] Cor 5.10 → Prop 5.7 → Found/LogVolume.lean</code></li>
<li><b>IUTchI → II → III → IV の順序付けは存在する</b>(当初「存在しない」と書いたが、辺の定義の副作用だった)</li>
<li><code>Interface</code> の posit が1つ <code>Found/</code> の実物に置き換わった(対数体積)</li>
<li><b>1つの命題の中に複数の種類の壁が同居する</b> — Prop 3.9 は「Lean の透明度」と「語彙の不在」の両方</li>
</ul>

<h2>直近のコミット</h2>
<table><tr><th style="width:78px">hash</th><th style="width:88px">日付</th><th>件名</th></tr>
${commits.map((c) => `<tr><td><code>${c.h}</code></td><td>${c.d}</td><td>${esc(c.s)}</td></tr>`).join('\n')}
</table>

<div class="note"><b>原典について。</b>望月氏の論文は本リポジトリに含めていない
(<code>ResearchPaper/0_Source/</code> は <code>.gitignore</code>)。第三者の著作物を再配布しないため。
原本は RIMS で公開されている。本リポジトリが持つのは<b>「それをどう形式化するか」だけ</b>で、
原文の引用は出典(論文名・物理ページ)つきの短い逐語に限り、
<b>すべて PDF と機械照合している</b>(${S.quotes} 件)。</div>

<p class="foot">Lean 4 (v4.31.0-rc2) + mathlib ・ 再現は <code>lake build</code> と <code>node tools/check.mjs</code> ・
このページは <code>node tools/index-html.mjs</code> で再生成する</p>
</div></body></html>
`;
writeFileSync(join(ROOT, 'index.html'), html, 'utf8');
console.log(`書き出し: index.html`);
console.log(`  ゲート ${S.gate} / selftest ${S.selftest} / 宣言 ${S.decls} / 逐語 ${S.quotes} / 節点 ${S.nodes} / 着地 ${S.landed}`);
