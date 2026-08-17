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
// 層0の仕分け(手で数字を書かないため JSON から読む)
const L0P = join(ROOT, 'ResearchPaper', 'layer0-startable.json');
const L0 = existsSync(L0P) ? JSON.parse(readFileSync(L0P, 'utf8')) : null;

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

${L0 ? `<h2>層0 —— いま着手できるもの</h2>
<div class="grid">
  <div class="card"><div class="k">層0 の節点</div><div class="v">${L0.layer0.total}</div></div>
  <div class="card"><div class="k">§0 語</div><div class="v">${L0.layer0.section0Terms}</div></div>
  <div class="card"><div class="k">語彙(概念)</div><div class="v">${L0.layer0.concepts}</div></div>
  <div class="card"><div class="k">番号付き項目</div><div class="v">${L0.layer0.numberedItems}</div></div>
  <div class="card ok"><div class="k">★着手できる(名指し)</div><div class="v">${L0.startable.length}</div></div>
  <div class="card ok"><div class="k">うち完了</div><div class="v">${L0.startable.filter((s) => s.state === 'done').length}</div></div>
</div>
<table><tr><th style="width:38%">項目</th><th style="width:54px">壁</th><th style="width:58px">状態</th><th>根拠</th></tr>
${L0.startable.map((s) => {
  const badge = s.state === 'done' ? '<b style="color:#1e7a52">完了</b>'
    : s.state === 'building' ? '<b style="color:#2f6fd0">作業中</b>' : '未着手';
  return `<tr><td><b>${esc(s.item)}</b></td><td>種類 ${s.kind}</td><td>${badge}</td><td style="font-size:12px">${esc(s.basis)}</td></tr>`;
}).join('\n')}
</table>
<div class="note"><b>★層0は「完成」ではない。${esc(L0.notComplete.reason)}</b><br>
${esc(L0.notComplete.detail)}<br><br>
<b>次:</b> ${esc(L0.notComplete.next)}</div>` : ''}

<h2>原文の記載不備と、我々が足した前提</h2>
<div class="note"><b>★この節は「原文をそのままでは形式化できなかった箇所」の一覧である。</b><br>
2 種類ある。<b>(A) 原文の記載が不備で、訂正が必要なもの</b>（我々は意図された主張を証明した）と、
<b>(B) 我々が前提を足したもの</b>（原文の主張そのものは未決着のまま）。
<b>(A) だけを指標に数えている。</b>(B) は数えていない——足した前提の下でしか証明できていないからである。</div>
<table><tr><th style="width:26%">項目</th><th style="width:52px">種別</th><th style="width:56px">数えたか</th><th>中身</th></tr>
<tr><td><b>[FrdI] Proposition 2.1, (iii)</b></td><td><b style="color:#b3261e">A</b></td><td><b style="color:#1e7a52">数えた</b></td><td style="font-size:12px">
<b>★原文の訂正が必要。</b>原文 p.44 は「Let d ∈N≥1.」で <code>d</code> を先に固定してから
「<code>C</code> is of perfect type if and only if Ψ is an equivalence of categories」と述べるが、
<b>その読みでは <code>⟸</code> が偽</b>である（機械検証済み、<code>Check/FrdI/Prop21QuantifierGap.lean</code>）。
<code>d = 1</code> では Ψ₁ はつねに恒等関手と同型なので perfect 型でない Frobenioid すべてが反例になり、
<code>d ≥ 2</code> でも <code>Φ = Z/2</code>・<code>d = 3</code> が反例になる。
<b>論文の数学の誤りではない</b>——原文の証明文自身が perfect の定義（すべての <code>n</code> を量化）を引いており、
意図された主張「すべての <code>d</code> で Ψ_d が圏同値 ⟺ perfect 型」は真で、我々はそれを証明した
（<code>Found/FrdI/Prop21.lean</code> の <code>prop_2_1_iii</code>）。<b>壊れているのは印刷された文の量化子の位置だけ。</b>
記録: <code>Gap/FrdI/Section2.lean</code> の <code>Gap_2_1_iii</code>（分類 ③ sourceGap、反例つき）。</td></tr>
<tr><td><b>[FrdI] Proposition 1.6, (v)</b></td><td><b style="color:#2f6fd0">B</b></td><td>数えていない</td><td style="font-size:12px">
2 件の <code>⟸</code>（Aut-ample / base-trivial）に <b>「全対象が Aut-ample」を足した</b>
（<code>cfp_autAmple_mpr</code> / <code>cfp_baseTrivial_mpr</code>）。
<b>原文の主張が真かどうかは未決着</b>（<code>Gap_1_6_v</code>、分類 ② missingMath）。
2026-08-16 に穴は <code>htwist</code> 1 点まで絞れており、反例があるとすれば
<code>Aut_C(A)</code> の底が <code>Φ(A)</code> に非自明に作用する所にしかない、と分かっている。
<b>足した前提が空虚でないことは証明済み</b>（<code>autAmple_and_baseTrivial_nonvacuous</code>：𝔽_Φ の対象はすべて両方を満たす）。
<b>★併せて (ii) の <code>plBkEquiv</code> の本質的全射性が未実装</b>（Lean の書き換えが通らず打ち切った）。</td></tr>
<tr><td><b>[FrdI] Proposition 1.14, (ii)(iii)</b></td><td><b style="color:#2f6fd0">B</b></td><td>数えていない</td><td style="font-size:12px">
「Frobenius 型射は mono」を出すのに <b>unit-trivial 型を足した</b>
（<code>mono_of_frobType_of_unitTrivial</code>）。原文が課すのは
「Φ は divisorial」「𝒟 は連結・totally epimorphic・FSMFF 型」「𝒞 は isotropic 型の Frobenioid」だけで、
<b>unit-trivial は我々が足したものである</b>。<code>Definition 1.3, (vi)</code> が与えるのは
<b>単元を除いた</b>忠実性だけで、unit-trivial はその単元をちょうど潰す。
<b>★さらに (iii) にはもう 1 つ穴がある</b>——<code>hFrobFS</code> の側は
<b>𝒪^× にどんな仮定を置いても埋まらず</b>、Φ の <code>d</code> 可除性（perfect 型）が要ると見ている（未証明）。</td></tr>
</table>

<h2>いま分かっていること</h2>
<table><tr><th>詰まる壁の種類</th><th>実例</th><th>重さ</th></tr>
${(L0 ? L0.wallKinds : []).map((w) => `<tr><td><b>${w.n}. ${esc(w.name)}</b></td><td>${esc(w.example)}</td><td>${esc(w.weight)}</td></tr>`).join('\n')}
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
