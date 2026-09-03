// 2026-09-01（第 994）—— 正規化は u = 1 で行える（c₄・Δ が動かない）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35NormalizeU1_20260901 = {
  measuredAt: '2026-09-01',
  block: '第 994',
  supersedes: 'lemma35SplitOrTwistBad20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★`exists_charNeTwoNF_u_one`（994）——' +
      '`2` が単元なら **`u = 1` の**変数変換で `a₁ = a₃ = 0` にできる。' +
      '☆mathlib の `toCharNeTwoNF` は `{ u := 1, r := 0, s := ⅟2·(−a₁), t := ⅟2·(−a₃) }` である。',
    '★★★★★★`variableChange_c4_Delta_of_u_one`（994）——' +
      '`u = 1` の変数変換は `c₄` と `Δ` を動かさない。',
  ],
  whyItMatters:
    '★第 993 が受ける `hc4`（`v_p(c₄) = 0`）と `hΔpos`（`0 < v_p(Δ)`）は、' +
    '正規化しても**そのまま移る**——`u = 1` だからである。' +
    '☆「正規化すると仮説が壊れるのでは」という懸念が消えた。',
  remaining:
    '☆残るのは (i) `primeSubring p` で `2` が単元であること（`p ∤ 2`、第 953 に `n = 2`）と、' +
    '正規化した曲線の `p` での整性を言う段、' +
    '(ii) 993 の 2 つの枝を 972 に流し `minDeltaExp_variableChange` で戻す段。' +
    '★`p ∣ 2` は第 991 の道（不分岐 2 次拡大＋943＋944）。',
  note:
    '★944-994 の 46 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35NormalizeU1_20260901 を書いた');
