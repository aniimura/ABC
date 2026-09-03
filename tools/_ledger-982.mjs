// 2026-09-01（第 982）—— 二者択一に要るのは c₄ が単元であることだけ。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SplitOnlyC420260901 = {
  measuredAt: '2026-09-01',
  block: '第 982',
  supersedes: 'lemma35IntegralNF20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`splits_or_exists_twist_splits\'\'`（982）——' +
      '`φK = 0` なら 2 次式は `c₄X²`（`a₁ = 0` だから）で `X = 0` を根にもつので分裂する。' +
      '☆したがって二者択一に要るのは **`φ c₄ ≠ 0`（＝乗法還元）と `IsCharNeTwoNF` の 2 つだけ**。',
  ],
  progression:
    '★分裂性の仮説はこの 3 ブロックで削れた: ' +
    '第 963（`ringChar k ≠ 2`・`hc`・`hK`）→ ' +
    '第 979（`ringChar` が落ちた）→ ' +
    '第 982（`hK` も落ちた）。' +
    '☆`φ c₄ ≠ 0` は `HasMultiplicativeReduction` の `multiplicativeReduction` そのもので、' +
    '第 976 が与える。`IsCharNeTwoNF` は第 981 が体の側の `a₁ = a₃ = 0` から与える。',
  remaining:
    '☆残るのは (i) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段' +
    '（`p ∤ 2` なら変数変換が `p` で整、`p ∣ 2` は別の手当て）、' +
    '(ii) 982 の二者択一を `hasSplitMultiplicativeReduction_baseChange` に流す段、' +
    '(iii) 972 に全部並べて `minDeltaExp_variableChange` で戻す段。',
  note:
    '★944-982 の 39 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SplitOnlyC420260901 を書いた');
