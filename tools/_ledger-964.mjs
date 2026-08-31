// 2026-09-01（第 964）—— (b1) が閉じた：完備化の付値の橋。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35ValuationBridge20260901 = {
  measuredAt: '2026-09-01',
  block: '第 964',
  supersedes: 'lemma35SplitDichotomy20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`le_exp_neg_one_of_lt_one` / `exp_neg_one_zpow`（964）——`ℤᵐ⁰` の小道具。',
    '★★★★★★★★★★★★`eq_of_isEquiv_of_surjective`（964）——' +
      '同値で両方全射な `ℤᵐ⁰`-付値は等しい。' +
      '☆mathlib に見当たらなかったので内製した。' +
      '一意化元 `π`（`v π = exp(-1)`）と `ρ`（`w ρ = exp(-1)`）を取り、' +
      '`w π = exp(-1)` を挿んでから任意の `x` を `π^{-n}` と比べる。',
    '★★★★★★★★★★★★★★★★`valuation_maximalIdeal_adicCompletion_eq`（964）——' +
      '`𝔪_R` の付値と `Valued.v` は一致する。' +
      '両方の単位球が `adicCompletionIntegers`、両方全射だから。',
    '★★★★★★★★★★★★★★★★★★★★`valuation_algebraMap_adicCompletion`（964）——' +
      '**局所の定理群の仮説 `hp` そのもの**。',
  ],
  allLocalDataDone:
    '★★★★(D3) が要求していた 3 つの局所データはすべて揃った: ' +
    '(a) 分裂性 → 963、(b1) 付値の橋 → 964、(d) Vélu の係数 → 956-962。' +
    '☆`hlu` は 953、`hu` は 951、極小モデルは 954 が与える。',
  remainingLayers: [
    '☆残るのは `isMuAtBadPrimes_of_veluQuotient` の本体——' +
      '各悪い素点で 954/953/951/963/964/961/962 を並べ、' +
      '950 →(955)→ 904 → 903 と流す長い呼び出しだけである。',
  ],
  note:
    '★944-964 の 21 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35ValuationBridge20260901 を書いた');
