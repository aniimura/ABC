// 2026-09-01（第 957）—— 訂正：veluV2 は反転で不変ではない。対ごとの偶性で受ける。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35PairEven20260901 = {
  measuredAt: '2026-09-01',
  block: '第 957',
  supersedes: 'lemma35SymmSum20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  correction:
    '★★第 956 の要約で「Vélu の和は反転で対称」と書いたが、**不正確だった**。' +
    'Tate 曲線は `a₁ = 1` なので `veluV2 = 3x² + a₄ − y` であり、`y ↦ −y−x` で変わる。' +
    '☆反転で不変なのは `veluU = (2y+x)²` の方だけである。' +
    '★正しくは「対の和が偶」——' +
    '`g(i) + g(l−i) = 2·veluU_i + 2·x_i·(veluV2_i + veluV2_{l−i})`。',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`two_mul_sum_eq_of_pair_even`（957）——' +
      '対ごとに `g i + g (l−i) = 2 * c i` なら `∑ g = 2 * ∑ c`。' +
      '☆証人 `c` を関数として受けるので選択公理も要らない。',
    '★★★★★★★★★★★★`exists_two_mul_of_pair_even`（957）——' +
      'これが Vélu の `hw : 2 * w = ∑ (…)` を満たす `w` の作り方である。',
  ],
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d) `tateXpair (ζ^{l−i}) … = tateXpair (ζ^i) …`（反転で x は変わらない）と' +
      '`veluU` の反転不変性を 957 に繋ぐ段、および `veluCurve` の楕円性。',
  ],
  note:
    '★944-957 の 14 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35PairEven20260901 を書いた');
