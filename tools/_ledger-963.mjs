// 2026-09-01（第 963）—— (a)：分裂するか、捻れば分裂するか。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SplitDichotomy20260901 = {
  measuredAt: '2026-09-01',
  block: '第 963',
  supersedes: 'lemma35VeluDelta20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`splits_or_exists_twist_splits`（963）——' +
      '剰余体で 2 次式が分裂するか、ある捻りで分裂するかの二者択一。' +
      '☆場合分けは `∃ x, c₄x² = K` かどうかだけ:' +
      ' あれば `a₁ = 0`（`IsCharNeTwoNF`）なので `c₄X² − K` が根 `x` を持ち分裂（第 937）、' +
      ' なければ非平方元を取り（`FiniteField.exists_nonsquare`）第 938 を当てる。',
  ],
  whyItMatters:
    '★mathlib の `HasSplitMultiplicativeReduction` は剰余体で 2 次式が分裂することを要求する。' +
    '☆本定理はその要求が**常に満たせる**ことを示す——分裂しなければ捻ればよい。' +
    '非分裂の側は `minDeltaExp_eq_mul_of_nonsplit`（929）が受ける。',
  remainingLayers: [
    '☆(D3) の残りは **(b1) `hp`（完備化の付値の橋）だけ**である。' +
      '★mathlib の `valuedAdicCompletion_eq_valuation` は `Valued.v` の言葉、' +
      '本プロジェクトの `hp` は `HeightOneSpectrum.valuation Lv 𝔪_R` の言葉。' +
      '☆道筋の測定: 両者は同じ単位球 `R` を持つ全射な ℤᵐ⁰-付値なので等しい。' +
      '`Valuation.isEquiv_of_val_le_one` で同値までは出るが、' +
      '「同値かつ両方全射 ⇒ 等しい」を ℤᵐ⁰ で言う補題が mathlib に見当たらない。' +
      '具体的には「一意化元 π で `w₁ π = w₂ π = ofAdd(-1)`、' +
      '`x = u·π^n`（`u` は `R` の単元）」で押すことになる。',
  ],
  note:
    '★944-963 の 20 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SplitDichotomy20260901 を書いた');
