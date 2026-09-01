// 2026-09-01（第 987）—— 剰余体の有限性を建てるための在庫の追加調査。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

const g = d.residueFieldFiniteGap20260901;
g.additionalStock = [
  '☆2026-09-01 の追加調査（`lean_check` で確認）:',
  '★あるもの: `IsDedekindDomain.exists_forall_sub_mem_ideal`——' +
    '有限個の素イデアルについての中国剰余型の近似（`∃ y, ∀ i, y − x i ∈ P i ^ e i`）。' +
    '☆これが「`𝓞 L` から `p` を法として代表元を取る」段の核になる。',
  '★あるもの: `IsDedekindDomain.HeightOneSpectrum.Support.finite`——' +
    '`k : K` の付値が `< 1` になる素点は有限個。' +
    '☆`y ∈ L` の分母に現れる素点が有限個であることを言うのに使う。',
  '★あるもの: `HeightOneSpectrum.valuation_lt_one_iff_dvd`（`R` の元について）',
  '✗ 無いもの: `HeightOneSpectrum.valuationSubring`／`mem_valuationSubring_iff`',
  '✗ 無いもの: `IsLocalization.AtPrime.quotientEquiv`／' +
    '`Ideal.quotientEquivQuotientLocalization`（局所化と剰余の同一視）',
];
g.refinedPlan = [
  '☆精密化した道筋:',
  '(a) `x : R` に稠密性を当てて `y : L` を取る（`v(x − y) < 1`）。',
  '(b) `y` の「分母」に現れる素点は有限個（`Support.finite`）で、いずれも `p` と異なる' +
    '（`v_p(y) ≤ 1` だから）。',
  '(c) `exists_forall_sub_mem_ideal` で、それらの素点で `y` の分母を払い、' +
    '`p` では `1` に合同な `s ∈ 𝓞 L` を取る。`a := s·y` が `𝓞 L` に入り `y ≡ a (mod p)`。',
  '(d) 以上で `𝓞 L → ResidueField R` が全射。核が `p` を含むので ' +
    '`𝓞 L / p → ResidueField R` が全射で、`Ideal.finiteQuotientOfFreeOfNeBot` から `Finite`。',
];

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: residueFieldFiniteGap に在庫調査を追記した');
