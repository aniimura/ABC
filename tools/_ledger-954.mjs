// 2026-09-01（第 954）—— (D3) の (b) の半分が閉じた：半安定性が極小モデルを直に与える。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35BadPrimeData20260901 = {
  measuredAt: '2026-09-01',
  block: '第 954',
  supersedes: 'lemma35UnitAtPrime20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`minDeltaExp_pos_of_jExp_neg`（954）——半安定で `jExp < 0` なら `Δ_min > 0`。',
    '★★★★★★★★★★★★★★★★`exists_minimal_c4_unit_of_jExp_neg`（954）——' +
      '悪い素点では `SemistableAt` の左の選択肢（`Δ_min = 0`）が落ちるので、' +
      '`C`・`hC`・`hc4ne`・`hc4` の 4 つがそのまま出る。' +
      '☆これが `minDeltaExp_eq_mul_of_veluMu`（904）の極小モデル引数の出どころである。',
  ],
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性（`hasMultiplicativeReduction_baseChange`(909) が' +
      '乗法還元まで与えるので、残るのは剰余体で 2 次式が分裂すること）、' +
      '(b1) `hp`（完備化の付値の両立）——' +
      'mathlib の `valuedAdicCompletion_eq_valuation` は `Valued.v` の言葉で述べており、' +
      '本プロジェクトの `hp` は `HeightOneSpectrum.valuation Lv 𝔪_R` の言葉である。' +
      'その 2 つを同一視する橋が要る。' +
      '(d) 各 `ζ` に対する `v`・`w` と `veluCurve` の楕円性。',
  ],
  note:
    '★944-954 の 11 ブロックで、(c) と `hu` と (b) の極小モデル側が閉じた。' +
    '`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35BadPrimeData20260901 を書いた');
