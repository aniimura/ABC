// 2026-08-31（第 897-898）—— mathlib-gap.json に「完備化の橋」の測定と埋め戻しを記録する。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.adicCompleteBridge20260831 = {
  measuredAt: '2026-08-31',
  block: '第 896-898',
  what:
    '★★★★★★★★★★**`Lemma 3.5` の「大域を各悪い素点で完備化に落とす」段**に必要な型は、' +
    '`R := v.adicCompletionIntegers K`・`Lv := v.adicCompletion K` で揃う。',
  measured: {
    present: [
      'IsDiscreteValuationRing R（mathlib: NumberField/Completion/FinitePlace.lean）',
      'IsDomain R / IsNoetherianRing R / IsLocalRing R',
      'IsFractionRing R Lv',
      'Algebra K Lv / CompleteSpace Lv',
      'IsHausdorff (maximalIdeal R) R',
    ],
    absentAtMeasurement: [
      '★★★IsPrecomplete (maximalIdeal R) R（したがって IsAdicComplete も）',
      'CompactSpace R / CompleteSpace R',
      'LocallyCompactSpace Lv / ValuativeRel Lv（IsNonarchimedeanLocalField 経由は使えない）',
    ],
    witness: 'lean/ABC3/Check/GenEll/AdicCompleteMissing.lean（第 896）',
  },
  closed: {
    at: '2026-08-31（第 897）',
    where: 'lean/ABC3/Found/GaloisRep/AdicCompleteIntegers.lean',
    statement:
      'instance isAdicComplete_adicCompletionIntegers : ' +
      'IsAdicComplete (IsLocalRing.maximalIdeal (v.adicCompletionIntegers K)) ' +
      '(v.adicCompletionIntegers K)',
    how:
      '3 段。(1) 整数環は開部分群なので閉集合、(2) 完備空間の閉集合は完備、' +
      '(3) 一次元子 ϖ に対し 𝔪^n = {y | v(y) ≤ v(ϖ)^n}（mathlib の ' +
      'maximalIdeal_pow_eq_setOf_le_v_algebraMap_pow）で位相が 𝔪 進位相であることを見る' +
      '（exists_pow_lt₀ が逆向きを与える）。あとは mathlib の IsAdic.isAdicComplete_iff。',
    nonvacuity:
      'lean/ABC3/Check/GenEll/CompletionBridgeWitness.lean（第 898）——' +
      'minDeltaExp_eq_mul_of_tateParamR（第 892）が数体の素点の完備化で実際に型付けできることを示した。',
  },
  remainingForLemma35: [
    '☆E ⊗ Lv が分裂乗法還元をもつこと（半安定 ＋ v_p(j) < 0 から）',
    '☆Tate 一意化 Φ の存在（tatePhiAddEquivAll の仮説 hloc・hI・hquad・hvalring・hvR・hvI）',
    '☆H の像が ⟨Φ(ζ)⟩ になること（Lemma 3.2, (i) の対偶——l が局所高さと互いに素）',
    '★★これ以外の数学（E_q/μ_l = E_{q^l}、j による母数の決定、Δ_min の関係、hDX）は' +
      'すべて証明済みである（第 853-895）。',
  ],
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: adicCompleteBridge20260831 を書いた');
