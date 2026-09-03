// 2026-09-01（第 989）—— 完備化の剰余体は有限：mathlib の穴を埋めた。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35ResidueFieldFinite20260901 = {
  measuredAt: '2026-09-01',
  block: '第 989',
  supersedes: 'lemma35IntegerCongr20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`finite_residueField_adicCompletionIntegers`（989）——' +
      '**`Finite (IsLocalRing.ResidueField (p.adicCompletionIntegers L))` を内製した**。' +
      '☆第 983・985 で mathlib に無いと実測した穴が埋まった。',
    '★`valued_algebraMap_eq`（989）——`Valued.v` と `p.valuation` の橋（`L` の元について）。',
    '★`mem_max_of_valued_lt_one`（989）——付値が `< 1` なら `𝔪` の元。',
    '★★★★★★★★★★★★`exists_integer_congr_completion`（989）——' +
      '完備化の整数環の元は `𝓞 L` の元と `𝔪` を法として合同。',
  ],
  howItWasDone: [
    '☆道は 3 段（第 985・987 で立てた道筋のとおり）:',
    '(1) `L` は `Lv` で稠密（`denseRange_algebraMap`、mathlib）',
    '(2) `v_p(y) ≤ 1` なら `y ≡ a (mod p)` なる `a ∈ 𝓞 L`（第 988）',
    '(3) よって `𝓞 L / p → ResidueField R` が全射、' +
      '`Ideal.finiteQuotientOfFreeOfNeBot`（mathlib）で `Finite`',
    '★第 897（`IsAdicComplete` の内製）と同じ形の作業だった。',
  ],
  plumbing:
    '★配管（3 度目の同じ穴）: `Valued.mem_nhds` は `Valued.v.restrict` と ' +
    '`ValueGroup₀` の言葉で述べられている。`γ := 1` を渡し、' +
    '`Valuation.restrict_lt_iff_lt_embedding` ＋ `simp [Units.val_one, map_one]` で落とす。' +
    '☆`tools/lean-idioms.md` に追記した。',
  remaining:
    '☆残るのは (i) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段、' +
    '(ii) 捻り `d` を `Lv` の整数環から `L` に降ろす段（989 の全射性が効く）、' +
    '(iii) 982 の二者択一を `hasSplitMultiplicativeReduction_baseChange` に流す段、' +
    '(iv) 972 に並べて `minDeltaExp_variableChange` で戻す段。',
  note:
    '★944-989 の 42 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35ResidueFieldFinite20260901 を書いた');
