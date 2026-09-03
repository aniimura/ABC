// 2026-09-01（第 983）—— 乗法還元なら剰余体で c₄ ≠ 0。mathlib の穴も記録。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35ResidueC4_20260901 = {
  measuredAt: '2026-09-01',
  block: '第 983',
  supersedes: 'lemma35SplitOnlyC420260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`residue_c4_ne_zero_of_multiplicativeReduction`（983）——' +
      '乗法還元なら整モデルの `c₄` は剰余体で `0` でない。' +
      '☆`v(c₄) = 1` なら `c₄ ∉ 𝔪`。★これが第 982 の唯一の実質的な仮説である。',
  ],
  mathlibGap:
    '★★★★**mathlib の穴（2026-09-01 実測）**: ' +
    '`Finite (IsLocalRing.ResidueField (p.adicCompletionIntegers L))` が' +
    '**インスタンスになっていない**（`exact?`・`infer_instance` とも失敗）。' +
    '☆数体の素点の完備化の剰余体は `𝓞 L / p` に等しく有限だが、その同一視が mathlib に無い。' +
    '★第 982 は `[Fintype k]` を要求する（`FiniteField.exists_nonsquare` と' +
    '`FiniteField.isSquare_of_char_two` を使うため）ので、' +
    'それを完備化の剰余体に当てるには**この有限性を内製する**必要がある。' +
    '☆これは第 897（`IsAdicComplete`）と同じ種類の穴である——' +
    '真だが mathlib に無いので自分で建てる。',
  remaining:
    '☆残るのは (i) 完備化の剰余体の有限性（上の mathlib の穴）、' +
    '(ii) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段、' +
    '(iii) 捻り `d` を `Lv` の整数環から `L` に降ろす段、' +
    '(iv) 972 に全部並べて `minDeltaExp_variableChange` で戻す段。',
  note:
    '★944-983 の 40 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35ResidueC4_20260901 を書いた');
