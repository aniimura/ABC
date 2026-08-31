// 2026-09-01（第 965）—— 各悪い素点での終点が捩れ点と j で受ける形になった。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35LocalTerminal20260901 = {
  measuredAt: '2026-09-01',
  block: '第 965',
  supersedes: 'lemma35ValuationBridge20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`minDeltaExp_eq_mul_of_torsion`（965）——' +
      '第 904 は `ζ`・`uζ`・`hu`・`v`・`w` を受け、さらに `W′` が Vélu の商に' +
      '**等しい**ことを要求していた。☆本定理はそれを' +
      ' (i) 位数 `l` の点 `P`（`P ≠ 0`）と `l ∤ v(q)`、' +
      ' (ii) `j` の一致（950 が与える）、' +
      ' (iii) Vélu の帳簿 `hvw`（961・962 が与える）' +
      ' に置き換える。★これが各悪い素点での**終点**である。',
  ],
  remainingLayers: [
    '☆残るのは `isMuAtBadPrimes_of_veluQuotient` の本体で、' +
      '各悪い素点について 965 の仮説を供給する段だけである:' +
      ' `h`・`h′`（963＋909＋929）、`hp`（964）、`C`・`hC`・`hc4`（954）、' +
      ' `hΔ`・`hcop`（932）、`hlu`（953）、`hql`・`h2`・`h2K`、' +
      ' `hvw`（961＋962）、`P`（`Q ⊗ Lv`）、`hW′j`（950）。',
  ],
  note:
    '★944-965 の 22 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35LocalTerminal20260901 を書いた');
