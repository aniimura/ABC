// 2026-09-01（第 976）—— 悪い素点で局所の還元型を出す。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35BadPrimeLocal20260901 = {
  measuredAt: '2026-09-01',
  block: '第 976',
  supersedes: 'lemma35JSwap20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`jExp_variableChange`（976）——`jExp` は変数変換で不変（`j` が不変だから）。',
    '★★★★★★★★`isMinimal_baseChange_at_bad_prime`（976）——' +
      '第 908 を `C • E` に当てた形。`IsIntegral` は `IsMinimal` から' +
      '（mathlib の `instIsIntegralOfIsMinimal`）自動で出る。',
    '★★★★★★★★★★★★`hasMultiplicativeReduction_at_bad_prime`（976）——' +
      '悪い素点では完備化で**乗法還元**をもつ。' +
      '☆`0 < v_p(Δ)` は第 973 が `jExp < 0` から与える。',
  ],
  remaining:
    '★残るのは (i) 分裂性（963 の二者択一＋整モデルの `IsCharNeTwoNF` 正規化、' +
    '非分裂は 938 → 929）、(ii) `hΔ`・`hcop`（932）・`hlu`（953）・`hql`・`h2`・`h2K`・' +
    '`hvw`（961＋962）、(iii) それらを 972 に並べて `minDeltaExp_variableChange` で戻す段。',
  note:
    '★944-976 の 33 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35BadPrimeLocal20260901 を書いた');
