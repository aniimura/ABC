// 2026-09-01（第 967）—— hW′j を作る段。在庫の引き忘れも記録。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35JMapVelu20260901 = {
  measuredAt: '2026-09-01',
  block: '第 967',
  supersedes: 'lemma35PointTransport20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`j_map_velu_vcPoint`（967）——' +
      '`L` の上で与えた Vélu の商の `j` を、底変換して変数変換した先の商の `j` に繋ぐ。' +
      '☆道は 2 段: `veluQuotientFull_baseChange`（`E′ ⊗ K` は `E ⊗ K` の `⟨Q ⊗ K⟩` による商）' +
      '＋ `j_veluQuotientFull_nsmul_variableChange`（950）。' +
      '★これが `isMuAtBadPrimes_of_veluQuotient` の `hW′j` を作る段である。',
  ],
  alreadyInStock:
    '★★`addOrderOf_rhPoint`・`addOrderOf_vcPoint` は' +
    '**すでに `Found/GenEll/PointVariableChange.lean:452, 1001` にあった**。' +
    '☆第 958 に続いて 1 セッションで 2 回目の在庫引き忘れ。' +
    '`tools/lean-idioms.md` に「その回に使う名前を全部まとめて grep する」と追記した。',
  remainingLayers: [
    '☆残るのは `isMuAtBadPrimes_of_veluQuotient` の本体。' +
      '967 が `hW′j` を `C • (E ⊗ Lv)` の形で与えるので、' +
      '`tateModel_baseChange`（944）で Tate モデルの形に書き換えてから 965 に渡す。' +
      '★点は `vcPoint C (E ⊗ Lv) (rhPoint φ E Q)` で、966 が `l • P = 0`・`P ≠ 0` を与える。',
  ],
  note:
    '★944-967 の 24 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35JMapVelu20260901 を書いた');
