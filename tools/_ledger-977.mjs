// 2026-09-01（第 977）—— 第 972 の残りの小さな仮説。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SmallHyps20260901 = {
  measuredAt: '2026-09-01',
  block: '第 977',
  supersedes: 'lemma35BadPrimeLocal20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`two_ne_zero_adicCompletion` / `two_ne_zero_adicCompletionIntegers`（977）——' +
      '標数 0（第 897 の `charZero_*`）から。',
    '★`pow_mem_of_mem_ideal`（977）——`q ∈ I` なら `q^l ∈ I`。',
    '★★★★★★★★`tateModel_map_Delta_ne_zero`（977）——' +
      'Tate モデルを `K` に上げた曲線の `Δ` は `0` でない。' +
      '☆`tateModel_baseChange`（944）により Tate モデルは `W` の変数変換だから。',
  ],
  supplied: [
    '☆第 972 の仮説のうち、これで供給済みになったもの:',
    '`h2`・`h2K`・`hql`・`hΔ`（977）',
    '`C`・`hC`・`hc4ne`・`hc4`（954）',
    '極小性・乗法還元（976、908＋909＋973 経由）',
    '`hp`（964）、`hlu`（953）、`hvw`（961＋962）',
    '`P`・`hP`・`hP0`・`hW′j`・`hellQ`（970＋947＋971、972 が内製）',
  ],
  remaining:
    '★残るのは (i) 分裂性（963 ＋ 整モデルの `IsCharNeTwoNF` 正規化、非分裂は 938 → 929）、' +
    '(ii) `hcop`（932 経由で `¬ l ∣ jExp p E` から）、' +
    '(iii) それらを 972 に並べて `minDeltaExp_variableChange` で `E`・`E′` に戻す段。',
  note:
    '★944-977 の 34 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SmallHyps20260901 を書いた');
