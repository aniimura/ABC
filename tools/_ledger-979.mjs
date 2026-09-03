// 2026-09-01（第 979）—— 剰余標数 2 でも二者択一は成り立つ。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SplitCharTwo20260901 = {
  measuredAt: '2026-09-01',
  block: '第 979',
  supersedes: 'lemma35Hcop20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`exists_sq_of_charTwo`（979）——標数 2 の有限体では `c·x² = K` が常に解ける' +
      '（`x ↦ x²` が全単射、mathlib の `FiniteField.isSquare_of_char_two`）。',
    '★★★★★★★★★★★★★★★★`splits_or_exists_twist_splits\'`（979）——' +
      '**第 963 から `ringChar k ≠ 2` が落ちた**。' +
      '☆標数 2 では非平方元が無い代わりに 2 次式が必ず分裂する' +
      '（`IsCharNeTwoNF` により `a₁ = 0` なので 2 次式は `c₄X² − K`）。',
  ],
  measurement:
    '★★測定: 剰余標数 2 の心配は**非平方元の不在**の側ではなく、' +
    '**整モデルを `IsCharNeTwoNF`（`a₁ = a₃ = 0`）に正規化できるか**の側にある。' +
    '☆その正規化には `2` が可逆であることが要るので、`p ∣ 2` では別の手当てが要る。' +
    '★原文が「after a quadratic extension」と書くのはこの点であり、' +
    '第 919 で「捻りで不分岐 2 次拡大は不要」と書いた範囲は' +
    '**剰余標数が奇数の素点**である（訂正）。',
  remaining:
    '☆残るのは (i) 整モデルの `IsCharNeTwoNF` 正規化（`p ∤ 2` では変数変換 1 回、' +
    '`p ∣ 2` では不分岐 2 次拡大か別の議論）、' +
    '(ii) 979 の二者択一を `hasSplitMultiplicativeReduction_baseChange` に流す段、' +
    '(iii) 972 に全部並べて `minDeltaExp_variableChange` で戻す段。',
  note:
    '★944-979 の 36 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SplitCharTwo20260901 を書いた');
