// 2026-09-01（第 973）—— 極小モデルでは v_p(Δ) = −jExp。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35DeltaFromJExp20260901 = {
  measuredAt: '2026-09-01',
  block: '第 973',
  supersedes: 'lemma35GlobalVelu20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`valAdd_Delta_eq_neg_jExp`（973）——' +
      '極小モデル（`v_p(c₄) = 0`）では `v_p(Δ) = −jExp`。' +
      '☆`j = Δ⁻¹c₄³` を単元の等式に直し、`valAdd` の乗法性で割るだけ。',
    '★★★★★★★★`valAdd_Delta_pos_of_jExp_neg`（973）——' +
      '悪い素点では `0 < v_p(Δ)`。★これが第 909' +
      '（`hasMultiplicativeReduction_baseChange`）の仮説である。',
  ],
  remaining:
    '★`isMuAtBadPrimes_of_veluQuotient` の本体に残るのは:' +
    ' (i) `C • E` から `E` への変数変換の転送' +
    '（908/909 は `(C • E) ⊗ Lv` について述べるが、972 は `E ⊗ Lv` について要求する）、' +
    ' (ii) 分裂性（963 を当てるには整モデルを `IsCharNeTwoNF` に正規化する段が要る）、' +
    ' (iii) `hΔ`・`hcop`（932 経由）・`hlu`（953）・`hql`・`h2`・`h2K`・`hvw`（961＋962）。',
  note:
    '★944-973 の 30 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35DeltaFromJExp20260901 を書いた');
