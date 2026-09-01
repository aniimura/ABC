// 2026-09-01（第 975）—— j が同じ相方に付け替える（第 974 の食い違いを吸収）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35JSwap20260901 = {
  measuredAt: '2026-09-01',
  block: '第 975',
  supersedes: 'lemma35VeluVcTransport20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`exists_point_j_tateModel\'`（975）——' +
      '`E″` が `E′` と同じ `j` を持てば `hW′j` は `E″` についても成り立つ。' +
      '☆これで第 974 で測った食い違い' +
      '（Vélu の関係は `C • E′` の側、極小性は `C′ • E′` の側）が吸収できる。' +
      '★route (α) の要石である。',
  ],
  routeAlpha: [
    '☆本体の道筋（route α、部品はすべて在庫）:',
    '(1) 954 で `C`・`hC`・`hc4`（`E`）と `C′`・`hC′`・`hc4′`（`E′`）。',
    '(2) 973 で `0 < v_p(Δ)`、908 で `(C • E) ⊗ Lv` と `(C′ • E′) ⊗ Lv` の極小性。',
    '(3) 909 ＋ 963 で分裂乗法還元（非分裂は 938 → 929）。',
    '(4) 974 で `C • E′ = (C • E)/⟨C • Q⟩`。',
    '(5) 975 で `hW′j` を `C′ • E′` に付け替える' +
      '（`j(C • E′) = j(C′ • E′)`——`variableChange_j` を 2 回）。',
    '(6) 964・953・932・961・962 で残りの局所データ。',
    '(7) 972 に渡し、`minDeltaExp_variableChange` で `E`・`E′` に戻す。',
  ],
  note:
    '★944-975 の 32 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35JSwap20260901 を書いた');
