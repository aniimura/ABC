// 2026-09-01（第 956）—— Vélu の w を作る鍵：対称な和は 2 で割れる。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SymmSum20260901 = {
  measuredAt: '2026-09-01',
  block: '第 956',
  supersedes: 'lemma35JForm20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`erase_zero_range_eq_union` / `disjoint_Icc_image_sub`（956）——' +
      '`l = 2m+1` のとき `(range l).erase 0` は `Icc 1 m` とその反転像の素分割。',
    '★★★★★★★★★★★★`sum_eq_two_nsmul_of_symm`（956）——' +
      '`i ↦ l-i` で対称な関数の和は `2 • (半分の和)`。' +
      '☆`l` が奇数だから不動点がなく、添字集合がちょうど 2 つずつの組に分かれる。',
    '★★★★★★★★`exists_two_mul_of_symm`（956）——環での `2 * w = ∑` の形。' +
      '☆これが Vélu の `w` を作る形である（`hw` は `2` で割らずに書いてある）。',
  ],
  whyItMatters:
    '★V\u00e9lu の `w` は `2 * w = ∑ (veluU + 2·veluV2·x)` という形で受けている。' +
    '`2` で割れない環でも `w` を作るには、右辺が 2 の倍数であることが要る。' +
    '☆`μ_l` の点は反転で対になり、`veluU = (2y + x)²`（Tate 曲線は `a₁ = 1`, `a₃ = 0`）は' +
    '反転 `y ↦ -y-x` で不変だから、和は対称である。' +
    '★956 はその「対称なら 2 で割れる」を楕円曲線と無関係な `Finset` の事実として取った。',
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d) `veluU` の反転不変性を 956 に繋ぐ段と `veluCurve` の楕円性。',
  ],
  note:
    '★944-956 の 13 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SymmSum20260901 を書いた');
