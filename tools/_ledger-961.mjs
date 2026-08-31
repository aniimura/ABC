// 2026-09-01（第 961）—— (d4) の hw が出た。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluWMu20260901 = {
  measuredAt: '2026-09-01',
  block: '第 961',
  supersedes: 'lemma35VeluW20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`exists_veluW_mu`（961）——' +
      '`μ_l` の Vélu の `w` は環の中で作れる。' +
      '☆`exists_veluW_of_inv`（960）に `tateXpair_mu_inv`（959）を当てるだけ。' +
      '★`l` を `2m+1` の形で受ける——奇数であることが対分けの本質だから。',
    '☆これで `tateParam_quot_velu_of_torsion`（948）の `hvw` のうち **`hw` が出た**。' +
      '`hv` は `v := ∑ …` と置けば定義的に成り立つ。',
  ],
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d5) `veluCurve (tateCurveAt q hq) v w` の楕円性。',
  ],
  note:
    '★944-961 の 18 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluWMu20260901 を書いた');
