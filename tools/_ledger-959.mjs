// 2026-09-01（第 959）—— (d2) が閉じた：添字 i ↦ l-i は点の反転である。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35MuInvolution20260901 = {
  measuredAt: '2026-09-01',
  block: '第 959',
  supersedes: 'lemma35VeluPairEven20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`map_negY_algebraMap`（959）——`negY` は底変換と可換。',
    '★`tateMu_pointCoords`（959）——添字 `i` の `μ_l` の点の座標。第 890 の中身を取り出した。',
    '★★★★★★★★★★★★★★★★`tateXpair_mu_inv`（959）——' +
      '**`X_{l−i} = X_i` かつ `Y_{l−i} = negY X_i Y_i`**。',
  ],
  howItWasDone:
    '★級数の側で取ろうとすると難しい——`tateXpair_symm` は `X(q/u) = X(u)` であり、' +
    '欲しいのは `X(1/u) = X(u)` だから、周期性 `u ∼ qu` を級数の側で通す必要がある。' +
    '☆**点を経由すれば 1 行**である: `(l−i) • P = −(i • P)`（`l • P = 0` だから、第 949）。' +
    '座標は `pointCoords_tatePtPair`（912）が `tateXpair`・`tateYpair` で与え、' +
    '反転の座標は `pointCoords_neg`（949）が `(x, negY x y)` で与える。' +
    '`algebraMap R K` は単射（`S.hinj`）なので `R` の等式に降りる。',
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d3) 959 と 958 と 957 を並べて `hvw` を作る段（`veluCurve` の楕円性を含む）。',
  ],
  note:
    '★944-959 の 16 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35MuInvolution20260901 を書いた');
