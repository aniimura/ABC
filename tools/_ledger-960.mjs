// 2026-09-01（第 960）—— (d3)：Vélu の w が環の中で作れる。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluW20260901 = {
  measuredAt: '2026-09-01',
  block: '第 960',
  supersedes: 'lemma35MuInvolution20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`exists_veluW_of_inv`（960）——' +
      '添字の反転が点の反転（`X_{l−i} = X_i`、`Y_{l−i} = negY X_i Y_i`）なら、' +
      '`∃ w, 2 * w = ∑ (u_Q + 2 v_Q x_Q)`。' +
      '☆三つを並べるだけ: 957（対ごとに偶なら和も偶）＋958（Vélu の項の対は偶）＋' +
      '仮説（959 が与える）。',
  ],
  dDone:
    '★(d) の 4 枚がすべて揃った: 957（Finset の対分割）・958（Vélu の項）・' +
    '959（添字と点の反転）・960（組み立て）。' +
    '☆残るのは 960 を具体的な `tateXpair`/`tateYpair` に当てて `hvw` を作る呼び出しと、' +
    '`veluCurve` の楕円性である。',
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d4) 960 を `tateXpair` に当てる呼び出しと `veluCurve` の楕円性。',
  ],
  note:
    '★944-960 の 17 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluW20260901 を書いた');
