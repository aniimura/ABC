// 2026-09-01（第 958）—— Vélu の項の対の和は偶。既存在庫の再発見も記録。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluPairEven20260901 = {
  measuredAt: '2026-09-01',
  block: '第 958',
  supersedes: 'lemma35PairEven20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`veluTerm_pair_even`（958）——' +
      '`(u_Q + 2v_Q x_Q) + (u_{-Q} + 2v_{-Q} x_Q) = 2(u_Q + x_Q(v_Q + v_{-Q}))`。' +
      '☆第 957 の `two_mul_sum_eq_of_pair_even` がそのまま消費する形である。',
  ],
  alreadyInStock:
    '★`veluGy_negY`・`veluU_negY`・`veluGx_negY`・`veluV_negY`・`veluW_negY` は' +
    '**すでに `Found/GenEll/Velu.lean:740-765` にあった**（第 958 で再発見）。' +
    '☆CLAUDE.md の「在庫」の通り、書く前に decl-index を引くべきだった。' +
    '`tools/lean-idioms.md` に失敗形として追記した。',
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d2) `tateXpair` の添字 `i ↦ l-i` と点の反転を繋ぐ段' +
      '（点の水準では `tatePhi_inv`(第 866) と `pointCoords_neg`(第 949) が揃っているが、' +
      '`hv`・`hw` は `tateXpair (ζ^i) (q(ζ^i)^{l-1})` の添字で書かれているので、' +
      'その対応を書く段が要る）、および `veluCurve` の楕円性。',
  ],
  note:
    '★944-958 の 15 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluPairEven20260901 を書いた');
