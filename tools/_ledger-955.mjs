// 2026-09-01（第 955）—— 950 の j 出力をそのまま受ける捩れ点版。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35JForm20260901 = {
  measuredAt: '2026-09-01',
  block: '第 955',
  supersedes: 'lemma35BadPrimeData20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`tateParam_quot_velu_j_of_torsion`（955）——' +
      '第 948 は `W′` が Vélu の商に**等しい**ことを要求していたが、' +
      '実際に得られるのは変数変換を除いた `j` の一致である。' +
      '☆第 950 がその `j` を与えるので、本定理がそのままの形で受けられる。',
    '★`ζ` は引数に現れない（947 が作る）。`hu` も現れない（951 が作る）。',
  ],
  scale: {
    note: '★`Lemma 3.5` の規模の実測（2026-09-01）',
    lemma35Decls: 256,
    lemma32Decls: 213,
    def33Decls: 293,
    files: 194,
    lines: 53811,
    comment:
      '原文では半ページだが、Tate 一意化（Definition 3.3）と Lemma 3.2 を' +
      '土台として抱えている。3 つ合わせて 762 宣言・194 ファイル・53,811 行。',
  },
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）、' +
      '(d) 各 `ζ` に対する `v`・`w` と `veluCurve` の楕円性。',
  ],
  note:
    '★944-955 の 12 ブロックで節点の界面は 3 つの局所データに絞られた。' +
    '`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35JForm20260901 を書いた');
