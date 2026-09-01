// 2026-09-01（第 970）—— Tate モデルの上の点と j の一致を大域データから作る。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35TateModelPoint20260901 = {
  measuredAt: '2026-09-01',
  block: '第 970',
  supersedes: 'lemma35VeluEllipticFromEprime20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`exists_point_j_tateModel`（970）——' +
      '大域の Vélu の商から、Tate モデルの上の位数 `l` の点 `P` と `j` の一致を作る。' +
      '☆これで `minDeltaExp_eq_mul_of_torsion`（965）の ' +
      '`P`・`hP`・`hP0`・`hW′j` がすべて揃う。',
  ],
  fiveSteps: [
    '(1) `tateParamR_spec` ＋ 944 で `(E_q) ⊗ Lv = C₀ • (E ⊗ Lv)`',
    '(2) 969 で Vélu の商の楕円性（循環しない）',
    '(3) 967 で `j` の一致（`C₀ • (E ⊗ Lv)` 側）',
    '(4) 在庫の `addOrderOf_rhPoint`・`addOrderOf_vcPoint` で点の位数を運ぶ',
    '(5) 968 で曲線の等式に沿って点を Tate モデル側へ運ぶ',
  ],
  plumbing:
    '★`.j` の上で `rw [himg, hbase]` をやると `IsElliptic` で motive が壊れる——' +
    '曲線の等式を先に作って `j_congr_curve`（913）で渡す（944 と同じ穴）。' +
    '★楕円性を `∧` の中で運ぶと後続の `.j` がインスタンスとして使えないので、' +
    '結論を `∀ _hell : …, j の等式` の形にした。' +
    '★`vcPoint` は `GaloisRep` と `GenEll` の両方にあり、' +
    '`namespace ABC3.Found.GaloisRep` の中では前者が選ばれる。名前空間を明示すること。',
  remainingLayers: [
    '☆残るのは 965 に 970 を繋ぎ、各悪い素点で (2)-(4)' +
      '（954・909・963・929・964・953・932・961・962）を供給する段。',
  ],
  note:
    '★944-970 の 27 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35TateModelPoint20260901 を書いた');
