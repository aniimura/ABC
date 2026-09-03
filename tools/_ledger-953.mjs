// 2026-09-01（第 953）—— (D3) の (c) が閉じた：hlu は v_p(l) = 0 から出る。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35UnitAtPrime20260901 = {
  measuredAt: '2026-09-01',
  block: '第 953',
  supersedes: 'lemma35VeluData20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★`isUnit_of_vAdd_eq_zero`（953）——付値が `0` の元は環の単元。' +
      '`x` も `x⁻¹` も環に入る（`dvr_mem_of_nonneg`、第 895）から。',
    '★★★★★★★★★★★★`isUnit_natCast_of_valAdd_eq_zero`（953）——' +
      '`v_p(n) = 0` なら `n` は完備化の整数環で単元。' +
      '`vAdd_algebraMap_eq_valAdd`（第 893）で `vAdd = valAdd` に落とすだけ。',
    '☆これが `tateParam_quot_velu_of_torsion`（948）の `hlu : IsUnit (l : R)` の出どころ——' +
      '原文の「l is prime to the residue characteristics」が `v_p(l) = 0` に当たる。',
  ],
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂乗法還元と極小性、' +
      '(b) `hp`（付値の両立）と `C`・`hC`・`hc4`、' +
      '(d) 各 `ζ` に対する `v`・`w` と `veluCurve` の楕円性。' +
      '★(c) は 953 で、`hu` は 951 で閉じた。',
  ],
  note:
    '★944-953 の 10 ブロックで節点の界面は局所データ 3 つに絞られた。' +
    '`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35UnitAtPrime20260901 を書いた');
