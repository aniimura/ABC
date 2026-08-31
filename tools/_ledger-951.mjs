// 2026-09-01（第 951）—— 分円の等式で仮説 hu が消えた。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Cyclotomic20260901 = {
  measuredAt: '2026-09-01',
  block: '第 951',
  supersedes: 'lemma35PointSet20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`range_erase_zero_eq_image_succ` / `prod_one_sub_pow_erase`（951）——' +
      'mathlib の `IsPrimitiveRoot.prod_one_sub_pow_eq_order` を ' +
      '`(range l).erase 0` の形に乗せ換え、`∏_{i=1}^{l-1} (1 - ζ^i) = l` を得る。',
    '★★★★★★★★★★★★`isUnit_one_sub_pow_of_isUnit_natCast`（951）——' +
      '積が `l` に等しいので、`l` が単元なら各因子 `1 - ζ^i` も単元。',
    '★`tateParam_quot_velu_of_torsion`（948）から**仮説 `hu` を落とした**。' +
      '残る Vélu の帳簿は `hv`・`hw` の 2 つだけである。',
  ],
  remainingLayers: [
    '☆(D3) 各悪い素点で 950 → 948 → 904 → 903 と流す段。' +
      '★組み立てが供給すべきものは' +
      ' (a) `E ⊗ Lv` の分裂乗法還元と極小性、' +
      ' (b) `hp`（付値の両立）と `C`・`hC`・`hc4`（極小モデルの `c₄` が単元）、' +
      ' (c) `hlu : IsUnit (l : R)`（残余標数が `l` でないこと）、' +
      ' (d) `v`・`w` を和として定義する段。' +
      '☆`hu` は 951 で不要になった。',
  ],
  note:
    '★第 943 の (A)(B)(C) は 945・944・945 で消え、Galois は 946 で消え、' +
    '局所の段は 947-948 で 1 つの定理になり、点集合の帳簿は 949-950 で閉じ、' +
    '分円の仮説は 951 で消えた。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Cyclotomic20260901 を書いた');
