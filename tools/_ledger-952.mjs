// 2026-09-01（第 952）—— 948 の hv・hw を存在量化に直した（充足可能性の訂正）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluData20260901 = {
  measuredAt: '2026-09-01',
  block: '第 952',
  supersedes: 'lemma35Cyclotomic20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  correction:
    '★★第 948 の `hv : ∀ ζ, IsPrimitiveRoot ζ l → v = ∑ …(ζ)` は**充足不能に近かった**——' +
    '1 つの `v` が全ての原始根について成り立つことを要求しており、' +
    '和が `ζ` に依らないことを別に証明しない限り満たせない。' +
    '☆第 952 で `hvw : ∀ ζ, IsPrimitiveRoot ζ l → ∃ v w, … ∧ … ∧ IsElliptic …` に直した。' +
    '★補助データ（`v`・`w`）は結論に現れないので、存在量化を内側に入れるのが正しい。',
  remainingLayers: [
    '☆(D3) 各悪い素点で 950 → 948 → 904 → 903 と流す段。組み立てが供給すべきもの:' +
      ' (a) `E ⊗ Lv` の分裂乗法還元と極小性、' +
      ' (b) `hp`（付値の両立）と `C`・`hC`・`hc4`（極小モデルの `c₄` が単元）、' +
      ' (c) `hlu : IsUnit (l : R)`（残余標数が `l` でないこと）、' +
      ' (d) 各 `ζ` に対する `v`・`w` と `veluCurve` の楕円性。' +
      '☆`hu` は 951 で不要になった。',
  ],
  note:
    '★944-952 の 9 ブロックで節点の周辺は片付いた。' +
    '`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluData20260901 を書いた');
