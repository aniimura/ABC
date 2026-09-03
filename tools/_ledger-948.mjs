// 2026-09-01（第 947-948）—— 局所の段が 1 つの定理にまとまり、hcop が本物になった。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35LocalStep20260901 = {
  measuredAt: '2026-09-01',
  block: '第 947-948',
  supersedes: 'lemma35Rational20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`isPrimitiveRoot_of_pow_eq_one_of_ne_one` / `pow_ne_one_of_isPrimitiveRoot`（947）',
    '★★★★★★★★★★★★★★★★★★★★`exists_primitiveRoot_of_torsion_point`（947）——' +
      '有理な `l`-捩れ点から `tateParam_quot_velu_dvr`（927）が欲しい 6 つ' +
      '（`ζ : R`・`IsPrimitiveRoot ζ l`・`uζ`・`algebraMap ζ = uζ`・`uζ^l = 1`・`hord`）' +
      'を一気に作る。',
    '★★★★★★★★★★★★★★★★★★★★`tateParam_quot_velu_of_torsion`（948）——' +
      '**`ζ` が引数から消えた**。`P` が位数 `l` の点で `l ∤ v(q)` なら `q_{E′} = q_E^l`。' +
      '☆Vélu の帳簿（`hu`・`hv`・`hw`）だけが `ζ` について全称で残る。',
    '★`isMuAtBadPrimes_of_veluQuotient` の `hcop : True` を本物にした——' +
      '`∀ p, jExp p E < 0 → ¬ (l ∣ jExp p E)`（原文の「l is prime to the local heights」）。',
  ],
  remainingLayers: [
    '☆(D2) `H = ⟨Q⟩` の点集合を `pointCoords` の水準で Tate モデルの側に揃える段。' +
      '変数変換分は `j_veluCurve_variableChange`（926）で移す。',
    '☆(D3) 各悪い素点で 948 → 904 → 903 と流す段。' +
      '非分裂は 938/925/940/926/929 で捻りに移してから同じ連鎖を当てる。',
  ],
  note:
    '★第 943 の (A)(B)(C) は 945・944・945 で消え、Galois そのものが 946 で消え、' +
    '947-948 で局所の段が 1 つの定理にまとまった。' +
    '☆残るのは `E ⊗ Lv` と Tate モデルの間の変数変換の帳簿だけである。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35LocalStep20260901 を書いた');
