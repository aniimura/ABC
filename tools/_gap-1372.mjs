import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
j.integralJacobson20260902 =
  '★★★★★★★★**第 1365-1366**——`Found/GenEll/IntegralJacobson.lean` と '
  + '`Found/GenEll/CompleteIntegralClosure.lean`。'
  + '`map_maximalIdeal_le_jacobson`（整拡大では `m_R·S ≤ jacobson ⊥`）と '
  + '`isDiscreteValuationRing_integralClosure`'
  + '（**完備 DVR の有限分離整閉包は DVR**）。'
  + '☆第 1356 で「mathlib に無い」と測った primitive がこれで塞がった。';
j.adicComparable20260902 =
  '★★★★★★**第 1367-1368**——`Found/GenEll/AdicComparable.lean`。'
  + '`isAdicComplete_of_pow_le`（`I ≤ J` かつ `J^e ≤ I` なら完備性が移る）で '
  + '`m_A C`-進完備を `m_C`-進完備に直す。'
  + '☆mathlib には「同じ位相を与えるイデアルの間で `IsAdicComplete` が移る」が無い（2026-09-02 実測）。';
j.ramificationExponent20260902 =
  '★★★★★★★★**第 1369**——`exists_pow_eq_map_maximalIdeal`（`m_A C = m_C^e`、`1 ≤ e`）と '
  + '`le_finrank_of_pow_eq_map_maximalIdeal`（`e ≤ [L:K]`、mathlib の '
  + '`Ideal.ramificationIdx_le_finrank` に橋渡し）。'
  + '☆`[L:K] < l` なら `l ∤ e` が出る。';
j.hpRamified20260902 =
  '★★★★★★★★★★★★**第 1370-1372**——`Found/GaloisRep/RamifiedValuationBridge.lean`・'
  + '`RamifiedBadPrime.lean`・`Found/GenEll/AlphaFromBadPrimeRam.lean`。'
  + '★測定: `exists_h2_h1_of_bad_prime`（第 1320）で不分岐の仮定 `hp` が使われるのは '
  + '`vAdd_algebraMap_eq_valAdd` の**ただ 1 箇所**だった'
  + '（他はインスタンス引数なので呼び出し側が別途与える）。'
  + '☆その 1 本を `e` 倍版 `vAdd_algebraMap_eq_mul_valAdd` にしただけで '
  + '`exists_h2_h1_of_bad_prime_ram` まで通った。'
  + '★★★これで `L_v` として **`L_p` の任意の有限拡大**が使える。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('gap entries recorded');
