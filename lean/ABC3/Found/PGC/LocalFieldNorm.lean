import ABC3.Skeleton.PGC.Setup
import Mathlib.Analysis.Normed.Unbundled.SpectralNorm
import Mathlib.NumberTheory.Padics.ProperSpace
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Topology.Algebra.Valued.NormedValued
import ABC3.Found.ResidueFieldFinite

/-!
# Track B — p進局所体の付値構造(スペクトルノルム経由)

`PAdicLocalField p`(= `ℚ_[p]` の有限次拡大)に、mathlib の道具だけで
ノルム体・付値体の構造を与える。**`sorry` 無し。**

これは `Interface/PGC/LocalFieldData.lean` の `ResidueCardinality` を
実装するための土台。連鎖:

```
spectralNorm.normedField        → NormedField K
isNonarchimedean_spectralNorm   → IsUltrametricDist K
spectralNorm.normedAlgebra      → NormedAlgebra ℚ_[p] K
‖p⁻¹‖ = p > 1                   → NontriviallyNormedField K
NormedField.toValued            → Valued K ℝ≥0
NormedField の RankOne インスタンス → (Valued.v).RankOne
FiniteDimensional.proper        → ProperSpace K      (ℚ_[p] が proper なので)
```

## 到達点(2026-08-14)

剰余体の有限性まで到達した(`residueField_finite`)。当初はノルムのダイヤモンドで
止まっていたが、**round-trip 補題は不要だった**——`CompactSpace 𝒪[K]` を経由すれば
問題はノルムでなく**位相**の一致に落ち、そこは定義的に一致する。
詳細は `ABC3/Found/ResidueFieldFinite.lean` の docstring。

## 到達点(2026-08-14、続き)

`isPrimePow`(q = p^f, f > 0)まで到達し、`Interface` の `ResidueCardinality` を
discharge した(`Found/PGC/ResidueCardinality.lean`)。要ったのは
**`‖(p : K)‖ = 1/p < 1`** ただ一つ——そこから

- p は `𝒪[K]` の単元でない → 極大イデアルに入る → `(p : 𝓀[K]) = 0`
- `CharP.charP_iff_prime_eq_zero` で `CharP 𝓀[K] p`
- `FiniteField.card` で `Nat.card 𝓀[K] = p ^ f`(`f > 0`)

と繋がる。一般部分は `Found/ResidueFieldFinite.lean` に置いた。

★**先行する記述の訂正**: ここには以前
「mathlib に `charP_of_prime_eq_zero` に相当する直接の補題は見当たらず(実測)」
と書いてあったが**誤り**だった。`CharP.charP_iff_prime_eq_zero`
(`Mathlib/Algebra/CharP/Basic.lean:103`)が存在する。経緯は
`Found/ResidueFieldFinite.lean` の docstring に記録した。

## 使い方

インスタンスは `scoped` にしてある。使う側で `open scoped ABC3.Found.PGC` すること
(無条件に大域インスタンスにすると、他所のノルム構造と衝突しうる)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] Algebra.IsAlgebraic.of_finite

/-- スペクトルノルムによる `K` のノルム体構造。 -/
@[implicit_reducible] noncomputable scoped instance normedField (K : PAdicLocalField p) :
    NormedField K.carrier :=
  spectralNorm.normedField ℚ_[p] K.carrier

/-- スペクトルノルムは非アルキメデス的。 -/
noncomputable scoped instance isUltrametric (K : PAdicLocalField p) :
    IsUltrametricDist K.carrier :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := ℚ_[p]) (L := K.carrier))

/-- `K` は `ℚ_[p]` 上のノルム代数。 -/
@[implicit_reducible] noncomputable scoped instance normedAlgebra (K : PAdicLocalField p) :
    NormedAlgebra ℚ_[p] K.carrier :=
  spectralNorm.normedAlgebra ℚ_[p] K.carrier

/-- スペクトルノルムは `ℚ_[p]` のノルムを延長する。 -/
theorem norm_algebraMap (K : PAdicLocalField p) (x : ℚ_[p]) :
    ‖algebraMap ℚ_[p] K.carrier x‖ = ‖x‖ :=
  spectralNorm_extends x

/-- `‖p⁻¹‖ = p > 1` なのでノルムは自明でない。 -/
@[implicit_reducible] noncomputable scoped instance nontriviallyNormedField
    (K : PAdicLocalField p) : NontriviallyNormedField K.carrier where
  toNormedField := normedField K
  non_trivial := by
    refine ⟨algebraMap ℚ_[p] K.carrier ((p : ℚ_[p])⁻¹), ?_⟩
    rw [norm_algebraMap, norm_inv, Padic.norm_p, inv_inv]
    exact_mod_cast (Fact.out : p.Prime).one_lt

/-- ノルムから定まる付値体構造。 -/
@[implicit_reducible] noncomputable scoped instance valued (K : PAdicLocalField p) :
    Valued K.carrier NNReal :=
  NormedField.toValued

/-- 付値は階数 1。 -/
@[implicit_reducible] noncomputable scoped instance rankOne (K : PAdicLocalField p) :
    (Valued.v : Valuation K.carrier NNReal).RankOne :=
  inferInstanceAs (NormedField.valuation (K := K.carrier)).RankOne

/-- `ℚ_[p]` は proper なので、その有限次拡大も proper。 -/
noncomputable scoped instance properSpace (K : PAdicLocalField p) :
    ProperSpace K.carrier :=
  FiniteDimensional.proper ℚ_[p] K.carrier

open scoped NormedField Valued in
/-- **p進局所体の剰余体は有限**。`Found/ResidueFieldFinite.lean` の一般結果を適用。 -/
theorem residueField_finite (K : PAdicLocalField p) : Finite 𝓀[K.carrier] :=
  finite_residueField (K := K.carrier)

open scoped NormedField Valued in
/-- 剰余体の元の個数 `q`。`Interface` の `ResidueCardinality.card` の実体。 -/
noncomputable def residueCard (K : PAdicLocalField p) : ℕ := Nat.card 𝓀[K.carrier]

/-- `‖p‖ = 1/p < 1`。剰余体の標数が `p` であることの唯一の入力。

`p` の像は `algebraMap ℚ_[p] K` を通るので、スペクトルノルムが `ℚ_[p]` の
ノルムを延長すること(`norm_algebraMap`)から直ちに従う。 -/
theorem norm_natCast_p_lt_one (K : PAdicLocalField p) : ‖((p : ℕ) : K.carrier)‖ < 1 := by
  rw [show ((p : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((p : ℕ) : ℚ_[p]) from
        (map_natCast _ p).symm,
      norm_algebraMap, Padic.norm_p]
  exact inv_lt_one_of_one_lt₀ (by exact_mod_cast (Fact.out : p.Prime).one_lt)

open scoped NormedField Valued in
/-- **原文 [pGC] p.3「k is the field of q = p^f elements」**。

`Interface` の `ResidueCardinality.isPrimePow` の実体。`0 < f` 込み。 -/
theorem residueCard_isPrimePow (K : PAdicLocalField p) :
    ∃ f : ℕ, 0 < f ∧ residueCard K = p ^ f :=
  card_residueField_eq_prime_pow (K := K.carrier) Fact.out (norm_natCast_p_lt_one K)

/-! ### `K.closure`(代数閉包)にもスペクトルノルムでノルム体構造を与える

`K.carrier` に `ℚ_[p]` から延長したのと同じ手順を、今度は `K.carrier` を
基点として `K.closure := AlgebraicClosure K.carrier` へ繰り返す。
`spectralNorm.normedField`/`normedAlgebra` は基点が完備でありさえすれば
有限次拡大に限らず**任意の代数拡大**に対して働くので、`K.closure` が
`K.carrier` 上有限次でなくても(実際、代数閉包は一般に無限次)問題ない。
これは `Λ_n` の元(`K.carrier` 上有限次拡大 `K.carrier⟮x⟯` を張る)を
`PowerSeries.aeval` で実際に評価するのに必要な `CompleteSpace` へ
繋がる一歩。 -/

/-- スペクトルノルムによる `K.closure` のノルム体構造。 -/
@[implicit_reducible] noncomputable scoped instance closureNormedField (K : PAdicLocalField p) :
    NormedField K.closure :=
  spectralNorm.normedField K.carrier K.closure

/-- `K.closure` は `K.carrier` 上のノルム代数。 -/
@[implicit_reducible] noncomputable scoped instance closureNormedAlgebra (K : PAdicLocalField p) :
    NormedAlgebra K.carrier K.closure :=
  spectralNorm.normedAlgebra K.carrier K.closure

/-- `K.closure` 上のスペクトルノルムも非アルキメデス的。 -/
noncomputable scoped instance closureIsUltrametric (K : PAdicLocalField p) :
    IsUltrametricDist K.closure :=
  IsUltrametricDist.isUltrametricDist_of_forall_norm_add_le_max_norm
    (isNonarchimedean_spectralNorm (K := K.carrier) (L := K.closure))

end ABC3.Found.PGC
