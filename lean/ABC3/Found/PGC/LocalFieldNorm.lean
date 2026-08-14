import ABC3.Skeleton.PGC.Setup
import Mathlib.Analysis.Normed.Unbundled.SpectralNorm
import Mathlib.NumberTheory.Padics.ProperSpace
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Topology.Algebra.Valued.NormedValued

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

## ★残っている一歩(2026-08-14 時点)

剰余体の有限性は mathlib の
`Valued.integer.properSpace_iff_completeSpace_and_isDiscreteValuationRing_integer_and_finite_residueField`
から出るはずだが、**ノルムのダイヤモンド**で止まっている——
この補題は `Valued.toNormedField`(付値から作ったノルム)を前提に述べられているのに対し、
こちらの `ProperSpace` はスペクトルノルムに関するもので、両者は
**命題としては等しいが定義的には等しくない**。

必要なのは round-trip 補題
`Valued.toNormedField (NormedField.toValued (K := K)) = ‹NormedField K›`
(あるいはノルムが一致することの補題)で、**mathlib に見当たらない**(実測)。
これを書くのが次の一歩。有限の・よく定義された作業である。

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

end ABC3.Found.PGC
