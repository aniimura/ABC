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

## ★残っている一歩

`Interface` の `ResidueCardinality` を discharge するには、あと
`isPrimePow`(q = p^f, f > 0)が要る。そのためには `CharP 𝓀[K] p` が必要:

- `(p : 𝓀[K]) = 0` — p は極大イデアルに入る(‖p‖ = 1/p < 1 なので `𝒪[K]` の単元でない)
- そこから `ringChar` が p を割り、p が素数かつ体が非自明なので `ringChar = p`
- 続いて `FiniteField.card` が `Nat.card = p ^ n` を与える

mathlib に `charP_of_prime_eq_zero` に相当する直接の補題は見当たらず(実測)、
`ringChar` 経由で数行書く必要がある。

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

end ABC3.Found.PGC
