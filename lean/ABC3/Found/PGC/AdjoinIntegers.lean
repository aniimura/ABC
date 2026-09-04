import ABC3.Found.PGC.LubinTateDistinguishedSeparable

/-!
# `K.carrier⟮x⟯` の整数環を `Valued` を経由せず直接構成する(`sorry` 無し)

`Found/PGC/AdjoinPAdicLocalField.lean`(`K.carrier⟮x⟯` 自身を新たな
`PAdicLocalField p` として `ℚ_[p]` から再構成する経路)とは別の、
**より直接的な経路**。`IntermediateField.adjoin K.carrier {x}` が
`K.closure` の部分体として mathlib の一般論から自動的に持つ
`NormedField` 構造(`Found/PGC/LocalFieldNorm.lean::closureNormedField`
の制限)だけを使い、`Valued` 型クラスを一切導入せずに「整数環」
(ノルム `≤ 1` の元からなる部分環)を素朴に構成する。

## なぜ `Valued` を避けるか(2026-09-04 の実測)

`IntermediateField.adjoin K.carrier {x}` に `Valued _ NNReal :=
NormedField.toValued` を導入し、`Valued.integer`・
`Valued.integer.mem_iff`・`isCompact_closedBall` を組み合わせようと
すると、**120秒を超えても終わらない**深刻な単一化の詰まりに繰り返し
遭遇した(`maxHeartbeats` を200万まで上げても解決せず)。原因は
未特定だが、`IntermediateField extends Subfield extends Subring
extends Submonoid ...` という何層にも重なった部分構造の上で、
新しく導入した `Valued` 由来の位相と、既存の `NormedField` 由来の位相
(の定義的な一致)を検査するコストが高いためと推測される。

**対処**: `Valued` を一切使わず、`Subring.mk` で `{y | ‖y‖ ≤ 1}` を
直接構成する。この素朴な構成なら、`CompactSpace`・`IsClosed`・
`CompleteSpace` のすべてが**高速に**(1秒未満で)確認できた——
`Valued` 特有の位相の競合が存在しないため。

## 現状と次の一歩

`adjoinIntegers K x` の `CompactSpace`(`compactSpace_adjoinIntegers`)
まで確立できた。`Valued.integer.isDiscreteValuationRing_of_
compactSpace` は `Valued` 型クラスを要求するため、上記の理由でこの
`Subring` には直接使えない——`IsDiscreteValuationRing`(ひいては
`Ideal.isLinearTopology` 経由の `IsLinearTopology`)を`Valued` を経由
せずに得る経路が次の課題として残る。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `K.carrier⟮x⟯` の「整数環」——ノルムが `1` 以下の元からなる部分環。
`Valued` を経由しない素朴な構成(閉じ方は非アルキメデス距離の三角
不等式`IsUltrametricDist.norm_add_le_max`から)。 -/
noncomputable def adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    Subring (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) where
  carrier := {y | ‖y‖ ≤ 1}
  mul_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    calc ‖a * b‖ = ‖a‖ * ‖b‖ := norm_mul a b
      _ ≤ 1 * 1 := mul_le_mul ha hb (norm_nonneg b) zero_le_one
      _ = 1 := by ring
  one_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    haveI : IsUltrametricDist (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
      infer_instance
    calc ‖a + b‖ ≤ max ‖a‖ ‖b‖ := IsUltrametricDist.norm_add_le_max a b
      _ ≤ 1 := max_le ha hb
  zero_mem' := by simp
  neg_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq, norm_neg] at *
    exact ha

/-- `𝒪[K.carrier]` の元は `algebraMap` を通して `adjoinIntegers K x` に
入る——`spectralNorm_extends`(基点の元のスペクトルノルムはもとの
ノルムそのもの)と `Valued.integer.norm_le_one` から。`𝒪[K.carrier]`
が `𝒪[K.carrier⟮x⟯]` へ埋め込まれることの核心部分。 -/
theorem algebraMap_mem_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (y : 𝒪[K.carrier]) :
    (⟨algebraMap K.carrier K.closure (y : K.carrier),
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).algebraMap_mem _⟩ :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ∈ adjoinIntegers K x := by
  show spectralNorm K.carrier K.closure (algebraMap K.carrier K.closure (y : K.carrier)) ≤ 1
  rw [spectralNorm_extends]
  exact Valued.integer.norm_le_one y

/-- `Λ_n` の元 `x` 自身も(`K.carrier⟮x⟯` の元として)`adjoinIntegers K x`
に入る——`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints`
から直接。 -/
theorem mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    (⟨x, hmem⟩ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ∈ adjoinIntegers K x := by
  show spectralNorm K.carrier K.closure x ≤ 1
  exact (spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx).le

/-- `adjoinIntegers K x`(集合として)は閉単位球そのもの、したがって
`K.carrier⟮x⟯` の中で**閉集合**。 -/
theorem isClosed_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    IsClosed ((↑(adjoinIntegers K x) : Set (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))) := by
  have heq : (↑(adjoinIntegers K x) : Set (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) =
      Metric.closedBall 0 1 := by
    ext y
    simp only [Metric.mem_closedBall, dist_eq_norm, sub_zero]
    rfl
  rw [heq]
  exact Metric.isClosed_closedBall

/-- ★★★★★★★★★**`adjoinIntegers K x` は完備**——`K.carrier⟮x⟯` 自身が
完備(`FiniteDimensional.complete`、`Found/PGC/LubinTateDistinguished
Separable.lean::completeSpace_adjoin_of_mem_iteratedLubinTateTorsion
Points` と同じ経路)であることと、`adjoinIntegers K x` がその中で
閉集合であることから(`IsClosed.completeSpace_coe`)。 -/
theorem completeSpace_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    CompleteSpace (adjoinIntegers K x) := by
  haveI : CompleteSpace (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    FiniteDimensional.complete K.carrier _
  haveI := isClosed_adjoinIntegers K x
  exact IsClosed.completeSpace_coe

/-- ★★★★★★★★★**`adjoinIntegers K x` はコンパクト**——`K.carrier⟮x⟯` が
`K.carrier` 上有限次元であることから `ProperSpace`
(`FiniteDimensional.proper`)、そこで `adjoinIntegers K x` = 閉単位球
がコンパクト(`isCompact_closedBall`)。`Valued` を経由しないので、
`Valued.integer` で同じことをしようとしたときに遭遇した単一化の
詰まりを避けられる。 -/
theorem compactSpace_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    CompactSpace (adjoinIntegers K x) := by
  haveI : ProperSpace (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    FiniteDimensional.proper K.carrier _
  have heq : (↑(adjoinIntegers K x) : Set (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) =
      Metric.closedBall 0 1 := by
    ext y
    simp only [Metric.mem_closedBall, dist_eq_norm, sub_zero]
    rfl
  have hcompact : IsCompact (↑(adjoinIntegers K x) :
      Set (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    heq ▸ isCompact_closedBall 0 1
  exact isCompact_iff_compactSpace.mp hcompact

end ABC3.Found.PGC
