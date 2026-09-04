import ABC3.Found.PGC.LubinTateDistinguishedSeparable
import ABC3.Found.PGC.LubinTateActionEndomorphism
import Mathlib.RingTheory.PowerSeries.PiTopology
import Mathlib.RingTheory.MvPowerSeries.PiTopology

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

## 現状(2026-09-04、★★★★★★★★★★★★★節目——`PowerSeries.aeval` が
`Λ_n` の元で実際に組み立てられた)

`IsDiscreteValuationRing`(そこから `Ideal.isLinearTopology` 経由で
`IsLinearTopology` を得る、という当初の計画)を`Valued` を経由せずに
得る経路の代わりに、**もっと直接的な2つの経路**を見つけた:

1. **`ValuationRing`・`IsLocalRing`**: `adjoinIntegers K x` の任意の
   2元 `a,b` について、ノルムの大小で場合分けし、小さい方を大きい方
   で割った商(ノルム`≤1`なので`adjoinIntegers K x`自身の元)を
   `PreValuationRing`の証人として直接与えるだけ——`Valued`・
   コンパクト性のどちらも不要。`ValuationRing.iff_local_bezout_domain`
   で`IsLocalRing`も従う。
2. **`IsLinearTopology`**: 半径`ε>0`の閉球`{y | ‖y‖≤ε}`が(非
   アルキメデス三角不等式+ノルム`≤1`のスカラー倍で)そのまま
   `adjoinIntegers K x`の**イデアル**になることを直接示し
   (`adjoinIntegersBall`)、`Metric.nhds_basis_closedBall`と
   `IsLinearTopology.mk_of_hasBasis`を組み合わせる——`Ideal.
   isLinearTopology`(adic位相経由、`IsDiscreteValuationRing`が要る)
   も`Valued`も一切経由しない。

これに`𝒪[K.carrier]`から`adjoinIntegers K x`への埋め込み
(`adjoinIntegersAlgebraMap`、`spectralNorm_extends`から等長)と、その
連続性(`continuousSMul_adjoinIntegers`)を組み合わせると、
`PowerSeries.aeval`の要求する条件が**すべて**揃い、`Λ_n`の元`x`で
実際に冪級数を評価する`lubinTateEvalAtTorsionPoint`が完成した。
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

/-! ### `adjoinIntegers K x` は付値環——`Valued`・コンパクト性を経由しない -/

/-- `adjoinIntegers K x` は `PreValuationRing`——任意の2元 `a,b` の
ノルムを比較し、小さい方を大きい方で割った商(ノルム`≤1`なので
`adjoinIntegers K x`自身の元)を証人として直接与える。 -/
theorem preValuationRing_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    PreValuationRing (adjoinIntegers K x) := by
  constructor
  intro a b
  set L := IntermediateField.adjoin K.carrier ({x} : Set K.closure)
  rcases eq_or_ne (a : L) 0 with ha0 | ha0
  · refine ⟨0, Or.inr ?_⟩
    apply Subtype.ext
    show (b : L) * (0 : L) = (a : L)
    rw [mul_zero, ha0]
  · have hanorm : 0 < ‖(a : L)‖ := norm_pos_iff.mpr ha0
    rcases le_total ‖(b : L)‖ ‖(a : L)‖ with hle | hle
    · have hcmem : ‖(b : L) / (a : L)‖ ≤ 1 := by
        rw [norm_div, div_le_one hanorm]; exact hle
      refine ⟨⟨(b : L) / (a : L), hcmem⟩, Or.inl ?_⟩
      apply Subtype.ext
      show (a : L) * ((b : L) / (a : L)) = (b : L)
      rw [mul_div_cancel₀ _ ha0]
    · have hbnorm : 0 < ‖(b : L)‖ := lt_of_lt_of_le hanorm hle
      have hbne0 : (b : L) ≠ 0 := norm_pos_iff.mp hbnorm
      have hcmem : ‖(a : L) / (b : L)‖ ≤ 1 := by
        rw [norm_div, div_le_one hbnorm]; exact hle
      refine ⟨⟨(a : L) / (b : L), hcmem⟩, Or.inr ?_⟩
      apply Subtype.ext
      show (b : L) * ((a : L) / (b : L)) = (a : L)
      rw [mul_div_cancel₀ _ hbne0]

/-- ★★★★★★★★★**`adjoinIntegers K x` は付値環**——`PreValuationRing`
から直ちに。 -/
theorem valuationRing_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    ValuationRing (adjoinIntegers K x) :=
  haveI := preValuationRing_adjoinIntegers K x
  { }

/-- `adjoinIntegers K x` は局所環——`ValuationRing.iff_local_
bezout_domain`から。 -/
theorem isLocalRing_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    IsLocalRing (adjoinIntegers K x) :=
  haveI := valuationRing_adjoinIntegers K x
  (ValuationRing.iff_local_bezout_domain.mp inferInstance).1

/-! ### `IsLinearTopology`——半径 `ε` の閉球が直接イデアルになる -/

/-- 半径 `ε` の「閉球」——`adjoinIntegers K x` のイデアル(`ε≤0` の
場合は退化して `{0}` 寄りになるが、`ε>0` でのみ使うので問題ない)。
非アルキメデス三角不等式(`add_mem'`)とノルム`≤1`のスカラー倍
(`smul_mem'`)から直接イデアルの公理を満たす——`Valued`・
`IsDiscreteValuationRing` を一切経由しない。 -/
noncomputable def adjoinIntegersBall {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (ε : ℝ) : Ideal (adjoinIntegers K x) where
  carrier := {y | ‖(y : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ max ε 0}
  zero_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    haveI : IsUltrametricDist (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
      infer_instance
    calc ‖(a : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) +
        (b : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖
        ≤ max ‖(a : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖
              ‖(b : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ :=
          IsUltrametricDist.norm_add_le_max _ _
      _ ≤ max ε 0 := max_le ha hb
  smul_mem' := by
    intro c a ha
    simp only [Set.mem_setOf_eq] at *
    show ‖(c : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) *
        (a : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ max ε 0
    rw [norm_mul]
    calc ‖(c : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ *
        ‖(a : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖
        ≤ 1 * max ε 0 := mul_le_mul c.2 ha (norm_nonneg _) zero_le_one
      _ = max ε 0 := one_mul _

/-- `ε>0` のとき、`adjoinIntegersBall K x ε`(集合として)はちょうど
半径`ε`の閉球そのもの。 -/
theorem coe_adjoinIntegersBall {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (ε : ℝ) (hε : 0 < ε) :
    (↑(adjoinIntegersBall K x ε) : Set (adjoinIntegers K x)) = Metric.closedBall 0 ε := by
  ext y
  simp only [SetLike.mem_coe, adjoinIntegersBall, Submodule.mem_mk, AddSubmonoid.mem_mk,
    AddSubsemigroup.mem_mk, Set.mem_setOf_eq, max_eq_left hε.le, Metric.mem_closedBall,
    dist_eq_norm, sub_zero]
  show ‖(y : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ ε ↔ ‖y‖ ≤ ε
  rfl

/-- `adjoinIntegers K x` における `0` の近傍フィルターは、
`adjoinIntegersBall K x ε`(`ε>0`)を基底に持つ——`Metric.nhds_basis_
closedBall`(距離空間の一般論)と `coe_adjoinIntegersBall` から。 -/
theorem hasBasis_nhds_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    (nhds (0 : adjoinIntegers K x)).HasBasis (fun ε : ℝ => 0 < ε)
      (fun ε => (↑(adjoinIntegersBall K x ε) : Set (adjoinIntegers K x))) := by
  have h := Metric.nhds_basis_closedBall (x := (0 : adjoinIntegers K x))
  exact h.congr (fun _ => Iff.rfl) (fun ε hε => (coe_adjoinIntegersBall K x ε hε).symm)

/-- ★★★★★★★★★★★**`adjoinIntegers K x` は線形位相**——`0` の近傍基底が
イデアル(`adjoinIntegersBall`)で与えられることから mathlib の
`IsLinearTopology.mk_of_hasBasis` で直接従う。`PowerSeries.aeval` が
要求する「評価先が線形位相であること」の核心部分——`Ideal.
isLinearTopology`(adic位相・`IsDiscreteValuationRing`が要る)も
`Valued`も一切経由しない、より直接的な経路。 -/
theorem isLinearTopology_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    IsLinearTopology (adjoinIntegers K x) (adjoinIntegers K x) :=
  IsLinearTopology.mk_of_hasBasis (adjoinIntegers K x) (hasBasis_nhds_adjoinIntegers K x)

/-! ### `𝒪[K.carrier]` から `adjoinIntegers K x` への埋め込み -/

/-- `𝒪[K.carrier]` から `adjoinIntegers K x` への環準同型——
`algebraMap K.carrier K.closure` の制限。`algebraMap_mem_
adjoinIntegers` が値の所属を保証する。 -/
noncomputable def adjoinIntegersAlgebraMap {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    𝒪[K.carrier] →+* adjoinIntegers K x where
  toFun y := ⟨⟨algebraMap K.carrier K.closure (y : K.carrier),
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).algebraMap_mem _⟩,
    algebraMap_mem_adjoinIntegers K x y⟩
  map_one' := by
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure ((1 : 𝒪[K.carrier]) : K.carrier) = 1
    simp
  map_mul' := by
    intro a b
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure ((a * b : 𝒪[K.carrier]) : K.carrier) =
      algebraMap K.carrier K.closure (a : K.carrier) * algebraMap K.carrier K.closure (b : K.carrier)
    rw [← map_mul]
    norm_cast
  map_zero' := by
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure ((0 : 𝒪[K.carrier]) : K.carrier) = 0
    simp
  map_add' := by
    intro a b
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure ((a + b : 𝒪[K.carrier]) : K.carrier) =
      algebraMap K.carrier K.closure (a : K.carrier) + algebraMap K.carrier K.closure (b : K.carrier)
    rw [← map_add]
    norm_cast

/-- `adjoinIntegers K x` は `𝒪[K.carrier]`-代数——
`adjoinIntegersAlgebraMap` から。★これで `𝒪[K.carrier]` が
`𝒪[K.carrier⟮x⟯]` へ埋め込まれる、という何段階も前から見通していた
課題が解決した。 -/
noncomputable instance adjoinIntegersAlgebra {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    Algebra (𝒪[K.carrier]) (adjoinIntegers K x) :=
  (adjoinIntegersAlgebraMap K x).toAlgebra

/-- `adjoinIntegersAlgebraMap` は連続——`algebraMap K.carrier
K.closure` 自体が連続(`NormedAlgebra`から)であることの制限。 -/
theorem continuous_algebraMap_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    Continuous (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) := by
  apply Continuous.subtype_mk
  apply Continuous.subtype_mk
  exact (continuous_algebraMap K.carrier K.closure).comp continuous_subtype_val

/-- `adjoinIntegers K x` へのスカラー倍は連続——`algebraMap` の連続性
(`continuous_algebraMap_adjoinIntegers`)から mathlib の
`continuousSMul_of_algebraMap` で従う。`PowerSeries.aeval` が要求する
最後の条件。 -/
theorem continuousSMul_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    ContinuousSMul (𝒪[K.carrier]) (adjoinIntegers K x) :=
  continuousSMul_of_algebraMap _ _ (continuous_algebraMap_adjoinIntegers K x)

/-! ### ★★★★★★★★★★★★★節目——`Λ_n` の元で冪級数を実際に評価する -/

/-- `Λ_n` の元 `x`(`adjoinIntegers K x` の元として)は位相的冪零——
`spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints` から
直接。 -/
theorem hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints {p : ℕ} [Fact p.Prime]
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
    PowerSeries.HasEval
      (⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
          K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem⟩ : adjoinIntegers K x) := by
  apply tendsto_pow_atTop_nhds_zero_of_norm_lt_one
  show spectralNorm K.carrier K.closure x < 1
  exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx

/-- ★★★★★★★★★★★★★**節目——冪級数 `f` を `Λ_n` の元 `x` で実際に
評価する**。`PowerSeries.aeval` が要求する全条件(`CompleteSpace`・
`IsLinearTopology`・`ContinuousSMul`、他は自動)が `adjoinIntegers K x`
について揃うことを、このファイルと `LubinTateDistinguishedSeparable.
lean` で積み上げた事実から組み立てるだけ。`f` に限らず任意の
`PowerSeries (𝒪[K.carrier])`(特に `[a]_f` を表す `iteratedLubinTate`
系列)を `x` へ適用できる、`𝒪[K.carrier]`-代数準同型そのもの。 -/
noncomputable def lubinTateEvalAtTorsionPoint {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    PowerSeries (𝒪[K.carrier]) →ₐ[𝒪[K.carrier]] adjoinIntegers K x :=
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  PowerSeries.aeval (R := 𝒪[K.carrier])
    (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem)

/-- ★★★★★★★★★★★★★**Lubin-Tate の `𝒪_K` 作用そのもの——`a·x := [a]_f(x)`
を実際の元として計算する**。`LubinTateAction`(`Found/PGC/
LubinTateActionEndomorphism.lean`、`[a]_f` を表す形式冪級数、既に
`sorry` 無しで確立済み)を、`lubinTateEvalAtTorsionPoint` で `Λ_n` の
元 `x` へ評価するだけ。原典 [pGC] の Lubin-Tate 理論が `𝒪_K` 加群
`Λ_n` を定義する、まさにその作用の**実装**。

残る仕事(節目(3)の完成に向けて): この作用が実際に `Λ_∞:=∪Λ_n` に
戻ること(`a·x ∈ Λ_∞`)、加法性(`(a+b)·x = F_f(a·x, b·x)`)・乗法性
(`a·(b·x)=(ab)·x`)——これらは `LubinTateAction` 自身の関数等式
(`LubinTateAction_functional_equation`)・準同型性(既に確立済みの
`LubinTateAction_add`・`LubinTateAction_mul` 系)を、`PowerSeries.subst`
(形式的代入)ベースの等式から `PowerSeries.aeval`(位相的評価)を
経由した等式へ橋渡しする必要があり、次の一歩として残っている。 -/
noncomputable def lubinTateActionAtTorsionPoint {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) : adjoinIntegers K x :=
  lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
    (LubinTateAction hq hπmax f hf0 hf1 hf a)

/-- `PowerSeries.aeval` は `X` を評価点そのものへ送る——`X`を
`Polynomial.X`の冪級数への係数込みの像として書き直し、
`PowerSeries.aeval_coe`・`Polynomial.aeval_X`と組み合わせるだけの、
`PowerSeries.aeval`の最も基本的な性質(汎用、`AdjoinIntegers`
固有ではない)。 -/
theorem aeval_X_eq_self {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    [UniformSpace R] [UniformSpace S] [IsUniformAddGroup R] [IsTopologicalSemiring R]
    [IsUniformAddGroup S] [T2Space S] [CompleteSpace S] [IsTopologicalRing S]
    [IsLinearTopology S S] [ContinuousSMul R S]
    {a : S} (ha : PowerSeries.HasEval a) :
    PowerSeries.aeval ha (PowerSeries.X : PowerSeries R) = a := by
  have h1 : (PowerSeries.X : PowerSeries R) = ((Polynomial.X : Polynomial R) : PowerSeries R) := by
    simp
  rw [h1, PowerSeries.aeval_coe]
  simp

/-- ★**`1·x = x`**——Lubin-Tate の `𝒪_K` 作用の単位律。`LubinTateAction_
one_eq_X`(`[1]_f = X`、既に確立済み)と`aeval_X_eq_self`(`aeval`は
`X`を評価点そのものへ送る)を組み合わせるだけ。加群作用の公理の
1つが実際に成り立つことを確認した、最初の具体例。 -/
theorem lubinTateActionAtTorsionPoint_one {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem 1 =
      ⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
          K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem⟩ := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  show lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
    (LubinTateAction hq hπmax f hf0 hf1 hf 1) = _
  rw [LubinTateAction_one_eq_X hq hπmax hπne0 f hf0 hf1 hf]
  unfold lubinTateEvalAtTorsionPoint
  apply aeval_X_eq_self

/-- ★**`0·x = 0`**——`LubinTateAction_zero_eq_zero`([0]_f=0、
`Found/PGC/LubinTateActionPiPow.lean`)と`PowerSeries.aeval`が
`AlgHom`(したがって`0`を`0`へ送る)であることを組み合わせるだけ。
加群作用のもう1つの公理(零元の吸収律)を確認した。 -/
theorem lubinTateActionAtTorsionPoint_zero {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem 0 = 0 := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  show lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
    (LubinTateAction hq hπmax f hf0 hf1 hf 0) = 0
  rw [LubinTateAction_zero_eq_zero hq hπmax hπne0 f hf0 hf1 hf]
  unfold lubinTateEvalAtTorsionPoint
  exact map_zero _

/-- ★★★★★★★★★★**作用はノルムを増やさない: `‖a·x‖ ≤ ‖x‖`**——
`[a]_f` の定数項が `0`(`constantCoeff_LubinTateAction`)であることから
`[a]_f = X * h` と分解し(`PowerSeries.X_dvd_iff`)、`PowerSeries.aeval`
が(積を保つ)環準同型であることから
`a·x = aeval x (X*h) = x * (aeval x h)`。`aeval x h` は
`adjoinIntegers K x` の元(ノルム`≤1`)なので `‖a·x‖=‖x‖*‖aeval x h‖
≤‖x‖*1=‖x‖`。`subst`/`aeval` の一般的な橋渡し理論を一切経由しない、
より直接的な経路。 -/
theorem norm_lubinTateActionAtTorsionPoint_le {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    ‖(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ ‖x‖ := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  obtain ⟨h, hh⟩ := PowerSeries.X_dvd_iff.mpr (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)
  show ‖(↑(lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
      (LubinTateAction hq hπmax f hf0 hf1 hf a)) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ ‖x‖
  rw [hh, map_mul]
  unfold lubinTateEvalAtTorsionPoint
  rw [aeval_X_eq_self]
  set y := PowerSeries.aeval (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem) h with hy_def
  set xmem := (⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem⟩ : adjoinIntegers K x) with hxmem_def
  show ‖(↑(xmem * y) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ ‖x‖
  rw [Subring.coe_mul]
  calc ‖(↑xmem : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) * (↑y : _)‖
      = ‖(↑xmem : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ * ‖(↑y : _)‖ :=
        norm_mul _ _
    _ ≤ ‖x‖ * 1 := by
        apply mul_le_mul_of_nonneg_left _ (norm_nonneg _)
        exact y.2
    _ = ‖x‖ := mul_one _

/-- `‖a·x‖ < 1`——`norm_lubinTateActionAtTorsionPoint_le`(`‖a·x‖≤‖x‖`)
と`x`自身が位相的冪零(`spectralNorm_lt_one_of_mem_
iteratedLubinTateTorsionPoints`、`‖x‖<1`)であることを組み合わせる
だけ。これで`a·x`自身も`PowerSeries.HasEval`を満たす見込みが立ち、
作用を反復適用できる土台が整う。 -/
theorem norm_lubinTateActionAtTorsionPoint_lt_one {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    ‖(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1 :=
  lt_of_le_of_lt (norm_lubinTateActionAtTorsionPoint_le K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a)
    (by show spectralNorm K.carrier K.closure x < 1
        exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hx)

/-! ### ★★★★★★★★★★★★★★★★`subst`(形式的代入)と`aeval`(位相的評価)を
繋ぐ連鎖律——mathlib に無かった橋渡しを直接構成する

これまで「一般の`a,b`についての加法性・乗法性には`subst`/`aeval`の
深い橋渡し理論が要る」と記録し、mathlib唯一の`subst`連続性補題
(`MvPowerSeries.continuous_subst`)が`DiscreteUniformity`を要求する
ため直接使えないと判明していた。ここでは**トランケーション(次数
ごとの多項式近似)を経由する別の経路**でこの橋渡しを直接構成する:

1. `coeff_pow_eq_zero_of_lt`: `g`の定数項が`0`なら`g=X*h`と分解でき、
   `g^d`の次数`e`未満の係数は`d>e`のとき`0`。
2. `coeff_subst_trunc_eq`: `subst g p`の次数`e`の係数は、`p`の次数
   `e`より大きい部分を切り捨てても変わらない(`coeff_subst'`の
   `finsum`公式が`d≤e`の項だけで決まることから、1.を使って)。
3. `tendsto_subst_trunc`: 2.から、`subst g (trunc N p)`は`N→∞`で
   `subst g p`へ収束する(`PowerSeries.WithPiTopology.tendsto_iff_
   coeff_tendsto`+各次数で最終的に一定という議論)。
4. `aeval_subst_eq_aeval_aeval`(★本体): `PowerSeries.continuous_
   aeval`(既存)を3.に適用して極限を`aeval`の中へ運び、`Polynomial.
   aeval_algHom_apply`(代数準同型は多項式評価と可換)で有限段
   (トランケーション)での等式を得て、もう一度極限に戻す
   (`PowerSeries.WithPiTopology.tendsto_trunc_atTop`)ことで結論する。

`subst`自身の連続性は一切使わない——`subst`の**係数ごとの有限性**
(`finsum`が有限和に潰れること)だけを使う、より直接的な経路。 -/

/-- `g`の定数項が`0`のとき、`g^d`の次数`e`未満(`e<d`)の係数は`0`。 -/
theorem coeff_pow_eq_zero_of_lt {A : Type*} [CommRing A] (g : PowerSeries A)
    (hg0 : PowerSeries.constantCoeff g = 0) (d e : ℕ) (hde : e < d) :
    PowerSeries.coeff e (g ^ d) = 0 := by
  obtain ⟨h, hh⟩ := PowerSeries.X_dvd_iff.mpr hg0
  rw [hh, mul_pow, PowerSeries.coeff_X_pow_mul']
  simp [Nat.not_le.mpr hde]

/-- `subst g p`の次数`e`の係数は、`p`を次数`N>e`でトランケートしても
変わらない——`coeff_subst'`の`finsum`公式(`∑_d coeff d p • coeff e
(g^d)`)が`coeff_pow_eq_zero_of_lt`により`d≤e`の項だけで決まる
ことから。 -/
theorem coeff_subst_trunc_eq {A : Type*} [CommRing A] {g : PowerSeries A}
    (hg : PowerSeries.HasSubst g) (hg0 : PowerSeries.constantCoeff g = 0)
    (p : PowerSeries A) (e N : ℕ) (hN : e < N) :
    PowerSeries.coeff e (PowerSeries.subst g ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A)) =
      PowerSeries.coeff e (PowerSeries.subst g p) := by
  rw [PowerSeries.coeff_subst' hg, PowerSeries.coeff_subst' hg]
  apply finsum_congr
  intro d
  by_cases hd : d < N
  · congr 1
    rw [Polynomial.coeff_coe, PowerSeries.coeff_trunc]
    simp [hd]
  · rw [show PowerSeries.coeff d ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A) = 0 from by
      rw [Polynomial.coeff_coe, PowerSeries.coeff_trunc]; simp [Nat.not_lt.mp hd], zero_smul]
    have hdN : N ≤ d := Nat.not_lt.mp hd
    have hed : e < d := lt_of_lt_of_le hN hdN
    rw [coeff_pow_eq_zero_of_lt g hg0 d e hed, smul_zero]

open scoped PowerSeries.WithPiTopology in
/-- `subst g (trunc N p)`は`N→∞`で`subst g p`へ収束する
(次数ごとの標準位相について)——`coeff_subst_trunc_eq`により各次数
`e`で`N>e`以降は値が一定になることから。 -/
theorem tendsto_subst_trunc {A : Type*} [CommRing A] [TopologicalSpace A] {g : PowerSeries A}
    (hg : PowerSeries.HasSubst g) (hg0 : PowerSeries.constantCoeff g = 0) (p : PowerSeries A) :
    Filter.Tendsto (fun N => PowerSeries.subst g ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A))
      Filter.atTop (nhds (PowerSeries.subst g p)) := by
  rw [PowerSeries.WithPiTopology.tendsto_iff_coeff_tendsto]
  intro e
  apply tendsto_atTop_of_eventually_const (i₀ := e + 1)
  intro N hN
  exact coeff_subst_trunc_eq hg hg0 p e N (by omega)

open scoped PowerSeries.WithPiTopology in
/-- ★★★★★★★★★★★★★**`subst`と`aeval`を繋ぐ連鎖律**——
`aeval x (subst g p) = aeval (aeval x g) p`。`subst`自身の連続性
(mathlibには`DiscreteUniformity`版しか無く直接使えない)を経由せず、
`tendsto_subst_trunc`(トランケーション経由の収束)・
`PowerSeries.continuous_aeval`(既存)・`Polynomial.aeval_algHom_
apply`(代数準同型は多項式評価と可換)だけを組み合わせて構成する。
Lubin-Tate作用の加法性・乗法性(`(a+b)·x=F_f(a·x,b·x)`・
`a·(b·x)=(ab)·x`)を`LubinTateAction_add`・`LubinTateAction_comp`
(既存)と組み合わせて示すための、探していた本物の橋渡し。 -/
theorem aeval_subst_eq_aeval_aeval {A S : Type*} [CommRing A] [CommRing S] [Algebra A S]
    [UniformSpace A] [UniformSpace S] [IsUniformAddGroup A] [IsTopologicalSemiring A]
    [IsUniformAddGroup S] [T2Space S] [CompleteSpace S] [IsTopologicalRing S]
    [IsLinearTopology S S] [ContinuousSMul A S]
    {g p : PowerSeries A} (hg : PowerSeries.HasSubst g) (hg0 : PowerSeries.constantCoeff g = 0)
    {x : S} (hx : PowerSeries.HasEval x)
    {y : S} (hy : PowerSeries.HasEval y) (hxy : PowerSeries.aeval hx g = y) :
    PowerSeries.aeval hx (PowerSeries.subst g p) = PowerSeries.aeval hy p := by
  have h1 : Filter.Tendsto (fun N => PowerSeries.aeval hx
      (PowerSeries.subst g ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A)))
      Filter.atTop (nhds (PowerSeries.aeval hx (PowerSeries.subst g p))) :=
    (PowerSeries.continuous_aeval hx).continuousAt.tendsto.comp (tendsto_subst_trunc hg hg0 p)
  have h2 : ∀ N, PowerSeries.aeval hx
      (PowerSeries.subst g ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A)) =
      PowerSeries.aeval hy ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A) := by
    intro N
    rw [PowerSeries.subst_coe hg, PowerSeries.aeval_coe, ← hxy]
    exact (Polynomial.aeval_algHom_apply (PowerSeries.aeval hx) g (PowerSeries.trunc N p)).symm
  rw [funext h2] at h1
  have h3 : Filter.Tendsto (fun N => PowerSeries.aeval hy
      ((PowerSeries.trunc N p : Polynomial A) : PowerSeries A))
      Filter.atTop (nhds (PowerSeries.aeval hy p)) :=
    (PowerSeries.continuous_aeval hy).continuousAt.tendsto.comp
      (PowerSeries.WithPiTopology.tendsto_trunc_atTop (R := A) p)
  exact tendsto_nhds_unique h1 h3

/-! ### ★★★★★★★★★★★★★★★★節目——Lubin-Tate 作用の乗法性
`a·(b·x)=(ab)·x`

`aeval_subst_eq_aeval_aeval`(連鎖律)と既存の`LubinTateAction_comp`
(`[ab]_f=subst([b]_f)([a]_f)`)を組み合わせて、`𝒪_K` 加群作用の
乗法性を確立する。`b·x`(`adjoinIntegers K x` の元)を新たな評価点
として、同じ`adjoinIntegers K x`(座標変換無し、`CompleteSpace`・
`IsLinearTopology`・`ContinuousSMul` はすべて `x` にのみ依存し
`b·x` には依存しない)へ`aeval`を組み立て直すのが鍵。 -/

/-- `adjoinIntegers K x` の元 `z` の位相的冪零性は、周囲の体での
位相的冪零性と同値——`adjoinIntegers K x` の位相は周囲の体からの
誘導位相(部分環)なので、`tendsto_subtype_rng`(部分型での収束は
包含写像を通した収束と同値、mathlib一般論)から直ちに従う。 -/
theorem hasEval_iff_coe {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (z : adjoinIntegers K x) :
    PowerSeries.HasEval z ↔
      PowerSeries.HasEval (↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  show Filter.Tendsto (fun n => z ^ n) Filter.atTop (nhds 0) ↔
    Filter.Tendsto (fun n => (↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ^ n)
      Filter.atTop (nhds 0)
  rw [tendsto_subtype_rng]
  have heq : (fun n => (↑(z ^ n) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))) =
      (fun n => (↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ^ n) := by
    funext n; exact SubmonoidClass.coe_pow z n
  rw [heq]
  simp

/-- `b·x`(`lubinTateActionAtTorsionPoint`の値)は位相的冪零——
`norm_lubinTateActionAtTorsionPoint_lt_one`と`hasEval_iff_coe`から
直ちに従う。これで`b·x`自身を評価点として`[a]_f`を再評価できる。 -/
theorem hasEval_lubinTateActionAtTorsionPoint {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (b : 𝒪[K.carrier]) :
    PowerSeries.HasEval (lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b) := by
  rw [hasEval_iff_coe]
  apply tendsto_pow_atTop_nhds_zero_of_norm_lt_one
  exact norm_lubinTateActionAtTorsionPoint_lt_one K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b

/-- `adjoinIntegers K x` の**任意の**位相的冪零な元 `z` を評価点として
冪級数を評価する——`x` 自身に特化していた`lubinTateEvalAtTorsionPoint`
の一般化。`CompleteSpace`・`IsLinearTopology`・`ContinuousSMul` は
`x`(座標系)のみに依存し `z` には依存しないので、そのまま流用できる。 -/
noncomputable def lubinTateEvalAtPoint {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z) :
    PowerSeries (𝒪[K.carrier]) →ₐ[𝒪[K.carrier]] adjoinIntegers K x :=
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  PowerSeries.aeval hz

/-- ★★★★★★★★★★★★★★★★**節目——Lubin-Tate 作用の乗法性**:
`a·(b·x) = (ab)·x`。既存の`LubinTateAction_comp`
(`[ab]_f=subst([b]_f)([a]_f)`)の両辺を`x`で評価し、
`aeval_subst_eq_aeval_aeval`(連鎖律)で右辺を
`aeval(aeval x [b]_f)([a]_f) = a·(b·x)` へ変形するだけ。`𝒪_K`
加群の構造公理のうち最も本質的な1つ(乗法性、`(ab)·x=a·(b·x)`)を
初めて確立した。 -/
theorem lubinTateAction_mul {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a b : 𝒪[K.carrier]) :
    lubinTateEvalAtPoint K x (lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b)
        (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b)
        (LubinTateAction hq hπmax f hf0 hf1 hf a) =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem (a * b) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  show PowerSeries.aeval (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b)
      (LubinTateAction hq hπmax f hf0 hf1 hf a) = _
  show PowerSeries.aeval (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b)
      (LubinTateAction hq hπmax f hf0 hf1 hf a) =
    lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
      (LubinTateAction hq hπmax f hf0 hf1 hf (a * b))
  rw [LubinTateAction_comp hq hπmax hπne0 f hf0 hf1 hf a b]
  unfold lubinTateEvalAtTorsionPoint
  symm
  apply aeval_subst_eq_aeval_aeval
  · show IsNilpotent (PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf b))
    rw [constantCoeff_LubinTateAction]; exact IsNilpotent.zero
  · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf b
  · rfl

/-! ### ★★★★★★★★★★★★★★★★★★★★節目——Lubin-Tate 作用の加法性
`(a+b)·x = F_f(a·x,b·x)`

`aeval_subst_eq_aeval_aeval`(1変数の連鎖律)の**2変数版**を同じ
トランケーション経由の手法で構成し、既存の`LubinTateAction_add`
(`[a+b]_f=subst(family)(F_f)`)と組み合わせて加法性を確立する。 -/

/-- `family i = X * h_i` と分解できるとき、`d.prod` の次数`e`未満の
係数は`0`(`coeff_pow_eq_zero_of_lt`の2変数版)——各成分を`X`で
くくり出し、積全体が`X^|d|`で割り切れることから。 -/
theorem coeff_prod_pow_eq_zero_of_lt {A : Type*} [CommRing A] (family : Fin 2 → PowerSeries A)
    (hfamily : ∀ i, PowerSeries.constantCoeff (family i) = 0)
    (d : Fin 2 →₀ ℕ) (e : ℕ) (hde : e < d.sum (fun _ ee => ee)) :
    MvPowerSeries.coeff (Finsupp.single () e) (d.prod (fun s ee => (family s) ^ ee)) = 0 := by
  show PowerSeries.coeff e (d.prod (fun s ee => (family s) ^ ee)) = 0
  have hfact : ∀ i, ∃ h : PowerSeries A, family i = PowerSeries.X * h :=
    fun i => PowerSeries.X_dvd_iff.mpr (hfamily i)
  choose h hh using hfact
  have heq : d.prod (fun s ee => (family s) ^ ee) =
      PowerSeries.X ^ (d.sum (fun _ ee => ee)) * d.prod (fun s ee => (h s) ^ ee) := by
    rw [Finsupp.prod, Finsupp.prod, Finsupp.sum]
    rw [← Finset.prod_pow_eq_pow_sum]
    rw [← Finset.prod_mul_distrib]
    apply Finset.prod_congr rfl
    intro s _
    rw [hh s, mul_pow]
  rw [heq, PowerSeries.coeff_X_pow_mul']
  simp [Nat.not_le.mpr hde]

/-- `subst family Φ`の次数`e`の係数は、`Φ`を多重次数`n`でトランケート
(`trunc'`)しても変わらない——`n`が`d`のすべての成分を覆っているか、
`d`の総次数が`e`を超えるかのいずれかであれば十分
(`coeff_prod_pow_eq_zero_of_lt`と`MvPowerSeries.coeff_subst`の
`finsum`公式から)。 -/
theorem coeff_subst_family_trunc_eq {A : Type*} [CommRing A] {family : Fin 2 → PowerSeries A}
    (hfamily : MvPowerSeries.HasSubst family) (hfamily0 : ∀ i, PowerSeries.constantCoeff (family i) = 0)
    (Φ : MvPowerSeries (Fin 2) A) (e : ℕ) (n : Fin 2 →₀ ℕ)
    (hn : ∀ d : Fin 2 →₀ ℕ, e < d.sum (fun _ ee => ee) ∨ d ≤ n) :
    PowerSeries.coeff e (MvPowerSeries.subst family
        ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) : MvPowerSeries (Fin 2) A)) =
      PowerSeries.coeff e (MvPowerSeries.subst family Φ) := by
  show MvPowerSeries.coeff (Finsupp.single () e) _ = MvPowerSeries.coeff (Finsupp.single () e) _
  rw [MvPowerSeries.coeff_subst hfamily, MvPowerSeries.coeff_subst hfamily]
  apply finsum_congr
  intro d
  rcases hn d with hd | hd
  · rw [coeff_prod_pow_eq_zero_of_lt family hfamily0 d e hd, smul_zero, smul_zero]
  · have hcoeff : MvPowerSeries.coeff d ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) :
        MvPowerSeries (Fin 2) A) = MvPowerSeries.coeff d Φ := by
      show (MvPowerSeries.trunc' A n Φ).coeff d = MvPowerSeries.coeff d Φ
      rw [MvPowerSeries.coeff_trunc']
      simp [hd]
    rw [hcoeff]

open scoped MvPowerSeries.WithPiTopology in
/-- `subst family(trunc' A n Φ)`は`n→∞`で`subst family Φ`へ収束する
(1変数版`tendsto_subst_trunc`の2変数版)——各出力次数`e`について、
`n:=`(定数`e`)を境に`coeff_subst_family_trunc_eq`が値の一定性を
保証することから。 -/
theorem tendsto_subst_family_trunc {A : Type*} [CommRing A] [TopologicalSpace A]
    {family : Fin 2 → PowerSeries A}
    (hfamily : MvPowerSeries.HasSubst family) (hfamily0 : ∀ i, PowerSeries.constantCoeff (family i) = 0)
    (Φ : MvPowerSeries (Fin 2) A) :
    Filter.Tendsto (fun n => MvPowerSeries.subst family
        ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) : MvPowerSeries (Fin 2) A))
      Filter.atTop (nhds (MvPowerSeries.subst family Φ)) := by
  rw [PowerSeries.WithPiTopology.tendsto_iff_coeff_tendsto]
  intro e
  apply tendsto_atTop_of_eventually_const (i₀ := Finsupp.equivFunOnFinite.symm (fun _ : Fin 2 => e))
  intro n hn
  apply coeff_subst_family_trunc_eq hfamily hfamily0 Φ e n
  intro d
  have hsum : d.sum (fun _ ee => ee) = d 0 + d 1 := by
    rw [Finsupp.sum_fintype]
    · rw [Fin.sum_univ_two]
    · intro; rfl
  have hn0 : e ≤ n 0 := by have := hn 0; simpa using this
  have hn1 : e ≤ n 1 := by have := hn 1; simpa using this
  by_cases hd0 : d 0 ≤ e
  · by_cases hd1 : d 1 ≤ e
    · right
      intro i
      fin_cases i
      · simpa using le_trans hd0 hn0
      · simpa using le_trans hd1 hn1
    · left; omega
  · left; omega

open scoped MvPowerSeries.WithPiTopology in
/-- ★★★★★★★★★★★★★★★★**`subst`と`aeval`を繋ぐ連鎖律の2変数版**——
`aeval x(subst family Φ) = aeval y Φ`(`y i := aeval x(family i)`)。
`aeval_subst_eq_aeval_aeval`と全く同じ手法(トランケーション経由の
収束・`PowerSeries.continuous_aeval`/`MvPowerSeries.continuous_aeval`
・`MvPolynomial.comp_aeval_apply`)を2変数へ拡張しただけ。 -/
theorem aeval_subst_family_eq_aeval_aeval {A S : Type*} [CommRing A] [CommRing S] [Algebra A S]
    [UniformSpace A] [UniformSpace S] [IsUniformAddGroup A] [IsTopologicalSemiring A]
    [IsUniformAddGroup S] [T2Space S] [CompleteSpace S] [IsTopologicalRing S]
    [IsLinearTopology S S] [ContinuousSMul A S]
    {family : Fin 2 → PowerSeries A} {Φ : MvPowerSeries (Fin 2) A}
    (hfamily : MvPowerSeries.HasSubst family) (hfamily0 : ∀ i, PowerSeries.constantCoeff (family i) = 0)
    {x : S} (hx : PowerSeries.HasEval x)
    {y : Fin 2 → S} (hy : MvPowerSeries.HasEval y)
    (hxy : ∀ i, PowerSeries.aeval hx (family i) = y i) :
    PowerSeries.aeval hx (MvPowerSeries.subst family Φ) = MvPowerSeries.aeval hy Φ := by
  have h1 : Filter.Tendsto (fun n => PowerSeries.aeval hx
      (MvPowerSeries.subst family ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) :
        MvPowerSeries (Fin 2) A)))
      Filter.atTop (nhds (PowerSeries.aeval hx (MvPowerSeries.subst family Φ))) :=
    (PowerSeries.continuous_aeval hx).continuousAt.tendsto.comp (tendsto_subst_family_trunc hfamily hfamily0 Φ)
  have h2 : ∀ n, PowerSeries.aeval hx
      (MvPowerSeries.subst family ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) :
        MvPowerSeries (Fin 2) A)) =
      MvPowerSeries.aeval hy ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) :
        MvPowerSeries (Fin 2) A) := by
    intro n
    rw [MvPowerSeries.subst_coe, MvPolynomial.comp_aeval_apply, MvPowerSeries.aeval_coe,
      funext hxy]
  rw [funext h2] at h1
  have h3 : Filter.Tendsto (fun n => MvPowerSeries.aeval hy
      ((MvPowerSeries.trunc' A n Φ : MvPolynomial (Fin 2) A) : MvPowerSeries (Fin 2) A))
      Filter.atTop (nhds (MvPowerSeries.aeval hy Φ)) :=
    (MvPowerSeries.continuous_aeval hy).continuousAt.tendsto.comp
      (MvPowerSeries.WithPiTopology.tendsto_trunc'_atTop Φ)
  exact tendsto_nhds_unique h1 h3

/-- `family:=fun i=>if i=0 then a·x else b·x`は位相的冪零な族——
`Fin 2`は有限型なので`cofinite`フィルターが`⊥`になり、
`Tendsto _ ⊥ _`は自動的に成り立つ(`MvPowerSeries.HasEval`のもう
片方の条件)。各成分の位相的冪零性は`hasEval_lubinTateActionAtTorsionPoint`
から。 -/
theorem hasEval_actionFam2 {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a b : 𝒪[K.carrier]) :
    MvPowerSeries.HasEval (fun i : Fin 2 =>
      if i = 0 then (lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a : adjoinIntegers K x)
      else (lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b : adjoinIntegers K x)) := by
  constructor
  · intro i
    split
    · exact hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a
    · exact hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b
  · rw [show (Filter.cofinite : Filter (Fin 2)) = ⊥ from by rw [Filter.cofinite_eq_bot_iff]; infer_instance]
    exact Filter.tendsto_bot

/-- `adjoinIntegers K x`の中で2変数の冪級数(特に`formalGroupLaw`)を
位相的冪零な点`y : Fin 2 → adjoinIntegers K x`で評価する——
`lubinTateEvalAtPoint`の2変数版、同じ`CompleteSpace`等をそのまま
流用する。 -/
noncomputable def lubinTateEvalFormalGroupAt {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (y : Fin 2 → adjoinIntegers K x) (hy : MvPowerSeries.HasEval y)
    (Φ : MvPowerSeries (Fin 2) (𝒪[K.carrier])) : adjoinIntegers K x :=
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  MvPowerSeries.aeval hy Φ

/-- ★★★★★★★★★★★★★★★★★★★★**節目——Lubin-Tate 作用の加法性**:
`(a+b)·x = F_f(a·x,b·x)`。既存の`LubinTateAction_add`
(`[a+b]_f=subst(family)(F_f)`、`family i:=if i=0 then[a]_f else[b]_f`)
の両辺を`x`で評価し、`aeval_subst_family_eq_aeval_aeval`(2変数の
連鎖律)で右辺を`aeval(fun i=>aeval x(family i))(F_f)
=F_f(a·x,b·x)`へ変形するだけ。これで`𝒪_K`加群の構造公理——単位律・
零元の吸収律・乗法性・加法性——が**すべて**確立された。 -/
theorem lubinTateAction_add {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a b : 𝒪[K.carrier]) :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem (a + b) =
      lubinTateEvalFormalGroupAt K x
        (fun i => if i = 0 then lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a
          else lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem b)
        (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a b)
        (formalGroupLaw hq hπmax f hf0 hf1 hf) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  show lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
    (LubinTateAction hq hπmax f hf0 hf1 hf (a + b)) = _
  unfold lubinTateEvalAtTorsionPoint lubinTateEvalFormalGroupAt
  rw [LubinTateAction_add hq hπmax hπne0 f hf0 hf1 hf a b]
  apply aeval_subst_family_eq_aeval_aeval (hasSubst_actionFam2 hq hπmax f hf0 hf1 hf a b)
  · intro i
    split
    · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
    · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf b
  · intro i
    split
    · rfl
    · rfl

/-! ### ★★★★★★★★★★★★★★★★★★★★★★★★節目——`Λ_n` は `𝒪_K` 作用で閉じている

古典的な Lubin-Tate 理論の中核事実——`Λ_n` が実際に `𝒪_K`-加群である
こと(`a·x ∈ Λ_n` for `x ∈ Λ_n`, `a ∈ 𝒪_K`)——を確立する。鍵となる
2段階: (1) `x∈Λ_n ⟹ π^n·x=0`(`D_n∣[π^n]_f` から)、(2) `π^n·x=0`
かつ乗法性・可換性から `π^n·(a·x)=a·(π^n·x)=a·0=0`、(3) 再び
`[π^n]_f=D_n·U_n`(`U_n` は単位)の分解で `π^n·(a·x)=0` から
`D_n(a·x)=0` を復元する——`D_n(x)=0⟹[π^n]_f(x)=0`(乗法のみ)は
簡単だが、逆向き`[π^n]_f(a·x)=0⟹D_n(a·x)=0`には`U_n(a·x)`が単位
であることを要する。 -/

/-- `x∈Λ_n`ならば`D_n(x)=0`——`iteratedLubinTateTorsionPoints`の定義
そのもの(`Multiset.mem_toFinset`+`Polynomial.mem_roots'`)を、
`adjoinIntegers K x`から`K.closure`への環準同型`g`(単射)を通して
`Polynomial.aeval`のレベルへ引き戻すだけ。 -/
theorem aeval_iteratedLubinTateDistinguished_eq_zero {p : ℕ} [Fact p.Prime]
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
    Polynomial.aeval (⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem⟩ : adjoinIntegers K x)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 := by
  set x' : adjoinIntegers K x := ⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem⟩ with hx'_def
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hx
  obtain ⟨_, hxroot⟩ := hx
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have hgx' : g x' = x := rfl
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext
    intro y
    rfl
  have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g x'
  rw [← Polynomial.aeval_def, hgcomp, hgx'] at key
  have haevalx : Polynomial.eval₂ (algebraMap (𝒪[K.carrier]) K.closure) x
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 := by
    rw [Polynomial.eval₂_eq_eval_map]
    exact hxroot
  rw [haevalx] at key
  have hginj : Function.Injective g := fun a b h => Subtype.ext (Subtype.ext h)
  apply hginj
  rw [key, map_zero]

/-- `x∈Λ_n`ならば`π^n·x=0`——`LubinTateAction_pi_pow`(`[π^n]_f=`
`n`回自己合成)・`iteratedLubinTate_eq_distinguished_mul_unit`
(`[π^n]_f=D_n*U_n`)・`aeval_iteratedLubinTateDistinguished_eq_zero`
(`D_n(x)=0`)を組み合わせるだけ。 -/
theorem pi_pow_action_eq_zero {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem (π ^ n) = 0 := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  show lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem
    (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)) = 0
  rw [LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf n,
    iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf n]
  unfold lubinTateEvalAtTorsionPoint
  rw [map_mul, PowerSeries.aeval_coe]
  rw [aeval_iteratedLubinTateDistinguished_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem]
  rw [zero_mul]

/-- `g`の定数項が0なら、`adjoinIntegers K x`の`0`での評価は`0`——
`X`分解トリック(`aeval_X_eq_self`と`map_mul`)だけ。 -/
theorem lubinTateEvalAtPoint_zero {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hz0 : PowerSeries.HasEval (0 : adjoinIntegers K x))
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.constantCoeff g = 0) :
    lubinTateEvalAtPoint K x 0 hz0 g = 0 := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  show PowerSeries.aeval hz0 g = 0
  obtain ⟨h, hh⟩ := PowerSeries.X_dvd_iff.mpr hg0
  rw [hh, map_mul]
  show PowerSeries.aeval hz0 PowerSeries.X * PowerSeries.aeval hz0 h = 0
  rw [aeval_X_eq_self]
  simp

/-- `lubinTateEvalAtPoint_zero`を、評価点が`0`と**等しい**(定義的に
`0`そのものでなくても)場合へ一般化——`subst`で評価点を`0`に
置き換えるだけ(`HasEval`が依存型なので、直接`rw`すると動機が
型として正しくならない、という定型の回避策)。 -/
theorem lubinTateEvalAtPoint_eq_zero_of_eq_zero {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {z : adjoinIntegers K x} (hz : PowerSeries.HasEval z) (hz0 : z = 0)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.constantCoeff g = 0) :
    lubinTateEvalAtPoint K x z hz g = 0 := by
  subst hz0
  exact lubinTateEvalAtPoint_zero K x hz g hg0

/-- `x∈Λ_n`ならば`(a*π^n)·x=0`——`π^n·(a·x)=a·(π^n·x)`
(`lubinTateAction_mul`を2回・可換律)と`π^n·x=0`
(`pi_pow_action_eq_zero`)から。 -/
theorem pi_pow_action_action_eq_zero {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem (a * π ^ n) = 0 := by
  rw [← lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a (π ^ n)]
  exact lubinTateEvalAtPoint_eq_zero_of_eq_zero K x _
    (pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem)
    (LubinTateAction hq hπmax f hf0 hf1 hf a)
    (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)

/-- `[π^n]_f(z)=0`ならば`D_n(z)=0`(`z`は`adjoinIntegers K x`の
**任意**の位相的冪零な元)——`[π^n]_f=D_n*U_n`と`U_n(z)`が単位
(単位を環準同型で送った像は単位)であることから、単位を約分するだけ
(`IsUnit.mul_left_eq_zero`)。`aeval_iteratedLubinTateDistinguished_
eq_zero`の逆向きにあたる。 -/
theorem eq_zero_of_pi_pow_action_eq_zero {p : ℕ} [Fact p.Prime]
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (z : adjoinIntegers K x) (hz : PowerSeries.HasEval z)
    (hz0 : lubinTateEvalAtPoint K x z hz (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)) = 0) :
    Polynomial.aeval z (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  rw [LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf n,
    iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf n] at hz0
  unfold lubinTateEvalAtPoint at hz0
  rw [map_mul, PowerSeries.aeval_coe] at hz0
  rwa [IsUnit.mul_left_eq_zero] at hz0
  exact (isUnit_iteratedLubinTateUnit hq hπmax hπne0 f hf0 hf1 hf n).map (PowerSeries.aeval hz)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**節目——`Λ_n` は `𝒪_K` 作用で閉じて
いる**: `x∈Λ_n`・`a∈𝒪_K`ならば`a·x∈Λ_n`(実際の元として)。
`pi_pow_action_action_eq_zero`(`(a*π^n)·x=0`、可換律で`π^n·(a·x)=0`
と同値)と`eq_zero_of_pi_pow_action_eq_zero`(`[π^n]_f(a·x)=0⟹
D_n(a·x)=0`)を組み合わせ、`D_n(a·x)=0`を`Polynomial.hom_eval₂`で
`K.closure`のレベルへ押し出して`iteratedLubinTateTorsionPoints`の
定義そのものに一致させる。これで古典的な Lubin-Tate 理論の
`Λ_n`が真に`𝒪_K`-加群であるという中核事実が、実在の任意のp進局所体
`K`について`sorry`無しで確立された。 -/
theorem lubinTateActionAtTorsionPoint_mem {p : ℕ} [Fact p.Prime]
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  set z := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a with hz_def
  have hz0 : lubinTateEvalAtPoint K x z
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a)
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ n)) = 0 := by
    rw [lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem (π ^ n) a, mul_comm]
    exact pi_pow_action_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem a
  have hDn0 : Polynomial.aeval z (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n) = 0 :=
    eq_zero_of_pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x z _ hz0
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
  refine ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
    hq hπmax hπne0 f hf0 hf1 hf n).monic.map _ |>.ne_zero, ?_⟩
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
  rw [← Polynomial.aeval_def] at key
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  rw [hgcomp, hDn0, map_zero] at key
  show (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf n)).eval (g z) = 0
  rw [Polynomial.eval_map]
  exact key.symm

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★節目——`a·x=0 ↔ π^n∣a`(核の確定)

古典的な Lubin-Tate 加群の同型 `𝒪_K/π^n≅Λ_n` の**核**にあたる部分——
`x` が「原始的な」π^n-捩れ点(`ψ_n` の根)のとき、写像 `a↦a·x`
(`𝒪_K→Λ_n`)の核がちょうど `π^n𝒪_K` であること——を確立する。
鍵になるのは(1) `π^{n-1}·x≠0`(`x` が `Λ_{n-1}` の元でないこと、
`iteratedLubinTateTorsionPoints_disjoint_iteratedLubinTatePsiTorsionPoints`
から)と(2) 付値環の一意分解(`a=u*π^k`、`u` は単位)を使った単位の
約分——`lubinTateAction_mul`(乗法性)を2回使うだけで、`F_f` 加法の
結合律・逆元は一切不要だった。 -/

/-- ★★★★★★★★★★**「原始的な」π^n-捩れ点は `π^{n-1}` で消えない**——
`x` が `Λ_{n-1}` に属さないこと(`ψ_n` の根と `Λ_{n-1}` は交わらない、
`iteratedLubinTateTorsionPoints_disjoint_iteratedLubinTatePsiTorsionPoints`)
の対偶を、`aeval_iteratedLubinTateDistinguished_eq_zero`と同じ
`Polynomial.hom_eval₂`ベースの橋渡しで`π^{n-1}·x=0⟹D_{n-1}(x)=0
⟹x∈Λ_{n-1}`という形にして使う。 -/
theorem lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ (n - 1)) ≠ 0 := by
  intro h0
  set z : adjoinIntegers K x := ⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem⟩ with hz_def
  have hz : PowerSeries.HasEval z :=
    hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
  have h0' : lubinTateEvalAtPoint K x z hz (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1))) = 0 := h0
  have hkey := eq_zero_of_pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf (n - 1) x z hz h0'
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have hgz : g z = x := rfl
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  have hkey2 := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1))
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
  rw [← Polynomial.aeval_def, hgcomp, hgz, hkey, map_zero] at hkey2
  have heval : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1))).eval x = 0 := by
    rw [Polynomial.eval_map]; exact hkey2.symm
  have hxn1 : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (n - 1) := by
    rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
    exact ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
      hq hπmax hπne0 f hf0 hf1 hf (n - 1)).monic.map _ |>.ne_zero, heval⟩
  exact Finset.disjoint_left.mp (iteratedLubinTateTorsionPoints_disjoint_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn) hxn1 hxψ

/-- `hπmax`(`maximalIdeal=span{π}`)から`π`が既約元であること——
`IsDiscreteValuationRing.exists_irreducible`で得た既約元`ϖ`と
`span{π}=maximalIdeal=span{ϖ}`から`π`・`ϖ`が同伴であること
(`Ideal.span_singleton_eq_span_singleton`)を経由する。 -/
theorem irreducible_of_maximalIdeal_eq_span {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    Irreducible π := by
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have heq : Ideal.span ({π} : Set (𝒪[K.carrier])) = Ideal.span {ϖ} :=
    hπmax ▸ hϖ.maximalIdeal_eq
  rw [Ideal.span_singleton_eq_span_singleton] at heq
  exact heq.symm.irreducible hϖ

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**核の確定(⟹)**: `a·x=0⟹π^n∣a`
(`x`が「原始的な」π^n-捩れ点のとき)。`a≠0`のとき付値環の一意分解
`a=u*π^k`(`u`単位)を取り、`k<n`と仮定して矛盾を導く: `u⁻¹`を
`lubinTateAction_mul`で作用させて`π^k·x=0`を復元し(単位の約分)、
さらに`π^{n-1-k}`を掛けて`π^{n-1}·x=0`を導き、
`lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_
iteratedLubinTatePsiTorsionPoints`に矛盾する。`F_f`加法の逆元・
結合律は一切使わず、乗法性(`lubinTateAction_mul`)のみで閉じる。 -/
theorem lubinTateActionAtTorsionPoint_eq_zero_imp_dvd_of_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier])
    (h0 : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a = 0) :
    π ^ n ∣ a := by
  haveI := valuationRing_isDVR K
  have hπirr : Irreducible π := irreducible_of_maximalIdeal_eq_span K hπmax
  rcases eq_or_ne a 0 with ha0 | hane0
  · exact ha0 ▸ dvd_zero _
  obtain ⟨k, u, hauk⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hane0 hπirr
  by_contra hndvd
  have hklt : k < n := by
    by_contra hge
    apply hndvd
    rw [hauk]
    exact Dvd.dvd.mul_left (pow_dvd_pow π (by omega)) (u : 𝒪[K.carrier])
  have hzero : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ k) = 0 := by
    have hmul := lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (↑u⁻¹ : 𝒪[K.carrier]) a
    rw [lubinTateEvalAtPoint_eq_zero_of_eq_zero K x
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a) h0
      (LubinTateAction hq hπmax f hf0 hf1 hf (↑u⁻¹ : 𝒪[K.carrier]))
      (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (↑u⁻¹ : 𝒪[K.carrier]))] at hmul
    have hcalc : (↑u⁻¹ : 𝒪[K.carrier]) * a = π ^ k := by
      rw [hauk]
      have hassoc : (↑u⁻¹ : 𝒪[K.carrier]) * (↑u * π ^ k) =
          ((↑u⁻¹ : 𝒪[K.carrier]) * ↑u) * π ^ k := by ring
      rw [hassoc, Units.inv_mul, one_mul]
    rw [hcalc] at hmul
    exact hmul.symm
  have hpred : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ (n - 1)) = 0 := by
    have hmul := lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ (n - 1 - k)) (π ^ k)
    rw [lubinTateEvalAtPoint_eq_zero_of_eq_zero K x
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ k)) hzero
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1 - k)))
      (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1 - k)))] at hmul
    rw [← pow_add] at hmul
    have heq : n - 1 - k + k = n - 1 := by omega
    rw [heq] at hmul
    exact hmul.symm
  exact lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hmem hxn hpred

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**核の確定(iff)**:
`a·x=0 ↔ π^n∣a`(`x`が「原始的な」π^n-捩れ点のとき)。⟸は
`pi_pow_action_action_eq_zero`(既出、`(c*π^n)·x=0`)そのもの、⟹は
直前の`lubinTateActionAtTorsionPoint_eq_zero_imp_dvd_of_mem_
iteratedLubinTatePsiTorsionPoints`。古典的な同型`𝒪_K/π^n≅Λ_n`の
核がちょうど`π^n𝒪_K`であることの確定——`|𝒪_K/π^n|=q^n=|Λ_n|`
(既出`card_iteratedLubinTateTorsionPoints`)と組み合わせれば、この
写像の単射性から全単射性が従う見通し(次の一歩)。 -/
theorem lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a = 0 ↔ π ^ n ∣ a := by
  constructor
  · exact lubinTateActionAtTorsionPoint_eq_zero_imp_dvd_of_mem_iteratedLubinTatePsiTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hmem hxn a
  · rintro ⟨c, rfl⟩
    rw [mul_comm]
    exact pi_pow_action_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem c

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**単位は原始的な捩れ点を原始的な
捩れ点へ送る**: `u∈(𝒪_K)^×`・`x`が`ψ_n`の根(原始的なπ^n-捩れ点)
ならば`u·x`も`ψ_n`の根。`(𝒪_K)^×`が`ψ_nの根`の集合に作用することの
確認——`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`へ向けて、`(𝒪_K/π^n)^×`の
軌道が「原始的な捩れ点全体」に収まることの布石。

証明の筋: `u·x∈Λ_n`(既出`lubinTateActionAtTorsionPoint_mem`)は
自動。`u·x∉Λ_{n-1}`は背理法——`u·x∈Λ_{n-1}`と仮定すると
(`aeval_iteratedLubinTateDistinguished_eq_zero`と同じ`Polynomial.
hom_eval₂`橋渡しで)`D_{n-1}(u·x)=0`(`adjoinIntegers K x`のレベル、
`x`自身の座標系のまま)が出て、`[π^{n-1}]_f=D_{n-1}*U_{n-1}`から
`(π^{n-1}*u)·x=0`、`lubinTateAction_mul`で単位`u`を約分(`u⁻¹`を
掛けるだけ、`F_f`加法は不要)して`π^{n-1}·x=0`を得て
`lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_
iteratedLubinTatePsiTorsionPoints`に矛盾する——**新しい点`u·x`自身の
`adjoinIntegers`を一切構築せず、`x`自身の座標系だけで閉じる**のが鍵
(異なる点どうしの`adjoinIntegers`インスタンスを橋渡しする必要が
無かった)。 -/
theorem unit_action_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (u : (𝒪[K.carrier])ˣ) :
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  set z := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (u : 𝒪[K.carrier])
    with hz_def
  have hyn : (↑(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n :=
    lubinTateActionAtTorsionPoint_mem K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (u : 𝒪[K.carrier])
  rw [← iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn, Finset.mem_sdiff]
  refine ⟨hyn, fun hyn1 => ?_⟩
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have hginj : Function.Injective g := fun a b h => Subtype.ext (Subtype.ext h)
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hyn1
  obtain ⟨_, hyroot⟩ := hyn1
  rw [Polynomial.IsRoot, Polynomial.eval_map] at hyroot
  have hgz : g z = ↑(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := rfl
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  have haevalz0 : Polynomial.aeval z (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1)) = 0 := by
    apply hginj
    have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (n - 1))
      (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
    rw [← Polynomial.aeval_def, hgcomp] at key
    rw [key, map_zero, hgz]
    exact hyroot
  have hzero : lubinTateEvalAtPoint K x z
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (u : 𝒪[K.carrier]))
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1))) = 0 := by
    show PowerSeries.aeval _ (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (n - 1))) = 0
    rw [LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf (n - 1),
      iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf (n - 1),
      map_mul, PowerSeries.aeval_coe, haevalz0, zero_mul]
  have hmul := lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ (n - 1)) (u : 𝒪[K.carrier])
  have hprod0 : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
      (π ^ (n - 1) * (u : 𝒪[K.carrier])) = 0 := hmul.symm.trans hzero
  have hkey : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (π ^ (n - 1)) = 0 := by
    have hmul2 := lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
      (↑u⁻¹ : 𝒪[K.carrier]) (π ^ (n - 1) * (u : 𝒪[K.carrier]))
    rw [lubinTateEvalAtPoint_eq_zero_of_eq_zero K x
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
        (π ^ (n - 1) * (u : 𝒪[K.carrier]))) hprod0
      (LubinTateAction hq hπmax f hf0 hf1 hf (↑u⁻¹ : 𝒪[K.carrier]))
      (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (↑u⁻¹ : 𝒪[K.carrier]))] at hmul2
    have hcalc : (↑u⁻¹ : 𝒪[K.carrier]) * (π ^ (n - 1) * (u : 𝒪[K.carrier])) = π ^ (n - 1) := by
      have heq2 : (↑u⁻¹ : 𝒪[K.carrier]) * (π ^ (n - 1) * (u : 𝒪[K.carrier])) =
          π ^ (n - 1) * ((↑u⁻¹ : 𝒪[K.carrier]) * ↑u) := by ring
      rw [heq2, Units.inv_mul, mul_one]
    rw [hcalc] at hmul2
    exact hmul2.symm
  exact lubinTateActionAtTorsionPoint_pi_pow_pred_ne_zero_of_mem_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hmem hxn hkey

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**主単数`1+cπ^n`は`x`を固定する**:
`(1+c*π^n)·x = x`(`x∈Λ_n`ならば任意の`n`・任意の`c∈𝒪_K`について)。
`lubinTateAction_add`(`(1+cπ^n)·x=F_f(1·x,(cπ^n)·x)`)・
`lubinTateActionAtTorsionPoint_one`(`1·x=x`)・
`pi_pow_action_action_eq_zero`(`(cπ^n)·x=0`)を代入し、
`aeval_formalGroupLaw_eq_of_snd_eq_zero`(`Found/PGC/LubinTateIdentityLaw.lean`、
`F_f(y₀,0)=y₀`)で結論する。`F_f`加法の逆元・結合律は一切使わず、
単位元則の評価版だけで閉じる——単射性一般には引き算が要る
(前々回の記録)が、**この特別な場合**(主単数が`x`を固定すること)
は避けられた。`Gal(K(Λ_n)/K)`の作用が`(𝒪_K/π^n)^×`を経由して
well-defined であることの一部——`(𝒪_K)^×`から`Λ_n`への作用が
`1+π^n𝒪_K`を法として不変であることの直接の確認。 -/
theorem one_add_mul_pi_pow_action_eq_self
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (c : 𝒪[K.carrier]) :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (1 + c * π ^ n) =
      ⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem⟩ := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hadd := lubinTateAction_add K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1 (c * π ^ n)
  have hcpi0 : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
      (c * π ^ n) = 0 := pi_pow_action_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem c
  simp only [lubinTateActionAtTorsionPoint_one, hcpi0] at hadd
  rw [hadd]
  show MvPowerSeries.aeval _ (formalGroupLaw hq hπmax f hf0 hf1 hf) = _
  rw [aeval_formalGroupLaw_eq_of_snd_eq_zero hq hπmax f hf0 hf1 hf hπne0 _ (by simp)]
  rfl

/-- `lubinTateEvalAtPoint`は評価点が(等式として)一致すれば同じ値を
返す——`lubinTateEvalAtPoint_eq_zero_of_eq_zero`の`w=0`限定を一般化
したもの。依存型(`HasEval`が評価点自身に依存する)への`rw`が
「motive is not type correct」で失敗する定型の罠を、`subst`で
回避する。 -/
theorem lubinTateEvalAtPoint_congr {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {z w : adjoinIntegers K x} (hz : PowerSeries.HasEval z) (hzw : z = w)
    (g : PowerSeries (𝒪[K.carrier])) :
    lubinTateEvalAtPoint K x z hz g = lubinTateEvalAtPoint K x w (hzw ▸ hz) g := by
  subst hzw
  rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`a↦a·x`は`(1+π^n𝒪_K)`を法として
乗法的に不変**: `(a*(1+cπ^n))·x = a·x`(任意の`a,c∈𝒪_K`)。
`lubinTateAction_mul`(`a·((1+cπ^n)·x)=(a*(1+cπ^n))·x`)と
`one_add_mul_pi_pow_action_eq_self`(`(1+cπ^n)·x=x`)を組み合わせる
だけ——`lubinTateEvalAtPoint_congr`で依存型の書き換えを処理する。
`(𝒪_K)^×`から`ψ_nの根`への作用(`unit_action_mem_iteratedLubinTate
PsiTorsionPoints`)が`(𝒪_K/π^n)^×`を経由してwell-definedであること
の直接の確認——`u`と`u*(1+cπ^n)`(同じ剰余類の代表元)が`x`に
同じ値を与える。 -/
theorem mul_one_add_pi_pow_action_eq
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a c : 𝒪[K.carrier]) :
    lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (a * (1 + c * π ^ n)) =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hmul := lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a (1 + c * π ^ n)
  have hone : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
      (1 + c * π ^ n) = ⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem⟩ :=
    one_add_mul_pi_pow_action_eq_self K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem c
  rw [← hmul, lubinTateEvalAtPoint_congr K x
    (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (1 + c * π ^ n))
    hone]
  rfl

end ABC3.Found.PGC
