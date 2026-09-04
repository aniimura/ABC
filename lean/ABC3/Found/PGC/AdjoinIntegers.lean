import ABC3.Found.PGC.LubinTateDistinguishedSeparable
import ABC3.Found.PGC.LubinTateActionEndomorphism

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

end ABC3.Found.PGC
