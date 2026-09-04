import ABC3.Found.PGC.ValuationRingDVR
import Mathlib.RingTheory.AdicCompletion.Topology

/-!
# `𝒪[K.carrier]` は π-進完備(`sorry` 無し)

`Found/PGC/ValuationRingDVR.lean::valuationRing_isDVR`(任意の p進局所体
`K` の付値環は離散付値環)に続き、`Found/PGC/LubinTateAction*.lean`
(このセッションで確立した Lubin-Tate 形式群の理論一式)が要求する
残りの仮定——`IsAdicComplete (maximalIdeal) (𝒪[K.carrier])`・
`UniqueFactorizationMonoid`・`CharZero`——を、`K` が**実在の p進局所体**
であることから確立する。

## なぜこれが要るか

`Found/PGC/LubinTateAction{Add,PiPow,Weierstrass,Psi,PsiField}.lean` の
一連の定理は、抽象的な環 `A`(`[CommRing A][IsLocalRing A][IsDomain A]`、
さらに節目ごとに `[IsAdicComplete (maximalIdeal A) A]`・
`[UniqueFactorizationMonoid A]`・`[CharZero A]` を追加)について証明した。
本ファイルは、この `A` を実際の `𝒪[K.carrier]`(`K:PAdicLocalField p`)に
具体化するとき、追加した仮定がすべて**自動的に成り立つ**ことを示す
——抽象的な Lubin-Tate の理論一式が、原典([pGC] が前提する p進局所体)
に**そのまま適用できる**ことの検証。

## 証明の筋(`IsAdicComplete`)

`𝒪[K.carrier]` はコンパクト(既出、`Metric.closedBall 0 1`)。
`Valued.integer.compactSpace_iff_completeSpace_and_isDiscreteValuationRing_
and_finite_residueField` でコンパクト性から `CompleteSpace` を取り出す。
残るは「距離位相 = `maximalIdeal`-進位相」(`IsAdic`)であることの確認
——一意化元 `ϖ`(`IsDiscreteValuationRing.exists_irreducible`)を1つ取り、
`maximalIdeal^n = closedBall 0 ‖ϖ‖^n`(`Irreducible.maximalIdeal_pow_
eq_closedBall_pow`)と非アルキメデス距離での閉球の開性
(`IsUltrametricDist.isOpen_closedBall`)・`‖ϖ‖<1` から来る近傍基底性
(`exists_pow_lt_of_lt_one`)を組み合わせて `isAdic_iff` の両条件を確認する。
最後に `IsAdic.isAdicComplete_iff`(`CompleteSpace`・`T2Space` から
`IsAdicComplete` を得る)で結論する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

/-- `𝒪[K.carrier]` の `maximalIdeal`-進位相は、もとの距離位相と一致する
(`IsAdic`)——一意化元による閉球の記述と非アルキメデス距離の性質から。 -/
theorem isAdic_maximalIdeal_valuationRing {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IsAdic (IsLocalRing.maximalIdeal (𝒪[K.carrier])) := by
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hϖnorm0 : ‖ϖ‖ ≠ 0 := (Valued.integer.norm_irreducible_pos hϖ).ne'
  rw [isAdic_iff]
  constructor
  · intro n
    rw [show (IsLocalRing.maximalIdeal (𝒪[K.carrier])) = Valued.maximalIdeal K.carrier from rfl,
      Irreducible.maximalIdeal_pow_eq_closedBall_pow hϖ n]
    exact IsUltrametricDist.isOpen_closedBall 0 (pow_ne_zero n hϖnorm0)
  · intro s hs
    rw [Metric.mem_nhds_iff] at hs
    obtain ⟨ε, hεpos, hεsub⟩ := hs
    have hϖlt1 : ‖ϖ‖ < 1 := Valued.integer.norm_irreducible_lt_one hϖ
    obtain ⟨n, hn⟩ := exists_pow_lt_of_lt_one hεpos hϖlt1
    refine ⟨n, ?_⟩
    rw [show (IsLocalRing.maximalIdeal (𝒪[K.carrier])) = Valued.maximalIdeal K.carrier from rfl,
      Irreducible.maximalIdeal_pow_eq_closedBall_pow hϖ n]
    calc Metric.closedBall (0 : 𝒪[K.carrier]) (‖ϖ‖ ^ n)
        ⊆ Metric.ball (0 : 𝒪[K.carrier]) ε := Metric.closedBall_subset_ball hn
      _ ⊆ s := hεsub

/-- ★★★★★★★★★**`𝒪[K.carrier]` は `maximalIdeal`-進完備**——
`Found/PGC/LubinTateAction*.lean` が要求する `[IsAdicComplete
(maximalIdeal A) A]` が、実際の p進局所体 `K` の整数環 `A:=𝒪[K.carrier]`
に対して**自動的に**成り立つ。コンパクト性(既出)から `CompleteSpace`
を取り出し、`isAdic_maximalIdeal_valuationRing` と組み合わせるだけ。 -/
theorem isAdicComplete_valuationRing {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier]) := by
  have hcs : CompactSpace (𝒪[K.carrier]) := by
    have hset : (𝒪[K.carrier] : Set K.carrier) = Metric.closedBall 0 1 := by
      ext x
      simp only [Metric.mem_closedBall, dist_eq_norm, sub_zero]
      exact Valued.integer.mem_iff
    have hcompact : IsCompact (𝒪[K.carrier] : Set K.carrier) := hset ▸ isCompact_closedBall 0 1
    exact isCompact_iff_compactSpace.mp hcompact
  have hcomplete : CompleteSpace (𝒪[K.carrier]) :=
    (Valued.integer.compactSpace_iff_completeSpace_and_isDiscreteValuationRing_and_finite_residueField.mp
      hcs).1
  exact (isAdic_maximalIdeal_valuationRing K).isAdicComplete_iff.mpr ⟨hcomplete, inferInstance⟩

/-- `𝒪[K.carrier]` は一意分解整域——`IsDiscreteValuationRing` が
`IsPrincipalIdealRing` を含意し、PID は UFD であることから。
`Found/PGC/LubinTateActionPsiField.lean` の `[UniqueFactorizationMonoid A]`
(Gauss の補題)を実際の `K` に対して満たす。 -/
theorem uniqueFactorizationMonoid_valuationRing {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    UniqueFactorizationMonoid (𝒪[K.carrier]) := by
  haveI := valuationRing_isDVR K
  infer_instance

/-- `𝒪[K.carrier]` は標数0——`ℚ_[p]` からの代数射が単射であることから。
`Found/PGC/LubinTateActionPsiField.lean` の `[CharZero A]`(混標数の場合)
を実際の `K` に対して満たす。 -/
theorem charZero_valuationRing {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    CharZero (𝒪[K.carrier]) := by
  haveI : CharZero K.carrier := charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  infer_instance

end ABC3.Found.PGC
