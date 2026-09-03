/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import Mathlib.RingTheory.AdicCompletion.Topology

/-!
# 第 897 ブロック —— **★★★★★★★★★★★★★★★★★★★★完備化の整数環は
`I` 進完備である**（`Found`）

## ★★★★★★★★★★これは何か

`Check/GenEll/AdicCompleteMissing.lean`（第 896）で測ったとおり、
`Lemma 3.5` の完備化の橋は**インスタンス 1 つ**に落ちていた:

    `IsAdicComplete (maximalIdeal (v.adicCompletionIntegers K)) (v.adicCompletionIntegers K)`

★本ブロックはそれを**作った**。mathlib には無いものである。

## ★道は 3 段

| 段 | 定理 | 内容 |
|---|---|---|
| 1 | `isClosed_adicCompletionIntegers` | 整数環は閉集合（開部分群は閉） |
| 2 | `completeSpace_adicCompletionIntegers` | 完備空間の閉集合は完備 |
| 3 | `isAdic_maximalIdeal_adicCompletionIntegers` | 位相は `𝔪` 進位相である |

★段 3 が鍵である。`IsDiscreteValuationRing` の一次元子 `ϖ` を取ると

    `𝔪^n = {y | v(y) ≤ v(ϖ)^n}`   （mathlib の `maximalIdeal_pow_eq_setOf_le_v_algebraMap_pow`）

であり、右辺は**開集合**（閉球は非アルキメデス位相で開）である。
逆に `v(ϖ) < 1` なので `exists_pow_lt₀` が任意の球を `𝔪^n` で押さえる。

☆あとは mathlib の `IsAdic.isAdicComplete_iff`（進位相なら
`IsAdicComplete ↔ CompleteSpace ∧ T2Space`）を当てるだけである。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField IsLocalRing Valued Valuation

/-! ## ★非アルキメデス位相の閉球は開集合 -/

/-- ★**付値の閉球は開集合**（非アルキメデスなので）。 -/
theorem isOpen_valued_le_embedding {R Γ₀ : Type} [Ring R]
    [LinearOrderedCommGroupWithZero Γ₀] [Valued R Γ₀]
    (γ : (MonoidWithZeroHom.ValueGroup₀ (.ofClass (Valued.v : Valuation R Γ₀)))ˣ) :
    IsOpen {x : R | Valued.v x ≤ MonoidWithZeroHom.ValueGroup₀.embedding γ.1} := by
  rw [isOpen_iff_mem_nhds]
  intro x hx
  rw [Valued.mem_nhds]
  refine ⟨γ, ?_⟩
  intro y hy
  simp only [Set.mem_setOf_eq] at *
  rw [Valuation.restrict_lt_iff_lt_embedding] at hy
  rw [← sub_add_cancel y x]
  exact le_trans (Valuation.map_add _ _ _) (max_le (le_of_lt hy) hx)

/-! ## ★★段 1-2——整数環は閉であり、したがって完備である -/

/-- ★★**完備化の整数環は閉集合**——開部分群は閉だから。 -/
theorem isClosed_adicCompletionIntegers (K : Type) [Field K] [NumberField K]
    (v : HeightOneSpectrum (𝓞 K)) :
    IsClosed ((v.adicCompletionIntegers K : Set (v.adicCompletion K))) := by
  have hset : ((v.adicCompletionIntegers K : Set (v.adicCompletion K)))
      = ((v.adicCompletionIntegers K).toSubring.toAddSubgroup
          : Set (v.adicCompletion K)) := rfl
  have hopen : IsOpen ((v.adicCompletionIntegers K).toSubring.toAddSubgroup
      : Set (v.adicCompletion K)) := by
    rw [isOpen_iff_mem_nhds]
    intro x hx
    rw [Valued.mem_nhds]
    refine ⟨1, ?_⟩
    intro y hy
    have hx' : Valued.v x ≤ 1 := (HeightOneSpectrum.mem_adicCompletionIntegers _ _ _).1 hx
    refine (HeightOneSpectrum.mem_adicCompletionIntegers _ _ _).2 ?_
    rw [← sub_add_cancel y x]
    exact le_trans (Valuation.map_add _ _ _)
      (max_le (le_of_lt (by simpa using hy)) hx')
  rw [hset]
  exact AddSubgroup.isClosed_of_isOpen _ hopen

/-- ★★**完備化の整数環は完備空間**——完備空間の閉集合だから。 -/
instance completeSpace_adicCompletionIntegers (K : Type) [Field K] [NumberField K]
    (v : HeightOneSpectrum (𝓞 K)) :
    CompleteSpace (v.adicCompletionIntegers K) :=
  (isClosed_adicCompletionIntegers K v).completeSpace_coe

/-! ## ★★★★★★★★段 3——位相は `𝔪` 進位相である -/

/-- ★★★★★★★★**完備化の整数環の位相は極大イデアル進位相である**。

☆一次元子 `ϖ` を取ると `𝔪^n = {y | v(y) ≤ v(ϖ)^n}` であり、
右辺は開集合。逆に `v(ϖ) < 1` なので `exists_pow_lt₀` が任意の球を押さえる。 -/
theorem isAdic_maximalIdeal_adicCompletionIntegers (K : Type) [Field K] [NumberField K]
    (v : HeightOneSpectrum (𝓞 K)) :
    IsAdic (IsLocalRing.maximalIdeal (v.adicCompletionIntegers K)) := by
  obtain ⟨ϖ, hϖ⟩ : ∃ ϖ : (v.adicCompletionIntegers K), Irreducible ϖ :=
    ⟨_, (IsDiscreteValuationRing.exists_irreducible _).choose_spec⟩
  have hints := Valuation.valuationSubring.integers
    (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
  have hlt1 : (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
      (algebraMap _ _ ϖ) < 1 := hints.valuation_irreducible_lt_one hϖ
  have hne0 : (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
      (algebraMap _ _ ϖ) ≠ 0 := by
    simp only [ne_eq, map_eq_zero]
    intro h
    exact hϖ.ne_zero (Subtype.ext h)
  have hrne : (Valued.v : Valuation (v.adicCompletion K)
      (WithZero (Multiplicative ℤ))).restrict (algebraMap _ _ ϖ) ≠ 0 := by
    intro hz
    exact hne0 (by rw [← Valuation.embedding_restrict, hz, map_zero])
  set g : (MonoidWithZeroHom.ValueGroup₀ (.ofClass
      (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))))ˣ :=
    Units.mk0 _ hrne with hgdef
  have hemb : ∀ n : ℕ, MonoidWithZeroHom.ValueGroup₀.embedding ((g ^ n).1)
      = (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
        (algebraMap _ _ ϖ) ^ n := by
    intro n
    rw [Units.val_pow_eq_pow_val, hgdef, Units.val_mk0, map_pow, Valuation.embedding_restrict]
  have hset : ∀ n : ℕ, ((IsLocalRing.maximalIdeal (v.adicCompletionIntegers K) ^ n :
      Ideal (v.adicCompletionIntegers K)) : Set (v.adicCompletionIntegers K))
      = {y : (v.adicCompletionIntegers K) |
          (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
              (algebraMap _ _ y)
            ≤ (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
              (algebraMap _ _ ϖ) ^ n} :=
    fun n => hints.maximalIdeal_pow_eq_setOf_le_v_algebraMap_pow hϖ n
  rw [isAdic_iff]
  refine ⟨fun n => ?_, fun s hs => ?_⟩
  · rw [hset n]
    have heq : {y : (v.adicCompletionIntegers K) |
          (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
              (algebraMap _ _ y)
            ≤ (Valued.v : Valuation (v.adicCompletion K) (WithZero (Multiplicative ℤ)))
              (algebraMap _ _ ϖ) ^ n}
        = (Subtype.val : (v.adicCompletionIntegers K) → v.adicCompletion K) ⁻¹'
          {x : v.adicCompletion K | Valued.v x
            ≤ MonoidWithZeroHom.ValueGroup₀.embedding ((g ^ n).1)} := by
      ext y
      simp only [Set.mem_setOf_eq, Set.mem_preimage, hemb n]
      rfl
    rw [heq]
    exact (isOpen_valued_le_embedding _).preimage continuous_subtype_val
  · rw [nhds_subtype_eq_comap, Filter.mem_comap] at hs
    obtain ⟨t, ht, hts⟩ := hs
    have ht0 : t ∈ nhds (0 : v.adicCompletion K) := ht
    rw [Valued.mem_nhds_zero] at ht0
    obtain ⟨γ, hγ⟩ := ht0
    have hγ0 : MonoidWithZeroHom.ValueGroup₀.embedding γ.1 ≠ 0 := by
      intro hz
      exact γ.ne_zero (MonoidWithZeroHom.ValueGroup₀.embedding_strictMono.injective
        (by rw [hz, map_zero]))
    obtain ⟨n, hn⟩ := exists_pow_lt₀ hlt1 (Units.mk0 _ hγ0)
    refine ⟨n, ?_⟩
    rw [hset n]
    intro y hy
    refine hts ?_
    refine hγ ?_
    simp only [Set.mem_setOf_eq] at hy ⊢
    rw [Valuation.restrict_lt_iff_lt_embedding]
    exact lt_of_le_of_lt hy hn

/-! ## ★★★★★★★★★★★★★★★★★★★★結論 -/

/-- ★★★★★★★★★★★★★★★★★★★★**完備化の整数環は `𝔪` 進完備である**。

★これが `Check/GenEll/AdicCompleteMissing.lean`（第 896）が名指しした
**欠けていたインスタンス**である。これで `Lemma 3.5` の局所の道具
（`tateCurveAt`・`tatePhi`・`minDeltaExp_eq_mul_of_tateParamR` 等）が
**数体の素点の完備化に実際に当てられる**ようになった。 -/
instance isAdicComplete_adicCompletionIntegers (K : Type) [Field K] [NumberField K]
    (v : HeightOneSpectrum (𝓞 K)) :
    IsAdicComplete (IsLocalRing.maximalIdeal (v.adicCompletionIntegers K))
      (v.adicCompletionIntegers K) :=
  (isAdic_maximalIdeal_adicCompletionIntegers K v).isAdicComplete_iff.2
    ⟨completeSpace_adicCompletionIntegers K v, inferInstance⟩

/-! ## ★標数は 0 である（インスタンスが mathlib に無い） -/

/-- ★数体の完備化は標数 `0`。 -/
instance charZero_adicCompletion (K : Type) [Field K] [NumberField K]
    (v : HeightOneSpectrum (𝓞 K)) : CharZero (v.adicCompletion K) :=
  charZero_of_injective_algebraMap (algebraMap K (v.adicCompletion K)).injective

/-- ★その整数環も標数 `0`。 -/
instance charZero_adicCompletionIntegers (K : Type) [Field K] [NumberField K]
    (v : HeightOneSpectrum (𝓞 K)) : CharZero (v.adicCompletionIntegers K) :=
  (algebraMap (v.adicCompletionIntegers K) (v.adicCompletion K)).charZero

end ABC3.Found.GaloisRep
