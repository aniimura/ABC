/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicCompleteIntegers

/-!
# 第 943 ブロック —— **★★★★★★★★★★★★完備な付値体の整数環は `𝔪` 進完備**（`Found`）

## ★★★★★★★★★★これは何か

第 897 は `v.adicCompletionIntegers K`（数体の素点の完備化）について
`IsAdicComplete` を作った。★本ブロックは**同じ証明を一般の完備な付値体へ**広げる。

☆動機: `Lemma 3.5` の組み立てには、`μ_l` や `l`-捉れ点を含む
**`Lv` の有限 Galois 拡大**が要る（第 942 の測定）。
その整数環でも同じ `IsAdicComplete` が要るので、一般の形で取っておく。

★道は第 897 と同じ 3 段である:

1. 付値環は開部分群なので閉集合
2. 完備空間の閉集合は完備
3. 一次元子 `ϖ` を取ると `𝔪^n = {y | v(y) ≤ v(ϖ)^n}` で位相は `𝔪` 進位相
-/

namespace ABC3.Found.GaloisRep

open IsLocalRing Valued Valuation

variable {K Γ₀ : Type} [Field K] [LinearOrderedCommGroupWithZero Γ₀] [MulArchimedean Γ₀]
  [Valued K Γ₀]

/-- ★★**付値環は閉集合**——開部分群は閉だから。 -/
theorem isClosed_valuationSubring :
    IsClosed (((Valued.v : Valuation K Γ₀).valuationSubring : Set K)) := by
  have hset : (((Valued.v : Valuation K Γ₀).valuationSubring : Set K))
      = (((Valued.v : Valuation K Γ₀).valuationSubring.toSubring.toAddSubgroup : Set K)) := rfl
  have hopen : IsOpen
      (((Valued.v : Valuation K Γ₀).valuationSubring.toSubring.toAddSubgroup : Set K)) := by
    rw [isOpen_iff_mem_nhds]
    intro x hx
    rw [Valued.mem_nhds]
    refine ⟨1, ?_⟩
    intro y hy
    have hx' : Valued.v x ≤ 1 := hx
    show Valued.v y ≤ 1
    rw [← sub_add_cancel y x]
    exact le_trans (Valuation.map_add _ _ _)
      (max_le (le_of_lt (by simpa using hy)) hx')
  rw [hset]
  exact AddSubgroup.isClosed_of_isOpen _ hopen

/-- ★★**完備な付値体の付値環は完備空間**。 -/
theorem completeSpace_valuationSubring [CompleteSpace K] :
    CompleteSpace ((Valued.v : Valuation K Γ₀).valuationSubring) :=
  (isClosed_valuationSubring (K := K) (Γ₀ := Γ₀)).completeSpace_coe

/-- ★★★★★★★★**付値環の位相は `𝔪` 進位相である**（離散付値環なら）。 -/
theorem isAdic_maximalIdeal_valuationSubring
    [IsDiscreteValuationRing ((Valued.v : Valuation K Γ₀).valuationSubring)]
    [IsTopologicalRing ((Valued.v : Valuation K Γ₀).valuationSubring)] :
    IsAdic (IsLocalRing.maximalIdeal ((Valued.v : Valuation K Γ₀).valuationSubring)) := by
  obtain ⟨ϖ, hϖ⟩ : ∃ ϖ : ((Valued.v : Valuation K Γ₀).valuationSubring), Irreducible ϖ :=
    ⟨_, (IsDiscreteValuationRing.exists_irreducible _).choose_spec⟩
  have hints := Valuation.valuationSubring.integers (Valued.v : Valuation K Γ₀)
  have hlt1 : (Valued.v : Valuation K Γ₀) (algebraMap _ _ ϖ) < 1 :=
    hints.valuation_irreducible_lt_one hϖ
  have hne0 : (Valued.v : Valuation K Γ₀) (algebraMap _ _ ϖ) ≠ 0 := by
    simp only [ne_eq, map_eq_zero]
    intro h
    exact hϖ.ne_zero (Subtype.ext h)
  have hrne : (Valued.v : Valuation K Γ₀).restrict (algebraMap _ _ ϖ) ≠ 0 := by
    intro hz
    exact hne0 (by rw [← Valuation.embedding_restrict, hz, map_zero])
  set g : (MonoidWithZeroHom.ValueGroup₀ (.ofClass (Valued.v : Valuation K Γ₀)))ˣ :=
    Units.mk0 _ hrne with hgdef
  have hemb : ∀ n : ℕ, MonoidWithZeroHom.ValueGroup₀.embedding ((g ^ n).1)
      = (Valued.v : Valuation K Γ₀) (algebraMap _ _ ϖ) ^ n := by
    intro n
    rw [Units.val_pow_eq_pow_val, hgdef, Units.val_mk0, map_pow, Valuation.embedding_restrict]
  have hset : ∀ n : ℕ,
      ((IsLocalRing.maximalIdeal ((Valued.v : Valuation K Γ₀).valuationSubring) ^ n :
        Ideal ((Valued.v : Valuation K Γ₀).valuationSubring))
          : Set ((Valued.v : Valuation K Γ₀).valuationSubring))
      = {y : ((Valued.v : Valuation K Γ₀).valuationSubring) |
          (Valued.v : Valuation K Γ₀) (algebraMap _ _ y)
            ≤ (Valued.v : Valuation K Γ₀) (algebraMap _ _ ϖ) ^ n} :=
    fun n => hints.maximalIdeal_pow_eq_setOf_le_v_algebraMap_pow hϖ n
  rw [isAdic_iff]
  refine ⟨fun n => ?_, fun s hs => ?_⟩
  · rw [hset n]
    have heq : {y : ((Valued.v : Valuation K Γ₀).valuationSubring) |
          (Valued.v : Valuation K Γ₀) (algebraMap _ _ y)
            ≤ (Valued.v : Valuation K Γ₀) (algebraMap _ _ ϖ) ^ n}
        = (Subtype.val : ((Valued.v : Valuation K Γ₀).valuationSubring) → K) ⁻¹'
          {x : K | Valued.v x
            ≤ MonoidWithZeroHom.ValueGroup₀.embedding ((g ^ n).1)} := by
      ext y
      simp only [Set.mem_setOf_eq, Set.mem_preimage, hemb n]
      rfl
    rw [heq]
    exact (isOpen_valued_le_embedding _).preimage continuous_subtype_val
  · rw [nhds_subtype_eq_comap, Filter.mem_comap] at hs
    obtain ⟨t, ht, hts⟩ := hs
    have ht0 : t ∈ nhds (0 : K) := ht
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

/-- ★★★★★★★★★★★★**完備な離散付値体の整数環は `𝔪` 進完備**。

★第 897 の一般形である。`Lv` の有限 Galois 拡大の整数環にも当てられる。 -/
instance isAdicComplete_valuationSubring [CompleteSpace K]
    [IsDiscreteValuationRing ((Valued.v : Valuation K Γ₀).valuationSubring)]
    [IsTopologicalRing ((Valued.v : Valuation K Γ₀).valuationSubring)]
    [IsUniformAddGroup ((Valued.v : Valuation K Γ₀).valuationSubring)] :
    IsAdicComplete (IsLocalRing.maximalIdeal ((Valued.v : Valuation K Γ₀).valuationSubring))
      ((Valued.v : Valuation K Γ₀).valuationSubring) :=
  (isAdic_maximalIdeal_valuationSubring (K := K) (Γ₀ := Γ₀)).isAdicComplete_iff.2
    ⟨completeSpace_valuationSubring, inferInstance⟩

end ABC3.Found.GaloisRep
