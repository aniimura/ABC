/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.DiscreteValuationRing.Basic
import ABC3.Found.GaloisRep.RamifiedValuationBridge
import ABC3.Meta.Claim

/-!
# 第 1373 ブロック —— **付値の延長公式 `v_C = (v_A)^e`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1372 が要る `hpe` の中身

第 1372（`exists_h2_h1_of_bad_prime_ram`）は

    hpe : ∀ x : L, v_{𝔪_C}(algebraMap L L_v′ x) = (v_p(x))^e

を受ける。☆`L → L_v`（完備化、第 964 で `hp` は等号）と
`L_v → L_v′`（有限拡大、**本ブロック**）を合成すればこれが出る。

★本ブロックの中身は
**「`𝔪_A C = 𝔪_C^e` なら `v_C(φ y) = (v_A y)^e`」**——
すなわち分岐指数の定義そのものを付値の言葉に直すことである。

☆道は 2 段:

1. `A` の素元 `π` について `v_C(φ π) = exp(−e)`
   ——`span{φ π} = 𝔪_A C = 𝔪_C^e = span{ϖ^e}` なので `φ π` と `ϖ^e` は同伴
2. `A` の任意の元は `単元 × π^n` なので乗法性で伸ばし、
   商体へは `y = r/s` で伸ばす
-/

namespace ABC3.Found.GenEll

open ABC3.Meta IsDedekindDomain

section IntValuation

variable {A : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]
variable {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C]
variable [Algebra A C]

/-- ★★★★**同伴な元の `intValuation` は等しい**（第 1373）。 -/
theorem intValuation_eq_of_associated {x y : C} (h : Associated x y) :
    (IsDiscreteValuationRing.maximalIdeal C).intValuation x
      = (IsDiscreteValuationRing.maximalIdeal C).intValuation y := by
  obtain ⟨u, hu⟩ := h
  have hunit : (IsDiscreteValuationRing.maximalIdeal C).intValuation (u : C) = 1 := by
    rw [HeightOneSpectrum.intValuation_eq_one_iff]
    intro hmem
    exact (IsDiscreteValuationRing.maximalIdeal C).isPrime.ne_top
      (Ideal.eq_top_of_isUnit_mem _ hmem u.isUnit)
  rw [← hu, map_mul, hunit, mul_one]

/-- ★★★★★★★★**素元の像の `intValuation` は `exp(−e)`**（第 1373）。 -/
theorem intValuation_algebraMap_uniformizer {e : ℕ} {π : A} (hπ : Irreducible π)
    (hIe : (IsLocalRing.maximalIdeal A).map (algebraMap A C)
      = (IsLocalRing.maximalIdeal C) ^ e) :
    (IsDiscreteValuationRing.maximalIdeal C).intValuation (algebraMap A C π)
      = WithZero.exp (-(e : ℤ)) := by
  obtain ⟨ϖ, hirr⟩ := IsDiscreteValuationRing.exists_irreducible C
  have hspan : Ideal.span {algebraMap A C π} = Ideal.span {ϖ ^ e} := by
    rw [← Ideal.span_singleton_pow, ← hirr.maximalIdeal_eq, ← hIe,
      hπ.maximalIdeal_eq, Ideal.map_span]
    simp
  have hassoc : Associated (algebraMap A C π) (ϖ ^ e) :=
    Ideal.span_singleton_eq_span_singleton.1 hspan
  rw [intValuation_eq_of_associated hassoc, map_pow,
    HeightOneSpectrum.intValuation_singleton _ hirr.ne_zero hirr.maximalIdeal_eq,
    ABC3.Found.GaloisRep.withZero_exp_pow]
  congr 1
  ring

/-- ★★★★★★★★★★★★
**`intValuation` の延長公式**——★**無条件**（第 1373）。 -/
theorem intValuation_algebraMap_pow {e : ℕ} (he : 1 ≤ e)
    (hIe : (IsLocalRing.maximalIdeal A).map (algebraMap A C)
      = (IsLocalRing.maximalIdeal C) ^ e) (a : A) :
    (IsDiscreteValuationRing.maximalIdeal C).intValuation (algebraMap A C a)
      = ((IsDiscreteValuationRing.maximalIdeal A).intValuation a) ^ e := by
  rcases eq_or_ne a 0 with rfl | ha
  · rw [map_zero, map_zero, map_zero, zero_pow (by omega)]
  obtain ⟨ϖA, hirrA⟩ := IsDiscreteValuationRing.exists_irreducible A
  obtain ⟨n, u, hnu⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible ha hirrA
  have hunitA : (IsDiscreteValuationRing.maximalIdeal A).intValuation (u : A) = 1 := by
    rw [HeightOneSpectrum.intValuation_eq_one_iff]
    intro hmem
    exact (IsDiscreteValuationRing.maximalIdeal A).isPrime.ne_top
      (Ideal.eq_top_of_isUnit_mem _ hmem u.isUnit)
  have hunitC : (IsDiscreteValuationRing.maximalIdeal C).intValuation
      (algebraMap A C (u : A)) = 1 := by
    rw [HeightOneSpectrum.intValuation_eq_one_iff]
    intro hmem
    exact (IsDiscreteValuationRing.maximalIdeal C).isPrime.ne_top
      (Ideal.eq_top_of_isUnit_mem _ hmem (u.isUnit.map (algebraMap A C)))
  rw [hnu, map_mul, map_mul, map_mul, hunitA, hunitC, one_mul, one_mul,
    map_pow, map_pow, map_pow,
    HeightOneSpectrum.intValuation_singleton _ hirrA.ne_zero hirrA.maximalIdeal_eq,
    intValuation_algebraMap_uniformizer hirrA hIe,
    ABC3.Found.GaloisRep.withZero_exp_pow, ABC3.Found.GaloisRep.withZero_exp_pow,
    ABC3.Found.GaloisRep.withZero_exp_pow]
  congr 1
  ring

end IntValuation

section Valuation

variable {A : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]
variable {K : Type*} [Field K] [Algebra A K] [IsFractionRing A K]
variable {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C]
variable {L : Type*} [Field L] [Algebra C L] [IsFractionRing C L]
variable [Algebra A C] [Algebra A L] [Algebra K L]
variable [IsScalarTower A C L] [IsScalarTower A K L]

/-- ★★★★★★★★★★★★★★★★
**付値の延長公式 `v_C(φ y) = (v_A y)^e`**——★**無条件**（第 1373）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが第 1372（`exists_h2_h1_of_bad_prime_ram`）の `hpe` の中身である。 -/
theorem valuation_algebraMap_pow {e : ℕ} (he : 1 ≤ e)
    (hIe : (IsLocalRing.maximalIdeal A).map (algebraMap A C)
      = (IsLocalRing.maximalIdeal C) ^ e) (y : K) :
    (HeightOneSpectrum.valuation L (IsDiscreteValuationRing.maximalIdeal C))
        (algebraMap K L y)
      = ((HeightOneSpectrum.valuation K (IsDiscreteValuationRing.maximalIdeal A)) y) ^ e := by
  obtain ⟨r, s, hs, rfl⟩ := IsFractionRing.div_surjective (A := A) y
  have htower : ∀ a : A, algebraMap K L (algebraMap A K a)
      = algebraMap C L (algebraMap A C a) := by
    intro a
    rw [← IsScalarTower.algebraMap_apply A K L, ← IsScalarTower.algebraMap_apply A C L]
  have hL : algebraMap K L (algebraMap A K r / algebraMap A K s)
      = algebraMap C L (algebraMap A C r) / algebraMap C L (algebraMap A C s) := by
    rw [map_div₀, htower, htower]
  rw [hL, map_div₀, map_div₀,
    HeightOneSpectrum.valuation_of_algebraMap, HeightOneSpectrum.valuation_of_algebraMap,
    HeightOneSpectrum.valuation_of_algebraMap, HeightOneSpectrum.valuation_of_algebraMap,
    div_pow, intValuation_algebraMap_pow he hIe, intValuation_algebraMap_pow he hIe]

end Valuation

/-! ## ★出典の紐付け(`.src`) -/

def intValuation_algebraMap_pow.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(intValuation の延長公式 v_C = (v_A)^e。★無条件)",
    sectionId := "genell-thm-3-8" }

def valuation_algebraMap_pow.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(付値の延長公式 v_C(φ y) = (v_A y)^e。★無条件)",
    sectionId := "genell-thm-3-8" }

def valuation_algebraMap_pow.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_pow_eq_map_maximalIdeal(第 1369、証明済み。𝔪_A C = 𝔪_C^e)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_eq_map_maximalIdeal") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1373）**——第 1372 が要る `hpe` の中身である。" ++
       "☆`L → L_v`（完備化、第 964 で等号）と `L_v → L_v′`（本ブロック）を合成する。") 19 ]

end ABC3.Found.GenEll
