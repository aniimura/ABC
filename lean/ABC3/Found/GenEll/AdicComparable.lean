/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.AdicCompletion.Basic
import Mathlib.RingTheory.DiscreteValuationRing.Basic
import Mathlib.LinearAlgebra.SModEq.Basic
import ABC3.Meta.Claim

/-!
# 第 1367 ブロック —— **`I`-進完備性は「同じ位相」を与えるイデアルに移る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1366 の出口を「使える形」にする

第 1366 は `C = 整閉包(A, L)` について **`m_A C`-進完備**を与えた。
☆ところが下流（`exists_h2_h1_of_bad_prime`、第 1320）が要求するのは
**`m_C`-進完備**（`IsAdicComplete (maximalIdeal C) C`）である。

★`C` は DVR なので `m_A C = m_C^e`（分岐指数 `e ≥ 1`）であり、
2 つのイデアルは**同じ進位相**を定める。☆本ブロックはその移送を作る:

> `I ≤ J` かつ `J^e ≤ I`（`e ≥ 1`）なら `IsAdicComplete I R → IsAdicComplete J R`

★Hausdorff 側は `J^(e·n) = (J^e)^n ≤ I^n`、
precomplete 側は `n ↦ f (e·n)` を `I`-Cauchy 列として拾い、
`n ≤ e·n` で戻す。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {R : Type*} [CommRing R]

/-- ★★★★**環自身では `I^n • ⊤` は `I^n`**（第 1367）。 -/
theorem mem_pow_smul_top_self (I : Ideal R) (n : ℕ) (x : R) :
    x ∈ (I ^ n • (⊤ : Submodule R R)) ↔ x ∈ I ^ n := by
  rw [← Ideal.one_eq_top, Ideal.smul_eq_mul, mul_one]

/-- ★★★★**イデアルの包含は冪に移る**（第 1367）。 -/
theorem pow_le_pow_of_le {I J : Ideal R} (h : I ≤ J) (n : ℕ) : I ^ n ≤ J ^ n := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [pow_succ, pow_succ]
    exact Ideal.mul_mono ih h

/-- ★★★★★★★★**`IsHausdorff` は「粗い方」から「細かい方」へ移る**（第 1367）。 -/
theorem isHausdorff_of_pow_le {I J : Ideal R} {e : ℕ} (hJI : J ^ e ≤ I)
    [IsHausdorff I R] : IsHausdorff J R where
  haus' := by
    intro x hx
    refine IsHausdorff.haus' (I := I) (M := R) x ?_
    intro n
    have h := hx (e * n)
    rw [SModEq.zero, mem_pow_smul_top_self] at h ⊢
    have hle : J ^ (e * n) ≤ I ^ n := by
      rw [pow_mul]
      exact pow_le_pow_of_le hJI n
    exact hle h

/-- ★★★★★★★★**`IsPrecomplete` も移る**（第 1367）。 -/
theorem isPrecomplete_of_pow_le {I J : Ideal R} (hIJ : I ≤ J) {e : ℕ} (he : 1 ≤ e)
    (hJI : J ^ e ≤ I) [IsPrecomplete I R] : IsPrecomplete J R where
  prec' := by
    intro f hf
    have hcauchy : ∀ {m n : ℕ}, m ≤ n →
        f (e * m) ≡ f (e * n) [SMOD (I ^ m • (⊤ : Submodule R R))] := by
      intro m n hmn
      have h := hf (Nat.mul_le_mul_left e hmn)
      rw [SModEq.sub_mem, mem_pow_smul_top_self] at h ⊢
      have hle : J ^ (e * m) ≤ I ^ m := by
        rw [pow_mul]
        exact pow_le_pow_of_le hJI m
      exact hle h
    obtain ⟨y, hy⟩ := IsPrecomplete.prec' (I := I) (M := R) (fun n => f (e * n)) hcauchy
    refine ⟨y, fun n => ?_⟩
    have h1 : f n ≡ f (e * n) [SMOD (J ^ n • (⊤ : Submodule R R))] :=
      hf (Nat.le_mul_of_pos_left n he)
    have h2 : f (e * n) ≡ y [SMOD (J ^ n • (⊤ : Submodule R R))] := by
      refine SModEq.mono ?_ (hy n)
      intro z hz
      rw [mem_pow_smul_top_self] at hz ⊢
      exact pow_le_pow_of_le hIJ n hz
    exact h1.trans h2

/-- ★★★★★★★★★★★★
**`I`-進完備性は「同じ位相」を与えるイデアルに移る**——★**無条件**（第 1367）。 -/
theorem isAdicComplete_of_pow_le {I J : Ideal R} (hIJ : I ≤ J) {e : ℕ} (he : 1 ≤ e)
    (hJI : J ^ e ≤ I) [IsAdicComplete I R] : IsAdicComplete J R where
  toIsHausdorff := isHausdorff_of_pow_le hJI
  toIsPrecomplete := isPrecomplete_of_pow_le hIJ he hJI

/-! ## ★★★★★★★★★★★★DVR での使い方 -/

section Dvr

variable {R : Type*} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]

/-- ★★★★★★★★★★★★★★★★
**DVR では `⊥` でも `⊤` でもないイデアルに沿う完備性は極大イデアルに沿う完備性**
——★**無条件**（第 1367）。

★★★これが第 1366 の出口（`m_A C`-進完備）を
下流が要求する形（`m_C`-進完備）に直す段である。 -/
theorem isAdicComplete_maximalIdeal_of_isAdicComplete (I : Ideal R) (hI0 : I ≠ ⊥)
    (hItop : I ≠ ⊤) [IsAdicComplete I R] :
    IsAdicComplete (IsLocalRing.maximalIdeal R) R := by
  obtain ⟨ϖ, hirr⟩ := IsDiscreteValuationRing.exists_irreducible R
  obtain ⟨e, he⟩ := IsDiscreteValuationRing.ideal_eq_span_pow_irreducible hI0 hirr
  have hIe : I = (IsLocalRing.maximalIdeal R) ^ e := by
    rw [he, hirr.maximalIdeal_eq, Ideal.span_singleton_pow]
  have he1 : 1 ≤ e := by
    rcases Nat.eq_zero_or_pos e with h | h
    · exact absurd (by rw [hIe, h, pow_zero, Ideal.one_eq_top]) hItop
    · exact h
  refine isAdicComplete_of_pow_le (I := I) ?_ he1 ?_
  · rw [hIe]
    exact Ideal.pow_le_self (by omega)
  · rw [hIe]

end Dvr

/-! ## ★出典の紐付け(`.src`) -/

def isAdicComplete_of_pow_le.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(I-進完備性は同じ位相を与えるイデアルに移る。★無条件)",
    sectionId := "genell-thm-3-8" }

def isAdicComplete_maximalIdeal_of_isAdicComplete.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(DVR では ⊥ でも ⊤ でもないイデアルの完備性は極大イデアルの完備性。★無条件)",
    sectionId := "genell-thm-3-8" }

def isAdicComplete_maximalIdeal_of_isAdicComplete.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1367）**——第 1366 の出口（`m_A C`-進完備）を" ++
       "下流（第 1320 `exists_h2_h1_of_bad_prime`）が要求する `m_C`-進完備に直す段である。" ++
       "☆`m_A C = m_C^e`（分岐指数）なので 2 つは同じ進位相を定める。") 19 ]

end ABC3.Found.GenEll
