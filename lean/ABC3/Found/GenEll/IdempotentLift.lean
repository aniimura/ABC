/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.AdicCompletion.Basic
import Mathlib.RingTheory.Idempotents
import ABC3.Meta.Claim

/-!
# 第 1357 ブロック —— **完備なイデアルに沿って冪等元が持ち上がる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——第 1356 で特定した primitive

`Skeleton/GenEll/AlphaLocalData.lean` の `exists_h2_h1_unipotent` を止めているのは
「完備離散付値環の有限拡大は再び完備離散付値環」であり、
★第 1356 でその根を **1 本の補題**に特定した:

> `IsAdicComplete I R` なら `R/I` の冪等元は `R` へ持ち上がる。

☆mathlib の `RingTheory/Idempotents.lean` は**冪零核に沿った持ち上げ**しか持たない
（`existsUnique_isIdempotentElem_eq_of_ker_isNilpotent` など）。

## ★★★★★★★★道——Newton の反復

`f(a) = 3a² − 2a³` と置くと**恒等式**

* `f(a) − a = −(a² − a)(2a − 1)`
* `f(a)² − f(a) = (a² − a)²(4a² − 4a − 3)`

が成り立つ（どちらも `ring` で閉じる）。★したがって `a² − a ∈ I^n` なら
`f(a)² − f(a) ∈ I^{2n}` かつ `f(a) ≡ a (mod I^n)` であり、
反復列は `I`-進 Cauchy 列になる。☆`IsPrecomplete` で極限を取り、
`IsHausdorff` で `e² = e` を出す。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {R : Type*} [CommRing R]

/-- ★★★★**Newton の 1 段** `f(a) = 3a² − 2a³`（第 1357）。 -/
def idemStep (a : R) : R := 3 * a ^ 2 - 2 * a ^ 3

/-- ★★★★★★`f(a) − a = −(a² − a)(2a − 1)`——★**恒等式**（第 1357）。 -/
theorem idemStep_sub_self (a : R) : idemStep a - a = -((a ^ 2 - a) * (2 * a - 1)) := by
  simp only [idemStep]
  ring

/-- ★★★★★★★★`f(a)² − f(a) = (a² − a)²(4a² − 4a − 3)`——★**恒等式**（第 1357）。 -/
theorem sq_idemStep_sub (a : R) :
    idemStep a ^ 2 - idemStep a = (a ^ 2 - a) ^ 2 * (4 * a ^ 2 - 4 * a - 3) := by
  simp only [idemStep]
  ring

/-- ★★★★**反復列**（第 1357）。 -/
def idemSeq (a : R) : ℕ → R
  | 0 => a
  | n + 1 => idemStep (idemSeq a n)

@[simp] theorem idemSeq_zero (a : R) : idemSeq a 0 = a := rfl

@[simp] theorem idemSeq_succ (a : R) (n : ℕ) :
    idemSeq a (n + 1) = idemStep (idemSeq a n) := rfl

/-- ★★★★★★★★**`a_n² − a_n ∈ I^{2^n}`**——★**無条件**（第 1357）。 -/
theorem idemSeq_sq_sub_mem (I : Ideal R) {a : R} (ha : a ^ 2 - a ∈ I) (n : ℕ) :
    idemSeq a n ^ 2 - idemSeq a n ∈ I ^ (2 ^ n) := by
  induction n with
  | zero => simpa using ha
  | succ n ih =>
      rw [idemSeq_succ, sq_idemStep_sub]
      have h2 : (idemSeq a n ^ 2 - idemSeq a n) ^ 2 ∈ I ^ (2 ^ (n + 1)) := by
        have h := Ideal.pow_mem_pow ih 2
        rwa [← pow_mul, ← pow_succ] at h
      exact Ideal.mul_mem_right _ _ h2

/-- ★★★★★★**`a_{n+1} − a_n ∈ I^{2^n}`**——★**無条件**（第 1357）。 -/
theorem idemSeq_succ_sub_mem (I : Ideal R) {a : R} (ha : a ^ 2 - a ∈ I) (n : ℕ) :
    idemSeq a (n + 1) - idemSeq a n ∈ I ^ (2 ^ n) := by
  rw [idemSeq_succ, idemStep_sub_self]
  exact (Ideal.neg_mem_iff _).2 (Ideal.mul_mem_right _ _ (idemSeq_sq_sub_mem I ha n))

/-- ★★★★★★★★**反復列は Cauchy**——★**無条件**（第 1357）。 -/
theorem idemSeq_sub_mem_of_le (I : Ideal R) {a : R} (ha : a ^ 2 - a ∈ I) {m n : ℕ}
    (hmn : m ≤ n) : idemSeq a n - idemSeq a m ∈ I ^ m := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hmn
  clear hmn
  induction k with
  | zero => simpa using Ideal.zero_mem _
  | succ k ih =>
      have hstep : idemSeq a (m + k + 1) - idemSeq a (m + k) ∈ I ^ m := by
        refine Ideal.pow_le_pow_right ?_ (idemSeq_succ_sub_mem I ha (m + k))
        exact le_trans (Nat.le_add_right m k) (Nat.le_of_lt (Nat.lt_two_pow_self))
      have := Ideal.add_mem _ hstep ih
      have heq : idemSeq a (m + k + 1) - idemSeq a m
          = (idemSeq a (m + k + 1) - idemSeq a (m + k)) + (idemSeq a (m + k) - idemSeq a m) := by
        ring
      rw [show m + (k + 1) = m + k + 1 from rfl, heq]
      exact this

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**完備なイデアルに沿って冪等元が持ち上がる**——★**無条件**（第 1357）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`a² − a ∈ I` なら、`a` と `mod I` で合同な冪等元 `e ∈ R` が存在する。

★★★これが第 1356 で特定した **mathlib に無い primitive** である
——`RingTheory/Idempotents.lean` は冪零核の場合しか持たない。 -/
theorem exists_isIdempotentElem_of_isAdicComplete (I : Ideal R) [IsAdicComplete I R]
    {a : R} (ha : a ^ 2 - a ∈ I) :
    ∃ e : R, e ^ 2 = e ∧ e - a ∈ I := by
  -- ★ Cauchy 列であることを `SModEq` の言葉に直す
  have aux : ∀ {m n : ℕ}, m ≤ n →
      idemSeq a m ≡ idemSeq a n [SMOD (I ^ m • ⊤ : Ideal R)] := by
    intro m n hmn
    rw [← Ideal.one_eq_top, Ideal.smul_eq_mul, mul_one, SModEq.sub_mem]
    have h := idemSeq_sub_mem_of_le I ha hmn
    simpa using (Ideal.neg_mem_iff _).2 h
  obtain ⟨e, he⟩ := IsPrecomplete.prec' (idemSeq a) aux
  -- ★ `e − a_n ∈ I^n`
  have hen : ∀ n : ℕ, e - idemSeq a n ∈ I ^ n := by
    intro n
    have h := he n
    rw [← Ideal.one_eq_top, Ideal.smul_eq_mul, mul_one, SModEq.sub_mem] at h
    simpa using (Ideal.neg_mem_iff _).2 h
  refine ⟨e, ?_, ?_⟩
  · -- ★ `e² = e`——Hausdorff
    have hzero : ∀ n : ℕ, e ^ 2 - e ≡ 0 [SMOD (I ^ n • ⊤ : Ideal R)] := by
      intro n
      rw [← Ideal.one_eq_top, Ideal.smul_eq_mul, mul_one, SModEq.zero]
      have h1 : (e - idemSeq a n) * (e + idemSeq a n - 1) ∈ I ^ n :=
        Ideal.mul_mem_right _ _ (hen n)
      have h2 : idemSeq a n ^ 2 - idemSeq a n ∈ I ^ n :=
        Ideal.pow_le_pow_right (Nat.le_of_lt Nat.lt_two_pow_self)
          (idemSeq_sq_sub_mem I ha n)
      have heq : e ^ 2 - e
          = (e - idemSeq a n) * (e + idemSeq a n - 1)
            + (idemSeq a n ^ 2 - idemSeq a n) := by ring
      rw [heq]
      exact Ideal.add_mem _ h1 h2
    have := IsHausdorff.haus' (I := I) (M := R) (e ^ 2 - e) hzero
    exact sub_eq_zero.mp this
  · -- ★ `e − a ∈ I`
    have h1 : e - idemSeq a 1 ∈ I := by
      have := hen 1
      simpa using this
    have h2 : idemSeq a 1 - a ∈ I := by
      have := idemSeq_succ_sub_mem I ha 0
      simpa using this
    have heq : e - a = (e - idemSeq a 1) + (idemSeq a 1 - a) := by ring
    rw [heq]
    exact Ideal.add_mem _ h1 h2

/-! ## ★出典の紐付け(`.src`) -/

def idemStep.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Newton の 1 段 f(a) = 3a² − 2a³)",
    sectionId := "genell-thm-3-8" }

def exists_isIdempotentElem_of_isAdicComplete.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(完備なイデアルに沿って冪等元が持ち上がる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_isIdempotentElem_of_isAdicComplete.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1357）**——第 1356 で特定した" ++
       "**mathlib に無い primitive** をここで作った。" ++
       "☆`RingTheory/Idempotents.lean` は冪零核に沿った持ち上げしか持たない。" ++
       "★これで「完備局所環の上の加群有限代数は局所環の直積」→" ++
       "「完備 DVR の有限整閉包は局所」→ DVR 構造、という道の第 1 段が通る。") 19 ]

end ABC3.Found.GenEll
